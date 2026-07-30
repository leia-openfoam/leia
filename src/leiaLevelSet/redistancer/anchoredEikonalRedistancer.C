/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "anchoredEikonalRedistancer.H"
#include "planeAnchors.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(anchoredEikonalRedistancer, false);
addToRunTimeSelectionTable(redistancer, anchoredEikonalRedistancer, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

anchoredEikonalRedistancer::anchoredEikonalRedistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    mesh_(mesh),
    epsilon_(redistDict_.getOrDefault<scalar>("epsilon", 0.1)),
    tolerance_(redistDict_.getOrDefault<scalar>("tolerance", 1e-3)),
    maxIter_(redistDict_.getOrDefault<label>("maxIter", 10))
{
    Info<< "anchoredEikonalRedistancer: epsilon = " << epsilon_
        << ", tolerance = " << tolerance_
        << ", maxIter = " << maxIter_ << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool anchoredEikonalRedistancer::doRedistance(volScalarField& psi)
{
    if (!mesh_.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "anchoredEikonal redistancing requires the registered "
            << "'NarrowBand' field: set levelSet.narrowBand.type to "
            << "signChange (or a type derived from it) in fvSolution."
            << exit(FatalError);
    }

    const volScalarField& narrowBand =
        mesh_.lookupObject<volScalarField>("NarrowBand");

    const planeAnchorData anchors =
        computePlaneAnchors
        (
            mesh_, psi, narrowBand,
            redistDict_.getOrDefault<bool>("curvatureCorrection", true)
        );

    if (returnReduce(anchors.cells.size(), sumOp<label>()) == 0)
    {
        Info<< "anchoredEikonalRedistancer: no interface anchors, skipping."
            << endl;
        return false;
    }

    // Re-sign source, captured before psi is touched.
    const scalarField psiOld(psi.primitiveField());

    // Unsigned working distance; copy-construct-with-rename inherits psi's
    // boundary conditions (zeroGradient/empty/processor). The names
    // "psiDist" and "yPhi" are load-bearing: fvSchemes/fvSolution lookups
    // are keyed on them.
    volScalarField y("psiDist", psi);
    y == mag(psi);

    const scalarField anchorAbs(mag(anchors.dist));

    // Explicit solver dict: a bare solve() would look up "psiDistFinal" on
    // final PIMPLE iterations (pimpleControl sets finalIteration already on
    // the first corrector when nOuterCorrectors is 1) and FatalError.
    const dictionary& solverDict =
        fvSolution_.subDict("solvers").subDict("psiDist");

    // Tucker regularized steady Eikonal, transport direction frozen per
    // iteration (mirrors OpenFOAM advectionDiffusionPatchDistMethod).
    label iter = 0;
    scalar initialResidual = 0;
    do
    {
        volVectorField ny(fvc::grad(y));
        ny /= (mag(ny) + SMALL);

        surfaceVectorField nf(fvc::interpolate(ny));
        nf /= (mag(nf) + SMALL);

        surfaceScalarField yPhi("yPhi", mesh_.Sf() & nf);

        fvScalarMatrix yEqn
        (
            fvm::div(yPhi, y)
          + fvm::SuSp(-fvc::div(yPhi), y)
          - epsilon_*y*fvm::laplacian(y)
         ==
            dimensionedScalar("1", dimless, 1.0)
        );

        // OF's constrain order: relax first, then impose the internal
        // Dirichlet rows. setValues also writes the constrained values into
        // y, so the anchors feed the next iterate's gradient.
        yEqn.relax();
        yEqn.setValues(anchors.cells, anchorAbs);

        initialResidual = yEqn.solve(solverDict).initialResidual();

    } while (initialResidual > tolerance_ && ++iter < maxIter_);

    // Guard against slightly negative distances flipping signs below.
    y.clamp_min(SMALL);

    // Re-sign: sign pattern from the pre-redistance field; anchored cells
    // receive their SIGNED plane distances exactly (the plane may
    // legitimately realign a cut-cell centre sub-cell -- callers refresh
    // the narrow band after an event).
    scalarField& psiIn = psi.primitiveFieldRef();
    const scalarField& yIn = y.primitiveField();
    forAll(psiIn, celli)
    {
        psiIn[celli] = (psiOld[celli] >= 0) ? yIn[celli] : -yIn[celli];
    }
    forAll(anchors.cells, i)
    {
        psiIn[anchors.cells[i]] = anchors.dist[i];
    }
    psi.correctBoundaryConditions();

    Info<< "anchoredEikonalRedistancer: " << anchors.cells.size()
        << " anchors (" << anchors.bandCells.size() << " band), "
        << iter + 1 << " iterations, final initial residual "
        << initialResidual << endl;

    return true;
}

} // End namespace Foam

// ************************************************************************* //
