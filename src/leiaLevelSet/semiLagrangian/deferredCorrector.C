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

#include "deferredCorrector.H"
#include "addToRunTimeSelectionTable.H"
#include "slReconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deferredCorrector, 0);
    addToRunTimeSelectionTable(slCorrector, deferredCorrector, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deferredCorrector::deferredCorrector
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    slCorrector(mesh, dict),
    nIters_(Foam::max(label(1), dict.getOrDefault<label>("nDefCorr", 3))),
    maxIters_(Foam::max(nIters_, dict.getOrDefault<label>("nDefCorrMax", 10))),
    relax_(Foam::min(scalar(1), dict.getOrDefault<scalar>("defCorrRelax", 1.0))),
    tol_(dict.getOrDefault<scalar>("defCorrTol", 0.0))
{
    Info<< "deferredCorrector: nDefCorr = " << nIters_
        << ", nDefCorrMax = " << maxIters_
        << ", defCorrRelax = " << relax_
        << ", defCorrTol = " << tol_ << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::deferredCorrector::correct
(
    volScalarField& psi,
    const pointField& feet,
    slReconstruction& recon
) const
{
    scalarField newPsi(mesh_.nCells());
    label nNonFinite = 0;
    scalar deltaRel = 0;
    label nPasses = 0;

    for (label it = 0; it < maxIters_; ++it)
    {
        ++nPasses;

        // Rebuild the reconstruction from the current iterate (pass 0: psi^n),
        // so any gradient it uses is consistent with the reconstructed field.
        recon.update(psi);
        const slReconstruction& R = recon;
        const bool clip = R.clipToStencilBounds();

        if (it == 0)
        {
            footRadiusGuard(R, feet);        // geometry only -> once
        }

        nNonFinite = 0;                      // report the delivered pass's quality
        scalar deltaMax = 0;
        scalar rangeMax = SMALL;

        forAll(feet, c)
        {
            const scalar v0 = robustEvaluate(R, c, feet[c], clip, nNonFinite);
            const scalar vOld = psi[c];
            const scalar v = vOld + relax_*(v0 - vOld);   // under-relaxed update
            scalar lo, hi;
            R.stencilRange(c, lo, hi);
            deltaMax = Foam::max(deltaMax, Foam::mag(v - vOld));
            rangeMax = Foam::max(rangeMax, hi - lo);
            newPsi[c] = v;
        }

        psi.primitiveFieldRef() = newPsi;
        psi.correctBoundaryConditions();

        reduce(deltaMax, maxOp<scalar>());
        reduce(rangeMax, maxOp<scalar>());
        deltaRel = deltaMax/Foam::max(rangeMax, SMALL);

        // Stop once the requested pass count is reached and (optionally) the
        // range-relative update has fallen below tol.
        if (it + 1 >= nIters_ && (tol_ <= 0 || deltaRel <= tol_))
        {
            break;
        }
    }

    warnNonFinite(nNonFinite);

    if (mesh_.time().writeTime())
    {
        Info<< "deferredCorrector: " << nPasses
            << " pass(es), final relative update = " << deltaRel << endl;
    }
}

// ************************************************************************* //
