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

#include "linearTaylorReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearTaylorReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        linearTaylorReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearTaylorReconstruction::linearTaylorReconstruction(const fvMesh& mesh)
:
    slReconstruction(mesh),
    gradPsi_
    (
        IOobject
        (
            "slGradPsi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("slGradPsi", dimless, vector::zero),
        "zeroGradient"
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearTaylorReconstruction::update(const volScalarField& psiOld)
{
    // Cell-point-cell LSQ gradient (linear-exact); the fvSchemes key "gradPsi"
    // resolves the pointCellsLeastSquares scheme. fvc::grad gathers the halo
    // internally, so gradPsi_[c] is correct even where c's stencil crosses a
    // processor boundary. Assign values only (primitiveFieldRef): psi may carry
    // any dimensions (it is a signed distance here) -- we keep gradPsi_ a plain
    // numeric field, so a dimensioned operator= would abort.
    const volVectorField g(fvc::grad(psiOld, "gradPsi"));
    gradPsi_.primitiveFieldRef() = g.primitiveField();
    gradPsi_.correctBoundaryConditions();

    // Collect the stencil (psiOld + centres) so slAdvection can apply the
    // foot-radius guard and the optional monotonicity clip; also sets
    // psiOldPtr_. (fvc::grad already did the gradient's own halo exchange.)
    collectStencil(psiOld);
    computeLimiters();
}


Foam::scalar Foam::linearTaylorReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const vector d = x - mesh_.C()[c];
    return (*psiOldPtr_)[c] + (gradPsi_[c] & d);
}

// ************************************************************************* //
