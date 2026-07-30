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

#include "quadraticTaylorReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadraticTaylorReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        quadraticTaylorReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadraticTaylorReconstruction::quadraticTaylorReconstruction(const fvMesh& mesh)
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
    ),
    hessPsi_
    (
        IOobject
        (
            "slHessPsi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("slHessPsi", dimless, tensor::zero),
        "zeroGradient"
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::quadraticTaylorReconstruction::update(const volScalarField& psiOld)
{
    // First LSQ pass: cell gradient (key "gradPsi"). Assign values only
    // (primitiveFieldRef) so psi's dimensions (it is a signed distance) do not
    // make the plain-numeric helper field's operator= abort.
    const volVectorField g(fvc::grad(psiOld, "gradPsi"));
    gradPsi_.primitiveFieldRef() = g.primitiveField();
    gradPsi_.correctBoundaryConditions();

    // Second LSQ pass on the gradient field: the Jacobian d(g_i)/dx_j (key
    // "gradGradPsi"), symmetrised to the Hessian. Both passes are LSQ =>
    // exact for quadratic psi. (This is why fvc::grad(fvc::grad(psi)) with
    // DEFAULT schemes is rejected: the default would resolve to Gauss.)
    const volTensorField T(fvc::grad(gradPsi_, "gradGradPsi"));
    const volTensorField Tsym(0.5*(T + T.T()));
    hessPsi_.primitiveFieldRef() = Tsym.primitiveField();
    hessPsi_.correctBoundaryConditions();

    // Collect the stencil (psiOld + centres) for slAdvection's foot-radius
    // guard and the optional monotonicity clip; also sets psiOldPtr_.
    collectStencil(psiOld);
    computeLimiters();
}


Foam::scalar Foam::quadraticTaylorReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const vector d = x - mesh_.C()[c];
    return
        (*psiOldPtr_)[c]
      + (gradPsi_[c] & d)
      + 0.5*(d & (hessPsi_[c] & d));
}

// ************************************************************************* //
