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

#include "steadyUpwindLinear.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyUpwindLinear, 0);
    addToRunTimeSelectionTable(velocityExtension, steadyUpwindLinear, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyUpwindLinear::steadyUpwindLinear(const fvMesh& mesh)
:
    steadyUpwind(mesh)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::steadyUpwindLinear::extend()
{
    // Stage 1: the normal derivative of the physical velocity, per cell, from
    // the least-squares gradient (fvSchemes: grad(velExtLSGrad) entry); named
    // "Uext" (UNREGISTERED) so the steady transport picks up the same
    // div(velExtPhiW,Uext) scheme -- the phiExt=="phi" naming trick.
    // Constructed with explicit zeroGradient BCs (an expression-built field
    // would carry non-solvable 'calculated' patches -- fvm::div needs
    // evaluable boundary coefficients).
    volVectorField V
    (
        IOobject
        (
            "Uext",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false   // do not register (the real Uext is registered)
        ),
        mesh_,
        dimensionedVector("0", dimVelocity/dimLength, vector::zero),
        "zeroGradient"
    );
    V.primitiveFieldRef() =
        (nHat_ & fvc::grad(U_, "velExtLSGrad"))().primitiveField();
    V.correctBoundaryConditions();

    // Seeds for BOTH stages come from the interface layer.
    DynamicList<label> fixedCells(seedCells_.size());
    DynamicList<vector> fixedValsV(seedCells_.size());
    DynamicList<vector> fixedValsU(seedCells_.size());
    forAll(seedCells_, i)
    {
        const label c = seedCells_[i];
        fixedCells.append(c);
        fixedValsV.append(V[c]);
        fixedValsU.append(Uext_[c]);
    }

    // Stage 2: extend V constant along normals.
    solveSteady(V, fixedCells, fixedValsV);

    // Stage 3: extend U with the extended normal derivative as source:
    //   sign(psi) (n.grad) Uext = sign(psi) V_ext
    // => Uext varies LINEARLY along each normal ray.
    volScalarField signPsi
    (
        IOobject
        (
            "velExtSign",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("s", dimless, 1)
    );
    forAll(signPsi, c)
    {
        signPsi[c] = (psi_[c] < 0) ? -1.0 : 1.0;
    }

    const volVectorField S("velExtSrc", signPsi*V);

    solveSteady(Uext_, fixedCells, fixedValsU, &S);
}

// ************************************************************************* //
