/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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

#include "anisotropicDiffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicDiffusion, 0);
    addToRunTimeSelectionTable(velocityExtension, anisotropicDiffusion, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicDiffusion::anisotropicDiffusion(const fvMesh& mesh)
:
    interfaceExtension(mesh),
    epsilon_(velExtDict_.getOrDefault<scalar>("epsilon", 1e-3))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::anisotropicDiffusion::extend()
{
    // Normal-aligned (anisotropic) diffusion tensor with isotropic
    // regularization: D = n (x) n + eps I.
    volSymmTensorField D
    (
        IOobject
        (
            "velExtD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(nHat_)
      + dimensionedSymmTensor("eps", dimless, epsilon_*symmTensor::I)
    );

    fvVectorMatrix UextEqn(fvm::laplacian(D, Uext_));

    // Whole-domain harmonic extension: Dirichlet only at the single interface
    // layer (the interpolated interface velocity); natural zeroGradient BCs at
    // the domain boundaries (Uext carries zeroGradient patches).
    DynamicList<label> fixedCells(seedCells_.size());
    DynamicList<vector> fixedVals(seedCells_.size());
    forAll(seedCells_, i)
    {
        fixedCells.append(seedCells_[i]);
        fixedVals.append(Uext_[seedCells_[i]]);
    }

    UextEqn.setValues(fixedCells, fixedVals);

    // Linear-solver controls come from the case dictionary (fvSolution ->
    // solvers -> "Uext.*"), not the source -- e.g. PBiCGStab/DIC with a
    // per-step RELATIVE tolerance (the previous step's Uext is the initial
    // guess); demanding a tight absolute tolerance on this near-rank-1
    // anisotropic Laplacian costs hundreds of iterations for no accuracy gain.
    UextEqn.solve();
}

// ************************************************************************* //
