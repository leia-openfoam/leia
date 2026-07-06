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

#include "pseudoTime.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "gaussConvectionScheme.H"
#include "upwind.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudoTime, 0);
    addToRunTimeSelectionTable(velocityExtension, pseudoTime, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudoTime::pseudoTime(const fvMesh& mesh)
:
    interfaceExtension(mesh),
    nIterations_(velExtDict_.getOrDefault<label>("nIterations", 5)),
    deltaTau_(velExtDict_.getOrDefault<scalar>("deltaTau", 0.5))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pseudoTime::extend()
{
    // Outward propagation direction w = sign(psi) n (points away from the
    // interface on both sides); dimensionless.
    volVectorField w
    (
        IOobject
        (
            "velExtW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimless, vector::zero),
        "zeroGradient"
    );
    forAll(w, c)
    {
        w[c] = ((psi_[c] < 0) ? -1.0 : 1.0)*nHat_[c];
    }
    w.correctBoundaryConditions();

    const surfaceScalarField phiW("velExtPhiW", fvc::interpolate(w) & mesh_.Sf());

    // Dirichlet only at the single interface layer (the interpolated interface
    // velocity); whole-domain solve with natural zeroGradient BCs elsewhere.
    // The pseudo-time reach is limited by nIterations, so the correction stays
    // near the interface even without band confinement.
    DynamicList<label> fixedCells(seedCells_.size());
    DynamicList<vector> fixedVals(seedCells_.size());
    forAll(seedCells_, i)
    {
        fixedCells.append(seedCells_[i]);
        fixedVals.append(Uext_[seedCells_[i]]);
    }

    // Inverse pseudo-time step 1/(deltaTau*cellSize) [1/length], so the
    // implicit Sp term and the convection matrix are dimensionally consistent.
    volScalarField rDtau
    (
        IOobject
        (
            "velExtRDtau",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rDtau", dimless/dimLength, 0)
    );
    forAll(rDtau, c)
    {
        rDtau[c] = 1.0/(deltaTau_*Foam::pow(mesh_.V()[c], 1.0/3.0));
    }

    // Implicit backward-Euler pseudo-time march of the NON-CONSERVATIVE normal
    // advection (Adalsteinsson & Sethian velocity extension):
    //   dUext/dtau + (w . grad)Uext = 0,   w = sign(psi) n  (upwind).
    // (w . grad)Uext is assembled as div(phiW,Uext) - Uext (div phiW): the
    // conservative convection matrix minus the flux divergence, giving the
    // advective (non-conservative) derivative -- which is what an extension
    // needs (it transports Uext along w, it does not conserve phiW*Uext, and
    // w = sign(psi) n is NOT divergence-free). The rDtau Sp term gives a
    // positive diagonal; seed + band exterior are matrix Dirichlet constraints.
    const volScalarField divW("velExtDivW", fvc::div(phiW));
    for (label it = 0; it < nIterations_; ++it)
    {
        fvVectorMatrix UextEqn
        (
            fvm::Sp(rDtau, Uext_)
          + fv::gaussConvectionScheme<vector>
            (
                mesh_,
                phiW,
                tmp<surfaceInterpolationScheme<vector>>
                (
                    new upwind<vector>(mesh_, phiW)
                )
            ).fvmDiv(phiW, Uext_)
          - fvm::Sp(divW, Uext_)
         ==
            rDtau*Uext_
        );

        UextEqn.setValues(fixedCells, fixedVals);
        // Linear-solver controls come from the case dictionary (fvSolution ->
        // solvers -> "Uext.*"), not the source.
        UextEqn.solve();
    }
}

// ************************************************************************* //
