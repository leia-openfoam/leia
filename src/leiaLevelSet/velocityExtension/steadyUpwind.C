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

#include "steadyUpwind.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyUpwind, 0);
    addToRunTimeSelectionTable(velocityExtension, steadyUpwind, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyUpwind::steadyUpwind(const fvMesh& mesh)
:
    interfaceExtension(mesh),
    nDefCorrExt_(velExtDict_.getOrDefault<label>("nDefCorrExt", 5)),
    zeroInflowGuard_(velExtDict_.getOrDefault<scalar>("zeroInflowGuard", 1e-6)),
    defCorrRelax_(velExtDict_.getOrDefault<scalar>("defCorrRelax", 0.7))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::steadyUpwind::solveSteady
(
    volVectorField& q,
    const labelUList& fixedCells,
    const UList<vector>& fixedVals,
    const volVectorField* source
)
{
    // Outward propagation direction w = sign(psi) n: NOT solenoidal, so the
    // advective derivative is assembled as the conservative flux divergence
    // minus the compression term -- the exact steady twin of the solver's
    // psi-equation treatment.
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

    const surfaceScalarField phiW
    (
        "velExtPhiW",
        fvc::interpolate(w) & mesh_.Sf()
    );
    const volScalarField divW("velExtDivW", fvc::div(phiW));

    // Zero-inflow guard: cells at the discrete skeleton (local maxima of
    // |psi|) can have no upwind neighbour -> a (near-)zero matrix row. A tiny
    // relaxation towards the current value keeps the row regular without
    // measurably perturbing the transported solution.
    volScalarField lambda
    (
        IOobject
        (
            "velExtLambda",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("lambda", dimless/dimLength, 0)
    );
    forAll(lambda, c)
    {
        lambda[c] = zeroInflowGuard_/cellSize_[c];
    }

    // Defect-correction loop: the divergence scheme comes from fvSchemes
    // (div(velExtPhiW,Uext), Gauss linearUpwind grad(Uext) by default); its
    // implicit part is the upwind matrix (block-triangularizable along the
    // acyclic characteristic graph -> one solve reaches the whole domain);
    // each pass re-assembles the explicit deferred correction with the
    // latest q (~2-3 passes -> formal second order).
    // Steady deferred correction needs under-relaxation between passes: the
    // explicit (higher-order minus upwind) flux of pass k feeds pass k+1, and
    // without damping the fixed-point iteration rings on a pure-advection
    // steady system (no ddt diagonal).
    vectorField qPrev(q.primitiveField());

    for (label corr = 0; corr < nDefCorrExt_; ++corr)
    {
        // The lambda*q RHS relaxes guarded (zero-inflow) rows towards their
        // current value instead of dragging them to zero.
        fvVectorMatrix qEqn
        (
            fvm::div(phiW, q)
          - fvm::Sp(divW, q)
          + fvm::Sp(lambda, q)
         ==
            lambda*q
        );

        if (source)
        {
            qEqn -= *source;   // explicit RHS (Aslam cascade stage)
        }

        qEqn.setValues(fixedCells, fixedVals);
        // Linear-solver controls come from fvSolution -> solvers -> "Uext.*".
        qEqn.solve();

        if (corr < nDefCorrExt_ - 1)
        {
            q.primitiveFieldRef() =
                defCorrRelax_*q.primitiveField()
              + (1.0 - defCorrRelax_)*qPrev;
            qPrev = q.primitiveField();
            q.correctBoundaryConditions();
        }
    }
}


void Foam::steadyUpwind::extend()
{
    DynamicList<label> fixedCells(seedCells_.size());
    DynamicList<vector> fixedVals(seedCells_.size());
    forAll(seedCells_, i)
    {
        fixedCells.append(seedCells_[i]);
        fixedVals.append(Uext_[seedCells_[i]]);
    }

    solveSteady(Uext_, fixedCells, fixedVals);
}

// ************************************************************************* //
