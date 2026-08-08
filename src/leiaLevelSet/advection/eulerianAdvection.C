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

#include "eulerianAdvection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(eulerianAdvection, false);
addToRunTimeSelectionTable(levelSetAdvection, eulerianAdvection, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eulerianAdvection::eulerianAdvection
(
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const velocityModel& velModel
)
:
    levelSetAdvection(mesh, U, phi, velModel),
    velExt_(velocityExtension::New(mesh)),
    source_(sdplsSource::New(mesh))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void eulerianAdvection::initialise(volScalarField& psi)
{
    // Preserved pre-loop side effects of the pre-unification solver:
    // register + write the SDPLS diagnostic fields and the initial
    // extension velocity.
    source_->update(psi, U_);
    source_->write();

    velExt_->correct();
    velExt_->Uext().write();
}


void eulerianAdvection::advance(volScalarField& psi)
{
    // Correct the advecting velocity to be interface-normal-constant in a
    // narrow band (non-invasive: U/phi unchanged; "none" = identity).
    velExt_->correct();

    // The level set obeys the ADVECTIVE (Hamilton-Jacobi) law Dpsi/Dt = 0,
    // not a conservation law: the conservative fvm::div(phiExt, psi) equals
    // (v.grad)psi only for solenoidal v. Assemble the advective derivative
    // exactly: (v.grad)psi = div(phiExt psi) - psi (div phiExt).
    const volScalarField divPhiExt("divPhiExt", fvc::div(velExt_->phi()));

    // Defect-correction loop for the deferred second-order (linearUpwind)
    // spatial term: each pass re-assembles with the latest psi so the
    // explicit (linearUpwind - upwind) correction converges (~2-3 passes ->
    // formal 2nd order). fvm::ddt reuses psi.oldTime(), fixed within the
    // step, so re-solving does not corrupt the time derivative.
    const label nDefCorr =
        mesh_.time().controlDict().getOrDefault<label>("nDefCorr", 3);
    for (label corr = 0; corr < nDefCorr; ++corr)
    {
        fvScalarMatrix psiEqn
        (
            fvm::ddt(psi)
            + fvm::div(velExt_->phi(), psi)
            - fvm::Sp(divPhiExt, psi)
        ==
            source_->fvmsdplsSource(psi, U_)
        );

        psiEqn.solve();
    }
}

} // End namespace Foam

// ************************************************************************* //
