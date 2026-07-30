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

#include "semiLagrangianAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(semiLagrangianAdvection, false);
addToRunTimeSelectionTable(levelSetAdvection, semiLagrangianAdvection, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

semiLagrangianAdvection::semiLagrangianAdvection
(
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const velocityModel& velModel
)
:
    levelSetAdvection(mesh, U, phi, velModel),
    slAdv_(slAdvection::New(mesh)),
    U0_("U0.advection", U),
    Uold_("Uold.advection", U)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void semiLagrangianAdvection::advance(volScalarField& psi)
{
    // u^n: the previous-step velocity the Taylor trajectory needs. For the
    // oscillating (reversed) field it is U0 scaled at t^n = t^{n+1} - dt;
    // for the steady field u^n = u^{n+1} = U0. (The solver has already
    // rescaled U/phi to u^{n+1} via oscillateVelocity.)
    if (velModel_.isOscillating())
    {
        const scalar tn =
            mesh_.time().value() - mesh_.time().deltaT().value();
        Uold_ == U0_*velModel_.oscillationFactor(tn, velModel_.tau());
    }
    else
    {
        Uold_ == U0_;
    }

    // Semi-Lagrangian update: psi holds psi^n on entry, psi^{n+1} on exit.
    slAdv_->advect(psi, U_, Uold_);
}

} // End namespace Foam

// ************************************************************************* //
