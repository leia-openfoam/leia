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

Class
    Foam::semiLagrangianAdvection

Description
    Characteristic (semi-Lagrangian) level-set advection:
    psi^{n+1}(x_c) = psi^n(x_d) with the second-order Taylor backward foot
    x_d and a runtime-selectable reconstruction of psi^n at the foot
    (levelSet.semiLagrangian sub-dictionary; see slAdvection). Owns the
    previous-step velocity u^n it reconstructs from the prescribed velocity
    model's oscillation factor (reversed tests) or the steady initial
    velocity.

SourceFiles
    semiLagrangianAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef semiLagrangianAdvection_H
#define semiLagrangianAdvection_H

#include "levelSetAdvection.H"
#include "slAdvection.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class semiLagrangianAdvection Declaration
\*---------------------------------------------------------------------------*/

class semiLagrangianAdvection
:
    public levelSetAdvection
{
private:

    autoPtr<slAdvection> slAdv_;

    //- Initial (t = 0) velocity: the oscillation-factor base for u^n
    //  (private copy -- the solver's U0 error-reference field is
    //  constructed later and is not this model's concern).
    const volVectorField U0_;

    //- Previous-step velocity u^n the Taylor trajectory needs
    volVectorField Uold_;

public:

    TypeName("semiLagrangian");

    semiLagrangianAdvection
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const surfaceScalarField& phi,
        const velocityModel& velModel
    );

    virtual ~semiLagrangianAdvection() = default;

    virtual void advance(volScalarField& psi);
};

} // End namespace Foam

#endif

// ************************************************************************* //
