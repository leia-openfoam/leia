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
    Foam::eulerianAdvection

Description
    Implicit finite-volume level-set transport in the ADVECTIVE
    (Hamilton-Jacobi) form,

    \f[
        \partial_t\psi + \nabla\cdot(\vec{v}\,\psi) - \psi\nabla\cdot\vec{v}
        = S_{SDPLS}(\psi),
    \f]

    with a DEFERRED-CORRECTION loop (controlDict nDefCorr, default 2)
    re-assembling the equation so the explicit part of a second-order
    deferred div scheme (linearUpwind) converges within the time step.
    Composes the runtime-selected velocityExtension (advecting flux
    phiExt; "none" = identity) and sdplsSource (RHS source vanishing on
    the zero level set; "noSource" = off) models.

SourceFiles
    eulerianAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef eulerianAdvection_H
#define eulerianAdvection_H

#include "levelSetAdvection.H"
#include "velocityExtension.H"
#include "sdplsSource.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class eulerianAdvection Declaration
\*---------------------------------------------------------------------------*/

class eulerianAdvection
:
    public levelSetAdvection
{
private:

    autoPtr<velocityExtension> velExt_;
    autoPtr<sdplsSource> source_;

public:

    TypeName("eulerian");

    eulerianAdvection
    (
        const fvMesh& mesh,
        const volVectorField& U,
        const surfaceScalarField& phi,
        const velocityModel& velModel
    );

    virtual ~eulerianAdvection() = default;

    virtual void initialise(volScalarField& psi);

    virtual void advance(volScalarField& psi);

    virtual const surfaceScalarField& advectingFlux() const
    {
        return velExt_->phi();
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
