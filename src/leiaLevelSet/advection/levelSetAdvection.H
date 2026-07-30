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
    Foam::levelSetAdvection

Description
    Run-time-selectable level-set advection method: how psi^n becomes
    psi^{n+1}. The unified leiaLevelSetFoam solver composes ONE advection
    model with the (independent) redistancer, phase indicator and narrow
    band models -- the complete level-set method is configured from the
    fvSolution levelSet dictionary:

    \verbatim
    levelSet
    {
        advection
        {
            type eulerian;      // eulerian | semiLagrangian
        }
        ...
    }
    \endverbatim

    Models:
      - eulerian: implicit FV transport in the advective
        (conservative-minus-compression) form with a deferred-correction
        loop for second-order (linearUpwind) div schemes; composes the
        velocityExtension and sdplsSource models (both dictionary-selected,
        both defaulting to inactive).
      - semiLagrangian: characteristic (backward foot point) update with a
        runtime-selectable reconstruction (levelSet.semiLagrangian).

SourceFiles
    levelSetAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetAdvection_H
#define levelSetAdvection_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "velocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class levelSetAdvection Declaration
\*---------------------------------------------------------------------------*/

class levelSetAdvection
{
protected:

    const fvMesh& mesh_;
    const volVectorField& U_;
    const surfaceScalarField& phi_;
    const velocityModel& velModel_;

public:

    TypeName("levelSetAdvection");

    declareRunTimeSelectionTable
    (
        autoPtr,
        levelSetAdvection,
        Mesh,
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const velocityModel& velModel
        ),
        (mesh, U, phi, velModel)
    );

    // Constructors

        levelSetAdvection
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const velocityModel& velModel
        );

    // Selectors

        //- Select from fvSolution levelSet.advection.type (default eulerian
        //  -- every pre-existing case dictionary remains valid).
        static autoPtr<levelSetAdvection> New
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const velocityModel& velModel
        );

    //- Destructor
    virtual ~levelSetAdvection() = default;

    // Member Functions

    //- One-time setup needing the level set (called once after
    //  construction, before the time loop): register/write derived fields.
    virtual void initialise(volScalarField& psi)
    {}

    //- Advance psi^n (on entry) to psi^{n+1} (on exit).
    virtual void advance(volScalarField& psi) = 0;

    //- The flux that actually advects psi (CFL/deltaT control source);
    //  base flux unless a model substitutes its own (velocity extension).
    virtual const surfaceScalarField& advectingFlux() const
    {
        return phi_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
