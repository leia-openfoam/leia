/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

Class
    Foam::velocityModel

Description
    Abstract base for runtime-selectable *prescribed* (verification) advection
    velocity fields, e.g. rigid rotation, deformation, vortex. Selected by
    `velocityModel { type ...; }` in fvSolution. Sets the cell-centred velocity
    U and the volumetric face flux phi from an analytic `velocity(x)`. This is
    the verification velocity; it will be replaced by a Navier-Stokes solve.

SourceFiles
    velocityModel.C
    prescribedVelocityModels.C

\*---------------------------------------------------------------------------*/

#ifndef velocityModel_H
#define velocityModel_H

#include "fvMesh.H"
#include "fvSolution.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class velocityModel Declaration
\*---------------------------------------------------------------------------*/

class velocityModel
{
protected:

    const fvSolution& fvSolution_;
    const dictionary& velocityDict_;
    Switch isOscillating_;
    const scalar tau_;

public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        velocityModel,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );

    // Constructors

        //- Construct from mesh (reads fvSolution velocityModel sub-dict)
        velocityModel(const fvMesh& mesh);

    // Selectors

        //- Select from fvSolution velocityModel.type
        static autoPtr<velocityModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~velocityModel() = default;

    // Member Functions

        Switch isOscillating() const
        {
            return isOscillating_;
        }

        scalar oscillationFactor(const scalar time, const scalar tau) const;

        //- Scale U and phi by cos(pi t / tau) for reversed (oscillating) tests
        void oscillateVelocity
        (
            volVectorField& U,
            const volVectorField& U0,
            surfaceScalarField& phi,
            const surfaceScalarField& phi0,
            const Time& runTime
        ) const;

        //- Analytic velocity at a point (implemented by concrete models)
        virtual vector velocity(const vector& x) const;

        //- Fill U from velocity(cellCentre) on internal + boundary fields
        virtual void setVelocity(volVectorField& U) const;

        //- Fill phi from velocity(faceCentre) & Sf on internal + boundary faces
        virtual void setVolumetricFlux(surfaceScalarField& phi) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
