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
    Foam::velocityExtension

Description
    Abstract base for runtime-selectable velocity *correction* models that, in
    a narrow band, replace the advecting velocity with one that is constant
    along the interface normal (an extension velocity). Selected by
    `levelSet { velocityExtension { type ...; } }` in fvSolution.

    The model is NON-INVASIVE: it reads the base velocity U (set by the
    velocityModel or, later, a Navier-Stokes solve) and emits a corrected field
    Uext (registered, AUTO_WRITE, written at U's frequency for visualization and
    error norms) plus an advection flux phiExt. The level-set equation advects
    with phiExt(). U and the base flux phi are never mutated.

    The default type "none" is the identity: Uext == U, phiExt == phi.

SourceFiles
    velocityExtension.C

\*---------------------------------------------------------------------------*/

#ifndef velocityExtension_H
#define velocityExtension_H

#include "fvMesh.H"
#include "fvSolution.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class velocityExtension Declaration
\*---------------------------------------------------------------------------*/

class velocityExtension
{
protected:

    const fvMesh& mesh_;
    const fvSolution& fvSolution_;
    const dictionary& levelSetDict_;
    //- Copy (not reference): subOrEmptyDict returns a temporary when the
    //  velocityExtension sub-dict is absent (default "none").
    dictionary velExtDict_;

    //- Base (uncorrected) velocity and flux, looked up from the registry.
    const volVectorField& U_;
    const surfaceScalarField& phi_;

    //- Corrected velocity (AUTO_WRITE) and the advection flux derived from it.
    volVectorField Uext_;
    surfaceScalarField phiExt_;

public:

    //- Runtime type information
    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        velocityExtension,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );

    // Constructors

        //- Construct from mesh; looks up U/phi and creates Uext/phiExt.
        velocityExtension(const fvMesh& mesh);

    // Selectors

        static autoPtr<velocityExtension> New(const fvMesh& mesh);

    //- Destructor
    virtual ~velocityExtension() = default;

    // Member Functions

        //- The corrected (extension) velocity field.
        const volVectorField& Uext() const { return Uext_; }

        //- The advecting flux for the level-set equation.
        const surfaceScalarField& phi() const { return phiExt_; }

        //- Update Uext/phiExt. Base "none": identity (Uext=U, phiExt=phi).
        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
