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
    Foam::slSource

Description
    A source step of the semi-Lagrangian level-set transport,
    D_t psi = psi F, applied by slAdvection after the characteristic update:
    psi holds the arrival value psi^A on entry and psi^{n+1} on exit.

    Selected by

        levelSet { semiLagrangian { source { type <model>; ... } } }

    The base class is the inert model `none`: active() is false, slAdvection
    does not call apply(), and every run without a source is byte-identical to
    a run before this family existed (2026-09-26). New source models are new
    classes in this family (CLAUDE.md "Extension without modification").

SourceFiles
    slSource.C

\*---------------------------------------------------------------------------*/

#ifndef slSource_H
#define slSource_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class slReconstruction;

class slSource
{
protected:

    const fvMesh& mesh_;

    //- levelSet.semiLagrangian.source (empty if the entry is absent).
    const dictionary dict_;

public:

    //- Runtime type information
    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slSource,
        Mesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        explicit slSource(const fvMesh& mesh);

    // Selectors

        static autoPtr<slSource> New(const fvMesh& mesh);

    //- Destructor
    virtual ~slSource() = default;

    // Member Functions

        //- False for `none`: slAdvection then skips apply().
        virtual bool active() const { return false; }

        //- psi holds the arrival value psi^A on entry and psi^{n+1} on exit.
        //  Utrajectory is the velocity the characteristics used at t^{n+1};
        //  geometryFit the fit of the geometry path (not yet refitted on psi^A);
        //  dt the time step this advect() covers (Stormer-Verlet halves included).
        virtual void apply
        (
            volScalarField& psi,
            const volVectorField& Utrajectory,
            slReconstruction& geometryFit,
            const scalar dt
        )
        {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
