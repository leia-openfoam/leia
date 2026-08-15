/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt
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
    Foam::narrowBand

Description
    Abstract base class which starts the model hierachy for marking cells which are near the interface aka the "narrowBand".
    Calculate a field with {0,1} which mark cells near the zero level-set interface.


SourceFiles
    narrowBand.C

\*---------------------------------------------------------------------------*/

#ifndef narrowBand_H
#define narrowBand_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class narrowBand Declaration
\*---------------------------------------------------------------------------*/

class narrowBand
{
protected:
    const fvMesh& mesh_;
    const volScalarField& psi_;

    //- psi with GUARANTEED-CURRENT coupled-patch values.
    //
    //  Every narrow-band criterion is a statement about psi across a face, and
    //  on a processor boundary the far side of that face comes from
    //  patchNeighbourField(), which returns whatever was last exchanged. If the
    //  band is built at a point where psi changed without a subsequent boundary
    //  update -- after an explicit source, a redistancing pass, an intermediate
    //  defect-correction iterate -- then the parallel run tests against a STALE
    //  neighbour value while the serial run always sees the current one. The
    //  band then depends on the DECOMPOSITION, and silently: no error, just a
    //  different set of cells on 16 ranks than on 1.
    //
    //  Deriving the criterion from this copy removes the dependence on whatever
    //  the caller happened to do beforehand. It is cheap next to the face loops
    //  that follow, and it is the only way to make the guarantee hold for a
    //  criterion nobody has written yet.
    tmp<volScalarField> psiSynced() const;


public:

    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        narrowBand, 
        Dictionary, 
        (
            const dictionary& dict, const volScalarField& psi
        ), 
        (dict, psi)
    );

    narrowBand(const dictionary& dict, const volScalarField& psi);

    static autoPtr<narrowBand> New(const fvMesh& mesh, const volScalarField& psi);

    virtual ~narrowBand() = default;


    // Member Functions

    const fvMesh& mesh() const;

    const volScalarField& psi() const;

    virtual void calc();

    virtual void write() const;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
