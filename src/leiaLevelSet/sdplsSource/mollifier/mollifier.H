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
    Foam::mollifier

Description
    Base class for implementing mollfier / cut-off fields for SDPLS source terms.

    This class implements the mollifier 'none' which does not scale the Level
    Set field.

SourceFiles
    mollifier.C

\*---------------------------------------------------------------------------*/

#ifndef mollifier_H
#define mollifier_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class mollifier1 Declaration
\*---------------------------------------------------------------------------*/

class mollifier
{
public:

    // Static Data Members
    
    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        mollifier, 
        Dictionary, 
        (
            const dictionary& dict
        ), 
        (dict)
    );

    // Constructors

        //- Construct from components
        mollifier(const dictionary& dict);

    // Selectors

        //- Select default constructed
        static autoPtr<mollifier> New(const dictionary& dict);

    //- Destructor
    virtual ~mollifier() = default;

    // Member functions

        //- Calculates the mollifier field
        //
        //  VIRTUAL because not every cut-off is a pointwise function of psi.
        //  The default implementation evaluates mollify(psi[cellI]) per cell,
        //  which presupposes that |psi| is a usable DISTANCE from the interface.
        //  That presupposition fails exactly where it matters: on the 3D LeVeque
        //  deformation |grad psi| departs far from 1, so a fixed |psi| threshold
        //  no longer marks a fixed distance and the cut-off clips the source
        //  where it is still needed -- measured, the band error then DIVERGES
        //  (order -0.54 against -0.02 unmollified). A topology-based band
        //  overrides this member instead.
        //
        //  @param psi is the Level Set field
        virtual tmp<volScalarField> field(const volScalarField& psi) const;

protected:

    // Protected Member functions
        
        //- mollifier formula for a cell
        //
        //  @param x is the Level Set value
        virtual double mollify(double x) const;


};

} // End namespace Foam

#endif

// ************************************************************************* //
