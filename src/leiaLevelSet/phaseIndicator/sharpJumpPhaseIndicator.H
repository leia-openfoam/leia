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
    Foam::sharpJumpPhaseIndicator

Description
    A phase indicator with an abrupt or sharp jump.

SourceFiles
    sharpJumpPhaseIndicator.C

\*---------------------------------------------------------------------------*/

#ifndef sharpJumpPhaseIndicator_H
#define sharpJumpPhaseIndicator_H

#include "phaseIndicator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sharpJumpPhaseIndicator Declaration
\*---------------------------------------------------------------------------*/

class sharpJumpPhaseIndicator
:
    public phaseIndicator
{
    label nAverages_; 

public:

    // Static Data Members

    TypeName ("sharpJump");

    // Constructors

    //- Construct from components
    sharpJumpPhaseIndicator(const fvMesh& mesh);

    //- Destructor
    virtual ~sharpJumpPhaseIndicator() = default;

    // Member Functions
    
    virtual void calcPhaseIndicator
    (
        volScalarField& alpha, 
        const volScalarField& psi 
    ); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
