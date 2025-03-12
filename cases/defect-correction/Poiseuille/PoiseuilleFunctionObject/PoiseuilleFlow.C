/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Tomislav Maric, TU Darmstadt 
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
    Foam::PoiseuilleFlow

Description

    A concrete class implementing the Poiseuille flow, i.e. the flow 
    between two noSlip walls.

SourceFiles
    PoiseuilleFlow.C

\*---------------------------------------------------------------------------*/

#include "PoiseuilleFlow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

PoiseuilleFlow::PoiseuilleFlow
(
    const scalar H, 
    const scalar Vis,
    const scalar L,
    const scalar& P
)
: 
    H_(H),
    Vis_(Vis),
    L_(L),
    P_(P)
{}

PoiseuilleFlow::PoiseuilleFlow(const dictionary& dict)
: 
    H_(dict.getOrDefault<scalar>("H", 1.0)),
    P_(dict.getOrDefault<scalar>("P", 100.0)),
    Vis_(dict.getOrDefault<scalar>("Vis", 10.0)),
    L_(dict.getOrDefault<scalar>("L", 1.0))
{}


vector PoiseuilleFlow::velocityCartesian(const vector& pos) const
{
    // Compute distance of y from the axis point.
    scalar r = mag(pos);
    if (r == 0) // then undefined transformation at the origin.
    {
       return vector(0,0,0); 
    }

    // Compute velocity at y point
    scalar Uypos = 1/(2*Vis_)*(P_/L_)*(pos[1]*H_-sqr(pos[1]));

    return vector(Uypos, 0, 0);
}

scalar PoiseuilleFlow::pressureCartesian(const vector& pos) const
{
    // TODO(TM)
    return 0; 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

