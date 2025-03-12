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
    Foam::CouetteFlow

Description

    A concrete class implementing the Couette flow, i.e. the flow 
    between two rotating cylinders.

SourceFiles
    CouetteFlow.C

\*---------------------------------------------------------------------------*/

#include "CouetteFlow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

CouetteFlow::CouetteFlow
(
    const scalar H, 
    const scalar& U
)
: 
    H_(H),
    U_(U)
{}

CouetteFlow::CouetteFlow(const dictionary& dict)
: 
    H_(dict.getOrDefault<scalar>("H", 1.0)),
    U_(dict.getOrDefault<scalar>("U", 1.0))
{}


vector CouetteFlow::velocityCartesian(const vector& pos) const
{
    // Compute distance of y from the axis point.
    scalar r = mag(pos);
    if (r == 0) // then undefined transformation at the origin.
    {
       return vector(0,0,0); 
    }

    // Compute velocity at y point
    scalar Uypos = U_*(pos[1] / H_);

    return vector(Uypos, 0, 0);
}

scalar CouetteFlow::pressureCartesian(const vector& pos) const
{
    // TODO(TM)
    return 0; 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

