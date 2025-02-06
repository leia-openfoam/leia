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
    Foam::taylorCouetteFlow

Description

    A concrete class implementing the Taylor-Couette flow, i.e. the flow 
    between two rotating cylinders.

SourceFiles
    taylorCouetteFlow.C

\*---------------------------------------------------------------------------*/

#include "taylorCouetteFlow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

taylorCouetteFlow::taylorCouetteFlow
(
    const scalar Ri, 
    const scalar Ro, 
    const scalar& omegai, 
    const scalar& omegao, 
    const vector& axisPoint
)
: 
    Ri_(Ri),
    sqrRi_(sqr(Ri_)),
    Ro_(Ro),
    sqrRo_(sqr(Ro_)),
    omegai_(omegai),
    omegao_(omegao),
    axisPoint_(axisPoint)
{}

taylorCouetteFlow::taylorCouetteFlow(const dictionary& dict)
: 
    Ri_(dict.getOrDefault<scalar>("Ri", 1.0)),
    sqrRi_(sqr(Ri_)),
    Ro_(dict.getOrDefault<scalar>("Ro", 2.0)),
    sqrRo_(sqr(Ro_)),
    omegai_(dict.getOrDefault<scalar>("omegai", 1)),
    omegao_(dict.getOrDefault<scalar>("omegao", 0)),
    axisPoint_(dict.getOrDefault<vector>("axisPoint", vector(0,0,0)))
{}

vector taylorCouetteFlow::velocityCylindrical(const scalar r) const
{
    // Brenn, G. (2016). Analytical solutions for transport processes. New York: Springer.
    scalar u_theta = r*(omegao_*sqrRo_ - omegai_*sqrRi_) / (sqrRo_ - sqrRi_) - 
    	((omegao_ - omegai_) * sqrRi_*sqrRo_)/(r * (sqrRo_ - sqrRi_));

    return vector(0, u_theta, 0); 
}

vector taylorCouetteFlow::velocityCartesian(const vector& x) const
{
    // Compute radial distance of x from the axis point.
    scalar r = mag(x - axisPoint_);

    if (r < SMALL) // undefined transformation at the origin
    {
       return vector(0,0,0); 
    }

    // Compute cylindrical velocity 
    const vector uCylindrical = velocityCylindrical(r);

    // Compute trigonometric terms
    scalar cosTheta = x[0] / r;
    scalar sinTheta = x[1] / r;

    // Convert to Cartesian components
    scalar u_x = -uCylindrical[1] * sinTheta;
    scalar u_y = uCylindrical[1] * cosTheta;

    return vector(u_x, u_y, 0);
}

scalar taylorCouetteFlow::pressure(const scalar r) const
{
    // TODO(TM)
    return 0; 
}

scalar taylorCouetteFlow::pressureCartesian(const vector& pos) const
{
    // TODO(TM)
    return 0; 
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

