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

#ifndef taylorCouetteFlow_H
#define taylorCouetteFlow_H

#include "scalar.H"
#include "vector.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class taylorCouetteFlow Declaration
\*---------------------------------------------------------------------------*/

class taylorCouetteFlow
{
    private:
	const scalar Ri_;
	const scalar sqrRi_;
	const scalar Ro_;
	const scalar sqrRo_;
	const scalar omegai_;
	const scalar omegao_;
	const vector axisPoint_;

    public:

	taylorCouetteFlow
	(
	    const scalar Ri, 
	    const scalar Ro, 
	    const scalar& omegai, 
	    const scalar& omegao, 
	    const vector& axisPoint
	);

	taylorCouetteFlow(const dictionary& dict);


	vector velocityCylindrical(const scalar r) const;

	vector velocityCartesian(const vector& pos) const;

	scalar pressure(const scalar r) const;

	scalar pressureCartesian(const vector& pos) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
