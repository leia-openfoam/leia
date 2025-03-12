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

    A concrete class implementing the Poiseuille flow, i.e. the flow between two plates.
SourceFiles
    PoiseuilleFlow.C

\*---------------------------------------------------------------------------*/

#ifndef PoiseuilleFlow_H
#define PoiseuilleFlow_H
// prevent double declaration

#include "scalar.H"
#include "vector.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class PoiseuilleFlow Declaration
\*---------------------------------------------------------------------------*/

class PoiseuilleFlow
{
    private:
	const scalar H_;
	const scalar P_;
	const scalar Vis_;
	const scalar L_;

    public:

	PoiseuilleFlow
	(
	    const scalar H,
	    const scalar Vis,
	    const scalar L,
	    const scalar& P
	);

	PoiseuilleFlow(const dictionary& dict);

	vector velocityCartesian(const vector& pos) const;

	scalar pressureCartesian(const vector& pos) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
