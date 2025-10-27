/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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
    Foam::pde3Redistancer

Description
    Implementation based on: 
    Mark Sussman, Peter Smereka, Stanley Osher,
    A Level Set Approach for Computing Solutions to Incompressible Two-Phase Flow,
    Journal of Computational Physics,Volume 114, Issue 1, 1994, Pages 146-159,
    https://doi.org/10.1006/jcph.1994.1155.

    Reference:
    \verbatim

    \endverbatim


SourceFiles
    pde3Redistancer.C

\*---------------------------------------------------------------------------*/

#ifndef pde3Redistancer_H
#define pde3Redistancer_H

#include "redistancer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pde3Redistancer Declaration
\*---------------------------------------------------------------------------*/

class pde3Redistancer
:
    public redistancer
{
private:
    // Data Members

    label nCellsForDiffusion_;

    scalar deltaX_;

    dimensionedScalar epsilon_;

    scalar deltaTau_;

    label nIter_;

    Switch write_;

    Switch  constrain_;


public:

    // Static Data Members

    TypeName ("PDE3");

    // Constructors

    //- Construct from components
    pde3Redistancer(const fvMesh& mesh);

    //- Destructor
    virtual ~pde3Redistancer() = default;

    // Member Functions
    virtual void doRedistance(volScalarField& psi); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
