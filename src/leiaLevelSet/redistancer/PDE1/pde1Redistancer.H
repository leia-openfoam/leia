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
    Foam::pde1Redistancer

Description
    Level-Set Redistancer using the reinitialization equation of Sussman, 
    Smereka, and Osher without constraint.

    Reference:
    \verbatim
        Sussman, M., Smereka, P., & Osher, S. (1994). A Level Set Approach for 
        Computing Solutions to Incompressible Two-Phase Flow. Journal of 
        Computational Physics, 114(1), 146–159.
    \endverbatim


SourceFiles
    pde1Redistancer.C

\*---------------------------------------------------------------------------*/

#ifndef pde1Redistancer_H
#define pde1Redistancer_H

#include "redistancer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pde1Redistancer Declaration
\*---------------------------------------------------------------------------*/

class pde1Redistancer
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

        bool write_;
        
public:

    // Static Data Members

    TypeName ("PDE1");

    // Constructors

    //- Construct from components
    pde1Redistancer(const fvMesh& mesh);

    //- Destructor
    virtual ~pde1Redistancer() = default;

    // Member Functions
    virtual void doRedistance(volScalarField& psi); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
