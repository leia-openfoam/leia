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
    Foam::pde2Redistancer

Description
    Level-Set Redistancer using the reinitialization routine of Ferro et al.

    Reference:
    \verbatim
        Ferro, P., Landel, P., Landrodie, C., Guillot, S., & Pescheux, M. (2025). 
        Enhanced Level-Set Method for free surface flow applications. 
        International Journal of Multiphase Flow, 185. 
        https://doi.org/10.1016/j.ijmultiphaseflow.2025.105126
    \endverbatim


SourceFiles
    pde2Redistancer.C

\*---------------------------------------------------------------------------*/

#ifndef pde2Redistancer_H
#define pde2Redistancer_H

#include "redistancer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pde2Redistancer Declaration
\*---------------------------------------------------------------------------*/

class pde2Redistancer
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

        scalar maxCo_;

        bool write_;

public:

    // Static Data Members

    TypeName ("PDE2");

    // Constructors

    //- Construct from components
    pde2Redistancer(const fvMesh& mesh);

    //- Destructor
    virtual ~pde2Redistancer() = default;

    // Member Functions
    virtual void doRedistance(volScalarField& psi); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
