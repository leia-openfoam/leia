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

\*---------------------------------------------------------------------------*/

#include "redistancer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(redistancer, false);
defineRunTimeSelectionTable(redistancer, Mesh);
addToRunTimeSelectionTable(redistancer, redistancer, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

redistancer::redistancer(const fvMesh& mesh)
    :
        fvSolution_(mesh),
        levelSetDict_(fvSolution_.subDict("levelSet")),
        redistDict_(levelSetDict_.subDict("redistancer")),
        redistanceInterval_(redistDict_.getOrDefault<label>("redistanceInterval", 1))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<Foam::redistancer> redistancer::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution (mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& redistDict = levelSetDict.subDict("redistancer");
    const word& modelType = redistDict.getOrDefault<word>("type", "noRedistancing");
    
    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = MeshConstructorTable(modelType);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr) 
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "redistancer",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<redistancer>(ctorPtr(mesh));
}

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
