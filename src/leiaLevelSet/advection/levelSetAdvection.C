/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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

#include "levelSetAdvection.H"
#include "fvMesh.H"
#include "fvSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(levelSetAdvection, false);
defineRunTimeSelectionTable(levelSetAdvection, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

levelSetAdvection::levelSetAdvection
(
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const velocityModel& velModel
)
:
    mesh_(mesh),
    U_(U),
    phi_(phi),
    velModel_(velModel)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<levelSetAdvection> levelSetAdvection::New
(
    const fvMesh& mesh,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const velocityModel& velModel
)
{
    const fvSolution& fvSolutionDict(mesh);
    const dictionary& levelSetDict = fvSolutionDict.subDict("levelSet");
    const dictionary advDict = levelSetDict.subOrEmptyDict("advection");
    const word modelType = advDict.getOrDefault<word>("type", "eulerian");

    Info<< "Selecting levelSetAdvection type: " << modelType << nl << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolutionDict,
            "levelSetAdvection",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<levelSetAdvection>(ctorPtr(mesh, U, phi, velModel));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
