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

#include "slScheme.H"
#include "fvSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slScheme, 0);
    defineRunTimeSelectionTable(slScheme, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slScheme::slScheme(const fvMesh& mesh)
:
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::slScheme> Foam::slScheme::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution(mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& slDict = levelSetDict.subOrEmptyDict("semiLagrangian");
    const word modelType = slDict.getOrDefault<word>("scheme", "pointValue");
    Info<< "Selecting slScheme " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "slScheme",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<slScheme>(ctorPtr(mesh));
}

// ************************************************************************* //
