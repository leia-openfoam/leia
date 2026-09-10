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

#include "slValueBound.H"
#include "slReconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slValueBound, 0);
    defineRunTimeSelectionTable(slValueBound, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slValueBound::slValueBound(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict)
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::slValueBound> Foam::slValueBound::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool legacyClipSwitch
)
{
    word modelType = dict.getOrDefault<word>("valueBound", "fromClipSwitch");

    // The sentinel keeps every case that predates this class bit-identical. The
    // pattern is the one already used for the slope limiter, where the legacy
    // Switch limitSlope drives the default of the newer word slopeLimiter
    // (slReconstruction.C).
    if (modelType == "fromClipSwitch")
    {
        modelType = legacyClipSwitch ? "stencilBounds" : "none";
        Info<< "Selecting slValueBound " << modelType
            << " (from clipToStencilBounds " << Switch(legacyClipSwitch) << ")"
            << endl;
    }
    else
    {
        Info<< "Selecting slValueBound " << modelType << endl;

        // An explicit valueBound OVERRIDES the legacy Switch. Saying so is the
        // whole point: a study that sets both and gets the bound it did not ask
        // for is a confounded arm, and nothing in the log would show it.
        if (legacyClipSwitch && modelType != "stencilBounds")
        {
            WarningInFunction
                << "valueBound " << modelType << " overrides"
                << " clipToStencilBounds true: the quasi-monotone clip is OFF"
                << " in this run." << endl;
        }
    }

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "slValueBound",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<slValueBound>(ctorPtr(mesh, dict));
}

// ************************************************************************* //
