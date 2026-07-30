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

#include "surfaceTensionForce.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(surfaceTensionForce, false);
defineRunTimeSelectionTable(surfaceTensionForce, Mesh);

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceTensionForce>
Foam::surfaceTensionForce::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution (mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& dict = levelSetDict.subDict("surfaceTensionForce");
    const word& modelType = dict.get<word>("type");
    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "surfaceTensionForce",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<surfaceTensionForce>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceTensionForce::surfaceTensionForce(const fvMesh& mesh)
:
    mesh_(mesh),
    runTime_(mesh.time()), 
    transportProperties_
    (
        IOobject
        (
            "transportProperties", 
            "constant", 
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )    
    ), 
    sigma_
    (
        "sigma",
        dimForce / dimLength, 
        transportProperties_.get<scalar>("sigma")
    )
{}


const Foam::surfaceScalarField*
Foam::surfaceTensionForce::registeredFaceCurvature
(
    const dictionary& modelDict,
    const bool legacyConnectedSelection
) const
{
    word source
    (
        modelDict.getOrDefault<word>
        (
            "faceCurvatureSource",
            legacyConnectedSelection ? "registered" : "model"
        )
    );

    if (source == "model")
    {
        return nullptr;
    }
    if (source == "connectedInterface")
    {
        source = "registered";
    }
    if (source != "registered")
    {
        FatalIOErrorInFunction(modelDict)
            << "Unknown faceCurvatureSource '" << source << "'. Valid: "
            << "model, registered (or connectedInterface as a synonym)."
            << exit(FatalIOError);
    }

    const word fieldName
    (
        modelDict.getOrDefault<word>
        (
            "faceCurvatureField",
            "kappaInterfaceFace"
        )
    );
    if (!mesh_.foundObject<surfaceScalarField>(fieldName))
    {
        FatalIOErrorInFunction(modelDict)
            << "faceCurvatureSource=" << source << " requires registered "
            << "surfaceScalarField '" << fieldName << "', but it is absent."
            << exit(FatalIOError);
    }

    return &mesh_.lookupObject<surfaceScalarField>(fieldName);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::surfaceTensionForce::integratedCSFFlux
(
    const surfaceScalarField& kappaFace,
    const volScalarField& forceWeight,
    const word& resultName
) const
{
    tmp<surfaceScalarField> tFlux
    (
        sigma_*kappaFace*fvc::snGrad(forceWeight)*mesh_.magSf()
    );
    tFlux.ref().rename(resultName);
    return tFlux;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::surfaceTensionForce::faceSurfaceTensionForceFlux() const
{
    tmp<surfaceScalarField> tFlux = calcFaceSurfaceTensionForceFlux();
    const dimensionSet expected(dimForce/dimLength);

    if (tFlux().dimensions() != expected)
    {
        FatalErrorInFunction
            << "Surface-tension model '" << type() << "' returned dimensions "
            << tFlux().dimensions() << " for " << tFlux().name() << nl
            << "The production contract requires the integrated scalar face "
            << "force flux G_sigma,f with dimensions " << expected << "."
            << abort(FatalError);
    }
    if (!tFlux().is_oriented())
    {
        FatalErrorInFunction
            << "Surface-tension model '" << type() << "' returned an "
            << "unoriented field " << tFlux().name() << "." << nl
            << "G_sigma,f must be owner-oriented so it can be combined with "
            << "snGrad(p_rgh)|Sf| in the pressure-flux space."
            << abort(FatalError);
    }

    return tFlux;
}

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
