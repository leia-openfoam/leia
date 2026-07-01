/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "velocityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityModel, 0);
    defineRunTimeSelectionTable(velocityModel, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityModel::velocityModel(const fvMesh& mesh)
:
    fvSolution_(mesh),
    velocityDict_(fvSolution_.subDict("velocityModel")),
    isOscillating_(velocityDict_.getOrDefault<Switch>("oscillation", "on")),
    tau_(velocityDict_.getOrDefault<scalar>("tau", mesh.time().endTime().value()))
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::velocityModel> Foam::velocityModel::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution(mesh);
    const dictionary& velocityDict = fvSolution.subDict("velocityModel");
    const word modelType = velocityDict.getOrDefault<word>("type", "none");
    Info<< "Selecting advection velocity " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "velocityModel",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<velocityModel>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::velocityModel::oscillationFactor
(
    const scalar time,
    const scalar tau
) const
{
    return Foam::cos(M_PI * time / tau);
}


void Foam::velocityModel::oscillateVelocity
(
    volVectorField& U,
    const volVectorField& U0,
    surfaceScalarField& phi,
    const surfaceScalarField& phi0,
    const Time& runTime
) const
{
    if (isOscillating())
    {
        const scalar osFactor =
            oscillationFactor(runTime.timeOutputValue(), tau_);

        phi == phi0 * osFactor;
        U == U0 * osFactor;
    }
}


Foam::vector Foam::velocityModel::velocity(const vector& /*x*/) const
{
    notImplemented("Foam::velocityModel::velocity");
    return vector(0, 0, 0);
}


void Foam::velocityModel::setVelocity(volVectorField& U) const
{
    const fvMesh& mesh = U.mesh();
    const auto& cellCenters = mesh.C();

    forAll(U, cellI)
    {
        U[cellI] = velocity(cellCenters[cellI]);
    }

    const auto& Cf = mesh.Cf();
    auto& UboundaryField = U.boundaryFieldRef();
    const auto& CfBoundaryField = Cf.boundaryField();
    const auto& meshBoundary = mesh.boundary();

    forAll(meshBoundary, patchI)
    {
        const auto& CfPatchField = CfBoundaryField[patchI];
        auto& UpatchField = UboundaryField[patchI];
        forAll(UpatchField, faceI)
        {
            UpatchField[faceI] = velocity(CfPatchField[faceI]);
        }
    }
}


void Foam::velocityModel::setVolumetricFlux(surfaceScalarField& phi) const
{
    const fvMesh& mesh = phi.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField& Cf = mesh.Cf();

    forAll(Cf, faceID)
    {
        phi[faceID] = (velocity(Cf[faceID]) & Sf[faceID]);
    }

    const auto& CfBoundaryField = Cf.boundaryField();
    const auto& SfBoundaryField = Sf.boundaryField();
    auto& phiBoundaryField = phi.boundaryFieldRef();
    const auto& meshBoundary = mesh.boundary();

    forAll(meshBoundary, patchI)
    {
        const auto& CfPatchField = CfBoundaryField[patchI];
        const auto& SfPatchField = SfBoundaryField[patchI];
        auto& phiPatchField = phiBoundaryField[patchI];
        forAll(phiPatchField, faceI)
        {
            phiPatchField[faceI] =
            (
                velocity(CfPatchField[faceI]) & SfPatchField[faceI]
            );
        }
    }
}

// ************************************************************************* //
