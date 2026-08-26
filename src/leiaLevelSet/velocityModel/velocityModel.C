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
    // A bool literal, not "on": Switch(const char*) is explicit, so a string
    // literal decays to a non-null pointer and converts to bool as TRUE. Here
    // that happened to coincide with the intended default, so this site was
    // benign -- but it is the same construct that made the `off` default in
    // uniaxialStrain mean ON, so it is written unambiguously.
    isOscillating_(velocityDict_.getOrDefault<Switch>("oscillation", true)),
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
        // COUPLED (processor, cyclic) PATCHES MUST BE SKIPPED.
        //
        // A coupled patch field does NOT hold the value at the face: it holds
        // the value in the NEIGHBOURING CELL CENTRE. That is what
        // processorFvPatchField::patchNeighbourField() returns, and it is what
        // every coupled-patch operator interpolates against --- Gauss gradient
        // forms w*U_own + (1 - w)*U_neighbour and expects the second term to be
        // a cell value.
        //
        // Writing velocity(Cf) here therefore stores a FACE value where a CELL
        // value is expected. On a uniform mesh with w = 1/2 the interpolation
        // returns 0.5*U(x_P) + 0.5*U(x_f) instead of U(x_f), an error of
        // 0.5*(U(x_f) - U(x_P)) = a*h/4 for a linear field of strain a. Divided
        // by the cell volume in the Gauss sum that is an O(a) --- NOT O(h) ---
        // error in grad(U) in every cell adjacent to a processor boundary, and
        // it does not vanish under refinement.
        //
        // MEASURED on cases/1Dstretch (v = a x e_x, a = 1 1/s, exact solution
        // known), band mean d(psi)/dx at t = 1, N = 64:
        //            serial      np=2       np=4       np=8
        //   noSource 0.367872  0.367872   0.367872   0.367872   (exact e^-1)
        //   R        1.000000  1.004681   1.005367   0.956834
        //   Rdiv     0.367872  0.322012   0.317157   0.312497
        // The sourceless arm is clean to 1.7e-07 at every processor count
        // because setVolumetricFlux writes phi PER FACE, where the face-centre
        // value IS the correct flux, so advection never sees this. Only
        // fvc::grad(U) does --- and the SDPLS source is its only consumer.
        // Rdiv is hit twice: w = (nhat & U) nhat inherits the error and
        // fvc::flux(w) then interpolates it, so phiW != phi at coupled faces
        // and the discrete cancellation the form depends on (exact to the last
        // digit in serial) fails.
        //
        // Leaving them to correctBoundaryConditions() below is also the more
        // accurate choice: the halo exchange returns the TRUE neighbour cell
        // value, not an analytic approximation to it, so a coupled face
        // interpolates exactly as an internal face does.
        if (meshBoundary[patchI].coupled())
        {
            continue;
        }

        const auto& CfPatchField = CfBoundaryField[patchI];
        auto& UpatchField = UboundaryField[patchI];
        forAll(UpatchField, faceI)
        {
            UpatchField[faceI] = velocity(CfPatchField[faceI]);
        }
    }

    // Halo exchange: fills every coupled patch with the neighbour cell values.
    U.correctBoundaryConditions();
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
