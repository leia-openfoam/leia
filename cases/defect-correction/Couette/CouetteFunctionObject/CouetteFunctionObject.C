/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

#include "CouetteFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CouetteFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, CouetteFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CouetteFunctionObject::CouetteFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    cFlow_(dict), 
    pCouette_ 
    (
        IOobject(
            "pCouette",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dimensionedPressure", dimPressure, 0.0)
    ),
    UCouette_
    (
        IOobject
        (
            "UCouette",
            runTime.timeName(), 
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("dimensionedVelocity", dimVelocity, vector::zero)
    ),
    Uerr_("Uerr", UCouette_)
{
    // Compute the Couette pressure and velocity fields 
    
    // Access the mesh cell centers
    const vectorField& C = mesh_.C();

    // Iterate over all cell-centered values 
    forAll(C, cellI)
    {
        // Get the position of the cell center
        const vector& xc = C[cellI];
        vector xcCylindrical (xc[0], xc[1], 0);

        // Compute and store the velocity and pressure for the cell
        UCouette_[cellI] = cFlow_.velocityCartesian(xcCylindrical);
        pCouette_[cellI] = cFlow_.pressureCartesian(xcCylindrical);
    }

    // Iterate over all boundary patches for face-centered values
    forAll(UCouette_.boundaryField(), patchI)
    {
        // Access the current boundary patch
        fvPatchVectorField& velocityPatch = UCouette_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pressurePatch = pCouette_.boundaryFieldRef()[patchI];

        // Access the face centers of the boundary patch
        const vectorField& faceCenters = velocityPatch.patch().Cf();

        // Iterate over all faces on the current patch
        forAll(faceCenters, faceI)
        {
            const vector& xf = faceCenters[faceI];

            vector xfCylindrical (xf[0], xf[1], 0);

            // Compute and set the velocity and pressure for the boundary face
            velocityPatch[faceI] = cFlow_.velocityCartesian(xfCylindrical);
            pressurePatch[faceI] = cFlow_.pressureCartesian(xfCylindrical);
        }
    }

    // Compute the velocity error. 
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    Uerr_ == U - UCouette_; 
    
    // Write fields
    UCouette_.write(); 
    pCouette_.write();
    Uerr_.write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CouetteFunctionObject::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::CouetteFunctionObject::execute()
{
    // Compute the velocity error. 
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    Uerr_ == U - UCouette_; 
    return true;
}


bool Foam::functionObjects::CouetteFunctionObject::end()
{
    volScalarField UerrMag ("UerrMag", mag(Uerr_));

    scalar UerrLinf = gMax(UerrMag); 
    scalar UerrL1 = gSum((UerrMag * mesh_.V())());
    scalar UerrL2 = sqrt(gSum((UerrMag * UerrMag * mesh_.V())()));
    scalar h = gMax(Foam::pow(mesh_.deltaCoeffs(), -1)());

    OFstream of("Couette.csv");
    of << "h,UerrLinf,UerrL1,UerrL2\n"
        << h << "," << UerrLinf << "," 
        << UerrL1 << "," << UerrL2 << "\n";
    
    return true;
}


bool Foam::functionObjects::CouetteFunctionObject::write()
{
    pCouette_.write();
    UCouette_.write();
    Uerr_.write();
    return true;
}


// ************************************************************************* //
