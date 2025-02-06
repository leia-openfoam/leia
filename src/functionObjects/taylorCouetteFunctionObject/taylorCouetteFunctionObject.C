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

#include "taylorCouetteFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(taylorCouetteFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, taylorCouetteFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::taylorCouetteFunctionObject::taylorCouetteFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    tcFlow_(dict), 
    pTaylorCouette_ 
    (
        IOobject(
            "pTaylorCouette",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dimensionedPressure", dimPressure, 0.0)
    ),
    UTaylorCouette_
    (
        "UTaylorCouette",
	mesh_.lookupObject<volVectorField>("U")
    ),
    Uerr_ 
    (
        "Uerr",
	mesh_.lookupObject<volVectorField>("U")
    )
{
    // Compute the Taylor-Couette pressure and velocity fields 
    
    // Access the mesh cell centers
    const vectorField& C = mesh_.C();


    // Iterate over all cell-centered values 
    forAll(C, cellI)
    {
        // Get the position of the cell center
        const vector& xc = C[cellI];
        vector xcCylindrical (xc[0], xc[1], 0);

        // Compute and store the velocity and pressure for the cell
        UTaylorCouette_[cellI] = tcFlow_.velocityCartesian(xcCylindrical);
        pTaylorCouette_[cellI] = tcFlow_.pressureCartesian(xcCylindrical);
    }

    volVectorField& U = const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("U"));

    // Iterate over all boundary patches for face-centered values
    forAll(UTaylorCouette_.boundaryField(), patchI)
    {
        // Access the current boundary patch
        fvPatchVectorField& UTaylorCouettePatch = UTaylorCouette_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pTaylorCouettePatch = pTaylorCouette_.boundaryFieldRef()[patchI];

	if (!isA<fixedValueFvPatchVectorField>(UTaylorCouettePatch))
	{
            // Access the face centers of the boundary patch
            const vectorField& faceCenters = UTaylorCouettePatch.patch().Cf();

            // Iterate over all faces on the current patch
            forAll(faceCenters, faceI)
            {
                const vector& xf = faceCenters[faceI];
		
                // Compute and set the velocity and pressure for the boundary face
                UTaylorCouettePatch[faceI] = tcFlow_.velocityCartesian(xf);
                pTaylorCouettePatch[faceI] = tcFlow_.pressureCartesian(xf);
            }
	}
    }

    // Compute the velocity error. 
    Uerr_ == U - UTaylorCouette_; 
    
    // Write fields
    UTaylorCouette_.write(); 
    pTaylorCouette_.write();
    Uerr_.write();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::taylorCouetteFunctionObject::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::taylorCouetteFunctionObject::execute()
{
    // Compute the velocity error. 
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    Uerr_ == U - UTaylorCouette_; 
    return true;
}


bool Foam::functionObjects::taylorCouetteFunctionObject::end()
{
    volScalarField UerrMag ("UerrMag", mag(Uerr_));

    scalar UerrLinf = gMax(UerrMag); 
    scalar UerrL1 = gSum((UerrMag * mesh_.V())()) / gSum(mesh_.V());
    scalar UerrL2 = sqrt(gSum((UerrMag * UerrMag * mesh_.V())()) / gSum(mesh_.V()));
    scalar h = gMax(Foam::pow(mesh_.deltaCoeffs(), -1)());
    
    const Time& runTime = mesh_.time();
    scalar elapsedCpuTime = runTime.elapsedCpuTime();

    // Local number of cells
    label nCells = mesh_.nCells();

    // Global number of cells (sum over all ranks)
    reduce(nCells, sumOp<label>());

    OFstream of("TaylorCouette.csv");
    of << "h,UerrLinf,UerrL1,UerrL2,N_CELLS,ELAPSED_CPU_TIME\n"
        << h << "," << UerrLinf << "," << UerrL1 << "," << UerrL2 << ","  
        << nCells << "," << elapsedCpuTime << "\n";
    return true;
}


bool Foam::functionObjects::taylorCouetteFunctionObject::write()
{
    pTaylorCouette_.write();
    UTaylorCouette_.write();
    Uerr_.write();
    return true;
}


// ************************************************************************* //
