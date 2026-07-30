/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

#include "constantCurvaturePressurePotential.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSnGrad.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(constantCurvaturePressurePotential, false);
addToRunTimeSelectionTable
(
    surfaceTensionForce,
    constantCurvaturePressurePotential,
    Mesh
);

constantCurvaturePressurePotential::constantCurvaturePressurePotential
(
    const fvMesh& mesh
)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    curvature_
    (
        "curvature",
        pow(dimLength, -1),
        surfTensionDict_.get<scalar>("curvature")
    ),
    alpha_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("alpha", "alpha.dispersed")
        )
    )
{}

tmp<surfaceScalarField>
constantCurvaturePressurePotential::calcFaceSurfaceTensionForceFlux() const
{
    tmp<volScalarField> tCapillaryPressure
    (
        sigma_*curvature_*alpha_
    );

    tmp<surfaceScalarField> tFlux
    (
        fvc::snGrad(tCapillaryPressure())*mesh_.magSf()
    );
    tFlux.ref().rename("GSigmaConstantCurvaturePressurePotential");
    return tFlux;
}

} // End namespace Foam

// ************************************************************************* //
