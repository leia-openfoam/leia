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

#include "curvaturePressurePotential.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(curvaturePressurePotential, false);
addToRunTimeSelectionTable
(
    surfaceTensionForce,
    curvaturePressurePotential,
    Mesh
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvaturePressurePotential::curvaturePressurePotential(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    curvature_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("kappa", "kappa")
        )
    ),
    alpha_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("alpha", "alpha.water")
        )
    )
{
    Info<< "curvaturePressurePotential: DIAGNOSTIC arm -- the capillary force "
        << "is the exact discrete gradient snGrad(sigma*kappa*alpha)*|Sf|, so "
        << "it produces NO velocity and CANNOT relax a deformed interface. "
        << "The discarded term -sigma*alpha*grad(kappa) is the shape-relaxation "
        << "driver." << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
curvaturePressurePotential::calcFaceSurfaceTensionForceFlux() const
{
    // sigma*kappa*alpha as a CELL field, then ONE snGrad: by construction in
    // the range of the same operator the pressure equation inverts.
    tmp<volScalarField> tCapillaryPressure(sigma_*curvature_*alpha_);

    tmp<surfaceScalarField> tFlux
    (
        fvc::snGrad(tCapillaryPressure())*mesh_.magSf()
    );
    tFlux.ref().rename("GSigmaCurvaturePressurePotential");
    return tFlux;
}

} // End namespace Foam

// ************************************************************************* //
