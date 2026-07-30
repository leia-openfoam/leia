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

#include "divGradPsiSnGradAlpha.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
using Foam::mag;

defineTypeNameAndDebug(divGradPsiSnGradAlpha, false);
addToRunTimeSelectionTable(surfaceTensionForce, divGradPsiSnGradAlpha, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
divGradPsiSnGradAlpha::divGradPsiSnGradAlpha(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    // NOTE: a "normals" (nc) field lookup used to live here but was never used by
    // the force -- removed so the model constructs in solvers that do not register
    // geometric normal fields (e.g. the semi-Lagrangian two-phase solver).
    alpha_(mesh_.lookupObject<volScalarField>(surfTensionDict_.getOrDefault<word>("alpha", "alpha.dispersed"))),
    psi_(mesh_.lookupObject<volScalarField>(surfTensionDict_.getOrDefault<word>("levelSet", "psi")))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
divGradPsiSnGradAlpha::calcFaceSurfaceTensionForceFlux() const
{
    if (const surfaceScalarField* sharedKappa =
        registeredFaceCurvature(surfTensionDict_))
    {
        return integratedCSFFlux
        (
            *sharedKappa, alpha_, "GSigmaDivGradPsiConnected"
        );
    }

    // Compute interface-normals using the gradient of the level set field. 
    tmp<volVectorField> nPsiTmp = fvc::grad(psi_);
    nPsiTmp->rename("nPsi");
    volVectorField& nPsi = nPsiTmp.ref();
    // Epsilon-guarded normalisation. Away from the interface grad(psi) can vanish
    // (a local extremum, or exactly at the algebraic-ellipsoid centre), where the
    // unguarded nPsi/mag(nPsi) is 0/0 -> SIGFPE. deltaN (carried in nPsi's own
    // dimensions so this is dimension-safe for any psi units) makes nPsi -> 0 there;
    // those cells carry snGrad(alpha) = 0 and so never enter the surface-tension
    // force -- the normalisation only matters in the interface band, where
    // mag(grad psi) >> deltaN and the guard is inert.
    const dimensionedScalar deltaN("deltaN", nPsi.dimensions(), SMALL);
    nPsi = nPsi / (mag(nPsi) + deltaN);

    // Face-centered curvature as a linear interpolation of the trace of the gradient
    // of the interface-normal-field.
    tmp<surfaceScalarField> tkf = fvc::interpolate(fvc::div(nPsi));
    return integratedCSFFlux(tkf(), alpha_, "GSigmaDivGradPsi");
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
