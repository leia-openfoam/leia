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

#include "divGradAlphaSnGradAlpha.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(divGradAlphaSnGradAlpha, false);
addToRunTimeSelectionTable(surfaceTensionForce, divGradAlphaSnGradAlpha, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

divGradAlphaSnGradAlpha::divGradAlphaSnGradAlpha(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    alpha_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("alpha", "alpha.water")
        )
    ),
    nSmooth_(surfTensionDict_.getOrDefault<label>("nSmooth", 2))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
divGradAlphaSnGradAlpha::calcFaceSurfaceTensionForceFlux() const
{
    // Interface normals from the volume-fraction gradient (interFoam style),
    // stabilised away from the interface where grad(alpha) -> 0 with the usual
    // mesh-scale deltaN (grad alpha ~ 0 there anyway, and the force weight
    // snGrad(alpha) vanishes with it).
    const dimensionedScalar deltaN
    (
        "deltaN",
        dimless/dimLength,
        1e-8/cbrt(average(mesh_.V()).value())
    );

    // Smooth alpha for the CURVATURE ONLY (smoothed-CSF): the geometric
    // indicator is a ~one-cell step, on which div(grad alpha/|grad alpha|)
    // picks up O(1/h) support-edge artifacts (measured: instant blow-up).
    // nSmooth cell->face->cell averaging passes widen the transition to the
    // 2-3 cells an algebraic-VoF alpha has; the smoothed field is still a
    // functional of the interface geometry only.
    volScalarField alphaS("alphaS", alpha_);
    for (label k = 0; k < nSmooth_; ++k)
    {
        alphaS = fvc::average(fvc::interpolate(alphaS));
    }

    if (const surfaceScalarField* sharedKappa =
        registeredFaceCurvature(surfTensionDict_))
    {
        // Retain this model's alpha-based force localisation, but consume the
        // same connected face curvature as the level-set CSF models.
        return integratedCSFFlux
        (
            *sharedKappa, alphaS, "GSigmaAlphaCSFConnected"
        );
    }

    tmp<volVectorField> gradAlphaTmp = fvc::grad(alphaS);
    const volVectorField& gradAlpha = gradAlphaTmp();

    volVectorField nAlpha("nAlpha", gradAlpha/(mag(gradAlpha) + deltaN));

    // Curvature from alpha alone: kappa = -div(nAlpha). Sign: alpha decreases
    // outward, so nAlpha points INWARD; -div(nAlpha) matches the outward-normal
    // convention of div(grad psi/|grad psi|) (positive for a drop, = 1/R).
    //
    // BALANCED-FORCE CONSISTENCY: the force weight must be the snGrad of the
    // SAME field the curvature was computed from (interFoam pairs kappa(alpha)
    // with snGrad(alpha) of one and the same alpha). Pairing kappa(alphaS) with
    // the SHARP snGrad(alpha_) mixes two discrete interfaces and is not
    // balanced even for constant kappa (measured: fast blow-up). With the
    // smoothed pair the discrete jump spreads over the smoothed transition,
    // exactly as in algebraic-VoF CSF; its integral (the Laplace jump) is
    // unchanged. The flow physics (rho, rhoPhi) keeps the sharp alpha.
    tmp<surfaceScalarField> tkf =
        fvc::interpolate(-fvc::div(nAlpha));
    return integratedCSFFlux(tkf(), alphaS, "GSigmaAlphaCSF");
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
