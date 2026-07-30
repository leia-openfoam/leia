/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Tomislav Maric, TU Darmstadt
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

#include "fvPatchFieldsFwd.H"
#include "primitiveFieldsFwd.H"
#include "surfaceTensionForce.H"
#include "traceGradGeoNormalSnGradAlpha.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "volFieldsFwd.H"
#include "levelSetPlaneReconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
using Foam::mag;

defineTypeNameAndDebug(traceGradGeoNormalSnGradAlpha, false);
addToRunTimeSelectionTable(surfaceTensionForce, traceGradGeoNormalSnGradAlpha, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
traceGradGeoNormalSnGradAlpha::traceGradGeoNormalSnGradAlpha(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    alpha_(mesh_.lookupObject<volScalarField>(surfTensionDict_.getOrDefault<word>("alpha", "alpha.dispersed"))),
    psi_(mesh_.lookupObject<volScalarField>(surfTensionDict_.getOrDefault<word>("levelSet", "psi"))),
    narrowBand_(mesh_.lookupObject<volScalarField>(surfTensionDict_.getOrDefault<word>("narrowBand", "NarrowBand")))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
traceGradGeoNormalSnGradAlpha::calcFaceSurfaceTensionForceFlux() const
{
    if (const surfaceScalarField* sharedKappa =
        registeredFaceCurvature(surfTensionDict_))
    {
        return integratedCSFFlux
        (
            *sharedKappa, alpha_, "GSigmaTraceGeoNormalConnected"
        );
    }

    // Compute interface-normals using the gradient of the level set field. 
    tmp<volVectorField> nPsiTmp = fvc::grad(psi_);
    nPsiTmp->rename("nPsi");
    volVectorField& nPsi = nPsiTmp.ref();
    nPsi = nPsi / (mag(nPsi) + dimensionedScalar("deltaN", nPsi.dimensions(), SMALL));
    Info << nPsi.name() << endl;

    // Recompute the same normalized linear-LSQ plane normals used by the
    // Detrixhe-Aslam phase indicator. Unlike the historical implementation,
    // this does not require the separate geometricPhaseIndicator model to
    // register an `nc` field, so the force works with Detrixhe-Aslam directly.
    const coupledFaceNeighbours coupledNei(mesh_, psi_);
    forAll(narrowBand_, cellI)
    {
        if (narrowBand_[cellI] == 1)    
        {
            const scalarList pc =
                leastSquaresPlaneCoeffs(mesh_, psi_, cellI, coupledNei);
            const vector nc(pc[0], pc[1], pc[2]);
            const scalar nmag = mag(nc);
            if (nmag > SMALL)
            {
                nPsi[cellI] = nc/nmag;
            }
        }
    }
    
    // Face-centered curvature as a linear interpolation of the trace of the gradient 
    // of the interface-normal-field. 
    tmp<surfaceScalarField> tkf =
        fvc::interpolate(tr(fvc::grad(nPsi)));
    return integratedCSFFlux(tkf(), alpha_, "GSigmaTraceGeoNormal");
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
