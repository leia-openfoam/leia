/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "isoCurvatureSurfaceTension.H"
#include "isoGeometricCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
using Foam::mag;

defineTypeNameAndDebug(isoCurvatureSurfaceTension, false);
addToRunTimeSelectionTable(surfaceTensionForce, isoCurvatureSurfaceTension, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

isoCurvatureSurfaceTension::isoCurvatureSurfaceTension(const fvMesh& mesh)
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
    field_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("field", "psi")
        )
    ),
    isoValue_(surfTensionDict_.getOrDefault<scalar>("isoValue", 0.0)),
    kappaName_(surfTensionDict_.getOrDefault<word>("kappa", "kappa")),
    sgn_(surfTensionDict_.getOrDefault<scalar>("sign", 1.0)),
    nExtend_(surfTensionDict_.getOrDefault<label>("nExtend", 2)),
    nSmooth_(surfTensionDict_.getOrDefault<label>("nSmooth", 0)),
    relax_(surfTensionDict_.getOrDefault<scalar>("relax", 1.0)),
    method_
    (
        surfTensionDict_.getOrDefault<word>("method", "conormal") == "parabola"    ? 1 :
        surfTensionDict_.getOrDefault<word>("method", "conormal") == "conormalCSF" ? 2 : 0
    )
{
    Info<< "isoCurvatureSurfaceTension: interface = {" << field_.name()
        << " = " << isoValue_ << "}, CSF weight snGrad(" << alpha_.name()
        << "), fills field '" << kappaName_ << "'" << nl;
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

tmp<surfaceScalarField>
isoCurvatureSurfaceTension::calcFaceSurfaceTensionForceFlux() const
{
    if (const surfaceScalarField* sharedKappa =
        registeredFaceCurvature(surfTensionDict_))
    {
        return integratedCSFFlux
        (
            *sharedKappa, alpha_, "GSigmaIsoCSFConnected"
        );
    }

    if (method_ == 0)
    {
        // DIRECT CONORMAL FORCE (integral surface tension). The per-cell surface-tension
        // body force per volume is sigma * g, g = -sgn*(sum_f c_f)/V_c -- divided by the
        // CELL VOLUME (never zero), so it is well conditioned for every cut geometry
        // (unlike the pointwise curvature, which needs /A_c -> 0 at corner cuts). The
        // solver contract wants the integrated FACE-NORMAL component of this body
        // force, so return interpolate(sigma*g) . Sf [N/m].
        volVectorField g
        (
            IOobject("gConormal", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_, dimensionedVector("g", dimless/dimArea, Zero), "zeroGradient"
        );
        const label nGood = computeIsoConormalForceDensity
        (
            mesh_, field_, isoValue_, sgn_, nExtend_, nSmooth_, g
        );
        const volVectorField fStVol(sigma_*g);       // body force [N/m^3]
        tmp<surfaceScalarField> tf
        (
            fvc::interpolate(fStVol) & mesh_.Sf()
        );

        Info<< "isoCurvatureSurfaceTension[conormal-force]: |fSigma_vol| max = "
            << gMax(mag(fStVol.primitiveField())) << " N/m^3, Gsigma in ["
            << gMin(tf().primitiveField()) << ", " << gMax(tf().primitiveField())
            << "] N/m; interface cells = " << returnReduce(nGood, sumOp<label>())
            << nl;
        return tf;
    }

    if (method_ == 2)
    {
        // CONORMAL-CSF (integral formulation INSIDE the CSF): the conormal integral
        // gives the cell force density g = -sgn*(sum_f n_f x L_f).n_c n_c / V_c, which
        // in the continuum is kappa*grad(alpha). Instead of dividing by the GEOMETRIC
        // interface area A_c (-> 0 at corner cuts, ill-posed) and multiplying by the
        // interface delta again in the CSF, normalise by the ALPHA-GRADIENT interface
        // measure (coarea: A_c ~ |grad alpha|_c V_c, never degenerate in the band):
        //     kappa_c = (g . grad alpha)_c / |grad alpha|_c^2      [1/m]
        // and apply it through the BALANCED CSF face structure sigma*interp(kappa)*
        // snGrad(alpha) -- the area weighting is carried by snGrad(alpha), which shares
        // the discrete operator with snGrad(p_rgh). Extension of kappa (unlike the
        // direct force) does NOT over-count: snGrad(alpha) self-normalises the integral.
        volVectorField g
        (
            IOobject("gConormal", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_, dimensionedVector("g", dimless/dimArea, Zero), "zeroGradient"
        );
        boolList filled;
        const label nGood = computeIsoConormalForceDensity
        (
            mesh_, field_, isoValue_, sgn_, /*nExtend*/ 0, /*nSmooth*/ 0, g, &filled
        );

        const volVectorField gradAlpha(fvc::grad(alpha_));
        volScalarField& kappa = mesh_.lookupObjectRef<volScalarField>(kappaName_);
        scalarField& k = kappa.primitiveFieldRef();
        k = 0.0;
        const vectorField& gI  = g.primitiveField();
        const vectorField& gaI = gradAlpha.primitiveField();
        forAll(k, c)
        {
            if (!filled[c]) { continue; }
            const scalar ga2 = magSqr(gaI[c]);
            if (ga2 < SMALL) { filled[c] = false; continue; }   // no alpha support
            k[c] = (gI[c] & gaI[c])/ga2;
        }
        extendSmoothCellScalar(mesh_, filled, nExtend_, nSmooth_, k);
        kappa.correctBoundaryConditions();

        tmp<surfaceScalarField> tkf = fvc::interpolate(kappa);
        tmp<surfaceScalarField> tf = integratedCSFFlux
        (
            tkf(), alpha_, "GSigmaIsoConormalCSF"
        );

        scalar kMin = GREAT, kMax = -GREAT;
        forAll(k, c) { if (filled[c]) { kMin = min(kMin, k[c]); kMax = max(kMax, k[c]); } }
        Info<< "isoCurvatureSurfaceTension[conormalCSF]: kappa(interface) in ["
            << returnReduce(kMin, minOp<scalar>()) << ", "
            << returnReduce(kMax, maxOp<scalar>()) << "] 1/m, fSigma in ["
            << gMin(tf().primitiveField()) << ", " << gMax(tf().primitiveField())
            << "] N/m; interface cells = " << returnReduce(nGood, sumOp<label>())
            << nl;
        return tf;
    }

    // PARABOLA method: local height-function curvature -> standard CSF face force.
    volScalarField& kappa = mesh_.lookupObjectRef<volScalarField>(kappaName_);
    const scalarField kOld(kappa.primitiveField());   // for time under-relaxation
    boolList filled;
    label nHole = 0;
    const label nGood = computeIsoGeometricCurvature
    (
        mesh_, field_, isoValue_, sgn_, nExtend_, nSmooth_, kappa, filled, nHole,
        method_
    );
    if (relax_ < 1.0 - SMALL)
    {
        scalarField& k = kappa.primitiveFieldRef();
        forAll(k, c) { if (filled[c]) { k[c] = relax_*k[c] + (1.0 - relax_)*kOld[c]; } }
    }
    kappa.correctBoundaryConditions();

    tmp<surfaceScalarField> tkf = fvc::interpolate(kappa);
    tmp<surfaceScalarField> tf = integratedCSFFlux
    (
        tkf(), alpha_, "GSigmaIsoParabolaCSF"
    );

    const scalarField& k = kappa.primitiveField();
    scalar kMin = GREAT, kMax = -GREAT;
    forAll(k, c) { if (filled[c]) { kMin = min(kMin, k[c]); kMax = max(kMax, k[c]); } }
    Info<< "isoCurvatureSurfaceTension[parabola]: kappa(interface) in ["
        << returnReduce(kMin, minOp<scalar>()) << ", "
        << returnReduce(kMax, maxOp<scalar>()) << "] 1/m, fSigma in ["
        << gMin(tf().primitiveField()) << ", " << gMax(tf().primitiveField())
        << "] N/m; interface cells: " << returnReduce(nGood, sumOp<label>())
        << " clean, " << returnReduce(nHole, sumOp<label>()) << " fit-failed" << nl;
    return tf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
