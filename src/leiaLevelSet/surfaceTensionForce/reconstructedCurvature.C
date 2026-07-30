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

#include "reconstructedCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(reconstructedCurvature, false);
addToRunTimeSelectionTable(surfaceTensionForce, reconstructedCurvature, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reconstructedCurvature::reconstructedCurvature(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    kappa_
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
    ),
    psi_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("levelSet", "psi")
        )
    ),
    faceInterpolation_
    (
        surfTensionDict_.getOrDefault<word>("faceInterpolation", "arithmetic")
    ),
    forceWeight_
    (
        surfTensionDict_.getOrDefault<word>("forceWeight", "alpha")
    )
{
    Info<< "reconstructedCurvature: faceInterpolation = " << faceInterpolation_
        << ", forceWeight = " << forceWeight_ << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField>
reconstructedCurvature::calcFaceSurfaceTensionForceFlux() const
{
    const surfaceScalarField* sharedKappa = registeredFaceCurvature
    (
        surfTensionDict_, faceInterpolation_ == "connectedInterface"
    );
    if (sharedKappa)
    {
        // The curvature is already attached to connected zero-surface
        // elements and extended to the CSF faces.  Do not interpolate the
        // diagnostic cell field back to faces.
        if (forceWeight_ == "sharpHeaviside")
        {
            volScalarField Hs
            (
                IOobject("Hs", mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                0.5*(1.0 - sign(psi_)),
                zeroGradientFvPatchScalarField::typeName
            );
            Hs.correctBoundaryConditions();
            return integratedCSFFlux
            (
                *sharedKappa, Hs, "GSigmaReconstructedConnected"
            );
        }
        return integratedCSFFlux
        (
            *sharedKappa, alpha_, "GSigmaReconstructedConnected"
        );
    }

    if (faceInterpolation_ == "interfaceWeighted")
    {
        // Kang/GFM interface-weighted face curvature: the two cell values are
        // the curvatures of the LOCAL psi-contours through the cell centres
        // (parallel curves: kappa_d = kappa/(1 + d kappa) != kappa), so their
        // arithmetic average retains an O(h) normal-offset error that varies
        // across the force support and makes Laplace equilibrium unreachable.
        // Inverse-distance interpolation TO the psi=0 crossing,
        //     kappa_f = (kappa_P |psi_N| + kappa_N |psi_P|)/(|psi_P|+|psi_N|),
        // evaluates the (near-linear in d) parallel-curve law at d = 0 -- the
        // interface curvature -- on every face. Kang, Fedkiw & Liu (2000);
        // Francois et al. (2006, SSF); Abadie et al. (2015): exact balance for
        // level sets, vs significant currents with the arithmetic average.
        tmp<surfaceScalarField> tkf
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "kappaf",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("kappaf", dimless/dimLength, Zero)
            )
        );
        surfaceScalarField& kf = tkf.ref();

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& kIn = kappa_.primitiveField();
        const scalarField& psiIn = psi_.primitiveField();

        scalarField& kfIn = kf.primitiveFieldRef();
        forAll(kfIn, facei)
        {
            const scalar wP = mag(psiIn[nei[facei]]);   // weight of OWNER value
            const scalar wN = mag(psiIn[own[facei]]);   // weight of NEIGHBOUR
            kfIn[facei] =
                (kIn[own[facei]]*wP + kIn[nei[facei]]*wN)/(wP + wN + VSMALL);
        }

        // Boundary faces: one-sided (patch-internal) curvature. Coupled patches
        // use the neighbour-side psi/kappa for the same weighting.
        surfaceScalarField::Boundary& kfbf = kf.boundaryFieldRef();
        forAll(kfbf, patchi)
        {
            const fvPatchScalarField& kp = kappa_.boundaryField()[patchi];
            const fvPatchScalarField& pp = psi_.boundaryField()[patchi];
            if (kp.coupled())
            {
                const scalarField kP(kp.patchInternalField());
                const scalarField kN(kp.patchNeighbourField());
                const scalarField pP(pp.patchInternalField());
                const scalarField pN(pp.patchNeighbourField());
                kfbf[patchi] =
                    (kP*mag(pN) + kN*mag(pP))/(mag(pP) + mag(pN) + VSMALL);
            }
            else
            {
                kfbf[patchi] = kp.patchInternalField();
            }
        }

        if (forceWeight_ == "sharpHeaviside")
        {
            // Faithful LS-SSF/GFM pairing (Abadie et al. 2015): the force weight
            // is the SHARP sign-based Heaviside of psi, so snGrad(H) is nonzero
            // ONLY on sign-change faces -- every force face straddles psi=0 and
            // the Kang formula is a true interpolation on each of them (no
            // same-side faces receiving non-interface curvature). The jump
            // integral sigma*kappa*[H] is unchanged; only the force
            // localisation sharpens. The flow physics (rho, mu, rhoPhi) keeps
            // the geometric alpha.
            volScalarField Hs
            (
                IOobject("Hs", mesh_.time().timeName(), mesh_,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                0.5*(1.0 - sign(psi_)),
                zeroGradientFvPatchScalarField::typeName
            );
            Hs.correctBoundaryConditions();
            return integratedCSFFlux
            (
                tkf(), Hs, "GSigmaReconstructedKangSharp"
            );
        }

        return integratedCSFFlux
        (
            tkf(), alpha_, "GSigmaReconstructedKang"
        );
    }

    // Canonical CSF face force. The curvature was evaluated symbolically from
    // the reconstruction (kappa = div(grad psi/|grad psi|)) and normal-extended
    // by the solver; here it is only interpolated linearly to the faces. Same
    // sign/definition as the trace(grad n) models, so the momentum wiring is
    // unchanged.
    tmp<surfaceScalarField> tkf = fvc::interpolate(kappa_);
    return integratedCSFFlux
    (
        tkf(), alpha_, "GSigmaReconstructedArithmetic"
    );
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
