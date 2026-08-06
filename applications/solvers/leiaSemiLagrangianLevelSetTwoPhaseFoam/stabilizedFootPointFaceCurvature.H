// Stabilized foot-point FACE curvature delivery -- the winning combination of
// the face-centered curvature gate (2026-08-06, docs/method-comparison):
// arithmetic face interpolation of the raw contour-referenced cell curvature,
// re-referenced to the interface by the parallel-curve inverse at the
// STABILIZED FOOT-POINT signed offset of the face centre
// (slReconstruction::footPointDistance, stable-foot-point-3d.md):
//
//     kappa_f^Gamma = kappa_f / (1 - d_f kappa_f),
//     d_f = (d_own |psi_nei| + d_nei |psi_own|) / (|psi_own| + |psi_nei|)
//
// (inverse-|psi| Kang weighting: the near-interface side's fit dominates).
// Guards: the inverse is skipped past half the local curvature radius
// (|1 - d kappa| <= 1/2), on faces whose foot search fails on both sides, and
// outside the CSF force support (|snGrad(alpha)| = 0 faces keep the plain
// interpolate() value the force never uses). Static-circle measurement:
// kappa_f converges at h^2.0 in the active-face L2 (0.105 1/m at N=512)
// versus h^1.1 for the uncorrected delivery (11.35 1/m) -- the offset, not
// the estimator, is the leading defect, and the correction MUST act on the
// face-combined value (per-side correction reaches only h^1.4).
//
// The result is written into a REGISTERED surface field consumed by the
// balanced-CSF reconstructedCurvature force model via
//     surfaceTensionForce { faceCurvatureSource registered;
//                           faceCurvatureField  kappaStableFootFace; }
// PRECONDITIONS (enforced in createSLFields.H): the cell kappa must be the
// raw contour-referenced fill (offsetCorrection none -- anything else
// double-corrects) and the reconstruction must provide footPointDistance
// (quadratic reconstructions; the base class FatalErrors otherwise).
//
// Parallel: footPointDistance is cell-local polynomial arithmetic. Internal
// faces use both adjacent fits; coupled (processor) patch faces use the
// owner-side foot only (the neighbour fit is rank-local) -- a one-sided but
// still O(h^2) offset estimate. Non-coupled physical patches keep the
// interpolate() value (snGrad(alpha) = 0 there for the droplet benchmarks).

#ifndef stabilizedFootPointFaceCurvature_H
#define stabilizedFootPointFaceCurvature_H

#include "slReconstruction.H"

namespace Foam
{

inline void computeStabilizedFootPointFaceCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& alpha,
    const volScalarField& kappa,
    const slReconstruction& recon,
    surfaceScalarField& kappaFace
)
{
    // The CSF arithmetic delivery everywhere; corrected on the force support.
    kappaFace = fvc::interpolate(kappa);

    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snIn = snA.primitiveField();
    const scalar snMax = max(gMax(mag(snIn)()), VSMALL);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& CfIn = mesh.Cf().primitiveField();
    const scalarField& psiIn = psi.primitiveField();

    scalarField& kfIn = kappaFace.primitiveFieldRef();

    label nActive = 0, nCorrected = 0, nFootFail = 0;

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        ++nActive;

        bool okO = false, okN = false;
        const scalar dO = recon.footPointDistance(own[f], CfIn[f], 0.0, okO);
        const scalar dN = recon.footPointDistance(nei[f], CfIn[f], 0.0, okN);

        scalar d = 0;
        if (okO && okN)
        {
            const scalar wP = mag(psiIn[nei[f]]);   // weight of the OWNER foot
            const scalar wN = mag(psiIn[own[f]]);
            d = (dO*wP + dN*wN)/(wP + wN + VSMALL);
        }
        else if (okO) { d = dO; }
        else if (okN) { d = dN; }
        else { ++nFootFail; continue; }

        const scalar denom = 1.0 - d*kfIn[f];
        if (mag(denom) > 0.5)
        {
            kfIn[f] /= denom;
            ++nCorrected;
        }
    }

    // Coupled (processor) patch faces: owner-side foot only.
    surfaceScalarField::Boundary& kfBf = kappaFace.boundaryFieldRef();
    const surfaceScalarField::Boundary& snBf = snA.boundaryField();

    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }

        fvsPatchScalarField& kfp = kfBf[patchI];
        const labelUList& fc = patch.faceCells();
        const vectorField& Cfp = patch.Cf();

        forAll(kfp, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            ++nActive;

            bool okO = false;
            const scalar dP = recon.footPointDistance(fc[i], Cfp[i], 0.0, okO);
            if (!okO) { ++nFootFail; continue; }

            const scalar denom = 1.0 - dP*kfp[i];
            if (mag(denom) > 0.5)
            {
                kfp[i] /= denom;
                ++nCorrected;
            }
        }
    }

    reduce(nActive, sumOp<label>());
    reduce(nCorrected, sumOp<label>());
    reduce(nFootFail, sumOp<label>());
    Info<< "stabilizedFootPointFace curvature: " << nCorrected << "/"
        << nActive << " active faces corrected ("
        << nFootFail << " foot searches failed)" << endl;
}

} // End namespace Foam

#endif

// ************************************************************************* //
