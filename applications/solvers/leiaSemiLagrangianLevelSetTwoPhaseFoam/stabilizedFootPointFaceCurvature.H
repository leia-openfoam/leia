// Stabilized foot-point FACE curvature delivery -- the winning combination of
// the face-centered curvature gate (2026-08-06, docs/method-comparison):
// arithmetic face interpolation of the raw contour-referenced cell curvature,
// re-referenced to the interface by the parallel-SURFACE inverse at the
// STABILIZED FOOT-POINT signed offset of the face centre
// (slReconstruction::footPointDistance, stable-foot-point-3d.md).
//
// DIMENSION-SEAMLESS (2D and 3D). With kappa = div(n) = k1 + k2 (sum
// convention) and K = k1 k2 the GAUSSIAN curvature, the offset surface at
// signed distance d obeys H_d = (H + K d)/Q, K_d = K/Q, Q = 1 + 2Hd + Kd^2
// (H = kappa/2). Inverting for the interface value:
//
//     kappa_f^Gamma = (kappa_f - 2 K_f d_f) / (1 - kappa_f d_f + K_f d_f^2),
//
//     d_f, K_f = inverse-|psi| (Kang) weighted owner/neighbour values:
//     d from footPointDistance(x_f), K = (g.(cof(H).g))/|g|^4 from
//     fitDerivatives (the implicit-surface Gaussian curvature; cof = adjugate
//     for the symmetric Hessian).
//
// Sphere check: kappa_d = 2/(R+d), K_d = 1/(R+d)^2 -> exactly 2/R. In 2D
// (pseudo-2D fits) K = +0 EXACTLY (every adjugate term carries a zero
// z-component factor), so the formula reduces bit-identically to the measured
// 2D form kappa_f/(1 - d_f kappa_f) -- committed 2D results stay valid. The
// denominator factors as (1 - d k1_d)(1 - d k2_d); the |denom| > 1/2 guard
// trips at d ~ 0.414R on a sphere (vs R in 2D), far outside any active band.
// NOTE the 2D scalar inverse WITHOUT K is first-order-wrong in 3D (it gives
// 2/(R - d) on a sphere): K is not an optional refinement.
//
// Guards: the inverse is skipped past the |1 - kappa_f d_f + K_f d_f^2| > 1/2
// trust bound, on faces whose foot search fails on both sides, and outside
// the CSF force support (|snGrad(alpha)| = 0 faces keep the plain
// interpolate() value the force never uses). K falls back per side to the
// other side's fit, then to 0 (the 2D form; counted) on order-deficient fits.
// Static-gate measurements (exact SDF): 2D circle kappa_f converges at h^2.04
// in the active-face L2 (0.105 1/m at N=512) versus h^1.13 uncorrected.
//
// The result is written into a REGISTERED surface field consumed by
// reconstructedCurvature via
//     surfaceTensionForce { faceCurvatureSource registered;
//                           faceCurvatureField  kappaStableFootFace; }
// PRECONDITIONS (enforced in createSLFields.H): the cell kappa must be the
// raw contour-referenced fill (offsetCorrection none -- anything else
// double-corrects) and the reconstruction must provide footPointDistance/
// fitDerivatives (the uncached quadratic; the base class FatalErrors).
//
// Parallel: footPointDistance/fitDerivatives are cell-local. Internal faces
// use both adjacent fits; coupled (processor) patch faces use the owner-side
// foot AND owner-side K only (the neighbour fit is rank-local) -- one-sided
// but still O(h^2) estimates of the same face quantities.

#ifndef stabilizedFootPointFaceCurvature_H
#define stabilizedFootPointFaceCurvature_H

#include "slReconstruction.H"
#include "syncTools.H"

namespace Foam
{

//- Implicit-surface Gaussian curvature of the cell's fit at its centre:
//  K = (g . (cof(H) . g))/|g|^4. ok = false (K = 0) when the fit is not a
//  full quadratic or the gradient is degenerate. Exactly +0 on pseudo-2D
//  fits (zero z-components make every adjugate contraction term vanish).
inline scalar fitGaussianCurvature
(
    const slReconstruction& recon,
    const label c,
    bool& ok
)
{
    vector g(Zero);
    symmTensor H(Zero);
    ok = false;

    if (recon.fitDerivatives(c, g, H) < 2) { return 0; }

    const scalar gm = mag(g);
    if (gm < SMALL) { return 0; }

    ok = true;
    return (g & (cof(H) & g))/(gm*gm*gm*gm);
}


//- The parallel-surface inverse with trust guard; falls back to the
//  uncorrected value past half the local curvature radius.
inline scalar parallelSurfaceInverse
(
    const scalar kf,
    const scalar d,
    const scalar K
)
{
    const scalar denom = 1.0 - d*kf + K*d*d;
    return (mag(denom) > 0.5) ? (kf - 2.0*K*d)/denom : kf;
}


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

    label nActive = 0, nCorrected = 0, nFootFail = 0, nGaussFallback = 0;

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        ++nActive;

        const scalar wP = mag(psiIn[nei[f]]);   // weight of the OWNER values
        const scalar wN = mag(psiIn[own[f]]);

        bool okO = false, okN = false;
        const scalar dO = recon.footPointDistance(own[f], CfIn[f], 0.0, okO);
        const scalar dN = recon.footPointDistance(nei[f], CfIn[f], 0.0, okN);

        scalar d = 0;
        if (okO && okN)      { d = (dO*wP + dN*wN)/(wP + wN + VSMALL); }
        else if (okO)        { d = dO; }
        else if (okN)        { d = dN; }
        else { ++nFootFail; continue; }

        bool okKO = false, okKN = false;
        const scalar KO = fitGaussianCurvature(recon, own[f], okKO);
        const scalar KN = fitGaussianCurvature(recon, nei[f], okKN);

        scalar K = 0;
        if (okKO && okKN)    { K = (KO*wP + KN*wN)/(wP + wN + VSMALL); }
        else if (okKO)       { K = KO; }
        else if (okKN)       { K = KN; }
        else                 { ++nGaussFallback; }   // K = 0: the 2D form

        const scalar corrected = parallelSurfaceInverse(kfIn[f], d, K);
        if (corrected != kfIn[f]) { ++nCorrected; }
        kfIn[f] = corrected;
    }

    // Coupled (processor) patch faces. Each rank can evaluate only its OWN
    // side's fit; correcting independently leaves the two ranks holding
    // DIFFERENT kappa_f for the same physical face -- an INCONSISTENT
    // capillary flux across processor seams (measured on the filtered
    // stationary droplet: np 8 runs ~7-10x noisier than np 4 at identical
    // (N, t) before this synchronization). Each rank computes its one-sided
    // correction, the values are SWAPPED across the seam, and the two
    // estimates averaged -- both copies then agree exactly (the base
    // interpolate() value is already seam-consistent, and 0.5*(a+b) is
    // commutative). Faces where only one side has a converged foot take that
    // side's value on both ranks.
    surfaceScalarField::Boundary& kfBf = kappaFace.boundaryFieldRef();
    const surfaceScalarField::Boundary& snBf = snA.boundaryField();

    const label nBndFaces = mesh.nBoundaryFaces();
    scalarField bCorr(nBndFaces, 0.0);
    scalarField bOk(nBndFaces, 0.0);
    label nSeamActive = 0;

    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }

        const labelUList& fc = patch.faceCells();
        const vectorField& Cfp = patch.Cf();
        const label bStart = patch.start() - mesh.nInternalFaces();

        forAll(fc, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            ++nSeamActive;   // each physical seam face counted on both ranks

            bool okD = false;
            const scalar dP = recon.footPointDistance(fc[i], Cfp[i], 0.0, okD);
            if (!okD) { continue; }

            bool okK = false;
            const scalar KP = fitGaussianCurvature(recon, fc[i], okK);
            if (!okK) { ++nGaussFallback; }

            bCorr[bStart + i] =
                parallelSurfaceInverse
                (
                    kfBf[patchI][i], dP, okK ? KP : scalar(0)
                );
            bOk[bStart + i] = 1.0;
        }
    }

    scalarField bCorrNei(bCorr);
    scalarField bOkNei(bOk);
    syncTools::swapBoundaryFaceList(mesh, bCorrNei);
    syncTools::swapBoundaryFaceList(mesh, bOkNei);

    label nSeamFootFail = 0;
    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }

        fvsPatchScalarField& kfp = kfBf[patchI];
        const label bStart = patch.start() - mesh.nInternalFaces();

        forAll(kfp, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }

            const label j = bStart + i;
            const bool ownOk = (bOk[j] > 0.5);
            const bool neiOk = (bOkNei[j] > 0.5);

            scalar corrected = kfp[i];
            if (ownOk && neiOk) { corrected = 0.5*(bCorr[j] + bCorrNei[j]); }
            else if (ownOk)     { corrected = bCorr[j]; }
            else if (neiOk)     { corrected = bCorrNei[j]; }
            else                { ++nSeamFootFail; }

            if (corrected != kfp[i]) { ++nCorrected; }
            kfp[i] = corrected;
        }
    }

    nActive += nSeamActive;
    nFootFail += nSeamFootFail;

    reduce(nActive, sumOp<label>());
    reduce(nCorrected, sumOp<label>());
    reduce(nFootFail, sumOp<label>());
    reduce(nGaussFallback, sumOp<label>());
    Info<< "stabilizedFootPointFace curvature: " << nCorrected << "/"
        << nActive << " active faces corrected ("
        << nFootFail << " foot searches failed, "
        << nGaussFallback << " Gaussian-curvature fallbacks; seam faces "
        << "synchronized, counted from both sides)" << endl;
}

} // End namespace Foam

#endif

// ************************************************************************* //
