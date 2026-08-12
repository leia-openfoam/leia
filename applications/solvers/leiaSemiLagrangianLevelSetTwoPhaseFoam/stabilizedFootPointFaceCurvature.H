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



//- CUT-CELL delivery (dimension-agnostic): invert ONCE PER CUT CELL and assign
//  that single interface value to every active face of the cell.
//
//  Per cut cell (0 < alpha < 1 -- the geometric indicator's own definition,
//  identical in 2D and 3D):
//      kappa_c, K_c  from the cell's own quadratic fit AT ITS CENTRE,
//      d_c           = stabilized foot-point distance of the cell centre,
//      kappa^Gamma_c = parallelSurfaceInverse(kappa_c, d_c, K_c).
//  Per active face: the value of its cut-cell side; the mean of the two if both
//  sides are cut cells; the face-based construction as a fallback if neither is
//  (counted -- with a sharp geometric indicator every active face has at least
//  one cut side, so this should stay at zero).
//
//  The point is structural, not accuracy: all active faces of one cut cell
//  carry the SAME number, so the difference of kappa_f taken across the force
//  support vanishes identically at every step, however far psi has drifted from
//  a signed distance. WP8.1 measured that across-support difference to be the
//  dominant part of the driver (4-6x the along-interface part) and to regrow to
//  parity within ~0.05 s when only the accuracy is improved.
//
//  Evaluating the fit's curvature ON the interface instead would NOT do this:
//  the fitted kappa is nearly constant along the normal (constant Hessian), so
//  it keeps the O(h) offset -- measured 11.3492 vs 11.3502 (2D) and 0.4184 vs
//  0.4191 (3D). The interface value has to be produced algebraically.
inline void computeCutCellFootPointFaceCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& alpha,
    const volScalarField& kappa,
    const slReconstruction& recon,
    surfaceScalarField& kappaFace
)
{
    // start from the classical delivery: only the ACTIVE faces are overwritten
    kappaFace = fvc::interpolate(kappa);

    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snIn = snA.primitiveField();
    const scalar snMax = max(gMax(mag(snIn)()), VSMALL);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& CIn = mesh.C().primitiveField();
    const scalarField& aIn = alpha.primitiveField();

    // --- one interface curvature per cut cell -----------------------------
    scalarField kappaGamma(mesh.nCells(), 0.0);
    scalarField isCut(mesh.nCells(), 0.0);
    label nCut = 0, nFootFail = 0, nGaussFallback = 0;

    forAll(aIn, c)
    {
        if (aIn[c] <= SMALL || aIn[c] >= 1 - SMALL) { continue; }

        vector g(Zero);
        symmTensor H(Zero);
        if (recon.fitDerivatives(c, g, H) < 2) { continue; }
        const scalar gm = mag(g);
        if (gm < SMALL) { continue; }

        const scalar kc = (tr(H)*gm*gm - (g & (H & g)))/(gm*gm*gm);
        bool okK = false;
        const scalar Kc = fitGaussianCurvature(recon, c, okK);
        if (!okK) { ++nGaussFallback; }

        bool okD = false;
        const scalar dc = recon.footPointDistance(c, CIn[c], 0.0, okD);
        if (!okD) { ++nFootFail; continue; }

        kappaGamma[c] = parallelSurfaceInverse(kc, dc, okK ? Kc : scalar(0));
        isCut[c] = 1.0;
        ++nCut;
    }

    // --- assign to the active faces ---------------------------------------
    scalarField& kfIn = kappaFace.primitiveFieldRef();
    label nActive = 0, nAssigned = 0, nNoCutSide = 0;

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        ++nActive;

        const bool cO = (isCut[own[f]] > 0.5);
        const bool cN = (isCut[nei[f]] > 0.5);
        if (cO && cN)
        {
            kfIn[f] = 0.5*(kappaGamma[own[f]] + kappaGamma[nei[f]]);
            ++nAssigned;
        }
        else if (cO) { kfIn[f] = kappaGamma[own[f]]; ++nAssigned; }
        else if (cN) { kfIn[f] = kappaGamma[nei[f]]; ++nAssigned; }
        else { ++nNoCutSide; }   // keeps the interpolated value
    }

    // --- coupled patches: the neighbour cell lives on another rank ---------
    scalarField kappaGammaNei, isCutNei;
    syncTools::swapBoundaryCellList(mesh, kappaGamma, kappaGammaNei);
    syncTools::swapBoundaryCellList(mesh, isCut, isCutNei);

    surfaceScalarField::Boundary& kfBf = kappaFace.boundaryFieldRef();
    const surfaceScalarField::Boundary& snBf = snA.boundaryField();

    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }

        fvsPatchScalarField& kfp = kfBf[patchI];
        const labelUList& fc = patch.faceCells();
        const label bStart = patch.start() - mesh.nInternalFaces();

        forAll(kfp, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            ++nActive;

            const bool cO = (isCut[fc[i]] > 0.5);
            const bool cN = (isCutNei[bStart + i] > 0.5);
            if (cO && cN)
            {
                // both ranks compute the same mean -> seam-consistent
                kfp[i] = 0.5*(kappaGamma[fc[i]] + kappaGammaNei[bStart + i]);
                ++nAssigned;
            }
            else if (cO) { kfp[i] = kappaGamma[fc[i]]; ++nAssigned; }
            else if (cN) { kfp[i] = kappaGammaNei[bStart + i]; ++nAssigned; }
            else { ++nNoCutSide; }
        }
    }

    reduce(nCut, sumOp<label>());
    reduce(nActive, sumOp<label>());
    reduce(nAssigned, sumOp<label>());
    reduce(nNoCutSide, sumOp<label>());
    reduce(nFootFail, sumOp<label>());
    reduce(nGaussFallback, sumOp<label>());
    Info<< "cutCellFootPointFace curvature: " << nCut << " cut cells, "
        << nAssigned << "/" << nActive << " active faces assigned ("
        << nNoCutSide << " with no cut side, " << nFootFail
        << " foot searches failed, " << nGaussFallback
        << " Gaussian fallbacks)" << endl;
}

//- CELL-MEAN delivery: keep the PER-FACE inversion, then average the corrected
//  values over the active faces that come from the same cell fit.
//
//  This is the per-face delivery followed by one averaging pass, NOT the
//  cut-cell delivery. The distinction is the whole point:
//
//    cutCellFootPointFace : ONE inversion per cell, from kappa_c, K_c and ONE
//                           foot distance d_c, applied to all n active faces.
//                           MEASURED G h^2 = 0.821, late growth rate 454 1/s.
//    cellMeanFootPointFace: n inversions, each with its OWN d_f, then averaged.
//                           n partly independent foot-distance errors are
//                           averaged instead of one being applied n times.
//
//  Both give every active face of a cut cell the same number, so both make the
//  across-support variation of kappa_f vanish within a cut cell. WP8.2 measured
//  that this structural property is NOT what governs the instability; what
//  governs it is the gain, how far kappa_f moves per unit perturbation of psi.
//  Averaging n estimates lowers the gain, concentrating on one raises it, which
//  is why these two structurally identical deliveries are expected to behave
//  very differently in the coupled run.
//
//  Accuracy cost: averaging over a cell's faces is a local smoothing over a
//  stencil of size h, whose error is O(h^2 d^2 kappa) -- second order, so it
//  does not degrade the h^2 convergence on an interface whose curvature varies.
//  A static circle cannot show this either way (its exact curvature is
//  constant), so the accuracy claim must be checked on a shape with
//  non-constant curvature.
//
//  Faces with no cut-cell side keep their per-face corrected value (counted).
inline void computeCellMeanFootPointFaceCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& alpha,
    const volScalarField& kappa,
    const slReconstruction& recon,
    surfaceScalarField& kappaFace
)
{
    // 1. the per-face inversion, seam-synchronized (this is the current
    //    production delivery, unchanged)
    computeStabilizedFootPointFaceCurvature
    (
        mesh, psi, alpha, kappa, recon, kappaFace
    );

    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snIn = snA.primitiveField();
    const surfaceScalarField::Boundary& snBf = snA.boundaryField();
    const scalar snMax = max(gMax(mag(snIn)()), VSMALL);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& aIn = alpha.primitiveField();

    scalarField& kfIn = kappaFace.primitiveFieldRef();
    surfaceScalarField::Boundary& kfBf = kappaFace.boundaryFieldRef();

    // 2. accumulate each cell's corrected ACTIVE-face values. Coupled-patch
    //    faces are included so a cut cell straddling a processor boundary
    //    averages over the same face set it would in serial.
    scalarField fSum(mesh.nCells(), 0.0);
    scalarField fCnt(mesh.nCells(), 0.0);

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        fSum[own[f]] += kfIn[f]; fCnt[own[f]] += 1.0;
        fSum[nei[f]] += kfIn[f]; fCnt[nei[f]] += 1.0;
    }
    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }
        const labelUList& fc = patch.faceCells();
        forAll(fc, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            fSum[fc[i]] += kfBf[patchI][i]; fCnt[fc[i]] += 1.0;
        }
    }

    // 3. one averaged interface curvature per CUT cell (0 < alpha < 1: the
    //    geometric indicator's own definition, identical in 2D and 3D)
    scalarField cellMean(mesh.nCells(), 0.0);
    scalarField isCut(mesh.nCells(), 0.0);
    label nCut = 0;

    forAll(aIn, c)
    {
        if (aIn[c] <= SMALL || aIn[c] >= 1 - SMALL) { continue; }
        if (fCnt[c] < 0.5) { continue; }
        cellMean[c] = fSum[c]/fCnt[c];
        isCut[c] = 1.0;
        ++nCut;
    }

    // 4. assign, exactly as the cut-cell delivery does
    label nActive = 0, nAssigned = 0, nNoCutSide = 0;

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        ++nActive;

        const bool cO = (isCut[own[f]] > 0.5);
        const bool cN = (isCut[nei[f]] > 0.5);
        if (cO && cN) { kfIn[f] = 0.5*(cellMean[own[f]] + cellMean[nei[f]]); ++nAssigned; }
        else if (cO)  { kfIn[f] = cellMean[own[f]]; ++nAssigned; }
        else if (cN)  { kfIn[f] = cellMean[nei[f]]; ++nAssigned; }
        else          { ++nNoCutSide; }   // keeps its per-face corrected value
    }

    scalarField cellMeanNei, isCutNei;
    syncTools::swapBoundaryCellList(mesh, cellMean, cellMeanNei);
    syncTools::swapBoundaryCellList(mesh, isCut, isCutNei);

    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }

        fvsPatchScalarField& kfp = kfBf[patchI];
        const labelUList& fc = patch.faceCells();
        const label bStart = patch.start() - mesh.nInternalFaces();

        forAll(kfp, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            ++nActive;

            const bool cO = (isCut[fc[i]] > 0.5);
            const bool cN = (isCutNei[bStart + i] > 0.5);
            if (cO && cN)
            {
                // both ranks form the same mean -> the seam stays consistent
                kfp[i] = 0.5*(cellMean[fc[i]] + cellMeanNei[bStart + i]);
                ++nAssigned;
            }
            else if (cO) { kfp[i] = cellMean[fc[i]]; ++nAssigned; }
            else if (cN) { kfp[i] = cellMeanNei[bStart + i]; ++nAssigned; }
            else { ++nNoCutSide; }
        }
    }

    reduce(nCut, sumOp<label>());
    reduce(nActive, sumOp<label>());
    reduce(nAssigned, sumOp<label>());
    reduce(nNoCutSide, sumOp<label>());
    Info<< "cellMeanFootPointFace curvature: " << nCut << " cut cells, "
        << nAssigned << "/" << nActive << " active faces assigned the mean of "
        << "their cell's per-face inversions (" << nNoCutSide
        << " with no cut side, kept per-face)" << endl;
}


//- SYMMETRIC FACE-MEAN delivery: the per-face inversion, then a smoothing pass
//  that is symmetric ABOUT EACH FACE.
//
//  Measured motivation (plan sec. 12 and 13): lowering the gain helps the coupled
//  run, but every construction that lowered it so far did so by giving all active
//  faces of a cut cell ONE value, and that lumping is first order on any interface
//  whose curvature varies -- a face at the edge of a cell receives a value centred
//  on the CELL, an O(h) offset. Measured on the 2:1 ellipse gate: order 1.02
//  (cut-cell) and 1.03 (cell-mean) against 1.98 for the per-face inversion.
//
//  The distinction that matters is not "averaging" versus "not averaging", it is
//  whether the average is CENTRED ON THE FACE. For face f between cells P and N,
//  the neighbour set
//
//      {active faces of P} + {active faces of N} - {f}
//
//  is symmetric about x_f by construction, because P and N are mirror images
//  through f. A mean over a symmetric set annihilates linear fields exactly on a
//  uniform mesh, so the smoothing error is O(h^2 d^2 kappa) -- second order
//  preserved -- while the noise is reduced by the averaging. That is precisely
//  what the cell-centred assignment cannot do.
//
//      kappa_f <- (1 - theta) kappa_f + theta * mean(neighbour set)
//
//  theta = 0 recovers the per-face inversion exactly; theta = 1 is the pure
//  symmetric neighbour mean.
//
//  Parallel: the neighbour set of an INTERNAL face is rank-local (both adjacent
//  cells are owned). For a coupled-patch face the far cell lives on another rank,
//  so each rank contributes its own side's (sum, count) and the two are combined
//  after a swap -- both ranks then form the same mean and the seam stays
//  consistent, the same discipline as the inversion itself.
inline void computeSymmetricFaceMeanCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& alpha,
    const volScalarField& kappa,
    const slReconstruction& recon,
    surfaceScalarField& kappaFace,
    const scalar theta = 0.5
)
{
    // 1. the per-face inversion, seam-synchronized (unchanged production path)
    computeStabilizedFootPointFaceCurvature
    (
        mesh, psi, alpha, kappa, recon, kappaFace
    );

    if (theta <= 0) { return; }

    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snIn = snA.primitiveField();
    const surfaceScalarField::Boundary& snBf = snA.boundaryField();
    const scalar snMax = max(gMax(mag(snIn)()), VSMALL);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const scalarField& kfIn = kappaFace.primitiveField();
    const surfaceScalarField::Boundary& kfBf0 = kappaFace.boundaryField();

    // 2. per-cell sum and count over that cell's ACTIVE faces
    scalarField cSum(mesh.nCells(), 0.0);
    scalarField cCnt(mesh.nCells(), 0.0);

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        cSum[own[f]] += kfIn[f]; cCnt[own[f]] += 1.0;
        cSum[nei[f]] += kfIn[f]; cCnt[nei[f]] += 1.0;
    }
    forAll(kfBf0, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }
        const labelUList& fc = patch.faceCells();
        forAll(fc, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            cSum[fc[i]] += kfBf0[patchI][i]; cCnt[fc[i]] += 1.0;
        }
    }

    scalarField cSumNei, cCntNei;
    syncTools::swapBoundaryCellList(mesh, cSum, cSumNei);
    syncTools::swapBoundaryCellList(mesh, cCnt, cCntNei);

    // 3. the symmetric neighbour mean, then the blend. Written into a copy so
    //    every face averages the UNSMOOTHED field (a single Jacobi sweep, not a
    //    sweep whose result depends on face ordering).
    scalarField kfNew(kfIn);
    label nActive = 0, nSmoothed = 0, nNoNeighbour = 0;

    forAll(kfIn, f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        ++nActive;

        const scalar s = (cSum[own[f]] - kfIn[f]) + (cSum[nei[f]] - kfIn[f]);
        const scalar n = (cCnt[own[f]] - 1.0) + (cCnt[nei[f]] - 1.0);
        if (n < 0.5) { ++nNoNeighbour; continue; }

        kfNew[f] = (1.0 - theta)*kfIn[f] + theta*(s/n);
        ++nSmoothed;
    }

    surfaceScalarField::Boundary& kfBf = kappaFace.boundaryFieldRef();
    scalarField bNew(mesh.nBoundaryFaces(), 0.0);
    boolList bSet(mesh.nBoundaryFaces(), false);

    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }
        const labelUList& fc = patch.faceCells();
        const label bStart = patch.start() - mesh.nInternalFaces();

        forAll(fc, i)
        {
            if (mag(snBf[patchI][i]) <= 1e-8*snMax) { continue; }
            ++nActive;

            const scalar kf = kfBf0[patchI][i];
            const scalar s = (cSum[fc[i]] - kf)
                           + (cSumNei[bStart + i] - kf);
            const scalar n = (cCnt[fc[i]] - 1.0)
                           + (cCntNei[bStart + i] - 1.0);
            if (n < 0.5) { ++nNoNeighbour; continue; }

            // both ranks form the same s and n -> the seam value agrees exactly
            bNew[bStart + i] = (1.0 - theta)*kf + theta*(s/n);
            bSet[bStart + i] = true;
            ++nSmoothed;
        }
    }

    kappaFace.primitiveFieldRef() = kfNew;
    forAll(kfBf, patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        if (!patch.coupled()) { continue; }
        const label bStart = patch.start() - mesh.nInternalFaces();
        fvsPatchScalarField& kfp = kfBf[patchI];
        forAll(kfp, i)
        {
            if (bSet[bStart + i]) { kfp[i] = bNew[bStart + i]; }
        }
    }

    reduce(nActive, sumOp<label>());
    reduce(nSmoothed, sumOp<label>());
    reduce(nNoNeighbour, sumOp<label>());
    Info<< "symmetricFaceMean curvature: theta = " << theta << ", "
        << nSmoothed << "/" << nActive << " active faces smoothed over their "
        << "owner+neighbour active-face set (" << nNoNeighbour
        << " with no neighbour, left per-face)" << endl;
}


} // End namespace Foam

#endif

// ************************************************************************* //
