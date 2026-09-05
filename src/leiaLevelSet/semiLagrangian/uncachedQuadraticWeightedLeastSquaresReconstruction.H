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

Class
    Foam::uncachedQuadraticWeightedLeastSquaresReconstruction

Description
    CPU-BOUND, PARALLEL-SCALING sibling of quadraticWeightedLeastSquares: the
    SAME IDW-weighted quadratic value fit (same 1/|d| weights, same constant-free
    5-coeff (2D) / 9-coeff (3D) basis, same under-determined linear/constant
    fallback), hence IDENTICAL accuracy (O(h^2..3)) and stability (pure value
    fitting -- the slope is a byproduct, never an injected input) BY CONSTRUCTION.

    The difference is purely the linear algebra and where the cost lives:

      * quadraticWeightedLeastSquares caches, per cell, the ncoeff x nNbr float
        pseudo-inverse map MW_ (compute-once from the static geometry) and per
        step does a contiguous mat-vec coeffs = MW_ * ds. That mat-vec is
        MEMORY-BANDWIDTH bound (~0.25 flop/byte) and the MW_ cache is the binding
        constraint at 128^3/256^3 (it forces single precision and the release of
        the stencil centres).

      * this model caches NOTHING per cell. Each step, per cell, it assembles the
        tiny SYMMETRIC weighted normal-equations matrix M = A^T W^2 A (ncoeff x
        ncoeff; 5x5 or 9x9) and rhs g = A^T W^2 ds in one sweep over the (already
        haloed) stencil, then solves M coeffs = g by a small in-place Cholesky.
        The work is register/L1-resident (~2 flop/byte) -> COMPUTE bound and
        embarrassingly parallel over the existing stencil halo (no extra
        communication, no per-cell cache). The stencil centres are RETAINED (the
        fit is reassembled every step), so it does not call releaseStencilCentres.

    Normal equations square the condition number vs the pseudo-inverse; for the
    9x9 SPD moment matrix over an O(1)-conditioned IDW stencil this is harmless.
    A degenerate/near-singular cell (Cholesky not SPD) falls back to a constant
    reconstruction (= psi_c) for that step; a small Tikhonov ridge `ridgeEps`
    (default 0) can be added to the diagonal if needed on distorted meshes.

    Selected with `reconstruction uncachedQuadraticWeightedLeastSquares`.
    Composes with the default `correction direct` (the fit is exact per step; no
    iteration). Trades quad's per-cell memory for on-the-fly compute: the win is
    bandwidth relief + removal of the hard large-N memory ceiling + strong
    scaling, not a large raw-byte saving.

SourceFiles
    uncachedQuadraticWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef uncachedQuadraticWeightedLeastSquaresReconstruction_H
#define uncachedQuadraticWeightedLeastSquaresReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class uncachedQuadraticWeightedLeastSquaresReconstruction
:
    public slReconstruction
{
    //- Active (non-empty) mesh directions (2 in 2D, 3 in 3D).
    labelList activeDirs_;

    //- Full quadratic coefficient count = 2*nd + nd*(nd-1)/2 (5 in 2D, 9 in 3D).
    label ncoeffFull_;

    //- Per-cell neighbour count (stencil size - 1) and its max (scratch sizing).
    labelList nNbr_;
    label maxNbr_;

    //- Per-cell active coefficient count (ncoeffFull_ | nd linear | 0 constant).
    labelList ncoeff_;

    //- Fitted coefficients, ONE flat contiguous buffer [c*ncoeffFull_ + k],
    //  refreshed every step (no per-cell List handle / heap alloc).
    scalarField coeffsFlat_;

    //- One-off basis/fallback setup (mesh static within a run).
    bool built_;

    //- Optional Tikhonov ridge on the normal-equations diagonal (default 0):
    //  M(k,k) += ridgeEps * trace(M)/ncoeff. Conditions distorted-mesh solves.
    scalar ridgeEps_;

    //- Geometric admissibility of the quadratic fit (2026-09-05): the smallest
    //  Cholesky pivot of the Jacobi-scaled (unit-diagonal) weighted normal matrix,
    //  built from stencil POSITIONS alone, must exceed this or the cell uses the
    //  linear fit (a linear stencil below it keeps the cell value). Default 0.3,
    //  read off the pivot census of the meshes in use: band cells >= 0.74 on every
    //  mesh, hex >= 0.64 everywhere, polyhedral interior >= 0.50, while cfMesh's
    //  boundary-layer / grading-transition stencils form a tail from 0 to ~0.3 with
    //  an almost empty 0.3-0.4 bin. Admits a scaled condition number of ~100, i.e.
    //  stable for U dt / h < 0.1 (a cell of condition c amplifies a neighbour's
    //  error by ~c (U dt / h)^2 per step; 1e-2 admitted c ~ 1e4 and diverged at
    //  U dt / h = 0.016). 0 = never (the behaviour before the fix). Dictionary
    //  entry quadraticPivotTol.
    scalar quadPivotTol_;

    //- Write the per-cell fit order (2 quadratic, 1 linear, 0 constant) as the
    //  volScalarField slFitOrder once at build time (dictionary writeFitOrder).
    bool writeFitOrder_;

    //- Least-squares solver: "normalEquations" (default; A^T W^2 A + Cholesky,
    //  fast but squares the condition number and drops to a constant fit on a
    //  non-SPD/degenerate cell) or "householderQR" (QR of the weighted design
    //  rows W A directly -- no squaring, a bounded LS solution even when the
    //  normal equations would be near-singular, at a modest extra cost).
    const word fit_;

    //- Foot-point Newton iterations for the curvature normal-extension (default 3).
    label curvatureNewtonIters_;

    //- Iteration cap for the closest-point (KKT) Newton foot
    label closestPointIters_;

    //- Parallel-curve (offset) correction of a cell-centred curvature.
    //
    //  Curvature evaluated at a cell centre is the curvature of the level
    //  contour THROUGH THAT CENTRE, not of the interface. For parallel curves at
    //  signed normal offset d, kappa_d = kappa/(1 + d kappa), so the exact
    //  inverse is kappa = kappa_d/(1 - d kappa_d). Over a force band about three
    //  cells thick this contributes mean|d| kappa ~ (3/4) h/R, i.e. a FIRST-ORDER
    //  error that dominates the measured band norms on a clean circle.
    //
    //  Modes:
    //    none            d = 0 (default; byte-identical to the previous behaviour)
    //    psi             d = psi_c. Correct only where psi is a true signed
    //                    distance, since psi_c is then the distance to the
    //                    interface. Reinitialisation-free transport does not
    //                    maintain that.
    //    psiOverGradPsi  d = psi_c/|grad P(0)|, the first-order distance to the
    //                    zero contour of the fit. Needs no |grad psi| = 1
    //                    assumption and uses only coefficients the fit already
    //                    produced, so it is the form that survives a drifted psi.
    //    quadraticRoot   d from the stable quadratic root along the normal ray
    //                    (offsetDistance below) -- second order in the offset,
    //                    and the SAME d the normal-projected SL transport uses,
    //                    so transport and curvature never disagree about a
    //                    cell's offset (nSL design note, Sec. 5.4).
    const word offsetCorrection_;

    //- Stabilized foot-point controls (footPointDistance; the nSL Variant (a)
    //  offset engine, stable-foot-point-3d.md). All runtime-configurable:
    //    footPointTolRel      convergence tolerance, relative to
    //                         stencilRadius(c) (default 1e-6, the note's
    //                         "1e-6 in units of the geometry");
    //    footPointAlphaMax    guard interval for the parabola extrapolation
    //                         (default 20);
    //    footPointCycles      cap on tangent-plane/surfacepoint cycles
    //                         (default 10);
    //    footPointSurfIters   cap on surfacepoint Newton iterations
    //                         (default 20).
    scalar fpTolRel_;
    scalar fpAlphaMax_;
    label fpCycles_;
    label fpSurfIters_;

    //- Acceptance threshold beta for the quadratic-root offset conversion:
    //  the root is used while D = |g|^2 - 2 psi_c h_nn >= beta |g|^2, i.e.
    //  comfortably before the local model's turning point, where the
    //  conversion becomes ill-conditioned (d(d_c)/d(psi_c) -> inf). Below the
    //  threshold the first-order psi_c/|g| fallback applies (nSL note, 5.3).
    scalar offsetBeta_;

    //- Signed offset from the cell centre to the zero set of the local model
    //  along the normal ray. Restricting the fit to the ray, R(s) = psi_c
    //  + s|g| + (1/2) s^2 h_nn, the signed distance is d_c = -s* with
    //  R(s*) = 0, taking the root that reduces to psi_c/|g| as h_nn -> 0; in
    //  cancellation-free form d_c = 2 psi_c/(|g| + sqrt(D)). Exactly
    //  zero-preserving (psi_c = 0 => d_c = 0 for any g, H). fallback reports
    //  the first-order branch. NOTE the sign convention: d_c has the sign of
    //  psi_c (distance measured FROM the interface TO the cell), matching
    //  psi_c/|g| and the psi/psiOverGradPsi curvature modes.
    inline scalar offsetDistance
    (
        const scalar psiC,
        const scalar gm,
        const scalar hnn,
        bool& fallback
    ) const
    {
        const scalar D = gm*gm - 2.0*psiC*hnn;

        if (D >= offsetBeta_*gm*gm)
        {
            fallback = false;
            return 2.0*psiC/(gm + Foam::sqrt(D));
        }

        fallback = true;
        return psiC/gm;
    }

    //- Apply the parallel-curve correction to a cell-centred curvature kd,
    //  given the signed offset d of the evaluation point. Guarded: the
    //  correction is singular where the offset reaches the local centre of
    //  curvature (1 - d kd -> 0), so it is skipped once the denominator
    //  leaves a safe band.
    inline scalar offsetCorrected(const scalar kd, const scalar d) const
    {
        const scalar denom = 1 - d*kd;

        // |denom| < 1/2 means the cell centre is past half the local radius of
        // curvature: the parallel-curve map is no longer invertible there in any
        // useful sense, so fall back to the uncorrected value rather than
        // amplify.
        return (Foam::mag(denom) > 0.5) ? kd/denom : kd;
    }

    //- Per-mode offset distance for the curvature correction (H = fit Hessian
    //  of the cell, needed only by quadraticRoot).
    inline scalar curvatureOffset
    (
        const scalar psiC,
        const vector& g,
        const scalar gm,
        const symmTensor& H
    ) const
    {
        if (offsetCorrection_ == "psi")
        {
            return psiC;
        }
        else if (offsetCorrection_ == "psiOverGradPsi")
        {
            return psiC/gm;
        }
        else if (offsetCorrection_ == "quadraticRoot")
        {
            const vector n = g/gm;
            const scalar hnn = n & (H & n);
            bool fb = false;
            return offsetDistance(psiC, gm, hnn, fb);
        }

        return 0;   // "none" (call sites shortcut this case anyway)
    }

    //- Constant-free quadratic basis at displacement d (identical to quad).
    void basis(const vector& d, const label ncoeff, scalar* b) const;

    //- grad P at displacement d from the fitted coefficients (analytic: g + H d).
    void gradFromCoeffs
    (
        const scalar* cf,
        const label nc,
        const vector& d,
        vector& g
    ) const;

    //- Symmetric Hessian of P from the fitted coefficients (0 for a linear fit).
    void hessianFromCoeffs(const scalar* cf, const label nc, symmTensor& H) const;

    //- Active directions + ncoeffFull_ from mesh geometricD().
    void setupBasis();

    //- Set per-cell ncoeff_/nNbr_/maxNbr_ from the stencil sizes (no MW_ cache,
    //  no releaseStencilCentres: the centres are needed every step).
    void build();

public:

    TypeName("uncachedQuadraticWeightedLeastSquares");

    explicit uncachedQuadraticWeightedLeastSquaresReconstruction(const fvMesh& mesh);

    virtual ~uncachedQuadraticWeightedLeastSquaresReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;

    //- Interface mean curvature kappa = div(grad psi/|grad psi|) evaluated
    //  symbolically from each cell's fitted gradient + Hessian and constant-
    //  extended along the normal via a foot-point Newton step. Requires a prior
    //  update() (the coefficients + stencil are reused). Only near-interface
    //  cells with a full quadratic fit are filled; the rest are left at 0.
    virtual void meanCurvature(volScalarField& kappa) const override;

    //- SDF-assuming curvature kappa = Laplacian(psi) = tr(H) (no |grad psi|
    //  normalisation, no normal projection); for comparison with meanCurvature.
    virtual void meanCurvatureLaplacian(volScalarField& kappa) const override;

    //- As meanCurvature but evaluated at the cell centre (no foot-point Newton
    //  projection, so NOT normal-extended); for the with/without-extension study.
    virtual void meanCurvatureNoExtension(volScalarField& kappa) const override;

    //- As meanCurvature but the foot from a trust-region CLOSEST-POINT
    //  Newton (KKT system with the fit Hessian) instead of the
    //  gradient-projection flow; per-cell fallback to the cell-centre value.
    virtual void meanCurvatureClosestPoint(volScalarField& kappa) const override;

    //- Fitted gradient and Hessian at the cell centre, read from the fit
    //  coefficients. Returns the fit order available: 2 (full quadratic),
    //  1 (linear fallback, H = 0), 0 (constant fallback, g = H = 0).
    virtual label fitDerivatives
    (
        const label c,
        vector& g,
        symmTensor& H
    ) const override;

    //- Signed offset from the cell centre to the fit's zero set along the
    //  normal ray (stable quadratic root, first-order fallback; see
    //  offsetDistance). Shared by the nSL transport and the curvature
    //  offsetCorrection quadraticRoot mode.
    virtual scalar signedOffset(const label c, bool& fallback) const override;

    //- Signed distance from p to the iso-surface {R_c = level} of the cell's
    //  own quadratic model by the stabilized foot-point algorithm (see the
    //  base declaration and stable-foot-point-3d.md). Polynomial arithmetic
    //  only; the foot must stay inside the stencil trust region.
    virtual scalar footPointDistance
    (
        const label c,
        const point& p,
        const scalar level,
        bool& ok
    ) const override;

    //- Fit curvature AT the foot point of p on the cell's own quadratic zero
    //  set (see the base declaration). The foot search DUPLICATES the
    //  footPointDistance algorithm deliberately: footPointDistance is inside
    //  every measured production path and stays byte-untouched.
    virtual scalar curvatureAtFootPoint
    (
        const label c,
        const point& p,
        bool& ok
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
