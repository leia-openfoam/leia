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
};

} // End namespace Foam

#endif

// ************************************************************************* //
