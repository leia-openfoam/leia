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
    Foam::defectCorrectedIDWReconstruction

Description
    Second-order departure-point reconstruction by ITERATIVELY defect-corrected
    inverse-distance-weighted (IDW) interpolation.

    A single application of the gradient (Taylor) correction --
    psi(x_d) = sum_C w_C[psi_C + grad(psi)_C.(x_d - x_C)] / sum_C w_C -- is
    linear-exact as a STATIC interpolant, but it composes a bounded Shepard
    smoother with an INDEPENDENTLY computed (pointCellsLeastSquares) gradient:
    the gradient lift is never constrained by the reconstructed values, so under
    the repeated feedback of semi-Lagrangian transport the composite operator
    amplifies grid-scale slope every step and |grad(psi)| runs away. The defect
    correction only lifts the order ROBUSTLY when it is applied iteratively, with
    the gradient fed through a stable averaging operator rather than once.

    Per arrival cell c, over its cell-point-cell (CPC) stencil S(c) = {c} U {point
    neighbours} (centres stencilC_[c], values stencilPsi_[c], cell-centred LSQ
    gradients grad(psi)_C = stencilGrad_[c]), we solve for EFFECTIVE NODAL VALUES
    u_C = psi_C + delta_C such that the gradient-lifted blend REPRODUCES the
    stencil data at every node,

        sum_C A_{DC} (u_C + grad(psi)_C.(x_D - x_C)) = psi_D    for all D in S(c),

    where A is a FINITE-weight, row-stochastic averaging matrix built from a
    nugget-regularised weight

        wreg(r) = ( h_c^2 / (|r|^2 + regEps^2 h_c^2) )^{p/2},   h_c = radius_[c],
        A_{DC}  = wreg(x_D - x_C) / sum_E wreg(x_D - x_E).

    Regularisation is essential: with the SINGULAR weight |r|^{-p} the blend
    interpolates its own nodes (A = I) and the nodal defect collapses to zero, so
    the iteration degenerates to the broken single pass. The regularised A has a
    strictly dominant diagonal (A_DD > 1/2 by construction on the CPC stencils),
    so the damped Gauss-Seidel sweep

        delta_D <- delta_D + omega ( psi_D - sum_C A_{DC}(psi_C + delta_C) - b_D )
                                    / A_DD,     b_D = sum_C A_{DC} grad_C.(x_D-x_C)

    is a contraction; the arrival-node defect delta_0 is pinned to 0. Driving the
    nodal defect down makes the reconstruction INTERPOLATORY (data-reproducing) and
    hence a bounded, non-amplifying operator, like the stable quadratic WLSQ fit.

    The FINAL evaluation at the foot x_d uses the ORIGINAL SINGULAR weights on the
    converged corrected Taylor values,

        psi(x_d) = sum_C W_C [ psi_C + delta_C + grad_C.(x_d - x_C) ] / sum_C W_C,
        W_C = |x_d - x_C|^{-p},

    which keeps linear-exactness (delta = 0 on a linear field) and the coincident-
    node contract evaluateRaw(c, x_c) == psi_c (the self node, i = 0, has
    delta_0 = 0). The blend matrix is geometry-only; it is recomputed into a reused
    per-cell buffer each step (NOT cached m*m per cell) to bound 3D memory.

    Selected with `levelSet.semiLagrangian.reconstruction defectCorrectedIDW`.
    Sub-dict keys (semiLagrangian): idwPower p (default 2), dcIters Gauss-Seidel
    sweeps (default 3), dcMaxIters hard cap (default 10), dcTol range-relative
    residual early-out (default 1e-3), dcRelax under-relaxation omega<=1 (default
    0.8), dcRegEps nugget regEps (default 0.1).

SourceFiles
    defectCorrectedIDWReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef defectCorrectedIDWReconstruction_H
#define defectCorrectedIDWReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class defectCorrectedIDWReconstruction
:
    public slReconstruction
{
    //- Cell-centred LSQ gradient of psiOld (fvSchemes key "gradPsi").
    volVectorField gradPsi_;

    //- grad(psiOld) gathered over the CPC stencil (incl. remote neighbours),
    //  aligned with stencilC_/stencilPsi_; refreshed each step in update().
    List<List<vector>> stencilGrad_;

    //- IDW power p (final-blend weight = |x_d - x_C|^{-p}); key "idwPower".
    scalar power_;

    //- Per-cell converged effective-nodal-value defects delta_C, aligned with
    //  stencilPsi_[c]; delta_[c][0] (arrival node) is pinned to 0. Filled each
    //  step in update(), consumed by evaluateRaw().
    List<scalarList> delta_;

    //- Per-cell stencil size m == stencilC_[c].size() (cached once).
    labelList m_;

    //- Largest stencil size (sizes the reused per-cell work buffers).
    label maxM_;

    //- One-off sizing of m_/delta_/maxM_ (mesh is static within a run).
    bool built_;

    //- Gauss-Seidel controls (read once from the semiLagrangian sub-dict).
    label nIters_;      // dcIters,    default 3
    label maxIters_;    // dcMaxIters, default 10
    scalar dcTol_;      // dcTol,      default 1e-3 (range-relative)
    scalar omega_;      // dcRelax,    default 0.8 (clamped <= 1)
    scalar regEps_;     // dcRegEps,   default 0.1 (nugget)

public:

    TypeName("defectCorrectedIDW");

    explicit defectCorrectedIDWReconstruction(const fvMesh& mesh);

    virtual ~defectCorrectedIDWReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
