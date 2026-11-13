/*---------------------------------------------------------------------------*\
Class
    Foam::footAnchoredQuadraticWeightedLeastSquaresReconstruction

Description
    uncachedQuadraticWeightedLeastSquares with the inverse-distance weights
    measured from the DEPARTURE FOOT POINT instead of the arrival cell centre.

    The production fit weights every stencil point by 1/|x_i - x_c| and then
    evaluates the polynomial at x_d = x_c - u dt, which is up to one cell away.
    A weighted least-squares fit is most accurate near the centroid of its
    weights, so the production fit is weighted to be best where it is NOT
    evaluated. This model moves the weights to x_d:

        w_i = 1/max(|x_i - x_d|, SMALL)        instead of  1/max(|x_i - x_c|, SMALL)

    ONLY the weights move. The polynomial stays anchored at the cell centre and
    is still fitted on the differences psi_i - psi_c over the constant-free
    basis, so R_c(x_c) = psi_c remains exact by construction, the zero-velocity
    identity of config/gates/slZeroVelocityNull2D.yaml is untouched, and
    leiaTestSLReconstruction gate (a) still holds.

    PREDICTED TO BE A NULL ON UNIFORM HEXAHEDRA, before it was written
    (docs/analysis-amplification-sl/ms.tex). The exact Fourier symbol gives
    rho = 1 for BOTH anchors at every Courant number from 0.0164 to 0.95, in 2D
    and 3D; Lambda differs (1.509 against 1.357 at Courant 0.5 along x) but rho
    does not move at all. At the boundary stencils the two anchors agree to six
    decimal places. So this model is expected to change nothing where the symbol
    applies; it exists to be measured where the symbol does NOT apply -- graded
    and polyhedral meshes, where the stencil is not symmetric about either point.

    KNOWN DEFECT, stated because it is structural and not a bug. The weight has a
    POLE where the foot coincides with a stencil point: |x_i - x_d| -> 0. The
    SMALL floor caps it at 1/SMALL, which makes the fit interpolate that one
    datum essentially exactly. On a uniform mesh that happens at Courant 1 (a
    cell centre) and, at a boundary cell, at Courant 0.5 (the boundary face at
    h/2). The production cell-centre anchor has no such pole.

    The per-cell fit therefore CANNOT be cached across the step: the weights
    depend on the foot point, which differs per cell and per step. The fit is
    assembled and solved inside evaluateRaw, once per cell per evaluation. The
    base class's update() still computes the cell-centre-anchored coefficients,
    and every geometry query that reads them -- fitDerivatives, meanCurvature,
    signedOffset, footPointDistance -- is deliberately left on that fit, because
    those are evaluated AT the cell centre and have no foot point.

    The quadraticPivotTol demotion stays a once-per-mesh CELL-CENTRE decision in
    build(), so the fitted rank is a property of the mesh and not of the flow.

SourceFiles
    footAnchoredQuadraticWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef footAnchoredQuadraticWeightedLeastSquaresReconstruction_H
#define footAnchoredQuadraticWeightedLeastSquaresReconstruction_H

#include "uncachedQuadraticWeightedLeastSquaresReconstruction.H"

namespace Foam
{

class footAnchoredQuadraticWeightedLeastSquaresReconstruction
:
    public uncachedQuadraticWeightedLeastSquaresReconstruction
{
public:

    TypeName("footAnchoredQuadraticWeightedLeastSquares");

    explicit footAnchoredQuadraticWeightedLeastSquaresReconstruction
    (
        const fvMesh& mesh
    );

    virtual ~footAnchoredQuadraticWeightedLeastSquaresReconstruction() = default;

    //- Assemble and solve the fit with the weights anchored at x, then evaluate
    //  the polynomial there. Falls back to the base (cell-centre) behaviour when
    //  the local solve fails, so a degenerate cell can never be worse than the
    //  production model.
    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
