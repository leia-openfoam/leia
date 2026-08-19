/*---------------------------------------------------------------------------*\
    CELL-CENTRE parallel-surface inverse of the curvature, applied IN PLACE to
    the cell field. Not a face delivery: the transformed cell curvature is then
    handed to whatever consumes the cell field, and with
    faceCurvatureSource = model that is plain arithmetic interpolation.

    Per cell:
        kappa_c^Gamma = parallelSurfaceInverse(kappa_c, d_c, K_c)
    with kappa_c and K_c from the cell's own quadratic fit AT ITS CENTRE and
    d_c the signed offset from the cell centre to the fit's zero set. The
    algebra is the SAME dimension-seamless expression the face delivery uses
    (K = 0 identically in 2D, active in 3D):
        kappa^Gamma = (kappa - 2 K d) / (1 - d kappa + K d^2)

    WHY THIS AND NOT cutCellFootPointFace. That delivery computes the identical
    inverse but then assigns the ONE cell value to every active face of the cell
    -- a lumping, which discards precisely the cell-to-cell smoothness this
    construction exists to provide, and which is first order wherever curvature
    varies along the interface. Here the inverted CELL field is interpolated
    normally, so its smoothness survives to the faces.

    WHY signedOffset AND NOT footPointDistance for d_c. footPointDistance is the
    stabilized foot-point search: iterative, trust-region guarded, and it
    correctly REFUSES cells whose zero set lies outside the fit's trusted
    stencil -- at N=512 that is 261232 of 262144 cells, i.e. everything but the
    interface band. Leaving those uninverted makes the cell field discontinuous
    at the band edge, and snGrad across that edge reports the jump rather than
    the curvature. signedOffset is the algebraic stable quadratic root along the
    normal ray, d_c = 2 psi_c/(|g| + sqrt(D)) with D = |g|^2 - 2 psi_c h_nn,
    falling back to first order when the discriminant is small; it always
    returns, so the inverted field is defined in EVERY cell and no edge exists.
    It is also the single offset definition the nSL transport and the curvature
    offsetCorrection already share, so transport and curvature cannot disagree
    about a cell's offset.

    MEASURED MOTIVATION (leiaTestRemainderTerm.csv, exact-SDF circle, where the
    exact grad(kappa) is identically zero so the whole term is error). The
    non-gradient content of the capillary force is alpha_f*snGrad(kappa_c), and
    over the FULL CSF force support:

        cell curvature                 N=64    N=128    N=256    N=512   order
        quadraticCellCentre           5.14e5   5.03e5   4.56e5   4.46e5   +0.09
        closestPointNewton            5.11e5   5.02e5   4.56e5   4.46e5   +0.09
        newtonFoot                    5.09e5   5.02e5   4.56e5   4.46e5   +0.09
        THIS (signedOffset)             9647     2273    539.7    139.3   +2.01
        THIS (stabilized foot)          9095     2192    520.7    133.1   +2.02

    Second order against ~zeroth for every other cell curvature, and 3200x
    smaller at N=512. Both offset definitions agree, so the convergence is a
    property of the algebraic inversion rather than of how d_c is obtained. With
    the GEOMETRIC face alpha (the rhoLENT mass-flux alpha_f) it reads 149.1 at
    order +2.016 -- an 8% cost, convergence intact.

    CAVEAT, stated because the coupled prediction rests on it: this is measured
    on constant curvature, where the exact remainder vanishes. On a
    varying-curvature interface the exact remainder is genuinely non-zero and
    the cut-cell construction is known to be FIRST order there (the ellipse
    gates score cutCellInverse at 1.02). So the second-order behaviour may be
    specific to constant curvature -- which is the stationary droplet, but not
    the oscillating one.
\*---------------------------------------------------------------------------*/

#ifndef cellCentreInverseCurvature_H
#define cellCentreInverseCurvature_H

namespace Foam
{

inline void applyCellCentreInverseCurvature
(
    const fvMesh& mesh,
    const slReconstruction& recon,
    volScalarField& kappa,
    const bool useGaussianCurvature = true
)
{
    scalarField& k = kappa.primitiveFieldRef();

    label nInv = 0, nFirstOrder = 0, nNoFit = 0, nUnfilled = 0;
    forAll(k, c)
    {
        // NEVER INVERT A CELL THE CURVATURE FILL SKIPPED.
        //
        // meanCurvatureNoExtension zeroes the whole field and then fills only
        // cells that have a full-rank quadratic fit, a non-degenerate gradient,
        // and satisfy the band gate |psi_c|/|g| <= 3*stencilRadius. Every other
        // cell keeps EXACTLY 0. Feeding such a cell to the inverse gives
        //     parallelSurfaceInverse(0, d, K) = -2 K d/(1 + K d^2),
        // which is identically 0 in 2D (K vanishes there) but in 3D manufactures
        // a curvature out of nothing: for a sphere K = 1/(R+d)^2, so the
        // magnitude peaks at 1/R -- HALF the exact 2/R, with the wrong sign --
        // around d ~ R. This header used to claim "it always returns, so the
        // inverted field is defined in EVERY cell and no edge exists"; that is
        // false in 3D. The fill's edge is still there, and inverting across it
        // decorated the edge with a large negative curvature instead of leaving
        // the zero the fill produced, i.e. it made the jump WORSE.
        //
        // It is currently benign for the force, because the Detrixhe-Aslam alpha
        // is exactly 0/1 that far out so snGrad(alpha) = 0 there. It is not
        // benign for the written kappa field, for any diagnostic that norms over
        // more than the active faces, or for any change that widens the force
        // support (forceWeight sharpHeaviside, or a mollified indicator).
        //
        // An unfilled cell is exactly 0.0; a fitted curvature is not, so testing
        // against VSMALL identifies the fill's own gate without duplicating it.
        if (mag(k[c]) <= VSMALL)
        {
            ++nUnfilled;
            continue;
        }

        // useGaussianCurvature = false forces K == 0, i.e. runs 3D through the
        // 2D algebra of the same expression. The fit is then skipped entirely,
        // so its noise cannot enter through any other path either.
        bool okK = false;
        const scalar Kc =
            useGaussianCurvature ? fitGaussianCurvature(recon, c, okK) : scalar(0);

        bool fallback = false;
        const scalar dc = recon.signedOffset(c, fallback);
        if (!std::isfinite(dc))
        {
            ++nNoFit;             // leave the cell-centre value untouched
            continue;
        }
        if (fallback) { ++nFirstOrder; }

        k[c] = parallelSurfaceInverse
        (
            k[c], dc, (useGaussianCurvature && okK) ? Kc : scalar(0)
        );
        ++nInv;
    }
    kappa.correctBoundaryConditions();

    reduce(nInv, sumOp<label>());
    reduce(nFirstOrder, sumOp<label>());
    reduce(nNoFit, sumOp<label>());
    reduce(nUnfilled, sumOp<label>());
    // nNoFit counts only non-finite offsets, which signedOffset essentially never
    // returns -- so it used to read 0 while a large unfilled population was being
    // inverted silently. nUnfilled is the number that matters and is now reported.
    Info<< "cellCentreInverse curvature: " << nInv << " cells inverted ("
        << nFirstOrder << " on the first-order offset branch, "
        << nUnfilled << " skipped as unfilled by the curvature gate, "
        << nNoFit << " skipped on a non-finite offset)" << endl;
}

} // End namespace Foam

#endif

// ************************************************************************* //
