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
    volScalarField& kappa
)
{
    scalarField& k = kappa.primitiveFieldRef();

    label nInv = 0, nFirstOrder = 0, nNoFit = 0;
    forAll(k, c)
    {
        bool okK = false;
        const scalar Kc = fitGaussianCurvature(recon, c, okK);

        bool fallback = false;
        const scalar dc = recon.signedOffset(c, fallback);
        if (!std::isfinite(dc))
        {
            ++nNoFit;             // leave the cell-centre value untouched
            continue;
        }
        if (fallback) { ++nFirstOrder; }

        k[c] = parallelSurfaceInverse(k[c], dc, okK ? Kc : scalar(0));
        ++nInv;
    }
    kappa.correctBoundaryConditions();

    reduce(nInv, sumOp<label>());
    reduce(nFirstOrder, sumOp<label>());
    reduce(nNoFit, sumOp<label>());
    Info<< "cellCentreInverse curvature: " << nInv << " cells inverted ("
        << nFirstOrder << " on the first-order offset branch, "
        << nNoFit << " left uninverted)" << endl;
}

} // End namespace Foam

#endif

// ************************************************************************* //
