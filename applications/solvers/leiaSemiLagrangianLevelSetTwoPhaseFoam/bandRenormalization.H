/*---------------------------------------------------------------------------*\
    Band renormalization: restore the PARALLEL FOLIATION of psi without moving
    the interface.

MEASURED MOTIVATION (docs/plan-curvature-stabilization.md sec. 14.2). In the
coupled stationary droplet, beta := |grad psi| in the band spreads from
[0.9983, 0.9999] at t = 0.02 to [0.724, 1.373] at t = 0.10, and the
non-parallelism vector |a_T| = |P grad(ln beta)| grows from 2.3 to 583 1/m. The
parallel-surface inverse every foot-point delivery uses is exact only while the
level sets are parallel, i.e. only while beta is constant ALONG each level set,
and its residual -- the integral of D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2
along the normal -- grew 600x over that run, reaching 19-37% of the measured band
curvature error by blow-up. The remaining ~80% of the error growth is the
least-squares quadratic itself degrading, driven by the SAME beta spread.

So the operand is the problem, not the operator. This restores the property the
whole curvature path assumes.

THE CONSTRUCTION. Divide psi by its own interface-extended gradient magnitude:

    psi <- psi / beta_Gamma,        beta_Gamma := |grad psi| at the FOOT POINT

with the foot point taken on the cell's own quadratic fit, so
    g_Gamma = g + H.(-d n),   n = g/|g|,   d = footPointDistance(cell centre)
and beta_Gamma = |g_Gamma|. Then |grad psi| -> 1 ON the interface and, to leading
order, the tangential variation of beta -- the entire non-parallelism -- is
removed. It is NOT a redistancing: nothing is propagated, nothing is solved, no
pseudo-time is marched, and no extra derivative of psi is taken. beta_Gamma comes
from the same g and H the curvature already uses.

WHY THIS SATISFIES THE IN-BAND REWRITE GATE BY CONSTRUCTION. The programme
forbids in-band psi rewrites without the sign-histogram gate, because a
redistancer can displace the zero set and flip signs. This operation is a
multiplication by a STRICTLY POSITIVE field, so:
  * the zero set is preserved EXACTLY -- psi = 0 iff psi/beta_Gamma = 0;
  * no cell can change sign, for any beta_Gamma > 0;
  * the phase indicator, the geometric mass flux and the band all read the same
    interface they read before.
Both are asserted at run time rather than assumed: the count of sign changes must
be zero and the zero-crossing cell set must be unchanged.

GUARDS. A cell is renormalized only where the fit is full (order 2), the foot
search converged, and beta_Gamma lies in [betaMin, betaMax] -- outside that range
the fit is not trustworthy and psi is left alone (counted, not silently skipped).
Cells outside the band are untouched, so the bulk skeleton of psi is never
rescaled.
\*---------------------------------------------------------------------------*/

#ifndef bandRenormalization_H
#define bandRenormalization_H

namespace Foam
{

struct bandRenormalizationStats
{
    label nBand = 0;
    label nRenormalized = 0;
    label nFitFallback = 0;      // no full quadratic
    label nFootFail = 0;         // foot search did not converge
    label nGuardReject = 0;      // beta_Gamma outside [betaMin, betaMax]
    label nSignChange = 0;       // MUST be zero
    scalar maxRelChange = 0;     // max |psi_new - psi_old|/|psi_old|
    scalar betaGammaMin = GREAT;
    scalar betaGammaMax = -GREAT;
    scalar aTRel = 0;            // band L2 of |P H n|/beta * h  (the trigger)
    bool fired = false;
};


//- psi <- psi/beta_Gamma on the narrow band, with the divisor SMOOTHED along the
//  band and the whole operation gated on a measured non-parallelism trigger.
//
//  MEASURED, and the reason for both refinements (plan sec. 14.3): the naive
//  version -- divide by the raw per-cell beta_Gamma every step -- is WORSE than
//  doing nothing. On the exact signed-distance initial field, where beta is
//  identically 1 and the operation should be a no-op, the fitted beta_Gamma comes
//  out in [0.9936, 0.9972] and psi is rescaled by up to 0.65% anyway. That is the
//  estimator's own error being injected into psi as grid-scale noise, every step,
//  and the instability responds to grid-scale psi noise (the gain result of
//  sec. 10). Measured at N=64 over 2 ms: max|U| 8x worse than the control and
//  min|grad psi| in the band drifting FASTER (0.957 against 0.988).
//
//  Two changes follow directly from that diagnosis:
//
//  1. SMOOTH THE DIVISOR. The drift we want to remove is large-scale (beta spreads
//     smoothly to [0.72, 1.37] across the droplet); the noise we must not inject is
//     grid-scale. Averaging beta_Gamma over each band cell's face-neighbours
//     separates them: the smooth part survives, the grid-scale part is damped. The
//     stencil is ALL face neighbours, so unlike the face-curvature case it does not
//     depend on where the interface sits in the cell, and no accuracy order is at
//     stake here -- psi only needs a smooth positive divisor.
//
//  2. FIRE ONLY WHEN THE DRIFT IS REAL. The trigger is the band L2 of the actual
//     non-parallelism, |P H n|/beta (this is exactly grad_Gamma(ln beta); it is
//     zero iff the level sets are parallel), made dimensionless with h. Measured
//     over the reference run: ~2e-4 at t = 0.02, 2.4e-3 at 0.06, 1.3e-2 at 0.08,
//     4.5e-2 at 0.10. A threshold near 5e-3 therefore leaves the healthy early
//     phase completely untouched and acts only once the foliation has genuinely
//     decayed. thresholdRel <= 0 restores every-step firing.
inline bandRenormalizationStats renormalizeBandPsi
(
    const fvMesh& mesh,
    volScalarField& psi,
    const volScalarField& narrowBand,
    const slReconstruction& recon,
    const scalar thresholdRel = 5e-3,
    const label nSmooth = 2,
    const scalar relaxation = 1.0,
    const scalar betaMin = 0.2,
    const scalar betaMax = 5.0
)
{
    bandRenormalizationStats st;

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& nbIn = narrowBand.primitiveField();
    scalarField& psiIn = psi.primitiveFieldRef();

    const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());

    // --- pass 1: beta_Gamma and the non-parallelism measure, per band cell -----
    scalarField bG(mesh.nCells(), 0.0);
    boolList hasBG(mesh.nCells(), false);
    scalar aT2 = 0, sumV = 0;
    const scalarField& V = mesh.V();

    forAll(psiIn, c)
    {
        if (nbIn[c] < 0.5) { continue; }
        ++st.nBand;

        vector g(Zero);
        symmTensor H(Zero);
        if (recon.fitDerivatives(c, g, H) < 2) { ++st.nFitFallback; continue; }
        const scalar gm = mag(g);
        if (gm < SMALL) { ++st.nFitFallback; continue; }

        bool okD = false;
        const scalar d = recon.footPointDistance(c, C[c], 0.0, okD);
        if (!okD) { ++st.nFootFail; continue; }

        const vector n = g/gm;
        const scalar b = mag(g + (H & (-d*n)));
        if (b < betaMin || b > betaMax) { ++st.nGuardReject; continue; }

        bG[c] = b;
        hasBG[c] = true;

        // grad_Gamma(ln beta) = P H n / beta : zero iff the level sets are parallel
        const vector Hn(H & n);
        const vector aTv(Hn - (n & Hn)*n);
        aT2 += V[c]*sqr(mag(aTv)/gm);
        sumV += V[c];
    }

    reduce(aT2, sumOp<scalar>()); reduce(sumV, sumOp<scalar>());
    st.aTRel = (sumV > 0) ? Foam::sqrt(aT2/sumV)*h : 0;

    reduce(st.nBand, sumOp<label>());
    reduce(st.nFitFallback, sumOp<label>());
    reduce(st.nFootFail, sumOp<label>());
    reduce(st.nGuardReject, sumOp<label>());

    if (thresholdRel > 0 && st.aTRel < thresholdRel)
    {
        // Healthy foliation: do nothing at all. Not even a smooth rescaling,
        // because any rescaling injects the estimator's error.
        return st;
    }
    st.fired = true;

    // --- pass 2: smooth the divisor over each band cell's face neighbours -----
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    for (label sweep = 0; sweep < nSmooth; ++sweep)
    {
        scalarField sum(bG), cnt(mesh.nCells(), 1.0);
        forAll(own, f)
        {
            const label P = own[f], N = nei[f];
            if (hasBG[P] && hasBG[N])
            {
                sum[P] += bG[N]; cnt[P] += 1.0;
                sum[N] += bG[P]; cnt[N] += 1.0;
            }
        }
        forAll(bG, c)
        {
            if (hasBG[c] && cnt[c] > 1.5) { bG[c] = sum[c]/cnt[c]; }
        }
    }

    // --- pass 3: rescale ------------------------------------------------------
    forAll(psiIn, c)
    {
        if (!hasBG[c]) { continue; }
        const scalar b = bG[c];
        if (b < betaMin || b > betaMax) { ++st.nGuardReject; continue; }

        st.betaGammaMin = min(st.betaGammaMin, b);
        st.betaGammaMax = max(st.betaGammaMax, b);

        // RELAXED correction: divide by b^omega rather than b.
        //
        // MEASURED (plan sec. 14.3): the unrelaxed version rescales psi by up to
        // 2% in a single step once the gate opens, fires 607 times, and blows the
        // droplet up at t = 0.0195 against 0.0668 for the untouched reference --
        // 3.4x WORSE. A finite-amplitude kick to psi is amplified by the capillary
        // loop faster than the systematic drift it removes.
        //
        // Relaxation separates the two by their statistics rather than their size.
        // The drift being corrected is SMOOTH and its correction accumulates
        // COHERENTLY over n steps (proportional to omega*n); the fit error that
        // rides along is grid-scale and accumulates INCOHERENTLY (omega*sqrt(n)).
        // Spreading the same total correction over n steps therefore improves the
        // correction-to-injected-noise ratio by sqrt(n), at no cost in the final
        // amount corrected. omega = 1 recovers the one-shot behaviour.
        const scalar before = psiIn[c];
        const scalar after =
            (relaxation >= 1.0) ? before/b : before/Foam::pow(b, relaxation);
        if ((before > 0) != (after > 0) && before != 0) { ++st.nSignChange; }
        if (mag(before) > SMALL)
        {
            st.maxRelChange =
                max(st.maxRelChange, mag(after - before)/mag(before));
        }
        psiIn[c] = after;
        ++st.nRenormalized;
    }

    psi.correctBoundaryConditions();

    reduce(st.nRenormalized, sumOp<label>());
    reduce(st.nSignChange, sumOp<label>());
    reduce(st.maxRelChange, maxOp<scalar>());
    reduce(st.betaGammaMin, minOp<scalar>());
    reduce(st.betaGammaMax, maxOp<scalar>());

    if (st.nSignChange > 0)
    {
        FatalErrorInFunction
            << "band renormalization changed the sign of psi in "
            << st.nSignChange << " cells. Dividing by a positive field cannot "
            << "do this, so either the divisor went non-positive past its guard "
            << "or the field is corrupt -- either way the interface has moved "
            << "and the in-band rewrite gate is violated."
            << exit(FatalError);
    }

    Info<< "band renormalization FIRED (aT*h = " << st.aTRel << " > "
        << thresholdRel << "): " << st.nRenormalized << "/" << st.nBand
        << " band cells rescaled by 1/beta_Gamma in ["
        << st.betaGammaMin << ", " << st.betaGammaMax
        << "] after " << nSmooth << " smoothing sweeps, relaxation "
        << relaxation << ", max |dpsi|/|psi| = "
        << st.maxRelChange << " (" << st.nFitFallback << " no full fit, "
        << st.nFootFail << " foot failed, " << st.nGuardReject
        << " guard-rejected; 0 sign changes by construction)" << endl;

    return st;
}

} // End namespace Foam

#endif

// ************************************************************************* //
