// Band mode-SPECTRUM diagnostic (WP0 of docs/plan-curvature-stabilization.md).
//
// minGradPsiBand measurably does NOT lead the parasitic instability
// (|grad psi| ~ 0.9 at the kappa spikes). The mechanism observable is the
// amplitude of the grid-scale NORMAL-PROFILE modes of psi in the band -- the
// 2h compression wave the parasitic strain pumps (and, after the filter
// falsification, the 4h/8h content the loop can migrate to). Per band cell:
//
//   n_hat  from the linear LLS interface plane (halo-only, fit-free);
//   a +/-n_hat face-neighbour chain, picking per step the neighbour that
//   maximizes the advance n_hat . (x_N - x_cur);
//   A(kh)_c = |psi_{+k} - 2 psi_0 + psi_{-k}| / 4   for k = 1, 2, 4
//   (a pure 2h mode of amplitude a gives A(h)_c = a; the wider spacings
//   respond to 4h/8h content -- per-mode calibration constants are fixed
//   against the WP1 injector, phase-averaging happens in the band norm).
//
// Band norms A{2h,4h,8h}_{max,2} (volume-weighted L2). PARALLEL CONTRACT:
// A2h uses only face neighbours + the one-layer coupled halo
// (coupledFaceNeighbours pattern) => serial/np-identical to round-off.
// A4h/A8h chains are rank-local; a chain that cannot advance (processor
// boundary, physical boundary, or no aligned neighbour) EXCLUDES that cell
// from that mode's norms, and the excluded count is reported -- their norms
// are np-independent only up to the seam-adjacent exclusion set (documented;
// the carrier-identification use case is a large-signal question).
//
// Cost: one LLS plane + <= 8 neighbour visits per band cell per call --
// O(band) work, no fits, no solves.

#ifndef bandModeSpectrum_H
#define bandModeSpectrum_H

#include "levelSetPlaneReconstruction.H"

namespace Foam
{

struct bandModeSpectrumNorms
{
    // [0] = A2h (spacing h), [1] = A4h (spacing 2h), [2] = A8h (spacing 4h)
    scalar aMax[3] = {0, 0, 0};
    scalar aL2[3]  = {0, 0, 0};
    label nExcluded[3] = {0, 0, 0};   // cells excluded per mode (chain short)
};


inline bandModeSpectrumNorms computeBandModeSpectrum
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand
)
{
    bandModeSpectrumNorms out;
    static const label spacing[3] = {1, 2, 4};   // in cells along n_hat
    const label maxDepth = 4;

    const coupledFaceNeighbours coupledNei(mesh, psi);
    const auto& C = mesh.C();
    const auto& cellCells = mesh.cellCells();
    const scalarField& V = mesh.V();

    scalar sumV[3] = {0, 0, 0};
    scalar s2[3] = {0, 0, 0};

    // One-layer coupled halo for the A2h contract: neighbour psi/centre
    // across processor patches, indexed [patch][patchFace].
    const auto& pbm = mesh.boundaryMesh();

    forAll(narrowBand, c)
    {
        if (narrowBand[c] != 1.0) { continue; }

        const scalarList pc = leastSquaresPlaneCoeffs(mesh, psi, c, coupledNei);
        vector nHat(pc[0], pc[1], pc[2]);
        const scalar nMag = mag(nHat);
        if (nMag < SMALL) { continue; }   // degenerate plane: skip cell
        nHat /= nMag;

        // psi at chain stations -maxDepth..+maxDepth along n_hat.
        // station[maxDepth] = the cell itself.
        scalar station[2*4 + 1];
        bool haveStation[2*4 + 1] = {false};
        station[maxDepth] = psi[c];
        haveStation[maxDepth] = true;

        // Minimum advance per chain step: without it, a chain blocked along
        // n_hat (processor seam, misaligned neighbours) WANDERS through
        // near-tangential cells instead of truncating -- measured as
        // decomposition-dependent A4h/A8h with zero exclusions. Half the
        // local cell spacing keeps stations near +/- k h.
        scalar hC = GREAT;
        for (const label nb : cellCells[c])
        {
            hC = min(hC, mag(C[nb] - C[c]));
        }
        const scalar minAdv = 0.5*hC;

        for (const label dir : {+1, -1})
        {
            label cur = c;
            vector xCur = C[c];
            for (label k = 1; k <= maxDepth; ++k)
            {
                // Best advancing neighbour among rank-local cell neighbours.
                label best = -1;
                scalar bestAdv = minAdv;
                for (const label nb : cellCells[cur])
                {
                    const scalar adv = dir*(nHat & (C[nb] - xCur));
                    if (adv > bestAdv) { bestAdv = adv; best = nb; }
                }

                // k = 1 only: the one-layer coupled halo keeps A2h exactly
                // np-consistent when the best advance is across a seam.
                if (k == 1 && cur == c)
                {
                    const cell& cFaces = mesh.cells()[c];
                    forAll(cFaces, i)
                    {
                        const label faceI = cFaces[i];
                        if (mesh.isInternalFace(faceI)) { continue; }
                        const label patchI = pbm.whichPatch(faceI);
                        if (patchI < 0 || !coupledNei.valid(patchI)) { continue; }
                        const label pf = faceI - pbm[patchI].start();
                        const vector& xN = coupledNei.C(patchI)[pf];
                        const scalar adv = dir*(nHat & (xN - xCur));
                        if (adv > bestAdv)
                        {
                            bestAdv = adv;
                            best = -2;   // halo hit
                            station[maxDepth + dir*1] = coupledNei.psi(patchI)[pf];
                        }
                    }
                }

                if (best == -2)
                {
                    haveStation[maxDepth + dir*1] = true;
                    break;   // cannot walk further into the halo
                }
                if (best < 0) { break; }   // chain truncated

                cur = best;
                xCur = C[cur];
                station[maxDepth + dir*k] = psi[cur];
                haveStation[maxDepth + dir*k] = true;
            }
        }

        for (label m = 0; m < 3; ++m)
        {
            const label k = spacing[m];
            if (!haveStation[maxDepth + k] || !haveStation[maxDepth - k])
            {
                ++out.nExcluded[m];
                continue;
            }
            const scalar A =
                mag(station[maxDepth + k] - 2.0*station[maxDepth]
                    + station[maxDepth - k])/4.0;
            out.aMax[m] = max(out.aMax[m], A);
            sumV[m] += V[c];
            s2[m] += V[c]*A*A;
        }
    }

    for (label m = 0; m < 3; ++m)
    {
        reduce(out.aMax[m], maxOp<scalar>());
        reduce(sumV[m], sumOp<scalar>());
        reduce(s2[m], sumOp<scalar>());
        reduce(out.nExcluded[m], sumOp<label>());
        out.aL2[m] = (sumV[m] > 0) ? Foam::sqrt(s2[m]/sumV[m]) : 0;
    }

    return out;
}

} // End namespace Foam

#endif

// ************************************************************************* //
