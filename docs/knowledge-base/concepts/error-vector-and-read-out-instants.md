---
title: "The error vector and the read-out instants"
description: "Shape, volume, band gradient, spurious current and pressure jump are reported together, in L2 and L1 only, at the instant the flow allows; a headline metric alone has misread candidates three times"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-28
date_settled: 2026-09-03
decided_by: [author decision 2026-09-03, config/gates/methodGate2D.yaml]
code: [workflow/scripts/make_gate_summary.py, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H, workflow/scripts/make_translating_matrix_table.py]
sources: [CLAUDE research loop steps 5 and 6, STATUS 0, SL article sec:translating, PHL 5.4, METHOD 7, STATUS 11.15]
---
# The error vector and the read-out instants

> Verdict (2026-09-28). A result for an interface method is a vector, never one number: the shape (geometric) error, the volume error, the interface-profile diagnostic (`|grad psi|` in the band, or the boundedness of alpha), the spurious-current level and the pressure-jump error, reported together, with the commit and the binaries ([CLAUDE.md, step 6](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L502-L508)). Every entry is an L2 or an L1 norm. An L_inf norm is never reported as a result: on a parasitic current the maximum is a single-cell extremum set by the sub-cell interface position, it does not converge under refinement, and a verdict built on it was retracted on 2026-09-03 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L196-L200)). Each entry is read at the instant the flow allows: in a reversed flow the gradient at T/2, the shape at T, the volume at both, because the reversal cancels errors at the endpoint; two runs are compared at equal step counts ([CLAUDE.md, step 5](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L463-L471)). The method gates implement this vector ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L1-L27), [[concepts/method-gates]]).

## What it is

The gate's vector per arm ([PHL 5.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L539-L569), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L54-L65)):

| entry | kinematic arm | droplet arms | instant |
|---|---|---|---|
| shape | `E_GEOM_ALPHA_REL` (relative L1 of alpha against the initial field) | `zeroSetRadialL2 / R`, the L2 radial distance of the zero-set crossings to the exact circle | T |
| band gradient | `E_NARROW_L2_GRAD_PSI` | `gradPsiL2ErrorBand`, L2 of `|grad psi| - 1` with the unlimited least-squares gradient ([STATUS 11.3](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3261-L3270)) | T/2 in a reversed flow, T for the droplets |
| volume | `E_VOL_ALPHA_REL` | `phaseVolumeRelError` | T/2 and T |
| boundedness | `E_BOUND_ALPHA`, the largest violation of `0 <= alpha <= 1` | same | maximum over t; scored, no order |
| spurious current | none | `meanMagUPrime` (L1) and `l2MagUPrime` (L2) of `U - U0`, stationary and translating only | T |
| pressure jump | none | `abs(pLaplace - sigma (d-1)/R) / (sigma (d-1)/R)`, the jump over the pure phases | T |
| curvature | none | `kErrL2Band R/(d-1)` | T |
| travel | none | the centroid displacement over `U0 T`, minus 1 | T |
| oscillation | none | the period and the damping rate of the mode-2 amplitude, from the zero crossings of `m2CosCoefficient` and a fit of `ln|peak|`; the Lamb period as a check | T |

The pressure jump is taken over the pure phases, `alpha > 1 - eps` against `alpha < eps` with `eps = 1e-2`; partitioning at `alpha = 1/2` puts the whole transition into the drop average and biases the jump low by about 3 % (70.4456 Pa against the true discrete 72.7394 Pa), and the old metric was non-monotone in solver quality ([METHOD 7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L321-L339)). The absolute lengths (`zeroSetRadialL2`, `centroidError`, the `m2*` columns) are in metres and are divided by `R` before any threshold; in SI units they are 400 times smaller than in Popinet's units, and a fixed threshold read every SI run as clean ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2226-L2231)).

## Why it matters

Three misreadings that the vector prevents. A candidate improved the volume order while failing the gradient; a single-metric view would have called it a win ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L502-L508)). On the 1D stretch the sourceless arm has an exact zero contour while its distance property has collapsed by a factor `e`; a shape-only or volume-only score calls the failed arm perfect ([`sdpls1Dstretch`](https://github.com/leia-openfoam/leia/blob/8867581/config/sdpls1Dstretch.yaml#L26-L30)). And the first curation of the repaired translating matrix reported an 18x difference between the density ratios; it had read `L_inf` at a common step set by the first arm to die, inside every water/air blow-up. Read at step 5000 in `L1`, the factor is 1.4 to 1.6 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L119-L125)).

The empirical basis for the L_inf rule is Popinet's benchmark: the half-order maximum against a near-first-order RMS reproduces across two unrelated discretisations, so it is a property of the benchmark, not of the scheme ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L160-L163), [[cases/popinet-translating-droplet]]). On the translating droplet the same data gives `L_inf` orders of 0.84 and -0.19 against `L2` orders of 1.04 and 0.63 over N = 256 to 512, with `L_inf/L2` between 10 and 34 ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1920-L1930)).

## Evidence

| claim | number | where |
|---|---|---|
| L_inf does not converge where L2 does | orders 0.84 and -0.19 (L_inf) against 1.04 and 0.63 (L2), N = 256 to 512; L_inf/L2 10 to 34 | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1920-L1930), MEASURED |
| Popinet's own split reproduces | L2 order 0.88 (R = 0.999), L_inf 0.49 (0.999) against his "close to first order" and "less than first order" | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L145-L163), MEASURED |
| the instant changes the verdict | 18x (L_inf at step 8427) against 1.4 to 1.6x (L1 at step 5000) | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L119-L125), MEASURED |
| the pressure-jump partition | 70.4456 Pa (-3.15 %) at alpha = 1/2 against 72.7394 Pa (-0.00 %) over the pure phases, exact constant curvature | [METHOD 7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L321-L333), MEASURED |
| the gate scored an L_inf norm once | `maxMagU` was scored and given an order in the first smoke; corrected on 2026-09-26 | [STATUS 11.11](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3495-L3499), [PHL](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L554-L557), MEASURED |
| the oscillating velocity is physics, not error | a candidate that damped the oscillation more scored as better until 4ef98db; the norms are now `oscL2MagU`, `oscMeanMagU`, reported | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3867-L3878), MEASURED |
| a value at T on the edge of a divergence carries no verdict | the baseline oscillating arm at N = 200: band gradient error 5.13 at T (0.10 and 0.14 at the coarser rungs), period 8.06 ms against 10.0 and 9.82 ms | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3955-L3971), MEASURED |

## Decisions

- Never L_inf; L2 and L1 only (author decision 2026-09-03, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L196-L200)). Where the record only has L_inf (Popinet's table), the note says so.
- `rhoClipFraction` is reported, not scored ([[concepts/bound-rho]]).

## Open questions

1. A centred band-gradient metric is blind to a cell-scale mode; a second-difference norm or a one-sided `q` is the proposed repair ([[concepts/method-gates]]).
2. The oscillating arm needs a read-out over time, not only at T ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3955-L3971)).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/method-gates]], [[concepts/richardson-ladders-and-orders]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/bit-identity-and-inertness-gates]].
- Cases: [[cases/popinet-translating-droplet]], [[cases/exact-1d-stretch]], [[cases/oscillating-droplet]], [[cases/kinematic-advection-cases]].
- Mechanism: [[concepts/density-ratio-amplifier]].

## Log

### 2026-09-28
Created.
