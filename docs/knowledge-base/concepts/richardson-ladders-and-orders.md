---
title: "Richardson ladders and observed orders"
description: "Three rungs at least, matched horizon and time-step law, the observed order next to every error and the resolution range with every number; one rung called the distance-cone bound a gain and two rungs called it a floor; the 2D orders of METHOD 8.3.7 were 3/2 too high"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-28
date_settled: 2026-09-10
decided_by: [author decision 2026-09-10, config/gates/methodGate2D.yaml]
code: [workflow/scripts/richardson.py, workflow/scripts/make_gate_summary.py, workflow/scripts/advection_convergence_table.py, workflow/scripts/value_bound_ladder_table.py]
sources: [CLAUDE mesh convergence section, CLAUDE method gates section, METHOD 8.3.5, METHOD 8.3.7, STATUS 11.2, PHL 5.4]
---
# Richardson ladders and observed orders

> Verdict (2026-09-28). Any change that touches the advection of the level set is measured on a mesh convergence study of at least three resolutions, with a matched horizon and a matched time-step law, and no statement is made from a single resolution; the observed order is reported next to the error for every entry of the vector, and the resolution range with every number ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L510-L538), [[concepts/error-vector-and-read-out-instants]]). The rule is written from one incident: the distance-cone bound lowered every interface metric by 25 to 68 % at N = 64, and at N = 128 its eikonal error moved from 7.119e-03 to 7.066e-03 (order 0.01, a floor) while its centroid error reversed to 119 % worse than no bound ([METHOD 8.3.5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L692-L730), [[retractions/distance-cone-bound-as-transport-bound]]). A second incident corrected the ladders themselves: every 2D advection order in METHOD 8.3.7 was 3/2 of the true value, because the table script used `h_eff = nCells^(-1/3)` for a 2D case ([STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3243), [[retractions/advection-orders-3-2-factor]]). The method gates run Richardson ladders by construction ([[concepts/method-gates]]).

## What it is

Two kinds of quantity, two procedures ([`richardson.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/richardson.py#L1-L33)):

1. An error with the exact value zero (shape, band gradient, volume, spurious current, pressure jump, curvature): the pairwise orders `p_ij = ln(e_i/e_j) / ln(h_i/h_j)` between consecutive rungs and the least-squares slope of `ln e` against `ln h` over all rungs. The error itself is the measurement; there is no extrapolation ([same](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/richardson.py#L41-L68)).
2. A quantity without an exact value (the oscillation period, the damping rate): the procedure of Celik et al. ([DOI 10.1115/1.2960953](https://doi.org/10.1115/1.2960953)) with the safety factor `Fs = 1.25`: the apparent order `p` by fixed-point iteration, so a non-constant refinement ratio is handled exactly; the extrapolated value; the GCI of the fine rung; the convergence type from `R = eps21/eps32` (monotone for `0 < R < 1`, oscillatory for `-1 < R < 0`, divergent for `|R| > 1`); and the asymptotic-range ratio `GCI_32 / (r21^p GCI_21)`, about 1 inside the asymptotic range ([same](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/richardson.py#L70-L98)). The self-test recovers `p = 1, 2, 3` to 1e-8 on the four integer gate ladders, 64 checks ([STATUS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3281-L3306)).

The ladder of the gates: in 2D the cell count doubles per rung (h ratio 1.414); in 3D the h ratio is at least 1.3 per rung, about 2.2 times the cells, because doubling the cells (h ratio 1.26) puts the rungs too close for a stable order estimate; three rungs at least; the first rung has R/h >= 10; all rungs have the same parity of N, so the droplet centre sits on a vertex at every rung ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L621-L629)). The fit quality is reported with the order where a least-squares fit is used: `p = 0.88 (R = 0.999)` for Popinet's L2 error ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L145-L158)), and the 36-arm viscosity ladder chose the one model with a positive order in all three metrics at both viscosity ratios ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L730-L752)).

## Why it matters

A single mesh gives an error, not a result. An error that is lower at one resolution can be a floor, a crossover or a coarse-mesh artefact, and the order is the only column that separates them ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L510-L516)). Two rungs are not a ladder: the interior growth of the translating droplet shrinks with h between N = 100 and 142 in the 40 mm box, and no order is stated until the third rung runs ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3932-L3948)). Three rungs are the minimum, not the target: a fourth point once falsified a trend that three had made look like clean second order ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L463-L471); the study is not named there).

## Evidence

| claim | number | where |
|---|---|---|
| one rung said gain, two rungs said floor | at N = 64 every metric 25 to 68 % lower with the cone bound; at N = 128 the eikonal error 7.119e-03 to 7.066e-03 (order 0.01), the spurious-current order 0.24 to 0.03, the centroid error +119 %; the unbounded run converges at 1.10 and reaches the floor at about N = 256 | [METHOD 8.3.5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L692-L730), MEASURED |
| the converged advection ladders reverse a single-rung reading | vortex, cone against none: 8.6x at N = 32 to 189.7x at N = 256, order 2.09 against 0.79 at the finest pair; volume 3.95e-05 against 1.302e-02, order 2.50 against 0.61 | [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L657-L679), MEASURED |
| the unbounded translation saturates | orders 3.80, 2.10, -0.54; the error rises from 7.947e-05 at N = 128 to 1.159e-04 at N = 256 | [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L640-L655), MEASURED, open |
| the 2D orders were 3/2 too high | translation `none` 3.80, 2.10, -0.54 (published 5.70, 3.15, -0.82); vortex `none` 2.76, 3.23, 2.09 (published 4.15, 4.85, 3.14); `h_eff = nCells^(-1/3)` | [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3243), [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L625-L636), MEASURED |
| a ladder whose entries are not monotone in h | the long-box translating droplet at N = 100, 142, 200: shape 1.38e-4, 7.08e-5, 3.99e-5 (pairwise orders 1.90, 1.68); spurious current 5.14e-4, 1.80e-3, 7.36e-4 (-3.6, 2.6); curvature error 15.4, 21.6, 20.2 1/m (-1.0, 0.2) | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3800-L3815), MEASURED |
| the third rung was not run when two rungs decided | a 0.8 % change across a 2x refinement locates the floor; N = 256 would only locate the crossover | [METHOD 8.3.5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L726-L730), author decision |

## Decisions

- The three requirements of the mesh-convergence rule ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L527-L538)); the gate ladders ([[decisions/process-gates-2d-first-and-no-best-yaml]]).
- The standing regression set runs the default configuration at three resolutions on three rungs ([[concepts/advection-regression-set]]).

## Open questions

1. The curated `advConv2D*_convergence.csv` tables are regenerated on Lichtenberg with the corrected script (open in [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3243)).
2. The saturation of the unbounded translation at N = 256 needs its own investigation before that rung scores anything ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L649-L655)).
3. The third rung of the 40 mm translating box ([[cases/translating-droplet]]).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/method-gates]], [[concepts/error-vector-and-read-out-instants]], [[concepts/advection-regression-set]], [[concepts/bit-identity-and-inertness-gates]].
- Retractions: [[retractions/advection-orders-3-2-factor]], [[retractions/distance-cone-bound-as-transport-bound]].
- Advection: [[concepts/value-bounds-and-clips]], [[models/sl-value-bound]], [[cases/kinematic-advection-cases]].

## Log

### 2026-09-28
Created.
