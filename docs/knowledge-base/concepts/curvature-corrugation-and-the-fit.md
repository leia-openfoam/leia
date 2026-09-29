---
title: "Curvature corrugation, aliasing and the psi filter as instrument"
description: "How a grid-scale perturbation of the level set reaches the force through the quadratic fit, what the band modes and the foliation residual say about the growth, and why the band filter stays a research instrument: 5.86 times better at R/h = 15.8 and 1.61 times worse at R/h = 10 (2026-09-28)."
aliases: [corrugation, aliasing of the fit, psi filter, biharmonicBand, band renormalisation, foliation residual]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-08-20
decided_by: [config/filterOffAmplifier3D.yaml, CLAUDE no-filtering rule, METHOD 8.1 row PSI_FILTER]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/psiFilterEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/bandModeSpectrum.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/bandRenormalization.H, applications/test/leiaTestFoliationResidual/leiaTestFoliationResidual.C]
sources: [PCS 0, PCS 6, PCS 7, PCS 14.2, PCS 15, PCS 18, PCT 0, METHOD 9.3, STATUS 4, PSH 0c, PSH 0d, CLAUDE no filtering, SL negative deck 6/1]
---
# Curvature corrugation, aliasing and the psi filter as instrument

> Verdict (2026-09-28). A grid-scale mode of the transported level set reaches the force through the quadratic fit: on diagonal stencils the `h/sqrt(2)` sampling is incommensurate with a `2h` profile mode, the least-squares misfit leaks into the tangential Hessian blocks that the curvature formula does not cancel, and the fit returns curvatures between -883 and 2744 1/m against the exact 953 at amplitude 0.5; an offline refit of the measured stencils reproduces every solver spike, so the spikes are the fit's true response ([PCS 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L84-L99)). The band modes are not what drives the growth phase: `max|U|` leads the `2h` amplitude by about 11 ms and the grid-scale explosion is the endgame ([PCS 6](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L456-L464)); the growth is carried by the residual smooth variation of the delivered curvature ([PCS 7](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L522-L529)), and the unstable content is the low-`m` shape while `m = 6` stays damped ([PCS 18.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1543-L1552)). The biharmonic band filter measures the size of that defect and is not a fix: at matched kick it is 5.86 times better at R/h = 15.8 and 1.61 times worse at R/h = 10.0, where it turns damping into growth ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1007-L1014)). Every filtered result before f83a1ab carries a decomposition-dependent seam bug ([[retractions/psi-filter-seam-bug]]), and production runs with `PSI_FILTER none` ([[decisions/psi-filter-none]]).

## What it is

Three objects are distinct here.

1. The corrugation of the transported zero set. In an exact rigid translation the `psi = 0` contour develops an RMS deviation of `0.475 h` (maximum `0.985 h`) by t = 0.05, growing at a rate fixed per unit time, about 0.28 per cell of displacement, weakly dependent on CFL over 0.007 to 0.5 ([METHOD 9.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L767-L776)). The transport operator amplifies on every mesh, spectral radius 1.0044 on hexahedra and 1.0111 on pMesh polyhedra ([METHOD 9.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L767-L769), [[concepts/polyhedral-fit-amplification]]). The `m > 4` corrugation of the `pointValue` scheme is `0.209 h` at `16 h` displacement ([PCT](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L144-L149)).
2. The aliasing of that corrugation into the curvature, the mechanism of the first paragraph. The delivered error on a clean interface is the mesh-locked `m = 4` pole bias ([PCS 18.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1495-L1510), [[concepts/curvature-from-the-fit]]).
3. The loss of the parallel foliation. Under transport `beta = |grad psi|` spreads from [0.998, 1.000] to [0.724, 1.373] in the band. The foliation residual `D` grows 600 times over a run while the total curvature error grows 16 times; it supplies about 21 percent of the growth, and the degrading fit the rest ([PCS 14.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1164-L1192), [[concepts/cell-centre-inverse-curvature]]).

The instruments. The band mode spectrum `A2h, A4h, A8h` reads the normal-direction profile ([`bandModeSpectrum.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/bandModeSpectrum.H)). The filter `psi <- psi - theta L(L(psi))` in the band plus one ring, `L = I - avg(interp(.))`, annihilates linear fields and damps the `2h` carrier ([negative deck 6/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/6/1), [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L533-L544)). The band renormalisation `psi <- psi/beta_Gamma` restores the parallel foliation ([`bandRenormalization.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/bandRenormalization.H), [PCS 15](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1213-L1218)). `leiaTestFoliationResidual` evaluates `D` on saved fields with a wide stencil that would never go in the solver ([`leiaTestFoliationResidual.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestFoliationResidual/leiaTestFoliationResidual.C#L5-L46)).

## Why it matters

The reinitialisation-free transport has zero damping for perturbations of `psi` that leave the zero set in place ([PCS 18.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1512-L1521)), so a filter that damps them is the cleanest falsifier of the damping theory, and it worked as a delay: `theta = 0.2` moves `t_blow` from 0.105 to 0.442 s at N = 128 ([negative deck 6/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/6/1)). But it carries its own weak source, which dominates wherever there is nothing to remove, so `theta` would have to vanish with the corrugation content, that is, it needs retuning per resolution ([PSH 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L338-L349)). That is the definition of a tuning knob, and the repository rule is that a filter's benefit measures the defect's size ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L382-L398)). The gradient-control candidates of [[hubs/gradient-control]] attack the third object, the drift of `|grad psi|`, without a filter.

## Where in the code

- The filter: [`psiFilterEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/psiFilterEqn.H); tokens `PSI_FILTER none`, `PSI_FILTER_THETA 0.05` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L483-L485)). The seam fix: `L(psi)` inherited `calculated` patch types from `fvc::average`, and the band dilation looped internal faces only; both fixed in f83a1ab with `config/seamConsistency3D{serial,par4}.yaml` as the regression ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L894-L917)).
- The renormalisation with its gate, smoothing and relaxation keys ([`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L358-L376)).

## Evidence

| claim | number | where |
|---|---|---|
| the fit aliases a 2h profile mode into the curvature | kappa in [-883, 2744] against 953 at eps = 0.5 on diagonal stencils; +-10 percent at grid-aligned poles; offline refit matches to 0.00 percent | [PCS 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L84-L99), MEASURED |
| the band modes do not lead the growth | onsets at N = 128: minGradPsiBand 0.028 s, max abs U 0.0655, A2h 0.0763, kErrL2Band 0.0784, A4h 0.0791, FPE 0.0803 | [PCS 6](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L444-L464), MEASURED |
| no wavelength carries the filtered growth phase | N = 256 filtered: max abs U onset 0.0677 s; A2h 0.1270, A4h 0.1417, A8h 0.1536; FPE 0.1672 | [PCS 7](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L494-L521), MEASURED |
| the 3D spectrum is flat in the unstable arm | A2h / A4h / A8h at 1.00x of t = 0 while max abs U grows 70x; delivered variation 182x across, 143x along | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L820-L825), MEASURED, filtered runs |
| the filter is a delay, not a closure | theta 0.05 survives t = 0.3 with a tail at 33 against 190 1/s; theta 0.2 t_blow 0.442 against 0.105 s; N = 256 filtered blows at 0.167 s | [negative deck 6/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/6/1), [PCS 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L90-L97), MEASURED |
| the filter's damping flips sign with resolution | 5.86x better at R/h = 15.8 (G 3.52 to 1.76); 1.61x worse at R/h = 10.0 (G -0.21 to +0.27), kicks matched to 4 digits | [PSH 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L326-L349), MEASURED |
| the filter is not the source | filter off, 3D: A = +7.68e-4 and +1.16e-3 per step at R/h = 15.8 and 20 | [PSH 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L259-L264), MEASURED |
| the seam bug | filter on, before / after / off: max abs U 2.43e-4 / 1.72e-6 / 1.55e-6; 53 to 205x worse before | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L906-L917), MEASURED |
| the filtered 2D convergence claim was the bug | post-fix max abs U 2.24e-5 / 2.81e-5 / 8.83e-5 at N = 64 / 128 / 256, orders -0.32 and -1.65 | [PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L121-L143), MEASURED |
| the foliation residual takes over late | bias 0.35 to 212 1/m while the error goes 70 to 1100; 21 percent of the growth | [PCS 14.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1164-L1192), MEASURED |
| band renormalisation is worse | t_blow 0.0195 against 0.0668 s (3.4x); relaxation at omega 1 / 0.2 / 0.05 gives min abs grad psi 0.9716 / 0.9751 / 0.9800 against 0.9877 untouched | [PCS 15](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1220-L1258), MEASURED |

## Why it failed, or why we think so

Grid-scale dissipation buys constant factors, not a sign flip of the eigenvalue; the feedback cascades to the shortest undamped scale, and `min|grad psi|` erodes from 0.94 to 0.33 with no `2h` content ([negative deck 6/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/6/1)). The renormalisation fails for a quantitative reason: `beta_Gamma` comes from the same fit, whose error at N = 64 is 0.3 to 0.6 percent, while the drift to correct is about 0.1 percent, so the correction is dominated by the error of its own estimator at any relaxation ([PCS 15.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1252-L1258)).

## Decisions

- `PSI_FILTER none` in production; every candidate is scored with the filter off ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L385), [CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L397-L398)); see [[decisions/psi-filter-none]].
- Filtered results before f83a1ab are not method properties ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922)); see [[retractions/psi-filter-seam-bug]].

## Open questions

1. A diagnostic for the two tangential directions of a 3D interface; the normal-direction spectrum is blind to them ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L885-L890)).
2. Whether the corrugation growth rate `r(A2h)` of PCS 16.3 and the `m = 2` mode rate of PCS 18 describe one object; PCS 18 supersedes the finite-intercept reading of 16.1 ([PCS 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1571-L1577)).

## Related

[[hubs/surface-tension]], [[hubs/gradient-control]], [[concepts/curvature-from-the-fit]], [[concepts/cell-centre-inverse-curvature]], [[concepts/face-curvature-deliveries]], [[concepts/parasitic-current-mechanism]], [[concepts/polyhedral-fit-amplification]], [[concepts/normal-projected-sl]], [[concepts/redistancing-geometric-grl]], [[models/sl-reconstruction]], [[decisions/psi-filter-none]], [[retractions/psi-filter-seam-bug]], [[cases/stationary-droplet]], [[studies/curvature-stabilization-campaign]].

## Log

### 2026-09-28
Created from PCS sections 0, 6, 7, 14, 15 and 18, METHOD 9.3, STATUS 4, the Shannon plan sections 0c and 0d and the negative-results deck.
