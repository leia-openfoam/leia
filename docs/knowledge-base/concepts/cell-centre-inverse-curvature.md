---
title: "cellCentreInverse with the K-aware inverse"
description: "The parallel-surface inverse applied per cell to the fit curvature; its non-gradient force content converges at order +2.01 against +0.09, it lowers the unabsorbed residual 4.60 to 1.49 times on the 2D stationary ladder, the Gaussian term is load-bearing in 3D, and its second order is measured on constant curvature only (2026-09-28)."
aliases: [cellCentreInverse, parallel-surface inverse, K-aware inverse, Gaussian-curvature-aware inverse]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: [METHOD 8.1 row CURVATURE_EXTENSION, config/oscillatingLadder2Dshared.yaml, config/faceCurvatureSphere3D.yaml, config/gates/methodGate2D.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H, applications/test/leiaTestParallelSurfaceInverse/leiaTestParallelSurfaceInverse.C]
sources: [METHOD 4.1 CORRECTED, METHOD 4.3 CORRECTED, METHOD 8.1 row CURVATURE_EXTENSION, STATUS 0, STATUS 4 K exonerated, RM K-aware inverse, PCS 11, SL deck 4/20, DP CURVATURE_INVERSE_GAUSSIAN]
---
# cellCentreInverse with the K-aware inverse

> Verdict (2026-09-28). Since c935883 (2026-09-01) the global default `curvatureExtension cellCentreInverse` applies the parallel-surface inverse `kappa^Gamma = (kappa - 2 K d)/(1 - d kappa + K d^2)` in every filled cell, with `kappa` and `K` from the cell's own fit and `d` the signed offset to the fit's zero set, and hands the corrected cell field to the arithmetic face interpolation ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167), [`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L1-L13)). Its basis is static and stationary: the non-gradient force content `alpha_f snGrad(kappa_c)` converges at order +2.01 against +0.09 for every other cell curvature and is 3200 times smaller at N = 512 ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L36-L48)); on the 2D stationary ladder the unabsorbed capillary residual is 4.60 / 3.81 / 1.55 / 1.49 times lower than with `none` at N = 32 / 64 / 128 / 256 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L242-L244)). The Gaussian term is load-bearing: the 3D sphere gate converges at `h^1.95` with `K` and at `h^1.02` without it ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1378-L1388)). The second order is measured on constant curvature only ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L54-L60)); the translating droplet is undecided and the oscillating evidence is one resolution on an algebraic level set ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)).

## What it is

A surface offset by the signed normal distance `d` has principal curvatures `k_i/(1 + d k_i)`. With `kappa = k_1 + k_2` and `K = k_1 k_2` the contour through the cell centre carries `kappa_d = (kappa + 2 K d)/Q`, `K_d = K/Q`, `Q = 1 + kappa d + K d^2`; eliminating `Q` gives the boxed inverse above. It is exact on a sphere for both signs of `d`, reduces to `kappa_d/(1 - d kappa_d)` in 2D because a pseudo-2D fit gives `K = +0` exactly, and is skipped when `|1 - kappa_d d + K_d d^2| <= 1/2` ([SL deck 4/20](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/20), [METHOD 4.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L222-L239)). Three design choices are recorded in the header:

1. Per cell, not per face. The cut-cell face delivery computes the same inverse and then assigns one value to every active face of the cell; that lumping is first order wherever the curvature varies along the interface. Here the corrected cell field keeps its cell-to-cell smoothness ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L15-L20), [[concepts/face-curvature-deliveries]]).
2. `signedOffset`, not `footPointDistance`. The foot-point search refuses cells whose zero set lies outside the trusted stencil (261232 of 262144 cells at N = 512), which would leave a discontinuous edge. The algebraic root `d_c = 2 psi_c/(|g| + sqrt(D))`, `D = |g|^2 - 2 psi_c h_nn`, always returns ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L22-L34)).
3. Never invert a cell the fill skipped. An unfilled cell holds exactly zero; the inverse of zero is `-2 K d/(1 + K d^2)`, zero in 2D but a curvature made from nothing in 3D. 93.5 percent of the cells (201888 of 216000 at 3D `N_L = 60`) were inverted that way before the guard; the coupled metrics moved by about 1e-7 because `snGrad(alpha)` vanishes there ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L82-L110), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L951-L959)).

The inverse does not assume `|grad psi| = 1`. Its weakest sufficient hypothesis is that the level sets in the band are the parallel offsets of the interface, equivalently `beta = |grad psi|` constant on each level set ([PCS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L771-L796)). Where that fails, the residual is `d D + O(d^2)` with `D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2`: first order in the offset, and larger than the uncorrected error where `|D| > |kappa^2 - 2K - D|` ([PCS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L806-L839)).

## Why it matters

The balanced force absorbs only the gradient part of the flux, and the driver is `alpha_f snGrad(kappa_c)` ([[concepts/balanced-force-csf-flux]]). Without the inverse each cell carries the parallel-surface value at its own offset, so the band varies in the interface-normal direction and the force is not a discrete gradient ([`config/oscillatingLadder2Dshared.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/oscillatingLadder2Dshared.yaml#L3-L12)). The inverse is the only cell field that makes the unabsorbable part converge ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L380-L384)). It lowers the floor of the parasitic current, not its growth: "the floor improved, the growth survived" ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L529-L535)), and a lever that improves the delivered force by two orders but leaves the per-step gain alone buys no stability ([PSH 5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L1015-L1016), [[concepts/parasitic-current-mechanism]]).

## Where in the code

- `applyCellCentreInverseCurvature` ([header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L69-L147)), called after the symbolic fill in [`slAlphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L270-L276) and at t = 0 in [`leiaSemiLagrangianLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/leiaSemiLagrangianLevelSetTwoPhaseFoam.C#L185-L191).
- The dispatch and the `gaussianCurvature` sub-key: [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L294-L332); tokens `CURVATURE_EXTENSION` and `CURVATURE_INVERSE_GAUSSIAN` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L292-L301)).
- `parallelSurfaceInverse` and `fitGaussianCurvature` live in [`stabilizedFootPointFaceCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H), shared with the face deliveries.
- The exact unit gate `leiaTestParallelSurfaceInverse` checks the sphere at both offset signs, a torus at `K < 0`, the guard branch and the `K` expression, and exits non-zero on a mismatch ([`leiaTestParallelSurfaceInverse.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestParallelSurfaceInverse/leiaTestParallelSurfaceInverse.C#L13-L33)).
- The Eulerian two-phase solver has no dispatch and cannot run it ([STATUS 11.16](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3704-L3707), [[concepts/eulerian-solver-mass-flux-port]]).

## Evidence

| claim | number | where |
|---|---|---|
| the non-gradient content of the force converges | 9647 / 2273 / 539.7 / 139.3 at N = 64 to 512, order +2.01 (signedOffset), +2.02 (stabilised foot); every other cell curvature 5.1e5 to 4.5e5 at +0.09; geometric face alpha 149.1 at +2.016 | [header](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L36-L52), MEASURED (`leiaTestRemainderTerm.csv`; the app is not in the tree at 8867581) |
| the 2D stationary ladder residual | 4.60 / 3.81 / 1.55 / 1.49 times lower than `none` at N = 32 / 64 / 128 / 256; the `none` values are 2.37e-5 / 3.15e-6 / 8.53e-7 / 6.43e-7 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L242-L244), [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1029-L1032), MEASURED |
| the first coupled result, unfiltered, N = 64 | max abs U 8.05e-6, volume 7.2e-6, shape 7.2e-7 at t = 0.1 | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L524-L531), MEASURED |
| the t_blow rows 0.100 / 0.078 / 0.036 s | do not reproduce on the current code | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1057-L1068), RETRACTED, see [[retractions/t-blow-baseline]] |
| the 3D sphere gate | h^1.95, 5.25e-3 1/m at N = 128 with K; h^1.02, 0.421 without K (equal to the raw delivery) | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1378-L1388), MEASURED |
| K off in the coupled 3D ladder | delivered content 962.6 / 964.1 / 959.2 at order +0.01 against 3.184 / 1.942 / 1.252 at +2.03; the R/h = 12.7 and 15.8 arms die at 0.0955 and 0.0772 s | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L853-L871), MEASURED on filtered runs that predate the seam fix ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922), [[retractions/psi-filter-seam-bug]]) |
| why K off is zeroth order | on a sphere the inverse gives 2/(R - d) with K and 2/(R - 2d) without; relative error d/R, content O(1) | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L873-L882), DERIVED |
| the constant is not reparametrisation-robust | 22.795 / 22.770 / 22.580 against 22.727 for phi, phi + 2 phi^2, phi + 20 phi^2 | [PCS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L798-L804), MEASURED |
| the inverse can be worse than no correction | 3.01 times worse at the major vertex of the 2:1 quadratic-form ellipse at d = 1e-3 | [PCS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L832-L839), DERIVED |
| the translating droplet after the inlet/outlet fix | `none` diverges 20 percent earlier (0.0695 against 0.0868 s), N = 100, np 4 | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3732-L3738), MEASURED, one resolution |
| the oscillating droplet | completed N = 128 where `none` failed at 0.0982 s; `none` had the lower volume error at N = 32 and 64; algebraic psi | [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), MEASURED, see [[cases/oscillating-droplet]] |

## Why it failed, or why we think so

It did not fail; it did what a delivery can do. The parallel-foliation hypothesis holds at t = 0 and stops holding as the run proceeds: the foliation residual grows 600 times over a run and supplies about 21 percent of the curvature-error growth, the degrading fit the rest ([PCS 14.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1178-L1192)). Every correction built from the same fit injects that fit's error ([PCS 15.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1271-L1277)). A cubic fit that could supply the third derivatives of `D` is rank-deficient on the current stencil ([PCS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L859-L861)).

## Decisions

- `CURVATURE_EXTENSION cellCentreInverse` as the global default; the gates set it on every droplet arm ([`config/gates/methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L112-L163)); see [[decisions/curvature-extension-cell-centre-inverse]].
- `CURVATURE_INVERSE_GAUSSIAN yes` ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L399)); see [[decisions/curvature-inverse-gaussian]].

## Open questions

1. Never scored on the varying-curvature ellipse gate ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), [[cases/curvature-static-gates]]).
2. A torus gate for the 3D ranking, exact signed distance with non-constant mean curvature ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2669-L2675)).
3. `none` against `cellCentreInverse` on a matched translating setup after 440107f ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3620-L3622), [[retractions/closed-box-translating-droplet]]).
4. The K-off coupled verdict rests on filtered runs before the seam fix; the static sphere gate stands ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922)).

## Related

[[hubs/surface-tension]], [[models/curvature-extension]], [[concepts/curvature-from-the-fit]], [[concepts/face-curvature-deliveries]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[concepts/curvature-corrugation-and-the-fit]], [[concepts/eulerian-solver-mass-flux-port]], [[cases/stationary-droplet]], [[cases/oscillating-droplet]], [[cases/curvature-static-gates]], [[decisions/curvature-extension-cell-centre-inverse]], [[decisions/curvature-inverse-gaussian]], [[retractions/psi-filter-seam-bug]], [[retractions/t-blow-baseline]], [[retractions/closed-box-translating-droplet]].

## Log

### 2026-09-28
Created from the header of `cellCentreInverseCurvature.H`, METHOD 4 and 8.1, STATUS 0 and 4, the roadmap's sphere gate and PCS 11.
