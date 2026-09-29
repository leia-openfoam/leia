---
title: "CURVATURE_EXTENSION: cellCentreInverse, case-dependent"
description: "CURVATURE_EXTENSION cellCentreInverse in the global default since c935883 (2026-09-01), decided by config/stationaryLadder2Dshared.yaml: the unabsorbed capillary residual is 4.60 / 3.81 / 1.55 / 1.49x lower than with none at N = 32 / 64 / 128 / 256; the translating droplet is undecided and the Popinet family runs none"
aliases: [CURVATURE_EXTENSION cellCentreInverse]
kind: decision
status: settled
part: surface-tension
tags: [decision, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: [config/stationaryLadder2Dshared.yaml, "author decision 2026-09-01"]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H, cases/default.parameter, config/gates/methodGate2D.yaml]
sources: ["METHOD 4.1 (L161-L167, CORRECTED)", "METHOD 8.1 row CURVATURE_EXTENSION (L391, CORRECTED 2026-09-27)", "STATUS 0 (L242-L244)", "STATUS 11.13 (L3573-L3574, L3609)", "STATUS 11.14 (L3750-L3760)", "DP L360-L400", "PHL L695-L704"]
---
# CURVATURE_EXTENSION: cellCentreInverse, case-dependent

> `CURVATURE_EXTENSION cellCentreInverse` in the global default since c935883 (2026-09-01), set on the author's instruction to run one shared configuration ([DP L360-L386](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L360-L386), [METHOD 4.1 L161-L167](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167)). Decided for the stationary droplet by `config/stationaryLadder2Dshared.yaml`: the unabsorbed capillary residual over the second half is 4.60 / 3.81 / 1.55 / 1.49x lower than with `none` at N = 32 / 64 / 128 / 256 ([STATUS L242-L244](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L242-L244), [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)). The value is case-dependent, and the gate sets it per arm ([methodGate2D.yaml L57-L63](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L57-L63)): for the translating droplet the choice is undecided, because the only `none`-against-`cellCentreInverse` comparison ran on the closed box and is void (corrected 2026-09-27, [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)); the Popinet family runs `none` at density ratio 1 by reproduction, not by measurement; no `cases/<case>.parameter` sets the token ([STATUS L3573-L3574](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3573-L3574)).

## The question

The balanced-force assembly `kappa_f snGrad(alpha)_f` is a discrete gradient only when `kappa` is constant over the CSF support. Without extension each cell centre carries the parallel-surface value at its own offset `d`, `kappa(d) = kappa0 / (1 + d kappa0)`, so the band varies in the interface-normal direction and the force is not absorbable ([DP L373-L380](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L373-L380), [[concepts/balanced-force-csf-flux]]). `cellCentreInverse` applies the parallel-surface inverse in place in every cell, with the offset from the stable quadratic root along the normal ray and, in 3D, the Gaussian curvature ([METHOD 4.1 L161-L167](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167), [[concepts/cell-centre-inverse-curvature]], [[decisions/curvature-inverse-gaussian]]). Its non-gradient content converges at order +2.01 against +0.09 for every other cell curvature model and is 3200x smaller at N = 512 ([DP L380-L385](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L380-L385)). The question of 2026-09-01 was whether that theory shows in the coupled residual on the case with no mean flow.

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| stationary 2D, `none` (reproduction control), N = 32 / 64 / 128 / 256 | relative L2 capillary residual, projectedFlux | must reproduce 2.37e-5 / 3.15e-6 / 8.53e-7 / 6.43e-7 | MEASURED, [config L25-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryLadder2Dshared.yaml#L25-L30), [DP L1036-L1040](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1036-L1040) |
| stationary 2D, `cellCentreInverse` against `none` | ratio of the residual over the second half | 4.60x / 3.81x / 1.55x / 1.49x lower at N = 32 / 64 / 128 / 256 | MEASURED, [STATUS L242-L244](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L242-L244) |
| the non-gradient force content `alpha_f snGrad(kappa_c)` over the CSF support | convergence order | +2.01 against +0.09 for every other cell curvature model; 3200x smaller at N = 512 | MEASURED, [DP L380-L385](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L380-L385) |
| translating 2D, `bestConfigTranslating2D` | — | VOID (closed box) | [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391) |
| translating 2D after the fix, N = 100, np 4, one change per run | divergence time | reference 0.0868 s; `none` 0.0695 s (step 6404), 20 % earlier; one resolution, an indicator only | MEASURED, [STATUS L3750-L3760](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3750-L3760) |
| translating 2D, `translatingRepaired2D` (`none`, N = 128, ratio 838.8) | completion | diverged in all 8 arms at t = 0.063 to 0.075 s; the gate baseline with `cellCentreInverse` diverged at 0.077 to 0.094 s of 0.1 s | MEASURED, different setups, [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391) |
| oscillating 2D, `oscillatingLadder2Dshared`, N = 128 | completion; volume error at N = 32 and 64 | `cellCentreInverse` completed where `none` failed at 0.0982 s; `none` had the lower volume error at the two coarse rungs; the case used an algebraic psi | MEASURED, void decision open, [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391) |

Pre-registered read-out ([config L25-L40](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryLadder2Dshared.yaml#L25-L40)): the `none` arms must reproduce the recorded ladder; PASS = a lower residual at every rung and no rung lost. The ladder ran with `MASS_FLUX geometricFaceDensity`, the default of that day ([config L52-L62](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryLadder2Dshared.yaml#L52-L62)).

## What it does not cover

1. The translating droplet, undecided: a matched `none`-against-`cellCentreInverse` ladder on the repaired case is the missing measurement ([STATUS L3620-L3622](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3620-L3622)). The gate runs `cellCentreInverse` to 0.05 s, before both onsets; a proposal to record `none` for the translating arm is open ([PHL L695-L704](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L695-L704)).
2. The oscillating droplet, one rung and an algebraic psi ([[cases/oscillating-droplet]]).
3. `cellCentreInverse` was never scored on the varying-curvature ellipse gate; its second order is measured on constant curvature only ([METHOD 4.3 L240-L241](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L240-L241), [[cases/curvature-static-gates]]).
4. The face deliveries (`stabilizedFootPointFace` and its variants) are closed as a lever; they change the gain and lose an order ([[concepts/face-curvature-deliveries]], [[retractions/cell-mean-delivery-adoption]]).
5. The Popinet `none` is a reproduction setting at density ratio 1, in each Popinet config ([popinet3D_La12000_poly_dump4_qr.yaml L45](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet3D_La12000_poly_dump4_qr.yaml#L45)).

## Related

[[hubs/surface-tension]] - [[models/curvature-extension]] - [[concepts/cell-centre-inverse-curvature]] - [[concepts/curvature-from-the-fit]] - [[concepts/balanced-force-csf-flux]] - [[concepts/face-curvature-deliveries]] - [[decisions/curvature-inverse-gaussian]] - [[decisions/surface-tension-reconstructed-curvature]] - [[decisions/process-gates-2d-first-and-no-best-yaml]] - [[retractions/closed-box-translating-droplet]] - [[retractions/cell-mean-delivery-adoption]] - [[cases/stationary-droplet]] - [[cases/translating-droplet]] - [[cases/oscillating-droplet]] - [[cases/popinet-translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED for the stationary droplet on 2026-09-01; case-dependent; the translating rationale was voided on 2026-09-27. Entered in [[decision-log#2026-09]].
