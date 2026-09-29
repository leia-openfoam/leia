---
title: "PHASE_INDICATOR: the Detrixhe-Aslam tetrahedral fill"
description: "PHASE_INDICATOR detrixheAslam in the global default since c935883 (2026-09-01): second order on a static circle and equal to the geometric clip to eight digits, tolerance-free; the switch from geometric is not inert"
aliases: [PHASE_INDICATOR detrixheAslam]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: [config/phaseIndicatorConvergence.yaml, "author decision 2026-09-01"]
code: [src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.H, src/leiaLevelSet/phaseIndicator/levelSetPlaneReconstruction.C, cases/default.parameter]
sources: ["METHOD 3 (L129-L140)", "METHOD 8.1 row PHASE_INDICATOR (L388)", "DP L5-L9", "SL article sec:indicator (L622-L715) and sec:indicator-verif (L1702-L1716)", "PCS 11.3 (L863-L878)", "RM L170-L177 and L1339-L1352"]
---
# PHASE_INDICATOR: the Detrixhe-Aslam tetrahedral fill

> `PHASE_INDICATOR detrixheAslam` in the global default since c935883 (2026-09-01) ([DP L5-L9](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L5-L9), [METHOD 8.1 L388](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L388)). The indicator fits a local signed-distance plane by linear least squares over the cell and its neighbours, fans the cell into tetrahedra and fills each with the closed form of Detrixhe and Aslam; it is second order and tolerance-free ([METHOD 3 L129-L140](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L129-L140)). On a static circle the indicator volume converges at order about 2.0, and the geometric clip agrees with it to eight digits because both act on the same plane ([SL article L1702-L1716](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1702-L1716)). The change of the default from `geometric` is not inert; it fixed what every serious study had pinned by hand, after a convergence order once came within one run of a fit across two indicators ([DP L5-L8](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L5-L8)).

## The question

The shape and volume metrics, the capillary force support `snGrad(alpha)` and the face density all read `alpha`. A pointwise Heaviside of psi is first order at the interface and assumes `abs(grad psi) = 1`, which advection does not preserve ([SL article L622-L633](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L622-L633)). The plane-based members drop that assumption. Between the two plane-based members, `geometric` clips the cell against the plane with geometric tolerances; `detrixheAslam` integrates the plane per tetrahedron analytically and is robust in one-cell-thick 2D ([phaseIndicatorConvergence.yaml L3-L9](https://github.com/leia-openfoam/leia/blob/8867581/config/phaseIndicatorConvergence.yaml#L3-L9), [[models/phase-indicator]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| static circle r = 0.15, no advection, N = 32 to 256 | relative indicator volume error | order about 2.0; `geometric` and `detrixheAslam` agree to eight digits | MEASURED, [SL article L1702-L1716](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1702-L1716) |
| 2D vortex, N = 32 / 64 / 128, T = 0.5 / 2 / 4 / 8, both indicators | shape and volume errors | the two coincide to about 1e-11 | MEASURED, [phaseIndicatorConvergence.yaml L26-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/phaseIndicatorConvergence.yaml#L26-L30) |
| N = 512 on the 0.01 m box, before the guard fix | the plane fit | every band cell fell back to the flat plane; the curvature row read 10.2 against the honest 17.4 | MEASURED, [RM L1339-L1352](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1339-L1352) |
| N = 32 translating droplet, no disturbance velocity | Detrixhe-Aslam phase volume; centroid | grows about 6 %; lags about 0.26 h | MEASURED, one resolution, [RM L170-L177](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L170-L177) |

Pre-registered read-out: the study isolates the Heaviside computation, because the psi advection is identical for both members and only `alpha` differs ([phaseIndicatorConvergence.yaml L6-L9](https://github.com/leia-openfoam/leia/blob/8867581/config/phaseIndicatorConvergence.yaml#L6-L9)).

The default change of 2026-09-01 was the author's instruction to run one shared configuration ([DP L360-L371](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L360-L371)); the kinematic cases sweep both members in their case layers ([3Dshear.parameter L3](https://github.com/leia-openfoam/leia/blob/8867581/cases/3Dshear.parameter#L3)).

## What it does not cover

1. The signed-distance assumption survives in the face area fraction: its offset `d = psi / abs(grad psi)` carries the error `e_d = -(1/2) b_n d^2`; the Hessian-corrected root exists as `offsetDistance` and is not measured ([PCS L863-L878](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L863-L878)).
2. The determinant guard of the plane fit now runs on a centred, variance-normalised matrix; before the fix the N = 512 rows of two curated tables were biased ([RM L1339-L1352](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1339-L1352)).
3. The oscillating droplet initialised an algebraic psi (`implicitEllipsoid`) until 2026-09-26, so its gradient columns had no meaning; the past studies are not yet voided ([STATUS L3271-L3279](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3271-L3279), [[cases/oscillating-droplet]]).
4. `heaviside` and `sharpJump` have no ladder on record ([[models/phase-indicator]]).

## Related

[[hubs/advection]] - [[models/phase-indicator]] - [[models/narrow-band]] - [[models/mass-flux]] - [[concepts/balanced-force-csf-flux]] - [[concepts/error-vector-and-read-out-instants]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[cases/oscillating-droplet]] - [[cases/kinematic-advection-cases]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED; the default changed on 2026-09-01. Entered in [[decision-log#2026-09]].
