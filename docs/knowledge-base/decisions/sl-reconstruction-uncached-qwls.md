---
title: "SL_RECONSTRUCTION: the uncached quadratic value fit"
description: "SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares in every droplet case layer, decided by the transport ladders: shape order 2.84 on the 2D vortex at CFL 0.5, 2.95 and 3.28 on the 3D shear case on hexahedra and polyhedra"
aliases: [SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-08-27
decided_by: [config/uncachedConv2Dvortex.yaml, config/advConv2Dvortex.yaml, config/advConv3DshearHex.yaml, config/advConv3DshearPoly.yaml]
code: [src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.H, src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C, cases/default.parameter, cases/stationaryDroplet2D.parameter]
sources: ["METHOD 8.1 row SL_RECONSTRUCTION (L373)", "METHOD 2.2 (L98-L125)", "METHOD 8 (L343-L347)", "PCS L52-L62", "DP L29-L32", "SL article sec:recon", "STATUS 11.19"]
---
# SL_RECONSTRUCTION: the uncached quadratic value fit

> `SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares` in the per-case layer of every two-phase droplet case since 3495a44 (2026-08-14) ([cases/stationaryDroplet2D.parameter L87](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L87)), and in the `lineTokens` of the 2D method gate for every semi-Lagrangian arm ([methodGate2D.yaml L44-L49](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L44-L49)). Decided by the transport ladders ([METHOD 8.1 L373](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373)): the shape error converges at order 2.84 on the 2D reversed vortex at CFL 0.5, re-established on 2026-08-27 with the gradU fix ([PCS L52-L61](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61)), and at 2.95 on hexahedra and 3.28 on cfMesh polyhedra on the 3D shear case ([METHOD 8 L343-L347](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347)). The cached member `quadraticWeightedLeastSquares` computes the same fit. The uncached member is the only member that exposes `footPointDistance` and `fitDerivatives` to the curvature path, and it holds no per-cell cache, so a 128^3 case fits in serial memory ([faceCurvatureSphere3D.yaml L15-L17](https://github.com/leia-openfoam/leia/blob/8867581/config/faceCurvatureSphere3D.yaml#L15-L17)).

## The question

Which reconstruction of `psi^n` at the departure foot is the production choice, and which layer carries it? The global default of the token is the cached member ([DP L29-L32](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L29-L32)). The two-phase cases need the uncached member, so each case sets it in its `values { }` block ([METHOD L15-L16](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L15-L16)). Two properties of the fit are load-bearing: the fit interpolates stencil values, and its degree is at least two ([METHOD 2.2 L118-L121](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L118-L121)). The fit is a constant-free quadratic with the weights `1/|d|`, solved by Cholesky on the normal system ([METHOD 2.2 L98-L116](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L98-L116)); the stencil is cell-point-cell on hexahedra and cell-face-cell on polyhedra ([METHOD 2.2 L123-L125](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L123-L125)).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| 2D vortex, hex, CFL 0.5 and 1.0 | shape order and volume order | 2.84 / 3.30 at CFL 0.5, 2.38 / 3.54 at CFL 1.0 (2026-08-27); before the gradU fix 2.97 / 2.98 and 2.59 / 2.94 | MEASURED, [PCS L52-L61](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61) |
| 3D shear, hex and poly | shape order | 2.95 / 3.28 | MEASURED, [METHOD 8 L343-L347](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347) |
| 3D deformation, hex and poly | shape order, filament-limited | 1.36 / 1.46 | MEASURED, [METHOD 8 L343-L347](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347) |
| against the best Eulerian line, 512^2, T = 8 | accuracy at equal wall clock | 20x more accurate at half the wall clock | MEASURED, [METHOD 8 L343-L347](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347) |
| polyhedral 3D shear, no value bound | completion | diverged at step 198 in HEAD and in the pre-change binary | MEASURED, [METHOD 8.3.8 L682-L693](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L682-L693) |
| polyhedral boundary-layer cells | admissibility test `quadraticPivotTol 0.3` | 7.1 % of the cells of a uniform cfMesh box demoted to the linear fit, none within 12 h of an interface; hexahedra bit-identical | MEASURED, [STATUS L1657-L1681](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1657-L1681), [DP L1075-L1081](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1075-L1081) |

Pre-registered read-out: the standing advection ladder config ([advConv2Dvortex.yaml L1-L40](https://github.com/leia-openfoam/leia/blob/8867581/config/advConv2Dvortex.yaml#L1-L40)). The `none` arm's order must match the recorded order within 0.3 at every rung pair, and `E_BOUND_ALPHA` must stay at round-off.

On 2026-08-27 the 3D shear, the 3D deformation and the polyhedral rows were still to be re-run after the gradU fix ([PCS L62](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L62)). The regression of 2026-09-10 ran the hexahedral and the polyhedral 3D shear bit-identical against the pre-change binary ([METHOD 8.3.8 L682-L693](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L682-L693)); it did not re-fit the 3D orders.

## What it does not cover

1. The linear members. The consistent-linear line converges at order about 1.1 and is a research line ([[concepts/linear-semi-lagrangian]]).
2. Stability. The reconstruct-and-evaluate operator amplifies on every mesh: `rho(B) = 1.00441` on production hexahedra ([METHOD 8.2 L409-L417](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L409-L417), [[concepts/polyhedral-fit-amplification]]).
3. Uniform translation. RETRACTED 2026-09-29: "the unbounded scheme saturates at N = 256, order -0.54" ([METHOD 8.3.7 L641-L653](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L641-L653)) was read from a reversed flow ([[retractions/reversed-2dtranslation]]). One way, the shape error converges at orders 2.90, 2.16, 1.86 over N = 32 to 256 at CFL 0.5 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697)). The METHOD 8.1 row of this token now says that its hex CFL axis came from the reversed `kinematicTranslation2D`; one way, the CFL-1 orders are irregular, 0.88 and 3.49 ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L405)).
4. The fit solver. `SL_FIT` is a separate decision ([[decisions/sl-fit-normal-equations]]).
5. The curated `advConv2D*_convergence.csv` orders of METHOD 8.3.7 were 3/2 too high ([STATUS L3237-L3242](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3237-L3242), [[retractions/advection-orders-3-2-factor]]). The orders in this note come from the SL article ladders, not from those files.

## Related

[[hubs/advection]] - [[models/sl-reconstruction]] - [[models/sl-scheme]] - [[models/sl-value-bound]] - [[concepts/departure-foot-ab2-centring]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/advection-regression-set]] - [[concepts/linear-semi-lagrangian]] - [[decisions/sl-fit-normal-equations]] - [[decisions/sl-clip-and-value-bound-off]] - [[retractions/gradu-coupled-patch-contamination]] - [[retractions/advection-orders-3-2-factor]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the measurement of 2026-08-27; the per-case layer carries the value since 2026-08-14. Entered in [[decision-log#2026-08]].

### 2026-09-29
What it does not cover, item 3, CORRECTED: the translation saturation came from a reversed flow (STATUS 11.19, [[retractions/reversed-2dtranslation]]); the one-way orders and the CORRECTED METHOD 8.1 row are added. The decision does not change: it rests on the vortex and 3D shear ladders.
