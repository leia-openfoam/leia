---
title: "The well-balanced gate with exact curvature"
description: "Prescribing the exact constant curvature isolates the force balance from the estimator: the Laplace jump is 145.470 Pa against 145.48 on refined and uniform meshes, the spurious velocity stays at the solver's round-off floor, and the exact-curvature arms of the amplifier gate stay bounded at every density ratio and translation speed (2026-09-28)."
aliases: [well-balanced gate, exact-curvature gate, constant-curvature gate, prescribed curvature]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-04
decided_by: [config/amplifierGate2D.yaml, config/kickOriginGate2D.yaml, config/transISTConstantCurvatureGate.yaml, workflow/scripts/run_pressure_compatibility_gate.py]
code: [src/leiaLevelSet/surfaceTensionForce/constantCurvatureSurfaceTension.H, src/leiaLevelSet/surfaceTensionForce/constantCurvaturePressurePotential.H, cases/stationaryDroplet2D/system/fvSolution.template]
sources: [CLAUDE research loop step 3, STATUS 0, STATUS 4 hanging-node gate, SL article well-balanced gate, SL article tab:amplifier, RM frozen circle, RM translating gate, PSH 0g]
---
# The well-balanced gate with exact curvature

> Verdict (2026-09-28). The gate replaces the reconstructed curvature by the exact constant `kappa = 1/R` (2D) or `2/R` (3D) and keeps everything else. It is the rung of the multiphase ladder between kinematic transport and the coupled 2D run, and it separates a force-balance bug from an estimator error; prescribed exact curvature gave identically zero velocity in all six arms ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L430-L453)). Measured: the one-step velocity on a frozen circle is `3.770e-9` to `1.039e-8` m/s at N = 32 to 128 on uniform meshes, and the product and potential forms agree to `1.5e-13` ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1046-L1057)); the exact-curvature arms of the amplifier gate stay bounded over 8000 steps at both density ratios and both translation speeds, at a plateau of `3.00e-5` of `U0` at ratio 838.8 and at round-off at ratio 1 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L200-L215)); and on the refined 3D meshes the Laplace jump is 145.470 Pa in every arm against the exact 145.48, with the spurious velocity at the solver's round-off floor, `3.8e-10` to `8.2e-10` m/s over time ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1207-L1219), [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1379-L1400)). The gate is not a method: a constant curvature cannot relax a deformed interface ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L266-L268)).

## What it is

Two oracle models exist. `constantCurvatureSurfaceTension` returns `sigma kappa snGrad(alpha) |S_f|`, the CSF product form with a constant `kappa`; `constantCurvaturePressurePotential` returns `snGrad(sigma kappa alpha) |S_f|`, the potential form. Both enter the unchanged `rAUf`-weighted pressure equation and the same `fvc::reconstruct` ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1030-L1038), [[models/surface-tension-force]]). The token `CONST_CURVATURE` sets the value in `fvSolution`; every other force model leaves the entry unread ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L248-L269)). Because the force is then an exact discrete gradient, any velocity that remains is delivery, density-interpolation or projection error, and not curvature error ([`config/kickOriginGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/kickOriginGate2D.yaml#L25-L35)).

The gate has been run in five forms:

1. One capillary step on a frozen analytic circle, uniform and 10 percent perturbed meshes, N = 32 to 128 ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1040-L1057)).
2. The translating droplet at N = 32 with `sigma = 0`, exact `kappa`, and the integral model, to t = 0.05 ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L157-L188)).
3. The step-1 kick of the repaired translating droplet, reconstructed against exact `kappa`, crossed with four translation speeds ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L76-L102)).
4. The full-horizon amplifier gate, 2x2x2 over the force model, `U0` and the density ratio ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L228)).
5. The hanging-node gate: `kappa = 2/R = 2000` on the refined N_fine = 60 meshes and on the uniform mesh, 921 steps to t = 0.01, np 4 and np 1 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1197-L1226), [[concepts/static-local-refinement]]).

## Why it matters

A coupled failure with the reconstructed curvature is ambiguous between the estimator and the assembly. With the exact curvature the assembly is tested alone, and the residual that remains under exact curvature is pure transport error: at t = 0.06 s the volume error (4.32e-3), the shape error (6.9e-6) and the travelled fraction (0.9975) are identical to three digits at both density ratios, and the net propulsion of the reconstructed arm (travel 1.0257) is gone ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L221-L224)). That is what decided the campaign: no amplifier independent of the source needs attacking, and the work is the curvature ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L217-L219), [[concepts/parasitic-current-mechanism]]). The gate is also the instrument for the seed-starvation question of the SDPLS sources ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L260-L264), [[studies/sdpls-pre-print]]).

## Where in the code

- [`constantCurvatureSurfaceTension.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/constantCurvatureSurfaceTension.H) and [`constantCurvaturePressurePotential.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/constantCurvaturePressurePotential.H); the frozen-circle runner `workflow/scripts/run_pressure_compatibility_gate.py` ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1040-L1044)).
- The refined-mesh gate arms G0 / G1 / G2 / GC of the static-refinement workflow ([`workflow/README.md`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L461)); the decomposition pair `stationaryDroplet3D{refined,uniform}Gate{1,4}` ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1414-L1419)).

## Evidence

| claim | number | where |
|---|---|---|
| a constant curvature on a frozen circle | 3.770e-9 / 7.635e-9 / 1.039e-8 m/s uniform; 8.584e-6 / 7.328e-4 / 1.392e-3 on 10 percent perturbed meshes, N = 32 / 64 / 128 | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1046-L1057), MEASURED |
| the product and potential forms coincide | max direct field difference 2.02e-14 to 1.50e-13 | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1046-L1053), MEASURED |
| the translating droplet at N = 32 | sigma = 0: 2.33e-15 m/s; exact kappa: 1.49e-7 m/s, completes t = 0.05; the integral model: 5.65e-2 after one step, FPE near 0.0375 s | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L164-L168), MEASURED |
| the step-1 kick under exact curvature | 2.15e-3 to 1.69e-9, factor 1.27e6 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102), MEASURED |
| the exact arms over 8000 steps | ratio 838.8, U0 = 0: 1.08e-8 to 1.46e-11; U0 = 0.05: 4.15e-9 to 3.00e-5; ratio 1: 4.13e-11 to 1.39e-12 and 7.04e-11 to 7.77e-11 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L200-L209), MEASURED |
| the hanging-node gate | 145.470 Pa in four arms (exact 145.48, -6.9e-5); max over time L2 velocity 3.8e-10 (band, np 4), 4.0e-10 (ball), 5.3e-10 (band, np 1), 8.2e-10 (uniform); endpoints 1.3 to 1.8e-10 against 6.0e-12 | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1207-L1219), MEASURED |
| the pre-registered pass line was mis-set | 1e-12 at every step is below the solver's floor; the decision rests on refined against uniform | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1220-L1222), MEASURED |
| a constant kappa holds a mode-2 perturbation | amplitude 1.0000 over 3200 steps, Courant 1.2e-10, continuity 1e-12 | [PSH 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L42-L46), MEASURED |
| the N = 64 translating gate | 1.09e-8 m/s on the first step, 3.94e-8 over the run; 1.003e-2 for `integralSurfaceTension` | [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L251-L258), MEASURED |
| the residual is the PISO-like level | SAAMPLE's residual-driven pressure iteration reaches 1e-13 from the first step where a fixed corrector count leaves 1e-6 to 1e-8 | [PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L755-L771), MEASURED in the reference, HYPOTHESIS for this solver |

## Why it failed, or why we think so

The gate does not fail on orthogonal meshes. On skewed meshes it does: with exact curvature and analytic geometry a residual of about `3.5e-5` m/s at N = 64 and `7.7e-5` at N = 128 remains after the momentum predictor and the strict solver have each removed their share, insensitive to the force form, the corrector count, `rAUf` and the tolerance. The one operator never removed is the `fvc::reconstruct` of the velocity correction ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1226-L1269), [METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L794-L799), [[concepts/pressure-projection-and-linear-solvers]]).

## Decisions

- The rung is mandatory before a coupled run on a new force or assembly ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L445-L453)); the two constant-curvature models are oracles, not production models ([[models/surface-tension-force]]).
- The hex refinement route proceeds on the strength of this gate; the octree objection of `config/stationaryDroplet3Dwide.yaml` is amended ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1223-L1226)).

## Open questions

1. The skew-mesh residual with exact curvature ([METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L794-L799)).
2. Whether residual-driven pressure iteration lowers the first-step floor from `1.09e-8` toward `1e-13` here as in SAAMPLE ([PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L764-L771)).

## Related

[[hubs/surface-tension]], [[hubs/verification]], [[models/surface-tension-force]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[concepts/pressure-projection-and-linear-solvers]], [[concepts/static-local-refinement]], [[concepts/method-gates]], [[cases/stationary-droplet]], [[cases/translating-droplet]], [[studies/poly3d-roadmap]], [[studies/sdpls-pre-print]].

## Log

### 2026-09-28
Created from CLAUDE's research loop, STATUS 0 and 4, the SL article's well-balanced gate and the roadmap's frozen-circle and translating gates.
