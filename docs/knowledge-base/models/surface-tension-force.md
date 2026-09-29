---
title: "surfaceTensionForce: the twelve models"
description: "The runtime-selectable capillary force family; reconstructedCurvature is production, the exact-curvature models are oracles, and every alternative estimator or integral model is closed (2026-09-28)."
aliases: [surface tension force models, SURFACE_TENSION_FORCE]
kind: model
status: settled
part: surface-tension
tags: [model, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-27
decided_by: [METHOD 8.1 row SURFACE_TENSION_FORCE, config/gates/methodGate2D.yaml]
code: [src/leiaLevelSet/surfaceTensionForce, applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H, applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H]
sources: [METHOD 5, METHOD 8.1, SL article sec:surften, SL article sec:droplet, RM face-flux matrix, SL negative-results deck]
---
# surfaceTensionForce: the twelve models

> Verdict (2026-09-28). The production model is `reconstructedCurvature` with the arithmetic face value and the `alpha` weight ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L389-L390)). It delivers the capillary force as one integrated face flux `G_sigma,f = sigma kappa_f snGrad(alpha) |S_f|` that the pressure equation absorbs exactly when `kappa_f` is constant ([METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L243-L272)). The two constant-curvature models are oracles, not methods ([RM frozen-circle gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1028-L1058)). Every other member is a closed research line: the FVM and trace estimators survive the N = 64 translating matrix with a disturbance larger than the translation speed, and the alpha-CSF, the two integral models and the iso-conormal model diverge before 4 ms ([`transISTN64ForceFluxModelMatrix.csv`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv)). The rule that orders the family: the better the static balance, the higher the dynamic gain ([SL article, delivery study](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1780-L1818)).

## What it is

`surfaceTensionForce` is the runtime-selectable family behind `levelSet.surfaceTensionForce.type` in `fvSolution`. Every member returns the owner-oriented scalar face flux `G_sigma,f` ([`surfaceTensionForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.H#L29-L106)). The solver adds that flux to `phig` in the pressure equation and corrects the velocity from the face-wise difference `phig - p_rghEqn.flux()` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L62-L102)). The force never enters the momentum matrix ([`UEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H#L8)). The library is `libleiaSurfaceTension` ([`Make/files`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/Make/files)); it also carries the semi-implicit fvOption, see [[models/semi-implicit-capillary-force]].

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| reconstructed curvature | `reconstructedCurvature` | settled, production | CSF flux from the registered cell curvature; a constant `kappa` is absorbed to 3e-11 m/s; the parasitic current is the tangential variation of `kappa_f` | [METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L243-L272), [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L73-L133) |
| constant curvature, CSF form | `constantCurvatureSurfaceTension` | settled, oracle | exact `kappa = 1/R` cuts the step-1 kick 2.15e-3 to 1.69e-9 (1.27e6x); bounded at every density ratio and U0 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102), [STATUS amplifierGate](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L228) |
| constant curvature, potential form | `constantCurvaturePressurePotential` | settled, oracle | `snGrad(sigma kappa alpha) magSf`; agrees with the CSF form to 2e-14 to 1.5e-13 on uniform and perturbed meshes | [RM frozen circle](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1028-L1058) |
| computed-curvature potential | `curvaturePressurePotential` | open, diagnostic | the remainder-off arm with the computed `kappa`; cannot relax a deformed interface; no result recorded | [header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/curvaturePressurePotential.H#L29-L64), [`config/remainderOffBalance.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/remainderOffBalance.yaml#L1-L45) |
| Kang face value on an FVM curvature | `correctionKang` | retracted as production | replay 1.91e-3 to 2.50e-2 m/s; N = 64 translating matrix reaches 0.05 s with max disturbance 0.167 m/s, above U0 = 0.05 m/s | [RM replay](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L724-L778), [matrix CSV](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv) |
| interFoam operator on the geometric alpha | `divGradAlphaSnGradAlpha` | retracted | 0.14 to 0.84 m/s after one step, 15 to 150 m/s, FPE within 148 to 493 steps; the operator alone bounds nothing | [PSH 0f](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L658-L679), [negative deck 3/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/3/1) |
| FVM curvature of psi | `divGradPsiSnGradAlpha` | settled, research line | face error h^1.16 (11.36 at N = 512); N = 64 standing current 8e-4 m/s, 60x the reconstruction; N = 128 diverges at 0.116 s | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1297-L1361), [negative deck 3/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/3/0) |
| trace of the gradient of the geometric normal | `traceGradGeoNormalSnGradAlpha` | open, no separate record | identical to `divGradPsiSnGradAlpha` in the N = 64 matrix (9.06e-2 m/s, reaches 0.05 s) | [matrix CSV](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv) |
| trace of the gradient of grad psi | `traceGradGradPsiSnGradAlpha` | open, no separate record | identical to `divGradPsiSnGradAlpha` in the N = 64 matrix (9.06e-2 m/s) | [matrix CSV](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv) |
| integral surface tension (CST) | `integralSurfaceTension` | retracted, 2D structured prototype | first reachable equilibrium at N = 64 (floor 5e-7 m/s), N = 128 diverges at 0.047 s; translating FPE at 9.3e-4 s | [SL article delivery study](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1789-L1818), [[concepts/integral-surface-tension-cst]] |
| integral conormal traction | `integralConormalSurfaceTension` | retracted, pseudo-2D prototype | curvature-free; pressure-range residual 0.4945; fails at t = 0.001124 s at N = 64 | [RM conormal](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L214-L257), [RM bridge](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L411-L443) |
| iso-value conormal curvature | `isoCurvature` | retracted, 2D prototype | static face error 208.7 at N = 512 at order 0.14 (`isoConormal`); translating FPE at 1.42e-3 s | [`face_curvature_orders.csv`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/face_curvature_orders.csv), [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L647-L649) |
| semi-implicit fvOption | `semiImplicitCapillaryForce` | candidate | see [[models/semi-implicit-capillary-force]] | [`semiImplicitCapillaryForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L26-L108) |

Note on the base class: `surfaceTensionForce.H` declares `TypeName("none")` but registers no model of that name ([`surfaceTensionForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.H#L106)). The entry is required in every case, and `none` is not a selectable word.

### The sub-options of `reconstructedCurvature`

| key | words | production | what it changes |
|---|---|---|---|
| `faceInterpolation` | `arithmetic`, `interfaceWeighted`, `connectedInterface` | `arithmetic` | the face value of `kappa`; `interfaceWeighted` is the Kang/GFM weighting, see [[concepts/kang-gfm-and-sharp-heaviside]] |
| `forceWeight` | `alpha`, `sharpHeaviside` | `alpha` | the field whose `snGrad` localises the force |
| `faceCurvatureSource` | `model`, `registered` | `model` (token `FACE_CURVATURE_SOURCE`) | `registered` consumes a solver-filled face field such as `kappaStableFootFace`, see [[concepts/face-curvature-deliveries]] |

Locators: [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L73-L133) and [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L469-L472).

## Why it matters

The parasitic current of a balanced-force coupling is set by the curvature error and by nothing else in the force chain. The exact-curvature arms of the amplifier gate stay bounded at every combination of density ratio and translation speed, and the reconstructed-curvature arms are the ones that destabilise ([STATUS amplifierGate](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L228)). The family therefore isolates two questions: does the force assembly balance (the oracles), and does the estimator feed the coupling a usable curvature (every other member). See [[concepts/parasitic-current-mechanism]].

## Where in the code

1. Selection: `levelSet { surfaceTensionForce { type ...; } }` in `fvSolution`, rendered from the token `SURFACE_TENSION_FORCE` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L472-L480), [template](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L232-L275)).
2. Assembly of the flux: [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L73-L191); the constant-curvature form is one line, `sigma*curvature*snGrad(alpha)*magSf` ([`constantCurvatureSurfTension.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/constantCurvatureSurfTension.C#L69)).
3. Consumption: `UEqn.H` line 8 calls `faceSurfaceTensionForceFlux()`; `pEqn.H` adds `GSigma` to `phig`, solves `laplacian(rAUf, p_rgh) == div(phiHbyA)` and reconstructs the velocity from the residual flux ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L62-L102)).
4. The flux-space residual `R_f = phig - p_rghEqn.flux()` on the active faces is written to `capillaryFluxResidual.csv` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L170-L201)).

## Evidence

| claim | number | where |
|---|---|---|
| a constant curvature gives zero spurious velocity on a frozen circle | max abs U about 3e-11 m/s | [METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L259-L262), MEASURED |
| the exact curvature removes the step-1 kick | 2.15e-3 to 1.69e-9, factor 1.27e6 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102), MEASURED |
| only the exact control stays balanced on the N = 64 translating matrix | exact 2.17e-7 m/s; FVM and trace 9.06e-2; Kang 0.167; the five others FPE before 4 ms | [matrix CSV](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv), [MC deck 8/1](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/8/1), MEASURED |
| the static face curvature of the production arithmetic delivery | order h^1.13, 11.35 1/m at N = 512; with the foot point h^2.04, 0.105 | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1297-L1361), MEASURED |
| the projection absorbs almost the whole flux | 99.85 to 99.98 percent; 30x tighter solve moves the residual by at most 2.4e-4 | [SL article sec:fluxresidual](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2200-L2273), MEASURED |
| the shared face service gives identical replays for six CSF members | same `abs(U - Utrans)` to the printed digits with the same connected curvature | [RM shared service](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L830-L851), MEASURED |
| the delivery ordering rule | the better the static balance, the higher the dynamic gain | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1810-L1818), [negative deck 4/3](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/3), MEASURED |

## Why it failed, or why we think so

Every alternative member reads the transported level-set profile or a sharp geometric alpha through a differentiating operator. The exact-curvature arms prove that the pressure-velocity path is balanced ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L211-L219)). The operator swap to interFoam's `-div(nHat)` on the sharp geometric alpha runs away at once, so no operator bounds the loop by itself; interFoam is bounded by the pairing of operator, transport and state ([PSH 0f](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L594-L679)). See [[concepts/variational-capillary-force]] for the proposed pairing.

## Decisions

- `SURFACE_TENSION_FORCE reconstructedCurvature`, `FACE_CURVATURE_SOURCE model`, per case ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L389-L390)); see [[decisions/surface-tension-reconstructed-curvature]].
- The `directCell` delivery of the conormal model was removed; only the face-flux path remains ([`integralConormalSurfaceTension.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/integralConormalSurfaceTension.C#L52-L71)).

## Open questions

1. `curvaturePressurePotential` has a pre-registered prediction (`RfFracL2` collapses from 1.7e-5 to about 1e-9) and no recorded result ([`config/remainderOffBalance.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/remainderOffBalance.yaml#L28-L33)).
2. The two trace models have no static or coupled record separate from `divGradPsiSnGradAlpha`.

## Related

[[hubs/surface-tension]], [[models/curvature-extension]], [[models/semi-implicit-capillary-force]], [[concepts/balanced-force-csf-flux]], [[concepts/curvature-from-the-fit]], [[concepts/face-curvature-deliveries]], [[concepts/kang-gfm-and-sharp-heaviside]], [[concepts/integral-surface-tension-cst]], [[concepts/variational-capillary-force]], [[concepts/well-balanced-exact-curvature-gate]], [[concepts/parasitic-current-mechanism]], [[decisions/surface-tension-reconstructed-curvature]], [[cases/stationary-droplet]], [[studies/poly3d-roadmap]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from METHOD 5 and 8.1, the SL article, the roadmap gates and the force-flux model matrix; dictionary words taken from the `TypeName` strings of `src/leiaLevelSet/surfaceTensionForce`.
