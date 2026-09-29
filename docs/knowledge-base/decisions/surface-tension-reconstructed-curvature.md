---
title: "SURFACE_TENSION_FORCE: the balanced-force flux from the reconstructed curvature"
description: "SURFACE_TENSION_FORCE reconstructedCurvature with FACE_CURVATURE_SOURCE model, decided by the balanced-force gates: a constant curvature is absorbed to 3e-11 m/s on a frozen circle and the Laplace jump is 145.470 Pa in every arm of the well-balanced gate; the parasitic current that remains is the variation of the discrete curvature"
aliases: [SURFACE_TENSION_FORCE reconstructedCurvature, FACE_CURVATURE_SOURCE model]
kind: decision
status: settled
part: surface-tension
tags: [decision, part/surface-tension]
date: 2026-09-28
date_settled: 2026-07-31
decided_by: [workflow/scripts/run_pressure_compatibility_gate.py, config/stationaryDroplet3DrefinedWB.yaml, "author decision 2026-07-31"]
code: [src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C, applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H, cases/default.parameter, cases/stationaryDroplet2D/system/fvSolution.template]
sources: ["METHOD 5 (L243-L272)", "METHOD 8.1 rows SURFACE_TENSION_FORCE and FACE_CURVATURE_SOURCE (L389-L390)", "METHOD 10 (L807-L837)", "SL article sec:surften (L1013-L1078)", "RM frozen-circle gate (L1028-L1058)", "STATUS 0 (L97-L102)", "DP L471-L484"]
---
# SURFACE_TENSION_FORCE: the balanced-force flux from the reconstructed curvature

> `SURFACE_TENSION_FORCE reconstructedCurvature` and `FACE_CURVATURE_SOURCE model` in the global default ([DP L471-L484](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L471-L484)); METHOD 8.1 lists the layer as per-case, but at 8867581 no droplet case sets the token, and only the two integral-surface-tension cases override it ([METHOD 8.1 L389-L390](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L389-L390), [transISTDroplet2D.parameter L58](https://github.com/leia-openfoam/leia/blob/8867581/cases/transISTDroplet2D.parameter#L58)). The model was production on 2026-07-31, the date of the METHOD.md dictionary ([METHOD 10 L834-L835](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L834-L835)). Decided by the balanced-force gates ([METHOD 5 L243-L272](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L243-L272)): the force is one integrated oriented face flux `G_sigma,f = sigma kappa_f snGrad(alpha) abs(S_f)` that enters the pressure equation on the faces where the pressure gradient acts, so a spatially constant curvature is a discrete gradient that the pressure absorbs exactly; measured at max abs U about 3e-11 m/s on a frozen circle. The frozen-circle gate of 2026-07-28 gives 3.770e-9 to 1.039e-8 m/s on uniform meshes at N = 32 to 128 ([RM L1028-L1058](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1028-L1058)), and the well-balanced refinement gate gives a Laplace jump of 145.470 Pa in every arm, 7e-5 below the exact value ([SL article L1550-L1560](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1550-L1560)). What remains is the tangential variation of the delivered curvature, which no pressure can balance ([SL article L1071-L1077](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1071-L1077)).

## The question

Which capillary force model, and which curvature enters it? The twelve models of the family return the same object, an owner-oriented face flux ([[models/surface-tension-force]]). `reconstructedCurvature` reads the curvature that the quadratic value fit of the level set supplies in closed form, `kappa = (tr(H) abs(g)^2 - g^T H g) / abs(g)^3`, with no re-differentiation of a finite-volume field ([SL article L1027-L1041](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1027-L1041), [[concepts/curvature-from-the-fit]]). `FACE_CURVATURE_SOURCE model` lets the model assemble its own face value by arithmetic interpolation; `registered` reads a solver-registered face field of a face delivery ([DP L471-L473](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L471-L473), [[models/curvature-extension]]). The alternatives are the finite-volume estimators `divGradPsiSnGradAlpha` and the trace forms, the Kang face value, the integral models and the iso-conormal model.

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| frozen 1 mm circle, exact constant curvature, one capillary step, uniform mesh, N = 32 / 64 / 128 | max cell velocity | 3.770e-9 / 7.635e-9 / 1.039e-8 m/s; the CSF form and the potential form agree to 2.02e-14 to 1.50e-13 | MEASURED, [RM L1028-L1058](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1028-L1058) |
| the same on 10 %-perturbed meshes | max cell velocity | 8.584e-6 / 7.328e-4 / 1.392e-3 m/s, growing with refinement | MEASURED, [RM L1028-L1058](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1028-L1058) |
| 3D well-balanced gate, kappa = 2/R prescribed, refined and uniform meshes, 921 steps | Laplace jump; spurious velocity | 145.470 Pa in every arm, 7e-5 below exact; `meanMagUPrime` and `l2MagUPrime` at or below 1e-12 m/s is the PASS criterion | MEASURED, [SL article L1550-L1560](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1550-L1560), [config L26-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3DrefinedWB.yaml#L26-L30) |
| 2D translating droplet, N = 128, exact kappa = 1/R against the reconstructed curvature | step-1 kick `L1(U - U0) / U_ref` | 2.15e-03 to 1.69e-09, a factor 1.27e+06; the curvature error is 7 % at R/h = 12.8 | MEASURED, [STATUS L97-L102](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102) |
| the N = 64 translating force-model matrix (transIST family; case and date not in the CSV) | completion; max abs (U - U_trans) | exact curvature reaches 0.05 s at 2.17e-07 m/s; `reconstructedCurvature` diverges at 3.78e-03 s; the three FVM estimators reach 0.05 s at 9.06e-02 m/s, 1.8x U0; Kang 0.167 m/s; the alpha-CSF, the two integral models and the iso-conormal model diverge before 1.5e-03 s | MEASURED, [transISTN64ForceFluxModelMatrix.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv) |
| delivery study, stationary droplet | Kang face value with a sharp Heaviside | static balance up to 58x better, but the sharp pairing diverges at t about 0.07 s where the arithmetic default survives to 0.44 s (N = 64) | MEASURED, [SL article L1951-L1960](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1951-L1960) |
| delivery study, CST integral model | equilibrium and divergence | floor 5e-7 m/s at N = 64 and t = 0.2 s; N = 128 diverges at 0.047 s against 0.105 s for the CSF default | MEASURED, [SL article L1970-L1985](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1970-L1985) |

Pre-registered read-out of the well-balanced gate: [stationaryDroplet3DrefinedWB.yaml L20-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3DrefinedWB.yaml#L20-L30). The frozen-circle gate scores the maximum written cell velocity only ([RM L1043-L1046](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1043-L1046)).

The rule that orders the family: the better the static balance, the higher the dynamic feedback gain; the arithmetic average and the smeared support act as low-pass filters on the curvature-error force, and delivery-side remedies cannot close the fine-mesh gap ([SL article L1985-L1990](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1985-L1990)).

## What it does not cover

1. The curvature estimator is the source of the parasitic current; the force delivery is not ([STATUS L74-L125](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L74-L125), [[concepts/parasitic-current-mechanism]]). The choice of the cell curvature and its extension is a separate decision ([[decisions/curvature-extension-cell-centre-inverse]]).
2. The perturbed-mesh residual with exact curvature grows under refinement and is insensitive to every lever but the `fvc::reconstruct` of the velocity correction ([METHOD 9 L795-L801](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L795-L801)).
3. The semi-implicit capillary fvOption is a candidate, not part of this decision; it is untested on the translating droplet since its collective fix ([[models/semi-implicit-capillary-force]]).
4. The variational force is a proposal without code ([[concepts/variational-capillary-force]]).
5. The N = 64 force-model matrix shows `reconstructedCurvature` itself diverging at 3.8 ms on that case; the CSV records no case name and no date, so the row is an indicator of the family ranking, not a ladder.

## Related

[[hubs/surface-tension]] - [[models/surface-tension-force]] - [[models/curvature-extension]] - [[concepts/balanced-force-csf-flux]] - [[concepts/curvature-from-the-fit]] - [[concepts/well-balanced-exact-curvature-gate]] - [[concepts/parasitic-current-mechanism]] - [[concepts/kang-gfm-and-sharp-heaviside]] - [[concepts/integral-surface-tension-cst]] - [[concepts/static-local-refinement]] - [[decisions/curvature-extension-cell-centre-inverse]] - [[decisions/curvature-inverse-gaussian]] - [[decisions/psi-filter-none]] - [[cases/stationary-droplet]] - [[studies/poly3d-roadmap]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED; production since 2026-07-31, the frozen-circle gate ran on 2026-07-28, the well-balanced refinement gate on 2026-09-04. Entered in [[decision-log#2026-07]].
