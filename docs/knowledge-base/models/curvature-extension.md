---
title: "curvatureExtension: the delivery words"
description: "The fourteen curvature-extension words of the SL two-phase solver, each with its order, its gain and its coupled verdict; cellCentreInverse is the global default, measured on the stationary droplet only (2026-09-28)."
aliases: [curvatureExtension, CURVATURE_EXTENSION, curvature delivery words]
kind: model
status: open
part: surface-tension
tags: [model, part/surface-tension]
date: 2026-09-28
date_settled:
decided_by:
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H]
sources: [METHOD 4, METHOD 8.1 row CURVATURE_EXTENSION, STATUS 4, STATUS 11.13, PCS sections 6-18, RM face gates, SL negative-results deck]
---
# curvatureExtension: the delivery words

> Verdict (2026-09-28). The global default is `cellCentreInverse` since c935883 ([METHOD 4.1, CORRECTED](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167)). Its basis is the stationary 2D ladder only: the unabsorbed capillary residual is 4.60 / 3.81 / 1.55 / 1.49x lower than with `none` at N = 32 / 64 / 128 / 256 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L242-L244)). For the translating droplet the choice is undecided, because the only `none`-against-`cellCentreInverse` comparison ran on the closed box and is void ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), [[retractions/closed-box-translating-droplet]]). The face deliveries are closed as a lever: the per-face inverse `stabilizedFootPointFace` is the only one that meets the acceptance criterion (order 1.98 and `G h^2` 0.647 on the ellipse), and every construction that lowered the gain lost an order ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1119-L1154)). The Popinet cases run `none` at density ratio 1 by reproduction, not by measurement ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)).

## What it is

`levelSet.curvatureExtension.type` in `fvSolution` selects how the SL two-phase solver turns the symbolic cell curvature of the quadratic fit ([[concepts/curvature-from-the-fit]]) into the curvature that the CSF flux applies. The word dispatch is a chain of `const bool` flags in [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L97-L104), not a runtime-selection table. The token is `CURVATURE_EXTENSION` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L384-L400)). Two families exist:

1. Cell-field words. They change the cell curvature and deliver it by arithmetic interpolation (`faceCurvatureSource model`).
2. Face-delivery words. They fill the registered face field `kappaStableFootFace`, consumed with `faceCurvatureSource registered` ([`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L67-L86)). The six face words require `semiLagrangian offsetCorrection none`; the solver refuses anything else, because the inverse would be applied twice ([`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L393-L410)).

## Members

Orders are the fitted active-face L2 orders on N >= 128; `G h^2` is the curvature noise gain of `leiaTestCurvatureNoiseGain` ([PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L731-L767)).

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| no extension, the cell-centre value | `none` | settled, the Popinet setting | the classical delivery: circle h^1.13 at 11.35 (N = 512), ellipse 0.97; the most stable of the three original words, stable window 0.06 to 0.11 s at N = 128 | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1297-L1361), [SL article stabilisation](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1754-L1778) |
| cell-centre parallel-surface inverse | `cellCentreInverse` | settled for the stationary droplet, open elsewhere | non-gradient force content converges at order +2.01 against +0.09 for every other cell field; residual 4.60 to 1.49x lower on the 2D ladder; second order on constant curvature only | [`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L36-L61), [[concepts/cell-centre-inverse-curvature]] |
| Gaussian term of that inverse | `gaussianCurvature yes` (sub-key) | settled | K-aware inverse h^1.95 against h^1.02 without K on the sphere; K off, two 3D arms diverge | [RM K-aware](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1362-L1414), [STATUS K](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L846-L890) |
| closest-point Newton foot | `closestPointNewton` | retracted | no static gain (h^1.08 / h^1.13); N = 64 floor 3.0e-6 at t = 0.14 then blow-up at 0.279 against 0.44; the projection amplifies band curvature 40 to 85x | [PCS dead ends](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L365-L366), [negative deck 3/3](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/3/3) |
| harmonic Laplace smoothing | `harmonicLaplace` (sub-key `relax`) | retracted | N = 64 blows at 0.051 s (tight solve) or 0.013 s (loose) against 0.44 s; coupling correlates the error force | [negative deck 2/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/2/1), [template](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L317-L327) |
| foot-point height function | `footPointHeightFunction` | retracted, 2D only | static h^1.75 (0.56 at N = 512); coupled blow-up at 0.006 s (N = 64) and 0.003 s (N = 128) | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1314), [negative deck 3/2](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/3/2) |
| connected zero-contour curvature | `connectedInterface` | retracted, 2D serial only | quietest replay 4.55e-3 m/s (Helmholtz w 5, lambda 16) but halves the m = 2 mode; static h^0.48; the N = 64 translating run reaches 0.0893 m/s and the oscillating run dies at 0.0209 s | [RM connected](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L779-L829), [RM shared service](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L881-L903) |
| interface mean, spatially constant | `interfaceMean` | settled, diagnostic only | replay quiet at 3.84e-9 m/s: spatial variation drives the current, not the mean; zeroth order on the ellipse (506 at N = 512) | [RM replay](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L735-L751), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400) |
| skip the symbolic fill | `fvm` | settled, research use | the force model computes its own curvature; required for `divGradPsiSnGradAlpha`, `isoCurvature` and linear fits; `kappa` stays 0 | [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L233-L239) |
| per-face parallel-surface inverse | `stabilizedFootPointFace` | settled, the only delivery that meets the criterion | circle h^2.04 (0.105, 108x), sphere h^1.95, ellipse 1.98 (0.2785); `G h^2` 0.647; coupled onset unchanged; t_blow 0.0668 s (N = 128), 0.0348 s (N = 256) | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1297-L1361), [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L990-L1063) |
| one inverse per cut cell | `cutCellFootPointFace` | retracted | best static (circle h^2.00, 0.0761) but `G h^2` 0.821 and blow-up 3.3x sooner (0.0202 s at N = 128); ellipse collapses to 1.02 | [PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L666-L768), [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L990-L1063) |
| cell mean of the per-face inversions | `cellMeanFootPointFace` | retracted | lowest gain (`G h^2` 0.402), longest survival at N = 128 (0.1049 s), but first order on the ellipse (1.03) and exponent -1.09: a prefactor win | [PCS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L902-L989), [[retractions/cell-mean-delivery-adoption]] |
| symmetric face mean of the inversions | `symmetricFaceMeanFootPointFace` (sub-key `faceSmoothing`, default 0.5) | retracted | theta 0.5: 2.014 at order 1.10, `G h^2` 0.445; theta 1.0: 4.009 at 1.07; the active mask is asymmetric | [PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1119-L1154) |
| fit curvature at the face foot point | `footPointEvaluatedFace` | retracted | exact where the fit reproduces psi (4.8e-4 on the quadratic-form ellipse) but first order on true distances; coupled: blows up earlier at every N (t_blow ratio 0.36 / 0.42 / 0.32), volume 5 to 15x worse | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L417-L460), [PCS 17.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1445-L1484) |
| the same, once per cell at its centre | `cellFootPointEvaluatedFace` | open | rationale: restore the cell-field structure the projection can absorb; prediction recorded in `config/driverSplitDeliveryProbe.yaml`; no result recorded | [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L276-L289), [config](https://github.com/leia-openfoam/leia/blob/8867581/config/driverSplitDeliveryProbe.yaml#L1-L45) |

## Why it matters

The pressure projection absorbs the part of the capillary flux that is `snGrad` of a cell field. Only the face-to-face variation of `kappa_f` reaches the velocity ([[concepts/balanced-force-csf-flux]]). The delivery word decides how much variation the estimator error produces at the faces, and how strongly a perturbation of psi moves `kappa_f` (the gain). Both numbers are measured on static gates in seconds ([[cases/curvature-static-gates]]). The campaign found that accuracy and gain are dissociated: 27x more accurate for 0.2 percent of gain change, and every gain reduction cost an order ([PCS 15.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1260-L1299)).

## Where in the code

- Dispatch and per-word comments: [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L97-L420).
- The cell-centre inverse: [`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L1-L61); the face deliveries: `stabilizedFootPointFaceCurvature.H`; the others: `footPointCurvature.H`, `connectedInterfaceCurvature.H`, `interfaceMeanCurvature.H`, `harmonicMeanCurvatureExtensionEqn.H`.
- The symbolic fill and the closest-point foot: `slAdvection::meanCurvatureNoExtension` and `meanCurvatureClosestPoint` ([`slAdvection.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slAdvection.C#L140-L158)).
- Temporal under-relaxation `curvatureExtension.relax` in [`slAlphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L314-L334).
- The Eulerian solver has no dispatch; it runs the older closest-point fill plus `stabilizedFootPointFace` ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3704-L3707)).

## Evidence

| claim | number | where |
|---|---|---|
| the varying-curvature ellipse gate ranks the deliveries | arithmetic 0.97, per-face inverse 1.98, cut-cell 1.02, cell-mean 1.03, interface mean -0.01 | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400), MEASURED |
| the acceptance criterion | `G h^2 <= 0.65` and order `>= 1.9` on the ellipse gate, not on the circle | [STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2685-L2688), [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1044-L1047), MEASURED |
| the gain orders the coupled growth rate | equal gain (0.644, 0.643) gives equal rate (186, 189 1/s) despite 27x accuracy; +28 percent gain gives 2.4x rate | [PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L731-L767), MEASURED |
| no delivery changed the t_blow exponent | per-face -0.94, cut-cell -0.48, cell-mean -1.09 | [PCS 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064-L1116), MEASURED |
| `cellCentreInverse` on the translating droplet after the fix | `none` diverges 20 percent earlier (0.0695 s against 0.0868 s), one resolution, np 4 | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3743), MEASURED, indicator only |
| the oscillating droplet at N = 128 | `cellCentreInverse` completed where `none` failed at 0.0982 s; the study used an algebraic psi | [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), MEASURED, see [[cases/oscillating-droplet]] |

## Why it failed, or why we think so

Assigning one value to every active face of a cut cell forces `kappa_f` to be constant across a cell, an O(h) offset wherever the exact curvature varies ([PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1041-L1049)). Symmetric averaging about a face fails for a different reason: the active-face mask depends on where the interface sits inside the cell, so the mean is offset at O(h) ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1137-L1142)). Every corrective built from the same quadratic fit injects that fit's error, and the loop amplifies the injection faster than the removed systematic error ([PCS 15.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1268-L1275)).

## Decisions

- `CURVATURE_EXTENSION cellCentreInverse` as the global default since c935883, case-dependent where a measurement says so; no `cases/<case>.parameter` sets it ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3573-L3574)). See [[decisions/curvature-extension-cell-centre-inverse]] and [[decisions/curvature-inverse-gaussian]].
- The gate's translating arm runs `cellCentreInverse` to 0.05 s; the `none` of 8a9b85a was reverted because its evidence is void ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3613-L3616)).

## Open questions

1. `none` against `cellCentreInverse` on a matched translating setup after the fix 440107f ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3620-L3622)).
2. `cellCentreInverse` was never scored on the ellipse gate ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)).
3. A 3D varying-curvature companion (a torus) to confirm the ranking in 3D ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2669-L2675)).

## Related

[[hubs/surface-tension]], [[models/surface-tension-force]], [[concepts/curvature-from-the-fit]], [[concepts/cell-centre-inverse-curvature]], [[concepts/face-curvature-deliveries]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[cases/curvature-static-gates]], [[decisions/curvature-extension-cell-centre-inverse]], [[decisions/curvature-inverse-gaussian]], [[retractions/cell-mean-delivery-adoption]], [[retractions/closed-box-translating-droplet]], [[studies/curvature-stabilization-campaign]].

## Log

### 2026-09-28
Created; dictionary words read from the flag chain of `createSLFields.H`, verdicts from PCS sections 6 to 18, STATUS 4 and METHOD 8.1.
