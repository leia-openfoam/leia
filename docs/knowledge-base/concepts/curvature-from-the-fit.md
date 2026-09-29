---
title: "Curvature from the quadratic fit"
description: "The curvature is read symbolically from the gradient and Hessian of the quadratic value fit at the cell centre; the parallel-curve offset correction turns the first-order 11.5 percent band error into a second-order 0.48 percent on a clean circle, and on the moving interface the estimator does not converge (2026-09-28)."
aliases: [symbolic curvature, parallel-curve offset correction, kappa from the fit]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-07-31
decided_by: [config/curvatureDroplet2D.yaml, config/faceCurvatureDroplet2D.yaml, METHOD 4]
code: [src/leiaLevelSet/semiLagrangian/slAdvection.C, src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C, applications/test/leiaTestMeanCurvature/leiaTestMeanCurvature.C, workflow/Snakefile.curvature]
sources: [METHOD 4, SL article sec:surften, SL article sec:curvature, PCS 18.1, STATUS 0, RM face gate, gcls article results]
---
# Curvature from the quadratic fit

> Verdict (2026-09-28). The cell curvature is the closed form `kappa = (tr(H) |g|^2 - g^T H g)/|g|^3`, evaluated with the gradient `g` and the Hessian `H` that the quadratic value fit already carries; no field is differentiated again ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L148-L159), [SL article `sec:surften`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L933-L946)). The value belongs to the level contour through the cell centre, which is a parallel curve of the interface at the signed offset `d`. That offset is a first-order error of about `(3/4) h/R` over the force band: 11.48 percent at N = 64 and 0.890 percent at N = 512 on a clean circle ([METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L187-L197)). The inverse `kappa = kappa_d/(1 - d kappa_d)` with `d = psi_c/|g|` gives 0.477 percent at N = 64 and 0.00703 percent at N = 512, with ratios 4.15, 4.00 and 3.91 per mesh halving: second order ([METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L203-L214)). The `O(h^1.2)` of the article was the offset error, not the accuracy of the fit's Hessian ([METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L214-L217)). Today the correction is applied by `curvatureExtension cellCentreInverse` ([[concepts/cell-centre-inverse-curvature]]), and the reconstruction's own `offsetCorrection` stays `none` ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167), [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L267-L268)). On the moving interface the estimator does not converge: 2.4 to 3.6 percent of `1/R` at all three rungs of the translating arm, against order 1.9 on the stationary droplet ([gcls article](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L636-L640)).

## What it is

In each band cell the fit is `R_c(x_c + d) = psi_c + g.d + 1/2 d^T H d`, so `g` and `H` are coefficients, not derivatives of a field ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L933-L936)). Only the uncached quadratic reconstruction implements this path; the base class throws `NotImplemented`, which is why a `GEOMETRY_FIT` token exists to keep the quadratic curvature while the transport fit varies ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L273-L291), [[models/sl-reconstruction]]). The fill is restricted: a cell gets a curvature only when its fit is full rank, `|g|` is not degenerate and `|psi_c|/|g|` is within three stencil radii of the interface; every other cell keeps exactly zero ([`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L82-L87)). The cell field is interpolated linearly to the faces and enters the balanced-force flux ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L965-L966), [[concepts/balanced-force-csf-flux]]).

Two evaluation points were tried. The article evaluates the formula at the Newton foot point of the fit's zero set and extends it constant along the normal ([SL article `eq:footNewton`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L948-L962)); that projection was measured to amplify the band curvature 40 to 85 times on a drifted level set and is off ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L157-L159), [PCS 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L365-L366)). The cell-centre value plus the algebraic inverse is what runs now. In 3D the inverse needs the Gaussian curvature `K = (g.cof(H).g)/|g|^4`, also from the fit ([METHOD 4.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L222-L239)).

## Why it matters

The curvature error is the source of the parasitic current ([[concepts/parasitic-current-mechanism]]): the band error `kErrL2Band = 70.6` against `kappa = 1000`, 7 percent at R/h = 12.8, is a property of the initial circle and identical in all eight arms of the kick-origin gate ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L100-L102)). The delivered face error collapses onto one function of the angle between the face normal and the interface normal: +3.33 at alignment falling to -0.99 at 45 to 60 degrees. The mesh offers two face orientations while the interface normal rotates, so the error is the `cos(4 theta)` pole bias; `m = 4` carries 99.7 percent of the smooth fluctuation, resolved at 20 cells per wave at N = 128, and mesh-locked, so refinement only shrinks its amplitude at `h^2` ([PCS 18.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1495-L1510)). See [[concepts/curvature-corrugation-and-the-fit]] for what the fit does to a perturbed level set.

## Where in the code

1. The symbolic fill: `slAdvection::meanCurvatureNoExtension` ([`slAdvection.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slAdvection.C#L140-L158)); the foot-point search `footPointDistance` in [`uncachedQuadraticWeightedLeastSquaresReconstruction.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L865-L983) ([PCS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L780-L784)).
2. The static measurement: `leiaTestMeanCurvature` compares the band curvature with `1/R` on the cells adjacent to a face with `snGrad(alpha) != 0` ([`leiaTestMeanCurvature.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestMeanCurvature/leiaTestMeanCurvature.C#L13-L27)); `make curvature` runs [`config/curvatureDroplet2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/curvatureDroplet2D.yaml) through [`workflow/Snakefile.curvature`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/Snakefile.curvature#L1-L58) and writes `curvature_error.csv` and `sl_curvature_error.png` into the SL theme.
3. The tokens: `SL_OFFSET_CORRECTION none`, `GEOMETRY_FIT none` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L267-L291)).

## Evidence

| claim | number | where |
|---|---|---|
| the band error of the raw cell value is the predicted offset error | predicted / measured 11.66 / 11.48, 6.02 / 6.03, 3.01 / 2.98, 0.88 / 0.89 percent at N = 64 / 128 / 256 / 512 | [METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L187-L197), DERIVED and MEASURED |
| the corrected value is second order on a clean circle | 2.11 / 0.477 / 0.111 / 0.0276 / 0.00703 percent at N = 32 to 512 with d = psi_c/abs(g); ratios 4.15, 4.00, 3.91 | [METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L203-L214), MEASURED |
| the article's uncorrected rate | L2 order 1.2, relative error 35 to 1 percent over N = 32 to 512; extension 1.22 against 1.29; face-averaged 1.18 | [SL article `sec:curvature`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1614-L1634), MEASURED |
| the full estimator against tr(H) | 35 against 48 percent at N = 32; equal limit on a signed distance | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1618-L1622), MEASURED |
| the face curvature of the arithmetic delivery, and with the foot point | h^1.13, 11.35 1/m at N = 512; h^2.04, 0.105 1/m | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1310-L1313), MEASURED |
| the Newton-foot and closest-point extensions leave the face error unchanged | h^1.08 and h^1.13, 11.35 1/m | [RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1317-L1328), MEASURED |
| the step-1 curvature error of the coupled runs | 70.6 against 1000 1/m at R/h = 12.8 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L100-L102), MEASURED |
| the error is the m = 4 pole bias | 99.7 percent of the smooth fluctuation; +3.33 aligned, -0.99 at 45 to 60 degrees | [PCS 18.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1497-L1510), MEASURED |
| the moving interface | 2.4 to 3.6 percent at N = 100, 142, 200 translating; order 1.9 stationary | [gcls article](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L636-L640), MEASURED |

## Why it failed, or why we think so

The fit's Hessian is second order; the first-order rate came from where the curvature was evaluated. Under transport the level sets stop being parallel offsets, `|grad psi|` spreads over 0.84 to 1.37 in the band ([METHOD 4.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L181-L183)), and both the fit and the inverse that follows it degrade; the foliation residual supplies about 21 percent of the curvature-error growth in a coupled run and the fit degradation the rest ([PCS 14.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1186-L1192), [[concepts/curvature-corrugation-and-the-fit]]).

## Decisions

- The production curvature is the symbolic value of the uncached quadratic fit, corrected per cell by `cellCentreInverse` ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), [[decisions/curvature-extension-cell-centre-inverse]]).
- The foot-point Newton projection stays off ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L166-L167)).

## Open questions

1. Why the curvature of the moving interface does not converge on the translating arm ([gcls article](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L923-L926), [[cases/translating-droplet]]).
2. The signed-distance assumption survives in the phase indicator, which locates the interface with a linear fit and the first-order offset `psi/|grad psi|`; the Hessian-corrected root exists as `offsetDistance` and is unmeasured ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2661-L2668), [PCS 11.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L863-L878), [[models/phase-indicator]]).
3. Whether any interface shape makes the delivered curvature exactly constant, the fixed-point question of Popinet's height function ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L518-L521)).

## Related

[[hubs/surface-tension]], [[models/surface-tension-force]], [[models/curvature-extension]], [[models/sl-reconstruction]], [[models/phase-indicator]], [[concepts/cell-centre-inverse-curvature]], [[concepts/face-curvature-deliveries]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[concepts/curvature-corrugation-and-the-fit]], [[cases/curvature-static-gates]], [[cases/translating-droplet]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from METHOD 4, the SL article sections on surface tension and curvature accuracy, PCS 18.1 and the gcls results.
