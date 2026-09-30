---
title: "The curvature stabilization campaign (plan v0.2, 2026-08)"
description: "The map of docs/plan-curvature-stabilization.md: the ground truth of section 0, the work packages, and the measured sections 6 to 18 (2026-08-07 to 08-17) that closed the delivery lever, the psi-side lever, and located the instability in the explicit capillary coupling of the m=2 mode."
aliases: []
kind: study
status: settled
part: surface-tension
tags: [study, part/surface-tension]
date: 2026-09-28
code: [docs/plan-curvature-stabilization.md, src/leiaLevelSet/surfaceTensionForce, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam, applications/test/leiaTestCurvatureNoiseGain, applications/test/leiaTestFoliationResidual]
sources: [PCS sections 0-18, STATUS 4 (08-19 to 08-31), kb-raw C2, kb-raw C8]
---
# The curvature stabilization campaign (plan v0.2, 2026-08)

> **Verdict (2026-09-28).** `docs/plan-curvature-stabilization.md` (1592 lines, first commit 2026-08-07, last 8693326 of 2026-08-27) is the record of the August campaign against the parasitic-current instability of the stationary droplet. Its measured sections 6 to 18 closed every operator lever: the curvature delivery is exhausted (the per-face inverse stays production, [section 14](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1117)), the psi-side renormalisation makes the run 3.4x worse ([section 15](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1211)), every coupling lever is exonerated ([section 16](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1300)), and the instability at N=128 is the explicit time coupling of the capillary force to the m=2 interface mode, with a Richardson intercept of zero ([section 18](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1485)). Three of its own conclusions were retracted inside the document ([Retracted or superseded](#retracted-or-superseded-inside-it)), and its transport ground truth carries a contamination notice ([section 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L34-L61)). The campaign's decisions survive in [[decisions/curvature-extension-cell-centre-inverse]], [[decisions/curvature-inverse-gaussian]], [[decisions/psi-filter-none]] and the acceptance criterion of [[concepts/face-curvature-deliveries]].

## What it is

A working document for a coding agent, version v0.2 (v0.1 plus the 2026-08-07 review amendments). Sections:

| section | title | line |
|---|---|---|
| 0 | Ground truth: transport orders (with the contamination notice of 2026-08-26 and the re-established 2D row of 08-27), curvature statics, the coupled blow-up times, the mechanism, the seam-consistency rule | [14](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L14) |
| 1 | Problem decomposition: A static seed (solved), B feedback exponent (open), C pinch-off | [117](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L117) |
| 2 | Work packages WP0 (band mode spectrum), WP1 (anisotropic geometry fit), WP2 (its integration), WP3 (the psi filter as a delay device), WP4 (3D Gaussian inverse, done), WP5 (normal constancy), WP6 (gated profile reset, with the binding rule v0.3), WP7 (evaluation matrix) | [135](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L135) |
| 3 | Twelve measured dead ends | [363](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L363) |
| 4, 5 | Conventions, open questions | [399](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L399), [422](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L422) |
| 6, 7 | WP0 retrodictions (a) N=128 arithmetic and (b) N=256 filtered, 2026-08-07 | [438](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L438), [487](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L487) |
| 8 | WP8.0: the projection is converged, the residual is structural, 2026-08-07 | [546](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L546) |
| 9 | WP8.1: the across-support variation regrows and dominates, 2026-08-08 | [604](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L604) |
| 10 | WP8.2: the cut-cell delivery falsifies the WP8.1 consequence; the noise gain G h^2, 2026-08-08 | [666](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L666) |
| 11 | WP8.3: what the inverse assumes, the foliation residual D, the cell-mean delivery, 2026-08-10 | [769](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L769) |
| 12 | WP8.4: the varying-curvature ellipse gate kills the cut-cell family, 2026-08-12 | [990](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L990) |
| 13 | WP7 refit: no delivery has changed the t_blow scaling, 2026-08-12 | [1064](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064) |
| 14 | WP8.5: symmetric face averaging fails; the foliation residual on the saved fields, 2026-08-12 | [1117](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1117) |
| 15 | WP9: band renormalisation 3.4x worse; the finding across four levers, 2026-08-12 | [1211](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1211) |
| 16 | The loop model, the capillary-dt sweep and the exonerations, 2026-08-13/14 | [1300](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1300) |
| 17 | The foliation gate on the ellipse; `footPointEvaluatedFace`, 2026-08-14 | [1393](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1393) |
| 18 | The source measured: a neutral loop destabilised by the explicit coupling of the m=2 mode, 2026-08-17 | [1485](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1485) |

The curated outputs are in the method-comparison data folder ([[studies/method-comparison]]): `face_curvature_orders*.csv`, `curvature_gain*.csv`, `face_curvature_orders_foliation.csv`, `mode_rate_vs_drift.csv`, `mode_rate_dt_series.csv`, `capillary_flux_residual.csv` (in the SL theme), and the figures `wp0_retrodiction_*.png`, `wp81_driver_split_N128.png`, `mode_rate_vs_drift.png`.

## Why it matters

The campaign replaced the intuition "a more accurate curvature stabilises the droplet" by four measured facts: static accuracy is anti-correlated with coupled stability (section 10); a delivery is scored by its noise gain and its order on a varying-curvature interface, not on a circle (sections 10, 12); the level set losing its parallel foliation drives the curvature-error growth (section 14.2); and the growth rate at N=128 is c dt with c raised by drift, so the lever is the time coupling of the force (section 18). These facts set the rules "no filtering in production" and "no partial solutions" ([[decisions/psi-filter-none]]) and the score of the later campaign ([[studies/shannon-parasitic-currents-campaign]]).

## Where in the code

`src/leiaLevelSet/surfaceTensionForce/` and the `curvatureExtension` switch of `leiaSemiLagrangianLevelSetTwoPhaseFoam` (`stabilizedFootPointFace`, `cutCellFootPointFace`, `cellMeanFootPointFace`, `symmetricFaceMeanFootPointFace`, `footPointEvaluatedFace`, `cellCentreInverse`), `applications/test/leiaTestCurvatureNoiseGain`, `applications/test/leiaTestFoliationResidual`, `applications/test/leiaTestDeparturePoint`, `workflow/scripts/loop_spectrum.py`, `interface_mode_trajectory.py`, `mode_rate_vs_drift.py`, `cases/ellipseDroplet2D`, `config/faceCurvatureEllipse2D.yaml`, `config/stationaryDropletDtSweep.yaml`.

## Evidence

| claim | number | where |
|---|---|---|
| Static delivery accuracy (circle, exact SDF) | cell-centred curvature h^1.07, relative L2 35.4 % to 1.7 %; stabilised foot-point face h^2.04, 11.35 to 0.105 1/m at N=512; 3D sphere with the Gaussian term h^1.95 (0.419 to 5.25e-3), without it h^1.02 | [section 0, lines 64-76](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L64-L76), MEASURED; [[concepts/curvature-from-the-fit]], [[decisions/curvature-inverse-gaussian]] |
| The projection is converged | a 30x tighter pressure solve moves the non-absorbable fraction, the velocity residual and max U by at most 2.4e-4 relative; the projection absorbs 99.85 to 99.98 % | [section 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L546), MEASURED; [[concepts/balanced-force-csf-flux]] |
| The across-support variation regrows | the foot-point delivery starts 57x below the arithmetic one (17.0 against 967 1/m) and reaches its magnitude within 0.05 s | [section 9](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L604), MEASURED |
| Accuracy anti-correlated with stability | static L2 45.8, 1.67, 1.18 1/m against t_blow 0.0803, 0.0668, 0.0202 s (arithmetic, per-face, cut-cell); the gain G h^2 0.644, 0.643, 0.821 orders them | [section 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L666), MEASURED |
| The inverse does not assume a signed distance; its residual is first order in the foliation defect D | torus check to 7 digits; on the 2:1 ellipse vertex the inverse is 3.01x worse than no correction | [section 11.1, 11.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L771), DERIVED and MEASURED |
| The cell-mean delivery survives longest on the circle | t_blow 0.1049 s at N=128 against 0.0668 (per-face); gain 0.402 | [section 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L902), MEASURED |
| The ellipse gate collapses one-value-per-cell deliveries | per-face inverse order 1.98; cut-cell 1.02; cell-mean 1.03; interface mean -0.01 | [section 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L990), MEASURED; [[cases/curvature-static-gates]] |
| No delivery changed the exponent | t_blow about N^p with p = -0.94, -0.48, -1.09 (per-face, cut-cell, cell-mean) | [section 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064), MEASURED |
| Symmetric face averaging | order 1.10, gain 0.445: fails the amended criterion | [section 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1119), MEASURED |
| The foliation residual on the saved coupled fields | negligible early (0.35 against 70 1/m), 212 to 411 against 1100 at blow-up; about 21 % of the growth | [section 14.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1155), MEASURED |
| Band renormalisation | t_blow 0.0668 to 0.0195 s (3.4x worse); the estimator error 0.3 to 0.6 % exceeds the 0.1 % drift it corrects | [section 15](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1211), MEASURED; [[concepts/curvature-corrugation-and-the-fit]] |
| The capillary-dt sweep and the exonerations | r = r0 + c dt with the dt term 90 % at N=128; adaptivity, GAMG, decomposition, the force time lag (psiOuterCorrectors, r 202.8 against 198.9) and the foot kernel (order 3.00) exonerated | [section 16](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1300), MEASURED; [[concepts/psi-outer-correctors]], [[concepts/departure-foot-ab2-centring]] |
| The production delivery is second order where the foliation is parallel | 0.278 (order 1.98) on the SDF ellipse; 8.69 (0.94) on the quadratic-form ellipse; `footPointEvaluatedFace` exact there (4.8e-4) and first order on true distances | [section 17](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1393), MEASURED |
| The m=2 mode and the explicit coupling | r_2 = 18.8, 12.7, 8.01, 4.02 1/s as dt halves three times; Richardson intercept +0.03 1/s; drift raises c by about 55 % | [section 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1554), MEASURED; [[concepts/parasitic-current-mechanism]] |

## Decisions

- The production delivery of August: the per-face stabilised foot-point inversion (order 1.98, gain 0.647) with `offsetCorrection none` (sections 12, 14); superseded in September by `cellCentreInverse` ([[decisions/curvature-extension-cell-centre-inverse]]), which the plan does not cover.
- The acceptance criterion for any delivery: G h^2 at most 0.65 and a fitted L2 order of at least 1.9 on the varying-curvature ellipse gate ([section 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1044)); an arm is never promoted on t_blow at one resolution ([section 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064)).
- The binding rule v0.3 (2026-08-07): psi advection is never modified or specialised for the static case; promotion needs the kinematic suite with the candidate active and a moving coupled gate ([section 2, WP6](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L323)).
- The psi filter is a delay device, not a closure (WP3, [section 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L90-L97)) ([[decisions/psi-filter-none]]).
- The Gaussian-curvature inverse in 3D (WP4, done) ([[decisions/curvature-inverse-gaussian]]).
- The semi-implicit capillary force becomes the principal licensed lever after section 18 ([[models/semi-implicit-capillary-force]]); its later history is in [[studies/shannon-parasitic-currents-campaign]].

## Retracted or superseded inside it

1. The WP8.1 consequence (build an element-DOF delivery) was falsified by WP8.2 within a day ([section 10, line 716](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L716)).
2. The cell-mean delivery, adopted in section 11.5 on coupled survival, was retired by the ellipse gate ([section 12, line 1040](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1040)) ([[retractions/cell-mean-delivery-adoption]]).
3. Section 16.1's "r0 does not vanish and rises toward fine grids" is superseded by the mode-level intercept of zero in section 18.4 ([line 1571](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1571)).
4. The first write-up of section 17 mislabelled the delivery rows; corrected 2026-08-14 in place. The "no stability penalty" note for `footPointEvaluatedFace` (17.4) predates its coupled falsification (blows earlier, volume 5 to 15x worse; [[concepts/face-curvature-deliveries]]).
5. The transport ground truth of section 0 is under the contamination notice of 2026-08-26; only the 2D vortex row is re-established (2.84 / 3.30) ([[retractions/gradu-coupled-patch-contamination]]).
6. The blow-up times of section 0 are no longer a score ([[retractions/t-blow-baseline]]); the plan itself says t_blow is not a valid proxy ([section 16.1, line 1311](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1311-L1313)).
7. Every filtered result predates the psi-filter seam bug found 2026-08-19 ([[retractions/psi-filter-seam-bug]]); every np 8 verdict before the seam-synchronisation fix of 2026-08-07 was discarded ([section 0, line 100](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L100-L110)).
8. The sign of the foliation residual was first written wrong (+Laplace of d); corrected in section 11.2 ([line 826](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L826-L831)).

## What it does not cover

The translating droplet (the plan is stationary-droplet only), the September deliveries `cellCentreInverse` and the K-aware inverse under coupling ([[concepts/cell-centre-inverse-curvature]]), the two-factor law and the order parameter g ([[studies/shannon-parasitic-currents-campaign]]), and the variational force ([[concepts/variational-capillary-force]]).

## Related

[[hubs/surface-tension]], [[hubs/verification]]. [[models/curvature-extension]], [[models/surface-tension-force]], [[models/semi-implicit-capillary-force]], [[concepts/parasitic-current-mechanism]], [[concepts/face-curvature-deliveries]], [[concepts/curvature-corrugation-and-the-fit]], [[concepts/balanced-force-csf-flux]], [[concepts/capillary-time-step]], [[concepts/psi-outer-correctors]], [[cases/stationary-droplet]], [[cases/curvature-static-gates]]. Siblings: [[studies/shannon-parasitic-currents-campaign]], [[studies/poly3d-roadmap]], [[studies/method-comparison]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from the plan document at 8867581.
