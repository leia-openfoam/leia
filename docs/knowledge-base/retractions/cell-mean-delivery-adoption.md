---
title: "The adopted cellMean delivery, retracted by the varying-curvature ellipse gate (2026-08-12)"
description: "RETRACTED 2026-08-12 - the cell-mean curvature delivery was adopted on 2026-08-10 for the lowest gain (G h^2 0.402) and the longest coupled survival (t_blow 0.1049 s against 0.0668 s at N = 128); on the 2:1 ellipse gate it is first order (1.03 against 1.98) and 21x less accurate, and its survival advantage is a prefactor that vanishes under refinement (exponent -1.09 against -0.94)"
aliases: [cellMean retraction, cell-mean delivery, cellMeanFootPointFace]
kind: retraction
status: retracted
part: surface-tension
tags: [retraction, part/surface-tension]
date: 2026-09-28
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H, config/stationaryDropletCellMean.yaml, config/faceCurvatureEllipse2D.yaml, cases/ellipseDroplet2D]
sources: [PCS 11.4 and 11.5, PCS 12, PCS 13, PCS 14.1, STATUS 4 (2026-08-10 and the ellipse table), STATUS 7 item 3, commit 2a1be92]
---
# The adopted cellMean delivery, retracted by the varying-curvature ellipse gate (2026-08-12)

> RETRACTED 2026-08-12. The claim, adopted on 2026-08-10, was "the cell-mean delivery (`curvatureExtension cellMeanFootPointFace`: the per-face parallel-surface inversions averaged over each cut cell's active faces) is the first arm to move the coupled behaviour: it has the lowest gain of every delivery measured (G h^2 0.402 in 2D and 0.534 in 3D against 0.639 and 0.826 for the per-face inverse) at equal or better static accuracy, and it survives longest on the coupled stationary droplet (t_blow 0.1049 s against 0.0668 s for the per-face inverse at N = 128, 0.0493 s against 0.0348 s at N = 256)" ([plan 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L902-L944), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L675-L716), commit [2a1be92](https://github.com/leia-openfoam/leia/commit/2a1be92)). The measurement: the varying-curvature gate `config/faceCurvatureEllipse2D.yaml` on the 2:1 signed-distance ellipse (curvature 250 to 2000 1/m). The cell-mean delivery collapses to first order, 1.03 against 1.98 for the per-face inverse (order fitted on N >= 128), and its L2 error at N = 512 is 5.851 1/m against 0.2785, 21x worse, while both are second order on the circle (2.05 and 2.04) ([plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1004-L1018), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400)). The same day the refit of t_blow against N gave the cell-mean arm the worst exponent, -1.09 against -0.94 (per-face) and -0.48 (cut-cell): a prefactor win that decays from 1.57x at N = 128 to 1.04x at N = 2048 ([plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064-L1093)). Scope: the delivery is retired from consideration as a production delivery; the per-face inverse stays (order 1.98, G h^2 0.647) ([plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1035-L1047)). No data is void.

## The claim, and where it lived

- [plan 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L902-L989), `docs/plan-curvature-stabilization.md`: the gain table, the coupled t_blow table, the matched-window growth rates and the mechanism read from one run. [Plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1035-L1042) records the retraction under CONSEQUENCE.
- [STATUS 4, 2026-08-10](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L675-L716), `STATUS.md`: the same tables. The [ellipse table](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400) at the head of section 4 and [section 7 item 3](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2669-L2673) ("what retracted the cell-mean adoption") carry the retraction.
- [`config/stationaryDropletCellMean.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDropletCellMean.yaml#L1-L30): the pre-registered header ("if the gain is the governing quantity, this is the arm that should survive longest").
- [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L391): `cellMeanFootPointFace` in the `CURVATURE_EXTENSION` list of words.
- The quadratic SL deck, [slide 4/22](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/22): the delivery table with the cell-mean row.
- The curated table `face_curvature_orders_ellipse.csv` ([blob](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/face_curvature_orders_ellipse.csv)) holds the retracting measurement.

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| Static accuracy of the delivery is equal or better. | 2:1 ellipse, L2 of kappa_f on the active faces: cell-mean 43.44 / 24.27 / 11.81 / 5.851 at N = 64 to 512, order 1.03; per-face 14.17 / 4.323 / 1.116 / 0.2785, order 1.98; cut-cell 1.02; interface mean -0.01. | MEASURED, [plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1004-L1012) |
| Why it is structural. | One value per cut cell puts a value centred on the CELL onto every active face; where the curvature varies the offset is O(h dkappa/ds), first order by construction. A symmetric average about each FACE would keep O(h^2), a cell-centred assignment cannot. | DERIVED, [plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1020-L1028) |
| The circle gate could rank the deliveries. | On the circle every delivery is second order (1.97 to 2.05); its constant curvature hides the lumping defect exactly as the stationary droplet does. | MEASURED, [plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1014-L1018), [plan 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1144-L1145) |
| The coupled survival is a change of scaling. | t_blow ~ N^p: per-face -0.94, cut-cell -0.48, cell-mean -1.09. Extrapolated ratio cell-mean over per-face: 1.57 / 1.42 / 1.28 / 1.15 / 1.04 at N = 128 to 2048. | MEASURED, [plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1070-L1093) |
| The gain orders the deliveries. | The gain is unchanged by the geometry (cell-mean 0.370 to 0.416, per-face 0.612 to 0.673, cut-cell 0.779 to 0.849 across the ellipse ladder), yet the ranking by accuracy inverts. Arithmetic and per-face have the same gain to 0.2 % and different growth. | MEASURED, [plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1030-L1033), [plan 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L967-L970) |
| The remaining variant, a symmetric face mean, keeps second order. | Order 1.10 at theta = 0.5 (G h^2 0.445) and 1.07 at theta = 1.0: the active-face mask depends on where the interface cuts the cell, an O(1) asymmetry. The delivery lever is closed. | MEASURED, [plan 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1127-L1153) |

Two caveats on the coupled numbers of the adoption: the t_blow values are not comparable across commits ([[retractions/t-blow-baseline]]), and the coupled arms ran at np 8 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L685)) before the coupled-face density fix of 2026-09-27 ([[concepts/coupled-face-density-defect]]).

## What survives

1. The gain-versus-accuracy dissociation: the cut-cell and cell-mean deliveries stay as the two extreme points that established it ([plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1040-L1042)). The mechanism of the gain is averaging against concentrating, not the derivative count ([plan 11.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L886-L900)).
2. The amended acceptance criterion: G h^2 <= 0.65 AND a fitted order >= 1.9 ON THE VARYING-CURVATURE GATE ([plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1044-L1047)), see [[cases/curvature-static-gates]] and [[concepts/face-curvature-deliveries]].
3. The scoring rule: no arm is promoted on t_blow at a single resolution; two resolutions give an exponent, three make it a claim ([plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1110-L1115)).
4. The methodological rule: a static test case must not be simpler than the interfaces the method runs on; stability on a circle does not predict accuracy on anything else ([plan 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1056-L1062)).
5. The production delivery of that date, the per-face inverse `stabilizedFootPointFace` (order 1.98, G h^2 0.647), later replaced by `cellCentreInverse` (commit c935883), see [[decisions/curvature-extension-cell-centre-inverse]] and [[concepts/cell-centre-inverse-curvature]].

## Propagation (checklist, same commit)

Done:

- [x] The plan document: [12 CONSEQUENCE](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1035-L1047), [13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1088-L1093) and [14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1147-L1153).
- [x] `STATUS.md`: the [ellipse table](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400) at the head of section 4 and [section 7 item 3](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2669-L2673).
- [x] The curated tables `face_curvature_orders_ellipse.csv` and `face_curvature_orders.csv`.
- [x] The knowledge-base row in [[models/curvature-extension]] (status retracted).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] [STATUS 4, 2026-08-10](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L675-L716) shows the adoption tables without a marker; the reader meets the retraction only at the head of the section.
- [ ] The 3D companion gate (a torus: exact signed distance with non-constant mean curvature) and the psi-transform arm ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2669-L2673)) are not run.
- [ ] The cell-field smoothing variant, whose stencil does not depend on the interface position, is recorded as low expected value and unmeasured ([plan 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1150-L1153)).

## Related

Hubs: [[hubs/surface-tension]]. Siblings: [[models/curvature-extension]], [[concepts/face-curvature-deliveries]], [[cases/curvature-static-gates]], [[concepts/cell-centre-inverse-curvature]], [[concepts/parasitic-current-mechanism]], [[concepts/coupled-face-density-defect]], [[studies/curvature-stabilization-campaign]], [[decisions/curvature-extension-cell-centre-inverse]], [[retractions/t-blow-baseline]], [[retractions/psi-filter-seam-bug]].

## Log

### 2026-09-28
Written from plan-curvature-stabilization 11.5, 12, 13 and 14.1 and STATUS 4. Retracted 2026-08-12. Entered in [[retraction-log#2026-08]].
