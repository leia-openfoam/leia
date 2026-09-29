---
title: "Face curvature deliveries and their gain"
description: "The face-level curvature deliveries of the stabilization campaign; the per-face inverse is second order on the ellipse gate with gain 0.647, every construction that lowered the gain lost an order, no delivery changed the blow-up exponent, and the delivery lever is closed (2026-09-28)."
aliases: [face deliveries, curvature noise gain, G h^2, stabilizedFootPointFace, cut-cell delivery, cell-mean delivery]
kind: concept
status: retracted
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled:
decided_by:
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/capillaryDriverSplit.H, applications/test/leiaTestCurvatureNoiseGain/leiaTestCurvatureNoiseGain.C, applications/test/leiaTestMeanCurvature/leiaTestMeanCurvature.C]
sources: [RM face gates, PCS 9 to 15, PCS 17, STATUS 1, STATUS 4, STATUS 7 acceptance criterion, SL deck 4/19 and 4/22]
---
# Face curvature deliveries and their gain

> Verdict (2026-09-28). A face delivery re-references the interpolated face curvature to the interface through the stabilised foot point of the face centre and the parallel-surface inverse. It restores second order on the circle (`h^2.04`, 0.105 1/m at N = 512, 108 times better) and on the sphere (`h^1.95`) ([RM face gate](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1310-L1313), [RM sphere](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1378-L1381)). The campaign then built variants to lower the curvature noise gain `G h^2`, and measured that accuracy and gain are dissociated: the cut-cell delivery is the most accurate static delivery and blows up 3.3 times sooner ([PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L702-L729)); the cell-mean delivery survives longest but collapses to first order on the varying-curvature ellipse (1.03 against 1.98) ([PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1004-L1018)); the symmetric face mean loses the order too (1.10) ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1127-L1142)); and no delivery changed the exponent of `t_blow(N)` ([PCS 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1070-L1093)). The acceptance criterion is `G h^2 <= 0.65` and a fitted order `>= 1.9` on the ellipse gate; only the per-face inverse passes ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2685-L2688)). The delivery lever is closed ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1147-L1153)), the cell-mean adoption is retracted ([[retractions/cell-mean-delivery-adoption]]), and production moved to the corrected cell field on 2026-09-01 ([[concepts/cell-centre-inverse-curvature]]).

## What it is

The CSF flux applies the face curvature `kappa_f`, so the static gates measure `kappa_f` on the active faces, assembled exactly as the force model does ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1299-L1308)). The face words of [[models/curvature-extension]] fill the registered field `kappaStableFootFace`, consumed with `faceCurvatureSource registered`, and they require `offsetCorrection none` ([PCS 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L107-L113)). The gain is `G = ||d kappa_f|| / ||d psi||` on the active faces for a pseudo-random perturbation of `psi` of amplitude `eps h`, with alpha and the active set held fixed; `G h^2` is the amplification relative to a plain second difference ([`leiaTestCurvatureNoiseGain.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestCurvatureNoiseGain/leiaTestCurvatureNoiseGain.C#L11-L57)). It is constant to about 3 percent over N = 64 to 256, linear in `eps` over a factor of 100, and unchanged by the interface shape ([PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L743-L745), [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1030-L1033)).

| delivery | circle order, L2 at N = 512 | ellipse order, L2 at N = 512 | G h^2 | t_blow at N = 128 / 256 [s] |
|---|---|---|---|---|
| arithmetic (no extension) | 1.13, 11.35 | 0.97, 13.66 | 0.647 | 0.0803 / not recorded |
| per-face inverse `stabilizedFootPointFace` | 2.04, 0.105 | 1.98, 0.2785 | 0.647 | 0.0668 / 0.0348 |
| cut-cell inverse `cutCellFootPointFace` | 2.00, 0.0761 | 1.02, 5.848 | 0.818 | 0.0202 / 0.0145 |
| cell-mean inverse `cellMeanFootPointFace` | 2.05, 0.08422 | 1.03, 5.851 | 0.395 | 0.1049 / 0.0493 |
| symmetric face mean, theta 0.5 / 1.0 | 1.97 to 2.00 | 1.10, 2.014 / 1.07, 4.009 | 0.445 / not recorded | not run |
| foot-evaluated `footPointEvaluatedFace` | about 1.0, 10.5 | about 1.0, 12.7 (SDF); exact 4.8e-4 (quadratic form) | equal to production within 3 percent | earlier than production at every N on the oscillating case |
| interface mean (control) | 2.03, 0.15 | -0.01, 505.6 | not recorded | diagnostic only |

Sources of the rows: [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L391-L400), [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1004-L1012), [PCS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L907-L924), [PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1127-L1135), [PCS 17.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1455-L1471), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L434-L448). The t_blow values are comparisons inside one campaign; the historical baseline they were read against is retracted ([[retractions/t-blow-baseline]]).

## Why it matters

Three hypotheses about the parasitic current were tested with these deliveries ([STATUS 1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L310-L334)). The across-support variation of `kappa_f` dominates the along-interface variation four to six times for both deliveries, and the foot-point advantage (17.0 against 967 1/m at t = 0) is gone within 0.05 s ([PCS 9](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L616-L642)). Removing that variation with the cut-cell delivery made the instability faster, so the variation is not the driver ([PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L716-L723)). The gain orders the extremes: equal gain gives equal growth rate despite a 27-fold accuracy difference, and 28 percent more gain gives 2.4 times the rate; but it does not order every pair, since arithmetic and per-face have the same gain and arithmetic outlives per-face ([PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L737-L751), [PCS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L967-L970)). The static gates are seconds per mesh, so a candidate delivery can be rejected before any coupled run ([[cases/curvature-static-gates]]).

## Where in the code

- All six face words are in [`stabilizedFootPointFaceCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H); the dispatch is in [`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L97-L104).
- The across/along split of the delivered curvature: [`capillaryDriverSplit.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/capillaryDriverSplit.H), per step in the droplet metrics ([PCS 9](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L606-L614)).
- The face section of [`leiaTestMeanCurvature.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestMeanCurvature/leiaTestMeanCurvature.C#L29-L51) and the gain app [`leiaTestCurvatureNoiseGain.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestCurvatureNoiseGain/leiaTestCurvatureNoiseGain.C).
- The processor-seam rule: a face value from rank-local fits is synchronised by swap-and-average; the one-sided correction made np 8 runs 7 to 10 times noisier than np 4 ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1397-L1409), [[concepts/seam-checks-and-decomposition-invariance]]).
- The Eulerian two-phase solver still delivers `stabilizedFootPointFace` ([STATUS 11.16](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3704-L3707)).

## Evidence

| claim | number | where |
|---|---|---|
| the per-face inverse restores second order on the circle | h^2.04, 11.35 to 0.105 1/m at N = 512; per-side variant h^1.40 | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1310-L1334), MEASURED |
| the K-aware inverse on the sphere | h^1.95, 5.25e-3 1/m at N = 128; per-side h^1.98 | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1378-L1384), MEASURED |
| the driver split at N = 128 | across / along 967 / 179 (arithmetic) and 17.0 / 12.8 (foot point) at t = 0; foot point 2969 / 730 at t = 0.06 | [PCS 9](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L616-L626), MEASURED |
| the cut-cell delivery blows up sooner | t_blow 0.0202 against 0.0668 s (N = 128), 0.0145 against 0.0348 (N = 256); late rate 454 against 189 1/s | [PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L702-L714), MEASURED |
| the gain orders the rate | G h^2 0.644 / 0.643 / 0.821 gives 186 / 189 / 454 1/s | [PCS 10](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L737-L751), MEASURED |
| the cell-mean delivery survives longest on the circle | 0.1049 s at N = 128, rate 59 to 100 1/s from t = 0.05 to 0.10 | [PCS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L926-L944), MEASURED |
| the ellipse gate collapses the lumped deliveries | orders 0.97 / 1.98 / 1.02 / 1.03 / -0.01; the lumped deliveries are 21 times less accurate at N = 512 | [PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1004-L1018), MEASURED |
| no delivery changed the exponent | t_blow ~ N^p with p = -0.94 / -0.48 / -1.09; the cell-mean advantage falls 1.57 to 1.04 by N = 2048 | [PCS 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1070-L1093), MEASURED and extrapolated |
| the symmetric face mean loses the order | 1.10 (theta 0.5) and 1.07 (theta 1.0); gain 0.445 | [PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1127-L1142), MEASURED |
| the foot-evaluated delivery is falsified coupled | oscillating t_blow ratio 0.36 / 0.42 / 0.32 at N = 64 / 128 / 256; volume 5 to 15 times worse | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L434-L453), MEASURED |
| production is second order only where the foliation is parallel | 0.278 at 1.98 (SDF ellipse) against 8.69 at 0.94 (quadratic form); `quadraticNewtonFoot` inverts to 0.187 at 1.95 | [PCS 17](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1405-L1442), MEASURED |

## Why it failed, or why we think so

Assigning one value to every active face of a cut cell centres each face's value on the cell, an `O(h dkappa/ds)` offset wherever the curvature varies: first order by construction, and better averaging cannot repair the lumping ([PCS 12](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1020-L1028)). A symmetric mean about the face fails for another reason: the active-face mask depends on where the interface cuts the cell, so the stencil is asymmetric at O(1) ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1137-L1142)). The pattern across four levers is one statement: the system is not limited by systematic error in any operator but by its sensitivity to perturbations of `psi`; every corrective built from the same fit injects the fit's error, and accuracy improvements are free because they do not change the gain ([PCS 15.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1260-L1277)).

## Decisions

- Production delivery until c935883: `stabilizedFootPointFace` (order 1.98, `G h^2` 0.647); retired: `cutCellFootPointFace` and `cellMeanFootPointFace`, kept in the code and the gates as the two extreme points ([STATUS 1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L336-L347)).
- The acceptance criterion for any future delivery: `G h^2 <= 0.65` and order `>= 1.9` on the ellipse gate, not on the circle ([STATUS 7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2685-L2688)).
- No arm is promoted on `t_blow` at a single resolution; two resolutions are the minimum, three are needed for a claim ([PCS 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1110-L1115)).

## Open questions

1. `cellFootPointEvaluatedFace` has a recorded prediction and no result ([[models/curvature-extension]]).
2. A per-face selection or blend between the production inverse and the foot-evaluated value, both from the same fit ([PCS 17.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1478-L1483)).
3. Smoothing a cell field, whose stencil does not depend on the interface position, was recorded as low value ([PCS 14.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1150-L1153)); the cell-centre inverse is the cell-field construction that was adopted instead.

## Related

[[hubs/surface-tension]], [[models/curvature-extension]], [[concepts/cell-centre-inverse-curvature]], [[concepts/curvature-from-the-fit]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[concepts/curvature-corrugation-and-the-fit]], [[concepts/seam-checks-and-decomposition-invariance]], [[cases/curvature-static-gates]], [[cases/stationary-droplet]], [[cases/oscillating-droplet]], [[retractions/cell-mean-delivery-adoption]], [[retractions/t-blow-baseline]], [[studies/curvature-stabilization-campaign]].

## Log

### 2026-09-28
Created from the roadmap's face gates, PCS sections 9 to 17 and STATUS sections 1, 4 and 7.
