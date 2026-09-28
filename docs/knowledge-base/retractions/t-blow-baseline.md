---
title: "The t_blow baseline as a proxy (2026-08-31)"
description: "RETRACTED 2026-08-31 - the recorded blow-up times (0.100/0.078/0.036 s for cellCentreInverse at N = 64/128/256) do not reproduce on the current code, and t_blow is not a proxy for the growth rate; every later ladder carries its own matched baseline"
aliases: [t_blow retraction, blow-up time baseline]
kind: retraction
status: retracted
part: surface-tension
tags: [retraction, part/surface-tension]
date: 2026-09-28
code: [config/fullHorizonStability2D.yaml, config/stationaryDropletDtSweep.yaml, config/projFluxStationary2D.yaml]
sources: [STATUS 4 (2026-08-31), STATUS 4 (2026-08-18), PCS 13 and 16.1, SL deck 3/13 and 4/22]
---
# The t_blow baseline as a proxy (2026-08-31)

> RETRACTED 2026-08-31. The claim was that the recorded blow-up times form a baseline: `cellCentreInverse` 0.100 / 0.078 / 0.036 s against production 0.100 / 0.063 / 0.033 s at N = 64 / 128 / 256 in 2D, and a new arm is scored against them ([STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L529-L535)). The measurement: the reproduction control of `fullHorizonStability2D` (off + cellCentred, N = 128, filters off, 13 334 steps) reached t = 0.1 s alive with max|U| = 4.65e-3 m/s and a fitted rate of +118 1/s. The extrapolated blow-up is near t = 0.146 s, an arrival-time shift of about 1.9x, outside the 5 to 38 % seed sensitivity that the campaign had measured for t_blow ([STATUS 4, 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1057-L1068)). Before that, the capillary-dt sweep had found that t_blow is not a proxy for the growth rate: the e-fold count K = r t_blow runs from 5.4 to 13.3 across arms ([plan 16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1311-L1313)). Scope: no result may be compared against the older t_blow table. A comparison within one study (one commit, one set of binaries, one mesh, one dt law) stands ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1065-L1068)). Every later ladder carries its own matched baseline ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1166-L1175)). No data is void.

## The claim, and where it lived

- [STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L529-L535), `STATUS.md`: the t_blow rows of `cell_centre_inverse_coupled.csv`, without a marker.
- [STATUS 4, 2026-08-10](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L685-L716) and [plan 11.5](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L914-L944): the cell-mean delivery adopted on t_blow 0.1049 s, see [[retractions/cell-mean-delivery-adoption]].
- [plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1064-L1093): t_blow exponents p fitted from two resolutions.
- The quadratic SL deck, slides [3/13](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/3/13) (blow-up times 0.44 / 0.105 / 0.035 s for the quadratic pipeline) and [4/22](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/22) (the t_blow column of the delivery table).
- [STATUS 4, 2026-08-15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L439-L459): the footPointEvaluated verdict rests on t_blow, within one study at three resolutions; that use is inside the surviving scope.
- The SL article, [sec:droplet](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1720-L1721) `docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`: one ladder, t_blow 0.47 / 0.44 / 0.11 / 0.03 s at N = 32 to 256.

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The 0.078 s blow-up time at N = 128 reproduces. | The control is alive at t = 0.1 s; max|U| 4.65e-3 m/s; rate +118 1/s over the last fifth; extrapolated blow-up near 0.146 s (1.9x). | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1057-L1064) |
| t_blow measures the growth rate. | K = r t_blow between 5.4 and 13.3 across the dt-sweep arms; incubation varies. | MEASURED, [plan 16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1311-L1313) |
| t_blow at one resolution ranks deliveries. | The cell-mean prefactor win vanishes under refinement: p = -1.09 against -0.94 for the per-face inverse. | MEASURED, [plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1070-L1093) |
| A reproducible step count is an instability. | 9337 / 9331 / 9367 / 9324 on four occasions, against 5 to 38 % scatter: a geometric event. | MEASURED, [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L528-L536) |
| The cause of the 1.9x shift. | Not found. "Until the cause is found, no result may be compared against the older t_blow table." | HYPOTHESIS, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1063-L1066) |

## What survives

1. The instability itself, qualitatively: the control's rate is +118 1/s and accelerating ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1060-L1062)).
2. The verdict of `fullHorizonStability2D` within its own matrix: off + projectedFlux is the best arm on every metric, rate -52.0 1/s ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1087)), see [[concepts/trace-velocity-projected-flux]].
3. The scoring on the growth rate: r(A2h) is the order parameter ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L466-L468), [STATUS 5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2555-L2559)), and G = A T/dt is the e-fold count at fixed physical time ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L991-L1005)).
4. The scoring rule: two resolutions are the minimum for an exponent, three for a claim ([plan 13](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1110-L1115)).

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: [RETRACTION FIRST](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1050-L1068); the refinement ladders carry [their own baselines](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1166-L1175); the [scoring instruction](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2555-L2559).
- [x] The plan document: [16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1311-L1313) and the [binding scoring note](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1110-L1115).
- [x] `CLAUDE.md`: ["Is the interface still inside the domain?"](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L528-L544).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The cause of the 1.9x shift in arrival time: no record.
- [ ] [STATUS 4 L529-L535](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L529-L535) and the deck slides 3/13 and 4/22 show the old values without a marker.
- [ ] The SL article ladder ([L1720-L1721](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1720-L1721)) has no pointer to the rule that its values are not comparable across commits.

## Related

Hubs: [[hubs/surface-tension]], [[hubs/verification]]. Siblings: [[concepts/parasitic-current-mechanism]], [[concepts/trace-velocity-projected-flux]], [[concepts/cell-centre-inverse-curvature]], [[concepts/error-vector-and-read-out-instants]], [[cases/stationary-droplet]], [[retractions/cell-mean-delivery-adoption]], [[retractions/psi-filter-seam-bug]].

## Log

### 2026-09-28
Written from STATUS 4 (2026-08-31) and plan sections 13 and 16.1. Retracted 2026-08-31. Entered in [[retraction-log#2026-08]].
