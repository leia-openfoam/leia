---
title: "The CFL 1.0 row of the SL 2D convergence table: alpha read a stale narrow band (2026-10-01)"
description: "RETRACTED 2026-10-01 - the published 2D vortex orders at CFL 1.0 (shape 2.378, volume 3.542): the kinematic SL solver computed alpha with the narrow band of the previous step; with the fix the same ladder gives 2.465 and 3.193, the shape error at T is 10.5 to 23.3 % lower at five of seven rungs; the CFL 0.5 row stands (orders within 0.005)"
aliases: [stale narrow band, stale band, CFL 1.0 orders retracted]
kind: retraction
status: retracted
part: advection
tags: [retraction, part/advection]
date: 2026-10-01
code: [applications/solvers/leiaSemiLagrangeLevelSetFoam/leiaSemiLagrangeLevelSetFoam.C, applications/solvers/leiaSemiLagrangeLevelSetFoam/advectionErrors.H, src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.C, config/uncachedConv2Dvortex.yaml, workflow/scripts/make_convergence_table.py]
sources: ["STATUS 11.22 item 5.3", "STATUS 11.23 items 5 to 7", "sl_convergence_orders.csv", "uncachedConv2Dvortex_errors.csv"]
---
# The CFL 1.0 row of the SL 2D convergence table: alpha read a stale narrow band (2026-10-01)

> RETRACTED 2026-10-01. The claim: the production semi-Lagrangian transport has, on the reversed 2D vortex at CFL 1.0, the shape order 2.378 and the volume order 3.542 ([sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/d2984c5e/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv), row 2, study `uncachedConv2Dvortex`). The measurement: `leiaSemiLagrangeLevelSetFoam` computed alpha with the narrow band of psi^n, and both phase indicators give a sign-based 0/1 alpha outside the band ([STATUS 11.23 item 1](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L4958-L4964)). The pre-fix and the fixed solver on identical copies of the published ladder: psi is byte-identical, the alpha metrics are not. The pre-fix runs reproduce the published orders to within 0.008; the fixed runs give 2.465 and 3.193 ([STATUS 11.23 item 7](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5009-L5036)). The CFL 0.5 row stands: its fitted orders move by at most 0.005 and the finest-rung errors do not change. The fix is 6f63418a. The curated table is not regenerated yet: author decision.

## The claim, and where it lived

- [sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/d2984c5e/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv) row 2, and the generated `convergence_orders.tex` and `convergence_orders_extended.tex` that the SL article inputs ([SL article L1488-L1499](https://github.com/leia-openfoam/leia/blob/d2984c5e/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1488-L1499)).
- The rows at CFL 1.0 of [uncachedConv2Dvortex_errors.csv](https://github.com/leia-openfoam/leia/blob/d2984c5e/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/uncachedConv2Dvortex_errors.csv) (`shapeError`, `volumeError`).
- The knowledge base: [[models/level-set-advection]], [[models/sl-reconstruction]], [[decisions/sl-reconstruction-uncached-qwls]], [[studies/sl-quadratic-pre-print]] (CORRECTED there on 2026-10-01).
- The curvature plan, [plan-curvature-stabilization L55](https://github.com/leia-openfoam/leia/blob/d2984c5e/docs/plan-curvature-stabilization.md#L55) ("2.38 / 3.54 at CFL 1.0").

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The fix leaves the transport alone. | psi byte-identical in all 27 pairs; every column of `gradPsiError.csv` identical at every step | MEASURED, [STATUS 11.23 item 4](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L4989-L4993) |
| The pre-fix runs are the published runs. | fitted orders within 0.008 of the table; `volumeError` within 0.8 %; `shapeError` within 1 % at N = 32 | MEASURED, [STATUS 11.23 item 7](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5009-L5036) |
| CFL 1.0, shape order | 2.378 published, 2.374 pre, 2.465 post | MEASURED, same |
| CFL 1.0, volume order | 3.542 published, 3.540 pre, 3.193 post | MEASURED, same |
| CFL 1.0, shape error at T | 10.5 / 14.8 / 23.3 / 11.0 / 15.5 % lower at N = 64 / 90 / 128 / 181 / 256 | MEASURED, same |
| CFL 0.5, the three fitted orders | move by at most 0.005; the finest-rung errors are unchanged | MEASURED, same |
| Per step on the production 2D ladder (CFL 0.5) | the median relative change of the volume error is 0.52 / 1.82 / 4.51 / 6.22 % at N = 32 / 64 / 128 / 256 | MEASURED, [STATUS 11.23 item 6](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5001-L5008) |

Why CFL 1.0 more than 0.5 (DERIVED): at CFL 1 the interface crosses up to one cell per step, so a cell that enters the band is more often cut by the interface. Why T/2 less than T (DERIVED): near T/2 the factor cos(pi t/T) is near zero, so few cells enter the band in the last step.

## What still stands

The transport and every psi metric: the gradient errors, the band gradient errors and their orders. The CFL 0.5 row. The production polyhedral rung and the published polyhedral configuration did not change at the read-out instants ([STATUS 11.23 item 7](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5009-L5036)).

## Not yet measured

18 more curated tables carry rows of this solver. Eight of the 19 have arms at CFL 0.8 or 1.0: `uncachedConv2Dvortex`, `npslConv2Dvortex`, `nslConv2Dvortex`, `sdCompare2D`, `linearConv2Dvortex`, `linearConv2DvortexClip`, `linearConv3Dshear`, `kinematicTranslation2D` ([STATUS 11.23 item 7](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5009-L5036)). Their alpha metrics at T are not re-measured; their psi metrics stand.

## Propagation

STATUS 11.22 item 5.3 and 11.23; the four knowledge-base notes above (CORRECTED lines); the curvature plan at L55 (a dated note); the SL hand-over, task T4 ([[sessions/sl-session-handover]]). The curated tables, the generated order tables and the SL deck wait for the regeneration.

## Related

[[decisions/kinematic-solver-per-flow-solver]] - [[models/level-set-advection]] - [[models/phase-indicator]] - [[models/narrow-band]] - [[studies/sl-quadratic-pre-print]] - [[retractions/gradu-coupled-patch-contamination]] - [[retraction-log]]

## Log

### 2026-10-01
RETRACTED on the gate of STATUS 11.23. Entered in [[retraction-log#2026-10]].
