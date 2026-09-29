---
title: "Filtered results predate the psi-filter seam bug (2026-08-19)"
description: "VOIDED 2026-08-19 - every coupled result with psiFilter biharmonicBand on more than one rank: the filter was decomposition-dependent (an uncoupled L(psi) on processor patches and a band dilation over internal faces only), 53 to 205x off the unfiltered control; the re-runs changed every filtered number by 84 to 98 %"
aliases: [psi filter seam bug, biharmonicBand seam defect, f83a1ab]
kind: retraction
status: voided
part: surface-tension
tags: [retraction, part/surface-tension]
date: 2026-09-28
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/psiFilterEqn.H, config/seamConsistency3Dserial.yaml, config/seamConsistency3Dpar4.yaml, config/capillaryEnvelope.yaml, config/stationaryDroplet3Dwide.yaml, config/stationaryDroplet3DwideNoK.yaml]
sources: [STATUS 4 (2026-08-18 and 2026-08-19), PSH 0c, PCS 0 and WP3, CLAUDE no-filtering rule, commit f83a1ab]
---
# Filtered results predate the psi-filter seam bug (2026-08-19)

> VOIDED 2026-08-19. The claims were the coupled results with `psiFilter biharmonicBand` on more than one rank: "THE COMBINATION WORKS", `cellCentreInverse` plus theta = 0.2 converging in current, volume and shape at N = 64 / 128 / 256 (max|U| 2.24e-3 / 4.15e-4 / 1.02e-4, orders +2.43, +2.02) ([STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L470-L497)); the flat transport-order axis read in filtered arms ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L503-L518)); the 2D wide ladder to N = 512 (order 1.01 at the fourth point) and the 3D wide ladder that "destabilises at R/h ~ 16" with r ~ 79 1/s ([STATUS 4, 2026-08-19](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L775-L844)); and the K-exoneration ladder, which ran filtered too ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L846-L890), [consequence](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922)). The measurement: `psiFilterEqn.H` formed `Lpsi = psi - fvc::average(fvc::interpolate(psi))`, which inherits `calculated` patch types, so `L(psi)` was uncoupled on processor patches, and its band dilation looped internal faces only, so the filtered cell set depended on the decomposition ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L905)). Serial against np 4 on the 3D stationary droplet (N_L = 60, L = 6R, 20 steps, filter off as the control): max|U| differed by 2.43e-04 before the fix and 1.72e-06 after it, against 1.55e-06 with the filter off; the band amplitude A2h by 8.38e-07 before and 3.72e-10 after; the volume error by 1.92e-05 before and 6.55e-07 after, so the filtered solver was 53 to 205x further from its serial twin than the unfiltered one ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L907-L916)). Fixed in commit [f83a1ab](https://github.com/leia-openfoam/leia/commit/f83a1ab). Scope of the void: the 2D ladder (np 8, 4714 to 106689 steps) and the 3D ladder (np 32, 9207 to 18344 steps), and every convergence order, onset and verdict read from them ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922)). The re-runs on the fixed binary changed every filtered number by 84 to 98 % ([plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L172-L175)).

## The claim, and where it lived

- [STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L470-L518), `STATUS.md`: "THE COMBINATION WORKS (28 coupled arms)", the theta = 0.05 arm, the 3D N = 64 result and the flat transport-order axis. Marked SUPERSEDED 2026-09-28 ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L475-L478), commit d1e3414).
- [STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L520-L550): "THE CELL-CENTRE INVERSE AND THE FILTER"; the filtered half is void, the unfiltered `cellCentreInverse` arms are not. Marked SUPERSEDED 2026-09-28 ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L529-L533)).
- [STATUS 4, 2026-08-19](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L775-L844): "The wide ladders: 2D completes, 3D destabilises at R/h ~ 16". Marked SUPERSEDED 2026-09-28 ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L791-L794)).
- [STATUS 4, 2026-08-19](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L846-L890): "K is EXONERATED", the same ladder with the Gaussian term off. Its K-verdict rests on filtered runs ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L919-L922)); no marker on the subsection.
- The curated tables of the method-comparison theme: `capillary_envelope_coupled.csv`, `cell_centre_inverse_coupled.csv` and `wide_ladder_coupled.csv` ([the tables folder](https://github.com/leia-openfoam/leia/tree/8867581/docs/method-comparison/method-comparison-article/data/tables)), curated before the fix.
- [plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L90-L97) and [WP3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L235-L252), `docs/plan-curvature-stabilization.md`: the filter reframed as a delay device (blow-up at t = 0.167 s post seam-fix at N = 256, a 5x fuse extension).

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The filtered solver is decomposition-invariant. | Serial against np 4, filter on, 20 steps: max\|U\| 2.43e-04 before the fix, 1.72e-06 after; A2hL2Band 8.38e-07 before, 3.72e-10 after; volumeRelError 1.92e-05 before, 6.55e-07 after. The filter-off control: 1.55e-06, 4.09e-09, 3.64e-07, bit-unchanged by the fix. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L907-L916) |
| Why: `L(psi)` had no halo. | `psi - fvc::average(...)` inherits the patch types of the second operand; `fvc::average` carries `calculated` patch fields, so the second application of L read a stored face value, not the neighbour cell's `L(psi)`. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L896-L901), [commit f83a1ab](https://github.com/leia-openfoam/leia/commit/f83a1ab) |
| Why: the band never crossed a seam. | The dilation looped `mesh.owner()` and `mesh.neighbour()`, internal faces only. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L902-L905) |
| The 2D combination converges in every metric. | Post-fix max\|U\| 2.2433e-05, 2.8073e-05, 8.8277e-05 at N = 64 / 128 / 256: the current grows, orders -0.32 and -1.65. The pre-fix values were 99.0, 93.2 and 13.5 % too high, because a finer mesh at fixed np 8 puts less interface on seams. Shape and volume still converge (+3.41 / +1.38, +5.00 / +4.09). | MEASURED, [plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L119-L138) |
| The 3D ladder blows up at R/h = 15.8. | Post-fix the finest arm ends at 5.82e-03, not 7.04e-02, with min\|grad psi\| 0.985, not 0.859: a mild amplification, not a blow-up. Every filtered number moved 84 to 98 %. | MEASURED, [plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L177-L196) |
| The per-step gain is monotone in resolution. | Post-fix gAvg changes sign: 2D -7.61e-04, -6.41e-05, +7.22e-05 at N = 64 / 128 / 256; 3D -1.91e-04, -8.35e-05, +2.70e-04 at R/h = 10.0 / 12.7 / 15.8. A stability boundary, not a trend. | MEASURED, [plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L197-L201) |

## What survives

1. The unfiltered `cellCentreInverse` arms of 2026-08-18 are outside this void ([STATUS marker](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L529-L533)). Every np > 1 coupled run before 2026-09-27 still carries the parallel defects of [[concepts/coupled-face-density-defect]].
2. The location of the 3D onset, between R/h = 12.7 and 15.8, survives; its severity does not ([plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L193-L196)).
3. The filter is not the source: with the filter removed the two fine 3D arms still grow, A = +7.68e-04 and +1.16e-03 per step ([STATUS 4, 2026-08-20](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L986-L989)), see [[concepts/parasitic-current-mechanism]].
4. The seam check itself: `config/seamConsistency3D{serial,par4}.yaml` is the committed one-command regression, and the class of defect is on the list of [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L216-L223), see [[concepts/seam-checks-and-decomposition-invariance]].
5. The rule: no filtering in the production method; a filter's benefit measures the size of the defect (5.86x better at R/h = 15.8, 1.61x worse at R/h = 10.0) ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L382-L398)), see [[decisions/psi-filter-none]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: the [INVALIDATION](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L922); the three SUPERSEDED 2026-09-28 markers at [2026-08-18](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L475-L478), [the cell-centre inverse and the filter](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L529-L533) and [the wide ladders](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L791-L794) (commit d1e3414).
- [x] The plan documents: [plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L117-L206) (the post-fix re-runs and what they retract); [plan-curvature-stabilization WP3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L235-L252) (the filter as a delay device).
- [x] `CLAUDE.md`: the [4-rank rule's list](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L216-L223) and ["No filtering in the production method"](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L382-L398).
- [x] The regression configs `config/seamConsistency3D{serial,par4}.yaml` (commit f83a1ab).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The curated tables `capillary_envelope_coupled.csv`, `cell_centre_inverse_coupled.csv` and `wide_ladder_coupled.csv` still hold the pre-fix numbers without a marker.
- [ ] The [K-exoneration subsection](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L846-L890) has no marker; whether the `NoK` ladder was re-run after the fix is not recorded.
- [ ] The post-fix 2D N = 512 rung was "still running" ([plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L140-L143)); the record has no entry of its landing.

## Related

Hubs: [[hubs/surface-tension]], [[hubs/verification]]. Siblings: [[decisions/psi-filter-none]], [[concepts/curvature-corrugation-and-the-fit]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/cell-centre-inverse-curvature]], [[concepts/parasitic-current-mechanism]], [[concepts/coupled-face-density-defect]], [[models/narrow-band]], [[retractions/gradu-coupled-patch-contamination]], [[retractions/t-blow-baseline]], [[retractions/cell-mean-delivery-adoption]].

## Log

### 2026-09-28
Written from STATUS 4 (2026-08-18 and 2026-08-19), plan-shannon 0c and the commit f83a1ab. Voided 2026-08-19. Entered in [[retraction-log#2026-08]].
