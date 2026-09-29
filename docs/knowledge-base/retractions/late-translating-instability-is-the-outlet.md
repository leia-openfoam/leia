---
title: "\"The late translating instability is the outlet\" (corrected 2026-09-27: an interior growth remains)"
description: "CORRECTED 2026-09-27 - the reading that the outlet alone causes the late divergence of the translating droplet: the 20 mm box also diverges, at t = 0.2127 s (N = 100) and 0.2194 s (N = 142) with the droplet more than 8 mm from the outlet, and the 40 mm box degrades without an outlet nearby; the outlet triggers the fast phase, a slower interior growth of about 30 1/s remains and shrinks with h between two rungs"
aliases: [outlet trigger correction, late translating instability, interior growth of the translating droplet]
kind: retraction
status: retracted
part: mass-flux
tags: [retraction, part/mass-flux]
date: 2026-09-28
code: [cases/translatingDroplet2D, config/gates/methodGate2D.yaml, config/translatingRepaired2D.yaml]
sources: [STATUS 11.14, STATUS 11.15, gcls pre-print sec:res-translating, SL article sec:translating-late, METHOD 6 and 8.1 row CURVATURE_EXTENSION, CLAUDE step 5]
---
# "The late translating instability is the outlet" (corrected 2026-09-27: an interior growth remains)

> CORRECTED 2026-09-27 at 03:40, the same night. The claim was the first heading of the box study in STATUS 11.15, "the late translating instability is the outlet" ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3769-L3771), corrected in place). Its basis: the baseline translating case at N = 100 in a box twice as long (20 mm, the same h, the outlet 10 mm further away) completes t = 0.1 s and its spurious current decays, where the 10 mm box grows from 1.79e-3 to 1.35e-1 between t = 0.06 and 0.08 s ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3773-L3782)). The measurements that corrected it, all with the fixed binaries: at N = 142 a slower growth remains in the 20 mm box after t = 0.07 s, 4.3e-4 to 1.8e-3 in 0.03 s (about 46 1/s), with the droplet more than 11 mm from the outlet ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3795-L3798)); the 20 mm box run to 0.25 s DIVERGES at step 19590, t = 0.2127 s, at N = 100, with a growth of about 30 1/s from t = 0.10 s that turns explosive after 0.16 s with the droplet more than 8 mm from the outlet, and at step 34187, t = 0.2194 s, at N = 142 ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3816-L3829)); the 40 mm box, with the outlet more than 17 mm away, completes 0.3 s but degrades: L2|U-U0| 3.6e-4 to 8.3e-3, a 10 % volume drift and a 3.9 mm lead of the droplet over the stream at N = 100 ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3917-L3930)). The corrected reading: a slow interior growth in time, amplified as the droplet approaches the outlet; the outlet triggers the fast phase ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3860-L3863)). Scope: the mechanism, not the data. Every run stands; the gate horizon of 0.05 s ends before both phases.

## The claim, and where it lived

- [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3769-L3784), `STATUS.md`: the heading, corrected in place with "I was wrong", and the 20 mm box table.
- [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3743): the four one-change discriminators (serial 0.0775 s, np 4 0.0868 s; `curvatureExtension none` 20 % earlier; `footIntegrator rk2` inside the scatter; `capillaryForceCentring midpoint` 32 % earlier; density ratio 1 completes) and the sentence "the divergence onset is when the leading edge is 3.5 mm (35 cells) from the outlet".
- [`config/gates/methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L137-L143): the translating arm's `END_TIME 0.05` rationale.
- The gcls pre-print, `docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex`, [sec:res-translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L789-L823) `sec:res-translating`: written after the correction, with the five observations (item 4 "the outlet triggers it", item 5 "a slower growth remains in the interior"). Moved on 2026-09-28 to the SL article ([d1e3414](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2278-L2315), `sec:translating-late`); the gcls article keeps a pointer ([d1e3414](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L720-L728)).

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| A longer box removes the divergence. | The 20 mm box to 0.25 s: N = 100 DIVERGED at step 19590 (t = 0.2127 s); N = 142 DIVERGED at step 34187 (t = 0.2194 s). L2\|U-U0\| 5.1e-4 / 7.7e-4 / 1.4e-3 / 3.4e-3 / 8.5e-2 / 5.9e-1 at t = 0.10 to 0.21 s; the centroid runs ahead of U0 t after 0.16 s. | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3816-L3829) |
| The growth needs the outlet. | The 40 mm box (outlet more than 17 mm away) completes 0.3 s at N = 100 and 142 and degrades: current from t = 0.10 s at about 26 1/s at first and 5 1/s at the end; volume drift 10 %; lead 3.9 mm. | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3917-L3930) |
| The interior growth is a resolution-independent instability. | At t = 0.3 s the degradation is 2.3x (current), 5.4x (volume) and 2.1x (lead) smaller at N = 142 than at N = 100: it shrinks with h. TWO rungs only, no order stated. | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3932-L3947) |
| The explosive phase starts when the trailing edge crosses the box centre (the geometric signature). | The start-position test (droplet at 5 mm instead of 2.5 mm, 20 mm box, N = 100): jump at t = 0.1372 s with the centroid at 12.15 mm; neither pre-registered prediction holds; the slow phase of the two runs is identical to t = 0.12 s (7.4e-4 to 7.7e-4). | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3831-L3863) |
| The fast phase is the outlet. | In the 10 mm box the current rises a factor 27 in 0.005 s at t = 0.06 s with the leading edge 3.5 mm (35 cells) from the outlet; the 20 mm box has no growth there at N = 100, 142 and 200. This part of the claim survives. | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3779-L3782), [L3826-L3829](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3826-L3829) |
| The late instability is a seam defect. | Serial diverges at 0.0775 s, np 4 at 0.0868 s with the fix and 0.0904 s without it: inside the decomposition scatter. | MEASURED, [METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L311-L312) |

## What survives

1. The outlet triggers the fast phase; the density contrast is necessary (the ratio-1 run completes); the phenomenon is decomposition-independent ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3741)). A benchmark that runs the translating droplet longer than 0.05 s needs a longer box ([SL article, d1e3414](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2702-L2712)).
2. The long-box ladder to 0.1 s at N = 100 / 142 / 200: all three rungs complete; shape and centroid converge near order 1.7 to 1.9, the curvature error does not (1.5 to 2.2 %) ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3800-L3814)).
3. The reading "a slow interior growth in time, amplified as the droplet approaches the outlet" fits every run so far, including the start-position test ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3860-L3863)).
4. The rule of [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L528-L544): before any conclusion from t_blow, compute where the interface is and what the nearest boundary is, see [[concepts/error-vector-and-read-out-instants]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: the [heading corrected in place](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3769-L3771), the [two-phenomena paragraph](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3826-L3829), the [geometric-signature test and its result](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3831-L3863), the [40 mm box](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3917-L3947).
- [x] `METHOD.md`: [section 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L311-L312) and the [`CURVATURE_EXTENSION` row](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391) (the gate baseline diverged at t = 0.077 to 0.094 s; `END_TIME 0.05` ends before both onsets).
- [x] The gcls pre-print: [sec:res-translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L789-L823) and the [baseline's limits](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L921-L935); on 2026-09-28 the paragraph moved to the SL article, [sec:translating-late](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2278-L2315) (commit d1e3414).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The third rung of the 40 mm box (N = 200, 160 000 cells, 78 000 steps) belongs on the cluster; no order of the interior growth can be stated from two rungs ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3944-L3947)).
- [ ] A longer translating box for the gates and ladders (a new length token with the current box as its default) is an OPEN author decision ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3782-L3784)).
- [ ] Not tested: an outlet condition other than fixed `p_rgh` with zero-gradient U ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3862-L3863)) and the semi-implicit capillary force ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3742-L3743)), see [[models/semi-implicit-capillary-force]].
- [ ] The mechanism of the interior growth at the water/air density ratio is the open defect of the translating droplet ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3928-L3930)), see [[concepts/density-ratio-amplifier]].

## Related

Hubs: [[hubs/mass-flux]], [[hubs/verification]]. Siblings: [[cases/translating-droplet]], [[concepts/density-ratio-amplifier]], [[concepts/method-gates]], [[concepts/parasitic-current-mechanism]], [[concepts/coupled-face-density-defect]], [[concepts/force-time-centring]], [[concepts/error-vector-and-read-out-instants]], [[models/semi-implicit-capillary-force]], [[retractions/closed-box-translating-droplet]], [[retractions/t-blow-baseline]], [[studies/gcls-pre-print]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Written from STATUS 11.14 and 11.15, the gcls pre-print and the SL article. Corrected 2026-09-27. Entered in [[retraction-log#2026-09]].
