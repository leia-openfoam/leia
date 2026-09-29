---
title: "The geometrically redistanced level set pre-print (draft)"
description: "The map of the GRL draft (484 lines, static gates written, advected sections TODO), its two decks and its data; the static results that stand, and why the line is closed."
aliases: []
kind: study
status: open
part: advection
tags: [study, part/advection]
date: 2026-09-28
code: [docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex, docs/geometrically-redistanced-levelset/grl-level-set-presentation/geometrically-redistanced-level-set.template.html, docs/geometrically-redistanced-levelset/grl-level-set-presentation/geometrically-redistanced-level-set-negative-results.template.html, src/leiaLevelSet/redistancer, applications/solvers/leiaRedistancedLevelSetFoam]
sources: [GRL article, GRL decks, MC article sec:frozen, PCS dead end 7, IMPROVEMENTS.md, kb-raw B6]
---
# The geometrically redistanced level set pre-print (draft)

> **Verdict (2026-09-28).** The GRL document is a draft: the method and the two static gates are written, the advected sections and the conclusions are marked TODO ([title, line 40](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L40); last commit d2b50f3, 2026-08-06). The static results stand: `planeFootWave` is second order in the band on a circle and machine-exact on a plane, `anchoredEikonal` is first order in the band, and PDE reinitialisation raises the band error at every resolution ([`sec:static`, line 355](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L355); [`sec:idempotency`, line 385](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L385)). The line is closed for transport: the frozen-band variant, built to have zero interface displacement, still inflates the volume error 0.017 to 1.77 at N=256 on the T=8 vortex ([MC `sec:frozen`, line 304](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)), the curvature plan lists plane-based band rewrites as dead end 7 ([PCS line 363](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L363)), and the improvement briefs skip GRL on purpose ([IMPROVEMENTS.md line 45](https://github.com/leia-openfoam/leia/blob/8867581/docs/IMPROVEMENTS.md#L34-L36)). See [[concepts/redistancing-geometric-grl]] and [[models/redistancer]].

## What it is

**The article.** `docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex`, 484 lines, `elsarticle`, marked "Draft skeleton". Pre-print PDF: https://leia-openfoam.github.io/leia/preprints/geometricallyRedistancedLevelSet.pdf. Sections: Introduction, TODO ([`sec:intro`, 75](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L75)); Method: transport ([86](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L86)), cell-local least-squares planes ([`sec:llsplanes`, 97](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L97)), anchors ([`sec:anchors`, 135](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L135)), the donor-plane wave ([`sec:wave`, 222](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L222)), algorithm ([`sec:algorithm`, 310](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L310)), trigger ([336](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L336)); Software, TODO ([343](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L343)); Verification: the static gate ([`sec:static`, 355](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L355)), one-step displacement ([`sec:idempotency`, 385](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L385)), reversed vortex, 3D and trigger ablation, all TODO ([433-447](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L433-L447)); Negative results ([`sec:negative`, 448](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L448)); Limitations and Conclusions, TODO ([471-478](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L471-L478)).

**The decks.** https://leia-openfoam.github.io/leia/decks/geometrically-redistanced-level-set.html (23 sections in 8 groups: `#/2` level set and indicator; `#/3` planeFootWave, 6 vertical slides; `#/4` the solver; `#/5` verification, 3 slides including "What the advected studies measured (so far)"; `#/6` conclusions, skeleton; `#/7` software design) and https://leia-openfoam.github.io/leia/decks/geometrically-redistanced-level-set-negative-results.html (29 sections in 8 groups: `#/2` the gate and the governing lesson; `#/3` foot-cloud scalloping; `#/4` the Eikonal fill; `#/5` PDE reinitialisation, two failures; `#/6` measurement pitfalls; `#/7` scoreboard and standing rules).

**The data.** `data/tables/` (10 entries: `bulkVortexGRL_errors.csv`, `redistanceCircle2D_errors.csv`, `redistanceStatic2D_errors.csv`, `vortexBoundsGRL_errors.csv`, `vortexThresholdGRL_errors.csv`, `vortexTriggerGRL_errors.csv`, `vortexUnclippedGRL_errors.csv`, `grl_convergence_orders.tex`, `redistanceStatic2D_orders.tex`, `static_redistance_orders.tex`) and `data/figures/` (6 PNG files). No `data/archive/` folder.

**How to build.** `make article-grl`; `make studies-grl` (the `GRL_STUDIES` configs) ([Makefile lines 138-139 and 267-269](https://github.com/leia-openfoam/leia/blob/8867581/Makefile#L138-L139)); `make decks`.

## Why it matters

The draft records the one geometric redistancing construction that passes a static gate, and the reason no redistancing can serve transport here: a redistancer rebuilds the distance to whatever zero set the field carries, including advection artefacts. That lesson closed the line and made the source-term route the only psi-maintenance candidate ([[hubs/gradient-control]]).

## Where in the code

`src/leiaLevelSet/redistancer/` (`noRedistancing`, `PDE`, `anchoredEikonal`, `planeFootWave`; trigger `interval`, `gradPsiThreshold`, `signedDistanceBounds`), `applications/solvers/leiaRedistancedLevelSetFoam`, `applications/test/leiaTestRedistance`, `cases/2DredistanceCircle`, `cases/2DredistanceStatic`.

## Evidence

The static-gate tables report the band maximum norm; the record has no L2 for those gates, so the numbers below are maximum-norm values and are quoted as such.

| claim | number | where |
|---|---|---|
| Static gate, one event on the tanh circle, band error after the event at h=1/256 | `planeFootWave` 4.48e-5 (orders 1.97, 2.00, 2.54 between rungs); `anchoredEikonal` 3.53e-4 (1.98, 2.06, 2.35); PDE 4.07e-2 against 2.03e-2 before the event | [`sec:static`, 355](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L355) and [`static_redistance_orders.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/static_redistance_orders.tex), MEASURED (maximum norm) |
| One step on an exact plane | band error 1.61e-13 at h=1/256, volume change 2.36e-16 | [`redistanceStatic2D_orders.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/redistanceStatic2D_orders.tex), MEASURED |
| One step on an exact circle, spurious volume change | 1.11e-5 at h=1/32 to 1.17e-7 at h=1/256 (table); the prose quotes 1.05e-4 to 1.72e-6 and 1.5 % to 0.02 % of the phase volume | [`sec:idempotency`, 385](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L385) and the same table, MEASURED; the prose and the table disagree |
| Shape error equals volume error (one-signed displacement) | equal at three of the four rungs in the table (1.20e-7 against 1.17e-7 at the finest) | same table, MEASURED |
| PDE reinitialisation on the exact circle | band error 1.97e-1 at h=1/32 | same table, MEASURED (maximum norm) |
| Advected 2D vortex `bulkVortexGRL`, shape error, h=1/32 to 1/128 | `noRedistancing` 6.97e-3 to 2.20e-3; `planeFootWave` 3.59e-2 to 5.09e-2 (grows); `anchoredEikonal` about 5.0e-2 flat; PDE 4.97e-2 to 4.43e-3 | [`grl_convergence_orders.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/grl_convergence_orders.tex), MEASURED; the article's advected sections are TODO |
| The trigger degenerates to every-step firing | the ablation is identical to 3 digits; thresholds 0.01 and 0.05 change nothing | deck https://leia-openfoam.github.io/leia/decks/geometrically-redistanced-level-set.html#/5/4, MEASURED |
| Frozen-band variant under Eulerian transport, T=8 | volume error 0.017 to 1.77 at N=256, 0.238 to 2.53 at N=128; the indicator inflates to 2.8x the droplet | [MC `sec:frozen`, 304](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304), MEASURED |

## Why it failed, or why we think so

Measured, in the negative-results deck and `sec:negative` ([448](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L448)): the foot-point-cloud fill scallops the gradient (replaced by the donor-plane evaluation); the anchored-Eikonal fill is first order in the band; a fixed pseudo-time step of the PDE fill blows up (7.5e6 on the gate) and the CFL-safe step still injures the field (band mean 0.537 to 1.28). Under transport the fill rebuilds distance to spurious bulk zero crossings of the advected field, so even a zero-displacement event injures the run ([MC `sec:frozen`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)).

## Decisions

The line is closed: [[models/redistancer]], [PCS dead end 7](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L363), [IMPROVEMENTS.md](https://github.com/leia-openfoam/leia/blob/8867581/docs/IMPROVEMENTS.md#L34-L36). The standing rules of the negative deck (`#/7/3`): an event must not injure; trigger policy is downstream of fill quality; assert each model's contract region; mesh-relative parameters only; distances come from continuous geometry.

## What it does not cover

3D, polyhedral and perturbed meshes, two-phase coupling, the trigger ablation, and every advected study named in the deck (the article sections are TODO).

## Related

[[hubs/advection]], [[hubs/gradient-control]], [[hubs/method-lines]]. [[models/redistancer]], [[models/phase-indicator]], [[concepts/redistancing-geometric-grl]], [[concepts/eulerian-fv-transport]]. Siblings: [[studies/method-comparison]], [[studies/sdpls-pre-print]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from the draft at 8867581, the two deck templates and the tables.
