---
title: "The normal-projected SL design note and pre-print"
description: "The map of the normal-projected semi-Lagrangian documents (design note, pre-print of 2026-08-03, twelve-slide deck): the clean trace, the divergent write-back, the falsified corrugation hypothesis, and the closed line."
aliases: []
kind: study
status: settled
part: advection
tags: [study, part/advection]
date: 2026-09-28
code: [docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex, docs/normal-projected-semi-lagrangian/normal-projected-semi-lagrangian.md, docs/normal-projected-semi-lagrangian/stable-foot-point-3d.md, docs/normal-projected-semi-lagrangian/npsl-presentation/normal-projected-semi-lagrangian-level-set.template.html, src/leiaLevelSet/semiLagrangian/normalProjectedScheme.C]
sources: [nPSL article, nPSL design note sections 11-13, RM lines 121-155 and 502-534, PCT lines 78-95 and 167-187, kb-raw B9, kb-raw A3]
---
# The normal-projected SL design note and pre-print

> **Verdict (2026-09-28).** The normal-projected semi-Lagrangian update (nSL) is a closed line, and its document set says so. The trace along the normal-projected velocity is clean (0.017 h in one step) and preserves a free stream to machine precision, but every per-step geometric write-back of a fit-derived offset diverges at an engine-independent rate (1.7x per ten steps), the corrugation hypothesis that motivated the line is falsified (0.209 h against 0.223 h at 16 h of displacement), and on the reversed vortex the two variants have orders -1.24 and -0.17 against 2.9 for the value path ([`sec:results`, line 393](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L393)). The roadmap's instruction stands: do not promote the projected trajectory ([RM lines 513-534](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L513-L534)); see [[concepts/normal-projected-sl]] and [[models/sl-scheme]].

## What it is

**The pre-print.** `docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex`, 505 lines, one author, committed once (cde884d, 2026-08-03). Pre-print PDF: https://leia-openfoam.github.io/leia/preprints/normalProjectedSemiLagrangian.pdf (the site compiles every `docs/*/*-article/*.tex`; the Makefile has no target for this file). Sections: Motivation ([`sec:motivation`, 78](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L78)); Mathematical formulation ([`sec:method`, 99](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L99)) with the tangential redundancy ([`sec:exact`, 102](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L102)), the departure-centred AB2 trace ([`sec:trace`, 132](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L132)), the geometric update and the two offset engines ([`sec:update`, 191](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L191)), the update ladder ([`sec:ladder`, 262](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L262)) and the strain renormalisation ([`sec:strain`, 286](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L286)); Implementation ([`sec:implementation`, 312](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L312)); Measured verification record ([`sec:results`, 393](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L393)).

**The design note.** `normal-projected-semi-lagrangian.md`, 521 lines, sections 1 to 10 the derivation and test sequence, sections 11 to 13 the first measured results (2026-08-03). Its status line still reads "design stage, not implemented, not measured" ([line 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/normal-projected-semi-lagrangian.md#L3)); the sections below it show that the line is stale. The companion `stable-foot-point-3d.md` (185 lines) documents the stabilised foot-point algorithm that engine (a) uses.

**The deck.** https://leia-openfoam.github.io/leia/decks/normal-projected-semi-lagrangian-level-set.html, 12 flat slides (`#/0` to `#/11`): why, the exact redundancy, the trace, the geometric update, the foot-point engine, the update ladder, implementation, the clean trace `#/8`, the falsified corrugation hypothesis `#/9`, the deforming flow `#/10`, status `#/11`.

**The data.** The theme has no data folder; the tables are inline in the article. The study CSVs sit under the semi-Lagrangian theme: `npslConv2Dvortex_errors.csv`, `nslConv2Dvortex_errors.csv`, `nslConv3Ddeformation_errors.csv`, `nslKinematicGate_errors.csv` in `docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/`. Configs `config/nsl*.yaml` and `config/npsl*.yaml`.

**How to build.** No Makefile target; `latexmk -pdf normalProjectedSemiLagrangian.tex` in the article folder; the deck through `make decks`.

## Why it matters

The line tested the hypothesis that the grid-scale corrugation of the translating interface comes from the per-step value resampling. The measurement located the amplifier in the fitted normals instead, which is a statement about every scheme that reads the fit ([[concepts/curvature-corrugation-and-the-fit]]). Two by-products survive: the departure-centred AB2 trace ([[concepts/departure-foot-ab2-centring]]) and the shared quadratic-root offset conversion of the curvature path.

## Where in the code

`src/leiaLevelSet/semiLagrangian/normalProjectedScheme.{H,C}` (`slScheme` type `normalProjected`; dictionary entries `renormalization`, `offsetEngine`, `bandRadii`, `minGradPsi`, `offsetBeta`, `footPoint*`), the reconstruction hooks `fitDerivatives`, `signedOffset`, `footPointDistance` in `uncachedQuadraticWeightedLeastSquaresReconstruction.{H,C}`, the study token `SL_SCHEME`.

## Evidence

All numbers from [`sec:results`, line 393](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L393) unless stated.

| claim | number | where |
|---|---|---|
| The trace is clean (sigma=0 rigid translation, N=64, R=6.4 h) | one-step error 0.017 h in the band, 0.027 h everywhere; over 1600 steps volume error 1.17e-2 against 1.02e-2 for `pointValue` | paragraph "The trace is clean", MEASURED |
| The geometric write-back diverges, engine-independently | band error after 100 steps 1.107 h (quadratic root), 1.11 h (first-order conversion), 1.113 h (stabilised foot point) against 0.012 h for raw transport; 1.7x per 10 steps, independent of dt | paragraph "The geometric write-back diverges", MEASURED; [PCT lines 167-187](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L167-L187) |
| The corrugation hypothesis is falsified | m > 4 corrugation at 16 h displacement: `pointValue` 0.209 h, nSL 0.223 h | paragraph "The corrugation hypothesis is falsified", MEASURED |
| Reversed vortex, T=2, CFL 1/2, shape error at N=256 | `pointValue` 7.85e-7 (order about 2.9); nSL strain 2.53e-2 (-1.24); nSL geometric plus foot point 2.95e-2 (-0.17) | paragraph "Deforming flow", MEASURED; [PCT lines 78-95](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L78-L95) |
| The strain factor applied to the whole domain | band gradient errors 1e6 to 1e34 across the ladder; the reason for the band gate | same paragraph, MEASURED |
| Coupled translating droplet, N=32, `normalProjection` trajectory | latest runaway of the roadmap matrix, 0.04365 s, but physically invalid by t=0.03 (max U 0.248 m/s against U0=0.05) | [RM lines 121-155](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L121-L155), MEASURED |
| Fixed-Courant refinement, zero force | `normalProjection` volume orders 1.45 and 0.81, zero-set orders 1.17 and 1.11 against 2.47, 2.55 and 2.55, 1.89 for the physical velocity | [RM line 302](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L302), MEASURED |

## Why it failed, or why we think so

Writing any fit-derived offset back into psi feeds neighbour noise into the next fit at order one per step; the offset carries no factor of dt, so the engine's accuracy does not change the loop gain ([`sec:results`](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L393)). The corrugation re-enters through the fitted normals: a grid-scale wiggle perturbs the gradient, the normal, the trace direction and the written increment with the same gain per displacement as the resampling it replaced. The strain mode, the only renormalisation with a factor dt in its gain, is anti-convergent on a deforming flow.

## Decisions

Not promoted ([RM "P1", lines 513-534](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L513-L534)); `pointValue` stays the production scheme ([[models/sl-scheme]]). The article's proposed next step, the strain factor on top of the unchanged value transport, has no measurement in the sources of this note.

## Retracted or superseded inside it

- The design note's status line (design stage) is stale; sections 11 to 13 of the same note and the article report the measurements.
- The study `npslConv2Dvortex` (np 4) and the `nslConv*` studies are on the gradU contamination list; the scheme reads `grad(U)` at `normalProjectedScheme.C:161` ([post-mortem section 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L102-L146); [[retractions/gradu-coupled-patch-contamination]]). The four-orders gap to the baseline is far above any plausible seam effect, but the fitted orders are provisional.
- The brief `improvement-drift-gate.md` that `IMPROVEMENTS.md` links for this theme does not exist in the repository ([IMPROVEMENTS.md line 40](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/IMPROVEMENTS.md#L31)).

## What it does not cover

3D runs, polyhedral meshes, and any coupled run with the nSL scheme beyond the roadmap's N=32 matrix.

## Related

[[hubs/advection]], [[hubs/method-lines]]. [[models/sl-scheme]], [[models/sl-reconstruction]], [[concepts/normal-projected-sl]], [[concepts/departure-foot-ab2-centring]], [[concepts/curvature-corrugation-and-the-fit]], [[concepts/trace-velocity-projected-flux]]. Siblings: [[studies/sl-quadratic-pre-print]], [[studies/sl-linear-pre-print]], [[studies/poly3d-roadmap]].

## Log

### 2026-09-28
Created from the article, the design notes and the deck at 8867581.
