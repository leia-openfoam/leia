---
title: "The velocity extension deck and the stub article"
description: "The map of the velocity-extension theme: a 162-section deck that is the document of record, an 88-line article stub, 60 figures and no tables; the per-model verdicts, and the later measurements that dispute or supersede them."
aliases: []
kind: study
status: open
part: gradient-control
tags: [study, part/gradient-control]
date: 2026-09-28
code: [docs/velocity-extension/velocity-extension-presentation/velocity-extension.template.html, docs/velocity-extension/velocity-extension-article/velocityExtension.tex, src/leiaLevelSet/velocityExtension, applications/test/leiaTestVelocityExtension]
sources: [VE deck, VE article stub, MC decision table, SDPLS article sec:coupledpatch-consequences, RM lines 121-134, gcls article sec:res-hl0, PCS section 15.3, kb-raw B7]
---
# The velocity extension deck and the stub article

> **Verdict (2026-09-28).** The document of record for the velocity-extension family is the deck, not the article: the article is an 88-line skeleton with a placeholder abstract and TODO sections ([title, line 23](https://github.com/leia-openfoam/leia/blob/8867581/docs/velocity-extension/velocity-extension-article/velocityExtension.tex#L23); committed once, 2026-07-11). The deck ranks the models on the reversed vortex at N=256 (`meshWave` the robust default, `closestPoint` the accurate one where psi stays near a signed distance, `steadyUpwindLinear` divergent; slide `#/15/2`). Three later records stand against the family: for pure advection it is dominated at 12 to 27x the cost ([[studies/method-comparison]]); in the 2026-09 method gate the `closestPoint` reference candidate FP0 fails, with a 45 % decomposition dependence of its kinematic serial check ([[studies/gcls-pre-print]]); and the "2185x worse" claim of the record is disputed after the coupled-patch fixes ([[studies/sdpls-pre-print]]). The status is open until the article is written and the family is re-measured on the fixed binaries.

## What it is

**The deck.** https://leia-openfoam.github.io/leia/decks/velocity-extension.html, 162 sections in 19 horizontal groups: `#/2` the level set and the phase indicator (9 vertical slides); `#/3` the models, why, the interface velocity, the non-invasive coupling, the map of the seven stacks (4); `#/4`, `#/5` the framework and its algorithm; `#/6` to `#/12` one group per model with motivation, model, discretisation, formulation, code and measured verdict: `none` (`#/6`), `anisotropicDiffusion` (`#/7`), `pseudoTime` (`#/8`), `steadyUpwind` (`#/9`), `steadyUpwindLinear` (`#/10`), `closestPoint` (`#/11`), `meshWave` (`#/12`); `#/13` static verification (7); `#/14` advected verification, reversed vortex, reversibility bias, the non-reversing steady vortex (17); `#/15` conclusions, ranked by error, which extension by problem (2); `#/16` software design (6); `#/17`, `#/18` the field atlases of alpha and of the gradient defect, six models times four horizons (24 each).

**The article.** `docs/velocity-extension/velocity-extension-article/velocityExtension.tex`, 88 lines: sections Introduction, Velocity-extension models, Verification and Conclusions are TODO comments; one placeholder figure `convergence_none.png`. The site compiles it as https://leia-openfoam.github.io/leia/preprints/velocityExtension.pdf. The improvement brief `velocity-extension/improvement-metric-footpoint.md` that `IMPROVEMENTS.md` lists ([line 37](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/IMPROVEMENTS.md#L29)) does not exist in the repository.

**The data.** `velocity-extension-article/data/figures/` (60 PNG files: `alpha_evo_<model>_T*.png` and `graderr_<model>_T*.png` for six models times four horizons, `convergence_{none,anisotropicDiffusion,pseudoTime}.png`, `static_convergence.png`, `static_e_anchor0.png`, `maxdef_convergence.png`, `volume_history_T8.png`, `steady_*.png`, `interface_grid.png`, `static_indicator_volume_convergence.png`). No `data/tables/` and no `data/archive/` folder. The deck regenerates from `config/staticExtension.yaml`, `config/bulkVortexHighRes.yaml` and `config/steadyVortex2D.yaml` (slide `#/15/3`).

**How to build.** `make studies-ve` (the `VE_STUDIES` configs, [Makefile line 136](https://github.com/leia-openfoam/leia/blob/8867581/Makefile#L136-L137)); `make decks`; the article has no Makefile target.

## Why it matters

A velocity extension is the second route to gradient control: transport psi with a velocity whose normal strain vanishes on the interface ([[hubs/gradient-control]]). The family of 2026-07 measured that route for the Eulerian solver; the halo-limited extension of 2026-09 ([[concepts/halo-limited-extension]]) and the `closestPoint` reference were then gated together in the semi-Lagrangian solver, and both failed ([[studies/gcls-pre-print]]).

## Where in the code

`src/leiaLevelSet/velocityExtension/` (base `velocityExtension`, intermediate `interfaceExtension` with `nLayers`, `nAnchorLayers`, `projectFlux`, `fadeMode`; members `none`, `anisotropicDiffusion`, `pseudoTime`, `steadyUpwind`, `steadyUpwindLinear`, `closestPoint`, `meshWave`, `haloLimited`), `applications/test/leiaTestVelocityExtension`, the kinematic solver's `traceFlux extension` switch in `leiaSemiLagrangeLevelSetFoam`.

## Evidence

The deck quotes the reversed-vortex numbers as "band gradient error / volume error" at N=256.

| claim | number | where |
|---|---|---|
| Static extension defect ratio to the raw velocity, one `correct()` at t=0 | `closestPoint` 0.27 to 0.0045 over h=1/32 to 1/256 (about O(h^2)); `meshWave` 0.16 to 0.19 flat; `steadyUpwind` 0.11 to 0.34 and degrading; `anisotropicDiffusion` 0.55 to 0.72; `pseudoTime` 0.47 flat; `steadyUpwindLinear` diverges to 1e1 to 1e2 | https://leia-openfoam.github.io/leia/decks/velocity-extension.html#/13/4, MEASURED |
| Reversed vortex, T=2, N=256 | `meshWave` 0.082 / 0.3 %; `closestPoint` 0.104 / 0.015 %; `pseudoTime` 0.085 / 0.2 %; `steadyUpwind` 0.105 / 0.4 %; `anisotropicDiffusion` 0.125 / 2.6 %; `none` 0.070 / 0.02 % | https://leia-openfoam.github.io/leia/decks/velocity-extension.html#/15/2, MEASURED |
| Reversed vortex, T=8, N=256 | `meshWave` 0.21 / 0.5 %; `pseudoTime` 0.36 / 0.5 %; `closestPoint` 0.61 / 55 % (the Newton walks fail when psi is far from a signed distance) | same slide, MEASURED |
| The reversibility bias | `none` band L2 1.67 at t=T/2 against 0.0013 at t=T: about 1300x cancellation | https://leia-openfoam.github.io/leia/decks/velocity-extension.html#/14/13, MEASURED; [[concepts/error-vector-and-read-out-instants]] |
| Steady (non-reversing) vortex, h=1/64 | `none` annihilates the interface at t about 2.65, `meshWave` at 2.85, the other four reach T=3 | https://leia-openfoam.github.io/leia/decks/velocity-extension.html#/14/16, MEASURED |
| `closestPoint` in parallel (steady vortex, np 4) | fallback share 5.32 % to 1.28 % at N=64 and 0.26 % at N=128 with the halo; parallel-serial discrepancy 1.19e-3 to 3.37e-4 | https://leia-openfoam.github.io/leia/decks/velocity-extension.html#/11/7, MEASURED; predates the `updateFlux` fix |
| Method comparison, T=8, N=512 | VE `closestPoint` shape error 1.44e-3 at 1.61e4 s against SL 1.49e-5 at 595 s; "dominated at every resolution" | [MC decision table](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex), MEASURED |
| Coupled translating droplet, N=32, roadmap matrix, time of runaway | `none` 0.0103 s, `meshWave` 0.0131 s, `closestPoint` 0.0258 s, `steadyUpwind` 0.0388 s | [RM lines 121-134](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L121-L134), MEASURED |
| Method gate 2026-09, FP0 (`closestPoint` in the SL solver) | shear shape error 34x the baseline, order 1.4; kinematic serial check 45 % between one and four ranks | [gcls `sec:res-hl0`, line 720](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L720), MEASURED |
| The 2185x claim, re-measured serially after both coupled-patch fixes | band gradient ratio 0.172 (N=128) and 0.267 (N=256), better; shape 9.1x and 5.7x worse, volume 12.6x and 7.6x worse | [SDPLS `sec:coupledpatch-consequences`, 2805](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2805), MEASURED; disputed, not replaced |

## Why it failed, or why we think so

Per model, from the verdict slides: `anisotropicDiffusion` has an intrinsic regularisation floor (the cross-diffusion of a near-rank-one tensor Laplacian is itself O(h)); `pseudoTime` covers two cells of a three-cell band, a reach failure; `steadyUpwind` transports the O(1) seed-staircase noise unchanged while the upwind diffusion that smeared it shrinks; `steadyUpwindLinear` differentiates noisy seeds into an undamped steady fixed point; `closestPoint` needs psi near a signed distance and breaks at T=8 without redistancing; `meshWave` quantises the closest-point map to seeds, O(h), never converges and never fails (`#/7/6`, `#/8/6`, `#/9/7`, `#/10/5`, `#/11/7`, `#/12/6`). For pure advection the whole family costs 12 to 27x for the worst accuracy ([MC line 208](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L208)).

## Retracted or superseded inside it

1. "`closestPoint` is 2185x worse than plain advection on the band gradient" is marked disputed pending an N=512 re-measurement ([SDPLS `sec:coupledpatch-consequences`](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2805)).
2. `interfaceExtension::updateFlux` wrote the raw flux over processor faces until the fix of 2026-08; every parallel extension result in the deck predates it ([SDPLS `sec:extfluxdefect`, 2776](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2776)). `closestPoint.C:432` and `steadyUpwindLinear.C:75` also consume the biased `grad(U)` of the `setVelocity` defect ([[retractions/gradu-coupled-patch-contamination]]).
3. "Doing nothing keeps psi more regular" holds only in the reversed reading at t=T; the steady vortex and the T/2 reading invert it (`#/14/14`, `#/14/16`).
4. The extension models aborted in the two-phase SL solver on 2026-08-12 with "failed lookup of alpha" (the phase field is `alpha.water` there) ([PCS section 15.3, line 1260](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1260-L1298)); the 2026-09 gate ran `closestPoint` in that solver, so the lookup was repaired in between; the fixing commit is not recorded in the sources of this note.
5. The roadmap's July rule "`steadyUpwind` remains the materialised capillary-case default" ([RM line 525](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L525-L527)) is superseded: the production trace velocity is `projectedFlux` with `traceFlux physical` and no extension ([[decisions/sl-trace-velocity-projected-flux]], [[concepts/trace-velocity-projected-flux]]).

## What it does not cover

The halo-limited extension (2026-09, [[concepts/halo-limited-extension]]), 3D, polyhedral meshes, the coupled solvers (only the roadmap and the gate measure those), and any re-run on the fixed binaries.

## Open questions

- The article: everything but the scaffold.
- The N=512 re-measurement that settles the disputed claim.
- Whether `meshWave` or `closestPoint` is worth a coupled gate arm after HL0 and FP0 failed ([[concepts/gradient-control-next-experiments]]).

## Related

[[hubs/gradient-control]], [[hubs/advection]], [[hubs/method-lines]]. [[models/velocity-extension]], [[concepts/halo-limited-extension]], [[concepts/closest-point-extension]], [[concepts/eulerian-fv-transport]], [[concepts/seam-checks-and-decomposition-invariance]]. Siblings: [[studies/method-comparison]], [[studies/sdpls-pre-print]], [[studies/gcls-pre-print]], [[studies/poly3d-roadmap]].

## Log

### 2026-09-28
Created from the deck template, the article stub and the figures folder at 8867581.
