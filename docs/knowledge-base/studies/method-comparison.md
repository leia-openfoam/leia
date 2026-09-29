---
title: "The method comparison article and deck"
description: "The map of the method-comparison theme: the 354-line decision article on the reversed vortex, the 79-section deck with its capillary-balance and time-centring tracks, and the 134-table data folder that also hosts the curvature campaign and the VOID closed-box tables."
aliases: []
kind: study
status: settled
part: advection
tags: [study, part/advection]
date: 2026-09-28
code: [docs/method-comparison/method-comparison-article/methodComparison.tex, docs/method-comparison/method-comparison-presentation/level-set-method-comparison.template.html, docs/method-comparison/method-comparison-article/REPRODUCE.md, workflow/Snakefile.comparison]
sources: [MC article, MC deck, MC REPRODUCE.md, VOID_closedBox_20260902 README, kb-raw B3]
---
# The method comparison article and deck

> **Verdict (2026-09-28).** The article answers one question: which level-set advection line to choose at a given resolution. Given sufficient resolution the quadratic semi-Lagrangian transport is the most accurate and the cheapest line at every measured resolution and horizon: at 512^2 and T=8 it is 20x more accurate than the best Eulerian variant at half the wall clock and 27x cheaper than velocity extension ([Verdict, line 194](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194); [decision table](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex)). That decision stands ([[decisions/sl-reconstruction-uncached-qwls]], [[concepts/eulerian-fv-transport]]). The deck goes further than the article: its capillary-balance track (`#/8`, 34 slides) and its time-centring track (`#/9`, 10 slides) record the July and August surface-tension gates, and one of its own slides retracts the explicit-force reading ([[retractions/force-at-n-not-n-plus-1]]). The article's limitations carry two SDPLS statements that the SDPLS pre-print later retracted, and the data folder holds the `VOID_closedBox_20260902` tables that must not be cited ([[retractions/closed-box-translating-droplet]]).

## What it is

**The article.** `docs/method-comparison/method-comparison-article/methodComparison.tex`, 354 lines, last commit d2b50f3 (2026-08-06). Pre-print PDF: https://leia-openfoam.github.io/leia/preprints/methodComparison.pdf. Sections: The decision question ([53](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L53)); Setup ([62](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L62)): the reversed single vortex, N=32 to 512, T in {2, 8}, CFL 0.5, np 4 everywhere, one solver `leiaLevelSetFoam`; Results ([92](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L92)) with convergence, cost against accuracy, the decision table, the detailed table, the flux-form volume loss ([`sec:fluxloss`, 147](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L147)) and the field atlas; Verdict ([194](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194)); Limitations and outlook per method ([208](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L208)) with the measured SL improvements ([`sec:slimp`, 258](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L258)) and the frozen-band redistancing ([`sec:frozen`, 304](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)); Reproducibility ([344](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L344)).

**The deck.** https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html, 79 sections in 11 groups: `#/2` the contenders; `#/3` convergence at T=2 and T=8 (7 vertical slides); `#/4` cost against accuracy; `#/5` the decision table and the detailed comparison; `#/6` the field atlas; `#/7` faults and improvements per method, the SL improvements, the flux-form loss and the frozen-band result (10); `#/8` "Surface tension now has one enforceable face-flux contract" and the capillary-balance track: the exact-curvature control, five models failing before 4 ms, the replay gate, the face-centred curvature gate, the connected interface, the shared face-curvature service, the manufactured mode gate, the oscillating benchmark, the oracles, the pressure gates (34); `#/9` where the capillary force sits in time, the amplification matrix, the step limit, what is retracted (10); `#/10` the final verdict.

**The data.** `data/tables/` (134 entries) holds three things: the benchmark CSVs and tables (`benchVortex*_errors.csv`, `benchVortex_decision.tex`, `benchVortex_detailed.tex`, `advConv2D*_convergence.csv`, `flux_volloss_timeseries.csv`); the curated tables of the curvature and viscosity campaigns (`face_curvature_orders*.{csv,tex}`, `curvature_gain*.{csv,tex}`, `curvatureModeTransferGate*.csv`, `mode_rate_*.csv`, `mufGrid2D_*_errors.csv`, `mufDecide2D_*`, `kickOriginGate2D_errors.csv`, `domain_size_control.csv`, `driverSplitDeliveryProbe_errors.csv`, `foot_evaluated_face_coupled.csv`, `cell_centre_inverse_coupled.csv`, `oscIST*`, `drop3dTraceGate_*`); and the folder `VOID_closedBox_20260902/` with a README that forbids citing anything in it. `data/figures/` (84 entries: the `benchVortex_*` and `bench_*` atlases, `face_curvature_convergence*.png`, `droplet_filtered_evolution.png`, `bestConfigTranslating.gif`) and `data/animations/` (one mp4). `REPRODUCE.md` lists every study behind the curvature tables with its cost. No `data/archive/` folder.

**How to build.** `make comparison` (`workflow/Snakefile.comparison` runs the six `config/benchVortex*.yaml` studies, harvests the figures and the decision table, rebuilds the deck and compiles the article, [Makefile lines 208-212](https://github.com/leia-openfoam/leia/blob/8867581/Makefile#L208-L212)); `methodComparison.tex` is not in `make articles`.

## Why it matters

It is the one place where all four advection lines run through one solver at one decomposition, so the wall clocks compare. Its verdict fixed the production transport, and its data folder became the curated home of the surface-tension campaign of 2026-07 and 2026-08 ([[studies/poly3d-roadmap]], [[studies/curvature-stabilization-campaign]]).

## Where in the code

`applications/solvers/leiaLevelSetFoam` (one dictionary selects the method), `workflow/Snakefile.comparison`, `config/benchVortex*.yaml`, `workflow/scripts/make_decision_table.py` and the figure scripts named in `REPRODUCE.md`.

## Evidence

| claim | number | where |
|---|---|---|
| Decision table, T=8, N=512, shape error and wall clock | SL quadratic 1.49e-5 at 595 s; Eulerian 3.04e-4 at 1170 s; Eulerian plus SDPLS R 5.83e-4; VE `closestPoint` 1.44e-3 at 1.61e4 s | [`benchVortex_decision.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex), MEASURED |
| T=2, N=512 | SL 1.72e-5 at 160 s (cached) and 190 s (uncached); Eulerian 1.68e-5 at 323 s | same table, MEASURED |
| Flux-form volume loss at N=128, T=8 | -17 % by t=6.5, monotone; band gradient 0.72 against 0.96 for the point-value scheme; the integral of psi conserved to machine precision | [`sec:fluxloss`, 147](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L147), MEASURED; [[models/sl-scheme]] |
| Limiters | Barth-Jespersen collapses the order 3.0 to 0.1; Venkatakrishnan to 0.9 | [`sec:slimp`, 258](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L258), MEASURED; [[models/sl-value-bound]] |
| Householder QR against Cholesky | bit parity on hexahedra and on a 10 %-perturbed mesh | same section, MEASURED; [[decisions/sl-fit-normal-equations]] |
| Flux-form SL | conserves the integral of psi to 2e-14; volume error 3 to 13x the point-value scheme at T=8 (0.82 against 0.062 at N=64) | same section, MEASURED |
| Frozen-band redistancing | volume error 0.017 to 1.77 at N=256 and 0.238 to 2.53 at N=128 | [`sec:frozen`, 304](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304), MEASURED; [[models/redistancer]] |
| The exact-curvature control is the only balanced translating run | 2.17e-7 m/s against 9.06e-2 for the computed curvatures at N=64 | deck https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/8/2, MEASURED; [[concepts/well-balanced-exact-curvature-gate]] |
| Spatial curvature variation is the active defect (one-step replay) | interface-mean curvature 3.84e-9 against quadratic cell centre 2.03 at the t=0.05 snapshot | deck `#/8/4`, MEASURED; [RM line 724](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L724) |
| Face-centred curvature gate | stabilised foot point h^2.04, L2 0.10 against 11.4 at N=512 | deck `#/8/5`, MEASURED; [[concepts/face-curvature-deliveries]] |
| Pressure gates on perturbed meshes | PCG removes the GAMG artefact 18.2x at N=64 and 28.8x at N=128 and leaves a plateau of about 4.8e-5 m/s | deck `#/8/33`, MEASURED; [[concepts/pressure-projection-and-linear-solvers]] |
| Time centring | the force is built from psi at n+1 (symplectic Euler, det M = 1); the n+1/2 centring is spectrally identical; the practical step limit is omega_grid dt about 1.0 to 1.3 | deck `#/9/5`, `#/9/6`, `#/9/9`, DERIVED and MEASURED; [[concepts/force-time-centring]], [[concepts/capillary-time-step]] |

## Retracted or superseded inside it

1. The limitations state that the SDPLS source `R` over-flattens near the resolution limit and that `beta` with an explicit discretisation diverges at every CFL ([line 208](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L208)). Both were measured with the reversed-sign source branch and are retracted in the SDPLS pre-print ([[studies/sdpls-pre-print]], [[models/sdpls-source]]).
2. The deck's time-centring track retracts its own first reading ("our force is explicit at t^n, anti-damping +22 to 88 1/s", slide `#/9/10`) ([[retractions/force-at-n-not-n-plus-1]]). The midpoint centring was later measured 32 % worse on the translating droplet (2026-09-27, [[concepts/force-time-centring]]).
3. Every table in `data/tables/VOID_closedBox_20260902/` came from the closed-box translating case; the README there says "Do not curate, cite or re-use anything here" ([[retractions/closed-box-translating-droplet]]).
4. The candidates of the capillary-balance track (`connectedInterface`, `helmholtzPreserveModes`, the integral models) are all failed or retired ([[models/curvature-extension]], [[models/surface-tension-force]], [[concepts/integral-surface-tension-cst]]); the production delivery is `cellCentreInverse` ([[decisions/curvature-extension-cell-centre-inverse]]).
5. The benchmark studies ran at np 4 with the SDPLS and VE arms consuming `grad(U)`, so those arms carry the coupled-patch defect ([[retractions/gradu-coupled-patch-contamination]]); the SL and plain Eulerian arms use the flux and are unaffected by that defect.

## What it does not cover

Two-phase coupling in the article (the deck's capillary track does), 3D, polyhedral meshes, the gradient-control laws and the halo-limited extension ([[studies/gcls-pre-print]]), and any measurement after 2026-08-06 in the article.

## Related

[[hubs/advection]], [[hubs/surface-tension]], [[hubs/verification]], [[hubs/method-lines]]. [[models/level-set-advection]], [[models/sl-reconstruction]], [[models/sdpls-source]], [[models/velocity-extension]], [[concepts/eulerian-fv-transport]], [[concepts/richardson-ladders-and-orders]], [[cases/kinematic-advection-cases]], [[cases/curvature-static-gates]]. Siblings: [[studies/sl-quadratic-pre-print]], [[studies/sdpls-pre-print]], [[studies/velocity-extension-pre-print]], [[studies/grl-pre-print]], [[studies/poly3d-roadmap]], [[studies/shannon-parasitic-currents-campaign]].

## Log

### 2026-09-28
Created from the article, the deck template, the data folder and REPRODUCE.md at 8867581.
