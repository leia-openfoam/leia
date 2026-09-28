---
title: "levelSetAdvection: eulerian vs semiLagrangian"
description: "The two transport lines of the kinematic solver leiaLevelSetFoam; semiLagrangian is 20x more accurate at half the wall clock at 512^2, T = 8, and is the production transport of the SL two-phase solver."
aliases: [levelSetAdvection, ADVECTION]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-08-27
decided_by: [config/benchVortex.yaml, config/uncachedConv2Dvortex.yaml]
code: [src/leiaLevelSet/advection/levelSetAdvection.H, src/leiaLevelSet/advection/eulerianAdvection.H, src/leiaLevelSet/advection/semiLagrangianAdvection.H]
sources: ["MC article sec:decision (L124-L132) and sec:verdict (L194-L205)", "SL article sec:verification (L1200-L1271)", "METHOD 8 (L343-L347)", "DP L250-L253"]
---
# levelSetAdvection: eulerian vs semiLagrangian

> Verdict (2026-09-28). `levelSetAdvection` selects how `psi^n` becomes `psi^{n+1}` in the unified kinematic solver `leiaLevelSetFoam` ([levelSetAdvection.H L30-L55](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/levelSetAdvection.H#L30-L55)). The code default is `eulerian` ([levelSetAdvection.C L69](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/levelSetAdvection.C#L69)), and so is the token `ADVECTION` ([DP L250](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L250)). The measured choice is `semiLagrangian`: at 512^2 and T = 8 it is 20x more accurate than the best Eulerian variant at half the wall clock ([MC article L195-L199](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194)), with shape orders 2.843 (2D vortex) and 2.955 / 3.285 (3D shear, hex / poly) ([sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv)). The two-phase solvers do not use this family: the SL solver builds `slAdvection` directly, the Eulerian solver has its own psi equation.

## What it is

One advection model composes with the independent redistancer, phase indicator and narrow band models; the whole level-set method is configured from the `fvSolution` `levelSet` dictionary ([levelSetAdvection.H L30-L45](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/levelSetAdvection.H#L30-L45)). The key is `levelSet.advection.type`.

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| Eulerian FV transport | `eulerian` | settled, robust second | Implicit transport in the advective (conservative-minus-compression) form with a deferred-correction loop (`nDefCorr`, default 3); composes `velocityExtension` and `sdplsSource`; 20x behind SL at 512^2, T = 8. | [eulerianAdvection.H L30-L43](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/eulerianAdvection.H#L30-L43), [MC article L220-L224](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L220-L224) |
| semi-Lagrangian transport | `semiLagrangian` | settled, the choice | Characteristic update with the second-order Taylor foot and the runtime-selected reconstruction; owns the previous-step velocity `u^n` from the prescribed velocity model. | [semiLagrangianAdvection.H L30-L36](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/semiLagrangianAdvection.H#L30-L36), [MC article L195-L199](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194) |

The dictionary words are the `TypeName` strings ([eulerian L75](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/eulerianAdvection.H#L75), [semiLagrangian L74](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/advection/semiLagrangianAdvection.H#L74)).

## Where each line is production

| solver | transport | composition root | evidence |
|---|---|---|---|
| `leiaLevelSetFoam` (kinematic, unified) | `levelSetAdvection::New` (`eulerian` or `semiLagrangian`) | [createFields.H L90](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetFoam/createFields.H#L90) | the only solver that builds this family |
| `leiaSemiLagrangeLevelSetFoam` (kinematic SL) | `slAdvection::New` | [createFields.H L85](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangeLevelSetFoam/createFields.H#L85) | the advection regression set runs here |
| `leiaSemiLagrangianLevelSetTwoPhaseFoam` (SL two-phase) | `slAdvection::New` | [createSLFields.H L90](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L90) | the production two-phase line |
| `leiaLevelSetTwoPhaseFoam` (Eulerian two-phase) | its own psi equation, interIsoFoam-based | [[concepts/eulerian-solver-mass-flux-port]] | the `baselineEulerian` gate candidate |
| `leiaRedistancedLevelSetFoam` | Eulerian advective form plus criterion-gated redistancing | [createFields.H L64-L72](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaRedistancedLevelSetFoam/createFields.H#L64-L72) | the GRL line |

## Why it matters

The two lines target the same problem: transport psi accurately while it stays usable as a signed distance ([MC article L54-L59](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L53)). The Eulerian line couples the interface accuracy to the flux limiter and the divergence scheme; the SL line has no flux, no divergence and no linear system ([SL article L225-L227](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L216)). Neither line maintains the signed-distance property through maximal stretching: the half-time gradient defect saturates toward O(1) on every mesh ([SL article L1213-L1217](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1200)).

## Evidence

| claim | number | where |
|---|---|---|
| Decision table, T = 8, N = 512: SL shape error and clock | 1.49e-05 at 595 s (cached), 722 s (uncached) | MEASURED, [benchVortex_decision.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex) |
| Decision table, T = 8, N = 512: Eulerian | 3.04e-04 at 1170 s | MEASURED, same table |
| Decision table, T = 8, N = 512: velocity extension closestPoint | 1.44e-03 at 1.61e4 s | MEASURED, same table |
| SL against Eulerian at 512^2, T = 8 | 20x more accurate at half the wall clock; 27x cheaper than velocity extension | MEASURED, [MC article L195-L199](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194) |
| SL shape order, 2D vortex, CFL 0.5 / 1.0 | 2.843 / 2.378 (re-established 2026-08-27 after the gradU fix) | MEASURED, [PCS L52-L61](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61), [sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv) |
| SL shape order, 3D shear hex / poly; 3D deformation hex / poly | 2.955 / 3.285; 1.360 / 1.464 | MEASURED, same table; the SL article quotes 2.95 / 3.28 and 1.36 / 1.46 at [L1253-L1259](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1250) |
| Eulerian shape order at T = 2 / 8 (finest ladder) | 1.52 / 1.18 (fitted slope sign convention of the table: -1.52 / -1.18) | MEASURED, [benchVortex_detailed.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_detailed.tex) |
| Filament death at T = 8 for N <= 64 | method-independent; refinement is the cure | MEASURED, [MC article L204-L205](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194) |

## Decisions

- `SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares` per case, on the transport ladders ([METHOD 8.1 L373](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373)): [[decisions/sl-reconstruction-uncached-qwls]].
- The kinematic token default `ADVECTION eulerian` is historical; the decision-benchmark configs set the line per arm ([DP L250-L253](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L250-L253)).

## Open questions

1. The 2D orders of METHOD 8.3.7 were 3/2 too high until 2026-09-26; the curated `advConv2D*` CSVs are not yet regenerated ([STATUS L3237-L3242](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3237-L3242)).
2. The method comparison ran at `np = 4` before the 2026-08-26 gradU fix; the 2D vortex row was re-established, the 3D rows remain pending ([PCS L52-L61](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61)).

## Related

[[hubs/advection]] - [[models/sl-reconstruction]] - [[models/sl-scheme]] - [[models/sdpls-source]] - [[models/velocity-extension]] - [[models/redistancer]] - [[concepts/eulerian-fv-transport]] - [[concepts/linear-semi-lagrangian]] - [[concepts/advection-regression-set]] - [[cases/kinematic-advection-cases]] - [[studies/method-comparison]] - [[studies/sl-quadratic-pre-print]] - [[retractions/gradu-coupled-patch-contamination]]

## Log

### 2026-09-28
Created from the family headers, the solver composition roots and the method comparison article.
