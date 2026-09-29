---
title: "viscosityFaceModel: the six face viscosities"
description: "The face dynamic viscosity of the viscous term; alg_lin decided 2026-09-03 on the 36-arm ladder as the only model with positive orders, after a frozen-muf bug had briefly pointed at geo_lin (2026-09-28)."
aliases: [viscosityFaceModel, VISCOSITY_FACE_MODEL, face viscosity]
kind: model
status: settled
part: viscosity
tags: [model, part/viscosity]
date: 2026-09-28
date_settled: 2026-09-03
decided_by: [config/mufGrid2D_jump1000.yaml, config/mufGrid2D_jump1000_N128.yaml, config/mufGrid2D_jump1000_N256.yaml, config/mufGrid2D_prod_N128.yaml, config/mufGrid2D_prod_N256.yaml, config/mufGrid2D_prod_N512.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceViscosity.H, applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H, applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H]
sources: [DP VISCOSITY_FACE_MODEL, METHOD 8.1 row VISCOSITY_FACE_MODEL, STATUS 0 Popinet, STATUS 11.13, STATUS 11.14, SL article sec:viscous]
---
# viscosityFaceModel: the six face viscosities

> Verdict (2026-09-28). `VISCOSITY_FACE_MODEL alg_lin` is the default, decided 2026-09-03 on a 36-arm ladder: six models, N = 128 / 256 / 512, viscosity ratios 54.83 and 1000 ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L729-L763)). `alg_lin` is the only model with a positive order in all three metrics at both ratios; at ratio 1000 its L1 is 1.400e-3 at N = 512 with order 1.10 ([SL article, choosing mu_f](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L830-L873)). The brief switch to `geo_lin` is retracted: the algebraic arms had run with a face viscosity frozen at t = 0 ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L754-L759)). Two items stay open: the 3D droplet templates carry no token and run the solver default `geo_lin`, and the ladder ran at np 4 to 16 before the coupled-face fixes of 2026-09-27 ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L397), [[concepts/viscosity-open-items]]).

## What it is

The viscous term is assembled explicitly per cell with one face coefficient `mu_f` in the implicit Laplacian and in the explicit transpose part ([SL article sec:viscous](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L759-L793)). `mu_f` is a model choice with two independent axes, spelled out as one token with six values ([`createFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L61-L111)):

| axis | words | meaning |
|---|---|---|
| where `alpha_f` comes from | `alg_` / `geo_` | linear interpolation of the cell indicator, or the geometric face area fraction that also builds `rho_f` |
| how `mu_f` is mixed | `_lin` / `_harm` / `_blend` | `alpha_f mu1 + (1 - alpha_f) mu2`; `1/mu_f = alpha_f/mu1 + (1 - alpha_f)/mu2`; `w mu_harm + (1 - w) mu_lin` with `w = abs(n . S_f)/abs(S_f)` |

The legacy words `interpolated`, `geometric` and `geometricHarmonic` still parse and map to `alg_lin`, `geo_lin` and `geo_harm` ([`createFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L93-L96)). `alg_lin` is the framework default, because the linear model commutes with interpolation ([`faceViscosity.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceViscosity.H#L9-L18)).

## Members

Ladder values at ratio 1000 from the token comment; `p` is the least-squares order of `log(error)` against `log h` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L729-L752)).

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| algebraic alpha_f, linear mixing | `alg_lin` | settled, production | L1 1.400e-3 at N = 512, p(L1) 1.10 (R 0.998), p(shape) +0.71; lowest error and travelled fraction closest to 1 at both ratios | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L752), [SL article tab:muf](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L854-L873) |
| algebraic alpha_f, harmonic mixing | `alg_harm` | retracted | 6.98e-3 at ratio 1000, p(L1) 0.08, p(shape) -0.50 | [SL article tab:muf](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L854-L873) |
| algebraic alpha_f, normal blend | `alg_blend` | retracted | 3.09e-3, p(L1) 0.83, p(shape) +0.60: second at the large ratio, negative shape order at ratio 54.8 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L747) |
| geometric alpha_f, linear mixing | `geo_lin` | retracted as default, solver default in code | 1.83e-3, p(L1) 0.46, p(shape) -0.37; the frozen-muf bug made `alg_lin` look worse than this | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L754-L759), [`createFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L92) |
| geometric alpha_f, harmonic mixing | `geo_harm` | retracted | 4.67e-3, p(L1) 0.15, p(shape) -0.56 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L747) |
| geometric alpha_f, normal blend | `geo_blend` | retracted | 5.67e-3 at ratio 54.8, p(shape) -0.65; 2.65e-3 at ratio 1000, order 0.33 (the SL article's table, [L871](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L871); the token comment does not carry it) | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L735-L740) |

## Why it matters

The explicit viscous form is not a refactoring: on the translating droplet at N = 128 the product-interpolated framework form and the explicit form separate by 1.7 percent in L1 after one step and by 36 percent after 2667 steps ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L794-L802)). The explicit form improves the disturbance by 0.64x in L1 and 0.81x in L2, the shape error by 3.3x, and moves the travelled fraction from 1.0125 to 0.9971 ([SL article tab:viscous](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L804-L828)). Popinet's face properties use a simple average of the cell values, which is `alg_lin` ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L141-L143)). Viscosity spreads the spurious current rather than damping it: removing it lowers L1 by 1.42x and raises L2 by 2.61x ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L875-L884)); in Popinet's constant-property case the inviscid limit is the worst arm ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L164-L169)).

## Where in the code

1. The dictionary word is read from `levelSet.massFlux.viscosityFaceModel` with the default `geo_lin` in [`createFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L89-L111); the SL solver includes this file from its sibling ([`leiaSemiLagrangianLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/leiaSemiLagrangianLevelSetTwoPhaseFoam.C#L136)).
2. `mu_f` is rebuilt every step for all six models in [`faceViscosity.H`](https://github.com/tmaric/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceViscosity.H#L30-L88).
3. The banner prints the default as `alg_lin (default)` while the code default is `geo_lin` ([`printMethodBanner.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H#L139), [`createFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L92)). A 3D case without the token prints a wrong default.
4. The 2D templates render the token ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L294-L299)); the three 3D droplet templates contain no `viscosityFaceModel` entry.
5. The Eulerian two-phase solver kept `muf` at its t = 0 value until c094bd8 ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3681-L3686)).

## Evidence

| claim | number | where |
|---|---|---|
| the 36-arm ladder at ratio 54.83, L1 at N = 512 | alg_lin 3.494e-3 (p 0.77); geo_lin 4.953e-3 (0.15); alg_harm 6.479e-3; geo_harm 6.027e-3; alg_blend 6.247e-3; geo_blend 5.669e-3 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L733-L740), MEASURED |
| the ladder at ratio 1000 | alg_lin 1.400e-3 (p 1.10, R 0.998); alg_blend 3.086e-3; geo_lin 1.834e-3; geo_harm 4.670e-3 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L747), MEASURED |
| every model except alg_lin, and alg_blend at the large ratio, has a negative shape order | see the two tables | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L837-L843), MEASURED |
| the retracted switch to geo_lin | muf constructed once and updated only in an `if (geometric)` branch; fixed in 39e59b3 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L754-L759), MEASURED |
| the ladder ran in parallel before the coupled-face fixes | np 8 (jump1000, N256), np 4 (N128), np 16 (prod_N512) | [`config/mufGrid2D_jump1000.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/mufGrid2D_jump1000.yaml#L36-L40), [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3668-L3671), MEASURED |
| the viscous-term formulation on the translating droplet, N = 128, t = 0.02 s | product form L1 1.098e-2; explicit form with interpolated mu_f 7.013e-3, geometric 6.113e-3, harmonic 6.728e-3 | [SL article tab:viscous](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L804-L823), MEASURED |

## Why it failed, or why we think so

The harmonic mixing reproduces the exact viscous flux across a layered medium normal to the interface, and the blend interpolates between the series and parallel limits without a free parameter. Both are among the worst at every resolution and both ratios. Neither argument describes what this discretisation requires; the plain arithmetic mixing of a linearly interpolated `alpha_f` outperforms both, and consistency of `mu_f` with `rho_f` is not what limits the term ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L845-L852)). Pure harmonic applies the series formula even on faces that lie nearly in the interface plane, which is a candidate reason ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L724-L728)).

## Decisions

- `VISCOSITY_FACE_MODEL alg_lin`, default layer, 2D templates; see [[decisions/viscosity-face-model-alg-lin]].
- The token was split from one three-valued word into the six-valued grid in 39e59b3 so that no combination can hide ([`config/mufGrid2D_jump1000.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/mufGrid2D_jump1000.yaml#L1-L28)).

## Open questions

See [[concepts/viscosity-open-items]]: the 3D token, the banner default, the parallel provenance of the ladder, and the geometric `alpha_f` consistency question.

## Related

[[hubs/viscosity]], [[concepts/viscosity-open-items]], [[decisions/viscosity-face-model-alg-lin]], [[concepts/rholent-mass-flux]], [[concepts/eulerian-solver-mass-flux-port]], [[concepts/coupled-face-density-defect]], [[cases/translating-droplet]], [[cases/popinet-translating-droplet]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from the token comment, METHOD 8.1, the SL article section on the viscous term and the solver sources.

### 2026-09-29
CORRECTED: the geo_blend value at ratio 1000 is recorded in the SL article's table (2.65e-3, order 0.33).
