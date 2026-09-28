---
title: "slScheme: pointValue, fluxForm, normalProjected"
description: "The production scheme is pointValue with the taylor foot and the input trajectory velocity; fluxForm conserves the integral of psi but loses volume on long horizons; normalProjected is closed."
aliases: [slScheme, SL_SCHEME, SL_FOOT_INTEGRATOR]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-27
decided_by: [config/benchVortexSLimproved.yaml, config/npslConv2Dvortex.yaml, config/transISTKinematicRefinementGate.yaml]
code: [src/leiaLevelSet/semiLagrangian/slScheme.H, src/leiaLevelSet/semiLagrangian/pointValueScheme.C, src/leiaLevelSet/semiLagrangian/fluxFormScheme.C, src/leiaLevelSet/semiLagrangian/normalProjectedScheme.C]
sources: ["METHOD 8.1 row SL_FOOT_INTEGRATOR (L376)", "MC article sec:slimp (L258-L295)", "MC article sec:fluxloss (L147-L166)", "nPSL article sec:results (L393-L493)", "RM L136-L155 and L302-L327", "STATUS 11.15 (L3732-L3738)"]
---
# slScheme: pointValue, fluxForm, normalProjected

> Verdict (2026-09-28). The scheme turns `psi^n` into `psi^{n+1}` once the reconstruction and the corrector are fixed ([slScheme.H L30-L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slScheme.H#L30-L49)). `pointValue` is the production scheme and the code default ([slScheme.C L53](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slScheme.C#L53)). Its foot integrator is `taylor`, decided without a gate ([METHOD 8.1 L376](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L376)); `rk2` moved the translating divergence time inside the decomposition scatter ([STATUS L3736](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3736)). `fluxForm` conserves the integral of psi to 2e-14 but has 3 to 13 times the volume error at T = 8 ([MC article L282-L295](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L282-L295)). `normalProjected` is closed: its orders are -1.24 and -0.17 against 2.9 for `pointValue` ([nPSL article L456-L468](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L456-L468)).

## What it is

Every scheme is handed the same per-cell quadratic fit and the same corrector by `slAdvection` ([slScheme.H L48-L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slScheme.H#L48-L49)). The dictionary is `levelSet { semiLagrangian { scheme ...; } }`; the token is `SL_SCHEME` ([DP L263-L266](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L263-L266)). The base class carries `dtScale_`, the fraction of the step the drift covers; the midpoint centring of the capillary force drifts a half step before and after the momentum solve ([slScheme.H L81-L94](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slScheme.H#L81-L94), [[concepts/force-time-centring]]).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| point value | `pointValue` | settled, production | `psi^{n+1}(x_c) = psi^n(x_d)`, the Lagrange interpolation of the reconstructed departure value; non-conservative; the only scheme the value bound reaches. | [pointValueScheme.H L30-L38](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L30-L38), [README L727-L729](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L727-L729) |
| flux form | `fluxForm` | settled, not production | Conserves the discrete integral of psi to 2e-14 for a divergence-free flux at CFL <= 1; competitive at T = 2; 3 to 13x the volume error of `pointValue` at T = 8. | [fluxFormScheme.H L58](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/fluxFormScheme.H#L58), [MC article L282-L295](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L282-L295) |
| normal projected | `normalProjected` | retracted (closed line) | Writes `d_c + delta` along `(u.n)n`; the trace is clean, the write-back diverges x1.7 per 10 steps, the corrugation hypothesis is falsified. | [normalProjectedScheme.H L214](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/normalProjectedScheme.H#L214), [[concepts/normal-projected-sl]] |

The dictionary words are the `TypeName` strings ([pointValue L99](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L99), [fluxForm L58](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/fluxFormScheme.H#L58), [normalProjected L214](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/normalProjectedScheme.H#L214)).

## The keys of pointValue

| key | values | default | what it decides | evidence |
|---|---|---|---|---|
| `footIntegrator` | `taylor`, `rk2` | `taylor` | `taylor` needs `grad(U)` for its dt^2 term; `rk2` uses a second velocity sample at the midpoint and never touches `grad(U)`. Same formal order. | [pointValueScheme.C L63](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L63), [pointValueScheme.H L68-L80](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L68-L80) |
| `trajectoryVelocity` | `input`, `normalProjection`, `normalClosestPoint` | `input` | The velocity inserted into the foot formula: the supplied field, the local projection `(U.n)n`, or the normal speed at a closest point of the zero set. | [pointValueScheme.C L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L60), [L79-L91](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L79-L91) |
| `trajectoryU` | word | `U` | The registry name of the physical velocity for the two normal modes. | [pointValueScheme.C L61](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L61) |

The keys of `normalProjected` are `renormalization` (`geometric`, `strain`, `none`; default `strain`) and `offsetEngine` (`rayRoot`, `footPoint`) ([normalProjectedScheme.C L85-L102](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/normalProjectedScheme.C#L85-L102), tokens [DP L269-L272](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L269-L272)).

## Why it matters

The scheme decides which property the transport keeps. `pointValue` keeps the characteristic interpolation reversible: its errors are dispersive, so the volume oscillates by 1 % and recovers on flow reversal ([MC article L159-L161](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L159-L161)). `fluxForm` evaluates the upwind reconstruction at the face; that dissipation is irreversible and smears the steepened psi in the filament ([MC article L152-L166](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L152-L166)).

## Where in the code

- `slScheme::New` reads `scheme` with default `pointValue` ([slScheme.C L48-L53](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slScheme.C#L48-L53)).
- The Taylor foot and the rk2 midpoint foot: [pointValueScheme.C L363-L429](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L363-L429).
- The coupled solver hands the scheme the AB2 estimate `UextStar` and the old level ([slAlphaEqn.H L116-L127](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L116-L127)); see [[concepts/departure-foot-ab2-centring]].

## Evidence

| claim | number | where |
|---|---|---|
| `fluxForm` conserves the integral of psi | to 2e-14 | MEASURED, [MC article L283-L284](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L283-L284) |
| `fluxForm` shape error at T = 2, N = 256, against `pointValue` | 1.5e-05 (flux) edges the point value, whose shape error stalled from 128^2 to 256^2 | MEASURED, [MC article L286-L288](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L286-L288) |
| `fluxForm` volume error at T = 8, N = 64 | 0.82 against 0.062 for `pointValue` | MEASURED, [MC article L291-L293](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L291-L293) |
| Flux-form volume loss at N = 128, T = 8 | -17 % by t = 6.5; band gradient 0.72 against 0.96 | MEASURED, [MC article L157-L163](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L157-L163) |
| Decision table, uncached + flux at N = 256 | E_geom 1.46e-05 (T = 2), 7.66e-04 (T = 8) | MEASURED, [benchVortex_decision.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex) |
| `normalProjected` shape orders, reversed vortex, Co 1/2 | strain -1.24, geometric + foot point -0.17, against 2.9 for `pointValue` | MEASURED, [nPSL article L456-L468](https://github.com/leia-openfoam/leia/blob/8867581/docs/normal-projected-semi-lagrangian/npsl-article/normalProjectedSemiLagrangian.tex#L456-L468) |
| `trajectoryVelocity normalProjection` runaway time, N = 32, translating droplet | 0.04365 s against 0.01309 s for `normalClosestPoint`; physically invalid by t = 0.03 (max abs U 0.248) | MEASURED, [RM L140-L155](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L140-L155) |
| Fixed-Courant refinement, `input` against `normalProjection` | volume orders 2.47, 2.55 against 1.45, 0.81; zero-set orders 2.55, 1.89 against 1.17, 1.11 | MEASURED, [RM L302-L322](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L302-L322) |
| `footIntegrator rk2` on the translating droplet, N = 100, np 4 | diverged at t = 0.0842 s, inside the scatter of the reference (0.0868 s) | MEASURED, one resolution, [STATUS L3732-L3738](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3732-L3738) |

## Decisions

- `SL_SCHEME pointValue`, `SL_FOOT_INTEGRATOR taylor` ([METHOD 8.1 L376](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L376), [DP L1028](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1028)). No gate decided the integrator; the one-resolution discriminator of 2026-09-27 found no effect.
- `normalProjection` is a diagnostic trajectory, not the production one ([RM L525-L531](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L525-L531)).

## Open questions

1. `rk2` consumes no `grad(U)`, the operator that was wrong on coupled patches for the life of the code ([pointValueScheme.H L72-L79](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L72-L79)). A ladder with `rk2` is not on record.
2. `fluxForm` and `normalProjected` bypass the value bound ([README L727-L729](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L727-L729)).

## Related

[[hubs/advection]] - [[models/sl-reconstruction]] - [[models/sl-value-bound]] - [[concepts/departure-foot-ab2-centring]] - [[concepts/normal-projected-sl]] - [[concepts/eulerian-fv-transport]] - [[concepts/trace-velocity-projected-flux]] - [[concepts/force-time-centring]] - [[studies/method-comparison]] - [[studies/npsl-design]]

## Log

### 2026-09-28
Created from the code, the method comparison article, the nPSL article and the roadmap.
