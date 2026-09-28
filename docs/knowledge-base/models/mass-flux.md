---
title: "massFlux: interpolatedDensity, geometricFaceDensity, rhoLENT"
description: "The two-phase mass-flux family of both coupled solvers and its four sub-switches; rhoLENT is production, measured on the stationary droplet only; every translating measurement before 2026-09-02 is void"
kind: model
status: settled
part: mass-flux
tags: [model, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-02
decided_by: [config/rhoLENTStationary2D.yaml, author decision 2026-09-02]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H, cases/default.parameter]
sources: [METHOD 6, METHOD 8.1 row MASS_FLUX, STATUS 0, STATUS 11.13, DP 624-869, G2 65-84]
---
# massFlux: interpolatedDensity, geometricFaceDensity, rhoLENT

> Verdict (2026-09-28). The dictionary word `levelSet.massFlux.type` selects one of three ways to build the mass flux `rhoPhi` of `div(rhoPhi,U)` ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L8-L31)). `rhoLENT` is the production value since 82ca995 (2026-09-02) ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L657-L685)). Its measured basis is the stationary droplet: `rhoLENTStationary2D` moved the unabsorbed capillary residual by +1.0 / -22 / +0.1 % at N = 32 / 64 / 128, with volume and shape equal to three digits ([config result](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L45-L80), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L364-L401)). Every translating measurement of this family before 440107f (2026-09-02) ran on a closed box and is void ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L73), [[retractions/closed-box-translating-droplet]]). After the fix, dec002f falsified mass-momentum consistency as the dominant term of the translating instability ([[concepts/density-ratio-amplifier]], [[retractions/mass-momentum-consistency-dominant-term]]). The sub-switch `boundRho` is active with rhoLENT and its only measured basis is void ([[concepts/bound-rho]]).

## What it is

The one-field momentum equation carries the mass flux `rhoPhi` in `div(rhoPhi,U)`. The family decides how `rhoPhi` is built from the level set and the volumetric flux `phi`:

- `interpolatedDensity`: `rho` is the geometric cell density, `rhoPhi = interpolate(rho) * phi`.
- `geometricFaceDensity`: `rho_f = alpha_f rho1 + (1 - alpha_f) rho2` from a face area fraction `alpha_f` of the reconstructed plane, `rhoPhi = rho_f * phi`. No density equation is solved.
- `rhoLENT`: the same geometric `rho_f * phi`, plus an auxiliary density equation on every PIMPLE outer corrector, and a reset of `rho` from the phase indicator after the loop ([[concepts/rholent-mass-flux]]).

The face area fraction comes from the linear least-squares plane of a cell, cut against the face and filled per triangle ([`faceAreaFraction.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H#L1-L18)). Since c094bd8 (2026-09-27) the Eulerian two-phase solver reads the same dictionary through the same three headers ([[concepts/eulerian-solver-mass-flux-port]]).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| interpolatedDensity | `interpolatedDensity` | open (unmeasured) | the solver's own code default; no `config/*.yaml` and no `cases/*.parameter` selects it at 8867581 (only the token comment names it) | [`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L17-L18) |
| geometricFaceDensity | `geometricFaceDensity` | settled as the ablation | relative mass residual 0.04 to 0.56 on the stationary droplet, no harm there; the shared default until 82ca995 | [`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L45-L80), [`default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L645-L657) |
| rhoLENT | `rhoLENT` | settled (production) | neutral to better on the stationary droplet, relative residual about 1e-13; no translating benefit shown after the fix | [METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L320), [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3578-L3595) |

## Sub-switches

| switch | dictionary word and values | default (code / token) | status | verdict in one line | evidence |
|---|---|---|---|---|---|
| alpha_f source | `alphaFSource`: `donorPlane`, `averagedPlanes`, `donorPlaneAdvected` (deprecated `upwind`, `central`) | `donorPlane` / `MASS_FLUX_ALPHAF_SOURCE donorPlane` | settled by author instruction, unmeasured | "only Gauss upwind worked in rhoLENT"; the alternatives were measured only on the closed box | [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L33-L92), [DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L686-L699), [[concepts/alphaf-source-donor-plane]] |
| density bound | `boundRho`: `true`, `false` | `false` / `MASS_FLUX_BOUND_RHO true` | open | active with rhoLENT; the only measured basis (rho to -72.28 without it) is void | [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L117-L135), [`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L19-L61), [[concepts/bound-rho]] |
| time level of alpha_f | `alphaFTimeLevel`: `new`, `trapezoid` | `new` / `MASS_FLUX_ALPHAF_TIME_LEVEL new` | voided (measurement) | trapezoid was "marginal" on the closed box: relative residual 0.17 to 0.29 against 0.17 to 0.22 | [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L137-L165), [DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L817-L837), [config](https://github.com/leia-openfoam/leia/blob/8867581/config/alphaFTimeLevelTranslating2D.yaml#L50-L85) |
| flux projection | `projectMassFlux`: `true`, `false` | `false` / `MASS_FLUX_PROJECT false` | voided (measurement) | on the closed box: residual down 24x, travelled fraction 0.725 to 0.049 | [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L167-L213), [[concepts/mass-flux-projection]] |
| residual diagnostic | `massResidualDiagnostic`: `true`, `false` | `false` / `MASS_RESIDUAL_DIAGNOSTIC false` | settled (diagnostic) | computes `R = ddt(rho) + div(rhoPhi)` for every model; changes no field | [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L127-L135), [`massResidualDiag.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massResidualDiag.H#L1-L24) |

Two readings of the code against the token comments, DERIVED at 8867581. The code default of the type is `interpolatedDensity` and the token default is `rhoLENT`. The token comment says the solver defaults rhoLENT to `central` ([DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L686-L690)); the code defaults `alphaFSource` to `donorPlane` ([code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L59-L63)). The comment is stale.

## Why it matters

For a droplet that translates at a uniform `U0`, the momentum transient plus convection collapses to `U0 * R` with `R = ddt(rho) + div(rhoPhi)`. That source is not a gradient, so the pressure projection cannot absorb it, and it vanishes at `U0 = 0` ([`massResidualDiag.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massResidualDiag.H#L11-L16), DERIVED). The family controls `R`. The measurement of 2026-09-02 then showed that `R` is not what sets the translating instability at the water/air ratio ([[concepts/density-ratio-amplifier]]).

## Where in the code

1. Selection, switches and the face fields `alphaf`, `alphafOld`, `rhof`: [`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L1-L279).
2. The face density from the new interface, and the cell density for the models without a density equation: [`updateFaceDensity.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H#L1-L78).
3. `rhoPhi`, the optional projection, the auxiliary equation or the diagnostic, on every outer corrector: [`updateMassFlux.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H#L1-L37).
4. The auxiliary equation, the clip and the residual: [`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L1-L123); the reset: [`resetRhoLENT.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H#L1-L44).
5. The SL solver includes them in [`slAlphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L470-L473) and [`createTransportFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L254); the Eulerian solver in [`leiaLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/leiaLevelSetTwoPhaseFoam.C#L103-L109) and [`alphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H#L95-L127).
6. The gate pins the family in the `twoPhaseCoupling` block: [`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L65-L84).

## Evidence

| claim | number | where |
|---|---|---|
| rhoLENT is neutral to better on the stationary droplet | +1.0 / -22 / +0.1 % on the unabsorbed capillary residual at N = 32 / 64 / 128; volume and shape equal to three digits; all six arms complete | [`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L45-L80), MEASURED |
| the relative mass residual of rhoLENT | about 1e-13 up to every crash | [METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L320), MEASURED |
| the residual of geometricFaceDensity on the stationary droplet, and no harm there | 5.072e-01 / 3.955e-02 / 5.590e-01 relative at N = 32 / 64 / 128 | [`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L45-L80), MEASURED |
| the clip of rhoLENT on the stationary droplet is at round-off | clipL1 2.13e-12 to 2.86e-12 | same, MEASURED |
| nine orders in the residual move the velocity excess by less than 2x | 1.68e-11 to 1.57e-01 (limitedLinearV) and 757x (upwind) in `R`; `L1(U-U0)` moves by less than a factor of two | [SL article, sec. translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2028) `sec:translating`, [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), MEASURED |
| the two models coincide at density ratio 1 | max relative difference 0.000e+00 in all five metrics of all four scheme pairs | [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), MEASURED |
| the translating advantage of rhoLENT (travelled fraction 0.874 against 0.725, volume error 0.0024 against 0.0598 at N = 128) | closed box before 440107f | [`massFluxComparison2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L49-L91), VOID |
| the oscillating droplet at N = 128: rhoLENT completes where geometricFaceDensity reaches 98.2 %, shape 1.26e-5 against 2.51e-4 | one rung; the case used an algebraic psi at that time | [DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L652-L657), [STATUS 11.4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3271-L3280), MEASURED with the void decision open |

## Decisions

- `MASS_FLUX rhoLENT`, author decision 2026-09-02 (82ca995): [[decisions/mass-flux-rholent]].
- `MASS_FLUX_ALPHAF_SOURCE donorPlane`, author instruction 2026-09-01: [[decisions/mass-flux-alphaf-donor-plane]].
- `MASS_FLUX_BOUND_RHO true`, open: [[decisions/mass-flux-bound-rho]].

## Open questions

1. No valid measurement shows a translating benefit of rhoLENT after 440107f ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3578-L3595)).
2. The basis of `boundRho` ([[concepts/bound-rho]]).
3. `interpolatedDensity` has never been measured.
4. Every SL two-phase result on more than one rank before 28d13f0 carries the coupled-face density error; which studies to re-run is an author decision ([[concepts/coupled-face-density-defect]]).
5. The stale token comment on the code default of `alphaFSource` (above).

## Related

- Hub: [[hubs/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/bound-rho]], [[concepts/alphaf-source-donor-plane]], [[concepts/mass-flux-projection]], [[concepts/ddt-scheme-pairing-bdf2]], [[concepts/eulerian-solver-mass-flux-port]], [[concepts/coupled-face-density-defect]], [[concepts/density-ratio-amplifier]].
- Cases: [[cases/translating-droplet]], [[cases/stationary-droplet]].
- Decisions and retractions: [[decisions/mass-flux-rholent]], [[decisions/mass-flux-bound-rho]], [[retractions/closed-box-translating-droplet]], [[retractions/mass-momentum-consistency-dominant-term]].

## Log

### 2026-09-28
Created.
