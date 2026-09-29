---
title: "MOMENTUM_DDT_SCHEME backward, RHO_DDT_SCHEME backward, MOMENTUM_DIV_SCHEME upwind"
description: "BDF2 on the momentum and on the rhoLENT density equation, upwind momentum convection, defaults since 2026-08-20: BDF2 against Euler moves the stationary-droplet gain +11.1 / +2.9 / -3.0 % (noise) at no cost, the density scheme must match by the free-stream argument, and upwind is inert on the stationary droplet to four significant figures"
aliases: [MOMENTUM_DDT_SCHEME backward, RHO_DDT_SCHEME backward, MOMENTUM_DIV_SCHEME upwind]
kind: decision
status: settled
part: mass-flux
tags: [decision, part/mass-flux]
date: 2026-09-28
date_settled: 2026-08-20
decided_by: [config/upwindConvection2D.yaml, config/upwindConvection3D.yaml, config/cellCentreInverseFilteredPostFix.yaml, "author decision 2026-08-20"]
code: [cases/default.parameter, config/gates/methodGate2D.yaml, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H]
sources: ["METHOD 8.1 rows MOMENTUM_DDT_SCHEME, RHO_DDT_SCHEME, MOMENTUM_DIV_SCHEME (L368-L369, L396)", "CLAUDE BDF2 section (L292-L306)", "STATUS 4 (L1034-L1046)", "STATUS 11.13 (L3599-L3600)", "PSH L354-L372", "DP L340-L356 and L940-L955"]
---
# MOMENTUM_DDT_SCHEME backward, RHO_DDT_SCHEME backward, MOMENTUM_DIV_SCHEME upwind

> Three tokens of the momentum discretisation, all in the global default since 2026-08-20 ([DP L340-L356](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L340-L356)) and in the gate's coupling block ([methodGate2D.yaml L73-L75](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L73-L75)). `MOMENTUM_DDT_SCHEME backward` (BDF2) is a repository rule: the semi-Lagrangian foot-point trace is second order in time and must not be fed a first-order velocity ([CLAUDE L292-L306](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L292-L306)). Its measured basis is a null result: on matched windows of 1179 / 3333 / 9428 steps BDF2 against Euler moves the stationary-droplet per-step gain by +11.1 / +2.9 / -3.0 %, sign-flipping, with volume and shape within 1.2 %; BDF2 costs nothing and is formally right ([STATUS L1039-L1046](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1039-L1046), [METHOD 8.1 L368](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L368)). `RHO_DDT_SCHEME backward` follows by the matching argument: with `U = U0` the momentum transient reduces to `U0` times the mass residual only if `ddt(rho, U)` and `ddt(rho)` use the same scheme; the pairing tables that once decided it ran on the closed box and are void ([DP L949-L954](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L949-L954), [METHOD 8.1 L369](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L369)). `MOMENTUM_DIV_SCHEME upwind` was set at the user's direction on 2026-08-20 and is inert on the stationary droplet: max abs U agrees with `linearUpwind gradU` to four significant figures on the 3D droplet at R/h = 12.7 and 15.8 and to 4 to 5 figures over the 2D ladder ([DP L343-L356](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L343-L356), [STATUS L1034-L1038](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1034-L1038)).

## The question

Three questions, one per token. Is the momentum time scheme the amplifier of the parasitic current? Which scheme must the rhoLENT density equation use? Does `linearUpwind`, which extrapolates the face velocity with the cell gradient of `U`, act as an amplifier at the interface, where `U` is continuous but its gradient is not ([upwindConvection2D.yaml L1-L17](https://github.com/leia-openfoam/leia/blob/8867581/config/upwindConvection2D.yaml#L1-L17))? A second question on BDF2 is placement: BDF2 is second order only when the force is evaluated at `t^{n+1}`, which needs the interface re-advected on every outer corrector ([ddtOrderGain3D.yaml L1-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/ddtOrderGain3D.yaml#L1-L30), [[concepts/psi-outer-correctors]], [[retractions/force-at-n-not-n-plus-1]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| 2D stationary ladder, N = 64 / 128 / 256, matched windows 1179 / 3333 / 9428 steps: BDF2 with upwind against Euler with `linearUpwind` | per-step gain `g`; delta g; delta max abs U; delta volume; delta shape | -2.312e-03 / +3.899e-04 / +2.661e-04 against -2.600e-03 / +3.790e-04 / +2.744e-04; +11.1 / +2.9 / -3.0 %; +40.5 / +3.7 / -7.5 %; +1.1 / +0.4 / +0.4 %; +1.2 / +0.4 / -0.2 % | MEASURED, [PSH L362-L372](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L362-L372) |
| an unmatched comparison of the same arms | gain at N = 512 | BDF2 looked 3.1x worse; that was the horizon: never compare `gAvg` across unequal step counts | RETRACTED, [STATUS L1043-L1046](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1043-L1046) |
| 2D ladder, N = 64 to 512, BDF2: `upwind` against `linearUpwind gradU` | `u0`; max abs U; volume and shape | identical to 5 digits; within 0.8 % at N = 64 and 0.06 % for N >= 128; within 0.01 % | MEASURED, [STATUS L1034-L1038](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1034-L1038) |
| 3D stationary droplet, R/h = 12.7 and 15.8, BDF2 | max abs U | 8.0915e-05 against 8.0944e-05; 2.4093e-04 against 2.4084e-04; per-step gain within 0.03 % | MEASURED, [DP L346-L349](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L346-L349) |
| `linearUpwind gradU` against the historical `linearUpwind grad(U)` | metrics CSV | byte-equal | MEASURED, [upwindConvection3D.yaml L54-L56](https://github.com/leia-openfoam/leia/blob/8867581/config/upwindConvection3D.yaml#L54-L56) |
| the ddt pairing tables (`matchedBDF2Translating2D`) and the densities -27.7 and -1786 without the bound | — | VOID (closed box) | [DP L949-L954](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L949-L954) |
| `translatingRepaired2D`, N = 128, four convection schemes | completion | every arm diverges; upwind is the longest-lived (9987 steps) | MEASURED, [STATUS L3600](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3600) |

Pre-registered read-outs: [upwindConvection2D.yaml L19-L38](https://github.com/leia-openfoam/leia/blob/8867581/config/upwindConvection2D.yaml#L19-L38) (outcome c occurred: no change, grad U is not the amplifier on the stationary droplet), [ddtOrderGain3D.yaml L58-L69](https://github.com/leia-openfoam/leia/blob/8867581/config/ddtOrderGain3D.yaml#L58-L69).

Both `upwindConvection` studies and `ddtOrderGain3D` ran with the psi filter at theta = 0.2 ([upwindConvection2D.yaml L56-L57](https://github.com/leia-openfoam/leia/blob/8867581/config/upwindConvection2D.yaml#L56-L57)); the verdicts are inertness verdicts on a filtered baseline ([[decisions/psi-filter-none]]).

## What it does not cover

1. The translating and oscillating droplets. The convective term is negligible on the stationary droplet, so the inertness of `upwind` says nothing there; every translating arm diverges late whatever the scheme ([DP L350-L355](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L350-L355), [STATUS L3600](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3600)). The repaired matrix found the convection scheme not to be the mechanism of the late instability ([SL article L2196](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2196), [[concepts/density-ratio-amplifier]]).
2. The density scheme. BDF2's homogeneous update on the density equation is an extrapolation, so the bound `MASS_FLUX_BOUND_RHO` is active with it; its evidence is void ([[decisions/mass-flux-bound-rho]]). The scheme flip-flopped three times (f7307b5, 28a1383, b60e3df) ([DP L940-L948](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L940-L948), [[concepts/ddt-scheme-pairing-bdf2]]).
3. The force centring is a separate token: `endStep`; `midpoint` is spectrally identical for the linear oscillator and diverges 32 % earlier on the translating droplet ([[concepts/force-time-centring]]).
4. BDF2 starts on Euler at the first step ([plan-shannon L372-L373](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-shannon-parasitic-currents.md#L372-L373)).

## Related

[[hubs/mass-flux]] - [[hubs/surface-tension]] - [[concepts/ddt-scheme-pairing-bdf2]] - [[concepts/psi-outer-correctors]] - [[concepts/force-time-centring]] - [[concepts/density-ratio-amplifier]] - [[concepts/parasitic-current-mechanism]] - [[decisions/mass-flux-rholent]] - [[decisions/mass-flux-bound-rho]] - [[decisions/mass-flux-alphaf-donor-plane]] - [[decisions/psi-filter-none]] - [[retractions/force-at-n-not-n-plus-1]] - [[retractions/closed-box-translating-droplet]] - [[cases/stationary-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED: BDF2 and upwind as defaults on 2026-08-20 (86cba7c, d96a9d0); the density scheme by the matching argument, its tables void since 2026-09-27. Entered in [[decision-log#2026-08]].
