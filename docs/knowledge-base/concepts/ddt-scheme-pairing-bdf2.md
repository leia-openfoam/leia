---
title: "Pairing the ddt schemes: BDF2 everywhere"
description: "ddt(rho,U) and ddt(rho) use backward (BDF2) in every study; the basis is the matching argument, because the pairing tables ran on the closed box; the token flipped three times in one day"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-08-20
decided_by: [author decision 2026-08-20, config/ddtOrderGain3D.yaml]
code: [cases/default.parameter, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massResidualDiag.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H]
sources: [CLAUDE BDF2 section, METHOD 8.1 rows MOMENTUM_DDT_SCHEME and RHO_DDT_SCHEME, DP 314-354, DP 871-948, STATUS 11.13]
---
# Pairing the ddt schemes: BDF2 everywhere

> Verdict (2026-09-28). Momentum uses OpenFOAM's `backward` (BDF2) in every leia two-phase study, never Euler, since 2026-08-20 on the author's direction ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L248-L262), [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L314-L338)). The rhoLENT density equation uses the same scheme, `RHO_DDT_SCHEME backward`, and the basis is an argument, not a table: with `U == U0` the momentum transient reduces to `U0` times the mass residual only when `ddt(rho,U)` and `ddt(rho)` expand with the same coefficients ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L874-L883), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L372), DERIVED). Both pairing tables (`matchedBDF2Translating2D`, `volumeCorrectionTranslating2D`) ran on the closed-box translating case and are void ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L942-L947), [[retractions/closed-box-translating-droplet]]). The measured basis of BDF2 itself is the stationary droplet: BDF2 against Euler on matched windows moved the gain by +11.1 / +2.9 / -3.0 %, a sign-flipping difference, with volume and shape within 1.2 % ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L371)). BDF2 costs nothing and is formally right, so it is mandatory ([[decisions/momentum-schemes-bdf2-upwind]]).

## What it is

Two `fvSchemes` keys carry a time scheme: `ddt(rho,U)` (token `MOMENTUM_DDT_SCHEME`) and `ddt(rho)` (token `RHO_DDT_SCHEME`). The rhoLENT auxiliary equation `fvm::ddt(rho) + fvc::div(rhoPhi) = 0` reads the second key ([`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L13-L17)), and the solver banner prints both so the pairing is never in doubt ([`massResidualDiag.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massResidualDiag.H#L20-L24)).

The matching argument. For a droplet that translates at uniform `U0`, the momentum transient plus convection collapses to `U0 * [ddt(rho) + div(rhoPhi)]` only when the two `ddt` operators agree. Otherwise momentum inherits `U0 * [ddt_mom(rho) - ddt_mass(rho)]`, which for a backward/Euler pairing is `U0 * (1/2)(rho^{n+1} - 2 rho^n + rho^{n-1})/dt`: half the second difference of `rho` over `dt`, of order `U0 * drho/dt` in every cell the interface sweeps ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L874-L883), DERIVED). Euler/Euler and backward/backward both satisfy it; the repository keeps second order.

Why BDF2 in momentum. The semi-Lagrangian foot trace is second order in time and must not be fed a first-order velocity; BDF2 evaluates fluxes and sources at `t^{n+1}`, which matches the capillary force built from `psi^{n+1}` ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L248-L262)). A genuinely second-order configuration also needs `PSI_OUTER_CORRECTORS yes`, so that the force sits at `t^{n+1}` and not at an extrapolated trajectory ([`ddtOrderGain3D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/ddtOrderGain3D.yaml#L1-L28), [[concepts/psi-outer-correctors]], [[retractions/force-at-n-not-n-plus-1]]). The rhoLENT paper solves its density update explicitly and uses Euler momentum; leia keeps `backward` on both because only the pairing matters for the consistency ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L910-L929)).

## Why it matters

The token `RHO_DDT_SCHEME` flipped three times on 2026-09-01 and 2026-09-02: it had no entry and was automatically matched; f7307b5 set it to `Euler` for boundedness and broke the matching; 28a1383 restored `backward`; b60e3df set `Euler` again on the paper's algorithm and broke it a second time ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L931-L936)). The record is kept because the flip-flop is itself the lesson: change `MOMENTUM_DDT_SCHEME` and `RHO_DDT_SCHEME` together, and read the pairing argument before either.

## Evidence

| claim | number | where |
|---|---|---|
| BDF2 against Euler on the stationary droplet | gain +11.1 / +2.9 / -3.0 % (noise); volume and shape within 1.2 % | [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L371), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L257-L262), MEASURED; the config is not named in either source |
| the 2x2 of ddt order and force placement, 3D R/h = 15.8 | `backward` no longer diverges, peak 6 % below Euler's, volume error 1.092e-05 to 8.744e-06 (20 % better), shape 0.8 % better, per-step gain unchanged (+1.1 %) | [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L321-L328), [`ddtOrderGain3D`](https://github.com/leia-openfoam/leia/blob/8867581/config/ddtOrderGain3D.yaml#L38-L46), MEASURED |
| the pairing table, translating N = 128 | Euler/Euler 13334 steps, travel 0.7667; Euler/backward 13334, 0.8827; backward/Euler DIVERGED at 10431, travel 1.3210; backward/backward 13334, 0.8737, volume error 0.0024 | [`matchedBDF2Translating2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/matchedBDF2Translating2D.yaml#L68-L92), VOID (closed box) |
| `MOMENTUM_DIV_SCHEME upwind` is inert on the stationary droplet | max U 8.0915e-05 against 8.0944e-05 (R/h = 12.7) and 2.4093e-04 against 2.4084e-04 (R/h = 15.8) with `linearUpwind gradU`; per-step gain to 0.03 % | [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L340-L354), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L396), MEASURED |
| the convection schemes on the repaired translating droplet | upwind, limitedLinearV, vanLeerV, linearUpwind within 30 % at step 5000, ordered by numerical diffusion; upwind is the longest-lived arm (9987 steps) | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2029), [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3578-L3595), MEASURED |
| the second-order foot integrator and the force centring on the late instability | `rk2` inside the scatter (0.0842 s against 0.0868 s); `midpoint` 32 % earlier (0.0593 s) | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3744), MEASURED, one resolution |

## Decisions

- `MOMENTUM_DDT_SCHEME backward`, `RHO_DDT_SCHEME backward`, `MOMENTUM_DIV_SCHEME upwind`: [[decisions/momentum-schemes-bdf2-upwind]]. Every new case template references the tokens; no template hardcodes a scheme ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L248-L262)).
- `CAPILLARY_FORCE_CENTRING endStep`: spectrally identical to `midpoint` for the linear capillary oscillator; `midpoint` diverged 32 % earlier on the translating droplet ([[concepts/force-time-centring]]).

## Open questions

1. The config behind the "+11.1 / +2.9 / -3.0 %" windows is not named in METHOD 8.1 or CLAUDE.md.
2. A valid pairing table on the repaired translating case does not exist. The argument is derived; its size on the repaired case is not measured.

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/bound-rho]], [[concepts/density-ratio-amplifier]], [[concepts/psi-outer-correctors]], [[concepts/force-time-centring]].
- Decisions and retractions: [[decisions/momentum-schemes-bdf2-upwind]], [[retractions/force-at-n-not-n-plus-1]], [[retractions/closed-box-translating-droplet]].

## Log

### 2026-09-28
Created.
