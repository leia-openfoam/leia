---
title: "Time centring of the capillary force"
description: "The capillary force is built from the level set at the end of the step, which makes the coupling symplectic Euler with unit amplification for omega dt below 2; the midpoint centring is spectrally identical for the linear oscillator and diverged 32 percent earlier on the translating droplet, and the earlier force-at-n anti-damping reading is retracted (2026-09-28)."
aliases: [capillaryForceCentring, endStep, midpoint centring, symplectic Euler coupling, force at n+1]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-27
decided_by: [METHOD 8.1 row CAPILLARY_FORCE_CENTRING, config/gates/methodGate2D.yaml, docs/plan-shannon-parasitic-currents.md section 0e]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/leiaSemiLagrangianLevelSetTwoPhaseFoam.C, applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H]
sources: [PSH 0e retraction, MC deck 9, DP CAPILLARY_FORCE_CENTRING, METHOD 8.1, STATUS 11.13, STATUS 11.14, SL negative deck 1/1, PSH 0d BDF2 matched window, PCS 18.4]
---
# Time centring of the capillary force

> Verdict (2026-09-28). `slAlphaEqn.H` advects `psi` in place and rebuilds the band, the indicator, `alpha` and the whole curvature chain from the advected field before `UEqn.H` calls the force, so the capillary force is built from `psi^{n+1}` ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L538-L542)). For the linearised capillary oscillator that is symplectic Euler: `det(M) = 1`, `|lambda| = 1` for `omega dt < 2` ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L543-L556)). The midpoint centring (Stormer-Verlet, Popinet's staggered arrangement) has the same determinant and the same trace, so it changes the accuracy of the coupling and not its linear stability, its frequency or its step limit ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L558-L564), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L398)). Measured once, on the translating droplet at N = 100 after the parallel fixes, `midpoint` diverged 32 percent earlier than `endStep` (0.0593 against 0.0868 s), an indicator at one resolution ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3732-L3738)). The reading that an explicit force at `t^n` anti-damps at `+omega^2 dt/2` is retracted, with every number derived from it ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530-L536), [[retractions/force-at-n-not-n-plus-1]]). The default is `CAPILLARY_FORCE_CENTRING endStep` with BDF2 momentum ([[decisions/momentum-schemes-bdf2-upwind]]).

## What it is

Write one azimuthal mode of the interface as `eta` and its velocity as `u`, with `omega^2 = sigma k^3/(rho_1 + rho_2)`. Three placements of the force exist ([MC deck 9/2](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/9/2), [MC deck 9/4](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/9/4)):

| scheme | update | force built from | amplification per step |
|---|---|---|---|
| (a) fully explicit | `eta^{n+1} = eta^n + dt u^n`, `u^{n+1} = u^n - dt omega^2 eta^n` | the old interface | `sqrt(1 + (omega dt)^2)`, always above 1 |
| (b) symplectic Euler, `endStep` | `eta^{n+1} = eta^n + dt u^n`, `u^{n+1} = u^n - dt omega^2 eta^{n+1}` | the advected interface | 1 while `omega dt < 2` |
| (c) Stormer-Verlet, `midpoint` | half kick, drift, half kick | the midpoint interface | 1 while `omega dt < 2`, same eigenvalues as (b) |

Computed per step at `omega dt = 0.500 / 1.047 / 1.571 / 2.000 / 3.142`: (a) 1.118034 / 1.447829 / 1.862268 / 2.236068 / 3.297296; (b) and (c) 1.000000 at the first four and 7.743015 at the last ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L550-L556)). Popinet staggers the volume fraction at `n +- 1/2` against the velocity at `n` and evaluates the force at `n + 1/2` ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L437-L440)). The practical explicit wall measured here is `omega_grid dt` of 1.0 to 1.3, against the linear bound 2 and Popinet's `pi` ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L581-L583), [[concepts/capillary-time-step]]).

A symplectic integrator conserves only for a force that derives from a potential. The delivered curvature carries an error that depends on where the interface sits relative to the mesh, so the force is not a gradient of any potential and does net work around each oscillation cycle; that work, not the time level, is what a fix has to drive to zero ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L571-L579), [[concepts/parasitic-current-mechanism]], [[concepts/variational-capillary-force]]).

## Why it matters

The momentum time scheme is part of the stabilisation, not only of the accuracy budget: with the Euler-era pipeline, `backward` destabilised the N = 64 stationary droplet (peak 48 m/s at 0.048 s with a frozen mass flux, `1e102` at 0.051 s with an active one) and was stable only at `dt <= dt_sigma/8`, because the L-stable Euler damped the grid-scale capillary mode every step ([negative deck 1/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/1/1)). The mode analysis explains why: at N = 64 the mode's physical and numerical damping exceeds `c dt`, at N = 128 it does not ([PCS 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1578-L1580)). On the current pipeline BDF2 against Euler on matched windows moves the gain by +11.1 / +2.9 / -3.0 percent, sign-flipping, with volume and shape within 1.2 percent, so BDF2 is kept on formal grounds at no measurable cost ([PSH 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L362-L384)); at 3D R/h = 15.8 `backward` no longer diverges and improves the volume error by 20 percent ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L317-L327)).

## Where in the code

- The switch `levelSet.capillaryForceCentring endStep | midpoint`, fatal on any other word ([`createSLFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H#L118-L155)); the midpoint drift of half a step ([`slAlphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L122-L125)) and the second pipeline in the main loop ([`leiaSemiLagrangianLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/leiaSemiLagrangianLevelSetTwoPhaseFoam.C#L418)).
- The token `CAPILLARY_FORCE_CENTRING endStep` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L328-L336)), set in the gate's `twoPhaseCoupling` block ([`config/gates/methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L77)).
- The force is read once per outer corrector at [`UEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H#L8), after the interface pipeline of `slAlphaEqn.H`.

## Evidence

| claim | number | where |
|---|---|---|
| the amplification matrix | (a) 1.118 to 3.297; (b) = (c) 1.000000 for omega dt of 0.5 to 2.0, 7.743 at 3.142 | [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L550-L564), DERIVED |
| the midpoint centring on the translating droplet | DIVERGED at step 5467, 0.0593 s, against 0.0868 s for endStep; serial reference 0.0775 s | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3738), MEASURED, one resolution |
| endStep had no coupled measurement before that | the linear analysis and a 20-step smoke | [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3589), MEASURED |
| the retracted anti-damping numbers | +88 1/s at R/h = 10 to +22 1/s at R/h = 25, and the prediction table of `config/viscousHorizon2D.yaml` | [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L486-L496), [`config/viscousHorizon2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/viscousHorizon2D.yaml#L19-L36), RETRACTED |
| the m = 8 wave the centring was meant to fix | period 1.56 ms, h-independent; envelope -97 to +212 1/s against viscous -128 | [MC deck 9/1](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/9/1), MEASURED |
| backward destabilised the Euler-era N = 64 case | peak 48 at 0.048 s (frozen rhoPhi), 1e102 at 0.051 s; stable at dt below dt_sigma/8; CN 0.9 diverges at 0.017 s | [negative deck 1/1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/1/1), MEASURED |
| BDF2 against Euler on the current pipeline | gain +11.1 / +2.9 / -3.0 percent at N = 64 / 128 / 256, matched step counts 1179 / 3333 / 9428; u0 identical to 4 digits | [PSH 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L362-L376), MEASURED |
| the outer-loop iteration of the force is not the lever | psiOuterCorrectors changes the rate by +2 percent at N = 128; the 3D gain by +0.4 / -0.1 percent | [PCS 16.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1352-L1363), [PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L209-L230), MEASURED, see [[concepts/psi-outer-correctors]] |
| the explicit step wall | omega_grid dt of 1.0 to 1.3 in practice, 2 in the linear analysis, pi in Popinet's limit | [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L581-L583), MEASURED |

## Why it failed, or why we think so

The centring hypothesis assumed the force at `t^n`. Reading the solver settled that it is at `t^{n+1}`, which is already neutrally stable, so the missing damping could not be a time-level effect. The same reading removed the explanation that the outer-corrector null result was "damping instead of anti-damping of the same size"; the null result stands, its explanation does not ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L493-L496), [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530-L536)).

## Decisions

- `CAPILLARY_FORCE_CENTRING endStep`, `MOMENTUM_DDT_SCHEME backward` ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L371), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L398)); see [[decisions/momentum-schemes-bdf2-upwind]] and [[concepts/ddt-scheme-pairing-bdf2]].
- `midpoint` stays selectable for research; it costs a second interface pipeline per step.

## Open questions

1. The midpoint centring on the stationary and oscillating droplet with the whole metric vector: not measured; the translating indicator is one resolution ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3729-L3730)).
2. The semi-Lagrangian foot uses an AB2-extrapolated velocity and the force uses `psi^{n+1}`; the joint order of the coupling is not stated in the record ([[concepts/departure-foot-ab2-centring]]).

## Related

[[hubs/surface-tension]], [[concepts/parasitic-current-mechanism]], [[concepts/capillary-time-step]], [[concepts/variational-capillary-force]], [[concepts/psi-outer-correctors]], [[concepts/departure-foot-ab2-centring]], [[concepts/ddt-scheme-pairing-bdf2]], [[models/semi-implicit-capillary-force]], [[decisions/momentum-schemes-bdf2-upwind]], [[retractions/force-at-n-not-n-plus-1]], [[cases/stationary-droplet]], [[cases/translating-droplet]], [[studies/method-comparison]].

## Log

### 2026-09-28
Created from the Shannon plan section 0e and its retraction, the method-comparison deck's time-centring track, the token comment and STATUS 11.13 and 11.14.
