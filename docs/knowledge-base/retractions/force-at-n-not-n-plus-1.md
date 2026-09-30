---
title: "The force built from psi^n anti-damps the capillary wave (retracted: the force is at n+1)"
description: "RETRACTED 2026-08-25 - the reading that an explicit capillary force at t^n anti-damps the m = 8 capillary wave by +omega^2 dt/2 (+88 to +22 1/s), which Popinet's n+1/2 staggering would remove: the solver builds the force from psi^{n+1}, the scheme is symplectic Euler with |lambda| = 1 for omega dt < 2, and the n+1/2 centring is spectrally identical to it"
aliases: [force time level retraction, symplectic Euler capillary force, capillaryForceCentring]
kind: retraction
status: retracted
part: surface-tension
tags: [retraction, part/surface-tension]
date: 2026-09-28
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/UEqn.H, config/viscousHorizon2D.yaml, cases/default.parameter]
sources: [PSH 0e, PSH 0e RETRACTION, DP CAPILLARY_FORCE_CENTRING, METHOD 8.1 row CAPILLARY_FORCE_CENTRING, CLAUDE BDF2 section, STATUS 11.14, MC deck Time centring, commit 3e231c8]
---
# The force built from psi^n anti-damps the capillary wave (retracted: the force is at n+1)

> RETRACTED 2026-08-25 (commit [3e231c8](https://github.com/leia-openfoam/leia/commit/3e231c8)). The claim, written on 2026-08-23, was "a capillary force evaluated explicitly at t^n gives the amplification |1 + i omega dt|, an anti-damping of +omega^2 dt/2 per unit time: +88 1/s at R/h = 10 falling to +22 1/s at R/h = 25, the right sign and order for the measured envelope, and exactly what Popinet's n+1/2 staggering removes; iterating the force towards n+1 with `psiOuterCorrectors` swaps anti-damping for damping of the same size, so neutrality needs the force at n+1/2" ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L485-L496)). The measurement is a reading of the solver: `slAlphaEqn.H` advects psi in place, psi^n to psi^{n+1}, and the narrow band, the phase indicator, alpha and the whole curvature chain are rebuilt from the advected field before `UEqn.H` calls the face force, so the force is built from psi^{n+1} ([plan-shannon 0e, RETRACTION](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L538-L542)). For the linearised capillary oscillator that is symplectic Euler, whose amplification matrix has det(M) = 1 exactly, so |lambda| = 1 for omega dt < 2. Computed per step at omega dt = 0.5 / 1.047 / 1.571 / 2.0: the force at t^n amplifies by 1.118 / 1.448 / 1.862 / 2.236, the force at t^{n+1} (ours) and at t^{n+1/2} (Popinet) both by 1.000000 ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L544-L556)). The two centrings are spectrally identical: det 1 and trace 2 - omega^2 dt^2, the same eigenvalues, amplification and phase error ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L558-L564)). Scope: the mechanism, the numbers derived from it (+22 to +88 1/s), the `psiOuterCorrectors` explanation and the prediction table of `config/viscousHorizon2D.yaml`. No data is void.

## The claim, and where it lived

- [plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L485-L501), `docs/plan-shannon-parasitic-currents.md`: "One contribution is quantifiable immediately" and the falling-under-refinement corollary; [item 3 of "what we are missing"](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L516-L517) proposed the n+1/2 stagger as "the one change with a derived reason to expect neutrality".
- [`config/viscousHorizon2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/viscousHorizon2D.yaml#L21-L41), the header: the old prediction table (anti-damping +135.2 / +94.9 / +66.8 1/s, net +7 / -33 / -61 1/s at R/h = 25.0 / 31.7 / 40.0) now sits under a NOTE that names the retraction.
- The retraction itself: [plan-shannon 0e, RETRACTION](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530-L583).

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The solver evaluates the force at t^n. | `slAlphaEqn.H` advects psi in place; band, indicator, alpha and curvature are rebuilt from psi^{n+1} before `UEqn.H:8` calls `fSigma->faceSurfaceTensionForceFlux()`. | MEASURED (read from the source), [plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L538-L542) |
| An explicit force at t^n amplifies the oscillator. | True, but not ours: at omega dt = 0.5 / 1.047 / 1.571 / 2.0 the t^n force amplifies 1.118 / 1.448 / 1.862 / 2.236 per step; the t^{n+1} force gives 1.000000, and so does t^{n+1/2}. At omega dt = 3.142 all three amplify (3.297, 7.743, 7.743). | DERIVED, [plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L544-L556) |
| Centring at n+1/2 changes the linear stability. | (b) and (c) have det 1 and trace 2 - omega^2 dt^2: the same eigenvalues, amplification and phase error per step. The centring changes the offset between the velocity and interface samples, the accuracy of the coupling, not the step limit. | DERIVED, [plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L558-L564), [DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L329-L336) |
| `psiOuterCorrectors` swaps anti-damping for damping. | Its gate moved the endpoint by +0.4 % and -0.1 % at R/h = 10.0 and 15.8: the lag is exonerated, and the retracted mechanism was the wrong explanation of that null result. | MEASURED, [plan-shannon 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L207-L230) |
| The midpoint centring, once built, helps the coupled case. | On the translating droplet (N = 100, np 4, fixed binaries) `capillaryForceCentring midpoint` diverged at step 5467 (t = 0.0593 s) against step 8000 (0.0868 s) for `endStep`: 32 % earlier, outside the decomposition scatter. One resolution. | MEASURED, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3738) |

## What survives

1. The measurements the mechanism was meant to explain ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L566-L569)): the oscillation is the m = 8 capillary wave (period and azimuthal spectrum agree), its physical envelope is -128 1/s and the measured envelope reaches +212 1/s, the mode is mesh-locked, 83 % of the curvature error is non-absorbable variation, and the horizon is 0.6 % of T_v, see [[concepts/parasitic-current-mechanism]].
2. The replacement mechanism: a symplectic integrator conserves only for a force derived from a potential; the delivered curvature error depends on where the interface sits relative to the mesh, so the force does net work per cycle ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L571-L579)). That is the quantity a fix has to remove, and it is the argument for [[concepts/variational-capillary-force]].
3. The step limit: the practical threshold omega_grid dt of about 1.0 to 1.3 against the linear bound omega dt < 2 ([plan-shannon 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L581-L583)), see [[concepts/capillary-time-step]].
4. The `endStep` default with `midpoint` selectable: the two are spectrally identical for the linear oscillator, and `midpoint` costs a second interface pipeline per step ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L398), [DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L329-L336)), see [[concepts/force-time-centring]].
5. The placement rule of [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L248-L261): BDF2 evaluates fluxes and sources at t^{n+1}, matching the force built from psi^{n+1}, see [[decisions/momentum-schemes-bdf2-upwind]].

## Propagation (checklist, same commit)

Done:

- [x] The plan document: the [RETRACTION subsection](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530-L583) of 0e.
- [x] `config/viscousHorizon2D.yaml`: the [NOTE](https://github.com/leia-openfoam/leia/blob/8867581/config/viscousHorizon2D.yaml#L21-L26) above the old prediction table.
- [x] `cases/default.parameter`: the [`CAPILLARY_FORCE_CENTRING` token](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L329-L336) (endStep default, midpoint added as Stormer-Verlet, "SPECTRALLY IDENTICAL").
- [x] `METHOD.md`: the [8.1 row](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L398); `CLAUDE.md`: the [BDF2 section](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L248-L261).
- [x] The decks: the method-comparison deck, track "Time centring", [slide 9](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/9) with the amplification matrix and [slide 9/9](https://leia-openfoam.github.io/leia/decks/level-set-method-comparison.html#/9/9) ("centring at n+1/2 restores neutrality" listed as a falsified reading); the quadratic SL negative-results deck, track "Time discretization", [slide 1](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/1).
- [x] `STATUS.md`: the [four discriminators of 2026-09-27](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3738) and the [best-settings row](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3589) ("the linear analysis; a 20-step smoke; no coupled measurement"); the SL article's late-instability paragraph ([d1e3414](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2295-L2299), "the midpoint force centring 32 % earlier").
- [x] The line in [[retraction-log]].

Still missing:

- [ ] `STATUS.md` has no line of its own for this retraction; it lives in the plan and in the config header.
- [ ] The work per cycle of the position-dependent curvature error was never measured directly; the proposed columns `kappaStdDevBand` and `kappaStdDevActiveFaces` ([config header](https://github.com/leia-openfoam/leia/blob/8867581/config/viscousHorizon2D.yaml#L43-L49)) have no recorded result.
- [ ] `endStep` against `midpoint` has one coupled measurement at one resolution (the N = 100 discriminator); no ladder.

## Related

Hubs: [[hubs/surface-tension]]. Siblings: [[concepts/force-time-centring]], [[concepts/parasitic-current-mechanism]], [[concepts/psi-outer-correctors]], [[concepts/variational-capillary-force]], [[concepts/balanced-force-csf-flux]], [[concepts/capillary-time-step]], [[decisions/momentum-schemes-bdf2-upwind]], [[studies/shannon-parasitic-currents-campaign]], [[studies/method-comparison]], [[cases/stationary-droplet]].

## Log

### 2026-09-28
Written from plan-shannon 0e and its RETRACTION, cases/default.parameter and STATUS 11.14. Retracted 2026-08-25. Entered in [[retraction-log#2026-08]].
