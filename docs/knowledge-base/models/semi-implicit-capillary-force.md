---
title: "fv::option semiImplicitCapillaryForce"
description: "The Hysing/Raessi/SAAMPLE semi-implicit capillary term as an fvOption; off by default, not needed once projectedFlux traces the foot, and untested on the translating droplet after the deadlock fix (2026-09-28)."
aliases: [semi-implicit capillary force, SEMI_IMPLICIT_CAPILLARY]
kind: model
status: candidate
part: surface-tension
tags: [model, part/surface-tension]
date: 2026-09-28
date_settled:
decided_by:
code: [src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H, src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.C, cases/stationaryDroplet2D/system/fvOptions.template]
sources: [CLAUDE 4-rank section, STATUS full-horizon gate 2026-08-31, STATUS 11.13, PSH 0c and 0g, PCS 18, SL deck capillary coupling roadmap]
---
# fv::option semiImplicitCapillaryForce

> Verdict (2026-09-28). The term is off in every production study (`SEMI_IMPLICIT_CAPILLARY off`, [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L583-L593)). On the full-horizon stationary gate the arm that wins carries none of it: `off + projectedFlux` decays at -52.0 1/s, `increment + projectedFlux` at -45.3 1/s, and `increment + cellCentred` still grows at +10.5 1/s ([STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1050-L1101)). Its first cluster run deadlocked for 76 minutes, because the diagnostics called collective reductions inside a `Pstream::master()` guard ([CLAUDE, 4-rank gate](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L207-L214)); the code now reduces on every rank ([`semiImplicitCapillaryForce.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.C#L105-L125)). It has not been tried on the translating droplet after the fixes of 2026-09-27 ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3620-L3622)).

## What it is

The explicit CSF force uses the interface position of the current step. Linearising the force about the position the interface will have at `t^{n+1} = x^n + dt u^{n+1}` adds `dt sigma Lap_Gamma(u^{n+1}) delta_Sigma`, an interface-concentrated Laplace-Beltrami diffusion of the new velocity ([header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L32-L45)). The fvOption assembles that term with the face coefficient `mu_f = coeff sigma dt abs(snGrad(alpha))_f`, the same discrete object the balanced flux differences ([header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L71-L82)). Two forms exist:

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| inert | `form off` | settled | `addSup` returns without touching the matrix; verified bit-identical | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L583-L587) |
| SAAMPLE / Raessi value form | `form value` | candidate | `fvm::laplacian(mu_f, U) - mu_c normalCorrections(U)`; the term remains in the converged equation and dissipates; pairs with `psiOuterCorrectors no`; on the m = 2 mode `c` = +5.6e6 against +1.7e6 for the baseline, a suspected sign bug | [header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L49-L54), [SL deck 6/5](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/6/5) |
| deferred-Jacobian increment form | `form increment` | candidate | `fvm::laplacian(mu_f, U) - fvc::laplacian(mu_f, U)` acts on the outer-loop increment and vanishes at convergence; `c` = -6.2e5 on the m = 2 mode; the lagged term still stands at 0.07 percent of its magnitude after 3 outer correctors | [header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L56-L69), [STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1090-L1096) |

Sub-keys: `coeff` multiplies `sigma dt` (1 = backward Euler linearisation, 2/3 = the BDF2 Jacobian scale); `laplaceBeltrami` subtracts the projector terms in the value form; `kappaName` and `nHatName` name the registered fit curvature and fit normal; `diagnostics` writes a per-outer-call CSV ([header](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L89-L104)). The tokens are `SEMI_IMPLICIT_CAPILLARY`, `SEMI_IMPLICIT_CAPILLARY_COEFF` and `SEMI_IMPLICIT_LB` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L591-L593)).

## Why it matters

The stationary-droplet growth at N = 128 was measured as the explicit time coupling of the capillary force to the m = 2 interface mode: halving the step halves the rate, and the Richardson intercept is +0.03 1/s against an operating rate of 18.8 1/s ([PCS 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1554-L1592)). That made the semi-implicit force the principal licensed lever ([PCS 18](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1582-L1586)). SAAMPLE runs stably at `omega_grid dt = pi/2`, above the explicit wall of 1.0 to 1.3 measured here, through this door ([PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L773-L790)); the price SAAMPLE records is overdamping. Under the no-filtering rule the term is scored with the damping off for the stability claim and on only as the documented price of the larger step ([PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L784-L790), [[decisions/psi-filter-none]]).

## Where in the code

1. The class: [`semiImplicitCapillaryForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.H#L131-L175) and [`.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.C#L43-L44), registered in the fvOption table so any solver that folds in `fvOptions(rho, U)` can select it ([`UEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H#L45-L50)).
2. The library: `libleiaSurfaceTension` ([`Make/files`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/Make/files)).
3. The case: `system/fvOptions` rendered from `cases/stationaryDroplet2D/system/fvOptions.template`.
4. The collective calls: `gAverage`, `gMax` and `gSum` are reached on every rank; only the CSV write sits behind `Pstream::master()` ([`.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.C#L105-L125), [`.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/fvOptions/semiImplicitCapillaryForce.C#L190-L205)).

## Evidence

| claim | number | where |
|---|---|---|
| the full-horizon stationary gate, N = 128, 13334 steps, filters off | control (off + cellCentred) +118.3 1/s; off + projectedFlux -52.0; increment + cellCentred +10.5; increment + projectedFlux -45.3 | [STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1076), MEASURED |
| the winner carries none of the term | `off + projectedFlux` is best on every metric of the matrix | [STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1080-L1087), MEASURED; see [[concepts/trace-velocity-projected-flux]] |
| the increment form does not reach round-off at 3 outer correctors | lagged term contracts about 5x per pass, 0.07 percent of its magnitude on the last pass | [STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1090-L1096), MEASURED |
| the m = 2 mode slope on the stationary droplet (single geometry, N = 128) | baseline `c` +1.7e6, value form +5.6e6, increment form -6.2e5 | [SL deck 6/5](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/6/5), MEASURED, not a verdict |
| the class of fix was first dropped, then revived | psiOuterCorrectorsGain3D moved the answer by at most 0.4 percent, so the implicit-force class was dropped 2026-08-19; PCS 18.4 revived it 2026-08-17 (r_0 = 0) | [PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L207-L243), [PCS 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1582-L1586), MEASURED |
| the cluster deadlock | four arms alive but silent for 76 minutes inside the first momentum assembly; gated in serial only | [CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L207-L214), MEASURED 2026-08-28 |

## Why it failed, or why we think so

The instability the term was built to close (a `c dt` growth of the m = 2 mode) is removed without it by tracing the foot point with the projected face flux; the reconstruct operator accounts for 70 percent of the cellCentred amplifier and solenoidality for 30 percent ([STATUS 2026-08-31](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1103-L1164)). The term therefore became a candidate for a larger capillary step, not a stability fix. The serial-only gate was the process failure: an `fvOption` touches an `fvMatrix`, and the parallel smoke test is part of its inertness gate ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L244-L246), [[concepts/seam-checks-and-decomposition-invariance]]).

## Decisions

- `SEMI_IMPLICIT_CAPILLARY off` by default, bit-identical ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L583-L593)).
- The `increment` candidate pairs with `PSI_OUTER_CORRECTORS yes`, raised `N_OUTER_CORRECTORS` and residual-controlled outer iterations (`OUTER_U_TOL`, `OUTER_P_TOL`) ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L595-L600)).

## Open questions

1. The translating droplet after the fixes of 2026-09-27: not tried ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3620-L3622), [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3742-L3743)).
2. The value form's suspected sign bug ([SL deck 6/5](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/6/5)).
3. Whether the increment form reaches round-off with residual-controlled outer iterations, and what capillary step it then allows ([[concepts/capillary-time-step]]).

## Related

[[hubs/surface-tension]], [[models/surface-tension-force]], [[concepts/parasitic-current-mechanism]], [[concepts/force-time-centring]], [[concepts/capillary-time-step]], [[concepts/trace-velocity-projected-flux]], [[concepts/psi-outer-correctors]], [[concepts/seam-checks-and-decomposition-invariance]], [[decisions/psi-filter-none]], [[retractions/t-blow-baseline]], [[studies/shannon-parasitic-currents-campaign]].

## Log

### 2026-09-28
Created from the fvOption header, the full-horizon gate of 2026-08-31, the CLAUDE 4-rank section and the SL deck roadmap.
