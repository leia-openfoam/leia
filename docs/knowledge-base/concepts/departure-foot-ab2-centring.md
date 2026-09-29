---
title: "The departure foot: AB2, departure-centred"
description: "The coupled solver traces the semi-Lagrangian foot with an Adams-Bashforth estimate of the half-step velocity, expanded about the departure time; the arrival form fed with the same two levels leaves a +dt^2 d_t u error; the taylor and rk2 integrators give the same divergence time (2026-09-27)."
aliases: [departure foot, AB2 foot, SL_FOOT_INTEGRATOR, departure point kernel]
kind: concept
status: settled
part: advection
tags: [concept, part/advection]
date: 2026-09-28
date_settled: 2026-08-14
decided_by: [applications/test/leiaTestDeparturePoint, config/psiOuterCorrectorsGain3D.yaml, "author decision 2026-08-28, cases/default.parameter L403-L427"]
code: [src/leiaLevelSet/semiLagrangian/pointValueScheme.C, src/leiaLevelSet/semiLagrangian/pointValueScheme.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/test/leiaTestDeparturePoint/leiaTestDeparturePoint.C]
sources: ["METHOD 2.1 (L69-L96)", "METHOD 8.1 row SL_FOOT_INTEGRATOR (L376)", "PCS 16.2 (L1352-L1372)", "STATUS 4 time axis (L461-L468)", "STATUS 4 gradU invalidation (L619-L640)", "STATUS 11.15 discriminators (L3725-L3738)", "PSH gate result (L207-L230)", "SL article sec:foot (L231-L266)", "DP L403-L427 and L1028"]
---
# The departure foot: AB2, departure-centred

> Verdict (2026-09-28). The semi-Lagrangian update needs the foot `x_d` of the characteristic that arrives at the cell centre. The kernel expands the trace to second order with a material acceleration term ([SL article `sec:foot` L231-L263](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L231-L263)). The coupled solver advances the interface before its momentum solve, so it holds only `u^n` and `u^{n-1}`. It therefore hands the kernel the Adams-Bashforth estimate `u* = u^n + (dt/dt0)(u^n - u^{n-1})`, which makes the trace departure-centred and second order from those two levels ([METHOD 2.1 L69-L90](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L69-L90), [slAlphaEqn.H L116-L126](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L116-L126)). The arrival form fed with the same two levels leaves a `+dt^2 d_t u` foot error: 2 to 4 % of the per-step displacement early in an oscillating droplet, 35 to 47 % by t = 0.02 ([METHOD L92-L96](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L92-L96)). The shipped kernel has per-step order 3.00 on an exact rotating characteristic ([PCS 16.2 L1364-L1372](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1364-L1372)). The integrator `taylor` is the default without a gate ([METHOD 8.1 L376](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L376)); `rk2` moved the translating divergence time from 0.0868 s to 0.0842 s, inside the decomposition scatter ([STATUS 11.15 L3732-L3738](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3732-L3738)).

## What it is

The update is `psi^{n+1}(x_c) = psi^n(x_d)` ([METHOD 2 L63-L67](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L63-L67)). The kernel locates the foot by a Taylor expansion about the arrival time:

    x_d = x_c - u_new dt + (dt^2/2) [ (u_new - u_old)/dt + (u_new . grad) u_new ]

([SL article eq:foot L254-L259](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L254-L259), [pointValueScheme.C L363-L369](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L363-L369)). Both approximations perturb an `O(dt^2)` term by `O(dt)`, so the foot is second order per step ([SL article L260-L263](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L260-L263)). The linear displacement of the kernel is the average of the two velocity slots, so the kernel wants `(u^{n+1}, u^n)` and forms `u^{n+1/2}` itself ([DP L413-L420](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L413-L420)).

The kinematic solvers supply the exact `u(t^{n+1})` and `u(t^n)`. The coupled solver cannot: the interface moves before the momentum solve. It fills the new slot with the Adams-Bashforth extrapolation `u^n + (dt/dt0)(u^n - u^{n-1})` and the old slot with `u^n`. The kernel then reduces algebraically to the departure-centred form

    x_d = x_c - (dt/2) (3 u^n - u^{n-1}) + (dt^2/2) (u^n . grad) u^n

([METHOD 2.1 L76-L90](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L76-L90)). Relative to the arrival form only the temporal-derivative term changes sign; the convective term keeps its sign under both centrings ([METHOD L82-L84](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L82-L84)). The ratio `dt/dt0` keeps the difference quotient correct under a varying step ([slAlphaEqn.H L96-L102](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L96-L102)).

Two integrators exist for the `dt^2` term ([pointValueScheme.H L68-L80](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L68-L80)):

| integrator | the `dt^2/2` term comes from | reads `grad(U)` | code |
|---|---|---|---|
| `taylor` (default) | the Gauss-linear cell gradient of the new-slot velocity | yes | [pointValueScheme.C L415-L432](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L415-L432) |
| `rk2` | a second velocity sample at the midpoint `x_m = x_c - (dt/2) u_h(x_c)`, `u_h = (u_new + u_old)/2` | no | [pointValueScheme.C L371-L414](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L371-L414) |

Both are second order in `dt`. They differ only in how the `dt^2` correction is obtained ([pointValueScheme.C L379-L384](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L379-L384)).

## Why it matters

The foot error is the temporal error of the transport. An arrival-centred trace fed with `(u^n, u^{n-1})` is first order in time ([METHOD L92-L94](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L92-L94)). This error is invisible in every steady or uniform test, because `d_t u = 0` makes the two forms identical ([METHOD L94-L95](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L94-L95)). The `taylor` integrator also reads `fvc::grad(U)`, the operator that was biased by `O(1)` in every processor-adjacent cell until 2026-08-26 ([STATUS L621-L625](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L621-L625), [[retractions/gradu-coupled-patch-contamination]]). The foot must stay inside the reconstruction stencil, so the transport is limited to `CFL <= 1`; the kernel warns when the guard trips ([SL article L263-L266](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L263-L266)).

## Where in the code

- The kernel: `src/leiaLevelSet/semiLagrangian/pointValueScheme.C`, feet computed once per advect call ([L363-L432](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L363-L432)); the integrator word `footIntegrator` ([pointValueScheme.C L63](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.C#L63)); the token `SL_FOOT_INTEGRATOR taylor` ([DP L1028](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1028)).
- The coupled call site: `slAlphaEqn.H`, snapshot of `psi^n` ([L90-L92](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L90-L92)), the AB2 estimate `UextStar` ([L116-L126](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L116-L126)), the old slot written on pass 1 only ([L129-L139](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L129-L139)).
- The unit gate: `applications/test/leiaTestDeparturePoint` drives the shipped `slAdvection::advect` with `u = cos(omega t) Omega e_z ^ (x - c)`, whose backward characteristic is exact; the feet are read off advected linear fields ([leiaTestDeparturePoint.C L5-L44](https://github.com/leia-openfoam/leia/blob/8867581/applications/test/leiaTestDeparturePoint/leiaTestDeparturePoint.C#L5-L44)). Curated output: [departure_point_errors.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/departure_point_errors.csv).

## Evidence

| claim | number | where |
|---|---|---|
| The arrival form with `(u^n, u^{n-1})` leaves `+dt^2 d_t u`; share of the per-step displacement in an oscillating droplet | 2 to 4 % early, 35 to 47 % by t = 0.02 (config not named) | MEASURED, [METHOD 2.1 L92-L96](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L92-L96) |
| Per-step foot order of the shipped kernel, exact rotating characteristic, omega = 0 | 3.00 | MEASURED, [PCS 16.2 L1364-L1367](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1364-L1367) |
| Relative foot error with the exact velocity supply / with the AB2 substitution | 0.07 (omega dt)^2 / 0.35 (omega dt)^2 | MEASURED, [PCS 16.2 L1367-L1369](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1367-L1369) |
| AB2 mislocation of the stiffest capillary mode's foot at omega_max dt = 0.52 | 9.5 % of its per-step displacement; the growth rate r is unchanged when passes >= 2 use the true iterate | MEASURED, [PCS 16.2 L1369-L1372](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1369-L1372), [STATUS L464-L466](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L464-L466) |
| Re-advecting with the true iterate on passes >= 2 (which also removes the AB2 extrapolation), 3D stationary droplet | max abs U +0.4 % at N_L = 60, -0.1 % at N_L = 95 | MEASURED, [PSH L209-L230](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L209-L230) |
| `rk2` against `taylor`, translating droplet, N = 100, np 4 | diverged at t = 0.0842 s against 0.0868 s; the scatter of the unchanged case spans 0.0775 to 0.0904 s | MEASURED, one resolution, [STATUS 11.15 L3725-L3738](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3738) |
| The `dt^2/2` term of `taylor` consumed a `grad(U)` biased by `O(1)` at processor patches | 31 parallel kinematic SL studies contaminated; 1D gate 4.3e-2 at np = 8 where serial was exact | MEASURED, [STATUS L619-L636](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L619-L636) |
| Foot-radius guard: the foot must stay inside the stencil | CFL <= 1, warn only | DERIVED, [SL article L263-L266](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L263-L266) |

## Decisions

- The coupled call site feeds the AB2 estimate on pass 1 and the true velocity iterate on later passes; the old slot is never refreshed on passes >= 2 ([DP L403-L427](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L403-L427), [[concepts/psi-outer-correctors]]).
- `SL_FOOT_INTEGRATOR taylor`, decided without a gate ([METHOD 8.1 L376](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L376)); the one-resolution discriminator of 2026-09-27 found no effect of `rk2` ([[models/sl-scheme]]).
- The x_d kernel is exonerated as a growth carrier of the parasitic current ([STATUS L461-L468](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L461-L468)).

## Open questions

1. The 2 to 4 % and 35 to 47 % figures of METHOD 2.1 name no config. The oscillating-droplet run that produced them is not recorded.
2. No ladder with `rk2` is on record; the integrator never reads `grad(U)`, which the seam class of defects makes attractive ([pointValueScheme.H L72-L79](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/pointValueScheme.H#L72-L79)).
3. The AB2 mislocation of 9.5 % at the capillary operating point is quantified but its effect on the period of the oscillating droplet is not measured ([[cases/oscillating-droplet]]).

## Related

[[hubs/advection]] - [[models/sl-scheme]] - [[models/sl-reconstruction]] - [[concepts/psi-outer-correctors]] - [[concepts/trace-velocity-projected-flux]] - [[concepts/normal-projected-sl]] - [[concepts/force-time-centring]] - [[concepts/seam-checks-and-decomposition-invariance]] - [[retractions/gradu-coupled-patch-contamination]] - [[cases/oscillating-droplet]] - [[studies/sl-quadratic-pre-print]]

## Log

### 2026-09-28
Created from METHOD 2.1, the SL article, the plan section 16.2, the coupled call site and the kernel.
