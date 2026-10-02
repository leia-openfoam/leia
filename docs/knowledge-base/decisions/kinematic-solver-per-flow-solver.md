---
title: "One kinematic solver per flow solver, with the same interface step"
description: "Settled 2026-10-01 by the author: every flow solver keeps one kinematic solver that runs its interface step with a prescribed velocity, so leiaSemiLagrangeLevelSetFoam stays for leiaSemiLagrangianLevelSetTwoPhaseFoam and leiaLevelSetFoam serves leiaLevelSetTwoPhaseFoam; the stale narrow band showed what happens when a pair drifts apart."
aliases: [kinematic counterpart, kinematic solver rule, one kinematic solver per flow solver]
kind: decision
status: settled
part: verification
tags: [decision, part/verification]
date: 2026-10-01
date_settled: 2026-10-01
decided_by: ["author decision 2026-10-01"]
code: [applications/solvers/leiaSemiLagrangeLevelSetFoam/leiaSemiLagrangeLevelSetFoam.C, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/solvers/leiaLevelSetFoam/leiaLevelSetFoam.C, applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H]
sources: ["STATUS 11.23 (L4947-L5055)", "STATUS 11.22 item 5", "commit 6f63418a"]
---
# One kinematic solver per flow solver, with the same interface step

> Settled 2026-10-01 by the author: "we need one kinematic solver for each flow solver if we use reversed prescribed analytic velocity" ([STATUS 11.23](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L4947-L4957)). A kinematic solver is the interface step of a flow solver, run with a prescribed velocity instead of the momentum solution. A reversed-flow test verifies the flow solver only if the two run the same interface code. So `leiaSemiLagrangeLevelSetFoam` stays as the kinematic solver of `leiaSemiLagrangianLevelSetTwoPhaseFoam`, and `leiaLevelSetFoam` (`eulerian`) is the kinematic solver of `leiaLevelSetTwoPhaseFoam`. My proposal of 2026-09-30 to retire the SL kinematic solver as well is withdrawn: the `semiLagrangian` model of `leiaLevelSetFoam` has no `projectedFlux` trace, so it does not run the interface step of the SL two-phase solver ([[models/level-set-advection]]). The retirement of `leiaRedistancedLevelSetFoam` agrees with the rule: it was a second Eulerian kinematic solver next to `leiaLevelSetFoam` ([[decisions/retire-redistanced-solver]]).

## The question

After the redistanced solver was retired, I asked whether `leiaSemiLagrangeLevelSetFoam` could go too, because `leiaLevelSetFoam` has a `semiLagrangian` advection model. The pair test of 2026-09-30 answered the narrow question (not equivalent), and the author answered the general one: the number of kinematic solvers follows the number of flow solvers, not the number of advection models.

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| the kinematic SL solver against its flow solver, before 2026-10-01 | order of the band and the phase indicator | alpha with the band of psi^n; the two-phase SL solver refreshes the band first | MEASURED (code), [slAlphaEqn.H L167-L172](https://github.com/leia-openfoam/leia/blob/d2984c5e/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L167-L172), [STATUS 11.23 item 1](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L4958-L4964) |
| the drift, 2D vortex N = 32 | relative volume error at step 53 | 6.5e-04 in the SL kinematic solver against 9.2e-05 with the current band (7.1x) | MEASURED, [STATUS 11.22 item 5](https://github.com/leia-openfoam/leia/blob/fd3e6e8b/STATUS.md#L4923-L4938) |
| the drift in a published table, 2D vortex, CFL 1.0 | fitted shape and volume order | 2.374 / 3.540 with the stale band, 2.465 / 3.193 with the fix | MEASURED, [STATUS 11.23 item 7](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5009-L5036), [[retractions/sl-stale-band-alpha-metrics]] |
| the Eulerian pair | order of the band and the phase indicator | the same in both | MEASURED (code), [alphaEqn.H L52-L92](https://github.com/leia-openfoam/leia/blob/d2984c5e/applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H#L52-L92), [leiaLevelSetFoam.C L174-L181](https://github.com/leia-openfoam/leia/blob/d2984c5e/applications/solvers/leiaLevelSetFoam/leiaLevelSetFoam.C#L174-L181) |

## What it does not cover

1. The two solvers of a pair share the interface step only by convention: they are separate code. The SL pair has two time loops around `slAdvection`; the Eulerian pair has two implementations of the psi equation (`eulerianAdvection.C` and `alphaEqn.H`). A shared library function for the interface step would make the match a property of the code. Not done; it is a refactor with a bit-identity gate ([STATUS 11.23 item 11](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L5051-L5055)).
2. Only the order of the band and the phase indicator was compared in both pairs. The rest of the interface step of the Eulerian pair was not compared end to end.

## Related

[[decisions/retire-redistanced-solver]] - [[models/level-set-advection]] - [[retractions/sl-stale-band-alpha-metrics]] - [[cases/kinematic-advection-cases]] - [[concepts/advection-regression-set]] - [[hubs/verification]] - [[decision-log]]

## Log

### 2026-10-01
SETTLED by the author; the stale band of the SL kinematic solver is fixed in 6f63418a. Entered in [[decision-log#2026-10]].
