---
title: "PSI_FILTER: none, because a filter is a research instrument"
description: "PSI_FILTER none in the global default, decided by the filter-off scoring rule of 2026-08-25: the theta = 0.2 band filter is 5.86x better at R/h = 15.8 and 1.61x worse at R/h = 10.0, so its benefit is a measurement of the defect, not a fix; every candidate is scored with all filters off"
aliases: [PSI_FILTER none]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-08-25
decided_by: [config/filterOffAmplifier3D.yaml, "author decision 2026-08-25"]
code: [cases/default.parameter, config/gates/methodGate2D.yaml]
sources: ["METHOD 8.1 row PSI_FILTER (L385)", "CLAUDE no-filtering section (L426-L452)", "STATUS 4 (L892-L922, L1016-L1032)", "PSH L328-L352", "PCS L90-L97", "DP L485-L487"]
---
# PSI_FILTER: none, because a filter is a research instrument

> `PSI_FILTER none` in the global default ([DP L485-L487](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L485-L487)) and in the gate's coupling block ([methodGate2D.yaml L84](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L84)). Decided by the filter-off scoring rule ([METHOD 8.1 L385](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L385)), written into AGENTS.md on 2026-08-25 (69323bb) and CLAUDE.md on 2026-08-28 (ffb6bd7) ([CLAUDE L426-L452](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L426-L452)): filtering is a research instrument only; production must be stable without it; every candidate is scored with all filters off. The measurement behind the rule ([STATUS L1016-L1032](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1016-L1032)): at matched initial kick the theta = 0.2 biharmonic band filter is 5.86x better at R/h = 15.8 and 1.61x worse at R/h = 10.0, where it turns a damped state into growth. That behaviour is a tuning knob that needs a resolution-dependent coefficient, not a model. The filter is also not the source of the parasitic growth: with it removed the two fine 3D arms still grow ([STATUS L1004-L1007](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1004-L1007)).

## The question

The psi filter `biharmonicBand` applies `psi -= theta L(L(psi))` in the band plus one ring once per step; it is the one explicit dissipation of the scheme ([DP L485-L487](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L485-L487), [filterOffAmplifier3D.yaml L8-L13](https://github.com/leia-openfoam/leia/blob/8867581/config/filterOffAmplifier3D.yaml#L8-L13)). It carried the 3D ladder through the horizon at N = 64 to 256 with the band gradient pinned. Two questions decided its status: is it the source of the growth or its damper, and does its benefit hold across resolutions?

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| 3D stationary droplet, filter off, R/h = 15.8 and 20.0 | per-step gain `A` | +7.68e-04 and +1.16e-03: the fine arms still grow | MEASURED, [STATUS L1019-L1022](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1019-L1022) |
| filter off against theta = 0.2, R/h = 10.0, matched `u0` (2.120e-04 against 2.118e-04) | max abs U at T; e-fold count G | 1.7249e-04 against 2.7722e-04; G -0.21 to +0.27: 1.61x worse with the filter | MEASURED, [PSH L328-L340](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L328-L340) |
| the same at R/h = 15.8, matched `u0` (4.085e-05 against 4.061e-05) | the same | 1.3838e-03 against 2.3627e-04; G +3.52 to +1.76: 5.86x better with the filter | MEASURED, [PSH L328-L340](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L328-L340) |
| the 2D theta sweep | gain at theta 0.05 against 0.2 | 5.93e-4 against 2.16e-4, 1.50e-4 against 1.37e-4, 1.96e-4 against 7.61e-5: more damping, less gain | MEASURED, [filterOffAmplifier3D.yaml L15-L17](https://github.com/leia-openfoam/leia/blob/8867581/config/filterOffAmplifier3D.yaml#L15-L17) |
| N = 256, theta = 0.05, extended horizon, np 8, after the seam fix | blow-up time | 0.167 s (onset 0.113 s) against 0.035 s unfiltered: a delay device, not a closure | MEASURED, [PCS L90-L97](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L90-L97) |
| every filtered result before 2026-08-19 | decomposition dependence | the filter operator was uncoupled on processor patches and its cell set decomposition-dependent | VOID before f83a1ab, [STATUS L892-L922](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L922), [[retractions/psi-filter-seam-bug]] |

Pre-registered read-out ([filterOffAmplifier3D.yaml L19-L34](https://github.com/leia-openfoam/leia/blob/8867581/config/filterOffAmplifier3D.yaml#L19-L34)): a. `A(h) > 0` at every resolution means every stable arm is stable only because the filter outruns the amplifier; b. the filter's damping `D` is h-independent per step. Outcome a occurred; the damping changed sign, so the filter carries its own source ([STATUS L1023-L1032](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1023-L1032)).

## What it does not cover

1. The rule covers every smoothing: the psi filter, curvature relaxation (`curvatureExtension.relax`), and any smoothing of the level set, the curvature or the force ([CLAUDE L428-L431](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L428-L431)). A parameter-free discretisation choice is not a filter: the `projectedFlux` trace passes this test ([STATUS L1148-L1164](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1148-L1164), [[decisions/sl-trace-velocity-projected-flux]]).
2. Three measurements that other decisions cite ran with the filter at theta = 0.2: `upwindConvection2D`, `upwindConvection3D` and `ddtOrderGain3D` ([upwindConvection2D.yaml L56-L57](https://github.com/leia-openfoam/leia/blob/8867581/config/upwindConvection2D.yaml#L56-L57), [ddtOrderGain3D.yaml L88-L89](https://github.com/leia-openfoam/leia/blob/8867581/config/ddtOrderGain3D.yaml#L88-L89)). Their verdicts are inertness verdicts and stand as such; see [[decisions/momentum-schemes-bdf2-upwind]].
3. The filter remains an instrument: its benefit at R/h = 15.8 measures the size of the corrugation defect ([[concepts/curvature-corrugation-and-the-fit]]).
4. The gradient-control sources of 2026-09 are not filters of psi in this sense; they are gated as candidates on the method gate ([[hubs/gradient-control]]).

## Related

[[hubs/advection]] - [[hubs/surface-tension]] - [[concepts/curvature-corrugation-and-the-fit]] - [[concepts/parasitic-current-mechanism]] - [[concepts/trace-velocity-projected-flux]] - [[decisions/sl-trace-velocity-projected-flux]] - [[decisions/momentum-schemes-bdf2-upwind]] - [[retractions/psi-filter-seam-bug]] - [[studies/shannon-parasitic-currents-campaign]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED by the rule of 2026-08-25, measured on 2026-08-20. Entered in [[decision-log#2026-08]].
