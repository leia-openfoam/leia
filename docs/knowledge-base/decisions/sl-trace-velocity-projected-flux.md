---
title: "SL_TRACE_VELOCITY: trace the foot with the reconstructed projected flux"
description: "SL_TRACE_VELOCITY projectedFlux, decided by config/fullHorizonStability2D.yaml on 2026-08-31: growth rate -52.0 1/s against +118.3 1/s for cellCentred at N = 128 over the full horizon; the reconstruct operator carries 70 % of the win, solenoidality 30 %"
aliases: [SL_TRACE_VELOCITY projectedFlux]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-08-31
decided_by: [config/fullHorizonStability2D.yaml, config/traceAmplifierDt2D.yaml, config/projFluxStationary2D.yaml, config/projFluxOscillating2D.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H, cases/default.parameter]
sources: ["METHOD 8.1 row SL_TRACE_VELOCITY (L375, CORRECTED 2026-09-28)", "STATUS 4 (L1050-L1164)", "METHOD 8.2 (L437-L455)", "DP L1036-L1047", "PSH L791-L805"]
---
# SL_TRACE_VELOCITY: trace the foot with the reconstructed projected flux

> `SL_TRACE_VELOCITY projectedFlux` in the global default since c935883 (2026-09-01) ([DP L1036-L1043](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1036-L1043)), and in the `lineTokens` of the 2D method gate for every semi-Lagrangian arm ([methodGate2D.yaml L46-L49](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L46-L49)). Decided by `config/fullHorizonStability2D.yaml` on 2026-08-31 ([STATUS L1050-L1101](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1050-L1101)): on the 2D stationary droplet at N = 128, filters off, over the full horizon of 13 334 steps, the arm `off + projectedFlux` decays at -52.0 1/s where the control `off + cellCentred` grows at +118.3 1/s, and it is the best arm of the matrix on every metric. The mechanism is measured ([STATUS L1103-L1164](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1103-L1164)): the `fvc::reconstruct` operator removes 70 % of the cell-centred amplifier, the velocity extension 0 %, the solenoidality of the projected flux 30 %. METHOD 8.1 names `stationaryDropletFootEval*` as the gate; that is wrong, and the row carries the correction of 2026-09-28: the evidence is STATUS 2026-08-31 and the token comment ([METHOD 8.1 L375](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L375)).

## The question

The semi-Lagrangian foot needs a cell velocity. `cellCentred` traces the velocity-extension field `Uext`, which is the cell velocity itself when no extension is selected. `projectedFlux` traces `fvc::reconstruct(phi)`, the cell field rebuilt from the face flux that the pressure projection made divergence-free. `reconstructedU` is the confound control: the same reconstruct operator on the unprojected face flux of `Uext` ([createTransportFields.H L28-L62](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L28-L62)). The question of 2026-08-31 had two parts: does the trace act on the amplifier `G` of `max|U|(T) = u0(h) exp G(h)`, and is the mechanism the solenoidality or the smoothing of the reconstruct operator ([traceAmplifierDt2D.yaml L1-L60](https://github.com/leia-openfoam/leia/blob/8867581/config/traceAmplifierDt2D.yaml#L1-L60))?

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| `off + cellCentred` (control), N = 128, t = 0.1 s | final max abs U; fitted rate over the last fifth; volume error; shape L2; band min abs grad psi | 4.65e-03 m/s; +118.3 1/s; 9.86e-05; 3.59e-06; 0.9170 | MEASURED, [STATUS L1070-L1078](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1078) |
| `off + projectedFlux`, same run | the same vector | 3.34e-06 m/s; -52.0 1/s; 1.14e-06; 1.96e-07; 0.9980 | MEASURED, [STATUS L1070-L1078](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1078) |
| `increment + cellCentred`, `increment + projectedFlux` (semi-implicit force) | rate | +10.5 1/s; -45.3 1/s | MEASURED, [STATUS L1070-L1078](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1078) |
| dt sweep at fixed mesh, window 0.02 to 0.1 s | numerical amplifier `dr` at dt = 7.5e-06 s | cellCentred 85.85 1/s; reconstructedU 26.05 1/s; projectedFlux 0 (reference) | MEASURED, [STATUS L1103-L1126](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1103-L1126) |
| the split of the cell-centred amplifier | share | reconstruct operator 70 %; velocity extension 0 % (null step, byte-identical); solenoidality 30 % | MEASURED, [STATUS L1130-L1146](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1130-L1146) |
| 2D stationary ladder, N = 32 / 64 / 128 / 256 | unabsorbed capillary residual, relative L2 | projectedFlux 2.37e-5 / 3.15e-6 / 8.53e-7 / 6.43e-7; cellCentred blows up at the two finest rungs (5.48e-4, 2.79e-2) | MEASURED, [DP L1036-L1040](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1036-L1040) (the config `projFluxStationary2D` holds the pre-registration, not the result) |
| 2D oscillating ladder, three rungs | completion | projectedFlux reaches the full horizon at all three rungs, cellCentred at none of the two finest | MEASURED, [DP L1040-L1041](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1040-L1041) |
| kinematic uniform translation, 18 arms | difference between the traces | identical to every printed digit: `fvc::reconstruct` is exact for a uniform field | MEASURED, [kinematicTranslation2D.yaml L56-L59](https://github.com/leia-openfoam/leia/blob/8867581/config/kinematicTranslation2D.yaml#L56-L59) |
| frozen uniform stream, pMesh, power iteration | `rho(B)` with projectedFlux against cellCentred | 1.0110805 against 1.0110805, identical to 8 digits | MEASURED, [METHOD 8.2 L437-L455](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L437-L455) |

Pre-registered read-outs: [fullHorizonStability2D.yaml L38-L47](https://github.com/leia-openfoam/leia/blob/8867581/config/fullHorizonStability2D.yaml#L38-L47) (no stability prediction for projectedFlux; a diverged arm is a result), [traceAmplifierDt2D.yaml L43-L60](https://github.com/leia-openfoam/leia/blob/8867581/config/traceAmplifierDt2D.yaml#L43-L60) (`dr` halves with dt; `rho = dr_recU / dr_cc` decides the mechanism), [projFluxStationary2D.yaml L23-L37](https://github.com/leia-openfoam/leia/blob/8867581/config/projFluxStationary2D.yaml#L23-L37), [projFluxOscillating2D.yaml L20-L28](https://github.com/leia-openfoam/leia/blob/8867581/config/projFluxOscillating2D.yaml#L20-L28) (period and decay rate, not max abs U).

The choice is not a filter. `fvc::reconstruct(linearInterpolate(U) & Sf)` carries no coefficient; it restricts the traced velocity to the discrete face-flux space in which the pressure projection enforces the divergence constraint ([STATUS L1148-L1164](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1148-L1164), [[decisions/psi-filter-none]]).

## What it does not cover

1. A prescribed uniform velocity cannot separate the traces; the advantage lives in the reconstruct operator on a solved flux that carries a pressure correction ([METHOD 8.2 L446-L455](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L446-L455)). The kinematic gates do not test this decision.
2. The reconstruct operator is second order only for C^1 fields and diverges at a gradient jump (SAAMPLE, [PSH L791-L805](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L791-L805)). The trace reads the interface velocity through that operator.
3. With `projectedFlux` a velocity extension does not enter the trace: `VELOCITY_EXTENSION closestPoint` was a no-op on the semi-Lagrangian line ([CLAUDE L681-L685](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L681-L685)). `SL_TRACE_FLUX extension` (2026-09-26) is the path that admits one ([DP L1044-L1047](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1044-L1047)).
4. The 3D translating trace studies carry void disturbance columns: the 3D case had no `dropletReferenceVelocity` until 2026-09-27 ([STATUS L3733-L3739](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3733-L3739)).
5. `reconstructedU` and `reconstructedMomentum` are controls, not candidates.

## Related

[[hubs/advection]] - [[concepts/trace-velocity-projected-flux]] - [[concepts/parasitic-current-mechanism]] - [[concepts/balanced-force-csf-flux]] - [[models/semi-implicit-capillary-force]] - [[models/velocity-extension]] - [[decisions/psi-filter-none]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[retractions/t-blow-baseline]] - [[cases/stationary-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the measurement of 2026-08-31; the default changed on 2026-09-01. The gate that METHOD 8.1 names is wrong and is marked corrected there. Entered in [[decision-log#2026-08]].
