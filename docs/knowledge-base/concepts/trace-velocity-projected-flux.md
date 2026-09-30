---
title: "The trace velocity: projectedFlux and the reconstruct operator"
description: "Tracing the semi-Lagrangian foot with fvc::reconstruct(phi) instead of the cell-centred velocity turned the stationary-droplet growth of +118 1/s into a decay of -52 1/s at N = 128; the mechanism is the reconstruct operator (70 %), not solenoidality (30 %), and not the extension (0 %) (2026-08-31)."
aliases: [projectedFlux, SL_TRACE_VELOCITY, traceVelocity, reconstructedU, SL_TRACE_FLUX, trace velocity]
kind: concept
status: settled
part: advection
tags: [concept, part/advection]
date: 2026-09-28
date_settled: 2026-08-31
decided_by: [config/fullHorizonStability2D.yaml, config/traceAmplifierDt2D.yaml, config/projFluxStationary2D.yaml, config/projFluxOscillating2D.yaml, config/kinematicTranslation2D.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/solvers/leiaSemiLagrangeLevelSetFoam/errorCalculation.H]
sources: ["STATUS 4 full-horizon gate (L1050-L1101)", "STATUS 4 reconstruct operator (L1103-L1164)", "STATUS 4 ladders launched (L1166-L1195)", "STATUS 11.5 traceFlux (L3300-L3303)", "STATUS 11.16 VOID columns (L3716-L3723)", "METHOD 8.1 row SL_TRACE_VELOCITY (L375)", "METHOD 8.2 controls (L437-L455)", "PSH reconstruct warning (L791-L805)", "DP L1029-L1040", "STATUS 11.19 (L4262-L4267)"]
---
# The trace velocity: projectedFlux and the reconstruct operator

> Verdict (2026-09-28). The coupled solver traces the semi-Lagrangian foot with `fvc::reconstruct(phi)`, the cell field of the pressure-projected face flux, instead of the cell-centred velocity ([createTransportFields.H L29-L41](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L29-L41)). On the 2D stationary droplet at N = 128 over the full 0.1 s horizon, this turns the growth of the parasitic current (+118.3 1/s with `cellCentred`) into a decay (-52.0 1/s), with no semi-implicit force and no tuned constant ([STATUS L1070-L1087](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1087)). The dt sweep with the confound control `reconstructedU` shows the mechanism: the reconstruct operator carries 70 % of the win, the velocity extension 0 % (a null step proved by bit-identity), and solenoidality 30 % ([STATUS L1121-L1137](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1121-L1137)). `projectedFlux` is the default since c935883 ([METHOD 8.1 L375](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L375), [DP L1029-L1036](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1029-L1036)). Under a frozen uniform stream the two traces are identical to eight digits, because `reconstruct(flux(U))` returns `U` exactly ([METHOD 8.2 L447-L455](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L447-L455)). Since 2026-09-26 the entry `traceFlux physical | extension` decides which face flux the trace reconstructs ([STATUS L3300-L3303](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3300-L3303)).

## What it is

The dictionary entry `levelSet.semiLagrangian.traceVelocity` selects the cell field handed to the departure-foot kernel ([[concepts/departure-foot-ab2-centring]]). The coupled solver knows four words, which form a ladder in which every step changes one thing ([createTransportFields.H L56-L67](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L56-L67), [L179-L193](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L179-L193)):

| word | traced field | reconstruct | face field | role |
|---|---|---|---|---|
| `cellCentred` | `Uext` (the extension field; identity when no extension model is selected) | no | none | the historical trace, the control |
| `reconstructedU` | `reconstruct(linearInterpolate(Uext) & Sf)` | yes | plain interpolated flux, not divergence-free | the confound control: same extension as `cellCentred`, same operator as `projectedFlux` |
| `reconstructedMomentum` | `reconstruct(linearInterpolate(U) & Sf)` | yes | plain interpolated flux of the raw momentum velocity | isolates the extension; bit-identical to `reconstructedU` when the extension is `none` |
| `projectedFlux` (production) | `reconstruct(phi)`, or `reconstruct(phiExt)` with `traceFlux extension` | yes | the pressure-projected flux, discretely divergence-free | isolates solenoidality |

The kinematic solver knows `cellCentred` and `projectedFlux` only, and a velocity extension enters it only through `traceFlux extension` ([errorCalculation.H L16-L39](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangeLevelSetFoam/errorCalculation.H#L16-L39), [L41-L82](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangeLevelSetFoam/errorCalculation.H#L41-L82)). Both non-default sources share the two-level machinery `Utrace` and `UtraceOld`, whose old slot holds the `t^n` level for the whole step ([createTransportFields.H L195-L209](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L195-L209), [slAlphaEqn.H L104-L114](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L104-L114), [L134-L139](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L134-L139)). Tokens: `SL_TRACE_VELOCITY projectedFlux`, `SL_TRACE_FLUX physical` ([DP L1036-L1040](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1036-L1040)).

## Why it matters

The hypothesis under test was that the cell-centred velocity carries the part of the velocity the pressure projection failed to absorb, the parasitic residual itself, and that tracing the interface with it closes a feedback loop ([traceAmplifierDt2D.yaml L7-L14](https://github.com/leia-openfoam/leia/blob/8867581/config/traceAmplifierDt2D.yaml#L7-L14)). The measurement kept the effect and changed the reading: `reconstructedU`, smoothed but not divergence-free, already removes 70 % of the `cellCentred` amplifier ([STATUS L1121-L1124](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1121-L1124)). `fvc::reconstruct(linearInterpolate(U) & Sf)` carries no tunable coefficient. It restricts the traced velocity to the discrete face-flux space, the space in which the pressure projection enforces the divergence constraint, so it is a discretisation choice and not a filter ([STATUS L1139-L1147](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1139-L1147), [[decisions/psi-filter-none]]). The statement the data supports: advect the interface with a velocity that lives in the space the pressure projection controls ([STATUS L1146-L1147](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1146-L1147)).

The same operator has a known limit. SAAMPLE proves `fvc::reconstruct` second order only for `C^1` fields; it diverges at the interface for a field with a gradient jump, and the velocity update of the pressure equation uses the same operator ([PSH L791-L805](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L791-L805)). The amplifier's output stage is therefore formally divergent across the interface for the non-`C^1` velocity the capillary jump creates ([PSH L801-L805](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L801-L805), [[concepts/parasitic-current-mechanism]]).

## Where in the code

- Selection and validation: [createTransportFields.H L85-L121](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L85-L121); `traceFlux` and its two warnings (an extension model that does not enter the trace is flagged): [L122-L157](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L122-L157).
- The trace field: [createTransportFields.H L179-L193](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L179-L193); refreshed per step at [slAlphaEqn.H L107-L110](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L107-L110).
- Kinematic solver: [errorCalculation.H L16-L82](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangeLevelSetFoam/errorCalculation.H#L16-L82).
- The banner prints the trace word; five droplet templates once had no `traceVelocity` entry, so every arm ran at the solver default ([printMethodBanner.H L12](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H#L12), [STATUS L1188-L1195](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1188-L1195)).

## Evidence

| claim | number | where |
|---|---|---|
| Full-horizon gate, 2D stationary droplet, N = 128, 13 334 steps: final max abs U and fitted rate over the last 20 % | `cellCentred` 4.65e-03 m/s, +118.3 1/s; `projectedFlux` 3.34e-06 m/s, -52.0 1/s; with the semi-implicit increment 7.53e-06 (+10.5) and 4.55e-06 (-45.3) | MEASURED, [STATUS L1070-L1075](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1070-L1075) |
| Same gate, volume error, shape L2, band min abs grad psi | `cellCentred` 9.86e-05, 3.59e-06, 0.9170; `projectedFlux` 1.14e-06, 1.96e-07, 0.9980 (shape starts at 1.91e-07) | MEASURED, same table and [L1077-L1078](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1077-L1078) |
| dt sweep (12 arms, steps 13 334 to 106 668): amplifier excess `dr` of `cellCentred` against `projectedFlux` at dt = 7.50e-06 | 85.85 1/s; `reconstructedU` 26.05 1/s | MEASURED, [STATUS L1109-L1119](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1109-L1119) |
| Split of the `cellCentred` amplifier | reconstruct operator 70 %, extension 0 % (null step, byte-identical), solenoidality 30 % | MEASURED, [STATUS L1126-L1137](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1126-L1137) |
| `dr_cellCentred` does not scale as dt | ratios per halving 1.63 / 2.60 / 4.13 against 2.00; `G_cellCentred` changes sign between dt/2 and dt/4 | MEASURED, [STATUS L1156-L1160](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1156-L1160) |
| `projectedFlux`'s own rate across the sweep | -46.85 to -37.91 1/s; the traces have not converged (spread 85.8 to 8.75) | MEASURED, [STATUS L1161-L1164](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1161-L1164) |
| Stationary 2D ladder, unabsorbed capillary residual, N = 32 / 64 / 128 / 256 | `projectedFlux` 2.37e-5 / 3.15e-6 / 8.53e-7 / 6.43e-7; `cellCentred` blows up at the two finest rungs, 5.48e-4 and 2.79e-2 | MEASURED, recorded only in [DP L1029-L1032](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1029-L1032) |
| Oscillating 2D ladder, N = 32 / 64 / 128 | `projectedFlux` reaches the full horizon at all three rungs; `cellCentred` at none of the two finest | MEASURED, [DP L1032-L1033](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1032-L1033) |
| Uniform translation, kinematic, one way (re-measured 2026-09-29) | `projectedFlux` and `cellCentred` agree to 5.0e-12 column-scaled in every column of every row. The reversed run of 2026-09-01, which the token comment calls BIT-IDENTICAL, is VOID ([[retractions/reversed-2dtranslation]]) | MEASURED, [STATUS 11.19](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4266-L4267); VOID, [DP L1034-L1035](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1034-L1035), [[cases/kinematic-advection-cases]] |
| Frozen uniform stream, power iteration on the 95 969-cell pMesh | `rho(B)` 1.0110805 for both traces; `max abs Utrace` = `max abs U` = 0.069282032 | MEASURED, [METHOD 8.2 L439-L455](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L439-L455) |
| `VELOCITY_EXTENSION closestPoint` on the SL line with the default trace | CSVs identical to the baseline in every arm: the extension did not enter | MEASURED, [STATUS L3300-L3303](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3300-L3303) |
| `fvc::reconstruct` on a field with a gradient jump | formally divergent at the interface (SAAMPLE sec. 3.4) | DERIVED (literature), [PSH L791-L805](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L791-L805) |

## Decisions

- `SL_TRACE_VELOCITY projectedFlux`, not inert, default since c935883: [[decisions/sl-trace-velocity-projected-flux]]. The METHOD 8.1 row names `stationaryDropletFootEval*` as the deciding gate; the evidence on record is the STATUS entries of 2026-08-31 and the ladders in `cases/default.parameter` ([METHOD 8.1 L375](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L375), [DP L1029-L1035](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1029-L1035)).
- `SL_TRACE_FLUX physical` (inert, 2026-09-26); `extension` lets a velocity extension enter the trace ([DP L1037-L1040](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1037-L1040), [[models/velocity-extension]], [[concepts/halo-limited-extension]]).
- Every comparison carries its own `cellCentred` arm on the same commit: the historical t_blow table is retracted ([STATUS L1057-L1068](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1057-L1068), [[retractions/t-blow-baseline]]).

## Open questions

1. The refinement ladders `projFluxStationary2D` and `projFluxOscillating2D` have no STATUS entry of their own; their result lives in the token comment only ([DP L1029-L1035](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1029-L1035)). The oscillating read-out that was pre-registered, period and decay against Lamb, is not recorded ([projFluxOscillating2D.yaml L20-L28](https://github.com/leia-openfoam/leia/blob/8867581/config/projFluxOscillating2D.yaml#L20-L28)).
2. The 2026-08-31 studies ran on more than one rank before the coupled-face density fix of 2026-09-27, so they carry that seam error ([METHOD L304-L312](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L304-L312), [[concepts/coupled-face-density-defect]]). Which of them to re-run is open.
3. The 3D trace ladders `traceTranslating3Dhex` and `traceTranslating3Dpoly_r10p0/r12p7/r15p8` are VOID in their disturbance columns (no `dropletReferenceVelocity` in the 3D template) ([STATUS L3716-L3723](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3716-L3723)).
4. The reconstruct operator keeps only face-normal components; whether it damps the physical capillary mode is the concern the oscillating gate was built for and has no recorded verdict ([projFluxOscillating2D.yaml L11-L18](https://github.com/leia-openfoam/leia/blob/8867581/config/projFluxOscillating2D.yaml#L11-L18)).

## Related

[[hubs/advection]] - [[hubs/surface-tension]] - [[models/sl-scheme]] - [[models/velocity-extension]] - [[models/semi-implicit-capillary-force]] - [[concepts/departure-foot-ab2-centring]] - [[concepts/psi-outer-correctors]] - [[concepts/parasitic-current-mechanism]] - [[concepts/pressure-projection-and-linear-solvers]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/halo-limited-extension]] - [[concepts/coupled-face-density-defect]] - [[decisions/sl-trace-velocity-projected-flux]] - [[decisions/psi-filter-none]] - [[retractions/t-blow-baseline]] - [[cases/stationary-droplet]] - [[cases/oscillating-droplet]] - [[cases/kinematic-advection-cases]]

## Log

### 2026-09-28
Created from STATUS section 4 (2026-08-31), METHOD 8.1 and 8.2, the Shannon plan, the token comments and the solver headers.

### 2026-09-29
The uniform-translation row is re-measured on the one-way `2Dtranslation`: 5.0e-12 column-scaled, not bit-identical (STATUS 11.19). The reversed run of 2026-09-01 is void ([[retractions/reversed-2dtranslation]]). The `SL_TRACE_VELOCITY` comment in `cases/default.parameter` still says BIT-IDENTICAL at aaa0a7dd.
