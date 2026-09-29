---
title: "MASS_FLUX_ALPHAF_SOURCE: the donor-cell plane, by author instruction"
description: "MASS_FLUX_ALPHAF_SOURCE donorPlane in the global default, on the author's instruction of 2026-09-01 (only Gauss upwind worked in rhoLENT); no measurement after the closed-box fix separates it from averagedPlanes"
aliases: [MASS_FLUX_ALPHAF_SOURCE donorPlane]
kind: decision
status: settled
part: mass-flux
tags: [decision, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: ["author decision 2026-09-01"]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, cases/default.parameter]
sources: ["METHOD 8.1 row MASS_FLUX_ALPHAF_SOURCE (L392)", "STATUS 0 (L55-L69)", "STATUS 11.13 (L3597)", "DP L690-L703 and L821-L823"]
---
# MASS_FLUX_ALPHAF_SOURCE: the donor-cell plane, by author instruction

> `MASS_FLUX_ALPHAF_SOURCE donorPlane` in the global default ([DP L690-L703](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L690-L703)). Decided by the author's instruction of 2026-09-01: "only Gauss upwind worked in rhoLENT" ([METHOD 8.1 L392](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L392)). No measurement decided it ([STATUS L3597](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3597)). The value restores what the eight droplet cases had before the token existed: every one of them hardcoded `faceInterpolation upwind` ([DP L695-L699](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L695-L699)). The studies that compared the alternatives ran on the closed-box translating case before 440107f and are void ([STATUS L55-L69](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L55-L69), [DP L821-L823](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L821-L823)).

## The question

The face density of the mass flux is `rho_f = alpha_f rho1 + (1 - alpha_f) rho2`, and `alpha_f` is the area fraction of the face cut by a reconstructed interface plane ([[models/mass-flux]]). Which cell supplies the plane? `donorPlane` takes the plane of the upwind (donor) cell; `averagedPlanes` takes the mean of the owner and neighbour fractions; `donorPlaneAdvected` displaces the donor plane by the same trajectory as the level set ([createMassFluxFields.H L41-L44](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L41-L44), [L89-L92](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L89-L92)). The legacy words `upwind` and `central` map to the first two ([L73](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L73)). The donor plane makes the face fraction agree with the upwind momentum convection `div(rhoPhi, U)` ([DP L701-L702](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L701-L702), [[decisions/momentum-schemes-bdf2-upwind]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| the author's instruction, 2026-09-01 | — | "only Gauss upwind worked in rhoLENT" | author decision, [DP L690-L693](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L690-L693) |
| `alphaFTest_donorPlane`, `alphaFTest_averagedPlanes` (closed box) | — | VOID; the curated tables are in `VOID_closedBox_20260902/` | [STATUS L55-L62](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L55-L62) |
| `alphaFTimeLevelTranslating2D`, `massFluxComparison2D`, `translatingClearOutlet2D` (closed box, 2026-09-02 before 17:25) | — | VOID; the token block carries the marker | [STATUS L65-L69](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L65-L69), [DP L821-L823](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L821-L823) |

Pre-registered read-out: none. The instruction replaced a measurement.

The token comment states that the solver defaults rhoLENT to `central`. At 8867581 the code default is `donorPlane` ([createMassFluxFields.H L56-L62](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L56-L62)); the comment is stale ([[models/mass-flux]]).

## What it does not cover

1. `donorPlane` against `averagedPlanes` on the repaired translating case: not measured.
2. The caution came from a plain higher-order scheme; the limited schemes were re-tested on the closed box only ([momentumDivScheme2D.yaml L60-L63](https://github.com/leia-openfoam/leia/blob/8867581/config/momentumDivScheme2D.yaml#L60-L63)), and the repaired matrix found the convection scheme not to be the mechanism ([SL article L2196](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2196)).
3. The time level of `alpha_f` (`trapezoid`) and the mass-flux projection (`projectMassFlux`) stay at their inert defaults; their only measurements are void ([[concepts/alphaf-source-donor-plane]], [[concepts/mass-flux-projection]]).

## Related

[[hubs/mass-flux]] - [[models/mass-flux]] - [[concepts/alphaf-source-donor-plane]] - [[concepts/mass-flux-projection]] - [[concepts/rholent-mass-flux]] - [[decisions/mass-flux-rholent]] - [[decisions/mass-flux-bound-rho]] - [[decisions/momentum-schemes-bdf2-upwind]] - [[retractions/closed-box-translating-droplet]] - [[cases/translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED by author instruction 2026-09-01, unmeasured. Entered in [[decision-log#2026-09]].
