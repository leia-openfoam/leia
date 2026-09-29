---
title: "alpha_f source: donorPlane and the alternatives"
description: "The face area fraction of the mass flux comes from the reconstructed plane of the donor cell, by author instruction of 2026-09-01 and unmeasured; the averaged, advected and trapezoid alternatives were measured only on the closed box and are void"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: [author decision 2026-09-01]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H, cases/default.parameter]
sources: [METHOD 8.1 row MASS_FLUX_ALPHAF_SOURCE, DP 686-699, DP 817-837, STATUS 0, STATUS 11.2, G2 67]
---
# alpha_f source: donorPlane and the alternatives

> Verdict (2026-09-28). The face liquid-area fraction `alpha_f`, which builds the face density `rho_f = alpha_f rho1 + (1 - alpha_f) rho2`, is cut from the linear least-squares plane of the donor cell: the owner if `phi >= 0`, else the neighbour ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L33-L63)). `MASS_FLUX_ALPHAF_SOURCE donorPlane` is the default on the author's instruction of 2026-09-01, "only Gauss upwind worked in rhoLENT", with no measurement behind it ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L393), [[decisions/mass-flux-alphaf-donor-plane]]). The three alternatives were measured on `translatingDroplet2D` on 2026-09-02 before the inlet/outlet fix 440107f, so every number is void ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L55-L62), [[retractions/closed-box-translating-droplet]]): `averagedPlanes` diverged at step 4634 ([STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3260)), the `trapezoid` time level was "marginal" ([`alphaFTimeLevelTranslating2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/alphaFTimeLevelTranslating2D.yaml#L50-L85)), and the flux projection stopped the droplet ([[concepts/mass-flux-projection]]).

## What it is

`alphaFSource` is not an interpolation scheme. Each internal face has two candidate area fractions, one from the owner's plane and one from the neighbour's, and the switch selects between them ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L33-L47)):

| value | `alpha_f` | status |
|---|---|---|
| `donorPlane` | the upwind cell's plane cut against the face: the phase composition of the fluid that the flux carries | production, unmeasured |
| `averagedPlanes` | `1/2 (a_owner + a_neighbour)` | void (diverged at step 4634 on the closed box) |
| `donorPlaneAdvected` | `donorPlane` with the plane displaced by the same semi-Lagrangian step as the cell values: `psi^n(x - u_f dt) = n.x + d - (n.u_f) dt`, one dot product per face | void (closed box); its number is not recorded in the sources at 8867581 |

The old values `upwind` and `central` are still accepted with a warning; they were renamed on 2026-09-01 because `central` is not OpenFOAM's `linear` ([same](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L49-L54)). The plane is the same reconstruction that the Detrixhe-Aslam indicator uses, and the face is filled per triangle ([`faceAreaFraction.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H#L1-L18)). On a coupled face the two sides exchange their local fractions and apply the donor rule with the local flux sign ([[concepts/coupled-face-density-defect]]).

A second switch, `alphaFTimeLevel`, sets the time level of the geometric face density: `new` (the plane of the new interface) or `trapezoid` (`1/2 (rho_f^n + rho_f^{n+1})`), one extra plane reconstruction per outer corrector ([`updateFaceDensity.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H#L38-L64)). The production value is `new`.

## Why it matters

The donor choice is what the eight droplet cases hardcoded as `faceInterpolation upwind` before the token existed, and the momentum convection is `Gauss upwind` too, so `div(rhoPhi,U)` and the face fraction agree ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L691-L698)). The counting argument against the time-level blend is DERIVED: the mass balance is one constraint per cell, a two-level blend supplies one value per face from a pointwise rule, and no local formula reaches an exact balance; only a projection does ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L829-L834)).

## Evidence

| claim | number | where |
|---|---|---|
| the trapezoid time level does not reduce the residual | relative mass residual 0.17 to 0.29 against 0.17 to 0.22; travelled fraction 0.276 against 0.333 at N = 64, 0.734 against 0.725 at N = 128; volume error 0.0194 against 0.0397 at N = 64 | [`alphaFTimeLevelTranslating2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/alphaFTimeLevelTranslating2D.yaml#L50-L85), commit 5cd2c98 of 2026-09-02 before 440107f, VOID |
| `averagedPlanes` diverges | DIVERGED at step 4634 (`alphaFTest_averagedPlanes`) | [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3260); the tables are in [`VOID_closedBox_20260902/`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/VOID_closedBox_20260902/README.md), VOID |
| `donorPlane` is the production value | author instruction, "only Gauss upwind worked in rhoLENT" | [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L686-L699), not measured |
| the token comment on the code default is stale | it says the solver defaults to `central`; the code defaults to `donorPlane` | [token](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L686-L689) against [code](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L59-L63), DERIVED (reading at 8867581) |

## Decisions

- `MASS_FLUX_ALPHAF_SOURCE donorPlane`, author instruction 2026-09-01: [[decisions/mass-flux-alphaf-donor-plane]].
- `MASS_FLUX_ALPHAF_TIME_LEVEL new`: the default that every run before 2026-09-02 used ([token](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L817-L837)).

## Open questions

1. No valid comparison of `donorPlane` against `averagedPlanes` or `donorPlaneAdvected` exists on the repaired case.
2. The `donorPlaneAdvected` result that the raw record quotes as "about 1 %" has no source at 8867581; it is not recorded here.

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/mass-flux-projection]], [[concepts/bound-rho]], [[concepts/coupled-face-density-defect]].
- Decisions and retractions: [[decisions/mass-flux-alphaf-donor-plane]], [[decisions/momentum-schemes-bdf2-upwind]], [[retractions/closed-box-translating-droplet]].

## Log

### 2026-09-28
Created.
