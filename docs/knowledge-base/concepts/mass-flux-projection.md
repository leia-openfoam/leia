---
title: "Projecting the mass flux (curl-free correction)"
description: "A curl-free correction of the geometric mass flux that zeroes the discrete mass residual without an auxiliary density; measured once, on the closed box, where it stopped the droplet; the construction survives, its numbers do not"
aliases: []
kind: concept
status: voided
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled:
decided_by:
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massFluxProjection.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H, cases/default.parameter]
sources: [DP 839-869, DP 680-684, STATUS 0, config/massFluxComparison2D.yaml]
---
# Projecting the mass flux (curl-free correction)

> Verdict (2026-09-28). `massFlux { projectMassFlux true; }` corrects the geometric mass flux by the curl-free field `rhoPhi' = snGrad(lambda) |S_f|` with `laplacian(lambda) = R`, so that `rho` stays geometric and bounded and the discrete mass residual `R = ddt(rho) + div(rhoPhi)` is removed by a Poisson solve instead of an auxiliary density ([`massFluxProjection.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massFluxProjection.H#L1-L43)). The token `MASS_FLUX_PROJECT` is `false` ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L839-L869)). Its one study, `massFluxComparison2D` (2026-09-02, commit 86f8333), ran on the closed-box `translatingDroplet2D` before 440107f and is void ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L680-L684), [[retractions/closed-box-translating-droplet]]). On that box the projection cut the residual 24x and stopped the droplet: travelled fraction 0.725 to 0.049 ([`massFluxComparison2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L49-L90)). Nothing has been measured on the repaired case.

## What it is

One equation per cell for a face field is under-determined. The minimal-L2 member of the family of corrections is the curl-free one: a scalar potential `lambda` with `laplacian(lambda) = R`, `[lambda] = kg/(m s)`, which is to the mass flux what the pressure is to the volumetric flux in the projection step ([`massFluxProjection.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massFluxProjection.H#L20-L31), DERIVED). The correction runs on every outer corrector before the auxiliary equation and `div(rhoPhi,U)` read the flux ([`updateMassFlux.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H#L19-L25)). The case needs a `massFluxPotential` entry in its solvers dictionary; the potential has zero-gradient patches everywhere, so the correction vanishes where the flux is already right ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L167-L213)).

Compatibility. The problem is all-Neumann, so it is solvable only when the volume integral of `R` is zero. That integral is the rate of change of the total mass. The mean of `R` is subtracted before the solve and reported as `massProjGlobalDrift`; the projection therefore fixes the distribution of the imbalance and leaves its global part alone ([same](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massFluxProjection.H#L33-L43), DERIVED). The four diagnostics are `massProjResidualBefore`, `massProjResidualAfter`, `massProjGlobalDrift` and `massProjCorrFraction`.

## Why it failed, or why we think so

On the closed box the projection enforced the balance against the geometric density of a level set that already moved too slowly (travelled fraction 0.725 before the projection). It then reduced the mass flux to match `ddt(rho_geom)`, which locked the wrong interface motion in: a positive feedback on the lag, not a correction of it ([`massFluxComparison2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L66-L72), HYPOTHESIS on a void run). The 47x reduction of the reducible residual that the 4-rank gate reported was retracted the same day as evidence of anything useful ([same](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L74-L76)). The slow interface motion itself was the annihilated free stream of the closed box ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L54)), so the mechanism named there was read off a wrong setup and is not a finding.

The inflation of +5 to +8 % that the token comment gives as the global drift ([token](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L859-L865)) came from `translatingClearOutlet2D`, also pre-fix, and was retracted on 2026-09-27 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L260-L269)).

## Evidence

| claim | number | where |
|---|---|---|
| the projection zeroes the reducible residual and stops the droplet (N = 128) | geometricFaceDensity: residual 2.158e-01 to 9.048e-03, velocity excess 0.6876 to 0.1396, volume error 0.05977 to 0.02448, travelled fraction 0.725 to 0.049; rhoLENT with projection: travelled fraction -0.149 | [`massFluxComparison2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L49-L64), VOID (closed box) |
| the 4-rank gate's 47x reduction for an 18 % flux correction | retracted as evidence | [same](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L74-L76), VOID |
| the compatibility condition and what the mean subtraction removes | the global mass drift | [`massFluxProjection.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/massFluxProjection.H#L33-L43), DERIVED |

## Decisions

- `MASS_FLUX_PROJECT false`; no decision note exists, because no valid measurement exists. The construction stays selectable and inert.

## Open questions

1. A measurement on the repaired case, at density ratio 838.8 and 1, with the residual, the drift and the whole error vector. After dec002f the residual is not the dominant term ([[concepts/density-ratio-amplifier]]), so the expected gain is small; the study is cheap and the question is closed only by a run.

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/bound-rho]], [[concepts/alphaf-source-donor-plane]], [[concepts/density-ratio-amplifier]].
- Retractions: [[retractions/closed-box-translating-droplet]], [[retractions/mass-momentum-consistency-dominant-term]].
- Method: [[concepts/wrong-setup-voids]].

## Log

### 2026-09-28
Created.
