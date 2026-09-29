---
title: "MASS_FLUX: rhoLENT, the mass-momentum-consistent flux"
description: "MASS_FLUX rhoLENT in the global default since 82ca995 (2026-09-02), decided by config/rhoLENTStationary2D.yaml: the unabsorbed capillary residual moves +1.0 / -22 / +0.1 % at N = 32 / 64 / 128 with volume and shape equal to three digits; every translating rationale before the closed-box fix is void"
aliases: [MASS_FLUX rhoLENT]
kind: decision
status: settled
part: mass-flux
tags: [decision, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-02
decided_by: [config/rhoLENTStationary2D.yaml, "author decision 2026-09-02"]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H, cases/default.parameter]
sources: ["METHOD 6 (L276-L320)", "METHOD 8.1 row MASS_FLUX (L386)", "STATUS 0 (L14-L73, L229-L247)", "STATUS 11.13 (L3596-L3600)", "DP L660-L689", "SL article sec:rholent and L2186-L2196"]
---
# MASS_FLUX: rhoLENT, the mass-momentum-consistent flux

> `MASS_FLUX rhoLENT` in the global default since 82ca995 (2026-09-02), on the author's decision ([DP L660-L689](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L660-L689)). Decided by `config/rhoLENTStationary2D.yaml` ([config L42-L62](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L42-L62)): on the 2D stationary droplet rhoLENT is neutral to better against `geometricFaceDensity`, the unabsorbed capillary residual moves +1.0 / -22 / +0.1 % at N = 32 / 64 / 128, volume and shape agree to three digits, all six arms reach the horizon, and the clip is at 2e-12 ([METHOD 8.1 L386](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L386)). The translating evidence the token comment quotes (travelled fraction 0.874 against 0.725) ran on the closed box and is void ([DP L684-L688](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L684-L688), [[retractions/closed-box-translating-droplet]]). After the fix, dec002f falsified mass-momentum consistency as the dominant term of the translating instability: nine orders of magnitude in the mass residual moved the velocity excess by less than 2x ([SL article L2186-L2196](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2186-L2196), [[retractions/mass-momentum-consistency-dominant-term]]).

## The question

The one-field momentum equation carries `rhoPhi` in `div(rhoPhi, U)`. For a droplet that translates at a uniform `U0`, the momentum transient plus convection collapses to `U0 R`, with `R = ddt(rho) + div(rhoPhi)` the discrete mass residual; that source is not a gradient and vanishes at `U0 = 0` ([METHOD 6 L276-L288](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L288), [[concepts/rholent-mass-flux]]). rhoLENT solves an auxiliary density equation with the same face flux on every outer corrector and resets `rho` from the phase indicator after the loop, so `R` is about 1e-13 ([METHOD 6 L300-L302](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L300-L302)). The question of 2026-09-02 was the regression: can rhoLENT be adopted without damage to the stationary droplet, the case the campaign had made clean ([config L1-L24](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L1-L24))?

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| `geometricFaceDensity`, N = 32 / 64 / 128 (control) | relative L2 capillary residual, second half | 2.3693e-05 / 3.1449e-06 / 8.5341e-07; the shared ladder reproduced to five digits | MEASURED, [config L44-L56](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L56) |
| `rhoLENT`, N = 32 / 64 / 128 | the same | 2.3939e-05 / 2.4505e-06 / 8.5436e-07: +1.0 / -22 / +0.1 % | MEASURED, [config L44-L62](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L62) |
| both, N = 32 / 64 / 128 | volume error; shape L2 | equal to three digits (rhoLENT 6.656e-04 / 1.676e-04 / 1.289e-05; 3.268e-06 / 2.375e-06 / 2.262e-07) | MEASURED, [config L44-L52](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L52) |
| both | relative mass residual | geometricFaceDensity 5.072e-01 / 3.955e-02 / 5.590e-01; rhoLENT 3.195e-08 / 9.829e-09 / 2.437e-05 | MEASURED, [config L44-L52](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L52) |
| `rhoLENT` | density clip L1 | 2.13e-12 / 2.56e-12 / 2.86e-12 | MEASURED, [config L44-L52](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L52) |
| repaired translating matrix (dec002f), ratio 838.8 | mass residual between the two fluxes; velocity excess | 757x (upwind), nine orders 1.7e-11 to 1.6e-1 (limited schemes); `L1(u - u0)` moves by less than 2x | MEASURED, [SL article L2186-L2196](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2186-L2196) |
| translating, closed box: travelled fraction 0.874 against 0.725, volume 0.0024 against 0.0598 | — | VOID | [DP L666-L669](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L666-L669), [STATUS L55-L62](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L55-L62) |

Pre-registered read-out ([config L26-L40](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L26-L40)): the control must reproduce the shared ladder; PASS = rhoLENT within a small factor at every rung and all arms complete; neutral is the success criterion, a large improvement would be as suspect as a regression.

The same table shows the mechanism: `geometricFaceDensity` carries a mass residual of 0.04 to 0.56 relative on the stationary droplet and does no harm there, because the source is `U0 R` and `U0 = 0` ([config L64-L70](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L64-L70), [STATUS L237-L241](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L237-L241)).

## What it does not cover

1. The translating droplet. No valid measurement shows a translating benefit of rhoLENT after 440107f ([STATUS L3596-L3600](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3596-L3600)); the density ratio is the amplifier, not the consistency ([[concepts/density-ratio-amplifier]]).
2. The density bound. rhoLENT's auxiliary density is unbounded, so `MASS_FLUX_BOUND_RHO true` is active with it; its measured basis is void ([[decisions/mass-flux-bound-rho]]).
3. The face density on coupled faces was rank-local until 28d13f0 (2026-09-27): every result on more than one rank before that commit carries a seam error of 1e-5 to 5e-4 relative in the velocity metrics ([METHOD 6 L304-L312](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L304-L312), [[concepts/coupled-face-density-defect]]). The deciding study ran at np 4.
4. `interpolatedDensity` has never been measured ([[models/mass-flux]]).
5. The Eulerian two-phase solver reads the same family only since c094bd8 ([[concepts/eulerian-solver-mass-flux-port]]).

## Related

[[hubs/mass-flux]] - [[models/mass-flux]] - [[concepts/rholent-mass-flux]] - [[concepts/density-ratio-amplifier]] - [[concepts/bound-rho]] - [[concepts/coupled-face-density-defect]] - [[concepts/eulerian-solver-mass-flux-port]] - [[decisions/mass-flux-alphaf-donor-plane]] - [[decisions/mass-flux-bound-rho]] - [[decisions/momentum-schemes-bdf2-upwind]] - [[retractions/closed-box-translating-droplet]] - [[retractions/mass-momentum-consistency-dominant-term]] - [[cases/stationary-droplet]] - [[cases/translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the stationary gate of 2026-09-02; the translating rationale is void since 2026-09-27. Entered in [[decision-log#2026-09]].
