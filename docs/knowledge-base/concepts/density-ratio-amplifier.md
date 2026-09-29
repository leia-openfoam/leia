---
title: "The density ratio is the amplifier, not the mass-flux consistency"
description: "dec002f (2026-09-02): nine orders of magnitude in the mass residual moved the velocity excess by less than 2x, all eight ratio-1 arms complete and all eight water/air arms diverge; the source is the curvature error and the density ratio scales both factors"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-02
decided_by: [config/translatingRepaired2D.yaml, config/translatingRepairedEqualRho2D.yaml, config/kickOriginGate2D.yaml, config/amplifierGate2D.yaml]
code: [cases/translatingDroplet2D, cases/default.parameter]
sources: [STATUS 0, SL article sec:translating, STATUS 11.14, dec002f]
---
# The density ratio is the amplifier, not the mass-flux consistency

> Verdict (2026-09-28). On the repaired translating droplet (N = 128, sixteen arms, `MASS_FLUX` crossed with `MOMENTUM_DIV_SCHEME` at density ratio 838.8 and at ratio 1 with the density sum held at 999.39) all eight ratio-1 arms complete 13334 steps and all eight water/air arms diverge between 8427 and 9987 steps. In the water/air half, where the ratio is fixed and only the mass residual varies, `geometricFaceDensity` raises the residual by 757x (upwind) to 9.4e9x (limitedLinearV) over rhoLENT, and the velocity excess moves by 1.19x, 0.57x and 1.13x: less than a factor of two, not consistently in sign ([dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [SL article, sec:translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2029)). Mass-momentum consistency is therefore not the dominant term ([[retractions/mass-momentum-consistency-dominant-term]]). The source of the current is the curvature error: the one-step kick is independent of `U0` to 0.6 % and an exact curvature removes it by six orders; the density ratio scales the source (factor 40.5 at matched viscosity) and the stability outright ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L74-L129), [[concepts/parasitic-current-mechanism]]).

## What it is

In the two-factor form `L1(U - U0)(T) = u0 exp(G T)` the source `u0` is the first-step disturbance and `G` the growth rate ([CLAUDE.md, research loop step 1](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L413-L420)). Three gates on the repaired case placed the two factors:

1. `kickOriginGate2D` (8 arms, 200 steps): `u0` does not depend on `U0` (2.3338e-04 to 2.3198e-04 of `U_ref` across a fourfold change), so it is not the consistency term `U0 * R`, which is zero at `U0 = 0` and would scale linearly. With `constantCurvatureSurfaceTension` at `kappa = 1/R` the kick falls from 2.15e-03 to 1.69e-09 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L79-L102)).
2. `amplifierGate2D` (8 arms, 8000 steps): the exact-curvature arms stay bounded at every combination of `U0` and ratio; the reconstructed-curvature arm at `U0 = 0.05` and ratio 838.8 reaches 3.0e-02 of `U0`, a factor of 1000 above its exact-curvature twin, and is the arm that destabilises ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L228), [[concepts/well-balanced-exact-curvature-gate]]).
3. The 16-arm matrix (dec002f): the ratio changes the disturbance level by a small factor and the stability outright.

The equal-density control holds `rho1 + rho2` at 999.39, so both halves sit at the same fraction of the capillary time-step limit ([`translatingRepairedEqualRho2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/translatingRepairedEqualRho2D.yaml#L24-L36)). At ratio 1 the two mass-flux models coincide bit for bit, which validates the implementation ([dec002f](https://github.com/leia-openfoam/leia/commit/dec002f)).

## Why it matters

The finding sets the target of the method: the capillary force, not the mass flux ([[concepts/variational-capillary-force]]). It also corrects an earlier reading: "density-ratio independent" was asserted from a matrix that varied only `U0` at one ratio and was retracted on 2026-09-03; the same curvature error divided by the light-phase density is a larger velocity, so the ratio acts on `u0` too ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L104-L118)). And a first curation reported an 18x ratio between the halves; it had read `L_inf` at the common step 8427, inside the water/air blow-ups. Read at step 5000, where every arm is healthy, the `L1` ratio is 1.4 to 1.6 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L119-L125), [[concepts/error-vector-and-read-out-instants]]).

## Evidence

| claim | number | where |
|---|---|---|
| the kick is independent of `U0` | `L1(U-U0)/U_ref` at step 1: 2.3338e-04 (0), 2.3305e-04 (0.0125), 2.3271e-04 (0.025), 2.3198e-04 (0.05 m/s) | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L79-L96), [SL article, tab. kickorigin](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1966-L1989), MEASURED |
| an exact curvature removes it | 2.15e-03 to 1.69e-09 (factor 1.27e+06); `kErrL2Band` 70.6 against 1000 (7 % at R/h = 12.8) | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102), MEASURED |
| the ratio scales the source | ratio 1 at matched dynamic viscosity: 5.769e-06 against 2.334e-04, factor 40.5 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L104-L112), MEASURED |
| the ratio decides the stability | 8 of 8 ratio-1 arms complete 13334 steps; 8 of 8 ratio-838.8 arms diverge at 8427 to 9987 | [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [`translatingRepairedMatrix.csv`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/translatingRepairedMatrix.csv), MEASURED |
| the two halves at a healthy instant | step 5000: `L1` 7.8e-03 to 9.8e-03 against 5.7e-03 to 6.3e-03 (1.4 to 1.6x); volume 2.6 to 3.0e-3 against 7.8 to 8.3e-4; shape 2.2 to 2.7e-5 against 4.2 to 5.9e-6; travelled fraction 1.017 to 1.021 against 0.9955 to 0.9968 | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1991-L2014), MEASURED |
| the residual is not the term | `R` 757x / 9.4e9x / 7.1e9x for upwind / limitedLinearV / vanLeerV; `max|U-U0|` 1.19x / 0.57x / 1.13x | [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), MEASURED |
| exact curvature over the full horizon | plateau 3e-05 of `U0` at ratio 838.8; round-off at ratio 1; reconstructed 3.0e-02 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L199-L222), MEASURED |
| with the fixed binaries at N = 100 | ratio 1 (499.695 both) COMPLETES 0.1 s; ratio 838.8 diverges at 0.0868 s | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3725-L3744), MEASURED |
| ratio 1 is stable, not clean | Popinet's benchmark at ratio 1: 4.86 % of `U` at N = 64; the ratio-1 matrix arms grow from 1.07e-02 at t = 0.063 s to 4.55e-01 at 0.075 s and the droplet reaches 1.66 `U0` before its leading edge passes the outlet | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L145-L152), [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), MEASURED |

## Decisions

- The retraction of the consistency hypothesis: [[retractions/mass-momentum-consistency-dominant-term]]. `MASS_FLUX rhoLENT` stays, on the stationary evidence ([[decisions/mass-flux-rholent]]).

## Open questions

1. The equal-density matrix holds the kinematic viscosities fixed, so its dynamic viscosity ratio moved from 54.8 to 15.3; the matched-mu control exists only in `amplifierGate2D` ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L127-L129), [L225-L227](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L225-L227)).
2. The ratio-1 arms are not clean at the full horizon; the box and horizon pairing is the open translating question ([[cases/translating-droplet]]).

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/ddt-scheme-pairing-bdf2]], [[concepts/mass-flux-projection]].
- Mechanism: [[concepts/parasitic-current-mechanism]], [[concepts/well-balanced-exact-curvature-gate]], [[concepts/variational-capillary-force]].
- Cases: [[cases/translating-droplet]], [[cases/popinet-translating-droplet]].
- Retractions: [[retractions/mass-momentum-consistency-dominant-term]], [[retractions/closed-box-translating-droplet]].

## Log

### 2026-09-28
Created.
