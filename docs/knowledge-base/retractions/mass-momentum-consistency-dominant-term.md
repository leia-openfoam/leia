---
title: "Mass-momentum consistency as the dominant term of the translating-droplet failure (falsified 2026-09-02)"
description: "FALSIFIED 2026-09-02 - the reading that the momentum source U0 times the mass residual ddt(rho) + div(rhoPhi) drives the late failure of the translating droplet: on the repaired mesh nine orders of magnitude in the residual (1.7e-11 to 1.6e-01) move the velocity excess by less than 2x and not consistently in sign; the density ratio is the mechanism (dec002f)"
aliases: [mass-momentum consistency falsified, dec002f, U0 times the mass residual]
kind: retraction
status: retracted
part: mass-flux
tags: [retraction, part/mass-flux]
date: 2026-09-28
code: [config/translatingRepaired2D.yaml, config/translatingRepairedEqualRho2D.yaml, config/kickOriginGate2D.yaml, cases/default.parameter]
sources: [STATUS 0, METHOD 6 and 8.1 row MASS_FLUX, DP MASS_FLUX and MASS_RESIDUAL_DIAGNOSTIC and RHO_DDT_SCHEME comments, SL article sec:translating, commit dec002f]
---
# Mass-momentum consistency as the dominant term of the translating-droplet failure (falsified 2026-09-02)

> FALSIFIED 2026-09-02 (commit [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), 36 minutes after the closed-box fix). The claim was "for a droplet translating at U0 the momentum transient collapses to U0 times [ddt(rho) + div(rhoPhi)], not a gradient, so the pressure cannot absorb it, and identically zero at U0 = 0; the mass residual must vanish (only `rhoLENT` does that, 2.3e-11 against 4.8e-02 for `geometricFaceDensity`) and `ddt(rho,U)` must use the scheme of `ddt(rho)`" ([the old STATUS section 0](https://github.com/leia-openfoam/leia/blob/3dece5a/STATUS.md#L31-L39), [DP `MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L668-L672)). The measurement: `translatingRepaired2D`, 16 arms on the repaired mesh, `MASS_FLUX` (rhoLENT, geometricFaceDensity) crossed with four `MOMENTUM_DIV_SCHEME` values at the density ratios 838.824 and 1. In the ratio-838.8 half, where the ratio is fixed and only the residual varies, the discrete mass residual moves from 7.85e-05 to 5.94e-02 (757x) with upwind, from 1.68e-11 to 1.57e-01 (9.4e+09x) with limitedLinearV and from 2.10e-11 to 1.49e-01 (7.1e+09x) with vanLeerV, while max|U-U0| moves 1.19x, 0.57x and 1.13x ([commit dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [SL article sec:translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2028) `sec:translating`). Nine orders of magnitude in the residual move the velocity excess by less than a factor of two, not consistently in sign; a source 757x to 9.4e9x larger would show, and it does not. All eight ratio-1 arms complete 13334 steps and all eight ratio-838.8 arms diverge at steps 8427 to 9987: the density ratio is the mechanism ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L67-L72), [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L113-L117)). Scope: the mechanism. The claim's own numbers were already void as closed-box results ([[retractions/closed-box-translating-droplet]]); `rhoLENT` stays the default on its stationary evidence.

## The claim, and where it lived

- [The old STATUS section 0](https://github.com/leia-openfoam/leia/blob/3dece5a/STATUS.md#L14-L75), `STATUS.md` at commit 3dece5a: "Why" and "both conditions are necessary" (closed-box numbers).
- [`cases/default.parameter`, `MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L668-L672): "WHY IT WORKS, in one line" (the 4.8e-02 residual, now marked history at [d1e3414](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L675-L676)); [`MASS_RESIDUAL_DIAGNOSTIC`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L805-L812): "the quantity the mass-momentum consistency argument turns on"; [`RHO_DDT_SCHEME`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L938-L947): the matching argument and the pairing tables.
- [METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L301), `METHOD.md`: the rhoLENT formulation; [METHOD 8.1, row `MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L387) carries the CORRECTED 2026-09-27 note that names dec002f.
- The SL article, [eq:massmomentum](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L895-L896) `eq:massmomentum`: the derivation of the term, kept as the reason the residual is a controlled quantity.
- [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L237-L241): the stationary residual of `geometricFaceDensity` (0.04 to 0.56) is "now the motivation for the repaired matrix rather than a result about it".

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The velocity excess follows the mass residual. | At fixed ratio 838.8: residual 7.85e-05 to 5.94e-02 (757x), excess 1.19x; 1.68e-11 to 1.57e-01 (9.4e+09x), excess 0.57x; 2.10e-11 to 1.49e-01 (7.1e+09x), excess 1.13x. | MEASURED, [commit dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2024) |
| The failure is a translation effect. | The step-1 kick, `L1(U-U0)/U_ref`, is 2.3338e-04 at U0 = 0 and 2.3198e-04 at U0 = 0.05 (0.4 % over a fourfold change); the term U0 times the residual is zero at U0 = 0 and would scale linearly. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L79-L95), [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1944-L1952) |
| The source of the disturbance is the mass flux. | Exact curvature (kappa = 1/R) makes the force an exact discrete gradient and cuts the kick from 2.15e-03 to 1.69e-09, a factor 1.27e+06; the curvature error is 7 % at R/h = 12.8. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102) |
| The density ratio is not the mechanism. | All eight ratio-1 arms reach 13334 steps; all eight ratio-838.8 arms diverge at 8427 to 9987; read at the common step 8427, max\|U-U0\| is 1.07e-02 to 1.33e-02 at ratio 1 against 2.26e-01 to 4.45e+00. The two `MASS_FLUX` models are bit-identical at ratio 1 (0.000e+00), as they must be. | MEASURED, [commit dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L113-L117) |
| The momentum convection scheme decides. | Upwind, limitedLinearV, vanLeerV and linearUpwind lie within 30 % of one another at step 5000, ordered by numerical diffusion. | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2024-L2028) |
| A first curation read 18x between the ratios. | It was L_inf, sampled at the common step 8427 inside the ratio-838.8 blow-ups; at step 5000 the ratio is 1.4 to 1.6 in L1, and the light-phase argument that predicted 420 is dead. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L119-L125) |

## What survives

1. `rhoLENT` as the `MASS_FLUX` default, on the stationary droplet: +1.0 / -22 / +0.1 % at N = 32 / 64 / 128, volume and shape to three digits ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L233-L236), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L387)), see [[decisions/mass-flux-rholent]] and [[concepts/rholent-mass-flux]].
2. Mass-momentum consistency "remains necessary for a defensible discretisation, and it is what makes the residual a controlled quantity rather than an unknown one" ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2022-L2024)); the measured residual is about 1e-13 up to every crash ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L300-L301)).
3. The matching argument for `RHO_DDT_SCHEME backward` with `MOMENTUM_DDT_SCHEME backward`: with U = U0 the transient reduces to U0 times the residual only if both `ddt` terms use the same scheme ([DP](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L945-L947)), see [[concepts/ddt-scheme-pairing-bdf2]].
4. The replacement: the curvature estimator is the source and the density ratio acts on both factors of max|U|(T) = u0(h) exp(G(h)) ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L104-L117)); with exact curvature the instability does not survive at any ratio ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L227)), see [[concepts/density-ratio-amplifier]] and [[concepts/parasitic-current-mechanism]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: [section 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L67-L133) (the repaired matrix, the kick origin, the decomposition) and the [amplifier gate](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L227).
- [x] `METHOD.md`: the [8.1 row `MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L387) (CORRECTED 2026-09-27: "dec002f falsified mass-momentum consistency as the dominant term").
- [x] `cases/default.parameter`: the [CORRECTED 2026-09-27 note under `MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L680-L684).
- [x] The SL article: ["Mass-momentum consistency is not the mechanism here"](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2028) and ["the source is the curvature error"](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1944-L1961).
- [x] The curated table `translatingRepairedMatrix.csv` ([blob](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/translatingRepairedMatrix.csv)).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The [`MASS_RESIDUAL_DIAGNOSTIC` comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L805-L812) still states the U0 times residual mechanism without a marker.
- [ ] Ratio 1 is not clean at the full horizon either: max|U-U0| grows from 1.07e-02 at t = 0.063 s to 4.55e-01 at 0.075 s, and past t = 0.081 s the leading edge passes the outlet at 10 mm ([commit dec002f](https://github.com/leia-openfoam/leia/commit/dec002f)); the box and horizon question is in [[retractions/late-translating-instability-is-the-outlet]].
- [ ] The residual of `geometricFaceDensity` on the repaired translating case has no curated table of its own; the numbers live in the commit message and in the SL article.

## Related

Hubs: [[hubs/mass-flux]]. Siblings: [[concepts/rholent-mass-flux]], [[concepts/density-ratio-amplifier]], [[concepts/parasitic-current-mechanism]], [[concepts/ddt-scheme-pairing-bdf2]], [[models/mass-flux]], [[decisions/mass-flux-rholent]], [[decisions/momentum-schemes-bdf2-upwind]], [[cases/translating-droplet]], [[retractions/closed-box-translating-droplet]], [[retractions/late-translating-instability-is-the-outlet]].

## Log

### 2026-09-28
Written from STATUS 0, the commit dec002f, METHOD 6 and 8.1, cases/default.parameter and the SL article. Falsified 2026-09-02. Entered in [[retraction-log#2026-09]].
