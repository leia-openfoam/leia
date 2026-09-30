---
title: "rhoLENT: the mass-momentum-consistent flux"
description: "An auxiliary density equation with the momentum face flux, reset from the indicator after the outer loop; production since 2026-09-02 on the stationary evidence only; the translating evidence is void and consistency is not the dominant term"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-02
decided_by: [config/rhoLENTStationary2D.yaml, author decision 2026-09-02]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, cases/default.parameter]
sources: [METHOD 6, METHOD 8.1 row MASS_FLUX, STATUS 0, STATUS 11.13, STATUS 11.14, SL article sec:ns, SL article sec:translating, DP 624-685]
---
# rhoLENT: the mass-momentum-consistent flux

> Verdict (2026-09-28). `rhoLENT` solves an auxiliary density equation with the same face flux that the momentum convection uses, on every outer corrector, and resets the density from the phase indicator after the loop ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L288)). It is the production `MASS_FLUX` since 82ca995 (2026-09-02), on the author's decision ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L656-L657), [[decisions/mass-flux-rholent]]). Its valid basis is one study: on the stationary droplet it moves the unabsorbed capillary residual by +1.0 / -22 / +0.1 % at N = 32 / 64 / 128 with volume and shape equal to three digits ([`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L43-L80)), and its relative mass residual is about 1e-13 up to every crash ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L300-L302)). Every translating number of the family before 440107f ran on a closed box and is void ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L73), [[retractions/closed-box-translating-droplet]]). After the fix, dec002f moved the mass residual by nine orders and the velocity excess by less than a factor of two ([SL article, sec:translating](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2029)): the consistency is necessary, but it is not what sets the translating instability ([[concepts/density-ratio-amplifier]], [[retractions/mass-momentum-consistency-dominant-term]]).

## What it is

The one-field momentum equation carries `div(rhoPhi,U)`. For a uniform stream `U == U0` the transient and the convection collapse to `U0 * R`, with `R = ddt(rho) + div(rhoPhi)` the discrete mass residual ([SL article, eq. massmomentum](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L890-L905)). `R` is not a gradient, so the pressure projection cannot absorb it, and it is zero at `U0 = 0`.

rhoLENT (Liu et al., JCP 493, 2023, eq. 40 and Algorithm 1, [DOI 10.1016/j.jcp.2023.112426](https://doi.org/10.1016/j.jcp.2023.112426)) drives `R` to round-off with three steps per outer corrector ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L288)):

1. The face density comes from the reconstructed plane of the donor cell: `rho_f = alpha_f rho1 + (1 - alpha_f) rho2` ([[concepts/alphaf-source-donor-plane]]).
2. The auxiliary equation `(rho_c^{n+1} - rho_c^n)/dt + (1/V_c) sum_f rho_f^{n+1} F_f = 0` is solved with the outer-iteration flux `F_f`, the same flux that forms `rhoPhi = rho_f * phi` in the momentum convection ([`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L1-L17)). The matrix is diagonal; the case needs `"rho.*" { solver diagonal; }` ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L640-L641)).
3. After the outer loop, `rho` is reset to `alpha rho1 + (1 - alpha) rho2` from the Detrixhe-Aslam indicator, so that `rho.oldTime()` of the next step is geometric ([`resetRhoLENT.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H#L1-L44)).

The auxiliary density is clipped to `[rho2, rho1]` right after its solve when `boundRho` is true, which is the production value; the clip is measured, never silent ([[concepts/bound-rho]]). The paper's algorithm is an explicit two-level update; leia uses `backward` for both `ddt(rho,U)` and `ddt(rho)`, because the pairing must match ([[concepts/ddt-scheme-pairing-bdf2]], [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L910-L929)).

## Why it matters

The model decides whether the translating droplet inherits a momentum source `U0 * R`. Without a density equation (`geometricFaceDensity`) `R` is whatever the geometry leaves: 0.04 to 0.56 relative on the stationary droplet, where it does no harm because `U0 = 0` ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L237-L241)). One recorded tension stays open: the level set is advected by the semi-Lagrangian trace velocity, and `rho` is transported by the face flux `phi`. These are two transport operators, and rhoLENT asks for one ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L649-L655), HYPOTHESIS).

## Where in the code

1. Selection and switches: [`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L8-L31).
2. The equation, the clip and the residual: [`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L1-L123); called on every outer corrector from [`updateMassFlux.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H#L27-L30).
3. The reset: [`resetRhoLENT.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H#L1-L44). The mean reset distance `rhoLENTResetL1` goes to the log.
4. Shared with the Eulerian two-phase solver since c094bd8 ([[concepts/eulerian-solver-mass-flux-port]]).

## Evidence

| claim | number | where |
|---|---|---|
| neutral to better on the stationary droplet | unabsorbed capillary residual +1.0 / -22 / +0.1 % at N = 32 / 64 / 128; volume and shape equal to three digits; six arms complete; the control reproduces the shared ladder to five digits | [`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L43-L80), MEASURED |
| the residual of the two models on the stationary droplet | rhoLENT 3.195e-08 / 9.829e-09 / 2.437e-05 relative; geometricFaceDensity 5.072e-01 / 3.955e-02 / 5.590e-01 | same, MEASURED |
| the clip on the stationary droplet | clipL1 2.13e-12 to 2.86e-12 | same, MEASURED |
| relative mass residual over whole runs | about 1e-13, up to every crash | [METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L300-L302), MEASURED |
| nine orders in `R` move the excess by less than 2x | upwind: `R` 757x, `max|U-U0|` 1.19x; limitedLinearV: 9.4e9x, 0.57x; vanLeerV: 7.1e9x, 1.13x | [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2015-L2029), MEASURED |
| the two models coincide at density ratio 1 | max relative difference 0.000e+00 in five metrics, four scheme pairs | [dec002f](https://github.com/leia-openfoam/leia/commit/dec002f), MEASURED |
| the Eulerian solver with the shared flux | travelled fraction 1.00012, relative residual 1.2e-10 at ratio 840 | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3681-L3711), MEASURED |
| translating gain of rhoLENT: travelled fraction 0.874 against 0.725, volume error 0.0024 against 0.0598 at N = 128 | closed box | [`massFluxComparison2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/massFluxComparison2D.yaml#L49-L90), [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L680-L684), VOID |
| the 4-rank gate: 2.296e-11 against 4.7944e-02 relative over 167 steps | closed box | [`rhoBoundGate2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoBoundGate2D.yaml#L32-L61), [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L55-L62), VOID |
| oscillating N = 128: rhoLENT completes where geometricFaceDensity reaches 98.2 %; shape 1.26e-5 against 2.51e-4 | one rung; the case used an algebraic psi at that time | [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L649-L651), [STATUS 11.4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3271-L3280), MEASURED with the void decision open |

## Decisions

- `MASS_FLUX rhoLENT` (82ca995, 2026-09-02): [[decisions/mass-flux-rholent]]. The stationary study is the gate; the translating rationale in the token comment is void ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L387)).
- The record was corrected on 2026-09-27: METHOD 6 said the auxiliary density is never clipped; it is ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L289-L299)).

## Open questions

1. No valid measurement shows a translating benefit after 440107f ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3578-L3595)).
2. The basis of `boundRho` is void ([[concepts/bound-rho]]).
3. The two-operator tension (trace velocity for psi, `phi` for rho) is recorded, not measured.
4. Every SL two-phase result on more than one rank before 28d13f0 carries the coupled-face error ([[concepts/coupled-face-density-defect]]).

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/bound-rho]], [[concepts/alphaf-source-donor-plane]], [[concepts/mass-flux-projection]], [[concepts/ddt-scheme-pairing-bdf2]], [[concepts/eulerian-solver-mass-flux-port]], [[concepts/coupled-face-density-defect]], [[concepts/density-ratio-amplifier]].
- Cases: [[cases/translating-droplet]], [[cases/stationary-droplet]].
- Decisions and retractions: [[decisions/mass-flux-rholent]], [[decisions/mass-flux-bound-rho]], [[retractions/closed-box-translating-droplet]], [[retractions/mass-momentum-consistency-dominant-term]].

## Log

### 2026-09-28
Created.
