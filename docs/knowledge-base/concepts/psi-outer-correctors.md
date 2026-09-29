---
title: "psiOuterCorrectors: re-advect psi in every outer iteration"
description: "Since 2026-08-28 the coupled solver restores psi^n and re-advects it with the current velocity iterate on every PIMPLE outer corrector; the frozen-force lag it removes is exonerated as the growth mechanism (r 202.8 against 198.9 1/s in 2D, +0.4 and -0.1 percent in 3D), so the default is a consistency decision, not a stability gain."
aliases: [psiOuterCorrectors, PSI_OUTER_CORRECTORS, iterated interface-force coupling, frozen-force lag]
kind: concept
status: settled
part: advection
tags: [concept, part/advection]
date: 2026-09-28
date_settled: 2026-08-28
decided_by: [config/stationaryDropletIteratedCoupling.yaml, config/stationaryDropletCouplingDepth.yaml, config/psiOuterCorrectorsGain3D.yaml, config/outerLoopControl2D.yaml, "author decision 2026-08-28, cases/default.parameter L403-L409"]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H]
sources: ["DP L402-L449", "METHOD 6 correction (L314-L317)", "METHOD 8.1 rows PSI_OUTER_CORRECTORS and N_OUTER_CORRECTORS (L394-L395)", "PCS 16.2 (L1352-L1363)", "PSH gate result (L207-L236)", "STATUS 4 full-horizon gate (L1089-L1096)", "STATUS 4 sequencing (L636-L640)", "SL article sec:limitations (L2491-L2495)"]
---
# psiOuterCorrectors: re-advect psi in every outer iteration

> Verdict (2026-09-28). With `PIMPLE.psiOuterCorrectors yes` every outer corrector restores `psi^n` from a snapshot and re-advects it with the current velocity iterate, then re-runs the band, the phase indicator, the curvature and the face delivery, so the capillary force converges to its `t^{n+1}` value within the step ([createTransportFields.H L211-L224](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L211-L224)). The switch is the default since 2026-08-28, on consistency grounds ([DP L403-L409](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L403-L409), [METHOD L314-L317](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L314-L317)). The mechanism it removes, the once-per-step frozen-force lag, is exonerated as the growth carrier: at N = 128 the 2D growth rate is 202.8 1/s iterated against 198.9 frozen, +2 %, and 1, 3 or 6 outer correctors, frozen or iterated, all give r in [185, 203] 1/s ([PCS 16.2 L1352-L1363](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1352-L1363)). In 3D the per-step gain moves by +0.4 % at R/h = 10.0 and -0.1 % at R/h = 15.8 ([PSH L209-L220](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L209-L220)). The m = 2 mode rate is inert to four digits: 3.906 against 3.908 at dt, -1.267 against -1.269 at dt/2 ([DP L429-L435](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L429-L435)). The cost is one transport, band, indicator, curvature and delivery chain per extra outer pass.

## What it is

The momentum left-hand side is `fvm::ddt`, implicit Euler with every term at `t^{n+1}`. With the switch off, the force `sigma kappa_f snGrad(alpha)` is evaluated once per step from the AB2-predicted trajectory and then frozen across the outer correctors ([createTransportFields.H L215-L219](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L215-L219)). With the switch on, pass 1 keeps the AB2 predictor and every later pass advects `psi^n` by the actual momentum iterate, so at outer-loop convergence the step is backward-Euler consistent in all terms ([L220-L223](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L220-L223)).

The consistency question was asked and answered on 2026-08-28: does re-advecting with the current iterate wrongly feed `u^{n+1}` to a trace that wants `u^n`? No. The kernel's linear displacement is the average of its two slots, so it wants `(u^{n+1}, u^n)`; the AB2 extrapolation on pass 1 is the fallback for `u^{n+1}` not existing yet, and a later pass carrying the real iterate is more consistent, not less ([DP L411-L423](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L411-L423), [[concepts/departure-foot-ab2-centring]]). The failure mode that would break it, refreshing the old slot on later passes, is guarded: the old slot is written on pass 1 only ([slAlphaEqn.H L129-L139](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L129-L139)).

## Why it matters

The switch was built as a lever against the `c dt` term of the growth rate, `r = r0 + c dt`, which carried 90 % of the parasitic growth at N = 128 ([DP L447-L449](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L447-L449), [[concepts/parasitic-current-mechanism]]). The measurements show the frozen-force lag is not that term. One test eliminated two mechanisms at once, because the switch also replaces the AB2 velocity extrapolation on passes >= 2; neither moves the answer by more than 0.4 % ([PSH L222-L230](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L222-L230)). The consequence for the campaign: the semi-implicit capillary force is the same class of fix, and the gate says that class buys nothing here ([PSH L232-L236](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L232-L236), [[models/semi-implicit-capillary-force]]).

## Where in the code

- The switch and the `psiN` snapshot (registered `NO_WRITE`): [createTransportFields.H L225-L240](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H#L225-L240).
- The gate of the interface pipeline, `pimple.firstIter() || psiOuterCorrectors || slSecondHalfDrift`: [slAlphaEqn.H L2-L17](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L2-L17); the snapshot on the first iteration: [L88-L92](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L88-L92); the restore `psi == psiN` on passes >= 2: [L143-L149](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L143-L149).
- The banner prints the setting: [printMethodBanner.H L201-L202](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H#L201-L202).
- Tokens: `PSI_OUTER_CORRECTORS yes`, `N_OUTER_CORRECTORS 3` ([DP L402-L403](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L402-L403)).

## Evidence

| claim | number | where |
|---|---|---|
| Switch off is bit-identical to the pre-switch fixed-dt matrix; switch on runs 3.00 interface pipelines per step | bit-identical through blow-up; 3.00 | MEASURED, [PCS 16.2 L1352-L1359](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1352-L1359) |
| 2D stationary droplet, N = 128, 3 outer correctors: growth rate frozen against iterated | 198.9 against 202.8 1/s (+2 %, noise) | MEASURED, [PCS L1359](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1359), [stationaryDropletCouplingDepth.yaml L3-L8](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDropletCouplingDepth.yaml#L3-L8) |
| Coupling-depth bracket, nOuter 1 / 3 / 6, frozen or iterated | r in [185, 203] 1/s; the outer-coupling axis is flat | MEASURED, [PCS L1359-L1363](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1359-L1363) |
| 3D stationary droplet, quarter horizon, max abs U at the end, no against yes | N_L = 60: 2.7722e-04 against 2.7844e-04 (+0.4 %); N_L = 95: 2.3627e-04 against 2.3600e-04 (-0.1 %) | MEASURED, [PSH L209-L216](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L209-L216) |
| Same gate, per-step gain and e-folds | gAvg +1.164e-04 / +1.807e-04 (N_L 60), +3.828e-04 / +4.082e-04 (N_L 95); volume, shape, min abs grad psi agree to 3 to 4 digits | MEASURED, [PSH L211-L220](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L211-L220) |
| m = 2 mode rate, N = 128, 12 outer correctors, yes against no | 3.906 against 3.908 at dt; -1.267 against -1.269 at dt/2 | MEASURED, [DP L429-L432](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L429-L432) |
| 3 frozen correctors against 12 re-advected, 2D, quarter horizon | peak, t_peak, final value, volume and shape agree to 3 to 4 digits; about 4x cheaper per arm | MEASURED, [fullHorizonStability2D.yaml L16-L20](https://github.com/leia-openfoam/leia/blob/8867581/config/fullHorizonStability2D.yaml#L16-L20), [outerLoopControl2D.yaml L1-L36](https://github.com/leia-openfoam/leia/blob/8867581/config/outerLoopControl2D.yaml#L1-L36) |
| The increment form's deferred correction at 3 outer correctors | contracts about 5x per pass; 0.07 % of the term on the last pass, not round-off | MEASURED, [STATUS L1089-L1096](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1089-L1096) |
| Cost of the switch | one transport, band, indicator, curvature and delivery chain per extra outer pass | DERIVED, [DP L433-L435](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L433-L435) |

## Decisions

- `PSI_OUTER_CORRECTORS yes`, default changed 2026-08-28: if the loop is worth iterating, the interface it iterates towards has to move with it ([DP L403-L409](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L403-L409), [METHOD 8.1 L394](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L394)). The change waited until the kinematic baseline was re-established after the gradU fix ([STATUS L636-L640](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L636-L640)).
- `N_OUTER_CORRECTORS 3`, the historical value ([METHOD 8.1 L395](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L395)).
- The semi-implicit force class was dropped from the Shannon plan on this gate ([PSH L232-L236](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L232-L236)).

## Open questions

1. The loop equivalence was measured with the semi-implicit increment term off; with the term on, the loop decides whether its deferred correction converges ([STATUS L1094-L1096](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1094-L1096)).
2. The SL pre-print's sixth limitation still says the interface is advanced once per step and held fixed across the correctors ([SL article L2491-L2495](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2491-L2495)); METHOD records the correction ([L314-L317](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L314-L317)).
3. The 3D gate's baselines predate the psi-filter seam fix, which is why its control sits inside the study ([psiOuterCorrectorsGain3D.yaml L26-L29](https://github.com/leia-openfoam/leia/blob/8867581/config/psiOuterCorrectorsGain3D.yaml#L26-L29), [[retractions/psi-filter-seam-bug]]).

## Related

[[hubs/advection]] - [[hubs/surface-tension]] - [[concepts/departure-foot-ab2-centring]] - [[concepts/force-time-centring]] - [[concepts/trace-velocity-projected-flux]] - [[concepts/parasitic-current-mechanism]] - [[concepts/capillary-time-step]] - [[concepts/bit-identity-and-inertness-gates]] - [[models/semi-implicit-capillary-force]] - [[cases/stationary-droplet]] - [[studies/curvature-stabilization-campaign]] - [[studies/shannon-parasitic-currents-campaign]] - [[retractions/psi-filter-seam-bug]]

## Log

### 2026-09-28
Created from the token comment, METHOD 6 and 8.1, the plan sections 16.2 (PCS) and the gate result (PSH), and the solver headers.
