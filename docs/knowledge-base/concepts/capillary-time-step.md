---
title: "The capillary time step"
description: "Every droplet runs at dt = 0.010861 h^1.5, 0.2323 of the Brackbill limit and 0.164 of Popinet's, tied to the cell size; the growth rate is about 90 percent dt-proportional at N = 128, a halving of the accumulated e-folds costs 7.2 times the steps, and the coefficient is not raised while the growth stands (2026-09-28)."
aliases: [CAPILLARY_DT_COEFF, Brackbill limit, capillary step, MAX_DELTA_T]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-08
decided_by: [cases/stationaryDroplet2D.parameter, config/stationaryDropletDtSweep.yaml, config/popinet3D_poly_sigma0_dtSweep.yaml, config/gates/methodGate2D.yaml]
code: [workflow/scripts/materialize.py, cases/stationaryDroplet2D/system/controlDict.template, cases/stationaryDroplet2D.parameter, cases/default.parameter]
sources: [stationaryDroplet2D.parameter CAPILLARY_DT_COEFF, STATUS 4 domain size, STATUS 4 polyhedral dt sweep, STATUS 4 polyhedral audit, STATUS 11.13, PCS 16.1, PCS 18.4, PSH 0c fixed dt, PSH 0e, PSH 0g, SL article sec:droplet, SL negative deck 1/2, CLAUDE no partial solutions, METHOD 10]
---
# The capillary time step

> Verdict (2026-09-28). The surface-tension force is explicit, so the step is capillary-limited. Every droplet case derives `MAX_DELTA_T = CAPILLARY_DT_COEFF / nRef^1.5` with `CAPILLARY_DT_COEFF = 0.010861` s and `nRef = CAPILLARY_REF_LENGTH/h`, so the step follows the cell size `h` as `h^1.5` and not the cell count ([`cases/stationaryDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L25-L35), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L767-L773)). That is 0.2323 of the Brackbill limit `sqrt((rho_1 + rho_2) h^3/(2 pi sigma))` at every rung ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3593)) and 0.164 of Popinet's `sqrt(rho h^3/(pi sigma))` ([PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L448-L453)). The step is not the cause of the instability and not its cure: at fixed N = 128 the growth rate is 90 percent proportional to `dt` (`r = 20.8 + 2.40e7 dt`), but the accumulated e-folds scale only as `dt^0.35`, so halving them costs 7.2 times the steps ([PCS 16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1315-L1327), [PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L160-L170)); on the polyhedral mesh the failure time does not move at all under `dt`, `dt/2` and `dt/4` ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2194-L2220)). The coefficient is not raised for a speed-up while the growth it would cover stands ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L377-L380)).

## What it is

The limit for an explicitly integrated capillary term is the period of the shortest capillary wave the mesh carries: `dt_sigma = sqrt((rho_1 + rho_2) h^3/(2 pi sigma))` in Brackbill's form, `sqrt(rho h^3/(pi sigma))` in Popinet's ([METHOD 10](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L832-L833), [`cases/stationaryDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L48-L58)). The sum of the densities enters, so a density-ratio sweep must hold `rho_1 + rho_2 = 999.39` fixed, that is 499.695 each at ratio 1, or the safety factor moves with the sweep ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L961-L967)). Writing `h = L_ref/nRef` makes the bound separable, and the coefficient is a property of the fluids and the reference box: for water/air (`rho_1 = 998.2`, `rho_2 = 1.19`, `sigma = 0.07274`, `L_ref = 0.01` m) Popinet's coefficient would be 0.0661311 s ([`cases/stationaryDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L50-L58)). The value in force is 0.010861 s, which the parameter file also records as `6e-5 * 32^1.5`, the step `dt_sigma/4` at N = 32 ([`cases/stationaryDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L60-L73)). Note: the comment above that value says the coefficient was set to Popinet's limit on 2026-08-23, but the value line is 0.010861 and STATUS 11.13 confirms 0.2323 of the Brackbill limit; the comment is stale. The ratios `0.010861/sqrt(999.39/(2 pi 0.07274)) = 0.2323` and `0.010861/sqrt(999.39/(pi 0.07274)) = 0.1642` are derived here and agree with both records.

The step at the reference box: `7.50e-6` s at N = 128 (13334 steps to 0.1 s), `6.0e-5` s at N = 32 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1052-L1054), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1111)); in 3D on the 6R box `1.0861e-5 / 7.6186e-6 / 5.4514e-6 / 3.8399e-6` s for 9207 / 13126 / 18344 / 26042 steps at R/h = 10 to 20 ([`config/stationaryDroplet3Dwide.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3Dwide.yaml#L5-L10)). The step is fixed: `deltaT = MAX_DELTA_T`, `adjustTimeStep no` ([`controlDict.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/controlDict.template#L13-L27)); `FIXED_DELTA_T > 0` overrides the law for a sweep at fixed mesh ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L303-L312)).

## Why it matters

Refinement moves the mesh and the step together, so a ladder cannot separate a spatial effect from a per-step one; that is why the dt sweeps at fixed N and the fixed-dt ladder exist ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L303-L312)). The cost of the explicit step is intrinsic: the step count grows as `N^1.5` and the 2D wall time as about `N^3` (1670 and 4716 steps, 9.3 and 75 s at N = 32 and 64) ([SL article `sec:droplet`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1852-L1859)). The route to a larger step is the semi-implicit term, which lets SAAMPLE run at `omega_grid dt = pi/2` above the explicit wall of 1.0 to 1.3 measured here ([PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L773-L790), [[models/semi-implicit-capillary-force]]).

## Where in the code

- `materialize.py` evaluates `MAX_DELTA_T` at `nRef = CAPILLARY_REF_LENGTH/h`; before that fix a truncated box got `4.29e-5` instead of `1.0861e-5` s at `L = 4e-3`, N = 40 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L767-L773)).
- The tokens `CAPILLARY_DT_COEFF`, `CAPILLARY_REF_LENGTH`, `DOMAIN_LENGTH` in the case layer ([`cases/stationaryDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D.parameter#L18-L35)); the gate states the law `dt = 10.861 h^1.5` ([`config/gates/methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L18-L19)).

## Evidence

| claim | number | where |
|---|---|---|
| the fraction of the limit | 0.2323 of the Brackbill limit at every rung; 0.164 of Popinet's at every R/h | [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3593), [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L453), MEASURED and DERIVED |
| the growth rate at fixed N is dt-proportional | N = 128: r = 202.9 / 104.4 / 70.0 1/s at dt, dt/2, dt/4, r0 = 20.8, c = 2.40e7 (90 percent); N = 256: r0 = 91.7, c = 6.07e7 (65 percent); N = 64: c about 0 | [PCS 16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1315-L1327), MEASURED |
| the mode rate halves with the step | r_2 = 18.8 / 12.7 / 8.01 / 4.02 1/s, intercept +0.03 | [PCS 18.4](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1556-L1563), MEASURED |
| a smaller step improves the interface metrics | dt/16 at N = 128: volume 73x, shape 33x better | [PCS 16.1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1325-L1327), MEASURED |
| r dt looked constant on the 3D ladder | 4.395e-4 / 3.986e-4 / 4.936e-4 per step at R/h = 10 / 12.7 / 15.8 (24 percent) | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L926-L938), MEASURED, superseded |
| the fixed-dt ladder falsified it | at dt = 5.451e-6 the R/h = 12.7 arm goes from +2.60 to -0.40 e-folds; gAvg about dt^1.35, e-folds about dt^0.35; 7.2x cost per halving | [PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L145-L170), MEASURED |
| the polyhedral far-field failure is not a step effect | steps to failure 528 / 1009 / 1974 (1 : 1.91 : 3.74); failure time 4.877 / 4.660 / 4.558 ms; rate 126 to 133 1/s | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2194-L2220), MEASURED; see [[concepts/polyhedral-fit-amplification]] |
| the polyhedral audit | interface cells at 0.236 (Popinet poly) and 0.198 (stationary poly) of the local limit; far-field slivers 1.476 and 1.889 above it; Co 0.0165, nu dt/h^2 0.0062 | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2504-L2523), MEASURED |
| a smaller or adaptive step only delays | t_blow 0.466 s at N = 32 against 0.032 s at N = 256; CST dt/2 moved 0.047 to 0.050 s | [negative deck 1/2](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/1/2), MEASURED |
| the explicit wall and the semi-implicit door | omega_grid dt 1.0 to 1.3 here; SAAMPLE stable at pi/2 with the Laplace-Beltrami term | [PSH 0e](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L581-L583), [PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L773-L783), MEASURED |

## Why it failed, or why we think so

The step is not a failure; the two readings that made it look like a lever failed. "r dt is constant" was `h` and `dt` covarying along the native ladder ([PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L160-L164)), and "smaller dt would cure the fine mesh" runs into `e-folds ~ dt^0.35` ([PSH 0c](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L166-L170)). On polyhedra the amplification bound is `Lambda^(T/dt) = exp((Lambda - 1) T/dt) = exp(c U T)`, which contains no `dt` ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2213-L2220)).

## Decisions

- `CAPILLARY_DT_COEFF 0.010861` in every droplet case layer; the gates carry it in `twoPhaseCoupling` ([`config/gates/methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L81)).
- The step follows `h`, not the cell count; bit-identical for every pre-existing token shape ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L767-L773)).
- Ladders are matched on the time-step law, and `gAvg` is compared only at equal step counts ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L458-L460), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1027-L1028), [[concepts/richardson-ladders-and-orders]]).

## Open questions

1. The stationary polyhedral rung runs 16 percent below its own law (`N_CELLS = 84` describes 7.143e-5 m while cfMesh made 7.937e-5 m); re-pin to 76 before a dt-matched comparison ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2518-L2521)).
2. A case whose interface can reach a wall must pin the step to the smallest band cell, not to the interface cell ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2514-L2517)).
3. The stale comment in `cases/stationaryDroplet2D.parameter` that describes a switch to Popinet's coefficient which the value line does not carry.

## Related

[[hubs/surface-tension]], [[concepts/parasitic-current-mechanism]], [[concepts/force-time-centring]], [[concepts/polyhedral-fit-amplification]], [[concepts/richardson-ladders-and-orders]], [[concepts/density-ratio-amplifier]], [[models/semi-implicit-capillary-force]], [[cases/stationary-droplet]], [[cases/popinet-translating-droplet]], [[decisions/mesh-family-hexahedral]], [[studies/shannon-parasitic-currents-campaign]].

## Log

### 2026-09-28
Created from the case parameter file, STATUS 4 and 11.13, PCS 16 and 18, the Shannon plan sections 0c, 0e and 0g and the negative-results deck.
