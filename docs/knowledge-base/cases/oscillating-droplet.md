---
title: "The oscillating droplet"
description: "Fixed 2D gate, 2026-09-27: read by period and damping, the mode-2 droplet's gradient drift turns unstable at N = 200 within ten periods."
aliases: [oscillating droplet, oscillatingDroplet2D, oscillatingDroplet3D, mode-2 droplet, Lamb oscillation, DROPLET_SURFACE]
kind: case
status: settled
part: surface-tension
tags: [case, part/surface-tension]
date: 2026-09-29
date_settled: 2026-09-27
decided_by: [config/gates/methodGate2D.yaml]
code: [cases/oscillatingDroplet2D, cases/oscillatingDroplet3D, cases/oscillatingDroplet2D.parameter, cases/oscillatingDroplet3D.parameter, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H, workflow/scripts/make_gate_summary.py, docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_result_figures.py]
sources: [STATUS 11.4, STATUS 11.5, STATUS 11.15, SL article sec:limitations, gcls article sec:res-baseline, methodGate2D.yaml arm oscillating, PHL 5.4, METHOD 8.1 row CURVATURE_EXTENSION]
---
# The oscillating droplet

> Verdict (2026-09-29). The droplet starts at rest as a mode-2 ellipse: `a = 1.1 R`, with the area of the circle of radius R = 1 mm. Surface tension makes it oscillate toward the circle ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet2D/system/fvSolution.template#L391-L402)). The gate reads it by the period and the damping rate of the mode-2 coefficient. The inviscid Lamb period of 9.51 ms is a check ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L114-L144)). On the fixed 2D gate the band gradient error at T = 0.1 s is 0.10 and 0.14 at N = 100 and 142. At N = 200 it grows exponentially after t = 0.05 s and reaches 5.13 ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3973-L3988)). That rung completes the horizon close to a divergence. Its period is 8.06 ms, against 10.0 and 9.82 ms on the coarser rungs. The Celik procedure classifies the period and damping ladders as divergent, so the rung cannot carry a verdict (same lines). The SL pre-print now states the drift as a defect of the method itself, with the figure `gcls_oscillating_drift.pdf` ([SL article `sec:limitations`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2742-L2765)). Every oscillating study before 2026-09-26 used an algebraic psi, whose gradient columns have no meaning. The void of those studies is an OPEN author decision ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3289-L3297)).

## What it is

**Geometry and fluids.** The 2D case uses the mesh, the walls, the fluids and the step of the stationary droplet; only the implicit surface differs ([[cases/stationary-droplet]]). The box is 10 mm with no-slip walls, the fluids are water and air, and the step is `dt = 10.861 h^1.5` s ([`oscillatingDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet2D.parameter#L22-L34)). The ellipse has the semi-axes 1.1e-3 m and 0.90909e-3 m, so `a b = R^2`. The metrics use the equivalent radius R for the rest-state curvature 1/R ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet2D/system/fvSolution.template#L398-L401)).

**The 3D twin.** A prolate spheroid, `a = 1.1 R` and `b = c = R/sqrt(1.1)`, so `a b c = R^3`, initialised as an exact signed distance on the 6R box. The ladder is N = 60 / 78 / 102 and T = 0.025 s, about 2.1 periods of Lamb's 3D mode, 11.65 ms ([`oscillatingDroplet3D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet3D.parameter#L1-L21), [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet3D/system/fvSolution.template#L373-L386)).

**The Lamb period.** For mode n = 2 the gate uses the inviscid frequency ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L139-L144)):

- 2D: `omega^2 = (n^3 - n) sigma / ((rho_d + rho_a) R^3)`, which gives `T_Lamb = 2 pi/omega = 9.508` ms and `T_REF = 1/omega = 1.51e-3` s ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L155-L166), [`methodGate2D_baseline_summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_summary.csv#L15-L17)).
- 3D: `omega^2 = n (n - 1)(n + 2) sigma / ((n rho_d + (n + 1) rho_a) R^3)`, 11.65 ms, `T_REF = 1.854e-3` s ([`methodGate3D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate3D.yaml#L112-L123)).

The 2D horizon of 0.1 s holds about 10.5 periods ([plan-halo-limited 5.2](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-halo-limited-gradient-control.md#L505)).

**The read-out.** The solver fits `r - R = mean + a cos(2 theta) + b sin(2 theta)` to the psi = 0 crossings every step ([`writeDropletMetrics.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L156-L183)). The gate then computes ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L114-L136)):

1. the period `T_osc = 2 (t_last - t_first)/(n_cross - 1)` from the zero crossings of `m2CosCoefficient`, with at least three crossings;
2. the damping rate `gamma`, the negative slope of a least-squares line through `ln(peak)` against t, one peak between each pair of crossings;
3. `periodVsLamb = T_osc/T_Lamb - 1`.

Period and damping have no exact value in the discrete problem, so the gate classifies their convergence by the procedure of Celik et al. ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L266-L283), [gcls article `sec:gate-metrics`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L519-L521)). The velocity norms are the physical oscillation: the gate reports them as `oscL2MagU` and `oscMeanMagU` and does not score them ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L214-L224)). It also scores no shape, pressure-jump or curvature error on this arm, because the radial distance to the circle is the oscillation itself ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L225-L233)).

**The level-set surface is a token.** `DROPLET_SURFACE implicitEllipsoid` is the algebraic `psi = sum (x_i - c_i)^2/a_i^2 - 1`, so the gradient magnitude `q = 2/a_i` is 1.8e3 to 2.2e3 at the interface. `signedDistanceEllipse` is a signed distance. The default is still `implicitEllipsoid`, for bit-identity; only the gates pin the signed distance ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1029-L1033), [`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L165)).

## Why it matters

In this coupled case the interface moves and deforms every step without a mean flow. So it tests the transport, which the stationary droplet cannot do ([`config/oscillatingDroplet2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/oscillatingDroplet2D.yaml#L1-L6)). The transport is reinitialisation-free, and this case shows what that costs in a coupled flow ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2679-L2684)). A redistancing partner or a gradient-control source must first show a gain on this case ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2752-L2756), [gcls article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L880-L882), [[concepts/gradient-control-overview]]). It is also the moving-interface gate of the trace velocity, with period and decay against Lamb as the discriminator ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1195-L1201), [[concepts/trace-velocity-projected-flux]]).

## Where in the code

- The cases: `cases/oscillatingDroplet2D` and `cases/oscillatingDroplet3D`, with [`oscillatingDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet2D.parameter#L1-L52) and [`oscillatingDroplet3D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/oscillatingDroplet3D.parameter#L1-L21). The header of the 2D file still carries the stationary-droplet comment (lines 1 to 4).
- The mode-2 columns: [`writeDropletMetrics.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L59-L116) (the crossings) and lines 156 to 183 (the fit). In 3D the angle uses the x and y components over the 3D radius ([`writeDropletMetrics.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L95-L100)).
- The gate arms: [`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L155-L166) and [`methodGate3D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate3D.yaml#L112-L123).
- The drift figure: [`make_result_figures.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_result_figures.py#L130-L144) writes [`gcls_oscillating_drift.pdf`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/data/figures/gcls_oscillating_drift.pdf) from `gate/histories_oscillating.csv` of the archive `shared-method-config-2026-09-01-192-g1150e68` ([[concepts/data-archive-per-version]]).

## Evidence

The gate numbers are the fixed 2D gate, commit 1150e68, 4 ranks, N = 100 / 142 / 200.

| claim | number | where |
|---|---|---|
| the band gradient error at T | 0.1035 / 0.1403 / 5.129; least-squares order -5.61 | [`methodGate2D_baseline_orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_orders.csv#L3), MEASURED |
| the drift over time | N = 142: 0.036 at t = 0.01 s to 0.140 at 0.10 s, about linear; N = 200: 0.054, 0.087, 0.159, 0.386, 0.581, 1.24, 5.13 at t = 0.01, 0.03, 0.05, 0.07, 0.08, 0.09, 0.10 s | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3976-L3982), MEASURED |
| the period against Lamb | 10.01 / 9.82 / 8.06 ms against 9.51 ms (+5.2 %, +3.3 %, -15.2 %); Celik: divergent | [`methodGate2D_baseline_summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_summary.csv#L15-L17), [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_orders.csv#L6), MEASURED |
| the damping rate | 8.11 / 4.98 / 21.28 1/s; Celik: divergent | [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_orders.csv#L7), MEASURED |
| the velocity jumps at T on the finest rung | L2 norm of U (the oscillation) 1.22e-3 / 6.74e-3 / 8.06e-2 m/s | [`methodGate2D_baseline_summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_summary.csv#L15-L17), MEASURED |
| the volume error has no clean order | 7.64e-3 / 7.31e-3 / 3.04e-3; pairwise 0.12 and 2.56, least squares 1.32 | [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_orders.csv#L4), MEASURED |
| the t = 0 gradient floor converges | 1.73e-3 / 8.34e-4 / 4.12e-4, order 2.07 | [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/methodGate2D_baseline_orders.csv#L8), MEASURED |
| the drift was there before the seam fixes | N = 200 band gradient error 3.41 before, 5.13 after (+51 %); volume 8.18e-3 to 3.04e-3; N = 100 velocity 2.41e-3 to 1.22e-3 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3905-L3919), MEASURED |
| a candidate completes with destroyed fields | HL0 band gradient error 7.5 / 2.2e6 / 2.0e7, volume error 0.74 / 10.4 / 6.7 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3994-L4002), MEASURED |
| the 3D case starts correctly | initial volume error -3.23 % at N = 26 (sphere -3.22 %); initial mode-2 coefficient 7.4e-5 m (sphere 2e-8 m) | [STATUS 11.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3310-L3314), MEASURED (smoke, 20 steps) |

**Results that rest on the algebraic psi.** Every `oscillatingDroplet2D` study before 2026-09-26 ran `implicitEllipsoid` ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3291-L3294)). The record holds these, all MEASURED on that setup:

| claim | number | where |
|---|---|---|
| cellCentreInverse against none, `oscillatingLadder2Dshared` | cellCentreInverse completed N = 128 where none failed at 0.0982 s; none had the lower volume error at N = 32 and 64 | [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393) |
| rhoLENT against geometricFaceDensity | rhoLENT completes N = 128 where the other reaches 98.2 % of the horizon; shape 1.26e-5 against 2.51e-4 | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L651-L653) |
| projectedFlux against cellCentred | projectedFlux reaches the horizon at all three rungs, cellCentred at neither of the two finest | [`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1039-L1040) |
| the foot-evaluated face delivery | blow-up time ratio 0.36 / 0.42 / 0.32 against production at N = 64 / 128 / 256; volume error 5 to 15 times worse | [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L439-L453) |
| the arrival form of the foot | 2 to 4 % of the per-step displacement early, 35 to 47 % by t = 0.02 s | [METHOD 2.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L92-L96), [[concepts/departure-foot-ab2-centring]] |

## Why it failed, or why we think so

The transport keeps no signed-distance property, so the gradient magnitude `q` of psi drifts from 1 over a long run ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2679-L2684)). At N = 100 and 142 the drift stays near 0.1 over ten periods; at N = 200 it turns exponential after t = 0.05 s. The record does not isolate why only the finest rung turns unstable; no run has measured that mechanism (OPEN). STATUS names the drift as the defect that the gradient-control candidates exist to remove. The candidates fail for other reasons, first the shear transport ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3984-L3986), [[concepts/why-the-candidates-failed]]). For HL0 the gcls pre-print reads its growth of about 480 1/s as a velocity fit across the interface jump. There the air boundary layer gives a normal strain of about 250 1/s (HYPOTHESIS, [gcls article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L847-L856)).

## Decisions

- The arm is scored by the period and the damping rate, not by the velocity norms. CORRECTED 2026-09-27: before, a candidate that damped the oscillation more scored as better ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3890-L3894)).
- The gates pin `DROPLET_SURFACE signedDistanceEllipse`; the token default stays algebraic for bit-identity ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3294-L3296)).
- The arm runs `cellCentreInverse`; no valid measurement decides the extension for this case ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393), [[decisions/curvature-extension-cell-centre-inverse]]).

## Open questions

1. Void the algebraic-psi studies (rename `_VOID_algebraicPsi_<date>`) or keep them with the caveat: author decision ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3296-L3297), [plan-halo-limited 9](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-halo-limited-gradient-control.md#L701-L702), [[concepts/wrong-setup-voids]]).
2. The horizon of the arm, and a baseline check that reads the drift over time and not only at T: author decision ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3986-L3988)).
3. The mechanism of the exponential drift at N = 200.
4. A new study on this case that does not set `DROPLET_SURFACE` still runs the algebraic psi ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1033)).
5. The record has no analytic reference for the damping rate; Lamb's inviscid period is a check only ([plan-halo-limited 5.4](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-halo-limited-gradient-control.md#L539-L566)). The effect of the foot mislocation on the period is not measured ([[concepts/departure-foot-ab2-centring]]).
6. The 3D arm never ran, because no candidate passed in 2D ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L4040)). Whether the 3D mode-2 columns measure the spheroid's mode is checked only by the 20-step smoke ([plan-halo-limited B6](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-halo-limited-gradient-control.md#L292-L295)).
7. The arm's curvature delivery, `cellCentreInverse`, is second order on constant curvature only. Its header says that this can be true for the stationary droplet but not for the oscillating one. Nobody scored it on the varying-curvature ellipse gate ([`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L54-L60), [[cases/curvature-static-gates]]). UPDATED 2026-09-29: the face curvature of `cellCentreInverse` is now scored on the signed-distance ellipse and ellipsoid and is second order on both (1.98 and 2.00 over N = 128 to 512; 2.10 over N = 50 to 128) ([METHOD 4.1](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L169-L193), [[concepts/cell-centre-inverse-curvature]]). The gate's oscillating arm sets `DROPLET_SURFACE signedDistanceEllipse` ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/config/gates/methodGate2D.yaml#L165)), the psi on which the static gate is second order. The earlier oscillating studies ran the default algebraic psi (`implicitEllipsoid`), on which every delivery is first order (cCI 0.91 on the implicit-psi ellipsoid). The remainder term is still measured on constant curvature only.

## Related

[[hubs/surface-tension]], [[hubs/gradient-control]], [[hubs/verification]], [[cases/stationary-droplet]], [[cases/translating-droplet]], [[cases/curvature-static-gates]], [[cases/benchmark-cases]], [[concepts/method-gates]], [[concepts/error-vector-and-read-out-instants]], [[concepts/richardson-ladders-and-orders]], [[concepts/wrong-setup-voids]], [[concepts/departure-foot-ab2-centring]], [[concepts/cell-centre-inverse-curvature]], [[concepts/face-curvature-deliveries]], [[concepts/trace-velocity-projected-flux]], [[concepts/halo-limited-extension]], [[concepts/gradient-control-overview]], [[concepts/why-the-candidates-failed]], [[concepts/data-archive-per-version]], [[models/phase-indicator]], [[models/redistancer]], [[decisions/phase-indicator-detrixhe-aslam]], [[decisions/curvature-extension-cell-centre-inverse]], [[retractions/gcls-coupled-loop-reading]], [[studies/gcls-pre-print]], [[studies/sl-quadratic-pre-print]], [[studies/method-gate-2d-campaign-2026-09]]. Pre-prints: [SL](https://leia-openfoam.github.io/leia/preprints/semiLagrangianLevelSet.pdf), [gcls](https://leia-openfoam.github.io/leia/preprints/gclsLevelSet.pdf).

## Log

### 2026-09-29
Created from STATUS 11.4, 11.5 and 11.15, the SL article's Limitations section, the gcls article's results, the 2D gate tables and the case files.
UPDATED item 7 (the static ellipse and ellipsoid scoring of cellCentreInverse, METHOD 4.1).
