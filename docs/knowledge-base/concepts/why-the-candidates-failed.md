---
title: "Why the candidates of the first campaign failed"
description: "One row per candidate: how it failed, the mechanism the pre-print named, the mechanism the technical report holds responsible, the mark (MEASURED, DERIVED, HYPOTHESIS) and the experiment that decides it"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [technical report, gcls article sec:discussion, STATUS 11.15]
---
# Why the candidates of the first campaign failed

> Open (2026-09-28): the readings below are the technical report's; each carries its mark and
> the experiment that decides it. Every number is from the archive
> `shared-method-config-2026-09-01-192-g1150e68` ([[concepts/data-archive-per-version]]) and from
> [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023).

## The candidates and how they failed (MEASURED)

| candidate | what | how it failed |
|---|---|---|
| S1 | soft wall alone, band 3h | shear: shape error 410x (order 0.06), band gradient ratio 0.71 (target met); translating N = 142 diverges; oscillating volume error 21; coupled seam 0.76 FAIL; stationary droplet identical to the baseline to four digits (the wall is off there) |
| HL0 | halo-limited extension alone, R = h | shear: 14x, order 1.39, ratio 0.83 (target missed); oscillating band gradient 2e7, volume 6.7; translating shape 64x; stationary volume 2e-7 against 2e-8; seam 4.7e-6 PASS |
| HL1q | HL0 plus linearQ, mu = 1/T_ref | shear: ratio 0.20 (target met), shape 173x (order 0.02); stationary N = 200 destroyed (shape 3.8 R, volume 76 %, current 0.088 m/s); translating shape 59x; seam 4.9e-7 PASS |
| HL1z | HL0 plus linearZ | shear: band gradient 10.5, shape 1.65, volume 1.6; stationary destroyed (3.7 R, 84 %); seam 5.3e-7 PASS |
| HL2 | HL0 plus soft wall | translating diverges at all three rungs; oscillating volume 14.6; shear 400x; seam NOT_COMPARABLE |
| FP0 | closest-point extension | shear 34x (order 1.38), ratio 0.94; kinematic seam FAIL 0.45; coupled seam FAIL 1.11; stationary and oscillating N = 200 timed out |

## The readings, before and after the report

| failure | the pre-print's reading (2026-09-27) | the report's reading (2026-09-28) | mark | decides it |
|---|---|---|---|---|
| HL1q/HL1z destroy the stationary droplet | the coupled loop psi -> kappa -> f -> u -> psi -> q -> F at rate mu, as in SDPLS | the source's own map is unstable with a centred q: the mode $(-1)^{i+j} d$ grows by $(1 + \mu\Delta t)$ per step; the laws read no velocity; the measured rates equal mu; the flow loop is the second stage | DERIVED (mechanism), MEASURED (rates about 290 and 320 to 350 1/s over 10 to 20 ms, against 270; E1.4 decides the scaling) | E1.1, E1.4 |
| the O(1) shape errors of S1, HL1q, HL1z in the shear flow | the crossing shift $\delta = h\,ab/(a+b)^2\,\Delta t\,(F_A - F_B)$ | the shift is O(h) (at most 0.2h summed); the O(1) errors come from the band edge under the accumulated factor $q_\mathrm{far}/q_\mathrm{target} \approx 2$ to 3 and, for the soft wall, from its 40σ slope | DERIVED | E1.2, E1.3, E1.7, E2.2 |
| S1's 76 % seam failure | round-off of the decomposition carried by $q_c$ into F | the bang-bang map of the soft wall (dt dG/dq ≈ 0.5 per step); HL1q/HL1z pass the same check | DERIVED, MEASURED (seam) | E1.7 |
| HL0 loses the transport order (1.39) | the reconstruct of a flux correction that varies within a cell | the extension relocates the normal strain into a shell at d = 1 to 2h (K = 1.27 at 1.5R); the reconstruct is one O(h grad u) contribution, neither necessary (FP0) nor sufficient (HL0's early volume error) | DERIVED, MEASURED (1D band mean 0.46 to 0.59 against 0.49) | E2.3, E0.1 |
| HL0 destroys the coupled arms | the sampler fits across the viscosity and density jump | translating: the sampler injects the spurious current unprojected (D = 0 there, no jump to fit across); oscillating: the boundary-layer kink reading is plausible | MEASURED (183 and 178 1/s growth from 0.02 s), HYPOTHESIS (oscillating) | E1.6, E1.5 |
| FP0 | the same reconstruct | no shell, same order loss: the tangential shear of a normal-constant extension and the kinks at the shear tail; the search is decomposition-dependent | DERIVED, MEASURED (seam) | Tier-I rotation |

## What the evidence for the primary reading is (MEASURED unless marked)

1. The curvature error of HL1q at $N = 200$ e-folds at about 290 1/s over 10 to 20 ms (faster,
   about 460 1/s, over 2 to 10 ms) and the spurious current at 320 to 350 1/s, against
   $\mu = 1/T_\mathrm{ref} = 270$ 1/s; the mesh-scale capillary rate ($\sim 10^5$ 1/s) is not
   seen. The linearised eigenmode predicts a rate of exactly $\mu$; the excess is not explained
   (the nonlinear regime, or the flow loop adding to it), and E1.4's scaling with $\mu$ decides.
2. HL1q and HL1z agree to three digits until $t = 0.02$ s: their linearisations coincide.
3. The curvature error is 10 times the baseline's at 4.3 ms while the velocity is 1.4 times: the
   curvature leads the velocity by about 5 ms.
4. At 4.3 ms the centred band-gradient metric is BELOW the baseline's (3.2e-4 against 4.4e-4)
   while the curvature error is 10x; at 10 ms it is 3.9x while the curvature error is 175x (28.3
   against 0.16 1/m): the metric cannot see the mode (DERIVED: a pure checkerboard is invisible to
   a centred difference).
5. The smooth part of the initial q error decays as the law intends: at 2.3 ms the band gradient
   error is 2.2e-4 against the baseline's 4.1e-4.
6. $\min q$ and $\max q$ leave [0.99, 1.01] symmetrically (0.83/1.17 at 20 ms, 0.33/2.29 at 30 ms):
   an oscillatory cell-scale mode, not a smooth drift.
7. The seed estimate: the cell-to-cell part of the initial q error, O($h^2/R^2$) ≈ 6e-4 at
   $N = 200$, reaches O(1) at $\ln(1/6\times10^{-4})/270 = 27$ ms; observed 25 to 35 ms (DERIVED).

## Related

[[hubs/gradient-control]], [[concepts/source-discretisation-defects]],
[[concepts/extension-strain-relocation]], [[concepts/gradient-control-next-experiments]],
[[retractions/gcls-coupled-loop-reading]], [[studies/method-gate-2d-campaign-2026-09]].

## Log

### 2026-09-28
Created with the technical report.
