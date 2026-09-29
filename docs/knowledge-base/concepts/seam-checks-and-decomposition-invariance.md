---
title: "Seam checks: code that is right on one rank and wrong on several"
description: "Eight defects of one class, from setVelocity (2026-08-26) to the face density and the droplet metrics (2026-09-27); the 4-rank gate before the cluster, the serial-against-np-4 configs and the gate's seam arms are the standing checks"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-28
date_settled: 2026-08-28
decided_by: [author decision 2026-08-28, config/seamConsistency3Dserial.yaml, config/seamConsistency3Dpar4.yaml, config/gates/methodGate2D.yaml]
code: [src/leiaLevelSet/velocityModel/velocityModel.C, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H, workflow/scripts/make_gate_summary.py, config/gates/methodGate2D.yaml]
sources: [CLAUDE 4-rank section, STATUS 4 gradU sections, STATUS 4 psi-filter section, STATUS 11.14, STATUS 11.15, docs/gradU-coupled-patch-contamination.md, gcls article sec:fv-parallel]
---
# Seam checks: code that is right on one rank and wrong on several

> Verdict (2026-09-28). No new or changed algorithm goes to the cluster until it has run on at least four MPI ranks on the laptop and the result has been read; a serial pass says nothing about processor-patch values, halo exchange and collective calls ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L199-L246)). The rule is written from one class of defect that this repository hit eight times: `setVelocity` wrote face values into coupled patches and biased `fvc::grad(U)` by O(1) in every processor-adjacent cell, which contaminated 31 parallel kinematic SL studies ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L619-L640), [[retractions/gradu-coupled-patch-contamination]]); the psi filter's `fvc::average` carried `calculated` patches, so `L(psi)` was uncoupled across seams and the filtered cell set depended on the decomposition ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L923), [[retractions/psi-filter-seam-bug]]); the semi-implicit capillary force called collective reductions inside a `Pstream::master()` guard and left four cluster arms alive but silent for 76 minutes; the face density of a processor face came from each rank's own cell, and the droplet metrics counted internal faces only ([[concepts/coupled-face-density-defect]]). A decomposition-dependent answer is a defect even when both runs complete. The standing checks are the 4-rank gate, the serial-against-np-4 configs, and the seam arms of the method gates ([[concepts/method-gates]]).

## What it is

The class ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L217-L228)):

| defect | what was wrong | fixed | effect |
|---|---|---|---|
| `velocityModel::setVelocity` | a face-centre value written into a coupled patch, where every operator expects the neighbour-cell value; `fvc::grad(U)` off by `a h/4` per face, O(a) after the Gauss sum | 30e6ba9, 2026-08-26 | the SL Taylor foot consumes `grad(U)`; 31 parallel kinematic SL studies contaminated, 10 with curated tables ([gradU note](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L18-L60)) |
| `interfaceExtension::updateFlux` | the raw flux assigned over the whole boundary field | before 2026-08-28 | listed in the class |
| the psi filter | `Lpsi = psi - fvc::average(...)` inherited `calculated` patches; the band dilation looped internal faces only | f83a1ab, 2026-08-19 | 53 to 205x worse than the unfiltered control before the fix; every filtered ladder ran with it ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L923)) |
| a narrow-band dilation | internal faces only; the filtered cell set decomposition-dependent | same | same |
| `semiImplicitCapillaryForce` | `gAverage`, `gMax`, `gSum` inside a `Pstream::master()` guard; rank 0 blocked in the reduction | 2026-08-28 | four cluster arms alive and silent for 76 minutes; invisible in every serial test ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L207-L215), [[models/semi-implicit-capillary-force]]) |
| `computeFaceAreaFractions` | a processor face filled from each rank's own cell; `rho_f` up to 90 % apart across a seam | 28d13f0, 2026-09-27 | np 4 against serial 1e-5 to 5e-4 in the velocity metrics before, 1e-8 to 1e-10 after ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3630-L3641)) |
| the droplet metrics | the zero set, the band and the face statistics on internal faces only | b1798c3, 2026-09-27 | shape, band and curvature 3 to 8 % apart between serial and np 4 on fields equal to 1e-8 ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3653-L3667)) |
| the closest-point extension | the closest-point search trusts halo data and falls back to a steady solve; decomposition-dependent | open | seam FAIL 1.4e-2 in the kinematic smoke, 0.45 in the gate, 1.11 in the coupled check ([STATUS 11.11](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3479-L3481), [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L4015-L4017)) |

The Eulerian SDPLS solver had two coupled-patch defects of the same kind, a face-centre value and a raw flux ([[studies/sdpls-pre-print]]).

The gate, on the laptop before any `sbatch` ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L230-L246)):

1. `decomposePar -force` with at least four subdomains; a serial case decomposes to one domain and proves nothing.
2. `mpirun -np 4 <solver> -parallel`.
3. Read the step count and the metrics CSV, not the exit code: a deadlock leaves the job alive with a truncated log and no nonzero return code anywhere ([[concepts/log-classifier-and-waiters]]).
4. Where the study is itself parallel, run the serial-against-np-4 equivalence check: `config/seamConsistency3D{serial,par4}.yaml` is the committed pattern, 20 steps, the identical case with the filter on and off, the off arm as the control that separates a filter seam bug from any other seam defect ([`seamConsistency3Dserial.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/seamConsistency3Dserial.yaml#L1-L58)).
5. Anything that touches an `fvMatrix` (an fvOption, a new term, a new solver) runs the parallel smoke as part of its inertness gate, not as an extra.

In the method gates the coarsest shear rung runs again at np 1 and np 8 next to np 4 and every column must agree to 1e-10 (column-scaled maximum over the rows); the coarsest translating rung runs again in serial and nine error-vector columns must agree to 1e-5; a reference that did not complete is NOT_COMPARABLE, because a divergence decorrelates the runs ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L110), [same](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L145-L152), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L292-L340)).

## Why it matters

A coupled seam tolerance has a floor. The pressure solve converges to a relative tolerance and GAMG's agglomeration depends on the decomposition, so np 1 and np 4 on the same refined 3D mesh agreed to 1.7e-4 in the mean current, 5.8e-5 in its L2 norm and 1.6e-8 in the Laplace jump, and a pre-registered line of 1e-10 was unattainable by construction; the yardstick is the same pair on a uniform mesh ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1220-L1246)). A seam difference far above that floor is a property of the candidate: the soft-wall source S1 differed by 76 % in `l2MagUPrime` between serial and np 4 where the baseline differed by 3.5e-7, which is the mesh-noise-floor diagnostic of the regression set applied to the decomposition ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3992-L4001), [[concepts/advection-regression-set]]).

## Evidence

| claim | number | where |
|---|---|---|
| the gradU defect grew with the rank count | 1D stretch, source R, band mean `d(psi)/dx` against the exact 1: np 2 1.004681, np 4 1.005367, np 8 0.956834; serial 1.000000; after the fix 1.000000 at every count | [gradU note](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L45-L60), MEASURED |
| its effect on a curated endpoint | 2D vortex N = 256: shape 1.105e-06 (serial, bug-free) against 7.851e-07 (curated, buggy parallel); volume 3.95e-05 against 8.01e-05; no order change | [STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L580-L617), MEASURED |
| the psi-filter seam bug | 3D N_L = 60, 20 steps, serial against np 4, filter on: max U 2.43e-04 to 1.72e-06 (control off 1.55e-06); A2h band 8.38e-07 to 3.72e-10; volume 1.92e-05 to 6.55e-07 | [STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L892-L923), MEASURED |
| the coupled-face fix | `rho_f` up to 90 % apart on 4 processor faces before, 0 after; np 4 against serial 1e-5 to 5e-4 before, 1e-8 to 1e-10 after | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3630-L3641), MEASURED |
| the metric fix | `zeroSetRadialL2` 3.42e-02 to 1.20e-07, `kErrL2Band` 7.59e-02 to 2.11e-08 (column-scaled, serial against np 4) | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3653-L3667), MEASURED |
| the coupled seam check of the campaign | baseline 3.5e-7, HL0 4.7e-6, HL1q 4.9e-7, HL1z 5.3e-7 PASS; S1 0.76, FP0 1.11 FAIL; HL2 NOT_COMPARABLE | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L4015-L4017), MEASURED |
| the kinematic seam check of the campaign | np 1 and np 8 against np 4 at most 1.9e-12 (baseline smoke); FP0 1.4e-2 | [STATUS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3281-L3306), [STATUS 11.11](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3479-L3481), MEASURED |
| a collective in a guard deadlocks silently | 4 arms alive and silent for 76 minutes; the 4-rank gate of `boundRho` and the residual diagnostic exists because both contain reductions | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L207-L215), [`rhoBoundGate2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoBoundGate2D.yaml#L1-L34), MEASURED |

## Decisions

- The 4-rank gate before the cluster (2026-08-28); the seam arms of the gates (2026-09-26, coupled check 2026-09-27, aaae627).
- A per-rank residual is not a conservation check across ranks: the rhoLENT residual of 1e-10 to 1e-13 did not see the face-density defect ([[concepts/coupled-face-density-defect]]).

## Open questions

1. Which parallel SL two-phase studies before 2026-09-27 to re-run ([[hubs/mass-flux]]).
2. The closest-point extension's decomposition dependence ([[models/velocity-extension]]).
3. The one-sided fallback of the face-curvature delivery at seams is named in the curvature plan and has no measurement of its own ([plan](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1453)).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/method-gates]], [[concepts/bit-identity-and-inertness-gates]], [[concepts/advection-regression-set]], [[concepts/log-classifier-and-waiters]], [[concepts/wrong-setup-voids]].
- Defects and retractions: [[concepts/coupled-face-density-defect]], [[retractions/gradu-coupled-patch-contamination]], [[retractions/psi-filter-seam-bug]], [[models/semi-implicit-capillary-force]], [[models/narrow-band]].
- Cases: [[cases/exact-1d-stretch]] (the 1D gate that found the gradU defect).

## Log

### 2026-09-28
Created.
