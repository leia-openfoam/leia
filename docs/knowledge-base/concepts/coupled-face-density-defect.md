---
title: "The coupled-face density defect (fixed 28d13f0)"
description: "Until 2026-09-27 each rank built the face density of a processor face from its own cell, so rho_f differed by up to 90 percent across a seam, and the droplet metrics counted internal faces only; both are fixed, the gate checks them, and which earlier parallel studies to re-run is open"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-27
decided_by: [commit 28d13f0, commit b1798c3, commit aaae627]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H, config/gates/methodGate2D.yaml]
sources: [STATUS 11.14, STATUS 11.15, METHOD 6, gcls article sec:fv-parallel, RM 510-511, CLAUDE 4-rank section]
---
# The coupled-face density defect (fixed 28d13f0)

> Verdict (2026-09-28). Two defects of the coupled SL solver were found and fixed on 2026-09-27, both of one kind: code that handled a processor face as a boundary face of one cell. First, `computeFaceAreaFractions` filled every boundary face from its local cell, so on a processor face each rank used its own plane whatever the flux direction; after one step at N = 64 on 4 ranks, 4 processor faces had `alpha_f` different on the two sides and `rho_f` up to 90 % apart while `phi` was equal; mass was not conserved across the seam where the interface crossed it (28d13f0, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3630-L3641)). Second, the droplet metrics sampled the zero set, the band and the face statistics on internal faces only, so serial and 4-rank runs differed by 3 to 8 % in shape, band and curvature on fields equal to 1e-8 (b1798c3, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3653-L3667)). After both fixes the four-rank translating droplet agrees with the serial run to 3.5e-7 or better in every error-vector column over 0.05 s ([gcls article, sec:fv-parallel](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L397-L455)). Every SL two-phase result on more than one rank before 2026-09-27 carries both errors; which studies to re-run or void is an open author decision ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3668-L3675)). The roadmap had flagged the owner-only face fraction in July ([roadmap](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L510-L511)).

## What it is

The face area fraction of the mass flux is cut from the reconstructed plane of the donor cell ([[concepts/alphaf-source-donor-plane]]). A coupled face (processor, cyclic) has a donor on each side, and both sides must use the same `alpha_f`, or `rho_f * phi` differs between the two sides of one face. In the corrected code each side computes its local-cell fraction, the values are exchanged with a collective `syncTools::swapBoundaryFaceList`, and the donor rule is applied with the local flux sign: `phi > 0` leaves the local cell (the local cell is the donor), `phi < 0` enters it; at `phi == 0` both sides take the owner side's value ([`faceAreaFraction.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H#L242-L304)). A processor face then gets the value an internal face gets in a serial run, up to the summation order.

The diagnostics: the zero-set crossings on processor faces are counted once, on the owner side of the patch ([`writeDropletMetrics.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L117-L153)); the band threshold and the band cell set use internal and processor faces ([same](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L229-L260)); the face statistics count processor faces once ([same](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L339-L350)). The research diagnostics `A4h*`, `A8h*` and `driver*` still loop internal faces only ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3712-L3715)).

The gate check: `translatingSeamNp1` runs the coarsest translating rung in serial and compares nine error-vector columns of the droplet CSV with the np 4 run at a column-scaled tolerance of 1e-5; a reference that did not complete makes the check NOT_COMPARABLE ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L145-L152), aaae627, [[concepts/seam-checks-and-decomposition-invariance]]).

## Why it matters

The per-rank mass residual of rhoLENT (1e-10 to 1e-13) did not see the defect: each rank's auxiliary density balances its own flux ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3640-L3641)). A residual that is local to a rank is not a conservation check across ranks. The defect had been named by inspection on 2026-08-19 in the seam-consistency configs ("faceAreaFraction.H:201-212 takes the face area fraction from the patch-internal cell on every boundary patch including processor", [`seamConsistency3Dserial.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/seamConsistency3Dserial.yaml#L19-L21)) and in July in the roadmap, and it stayed in the code for five weeks.

## Evidence

| claim | number | where |
|---|---|---|
| the two sides of a processor face disagreed | 4 of 37 processor faces, `rho_f` worst relative difference 9.04e-01 at step 1; 3 faces, 9.68e-01 at step 5; 0 after the fix | [`gcls_face_density.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/gcls_face_density.tex), MEASURED |
| the solution error of the decomposition | np 4 against serial, velocity metrics 1e-5 to 5e-4 by t = 0.033 s before; 1e-8 to 1e-10 through step 5000 after | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3636-L3639), MEASURED |
| the diagnostic bias | `zeroSetRadialL2` 3.42e-02 to 1.20e-07; `gradPsiL2ErrorBand` 3.82e-02 to 1.24e-08; `kErrL2Band` 7.59e-02 to 2.11e-08; `m2Amplitude` 6.31e-01 to 1.92e-07 (column-scaled, to t = 0.05 s) | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3653-L3667), [`gcls_seam.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/gcls_seam.tex), MEASURED |
| serial runs are unchanged | bit-identical over 4604 steps, every column, tolerance 0 | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3665-L3667), MEASURED |
| the late instability is not this defect | divergence at 0.0775 s (serial), 0.0904 s (np 4 before), 0.0868 s (np 4 after) | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3642-L3652), MEASURED |
| the cluster reproduces the laptop | `translatingSeamNp1`: 3.47e-7 (`meanMagUPrime`), 1.2e-7, 2.1e-8, the same digits as the laptop | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3879-L3886), MEASURED |
| what the defects changed in the gate baseline | stationary N = 100 shape +9.4 %; oscillating N = 100 `l2MagUPrime` -50 %; oscillating N = 200 band gradient +51 %, volume -63 %; the kinematic shear arm 0.0 % | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3887-L3916), MEASURED |

## Decisions

- Both fixes shipped with a serial bit-identity gate and the np 4 check; the coupled decomposition check is part of the 2D gate since aaae627 ([[concepts/method-gates]]).
- The fixed gate re-ran the whole campaign with `PRESERVE=1` and was compared with the pre-fix record ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3887-L3916)).

## Open questions

1. Which earlier parallel SL two-phase studies to re-run or void (author decision; [[hubs/mass-flux]], open item 2). The record before 28d13f0 carries a decomposition error of about 1e-4 in the solution and a bias of several percent in the shape, band, curvature and mode-2 columns.
2. The research diagnostics that still loop internal faces only.

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/alphaf-source-donor-plane]], [[concepts/eulerian-solver-mass-flux-port]].
- Method: [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/method-gates]], [[concepts/error-vector-and-read-out-instants]].
- Retractions of the same class: [[retractions/gradu-coupled-patch-contamination]], [[retractions/psi-filter-seam-bug]].
- Case: [[cases/translating-droplet]].

## Log

### 2026-09-28
Created.
