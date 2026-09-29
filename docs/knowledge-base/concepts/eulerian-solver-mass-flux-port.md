---
title: "The Eulerian two-phase solver: frozen rho and the rhoLENT port"
description: "leiaLevelSetTwoPhaseFoam kept rho, rhoPhi and the face viscosity at their t = 0 values until c094bd8 (2026-09-27) and moved the heavy droplet at 42 percent of the stream; it now shares the mass-flux headers of the SL solver; the five frozen-rho studies wait for the void decision"
aliases: []
kind: concept
status: settled
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled: 2026-09-27
decided_by: [commit c094bd8]
code: [applications/solvers/leiaLevelSetTwoPhaseFoam/leiaLevelSetTwoPhaseFoam.C, applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H]
sources: [STATUS 11.14, gcls article sec:fv-parallel, PHL 98, PHL 434-440, PHL 703, c094bd8]
---
# The Eulerian two-phase solver: frozen rho and the rhoLENT port

> Verdict (2026-09-28). Until c094bd8 (2026-09-27) the Eulerian two-phase solver `leiaLevelSetTwoPhaseFoam` kept `rho`, `rhoPhi` and the face viscosity `muf` at their `t = 0` values for the whole run; `rho` was reassigned only on a dynamic mesh refinement ([`leiaLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/leiaLevelSetTwoPhaseFoam.C#L103-L109)). With the frozen density the heavy translating droplet travelled at 42 % of the stream (travelled fraction 0.420, `L2|U-U0|` 0.155 m/s); with the shared rhoLENT flux it travels at the stream speed to 1e-4 (1.00012, 4.85e-4 m/s) after 1000 steps at N = 100 ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3681-L3711), [gcls article, sec:fv-parallel](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L456-L476)). Every Eulerian two-phase result with a moving interface and a density contrast before that commit is wrong physics. Whether the five frozen-rho studies are renamed `_VOID_frozenRho_<date>` is an open author decision (plan item D-m, [PHL](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L434-L440)).

## What it is

The SL two-phase solver's mass-flux code was moved verbatim into three headers that both solvers include: the selection and switches ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L1-L13)), the face density from the new interface ([`updateFaceDensity.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateFaceDensity.H#L1-L9)) and the mass flux with the auxiliary equation on every outer corrector ([`updateMassFlux.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/updateMassFlux.H#L1-L8)). The SL files reconstruct byte for byte from them ([c094bd8](https://github.com/leia-openfoam/leia/commit/c094bd8)). The Eulerian solver now reads the same `levelSet.massFlux` dictionary with the same defaults, rebuilds `alphaf`, `rhof` and `muf` from the new interface, solves the rhoLENT density on every outer corrector, resets `rho` after the PIMPLE loop and prints the density-bound numbers to its log ([`alphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H#L88-L127)). `PIMPLE.psiOuterCorrectors yes` rebuilds the interface on every outer corrector from the current psi iterate ([`leiaLevelSetTwoPhaseFoam.C`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/leiaLevelSetTwoPhaseFoam.C#L111-L127)).

Two steps remain before the Eulerian line can join the coupled gate arms ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3702-L3711)):

1. The curvature pipeline. The Eulerian production curvature is the older `meanCurvatureClosestPoint` plus the `stabilizedFootPointFace` delivery; it does not read `curvatureExtension`, so it cannot run `cellCentreInverse`.
2. The droplet CSV (plan F2). `writeDropletMetrics.H` has the SL CSV name hard-coded and reads one SL fit object; everything else it needs exists in the Eulerian solver.

## Why it matters

The Eulerian line is the second solver line of the method gates ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L34-L36)); its candidate file still says the density is frozen and only its kinematic arms are meaningful ([`baselineEulerian.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/candidates/baselineEulerian.yaml#L1-L9), stale at 8867581). The SDPLS coupled results of the pre-print of that line were produced with the frozen density ([[studies/sdpls-pre-print]]).

## Evidence

| claim | number | where |
|---|---|---|
| the frozen density moves the heavy droplet at 42 % of the stream | travelled fraction 0.420, `L2|U-U0|` 1.55e-1 m/s, `L1` 3.05e-2 m/s at ratio 840; 1.00339 and 2.18e-4 at ratio 1 | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3690-L3696), [`gcls_eulerian.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/data/tables/gcls_eulerian.tex), MEASURED |
| the shared flux restores the stream | 1.00012, 4.85e-4 m/s, 2.54e-4 m/s at ratio 840 (np 4 and np 1); 1.00078, 1.04e-4 at ratio 1; relative mass residual 1.2e-10 | same, MEASURED |
| the port is decomposition-invariant | np 4 against serial: `U`, `alpha`, `psi`, `p_rgh` equal to 1e-8 relative | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3697-L3699), MEASURED |
| the stationary droplet with the port | displacement 2.5e-13 m, `L2|U|` 8.0e-5 m/s at np 4 | same, MEASURED |
| the extraction did not change the SL solver | bit-identical in serial; bit-identical on 4 ranks over 300 steps (`compare_metrics_csv --tol 0`) | [c094bd8](https://github.com/leia-openfoam/leia/commit/c094bd8), MEASURED |

## Decisions

- The port itself (c094bd8) is not inert by design and was gated on the translating droplet at both density ratios, on the stationary droplet and on 4 ranks ([PHL](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L434-L440)).

## Open questions

1. D-m: void `sdplsDropletNS2D`, `sdplsDropletMechanism2D`, `sdplsDropletBdf2Droplet2D`, `sdplsPsiBudgetDroplet2D` and `sdplsRdivDroplet2D`, or keep them with a caveat (author decision, [PHL](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L703)).
2. The curvature-extension dispatch and the droplet CSV of the Eulerian solver (F2).
3. The face viscosity `muf` was frozen in the same way; the viscosity decision was measured on the SL solver ([[models/viscosity-face-model]]).

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]], [[models/level-set-advection]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/bound-rho]], [[concepts/coupled-face-density-defect]], [[concepts/wrong-setup-voids]].
- Lines and studies: [[concepts/sdpls-source-eulerian]], [[studies/sdpls-pre-print]], [[concepts/method-gates]].

## Log

### 2026-09-28
Created.
