---
title: "iDEC defect correction diverges"
description: "The gradient-injecting defectCorrectedIDW reconstruction and its two iterative variants are unstable under repeated semi-Lagrangian transport: the band gradient defect grows 39 to 1.4e4 to 4.7e9 under refinement because the operator has spectral radius above one, and iterating an expansive operator cannot stabilise it (report of 2026-07-10)."
aliases: [iDEC, defectCorrectedIDW, deferred correction, deferredCorrection]
kind: concept
status: retracted
part: advection
tags: [concept, part/advection]
date: 2026-09-28
date_settled: 2026-07-10
decided_by: [config/idwCompare2D.yaml, config/idwCompare3Dshear.yaml, config/idwCompare3Ddeformation.yaml]
code: [src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.H, src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.C, src/leiaLevelSet/semiLagrangian/deferredCorrector.H, src/leiaLevelSet/semiLagrangian/slCorrector.H]
sources: ["iDEC failure report (supplementary, L119-L237)", "METHOD 2.2 (L118-L121)", "METHOD 8.1 row SL_CORRECTION (L374)", "DP L33-L36"]
---
# iDEC defect correction diverges

> Verdict (2026-09-28, report dated 2026-07-10). A departure-foot reconstruction that injects a separately computed cell gradient into the interpolation is unstable under repeated semi-Lagrangian transport. The single-pass `defectCorrectedIDW` drives the band gradient defect on the reversed 2D vortex from 39 to 1.4e4 to 4.7e9 at N = 32, 64, 128, with shape order -0.25 ([iDEC report L125-L136](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L125-L136)). The iterative nodal-defect variant lowers the constant by one to two orders and still grows under refinement ([L141-L155](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L141-L155)). The advection-level deferred correction diverges faster: `nDefCorr` 1, 3, 5, 10 gives 1.66e3, 4.33e12, 2.93e22, 1.25e47 ([L174-L191](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L174-L191)). The diagnosis is that the reconstruction operator has spectral radius above one; a defect correction stabilises only a contracting base operator ([L211-L237](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L211-L237)). The value-fitting `quadraticWeightedLeastSquares` is the stable control on the same flows, shape order 3.01 in 2D ([L193-L209](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L193-L209)). The member stays in the family as a documented failure; production uses the value fit with the `direct` corrector ([METHOD 8.1 L373-L374](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373-L374)).

## What it is

The single pass is an inverse-distance-weighted (Shepard) blend of per-node Taylor values,

    psi(x_d) = sum_C w_C [ psi_C + g_C . (x_d - x_C) ] / sum_C w_C,   w_C = abs(x_d - x_C)^(-p),

with `g_C` the cell-centred least-squares gradient of `psi^n` (`pointCellsLeastSquares`) ([iDEC report L73-L86](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L73-L86)). As a static interpolant it is linear-exact to 3e-16 ([L84-L85](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L84-L85)). Two iterative variants exist:

| variant | what is iterated | controls | code |
|---|---|---|---|
| A, nodal-defect correction | effective nodal values `u_C = psi_C + delta_C` so that a regularised row-stochastic blend reproduces the stencil data; damped Gauss-Seidel, arrival defect pinned to 0 | `dcIters 3`, `dcRelax 0.8`, `dcRegEps 0.1` | [defectCorrectedIDWReconstruction.H L43-L85](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.H#L43-L85) |
| B, advection-level deferred correction | the SL map itself over fixed feet, `psi^{(k+1)}(x_c) = R[psi^{(k)}](x_d)`, gradient recomputed from the current iterate | `nDefCorr`, `defCorrRelax`, `defCorrTol` | [deferredCorrector.H L29-L52](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/deferredCorrector.H#L29-L52) |

The dictionary words are `reconstruction defectCorrectedIDW` ([TypeName L137](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.H#L137)) and `correction deferredCorrection` ([deferredCorrector.H L41-L42](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/deferredCorrector.H#L41-L42)); the token `SL_CORRECTION direct` is the default ([DP L33-L36](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L33-L36)).

## Why it matters

The failure fixes one of the two load-bearing properties of the production reconstruction: the fit must interpolate stencil values, because gradient-injecting variants amplify grid-scale error and diverge ([METHOD 2.2 L118-L121](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L118-L121)). The same family holds `linearTaylor` and `quadraticTaylor`, which also inject a least-squares gradient ([iDEC report L220-L221](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L220-L221), [[concepts/linear-semi-lagrangian]]). It also settles a process point: more iterations of an expansive operator make the divergence faster, not slower ([L232-L234](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L232-L234)).

## Where in the code

- `src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.{H,C}`: the blend matrix is geometry-only and rebuilt into a reused buffer each step to bound 3D memory ([H L78-L79](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.H#L78-L79)).
- `deferredCorrector.{H,C}`: retained as the seam for a contracting deferred correction; its header records the divergence ([H L44-L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/deferredCorrector.H#L44-L49)).
- Studies: `config/idwCompare2D.yaml` (N = 32, 64, 128; CFL 0.5 and 1.0; np 4) ([L1-L26](https://github.com/leia-openfoam/leia/blob/8867581/config/idwCompare2D.yaml#L1-L26)), `idwCompare3Dshear`, `idwCompare3Ddeformation` and their `128` variants.

## Evidence

| claim | number | where |
|---|---|---|
| Single pass, 2D reversed vortex, CFL 0.5, N = 32 to 128: band gradient defect and shape error | 39 to 1.4e4 to 4.7e9; shape 1.6e-2 to 3.0e-2, order -0.25 | MEASURED, [iDEC report L125-L136](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L125-L136) |
| Single pass, 3D shear and 3D deformation | gradient 8.5e3 to 1.2e9 and 8.8e4 to 4.5e12; shape 0.83 to 0.92 and 0.45 to 0.48 (interface destroyed) | MEASURED, same table |
| Variant A (nodal defect), 2D, N = 32 to 128 | gradient 15.2 to 1.66e3 to 6.3e7; shape order -0.33 | MEASURED, [L141-L155](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L141-L155) |
| Variant A is converged, not under-iterated (2D, N = 64) | `dcIters` 3: 1657; 30: 1649; `limitSlope true`: 1657 | MEASURED, [L157-L172](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L157-L172) |
| Variant B (deferred correction), 2D, N = 64: gradient defect at T against `nDefCorr` | 1: 1.66e3; 3: 4.33e12; 5: 2.93e22; 10: 1.25e47; 10 with relax 0.5: 5.62e18 | MEASURED, [L174-L191](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L174-L191) |
| Control: value fit `quadraticWeightedLeastSquares`, same flows | 2D shape 3.0e-4 to 4.7e-6, order 3.01; 3D shear order 2.43; 3D deformation 1.34; gradient bounded at 0.1 to 0.5 | MEASURED, [L193-L209](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L193-L209) |
| The static interpolant is linear-exact | error about 3e-16 | MEASURED, [L84-L85](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L84-L85) |
| The operator is (bounded smoother) composed with (unbounded differentiator), so `rho(M) > 1` on grid-scale modes | qualitative | DERIVED, [L211-L221](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L211-L221) |

The gradient columns are L2 norms of `abs(grad psi) - 1` in the narrow band ([L66-L69](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L66-L69)).

## Why it failed, or why we think so

The lift term `g_C . (x_d - x_C)` feeds the gradient of `psi^n` back into `psi^{n+1}` without ever being constrained by the reconstructed values. Differentiation amplifies grid-scale content and the Shepard smoother does not remove it, so `abs(grad psi)` diverges geometrically in the step index ([L211-L221](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L211-L221)). Variant A drives the nodal defect to zero but still evaluates the singular-weight blend of `psi_C + g_C . r`: the gradient path is not severed. Variant B is Picard iteration on a divergent fixed point and amplifies faster ([L223-L237](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L223-L237)). The always-on value cap keeps the values bounded, so the shape error saturates instead of producing NaN ([L137-L139](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L137-L139)).

## Decisions

- `SL_RECONSTRUCTION` stays a value fit ([[decisions/sl-reconstruction-uncached-qwls]]); `SL_CORRECTION direct`, `deferredCorrection` is a research path no study selects ([METHOD 8.1 L374](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L374)).
- Untested directions for a contracting correction are listed in the report: a consistent-gradient lift, a bounded lift, a Khosla-Rubin split, BFECC, a structure-preserving gradient ([L239-L266](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L239-L266)).

## Open questions

1. The studies ran at np 4 with the kinematic solver on 2026-07-10, before the gradU coupled-patch fix of 2026-08-26 ([idwCompare2D.yaml L10-L11](https://github.com/leia-openfoam/leia/blob/8867581/config/idwCompare2D.yaml#L10-L11), [STATUS L619-L636](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L619-L636)). The numbers carry that defect; the divergence by nine orders of magnitude is far outside its measured effect (a factor of at most 2 at the endpoint, [PCS L58-L59](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L58-L59)), so the verdict stands and the digits do not.
2. None of the five directions has been tried.

## Related

[[hubs/advection]] - [[models/sl-reconstruction]] - [[concepts/linear-semi-lagrangian]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/eulerian-fv-transport]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[retractions/gradu-coupled-patch-contamination]] - [[studies/sl-quadratic-pre-print]]

## Log

### 2026-09-28
Created from the supplementary iDEC failure report, the reconstruction and corrector headers and the study configs.
