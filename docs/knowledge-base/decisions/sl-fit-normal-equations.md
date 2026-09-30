---
title: "SL_FIT: normal equations, because QR blows up identically"
description: "SL_FIT normalEquations, decided by config/popinet3D_La12000_poly_dump4_qr.yaml on 2026-09-05: householderQR diverges on the same step with the same phase volume to six figures, and the amplifying polyhedral cells are well conditioned (pivot 0.757)"
aliases: [SL_FIT normalEquations]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-09-05
decided_by: [config/popinet3D_La12000_poly_dump4_qr.yaml]
code: [src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C, cases/default.parameter]
sources: ["METHOD 8.1 row SL_FIT (L377)", "STATUS 4 (L2269-L2286)", "STATUS 4 (L1566-L1575)", "MC article L275-L280", "DP L263-L264"]
---
# SL_FIT: normal equations, because QR blows up identically

> `SL_FIT normalEquations` in the global default since 0f02aea (2026-07-30) ([DP L263-L264](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L263-L264)). Decided by `config/popinet3D_La12000_poly_dump4_qr.yaml` on 2026-09-05 ([METHOD 8.1 L377](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L377)): with `householderQR` the polyhedral Popinet 3D case blows up identically, with a step-3 phase volume of 0.017512193 against 0.017512208 for the normal equations, about six significant figures. The pivot census of 2026-09-09 explains it: the worst amplifying cell has `Lambda = 1.2608` and a scaled pivot of 0.757, six times the admissibility tolerance 0.3, so better arithmetic has nothing to repair ([STATUS L2269-L2278](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2269-L2278)). On hexahedra and on a 10 %-perturbed 2D mesh the two solvers are at bit parity ([MC article L275-L280](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L275-L280)).

## The question

The quadratic fit is solved as the normal system `M a = r`, `M = A^T W^2 A`, by in-place Cholesky ([METHOD 2.2 L108-L111](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L108-L111)). The reconstruction warns of near-degenerate stencils on cfMesh's near-wall polyhedra, with condition numbers of 3e7 to 6e12 and coefficients of 3e2 to 3e8 ([STATUS L1566-L1575](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1566-L1575)). The hypothesis of 2026-09-05: the normal equations square the condition number, so a Householder QR of the weighted design rows removes the polyhedral divergence ([config header L11-L16](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet3D_La12000_poly_dump4_qr.yaml#L11-L16)). The key is `fit` in `levelSet.semiLagrangian`, with the values `normalEquations` and `householderQR` ([uncached C L187-L213](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L187-L213), [L655-L692](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L655-L692)).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| Popinet 3D, pMesh, N = 64, np 4, `householderQR` | phase volume at step 3 | 0.017512193 against 0.017512208 with `normalEquations`: the same blow-up | MEASURED, [METHOD 8.1 L377](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L377) |
| the same mesh, one fit pass | pivot of the worst amplifier (cell 14736, `Lambda` 1.2608, size 0.333 h) | 0.757 (tolerance 0.3) | MEASURED, [STATUS L2269-L2278](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2269-L2278) |
| the same mesh | the top 0.1 % of `Lambda` (674 cells) | all still quadratic; the 47 622 cells the pivot test demotes reach `Lambda` 1.0511 only | MEASURED, [STATUS L2269-L2278](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2269-L2278) |
| 2D hex and 10 %-perturbed mesh | error against the Cholesky solve | identical to the last digit, zero degenerate fallbacks | MEASURED, [MC article L275-L280](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L275-L280) |

Pre-registered read-out ([config header L14-L16](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet3D_La12000_poly_dump4_qr.yaml#L14-L16)): PASS = phase volume constant to 1e-7 over four steps and Courant about 0.1; a blow-up as before exonerates the fit's conditioning. The blow-up occurred.

## What it does not cover

1. Conditioning is handled by the admissibility test `quadraticPivotTol 0.3`, which demotes a cell to the linear fit ([STATUS L1657-L1681](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1657-L1681), [[models/sl-reconstruction]]). The pivot test measures conditioning; `Lambda` measures amplification; they select different cells.
2. The amplification itself is open. A rank reduction keyed on `Lambda <= 1`, not on the singular values, is the route the record names ([STATUS L2280-L2286](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2280-L2286), [[concepts/polyhedral-fit-amplification]]).
3. The cost of QR against Cholesky on hexahedra is not recorded.

## Related

[[hubs/advection]] - [[models/sl-reconstruction]] - [[concepts/polyhedral-fit-amplification]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[decisions/mesh-family-hexahedral]] - [[cases/popinet-translating-droplet]] - [[retractions/polyhedral-popinet-3d-mesh-defect]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the gate of 2026-09-05, explained by the census of 2026-09-09. Entered in [[decision-log#2026-09]].
