---
title: "slReconstruction: the semi-Lagrangian value fit"
description: "The production member is uncachedQuadraticWeightedLeastSquares, a constant-free quadratic value fit; value fitting and degree 2 are load-bearing, measured 2026-08-27 to 2026-09-09."
aliases: [slReconstruction, SL_RECONSTRUCTION, UQWLSR]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-05
decided_by: [config/uncachedConv2Dvortex.yaml, config/uncachedConv3Dshear.yaml, config/uncachedConv3DshearPoly.yaml, config/popinet3D_La12000_poly_dump4_qr.yaml]
code: [src/leiaLevelSet/semiLagrangian/slReconstruction.H, src/leiaLevelSet/semiLagrangian/slReconstruction.C, src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C]
sources: ["METHOD 2.2 (L98-L125)", "METHOD 8.1 rows SL_RECONSTRUCTION and SL_FIT (L373-L377)", "STATUS 4 pivot tolerance (L1550-L1722)", "STATUS 4 boundary faces (L1874-L1990)", "STATUS 4 QR (L2269-L2278)", "SL article sec:recon"]
---
# slReconstruction: the semi-Lagrangian value fit

> Verdict (2026-09-28). The production reconstruction of the level set at the departure foot is `uncachedQuadraticWeightedLeastSquares` ([METHOD 8.1 L373](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373)). It fits a constant-free quadratic to the stencil values with weights `1/|d|` and solves the normal system by Cholesky ([METHOD 2.2 L98-L116](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L98-L116)). Two properties are load-bearing: the fit interpolates values, and the degree is at least two ([METHOD L118-L121](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L118-L121)). The shape error converges at order 2.843 (2D vortex, CFL 0.5), 2.955 (3D shear, hex) and 3.285 (3D shear, poly) ([sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv)). The Taylor members and `defectCorrectedIDW` diverge ([iDEC report L119-L139](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/supplementary/iDEC_failure_report.tex#L119-L139)). The fit needs a geometric admissibility test on polyhedra: `quadraticPivotTol 0.3` ([STATUS L1657-L1681](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1657-L1681)).

## What it is

The semi-Lagrangian update is `psi^{n+1}(x_c) = psi^n(x_d)` ([METHOD 2 L63-L67](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L63-L67)). The reconstruction supplies `psi^n(x_d)`. Every member is centred on the arrival cell `c` and reproduces the cell value there, `evaluate(c, x_c) == psiOld[c]` ([slReconstruction.H L27-L37](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.H#L27-L37)). The family is selected by `levelSet { semiLagrangian { reconstruction ...; } }`; the code default is `quadraticWeightedLeastSquares` ([slReconstruction.C L166](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.C#L166)), the global token default is the same ([DP L29](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L29)), and the per-case layer sets the uncached member ([METHOD L15-L16](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L15-L16)).

The production fit, in the words of the SL article ([`sec:recon`, L304-L385](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L304-L305)): basis `b(d) = (d_a, d_a^2/2, d_a d_b)` with `m = 5` terms in 2D and `9` in 3D, weights `w_i = 1/|d_i|`, normal system `M a = r` with `M = A^T W^2 A`, solved by in-place Cholesky. The coefficient vector splits into the gradient `g` and the Hessian `H`, which the curvature path reads directly ([METHOD L114-L116](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L114-L116)).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| linear Taylor | `linearTaylor` | retracted for kinematics | The band gradient defect grows to 1.4e9 at N = 128 and 2.0e21 at N = 256 on the reversed vortex. It serves only the two-phase consistent-linear pipeline. | [LSL article L279-L282](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L279-L282), [SL article L327-L337](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L327-L337) |
| linear value fit | `linearWeightedLeastSquares` | open | The header calls it stable by construction; the study header says it is unstable in advection without the clip. See the open question. | [header L43-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/linearWeightedLeastSquaresReconstruction.H#L43-L51), [config L34-L36](https://github.com/leia-openfoam/leia/blob/8867581/config/linearConv2Dvortex.yaml#L34-L36) |
| signed-distance linear fit | `signedDistanceLinearWeightedLeastSquares` | candidate | Returns `P(x_foot)/|g|`, which pins the gradient magnitude to 1. No ladder is on record. | [header L36-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/signedDistanceLinearWeightedLeastSquaresReconstruction.H#L36-L51) |
| quadratic Taylor | `quadraticTaylor` (formerly `nestedLSQ`) | retracted for production | Builds the quadratic from twice-differentiated psi and needs the stencil clip to stay bounded. | [RM L24-L28](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L24-L28), [header L42-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/quadraticTaylorReconstruction.H#L42-L51) |
| cached quadratic value fit | `quadraticWeightedLeastSquares` (alias `quadraticWLSQ`) | settled, not production | Same arithmetic as the uncached member; the cached pseudo-inverse forces single precision and a memory ceiling in 3D. | [alias L46-L52](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/quadraticWeightedLeastSquaresReconstruction.C#L46-L52), [SL article L550-L560](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L550-L551) |
| uncached quadratic value fit | `uncachedQuadraticWeightedLeastSquares` | settled, production | Second to third order on hex and poly; no per-cell cache; the pivot test lives here. | [METHOD 8.1 L373](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L373), [header L30-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.H#L30-L51) |
| signed-distance quadratic fit | `signedDistanceQuadraticWeightedLeastSquares` (alias `sdQuadraticWLSQ`) | candidate | Returns the signed distance to the fitted zero set and a normal-strain rescale. Its kinematic study `sdCompare2D` predates the gradU fix. | [header L30-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/signedDistanceQuadraticWeightedLeastSquaresReconstruction.H#L30-L51), [STATUS L626-L629](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L626-L629) |
| band quadratic value fit | `bandQuadraticWeightedLeastSquares` (alias `bandQuadraticWLSQ`) | candidate | Fits the band only and re-extends the far field by Eikonal sweeps. No ladder is on record. | [header L29-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/bandQuadraticWeightedLeastSquaresReconstruction.H#L29-L51) |
| defect-corrected IDW | `defectCorrectedIDW` | retracted | The gradient defect grows 39 to 1.4e4 to 4.7e9 under refinement; the operator has spectral radius above one. | [[concepts/idec-defect-correction-failure]] |

The dictionary words are the `TypeName` strings of the headers ([linearTaylor L61](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/linearTaylorReconstruction.H#L61), [uncached L300](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.H#L300), [defectCorrectedIDW L137](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/defectCorrectedIDWReconstruction.H#L137)).

## The keys of the family

| key | values | default | what it decides | evidence |
|---|---|---|---|---|
| `stencil` | `point`, `face` | `point` | Cell-point-cell on hexahedra (6 face neighbours are too few for 9 coefficients); cell-face-cell on polyhedra. | [slReconstruction.C L94-L98](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.C#L94-L98), [METHOD L123-L125](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L123-L125) |
| `stencilBoundaryFaces` | `include`, `exclude`, `inflowOnly` | `include` | A boundary face is a data point at h/2 with the cell's own value. `exclude` diverged at step 716 in 2D. `inflowOnly` fixes the outlet layer but not the polyhedral far field. | [STATUS L1883-L1930](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1883-L1930), [STATUS L1965-L1989](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1965-L1989) |
| `quadraticPivotTol` | scalar | `0.3` | The smallest scaled Cholesky pivot a quadratic stencil must reach; below it the cell uses the linear fit. | [uncached C L184](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L184), [STATUS L1657-L1681](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1657-L1681) |
| `fit` | `normalEquations`, `householderQR` | `normalEquations` | QR blows up identically: step-3 phase volume 0.017512193 against 0.017512208. | [METHOD 8.1 L377](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L377), [STATUS L2269-L2278](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2269-L2278) |
| `ridgeEps` | scalar | `0` | Optional Tikhonov ridge on the diagonal. | [uncached C L183](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L183) |
| `fitProbeDisplacement`, `writeFitOrder` | vector, bool | `(0 0 0)`, `false` | Diagnostics: the amplification bound `Lambda` and the per-cell fit order. Inert at the defaults. | [uncached C L185-L186](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C#L185-L186), [[concepts/polyhedral-fit-amplification]] |
| `clipToStencilBounds`, `clipRegion`, `clipKeepExtrema` | legacy clip | `false`, `all`, `false` | Superseded by the `slValueBound` family. | [[models/sl-value-bound]] |

## Why it matters

The reconstruction is the only place where the transport can create error. The foot is exact for a uniform velocity ([kinematicTranslation2D header L19-L23](https://github.com/leia-openfoam/leia/blob/8867581/config/kinematicTranslation2D.yaml#L19-L23)). A value fit stays bounded within the stencil data; a Taylor expansion injects a differentiated field and amplifies grid-scale error ([SL article L418-L427](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L418)). The same fit supplies the curvature, so a defect here reaches the capillary force ([[concepts/curvature-from-the-fit]]).

## Where in the code

- Base class and stencil: `src/leiaLevelSet/semiLagrangian/slReconstruction.{H,C}` ([L27-L37](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.H#L27-L37)).
- Production member: `uncachedQuadraticWeightedLeastSquaresReconstruction.{H,C}`; the admissibility test is in `build` ([STATUS L1577-L1589](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1577-L1589)).
- Unit gate: `applications/test/leiaTestSLReconstruction`.
- Tokens: `SL_RECONSTRUCTION`, `SL_STENCIL`, `SL_FIT`, `SL_QUAD_PIVOT_TOL`, `SL_WRITE_FIT_ORDER`, `SL_STENCIL_BOUNDARY_FACES` ([DP L29-L41](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L29-L41), [DP L1074-L1093](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1074-L1093)).

## Evidence

| claim | number | where |
|---|---|---|
| Shape order, 2D reversed vortex, hex, CFL 0.5 / 1.0 | 2.843 / 2.378 | MEASURED, [sl_convergence_orders.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_convergence_orders.csv), re-established after the gradU fix [PCS L52-L61](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61) |
| Shape order, 3D shear, hex / poly | 2.955 / 3.285 | MEASURED, same table; [SL article L1253-L1255](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1253-L1255) |
| Shape order, 3D deformation, hex / poly (filament-limited) | 1.360 / 1.464 | MEASURED, same table |
| Accuracy against the best Eulerian line at 512^2, T = 8 | 20x more accurate at half the wall clock | MEASURED, [MC article L195-L199](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L194) |
| Ill-conditioned polyhedral stencils: neighbours, condition, coefficients | 10 point-neighbours for 9 coefficients; condition 3e7 to 6e12; coefficients 3e2 to 3e8 | MEASURED, [STATUS L1566-L1575](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1566-L1575) |
| Tolerance history: divergence step of the 78-step polyhedral smoke | 8 / 12 / 16 at tolerance 0 / 1e-3 / 1e-2 | MEASURED, [STATUS L1642-L1643](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1642-L1643) |
| Pivot census: band cells, hex minimum, polyhedral interior | band >= 0.74 (poly >= 0.94); hex >= 0.64; poly interior >= 0.50 | MEASURED, [sl_fit_pivot_census.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_fit_pivot_census.csv), [STATUS L1670-L1676](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1670-L1676) |
| Cells demoted at 0.3 on the uniform polyhedral boxes | 7.1 %, none within 12h of the interface | MEASURED, [STATUS L1662-L1668](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1662-L1668) |
| Hex bit-identity with the pivot test | cmp-identical over 65 steps, 0 demoted cells | MEASURED, [STATUS L1695-L1696](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1695-L1696) |
| Boundary-normal transport fraction with `include`, poly / hex | inlet 0.082 / 0.198; outlet 0.240 / 0.491 | MEASURED, [STATUS L1877-L1884](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1877-L1884) |
| `inflowOnly` on the 2D Popinet horizon | every droplet metric identical to `include` to all printed digits | MEASURED, [STATUS L1924-L1930](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1924-L1930) |
| QR against normal equations on the diverging polyhedral case | phase volume 0.017512193 vs 0.017512208 at step 3; the worst amplifier has pivot 0.757 | MEASURED, [STATUS L2269-L2278](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2269-L2278) |
| Non-orthogonal correctors 1 / 3 / 6 on the polyhedral cases | metrics identical to every printed digit | MEASURED, [STATUS L1941-L1957](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1941-L1957) |

## Decisions

- `SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares` per case: [[decisions/sl-reconstruction-uncached-qwls]].
- `SL_FIT normalEquations`: [[decisions/sl-fit-normal-equations]].
- `SL_QUAD_PIVOT_TOL 0.3`, decided 2026-09-08 from the pivot census ([STATUS L1680-L1681](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1680-L1681)).
- `SL_STENCIL_BOUNDARY_FACES include` stays the default; `inflowOnly` is selectable ([STATUS L1987-L1989](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1987-L1989)).

## Open questions

1. Which linear member is unstable? METHOD says a linear value fit drives the gradient error to 1e21 ([L120-L121](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L120-L121)). The LSL article attributes 2.0e21 at N = 256 to `linearTaylor` ([L279-L282](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L279-L282)). The `linearConv2Dvortex` study sweeps `linearWeightedLeastSquares` and its header calls that member unstable ([L28-L36](https://github.com/leia-openfoam/leia/blob/8867581/config/linearConv2Dvortex.yaml#L28-L36)). See [[concepts/linear-semi-lagrangian]].
2. The unbounded scheme saturates on uniform translation at N = 256: the error rises from 7.947e-05 to 1.159e-04, order -0.54 ([METHOD 8.3.7 L641-L653](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L641-L653)). The cause is not found.
3. The transport operator amplifies on every mesh: `rho(B) = 1.00441` on production hexahedra ([METHOD 8.2 L409-L417](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L409-L417)). See [[concepts/polyhedral-fit-amplification]].
4. No ladder is on record for `signedDistanceLinearWeightedLeastSquares` and `bandQuadraticWeightedLeastSquares`.

## Related

[[hubs/advection]] - [[models/sl-scheme]] - [[models/sl-value-bound]] - [[models/level-set-advection]] - [[concepts/departure-foot-ab2-centring]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/idec-defect-correction-failure]] - [[concepts/linear-semi-lagrangian]] - [[concepts/value-bounds-and-clips]] - [[concepts/curvature-from-the-fit]] - [[concepts/advection-regression-set]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[decisions/sl-fit-normal-equations]] - [[retractions/gradu-coupled-patch-contamination]] - [[retractions/polyhedral-popinet-3d-mesh-defect]] - [[studies/sl-quadratic-pre-print]]

## Log

### 2026-09-28
Created from METHOD 2.2 and 8.1, STATUS section 4, the SL article and the code headers.
