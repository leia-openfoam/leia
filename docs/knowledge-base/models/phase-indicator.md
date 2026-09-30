---
title: "phaseIndicator: heaviside, sharpJump, geometric, detrixheAslam"
description: "detrixheAslam is production: order about 2.0 on a static circle and equal to the geometric clip to eight digits; the local plane needs no signed distance, but its first-order offset still assumes one (open)."
aliases: [phaseIndicator, PHASE_INDICATOR, Detrixhe-Aslam]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-01
decided_by: [config/phaseIndicatorConvergence.yaml, "author decision recorded in cases/default.parameter L5-L9"]
code: [src/leiaLevelSet/phaseIndicator/phaseIndicator.C, src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.H, src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.H, src/leiaLevelSet/phaseIndicator/levelSetPlaneReconstruction.C]
sources: ["METHOD 3 (L129-L140)", "METHOD 8.1 row PHASE_INDICATOR (L388)", "DP L5-L9", "SL article sec:indicator (L621-L680)", "SL article sec:indicator-verif (L1531-L1544)", "PCS 11.3 (L863-L878)", "RM L170-L177 and L1339-L1352", "STATUS 11.3 (L3261-L3269)"]
---
# phaseIndicator: heaviside, sharpJump, geometric, detrixheAslam

> Verdict (2026-09-28). The phase indicator turns the level set into a volume fraction `alpha`. The production member is `detrixheAslam` ([METHOD 8.1 L388](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L388), [DP L9](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L9)). It fits a local signed-distance plane by linear least squares, fans the cell into tetrahedra and fills each with the closed form of Detrixhe and Aslam ([METHOD 3 L131-L140](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L131-L140)). On a static circle the volume error converges at order about 2.0, and the geometric clip agrees with it to eight digits ([SL article L1537-L1540](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1531)). The switch from `geometric` to `detrixheAslam` is NOT inert ([DP L5-L8](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L5-L8)). Open: the face fraction still uses the first-order offset `psi/abs(grad psi)`, which assumes a signed distance ([PCS L863-L878](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L863-L878)).

## What it is

`alpha_c` is the cell average of `H(-psi)`, and `alpha = 1` in `{psi < 0}` ([SL article L623-L629](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621), [volumeCorrection.H L56-L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L56-L60)). The plane-based members recover the interface from psi without the assumption `abs(grad psi) = 1`: the plane `(n_c, d_c)` minimises the least-squares residual over the cell and its neighbours, and its normalised form is a signed distance reconstructed from psi ([SL article L633-L649](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621)). The key is `levelSet.phaseIndicator.type`; the code default is `geometric` ([phaseIndicator.C L55](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/phaseIndicator.C#L55)) and the token default is `detrixheAslam` ([DP L9](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L9)). The indicator is evaluated in the narrow band only ([geometricPhaseIndicator.C L123-L126](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C#L123-L126)).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| smoothed Heaviside | `heaviside` | settled, not production | A smoothed Heaviside of psi over `nCells`; first order at the interface and it presumes a signed distance. | [heavisidePhaseIndicator.H L30-L33](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/heavisidePhaseIndicator.H#L30-L33), [SL article L630-L633](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621) |
| sharp jump | `sharpJump` | settled, not production | An abrupt jump with `nAverages` averaging passes. No ladder is on record. | [sharpJumpPhaseIndicator.H L30-L55](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/sharpJumpPhaseIndicator.H#L30-L55) |
| geometric clip | `geometric` | settled, the code default | The same least-squares plane truncates the cell by polygon clipping; sensitive to geometric tolerances and known issues in 2D. | [geometricPhaseIndicator.H L30-L37](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.H#L30-L37), [detrixheAslam H L37-L39](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.H#L37-L39) |
| Detrixhe-Aslam tetrahedral fill | `detrixheAslam` | settled, production | Same plane, analytic per-tetrahedron fraction, tolerance-free, identical in 3D and one-cell-thick 2D. | [detrixheAslam H L30-L50](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.H#L30-L50), [SL article L651-L680](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621) |

The dictionary words are the `TypeName` strings ([heaviside L64](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/heavisidePhaseIndicator.H#L64), [sharpJump L61](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/sharpJumpPhaseIndicator.H#L61)) and the registrations ([geometric](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C#L51), [detrixheAslam](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.C#L41)). The base class declares `TypeName("none")` but never registers it ([phaseIndicator.H L69](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/phaseIndicator.H#L69)).

`detrixheAslam` reads `geometrySource levelSetField | analyticImplicitSurface`; the second evaluates the exact implicit surface instead of the transported psi ([detrixheAslam C L97-L130](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/detrixheAslamPhaseIndicator.C#L97-L130)).

## Why it matters

Every shape and volume metric of the verification reads this indicator ([SL article L1533-L1534](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1531)), and so does the capillary force support `snGrad(alpha)` ([[concepts/balanced-force-csf-flux]]). The volume correction restores the volume this indicator measures ([[models/volume-correction]]). A mixture of two indicators in one ladder once came within one run of a fitted order ([DP L5-L8](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L5-L8)).

## Where in the code

- Family: `src/leiaLevelSet/phaseIndicator/`; the shared plane fit is `leastSquaresPlaneCoeffs` in `levelSetPlaneReconstruction.C`.
- The determinant guard runs on a centred, variance-normalised copy of the 4x4 matrix; the solve uses the absolute system ([levelSetPlaneReconstruction.C L205-L225](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/levelSetPlaneReconstruction.C#L205-L225), guard at [L321-L324](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/levelSetPlaneReconstruction.C#L321-L324)).
- The first-order offset of the face fraction: [geometricPhaseIndicator.C L130-L141](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C#L130-L141), `faceAreaFraction.H` ([PCS L869-L870](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L869-L870)).
- The same planes seed the GRL redistancer ([[models/redistancer]]) and the geometric mass flux ([[models/mass-flux]]).

## Evidence

| claim | number | where |
|---|---|---|
| Static circle r = 0.15, N = 32 to 256: relative volume error order | about 2.0; geometric and Detrixhe-Aslam agree to eight digits | MEASURED, [SL article L1537-L1540](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1531) |
| The closed form is tolerance-free | every denominator is a product of differences that straddle the interface | DERIVED, [SL article L673-L679](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621) |
| The middle-case numerator of the published formula returns the complementary fraction | corrected in the implementation | DERIVED, [SL article L676-L679](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L621) |
| Determinant guard at N = 512 on the 0.01 m box | every band cell returned the flat-plane fallback before the fix; the N = 512 curvature row was biased low, 10.2 against 17.4 | MEASURED, [RM L1339-L1352](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1339-L1352) |
| Transported level set, N = 32 translating droplet, no disturbance velocity | Detrixhe-Aslam phase volume grows about 6 %; centroid lags about 0.26 h | MEASURED, one resolution, [RM L175-L177](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L175-L177) |
| Chord bias of the plane | O(h^2 kappa), one-signed; caps alpha at second order | DERIVED, [MC article L251-L255](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L251-L255) |
| Band metrics read the unlimited gradient | `gradPsiMetric` (unlimited leastSquares) since 2026-09-26; identical at tolerance 0 on a smooth band | MEASURED, [STATUS L3261-L3269](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3261-L3269) |

## Decisions

- `PHASE_INDICATOR detrixheAslam`: [[decisions/phase-indicator-detrixhe-aslam]]. The comparison config is `config/phaseIndicatorConvergence.yaml` ([README L21](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L21)).

## Open questions

1. The face fraction's offset `d = psi/abs(grad psi)` carries the error `e_d = -(1/2) b_n d^2`; the Hessian-corrected root exists as `offsetDistance` and would cost no extra derivative order. It moves the force support, so it must be re-gated on transport and volume ([PCS L863-L878](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L863-L878)). Not measured.
2. `cases/oscillatingDroplet2D` initialised an algebraic psi (`implicitEllipsoid`, `abs(grad psi)` 1.8e3 to 2.2e3 at the interface), so its gradient columns had no meaning; the token `DROPLET_SURFACE` now exists and the past studies are not yet voided ([STATUS L3271-L3279](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3271-L3279), [[cases/oscillating-droplet]]).
3. The sign-change band can miss a cell clipped at a corner whose neighbours do not change sign ([signChangeNarrowBand.C L58-L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/signChangeNarrowBand.C#L58-L60)); see [[models/narrow-band]].

## Related

[[hubs/advection]] - [[models/narrow-band]] - [[models/volume-correction]] - [[models/redistancer]] - [[models/mass-flux]] - [[concepts/balanced-force-csf-flux]] - [[concepts/curvature-from-the-fit]] - [[concepts/redistancing-geometric-grl]] - [[decisions/phase-indicator-detrixhe-aslam]] - [[cases/oscillating-droplet]] - [[studies/sl-quadratic-pre-print]]

## Log

### 2026-09-28
Created from METHOD 3, the SL article, the plan and roadmap entries and the code.
