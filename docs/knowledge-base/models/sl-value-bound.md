---
title: "slValueBound and the clips"
description: "No bound is in the best configuration: stencilBounds was falsified by gate G4 (2026-09-09) and lipschitzCone by the advection ladders (2026-09-10); the sentinel fromClipSwitch resolves to none."
aliases: [slValueBound, SL_VALUE_BOUND, SL_CLIP, clipToStencilBounds]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-10
decided_by: [config/popinet3D_poly_sigma0_clipGate.yaml, config/popinet2D_clipRegionGate.yaml, config/popinet2D_coneBoundGate.yaml, config/popinet2D_coneBoundLadderN128.yaml, config/advConv2Dtranslation.yaml, config/advConv2Dvortex.yaml]
code: [src/leiaLevelSet/semiLagrangian/slValueBound.H, src/leiaLevelSet/semiLagrangian/noValueBound.H, src/leiaLevelSet/semiLagrangian/stencilBoundsValueBound.H, src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H]
sources: ["METHOD 8.1 rows SL_CLIP to SL_CONE_INADMISSIBLE (L378-L384)", "METHOD 8.3 (L488-L749)", "METHOD 9 item 4 (L777-L790)", "STATUS 4 sigma-0 control (L1822-L1844)", "STATUS 4 clip decided (L2008-L2039)", "STATUS 4 G4 and the retraction (L2312-L2495)", "DP L57-L158", "STATUS 11.19 (L4234-L4252)"]
---
# slValueBound and the clips

> Verdict (2026-09-28). The bound on the reconstructed departure value is a runtime-selected family with three members: `none`, `stencilBounds` and `lipschitzCone` ([slValueBound.H L30-L52](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.H#L30-L52)). No bound is in the best configuration ([METHOD 8.1 L381](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L381)). `stencilBounds`, the quasi-monotone clip, removes the polyhedral far-field failure but costs +30.4 % volume error on the 2D hexahedral translating droplet; with the extremum exemption it fails at step 506 against the control's 527 ([METHOD 8.1 L378](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L378), [STATUS L2319-L2338](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2319-L2338)). `lipschitzCone` passes the exact-field gate and improves every coupled metric at N = 64, but the advection ladders falsify it: 189.7x worse on the vortex at N = 256 ([METHOD 8.3.7 L639-L667](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L639-L667)); the retraction is [[retractions/distance-cone-bound-as-transport-bound]]. CORRECTED 2026-09-29: the "3.6x worse on uniform translation at N = 64" of that ladder came from a reversed flow and is void ([[retractions/reversed-2dtranslation]]). One way, the cone bound is 14.1x worse at N = 64 and 169.6x at N = 256, with orders 0.38, 0.08, 0.36 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697)). Hard limiters collapse the order: Barth-Jespersen 3.0 to 0.1, Venkatakrishnan 3.0 to 0.9 ([METHOD 9 L786-L788](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L786-L788)).

## What it is

The bound decides one interval `[lo, hi]` per arrival cell; the value `psi^{n+1}` must lie in it. It never touches the reconstruction, the feet or the field ([slValueBound.H L31-L34](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.H#L31-L34)). It acts in `slCorrector::robustEvaluate`, so it reaches the `pointValue` scheme only ([slValueBound.H L64-L67](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.H#L64-L67)). The family exists because the quadratic fit is not a convex combination of the stencil values, so it can create a new extremum, and the frozen-velocity power iteration measures a spectral radius above one on every mesh ([slValueBound.H L54-L62](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.H#L54-L62), [[concepts/polyhedral-fit-amplification]]).

The dictionary key is `valueBound`; the default is the sentinel `fromClipSwitch` ([slValueBound.C L57](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.C#L57)). The sentinel resolves to `stencilBounds` when the legacy `clipToStencilBounds` is true and to `none` otherwise ([slValueBound.H L49-L52](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slValueBound.H#L49-L52)). `SL_CLIP` is `false`, so the sentinel resolves to `none` ([DP L57](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L57), [DP L104](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L104)).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| no bound | `none` | settled, production | Only the runaway cap (stencil mid +- 10 stencil ranges) that predates the family; not stable in the `rho <= 1` sense, `rho(B)` = 1.00441 on production hexahedra. | [noValueBound.H L68](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/noValueBound.H#L68), [DP L106-L111](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L106-L111) |
| quasi-monotone clip | `stencilBounds` | retracted as a fix | Clips to the stencil [min, max]; 59.2 % of the cells it must act on are stencil extrema; the exemption returns the failure at step 506. | [stencilBoundsValueBound.H L36-L48](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/stencilBoundsValueBound.H#L36-L48), [STATUS L2340-L2353](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2340-L2353) |
| distance-cone bound | `lipschitzCone` | retracted as a transport bound | Intersection of the intervals `psi_j +- L abs(x_d - x_j)`; exact at a distance cusp; worse than no bound on every advection gate. | [lipschitzConeValueBound.H L30-L45](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H#L30-L45), [METHOD 8.3.6 L739-L742](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L739-L742) |

The dictionary words are the `TypeName` strings ([none L68](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/noValueBound.H#L68), [stencilBounds L74](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/stencilBoundsValueBound.H#L74), [lipschitzCone L183](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H#L183)).

## The keys of the family

| key | values | default | read by | evidence |
|---|---|---|---|---|
| `clipRegion` (`SL_CLIP_REGION`) | `all`, `outsideBand` | `all` | `stencilBounds` only; the band exclusion removes zero firings on hex | [slReconstruction.C L57](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.C#L57), [METHOD 8.1 L379](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L379) |
| `clipKeepExtrema` (`SL_CLIP_KEEP_EXTREMA`) | bool | `false` | `stencilBounds` only; exempts a cell that is its stencil's extremum | [slReconstruction.C L63](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slReconstruction.C#L63), [METHOD 8.1 L380](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L380) |
| `lipschitzMode` (`SL_CONE_L_MODE`) | `unity`, `stencil` | `unity` | `lipschitzCone` only; `stencil` is a diagnostic arm and loses the phase on 3D shear | [lipschitzConeValueBound.C L50](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.C#L50), [METHOD 8.1 L382](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L382) |
| `lipschitzConstant` (`SL_CONE_L`) | scalar | `1` | the eikonal value; any other value is a tuned coefficient and the solver warns | [lipschitzConeValueBound.C L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.C#L51), [METHOD 8.1 L383](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L383) |
| `onInadmissible` (`SL_CONE_INADMISSIBLE`) | `cellOnly`, `none` | `cellOnly` | what to do with an empty interval; 0.58 % of cell-steps were empty and the choice moved the volume error from -67.8 % to -5.3 % | [lipschitzConeValueBound.C L52](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.C#L52), [METHOD 8.1 L384](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L384) |
| `coneBoundaryFaces` | bool | `false` | a zeroGradient face value is not a sample of the distance function | [lipschitzConeValueBound.H L117-L125](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H#L117-L125) |

Diagnostics written at write time: `slClipFired`, `slClipFiredEver`, `slBoundDelta`, `slBoundSlack`, `slBoundInadmissible` ([README L735-L739](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L735-L739)). Read the cumulative firing counter, never the per-step sample ([STATUS L2475-L2481](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2475-L2481)).

## Why it matters

The polyhedral far field fails without a bound: the sigma = 0 passive control on the Popinet 3D pMesh grows a fake zero set at step 508 while the velocity stays an exact uniform stream ([STATUS L1822-L1827](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1822-L1827)). A bound that fixes this on polyhedra must be inert on the hexahedral interface. No formulation on record does both ([STATUS L2376-L2381](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2376-L2381)). The reason is structural: a monotone bound cannot represent an extremum, and a signed distance has genuine extrema on its medial axis ([STATUS L2417-L2421](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2417-L2421)). There is also an impossibility result: a nonnegative, partition-of-unity, quadratic-exact rule applied to `p(x) = abs(x - x_d)^2` forces a sum of positive terms to zero, so quadratic exactness and `Lambda = 1` exclude each other at an off-node point ([METHOD 9 L781-L785](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L781-L785)).

## Where in the code

- Family: `src/leiaLevelSet/semiLagrangian/slValueBound.{H,C}`, members `noValueBound`, `stencilBoundsValueBound`, `lipschitzConeValueBound`.
- The cone interval itself: `slReconstruction::stencilConeRange` ([METHOD 8.3.1 L507-L508](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L507-L508)).
- Unit gate: `applications/test/leiaTestSLReconstruction`, section (d) ([METHOD 8.3.1 L510-L518](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L510-L518)).
- Nonlinear amplification instrument: `leiaTestTransportSpectrum -mode growth` ([METHOD 8.3.3 L544-L560](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L544-L560)).

## Evidence

| claim | number | where |
|---|---|---|
| Global clip on the sigma = 0 polyhedral control | fake zero set NONE against step 508; zero-set error 3.98e-04 against 8.66e-01; volume error 1.08e-06 against 3.79 | MEASURED, [STATUS L2017-L2021](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2017-L2021) |
| Global clip on Popinet 2D hex, N = 64 | volume error +30.4 %, shape +4.9 %, centroid +12.5 %, L1 abs u' +8.9 % | MEASURED, [STATUS L2403-L2407](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2403-L2407) |
| Band-aware clip (`outsideBand`) | identical damage to `all` to every printed digit; the six firing cells are the four box corners and the two apex cells | MEASURED, [STATUS L2398-L2415](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2398-L2415) |
| Extremum exemption on hex | damage cut 30x: volume +1.03 %, 3 cell-steps of clipping in total | MEASURED, [STATUS L2431-L2438](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2431-L2438) |
| G4: clip + exemption on the polyhedral control | failure at step 506 against 527 without a clip; the exemption withheld 59.2 % of the firings (2 680 917 of 6 578 187 cell-steps) | MEASURED, [STATUS L2319-L2342](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2319-L2342) |
| Family inertness | 8 arms of `popinet2D_clipRegionGate` byte-identical over 1563 steps at np = 4 | MEASURED, [METHOD 8.1 L381](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L381), [STATUS L280-L286](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L280-L286) |
| Cone bound, exact-field gate | 0 empty intervals; truth inside the interval to 1.1e-16; monotone clip errs by one cell width at 96 apex cells, the cone by 0 | MEASURED, [METHOD 8.3.1 L510-L518](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L510-L518) |
| Cone bound, coupled 2D hex, N = 64 | every interface metric improves 25 to 68 %; grad-psi band error -75 to -81 % | MEASURED, one resolution, [METHOD 8.3.2 L525-L540](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L525-L540) |
| Cone bound, nonlinear growth per step | 3.48e-04 against 1.456e-04 for `none` at amp 0.1 h: 1.5 to 2.4x more amplification | MEASURED, [METHOD 8.3.3 L550-L567](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L550-L567) |
| Cone bound, converged 2D ladders | translation, one way (2026-09-29): 2.5x / 14.1x / 60.1x / 169.6x worse at N = 32 / 64 / 128 / 256, cone orders 0.38, 0.08, 0.36 (no convergence); the reversed 0.85x / 3.6x / 4.3x / 1.9x of 2026-09-10 are VOID; vortex: 8.6x to 189.7x worse, order 0.79 against 2.09 | MEASURED, [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697), [METHOD 8.3.7 L639-L667](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L639-L667) |
| Cone bound, one-way translation, four arms on one N = 64 mesh (digest fd13ce29ac8e) | `none` 7.21e-03, cone unity 1.02e-01 (14.1x), cone stencil 6.82e-02 (9.5x), `stencilBounds` 7.21e-03; the cone arm gives 3.03e-02 with the exact outflow value and 1.02e-01 with zeroGradient | MEASURED, [STATUS 11.19](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4248-L4252) |
| Cone bound, resolution ladder N = 64 to 128 | grad-psi error 7.119e-03 to 7.066e-03, order 0.01 (a floor); centroid +119.4 % worse at N = 128 | MEASURED, [METHOD 8.3.5 L698-L724](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L698-L724) |
| Cone bound, 3D shear hex, `stencil` mode | E_VOL_ALPHA_REL = 1.0000, the phase is gone | MEASURED, [METHOD 8.3.4 L592-L597](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L592-L597) |
| Polyhedral 3D shear, unbounded | diverged at step 198 in both HEAD and the pre-change binary | MEASURED, [METHOD 8.3.8 L680-L686](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L680-L686) |
| Mesh-noise floor on the polyhedral rung | volume spread 0.0 / 0.0 / 14.2 % for none / stencilBounds / lipschitzCone | MEASURED, [CLAUDE L617-L623](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L617-L623) |

## Why it failed, or why we think so

1. `stencilBounds`: the fit undershoots at the apex of the distance cone because a smooth quadratic cannot follow a non-differentiable minimum; the clip then flattens the extremum at every step ([STATUS L2417-L2421](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2417-L2421)). The exemption hands the growing checkerboard exactly the cells it needs ([STATUS L2340-L2346](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2340-L2346)).
2. `lipschitzCone`: in strained flow the true Lipschitz constant grows as `L(t) <= L(0) exp(int abs(grad u) dt)`; the implementation clamps `L = 1` and destroys the field ([METHOD 8.3.4 L606-L615](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L606-L615)). The flow-map `L` needs a collective per step and is deliberately not implemented ([lipschitzConeValueBound.H L87-L97](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H#L87-L97)).

## Decisions

- `SL_CLIP false`, `SL_VALUE_BOUND fromClipSwitch`: [[decisions/sl-clip-and-value-bound-off]].
- The family survives as one study axis; the review's Rank 1 (monotone transporting base plus filtered correction) is the next candidate ([METHOD 8.3.6 L733-L737](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L733-L737)).

## Open questions

1. What separates a genuine extremum from a checkerboard peak is scale, not being an extremum; two threshold-free properties are proposed and untested ([STATUS L2362-L2374](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2362-L2374)).
2. Arm 00002 of G4 (global clip) still grows slowly: 4.46e-04 at step 1 to 8.97e-04 at step 649, just below the 1e-3 threshold ([STATUS L2383-L2387](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2383-L2387)).
3. The curated `advConv2D*_convergence.csv` still carry the wrong `hEff` column and must be regenerated ([STATUS L3241-L3242](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3241-L3242), [[retractions/advection-orders-3-2-factor]]). CORRECTED 2026-09-29: the one-way [translation table](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/docs/method-comparison/method-comparison-article/data/tables/advConv2Dtranslation_convergence.csv) is regenerated with h = 1/N; the [vortex table](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/docs/method-comparison/method-comparison-article/data/tables/advConv2Dvortex_convergence.csv) still carries the old column (0.0992 at N = 32).
4. Why the cone arm depends on the psi boundary value 0.1 away from the interface is not measured; the bound is falsified anyway ([STATUS 11.19 item 8](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4315-L4316)).

## Related

[[hubs/advection]] - [[concepts/value-bounds-and-clips]] - [[concepts/polyhedral-fit-amplification]] - [[models/sl-reconstruction]] - [[models/sl-scheme]] - [[concepts/advection-regression-set]] - [[decisions/sl-clip-and-value-bound-off]] - [[retractions/distance-cone-bound-as-transport-bound]] - [[retractions/clip-damage-is-the-narrow-band]] - [[retractions/advection-orders-3-2-factor]] - [[retractions/reversed-2dtranslation]] - [[cases/popinet-translating-droplet]]

## Log

### 2026-09-28
Created from METHOD 8.1, 8.3 and 9, STATUS section 4 and the family headers.

### 2026-09-29
The translation numbers of the cone bound are re-measured on the one-way `2Dtranslation` (STATUS 11.19, [[retractions/reversed-2dtranslation]]): the verdict and the ladder row CORRECTED (14.1x at N = 64 and no convergence, where the void reversed run gave 3.6x), a row for the four shared-mesh arms, open question 3 CORRECTED and 4 added. The falsification stands and is stronger.
