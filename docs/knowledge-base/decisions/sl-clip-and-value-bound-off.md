---
title: "SL_CLIP false and SL_VALUE_BOUND none: no bound on the reconstructed value"
description: "SL_CLIP false and SL_VALUE_BOUND fromClipSwitch (resolves to none), decided by gate G4 on 2026-09-09 and the advection ladders on 2026-09-10: the clip with the extremum exemption fails at step 506 against 527 without it, the global clip costs +30.4 % volume, and the cone bound is 3.6x to 189.7x worse on pure advection"
aliases: [SL_CLIP false, SL_VALUE_BOUND none]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-28
date_settled: 2026-09-10
decided_by: [config/popinet3D_poly_sigma0_clipGate.yaml, config/popinet2D_clipRegionGate.yaml, config/popinet2D_coneBoundGate.yaml, config/popinet2D_coneBoundLadderN128.yaml, config/advConv2Dtranslation.yaml, config/advConv2Dvortex.yaml, config/advConv3DshearHex.yaml, config/advConv3DshearPoly.yaml]
code: [src/leiaLevelSet/semiLagrangian/slValueBound.H, src/leiaLevelSet/semiLagrangian/stencilBoundsValueBound.H, src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.H, cases/default.parameter]
sources: ["METHOD 8.1 rows SL_CLIP to SL_CONE_INADMISSIBLE (L378-L384)", "METHOD 8.3 (L490-L749)", "METHOD 9 item 4 (L777-L790)", "STATUS 4 (L2008-L2039, L2312-L2495)", "DP L57-L158"]
---
# SL_CLIP false and SL_VALUE_BOUND none: no bound on the reconstructed value

> `SL_CLIP false`, `SL_CLIP_REGION all`, `SL_CLIP_KEEP_EXTREMA false` and `SL_VALUE_BOUND fromClipSwitch` in the global default ([DP L57-L106](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L57-L106)). The sentinel follows `SL_CLIP`, so the production bound is `none` ([METHOD 8.1 L381](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L381)). Decided in two steps. Gate G4, 2026-09-09 ([STATUS L2312-L2342](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2312-L2342)): on the polyhedral sigma = 0 control the quasi-monotone clip with the extremum exemption fails at step 506 against 527 for no clip, because 59.2 % of the cells the clip must act on are stencil extrema; the global clip removes the polyhedral failure but costs +30.4 % volume error on the 2D hexahedral translating droplet at N = 64 ([METHOD 8.1 L378](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L378)). The advection ladders, 2026-09-10 ([METHOD 8.3.7 L639-L667](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L639-L667)): the distance-cone bound is 3.6x worse than no bound on uniform translation at N = 64 and 189.7x worse on the vortex at N = 256; the earlier coupled gain at one resolution is retracted the same day ([METHOD 8.3 L490-L498](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L490-L498), [[retractions/distance-cone-bound-as-transport-bound]]). The token comment of `SL_CLIP` still calls the clip required on polyhedra; the correction of 2026-09-28 marks it refuted ([DP L64-L65](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L64-L65)).

## The question

The quadratic value fit is not a convex combination of its stencil values, so it can create a new extremum, and on cfMesh's small one-sided cells the far field grows a false zero set: the sigma = 0 passive control fails at step 508 while the velocity stays an exact uniform stream ([STATUS L1822-L1827](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1822-L1827)). Can a coefficient-free bound on the reconstructed value remove that failure without a cost at the interface? Two candidates were gated: the clip to the stencil interval (`stencilBounds`), and the intersection of the Lipschitz intervals `psi_j +- L abs(x_d - x_j)` with `L = 1` (`lipschitzCone`) ([DP L106-L126](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L106-L126), [[models/sl-value-bound]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| G4, poly sigma = 0, N = 64, np 32, 650 steps: clip off (two arms) | first step with a fake zero set | 527 in both, bit-identical: the tokens are unread when the clip is off | MEASURED, [STATUS L2319-L2338](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2319-L2338) |
| G4: global clip, no exemption | fake zero set; zero-set error at T; volume error at T | none in 650 steps; 8.97e-04; 3.93e-06 | MEASURED, [STATUS L2319-L2338](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2319-L2338) |
| G4: clip with the extremum exemption (the candidate) | first step with a fake zero set | 506, 4 % before the control, inside the 5 to 38 % scatter of genuine instabilities | MEASURED, [STATUS L2319-L2342](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2319-L2342) |
| G4: the exemption | share of the firings withheld | 59.2 % (2 680 917 of 6 578 187 cell-steps over 1950 corrector calls) | MEASURED, [STATUS L2340-L2346](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2340-L2346) |
| G5, Popinet 2D hex, N = 64, 1563 steps: global clip | volume error; shape; centroid; L1 abs u' | +30.4 %; +4.9 %; +12.5 %; +8.9 % | MEASURED, [STATUS L2403-L2407](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2403-L2407) |
| G5: band-aware clip (`outsideBand`) | the same vector | identical to the global clip to every printed digit; the firing cells are the four box corners and the two apex cells | MEASURED, [STATUS L2398-L2415](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2398-L2415) |
| G5: exemption on hex | volume error; firings | +1.03 %; 3 cell-steps in the run | MEASURED, [STATUS L2431-L2438](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2431-L2438) |
| B1: the family sentinel | 8 arms of `popinet2D_clipRegionGate`, 1563 steps, np 4 | byte-identical to the pre-family study | MEASURED, [METHOD 8.1 L381](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L381) |
| B7a, cone bound, coupled 2D hex, N = 64 | interface metrics | every metric improves 25 to 68 %, one resolution only | MEASURED, [METHOD 8.3.2 L525-L540](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L525-L540) |
| converged 2D advection ladders, cone bound against none | geometric error ratio | translation 3.6x / 4.3x / 1.9x worse at N = 64 / 128 / 256 (0.85x better at N = 32); vortex 8.6x to 189.7x worse, order 0.79 against 2.09 | MEASURED, [METHOD 8.3.7 L639-L667](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L639-L667) |
| resolution ladder N = 64 to 128, cone bound | band grad-psi error; centroid error | 7.119e-03 to 7.066e-03, order 0.01 (a floor); +119.4 % worse than none at N = 128 | MEASURED, [METHOD 8.3.5 L698-L724](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L698-L724) |
| 3D shear hex, cone bound in `stencil` mode | `E_VOL_ALPHA_REL` | 1.0000, the phase is gone | MEASURED, [METHOD 8.3.4 L592-L597](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L592-L597) |
| hard limiters on the fit | convergence order | Barth-Jespersen 3.0 to 0.1; Venkatakrishnan 3.0 to 0.9 | MEASURED, [METHOD 9 L786-L788](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L786-L788) |

Pre-registered read-outs: [popinet3D_poly_sigma0_clipGate.yaml L52-L70](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet3D_poly_sigma0_clipGate.yaml#L52-L70) (FALSIFIED if the candidate fails at any step: the exemption returns the growth the bound removed), [popinet2D_clipRegionGate.yaml L60-L70](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet2D_clipRegionGate.yaml#L60-L70), [popinet2D_coneBoundGate.yaml L41-L60](https://github.com/leia-openfoam/leia/blob/8867581/config/popinet2D_coneBoundGate.yaml#L41-L60) (FALSIFIED if an interface metric moves by more than 1 %), [advConv2Dvortex.yaml L25-L36](https://github.com/leia-openfoam/leia/blob/8867581/config/advConv2Dvortex.yaml#L25-L36) (a bound arm passes only if its order is within 0.3 of `none`).

Why no formulation does both: a monotone bound cannot represent an extremum, and the signed distance has genuine extrema on its medial axis; a growing checkerboard has extrema at its own peaks, so any exemption of "a cell that is its stencil's extremum" hands the defect the cells it needs ([STATUS L2340-L2346](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2340-L2346), [STATUS L2417-L2421](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2417-L2421)). A nonnegative, partition-of-unity, quadratic-exact rule is impossible at an off-node point ([METHOD 9 L781-L785](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L781-L785)).

## What it does not cover

1. The polyhedral far-field failure stays. The production scheme has no maximum principle and `rho(B) > 1` on every mesh ([METHOD 8.2 L404-L420](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L404-L420)); the mesh-family decision keeps polyhedra off the production path ([[decisions/mesh-family-hexahedral]]).
2. The family stays as one study axis; the next candidate is the review's Rank 1, a monotone transporting base plus a filtered correction ([METHOD 8.3.6 L733-L737](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L733-L737)).
3. The flow-map Lipschitz constant `L(t) <= L(0) exp(int abs(grad u) dt)` is deliberately not implemented; it needs a collective per step ([METHOD 8.3.6 L742-L749](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L742-L749)).
4. Two threshold-free properties that separate a genuine extremum from a checkerboard peak (grad psi towards 0 at a medial-axis extremum; survival under a wider stencil) are proposed and untested ([STATUS L2362-L2374](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2362-L2374)).
5. The mesh-noise floor of the polyhedral rung is a property of the candidate: volume spread 0.0 / 0.0 / 14.2 % for none / stencilBounds / lipschitzCone ([CLAUDE L617-L623](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L617-L623), [[concepts/advection-regression-set]]).

## Related

[[hubs/advection]] - [[models/sl-value-bound]] - [[models/sl-reconstruction]] - [[concepts/value-bounds-and-clips]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/advection-regression-set]] - [[concepts/richardson-ladders-and-orders]] - [[decisions/sl-reconstruction-uncached-qwls]] - [[decisions/mesh-family-hexahedral]] - [[retractions/distance-cone-bound-as-transport-bound]] - [[retractions/clip-damage-is-the-narrow-band]] - [[retractions/advection-orders-3-2-factor]] - [[cases/popinet-translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED: the clip by G4 on 2026-09-09, the cone bound by the advection ladders on 2026-09-10. Entered in [[decision-log#2026-09]].
