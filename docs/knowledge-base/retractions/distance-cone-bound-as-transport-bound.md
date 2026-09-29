---
title: "The distance-cone (lipschitzCone) bound as a transport gain (2026-09-10)"
description: "RETRACTED 2026-09-10, the same day - the cone bound improved every interface metric by 25 to 68 % on one coupled case at one resolution; the advection ladder and the resolution ladder falsify it as a transport bound"
aliases: [lipschitzCone retraction, cone bound retraction]
kind: retraction
status: retracted
part: advection
tags: [retraction, part/advection]
date: 2026-09-28
code: [src/leiaLevelSet/semiLagrangian/lipschitzConeValueBound.C, workflow/scripts/value_bound_ladder_table.py, config/popinet2D_coneBoundGate.yaml, config/popinet2D_coneBoundLadderN128.yaml, config/advConv2Dvortex.yaml, config/advConv2Dtranslation.yaml]
sources: [METHOD 8.3, STATUS header line, CLAUDE mesh-convergence rule, DP SL_VALUE_BOUND]
---
# The distance-cone (lipschitzCone) bound as a transport gain (2026-09-10)

> RETRACTED 2026-09-10, the same day. The claim was "the distance-cone bound improves the interface and worsens the amplifier: every interface metric improves by 25 to 68 %", made from one coupled case at one resolution ([METHOD 8.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L488-L496), the gain table at [METHOD 8.3.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L520-L540)). The measurement: on pure advection with every arm on the identical mesh, the bound is 3.6x worse than no bound on uniform translation (3.41e-04 against 1.24e-03), 19x worse on the reversed vortex and 9x worse on 3D shear, and its `stencil` mode loses the phase (E_VOL_ALPHA_REL = 1.0000) ([METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L592-L597)). The converged vortex ladder puts the disadvantage at 8.6x at N = 32 and 189.7x at N = 256 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L657-L667)). The resolution ladder of the coupled case shows a floor: the eikonal error moves 7.119e-03 to 7.066e-03 from N = 64 to 128 (order 0.01), the spurious-current order falls from 0.24 to 0.03, and the centroid error reverses to +119.4 % worse than no bound ([METHOD 8.3.5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L698-L729)). Scope of the retraction: the general claim. The N = 64 coupled numbers stand as one rung. `SL_VALUE_BOUND` stays at the sentinel that resolves to `none` ([METHOD 8.3.6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L739-L742)). No data is void.

## The claim, and where it lived

- [METHOD 8.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L488-L496), `METHOD.md`: the section was rewritten in place. Its first paragraph now records the retraction.
- [STATUS header](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L6), `STATUS.md`: the "Last updated 2026-09-10" line carries the falsification and the retraction.
- [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L554-L581): the rule "Touch advection, run a mesh convergence study" was written from this case.
- [cases/default.parameter](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L119-L146), `cases/default.parameter`: the `lipschitzCone` comment describes the bound and its modes.
- The curated ladder table `docs/method-comparison/method-comparison-article/data/tables/value_bound_ladder_popinet2D.csv` ([blob](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/value_bound_ladder_popinet2D.csv)) holds the N = 64 and N = 128 coupled numbers.

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The bound is inert where L = 1 is exactly valid (grad u = 0). | Uniform translation, one rung: 2.3 to 3.6x worse than `none`. Converged ladder: better at N = 32 (ratio 0.85), then 3.6x, 4.3x and 1.9x worse at N = 64, 128 and 256. | MEASURED, [METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L601-L604), [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L641-L653) |
| The bound is safe in strained flow. | Vortex: 8.6x to 189.7x worse, order 2.09 against 0.79 at the finest pair; volume error 3.95e-05 against 1.302e-02 (factor 330). | MEASURED, [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L657-L667) |
| The bound does not worsen the amplifier. | Per-step growth minus one: `none` 1.456e-04; `lipschitzCone` unity cellOnly 3.48e-04, 2.92e-04, 2.24e-04 at amplitudes 0.1 h, 1 h, 10 h (1.5 to 2.4x more). | MEASURED, [METHOD 8.3.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L550-L567) |
| The coupled gain is an asymptotic result. | Eikonal error order 0.01 (a floor at about 7e-03); the unbounded run converges at order 1.10 and reaches that floor near N = 256; centroid error +119.4 % worse at N = 128. | MEASURED, [METHOD 8.3.5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L698-L729) |
| Why the bound fails in strained flow. | The true Lipschitz constant grows as L(t) <= L(0) exp(int ||grad u||_inf ds); the implementation clamps L = 1 and fights the physics. | DERIVED, [METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L606-L615) |

The one case where the bound looked good, Popinet's translating droplet, has a nearly uniform velocity. That is the single regime where L = 1 is defensible ([METHOD 8.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L492-L496)).

## What survives

1. The `slValueBound` family (`none`, `stencilBounds`, `lipschitzCone`) as one study axis, gated byte-inert on 8 arms over 1563 steps at np 4 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L280-L286), [METHOD 8.3.6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L731-L737)).
2. The `-mode growth` instrument, validated against the power iteration to 3.6 % of (rho - 1) ([METHOD 8.3.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L558-L560)).
3. The exact-field unit gate: 0 empty intervals, and at 96 apex cells the cone bound errs by 0 where the monotone clip errs by one cell width (0.0157) ([METHOD 8.3.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L510-L518)).
4. Every bound prevents the polyhedral divergence at step 198, and the cheaper `stencilBounds` gives the best volume error there (4.4e-03 against 1.0e-01 and 3.2e-01) ([METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L617-L623)).
5. A finding about the baseline: the unbounded translation saturates at N = 256, order -0.54, error rising from 7.95e-05 to 1.16e-04 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L646-L653)). It is open.

## Propagation (checklist, same commit)

Done:

- [x] `METHOD.md`: [8.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L488-L496) rewritten with the retraction first; [8.3.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L578-L584) marked as single-rung numbers.
- [x] `STATUS.md`: the [header line](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L6).
- [x] `CLAUDE.md`: the [mesh-convergence rule](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L554-L581) and the [advection regression set](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L583-L604).
- [x] `cases/default.parameter`: the [token comments](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L104-L158).
- [x] The curated tables: `value_bound_ladder_popinet2D.csv` and the two `advConv2D*_convergence.csv` (commit 32a0d9b, 2026-09-10). The orders in those two CSVs were then found 3/2 too high, see [[retractions/advection-orders-3-2-factor]].
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The 3D hexahedral and polyhedral rungs of the converged ladders were "RUNNING, not reported" ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L677-L678)). The record has no entry of their landing.
- [ ] The flow-map Lipschitz constant is recorded as a poor bet and was not built ([METHOD 8.3.6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L744-L749)).

## Related

Hubs: [[hubs/advection]], [[hubs/verification]]. Siblings: [[models/sl-value-bound]], [[concepts/value-bounds-and-clips]], [[concepts/richardson-ladders-and-orders]], [[concepts/advection-regression-set]], [[decisions/sl-clip-and-value-bound-off]], [[retractions/clip-damage-is-the-narrow-band]], [[retractions/advection-orders-3-2-factor]].

## Log

### 2026-09-28
Written from METHOD 8.3, the STATUS header line and the CLAUDE.md rule. Retracted 2026-09-10. Entered in [[retraction-log#2026-09]].
