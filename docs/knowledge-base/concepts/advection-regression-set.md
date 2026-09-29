---
title: "The standing advection regression set"
description: "Any change to the advection path re-runs three coarse rungs (hex 2D 2Dvortex or 2Dtranslation, hex 3D 3Dshear, polyhedral 3D 3Dshear) at three resolutions against the preserved pre-change study; an inert change must be byte-identical, and a polyhedral rung reads its differences against a measured mesh-noise floor (2026-09-10)"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-29
date_settled: 2026-09-10
decided_by: [author decision 2026-09-10, config/advConv2Dtranslation.yaml, config/advConv2Dvortex.yaml, config/advConv3DshearHex.yaml, config/advConv3DshearPoly.yaml]
code: [config/advConv2Dtranslation.yaml, config/advConv2Dvortex.yaml, config/advConv3DshearHex.yaml, config/advConv3DshearPoly.yaml, workflow/scripts/advect_bound_arm.sh, workflow/scripts/compare_metrics_csv.py, workflow/scripts/advection_convergence_table.py]
sources: [CLAUDE regression set section, CLAUDE mesh convergence section, METHOD 8.3.4, METHOD 8.3.7, METHOD 8.3.8, STATUS 10.3, STATUS 11.7, STATUS 11.15]
---
# The standing advection regression set

> Verdict (2026-09-29). Any change to the advection path of the level set re-runs the standing regression set before any conclusion is drawn from any other study ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L583-L598)). The set has three rungs: hex 2D (`2Dvortex` or `2Dtranslation`), hex 3D (`3Dshear`) and polyhedral 3D (`3Dshear`). Each rung runs the default configuration at three resolutions against the preserved pre-change study. An inert change must be byte-identical there; a change that is not inert must show its order. The rule exists because a coupled two-phase gate does not exercise the transport the same way. The `slValueBound` refactor was bit-identical on the coupled solver over 1563 steps and 8 arms. Nobody ran the advection solver, although it reaches the same rewritten function (2026-09-10, [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L600-L604)). The regression then passed: all four cases were bit-identical in every physical column at every step at tolerance 0 ([METHOD 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L682-L692)). A polyhedral rung needs no bit-identical meshes, but it needs a measured mesh-noise floor. On 3D shear, `none` and `stencilBounds` moved 0.0 % between two independent meshes, and `lipschitzCone` moved its volume error by 14.2 % ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L606-L630)). A spread of that size says that the candidate is unstable to round-off.

## What it is

The rungs ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L590-L594)) and their committed studies:

| rung | case | mesh | why it is in the set | committed study |
|---|---|---|---|---|
| hex 2D | `2Dvortex` or `2Dtranslation` | `hex` | the cheapest transport order; `2Dtranslation` is the only O(1)-displacement gate | [`advConv2Dvortex.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv2Dvortex.yaml#L40-L66), [`advConv2Dtranslation.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv2Dtranslation.yaml#L40-L66): N = 32, 64, 128, 256; np 8 |
| hex 3D | `3Dshear` | `hex` | the 3D stencils and the cross terms of the quadratic fit | [`advConv3DshearHex.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearHex.yaml#L40-L66): N = 32, 64, 128; np 48 |
| polyhedral 3D | `3Dshear` | `poly` | the only rung where the polyhedral amplification defect appears | [`advConv3DshearPoly.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearPoly.yaml#L58-L84): `MAX_CELL_SIZE` 0.04, 0.02, 0.01; np 48 |

The procedure:

1. Preserve the pre-change study, or build the pre-change binaries into a separate prefix ([METHOD 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L684-L686)).
2. Run the default configuration at three resolutions on each rung with `leiaSemiLagrangeLevelSetFoam`. The velocity is prescribed and no force is in the loop, so a failure is a transport failure ([`advect_bound_arm.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/advect_bound_arm.sh#L15-L18)).
3. For a change that claims to be inert, run the pre-change and the new binary on the same mesh and the same `0/` ([METHOD 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L684-L686)). Compare every metric CSV with `compare_metrics_csv.py --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`. PASS means byte-identical ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L596-L598), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L637-L644)).
4. For a change that is not inert, report the error and the observed order of every metric, with the resolution range ([[concepts/richardson-ladders-and-orders]]). The committed studies pre-register the read-out. The order of the `none` arm in `E_GEOM_ALPHA_REL` lies within 0.3 of the recorded order at every rung pair. `E_BOUND_ALPHA` stays at round-off ([`advConv2Dtranslation.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv2Dtranslation.yaml#L25-L35)).
5. On the polyhedral rung, measure the mesh-noise floor with a null arm: the same configuration on two independently built meshes. Read every difference against that floor ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L615-L630)).

The mesh-noise floor. The workflow's `mesh` rule runs once per arm, and cfMesh is not bit-reproducible. Two runs of `pMesh` from the identical rendered case gave the same point count. But 2658 of 206528 points (1.3 %) differed in the last digit, at about 3e-12 relative, and some faces had another order ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L606-L613)). The two mesh files had different digests and sizes, 8277379 against 8277373 bytes ([`advConv3DshearPoly.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearPoly.yaml#L37-L48)). The null arm on the 3D shear polyhedral rung ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L617-L623)):

| arm | geometric error spread | volume error spread |
|---|---|---|
| `valueBound none` | 0.0 % | 0.0 % |
| `valueBound stencilBounds` | 0.0 % | 0.0 % |
| `valueBound lipschitzCone` | 0.3 % | 14.2 % |

So the mesh noise is not a general confound. It is a diagnostic. Two of the three arms are insensitive to it at four significant figures. The third turns a 3e-12 coordinate change into a 14 % volume difference ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L625-L630)).

One shared mesh across arms is the exception, for an effect that is expected to sit at the floor ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L632-L635)). The workflow builds the mesh once (`--until preprocess` on a one-arm config). `advect_bound_arm.sh` copies `constant/polyMesh` into every arm, regenerates `0/` from `0.org` and `leiaSetFields`, and prints the md5 of the mesh points. Arms are comparable only when the digests match ([`advect_bound_arm.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/advect_bound_arm.sh#L4-L13), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/advect_bound_arm.sh#L52-L62)).

The same diagnostic applies to the decomposition. The soft-wall source S1 differed by 76 % in `l2MagUPrime` between serial and np 4, where the baseline differed by 3.5e-7 ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L4010-L4014), [[concepts/seam-checks-and-decomposition-invariance]]).

## Why it matters

- A coupled gate cannot stand in for the advection gate. The function that the `slValueBound` refactor rewrote, `robustEvaluate`, is reached from `pointValueScheme` in both solvers ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L600-L604)).
- The method gates render every arm on a hexahedral mesh ([`render_gate_configs.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/render_gate_configs.py#L182)). The polyhedral rung therefore exists only in this set (DERIVED).
- A single rung gave a false gain. The distance-cone bound lowered every interface metric by 25 to 68 % at N = 64 on the coupled 2D case. At N = 128 the eikonal error moved from 7.119e-03 to 7.066e-03, an order of 0.01 ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L562-L569), [[retractions/distance-cone-bound-as-transport-bound]]).
- The polyhedral rung shows the amplification defect. The unbounded arm stopped with a floating-point exception at step 198 on 49 911 cells, and every bound prevented that ([METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L619-L625), [[concepts/polyhedral-fit-amplification]]).
- The order column needs its own check. The 2D orders of METHOD 8.3.7 were first published 3/2 too high, because the table script used h_eff = nCells^(-1/3) for a 2D case ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L629-L636), [[retractions/advection-orders-3-2-factor]]).

## Where in the code

- The four committed studies, each with its pre-registered read-out in the header and three `SL_VALUE_BOUND` arms: [`advConv2Dtranslation.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv2Dtranslation.yaml#L1-L39), [`advConv2Dvortex.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv2Dvortex.yaml#L1-L39), [`advConv3DshearHex.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearHex.yaml#L1-L39), [`advConv3DshearPoly.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearPoly.yaml#L1-L57).
- `workflow/scripts/advect_bound_arm.sh`: one arm on a given shared mesh ([L1-L62](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/advect_bound_arm.sh#L1-L62)), first commit [05a5447](https://github.com/leia-openfoam/leia/commit/05a5447) (2026-09-10).
- `workflow/scripts/compare_metrics_csv.py`: the comparator ([[concepts/bit-identity-and-inertness-gates]]).
- `workflow/scripts/advection_convergence_table.py`, `value_bound_ladder_table.py` and `value_bound_advection_census.py`: the read-out scripts that the config headers name.

## Evidence

| claim | number | where |
|---|---|---|
| the coupled gate did not cover the advection solver | 8 arms over 1563 steps, coupled solver only | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L600-L604), MEASURED |
| the regression passed | hex 2D translation, hex 2D vortex, hex 3D shear, polyhedral 3D shear: bit-identical in every physical column at every step, HEAD against 83db548 on an identical mesh and `0/`; the polyhedral rung stops at step 198 in both | [METHOD 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L682-L688), MEASURED |
| a whole-file `cmp` gives a false difference | DIFFERS on all four cases, from the two clock columns | [METHOD 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L690-L692), MEASURED |
| two `pMesh` runs of one case differ | 2658 of 206528 points in the last digit, about 3e-12 relative; 8277379 against 8277373 bytes | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L608-L613), [`advConv3DshearPoly.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/advConv3DshearPoly.yaml#L39-L43), MEASURED |
| the floor depends on the candidate | geometric and volume spread: `none` 0.0 % and 0.0 %, `stencilBounds` 0.0 % and 0.0 %, `lipschitzCone` 0.3 % and 14.2 % | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L617-L623), MEASURED |
| the unbounded polyhedral arm diverges | floating-point exception at step 198; geometric error 39.3 at step 186 | [METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L619-L621), MEASURED |
| the set in the library-split gate | 42 cases and 126 CSV pairs (the sum of the four rows), 0 differences; one mesh per polyhedral resolution, 49 911, 347 073 and 2 389 233 cells | [STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3084-L3093), MEASURED |
| the set in the Phase C gate | 15 cases at laptop sizes (the sum of the three rows), identical; polyhedral meshes of 18 082, 40 001 and 74 234 cells, shared | [STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3379-L3387), MEASURED |

## Decisions

- 2026-09-10: the standing regression set, the mesh-noise floor as a property of the candidate, and the committed comparator ([30cb989](https://github.com/leia-openfoam/leia/commit/30cb989), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L583-L645)).
- 2026-09-10: `SL_VALUE_BOUND` stays at the sentinel that resolves to `none` ([METHOD 8.3.6](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L741-L744), [[decisions/sl-clip-and-value-bound-off]]).

## Open questions

1. The unbounded `2Dtranslation` arm saturates: its order goes 3.80, 2.10, then -0.54 at N = 256. METHOD says that this needs its own investigation before the finest rung scores anything ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L641-L655)).
2. The curated `advConv2D*_convergence.csv` tables are to be regenerated on Lichtenberg with the corrected order script; STATUS records this as OPEN ([STATUS 11.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3255-L3260)).
3. METHOD 8.3.7 records the 3D hex and polyhedral ladders of 2026-09-10 as running and reports no order for them ([METHOD](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L679-L680)). The gates of 2026-09-23 and 2026-09-26 compared these rungs for bit identity ([STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3084-L3093), [STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3379-L3387)).

## Related

- Hubs: [[hubs/verification]], [[hubs/advection]].
- Siblings: [[concepts/richardson-ladders-and-orders]], [[concepts/value-bounds-and-clips]], [[concepts/bit-identity-and-inertness-gates]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/method-gates]].
- Models and decisions: [[models/sl-value-bound]], [[decisions/sl-clip-and-value-bound-off]], [[decisions/mesh-family-hexahedral]].
- Retractions: [[retractions/distance-cone-bound-as-transport-bound]], [[retractions/advection-orders-3-2-factor]].
- Cases and mechanism: [[cases/kinematic-advection-cases]], [[concepts/polyhedral-fit-amplification]].

## Log

### 2026-09-29
Created.
