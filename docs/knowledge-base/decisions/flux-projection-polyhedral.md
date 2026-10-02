---
title: "FLUX_PROJECTION helmholtz on polyhedral meshes"
description: "Settled 2026-10-02: the prescribed face flux U(x_f) & S_f, far from discretely divergence-free on cfMesh polyhedra (max |div| 2.3 to 4.0 1/s at every rung), is projected once at the start; two defects of correctFlux were fixed first; on the 3D shear to 4.8M cells no transport metric degrades by more than 1.1 %."
aliases: [FLUX_PROJECTION, fluxProjection, Helmholtz projection of the prescribed flux, correctFlux]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-10-02
date_settled: 2026-10-02
decided_by: [config/gates/fluxProjectionGate3DshearPolyClip.yaml, config/gates/fluxProjectionGate3DshearPoly.yaml, "author decision 2026-10-01"]
code: [src/leiaLevelSet/velocityModel/velocityModel.C, src/leiaLevelSet/velocityModel/fluxCorrection.C, cases/3Dshear_poly.parameter, workflow/scripts/flux_projection_gate.sh, workflow/scripts/flux_projection_gate_table.py]
sources: ["STATUS 11.25 (L5095-L5172)", "METHOD 8.1 row FLUX_PROJECTION (L414)"]
---
# FLUX_PROJECTION helmholtz on polyhedral meshes

> Settled 2026-10-02, on the author's instruction of 2026-10-01: on polyhedral meshes the face flux `U(x_f) & S_f` of the prescribed velocity is far from discretely divergence-free, and the semi-Lagrangian `projectedFlux` trace reconstructs its cell velocity from that flux, so the flux is projected once at the start ([STATUS 11.25](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5095-L5172)). The token `FLUX_PROJECTION` is `helmholtz` in the polyhedral parameter layer only; hex studies are byte-identical. Measured on the polyhedral 3D shear up to 4,832,366 cells: the shape error at T moves by at most 1.1 %, the T/2 volume error and the band gradient error do not move, the T volume residual moves toward zero ([STATUS 11.25 item 5](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5134-L5159)). The projection does not cure the failure of `SL_CLIP false` on polyhedra ([[decisions/sl-clip-and-value-bound-off]]).

## The question

How large is the divergence of the prescribed flux on polyhedra, and does removing it change the transport? On uniform hex meshes the face-centre flux of the 3D shear field is divergence-free to round-off, because the three face differences of a hex cell cancel exactly; on cfMesh polyhedra they do not ([STATUS 11.25 item 1](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5103-L5109)).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| prescribed flux, cfMesh 3D shear, 105k to 4.8M cells | max\|div(phi)\| | 3.97 / 3.98 / 3.99 / 4.00 1/s (does not fall with refinement); hex: 1.8e-14 to 2.7e-13 | MEASURED, [STATUS 11.25 item 1](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5103-L5109) |
| `correctFlux` before the fixes | divergence left in the pinned cell | 5 to 56 1/s: the net boundary flux (up to 4.8e-05 m3/s) of an incompatible Neumann problem; and one pinned cell per rank in parallel | MEASURED, [STATUS 11.25 item 2](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5110-L5124) |
| `correctFlux` after the fixes | max\|div(phi)\|; serial against np 4 | <= 2.8e-09 1/s; 1.0e-15 m3/s (2.3e-12 of the largest face flux) | MEASURED, same |
| hex, `fluxProjection none`, old against new build | every CSV column and field | byte-identical | MEASURED, [STATUS 11.25 item 4](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5131-L5133) |
| published polyhedral configuration, 4 rungs | shape error at T; its orders | +0.1 / +0.3 / -0.5 / +1.1 %; 2.59 / 2.97 / 2.93 against 2.59 / 2.98 / 2.89 | MEASURED, [STATUS 11.25 item 5](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5134-L5159) |
| the same | volume error at T/2; band gradient error at T/2 | unchanged to four and three digits | MEASURED, same |
| the same | volume error at T | -2.4 / -6.7 / -2.4 / -78 %; the signed value changes sign across the rungs, so the -78 % is an absolute shift of 2.6e-05 | MEASURED, same |
| production configuration (`SL_CLIP false`), 3 rungs | volume error at T | 7.6 / 20.1 / 25.8 without, +0.1 to +9.4 % with the projection: no rescue | MEASURED, same |

## What it does not cover

1. 3Ddeformation on polyhedra is not measured with the projection ([STATUS 11.25 item 8](https://github.com/leia-openfoam/leia/blob/2793c9a7/STATUS.md#L5170-L5172)).
2. The curated polyhedral tables were made without the projection; they stay valid records of that configuration, and a regeneration includes it.
3. `benchVortexSLimprovedPerturbed` ran the old `correctFlux` (`-fluxCorrection`, np 4) with one pinned cell per rank; it feeds no curated table.

## Related

[[models/level-set-advection]] - [[decisions/sl-trace-velocity-projected-flux]] - [[decisions/sl-clip-and-value-bound-off]] - [[decisions/mesh-family-hexahedral]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/seam-checks-and-decomposition-invariance]] - [[decision-log]]

## Log

### 2026-10-02
SETTLED on the cluster gate of STATUS 11.25. Entered in [[decision-log#2026-10]].
