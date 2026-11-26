---
title: "Boundary points as advected unknowns: the wall treatment of the semi-Lagrangian fit"
description: "Decided 2026-11-26: on cubic cells, boundary-face, edge and corner values are advected unknowns and a wall cell uses the 27-point tensor interpolant; exact for a quadratic, no growth at walls in the symbol, and with a cell-value velocity gradient the rotation near the walls is second order in dt."
aliases: [wall tensor, cubicLatticeWallTensor, faceUnknownsTensorI2, boundary face unknowns, taylorCellGradient, advanceBoundary]
kind: concept
status: open
part: advection
tags: [concept, part/advection]
date: 2026-11-26
decided_by: [docs/analysis-amplification-sl/sympy/wall_face_unknowns.py, docs/analysis-amplification-sl/sympy/wall_edges_corners.py, config/slAmr/wallTensorExact3Dserial.yaml, config/slAmr/wallTensorRotationCellGradDt3D.yaml]
code: [src/leiaLevelSet/semiLagrangian/cubicLatticeWallTensorReconstruction.C, src/leiaLevelSet/semiLagrangian/foot/slTaylorCellGradientFootIntegrator.C, src/leiaLevelSet/semiLagrangian/slReconstruction.H]
sources: ["STATUS 11.45", "docs/plan-sl-amr-preprint.md section 6 items 8 and 9"]
---
# Boundary points as advected unknowns: the wall treatment of the semi-Lagrangian fit

> Verdict (2026-11-26, MEASURED). The production boundary data (`stencilBoundaryFaces inflowOnly`) make a wall cell fall back to a linear fit wherever a face is dropped, and that error stays in the domain when the flow is tangential to the wall: the rigid rotation has 9.3e-3 h in the wall cells and an interface error that does not fall with dt ([[concepts/sl-amplification-symbol]], STATUS 11.45 item 3). Boundary-face, edge and corner values as advected unknowns with the wall cell on the 27-point tensor interpolant are exact for a quadratic and have no growth at walls in the symbol. The class `cubicLatticeWallTensor` with the foot integrator `taylorCellGradient` gives second order in dt near the walls (E 2.59e-4, 6.48e-5, 1.62e-5 h). Supervisor's decision for the AMR line on cubic cells; the unstructured form is a roadmap item.

## What it is

- Each wall cell owns the psi value of its boundary points: one face centre per wall, one point on each boundary edge (at the height of the cell centre), the corner. The values are cell fields exchanged through the stencil map, so the treatment works across processor boundaries.
- A wall cell is evaluated with the tensor-product quadratic on the nodes {wall, 0, 1} along each wall axis and {-1, 0, 1} along the others. Every other cell keeps its fit; an interior cell reads no face data.
- `advanceBoundary` (a hook of `slReconstruction`, inert by default) moves every boundary point from its own departure point in the interpolant of psi^n. A departure point outside the domain takes the patch value at the new time; an inflow edge or corner takes the quadratic in the wall plane through patch and cell values.
- The velocity at a boundary point and the dt^2 term of the cell foot must come from cell values: the Gauss gradient of `reconstruct(phi)` reads the adjacent cell value on a wall face.

## Evidence

| claim | number | where |
|---|---|---|
| production inflowOnly is inexact and grows at the inflow wall | linear fallback in 32,256 wall cells; symbol growth 1.2e-3 to 3.1e-3 at every slab thickness | STATUS 11.45 items 3 and 6 |
| face unknowns with a least-squares wall fit amplify | 4.0e-2 per step at the outflow wall, mode (pi, pi) along the wall | `wall_face_unknowns.py` |
| tensor interpolant + inner cells without face data + inflow from data | no growth at Co >= 0.1 at walls, edges, corners; 2.6e-4 at Co 0.0164, falling with thickness | `wall_face_unknowns.py`, `wall_edges_corners.py` |
| exactness in the solver, sphere leaving through the outflow wall | max error 2.4e-9 h over all cells, serial = 4 ranks | `config/slAmr/wallTensorExact3Dserial.yaml` |
| rotation near the walls | E 2.59e-4 / 6.48e-5 / 1.62e-5 h, wall cells 2.81e-4 / 7.05e-5 / 1.77e-5 h at dt, dt/2, dt/4 | `config/slAmr/wallTensorRotationCellGradDt3D.yaml` |

## Open questions

- Stage B: face values across level transitions along a wall (AMR), and the wall rule of the refinement band.
- The unstructured formulation (plan section 6 item 9).
- A mesh ladder of the rotation near the walls; the corner power bound.
- A contact-angle condition as a normal-derivative row of the same construction.

## Related

[[concepts/sl-amplification-symbol]] · [[concepts/sl-dynamic-amr]] · [[concepts/sl-ghost-cell-reconstruction]] · [[models/sl-reconstruction]] · [[concepts/departure-foot-ab2-centring]]

## Log

### 2026-11-26
Created. Symbol analysis (three candidates, the third decided), the class, the exactness gate and the rotation gates of STATUS 11.45 items 6 and 7.
