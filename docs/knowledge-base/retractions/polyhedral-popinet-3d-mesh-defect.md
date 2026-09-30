---
title: "The polyhedral Popinet 3D results: tilted wall faces (2026-09-05)"
description: "VOIDED 2026-09-05 - every polyhedral Popinet-3D result on the STL mesh: four box edges were not feature edges, 4.8 % of the wall faces were tilted into the flow, and a uniform stream was not a discrete solution of the mesh"
aliases: [tilted wall faces void, popinet3D poly VOID]
kind: retraction
status: voided
part: advection
tags: [retraction, part/advection]
date: 2026-09-28
code: [cases/popinetTranslating3D_poly.parameter, config/popinet3D_La12000_poly_r12p8.yaml]
sources: [STATUS 4 (2026-09-05 to 2026-09-08), commit 0da60f3, sl_fit_pivot_census.csv]
---
# The polyhedral Popinet 3D results: tilted wall faces (2026-09-05)

> VOIDED 2026-09-05. The claims were the polyhedral Popinet-3D results on the STL mesh: the two ladder rungs r12p8 and r19p2 that diverged at steps 1151 and 1301 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1846-L1864)), the 78-step smokes and their exports, the field dumps, the translating arms of the corrector sweep, the census row "poly uniform Popinet-3D N = 64", and the two `inflowOnly` horizon attempts ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1758-L1762)). The measurement: OpenFOAM's `simpleFoam` on the same pMesh does not keep the uniform stream. After 40 iterations from the exact state |U - U0| is 8.3 % in the eight concave inlet-corner cells, 3 % in the wall layer and 1 % in the interior. 1 984 of the 40 960 wall faces (4.8 %) have normals tilted into the flow, up to 8 degrees at the box corners, because `box2x1x1.stl` put the four side walls into one solid and cfMesh's dual wrapped faces around the edges that were never feature edges. The feature-edge surface gives 0 tilted faces, and `simpleFoam` then holds the stream to |U - U0| <= 5.7e-16 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1740-L1756)). Scope of the void: every polyhedral Popinet-3D result on the STL mesh. The directories carry the suffix `_VOID_tiltedWallFaces_20260905`, and the committed `popinet3D_La12000_poly_smoke4_errors.csv` was removed (commit [0da60f3](https://github.com/leia-openfoam/leia/commit/0da60f3)). The stationary polyhedral ladder is not affected: its surface has one solid per plane and 0 tilted faces ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1750-L1751)).

## The claim, and where it lived

- [STATUS 4, RUNNING 2026-09-05](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1514-L1548), `STATUS.md`: the smoke that diverged at step 3, the ladder plan and the two controls.
- [STATUS 4, ROOT CAUSE](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1550-L1575): the ill-conditioned quadratic stencils on the boundary-layer polyhedra. This finding survives (below).
- [STATUS 4, RESULT 2026-09-05 evening](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1846-L1907): "the Popinet-3D polyhedral ladder DIVERGES late". Every number in that subsection is from the tilted mesh. The subsection was written before the void and carries no marker.
- The census table `docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_fit_pivot_census.csv` ([blob](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_fit_pivot_census.csv)): the row "poly uniform Popinet-3D N=64" was regenerated on the feature-edge mesh (commit 1583245).

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| A uniform stream is a discrete solution on the mesh. | `simpleFoam`: 8.3 % in the corner cells, 3 % in the wall layer, 1 % in the interior; pressure -0.19 to +0.22 where 0 is exact; a converged 302-iteration solve lands on the same field. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1740-L1745) |
| Every wall face is perpendicular to the flow. | 1 984 of 40 960 wall faces (4.8 %) tilted, up to 8 degrees. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1746-L1748) |
| The cause. | One solid in the STL, so no feature edges. An explicit `boundaryLayers` block does not cure it (8.6 %, p +-0.18). | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1748-L1752) |
| The cure. | `surfaceFeatureEdges -angle 45` (12 edges): 0 tilted faces on walls, inlet and outlet; |U - U0| <= 5.7e-16 and p within 1e-14. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1752-L1756) |
| The solver's constant boundary-cell velocity 0.0945 was a level-set defect. | It is the single-phase mesh defect (0.09451 constant from step 20 in r12p8). | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1743), [L1853](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1853) |

## What survives

Measured on hexahedra or on the stationary box as well ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1762-L1766)):

1. The quadratic-fit admissibility rule, pivot 0.3, bit-inert everywhere ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1550-L1575)).
2. The boundary-face effect on the semi-Lagrangian stencil: normal to an inflow or outflow boundary the level set moves at 8 to 50 % of the exact rate (inlet 0.082 polyhedral, 0.198 hexahedral; outlet 0.240 and 0.491) ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1874-L1891)). `stencilBoundaryFaces inflowOnly` makes the outlet exact in 2D ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1924-L1930)).
3. The non-orthogonal correction converges in one pass (stationary r13p8; 1, 3 and 6 correctors agree to 5e-6 relative) ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1941-L1957)).

What replaced the voided results, on the corrected meshes:

1. The 4-rank gate on the feature-edge mesh passed: max|u'| 8.6e-3 against 9.5e-2 on the tilted mesh and 1.3e-2 on hexahedra, now at the droplet surface ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1772-L1781)).
2. Both rungs diverged again, at steps 1154 and 1308, by a different route: one sliver cell of 0.153 h on a box edge ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1783-L1807)).
3. The 0.0195 mesh (no flagged defects) passed the gate and diverged at step 1026; the fake zero set appears at step 508, and the sigma = 0 control fails at the same step. The transport alone fails ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1809-L1830)).
4. Three pMesh variants, three seeds, one class: the semi-Lagrangian far-field transport grows an alternating mode in cfMesh's small, asymmetric boundary cells near the outlet ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1817-L1821)). That led to the clip work and the amplification bound, see [[retractions/clip-damage-is-the-narrow-band]] and [[concepts/polyhedral-fit-amplification]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: the [VOID paragraph](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1738-L1770).
- [x] The study directories renamed on the laptop and on the cluster; the committed smoke CSV removed (commit 0da60f3).
- [x] The census row regenerated on the feature-edge mesh (commit 1583245).
- [x] The case ships the feature-edge surface: `SURFACE_FILE ( box5x2p5x2p5mm.fms )` in [`cases/popinetTranslating3D_poly.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/popinetTranslating3D_poly.parameter#L41-L45) (now in SI units), and all polyhedral configs switched ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1768-L1769)).
- [x] The rule: "tilted wall faces" is in the list of wrong setups, see [[concepts/wrong-setup-voids]].
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The subsection [STATUS 4 L1846-L1907](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1846-L1907) reads as a result and has no VOID marker in its heading.
- [ ] The polyhedral Popinet benchmark has no completed horizon on any pMesh variant; the pre-registered verdict is "cannot be rescued by a cell size" ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1820-L1821)). A different boundary mesh is proposed and not tried ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1805-L1807)).

## Related

Hubs: [[hubs/advection]], [[hubs/verification]]. Siblings: [[cases/popinet-translating-droplet]], [[concepts/polyhedral-fit-amplification]], [[concepts/wrong-setup-voids]], [[models/sl-reconstruction]], [[decisions/mesh-family-hexahedral]], [[retractions/clip-damage-is-the-narrow-band]], [[retractions/closed-box-translating-droplet]].

## Log

### 2026-09-28
Written from STATUS 4 (2026-09-05 to 2026-09-08). Voided 2026-09-05. Entered in [[retraction-log#2026-09]].
