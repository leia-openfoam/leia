---
title: "Popinet's translating droplet benchmark"
description: "The literature anchor of the translating droplet (Popinet 2009, Sec. 6.2.2) at density and viscosity ratio 1: the 2D ladder reproduces his two convergence statements, in SI units since 2026-09-08; the polyhedral 3D results on the tilted-face mesh are void, and no polyhedral run completes the horizon without the clip"
aliases: [popinetTranslating2D, popinetTranslating3D, Popinet benchmark, Popinet 2009 translating droplet]
kind: case
status: settled
part: mass-flux
tags: [case, part/mass-flux]
date: 2026-09-29
date_settled: 2026-09-08
decided_by: [config/popinet2D_La12000_N64.yaml, config/popinet2D_La12000_N128.yaml, config/popinet2D_La12000_N256.yaml, config/popinet2D_LaSweep_N64.yaml]
code: [cases/popinetTranslating2D, cases/popinetTranslating2D.parameter, cases/popinetTranslating3D, cases/popinetTranslating3D.parameter, cases/popinetTranslating3D_poly.parameter, workflow/scripts/make_popinet_table.py, workflow/scripts/popinet_si.py]
sources: [STATUS 0 ANCHORED, STATUS 4 Popinet-3D polyhedral (2026-09-05 to 2026-09-09), STATUS 4 SI units, SL article sec:popinet, METHOD 8.1 rows CURVATURE_EXTENSION and SL_CLIP, popinetTranslating.csv]
---
# Popinet's translating droplet benchmark

> Verdict (2026-09-29). Popinet introduced the translating droplet for the coupling of surface tension with the advection of a circular interface (JCP 228 (2009) 5838-5866, Sec. 6.2.2, DOI [10.1016/j.jcp.2009.04.042](https://doi.org/10.1016/j.jcp.2009.04.042); [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L135-L140)). The 2D reproduction ran on 2026-09-04 and again in SI units on 2026-09-08 at his dimensionless groups ([STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L186-L193)). The record uses his convention: the maximum over time of the spatial norm, relative to U. The L2 velocity error is 4.43e-3, 2.49e-3 and 1.26e-3 at N = 64, 128, 256, order 0.91 ([SL article, sec:popinet](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2382-L2404), `sec:popinet`). The shape error converges at order 1.71 ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2406-L2412)). Popinet states his maximum error only in L_inf ("of the order of 5 % of U"). The record compares with that number only in L_inf, and this note does not use that comparison as a result. The half-order maximum against the near-first-order L2 reproduces his own statements, and it is the basis of the L2/L1 rule ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2414-L2418), [[concepts/error-vector-and-read-out-instants]]). Every polyhedral 3D result on the first cfMesh surface is VOID: 4.8 % of the wall faces were tilted into the flow ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1756-L1786), [[retractions/polyhedral-popinet-3d-mesh-defect]]). On the corrected meshes the level-set transport alone fails in the far field, at step 508 with sigma = 0 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1840-L1848)). The quasi-monotone clip removes that failure but is not in the production method ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L380)).

## What it is

| item | Popinet 2009 | this repository (SI, since 2026-09-08) |
|---|---|---|
| droplet | D = 0.4 in the unit square | D = 1 mm, 0.4 of the box height 2.5 mm |
| box, boundaries | periodic in x, symmetry top and bottom | 5 mm x 2.5 mm, inlet left and outlet right, slip top and bottom |
| properties | rho and nu constant: both ratios 1 | rho = 1000 kg/m^3, nu = 1e-6 m^2/s in both phases |
| groups | We = 0.4; La = 120, 1200, 12 000, infinity | sigma = La rho nu^2/D = 0.012 N/m; U = sqrt(We sigma/(rho D)) = 0.0693 m/s |
| horizon | one t/T_U | T_U = D/U = 14.4 ms |
| method | adaptive quadtree VOF, height-function curvature | uniform-mesh semi-Lagrangian level set, least-squares curvature |

Sources: [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L142-L148), [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2345-L2376), [`popinetTranslating2D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/popinetTranslating2D.parameter#L1-L45). The face properties of Popinet are the `alg_lin` mixture, which the viscosity ladder of this repository selected independently ([[models/viscosity-face-model]]).

Two differences stay, and the record states them. First, the SL transport has no periodic support, so inflow and outflow replace periodicity. The box is 2 x 1. After one t/T_U the droplet is still 57 cells clear of the outlet at N = 64. Second, the discretisations differ, and that difference is the object of the comparison ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2369-L2376), [config header](https://github.com/leia-openfoam/leia/blob/d1e3414/config/popinet2D_La12000_N64.yaml#L38-L51)).

The 2D case has the patches `inlet`, `outlet` and `walls`, the same conditions as [[cases/translating-droplet]], and a box length of `POPINET_XLEN` box heights ([`blockMeshDict.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/popinetTranslating2D/system/blockMeshDict.template#L36-L41)). The 3D case is 5 x 2.5 x 2.5 mm at N = 64, 96, 128 ([`popinetTranslating3D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/popinetTranslating3D.parameter#L1-L32)). Its polyhedral twin uses the feature-edge surface `box5x2p5x2p5mm.fms` ([`popinetTranslating3D_poly.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/popinetTranslating3D_poly.parameter#L31-L45)). The time step is 0.2323 of the Brackbill limit, as in every production run ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/popinetTranslating2D.parameter#L24-L29)).

`CURVATURE_EXTENSION none`: every Popinet config sets it in its `axes_override` ([`popinet2D_La12000_N64.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/popinet2D_La12000_N64.yaml#L87)). It is the setting of the reproduction, not a measured preference; no case layer sets the token ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393), [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3591-L3592), [[models/curvature-extension]]).

## Why it matters

It places the SL two-phase method against an independent method on an identical problem ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2338-L2343)). At density ratio 1 it isolates the curvature-error source from the density-ratio amplifier of [[cases/translating-droplet]] ([[concepts/density-ratio-amplifier]]). Its velocity is nearly uniform, so it is a weak test of a transport bound: the distance-cone bound looked good only here ([METHOD 8.3](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L492-L498), [[retractions/distance-cone-bound-as-transport-bound]]). The 3D polyhedral twin was the first translating polyhedral case, and it found the fit defect that no kinematic gate had seen ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L516-L526)).

## Where in the code

- Cases and configs: `cases/popinetTranslating{2D,3D}`, the three `.parameter` files, and `config/popinet2D_*`, `config/popinet3D_*`.
- The SI set: `workflow/scripts/popinet_si.py` prints the set and checks a config against it ([`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L995-L1007)).
- The curation: `workflow/scripts/make_popinet_table.py` writes [`popinetTranslating.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/popinetTranslating.csv). It keeps L_inf only because Popinet tabulates it ([`make_popinet_table.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_popinet_table.py#L11-L17)).

## Evidence

| claim | number | where |
|---|---|---|
| the 2D ladder in L2 (max over time, relative to U) | 4.43e-3, 2.49e-3, 1.26e-3 at N = 64, 128, 256; order 0.91 (R = 0.999) | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2393-L2402), [CSV](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/popinetTranslating.csv#L2-L4) |
| the 2D ladder in L1 | 1.86e-3, 1.07e-3, 4.75e-4; the record states no L1 order; the least-squares slope of these values is 0.98 | MEASURED, [CSV](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/popinetTranslating.csv#L2-L4); the order is DERIVED (2026-09-29) |
| the shape error, relative to R | 4.08e-3, 1.04e-3, 3.78e-4; order 1.71 (R = 0.996) against his "roughly first order" | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2406-L2412) |
| the L2 comparison with Popinet | 4.43e-3 at N = 64 against a peak of about 3e-3 in his figure | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2406-L2408) |
| the comparison with his quoted maximum | L_inf only (his convention): 4.86 % of U against "of the order of 5 %"; order 0.49 against his "less than first order". Not a result of this note. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L159-L167) |
| the Laplace sweep, N = 64 | the article gives L_inf only (2.93e-2 to 5.88e-2). The CSV holds L2 2.50e-3, 3.49e-3, 4.43e-3, 5.40e-3 and L1 1.22e-3, 1.58e-3, 1.86e-3, 2.18e-3 at La = 120, 1200, 12 000, infinity: the inviscid case is the worst | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2420-L2427), [CSV](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/popinetTranslating.csv#L5-L8) |
| the SI conversion is exact | 28 curated entries within 0.051 %; the step counts match (1563, 4420, 12 501); the N = 64 similarity gate differs by at most 0.067 % | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L186-L193), [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2081-L2098) |
| the 3D hexahedral twin completes the horizon (N = 64, his units) | L1 7.4e-4, L2 2.0e-3 (2D: 1.9e-3, 4.4e-3); final volume 1.4e-3, shape 3.0e-4, Laplace jump 10.055 against 10 | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1928-L1933) |
| the first polyhedral smoke | the far field 1.4 h wrong after one step, divergence at step 3; stencil condition numbers 3e7 to 6e12 | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1544-L1597) |
| the tilted wall faces | 1 984 of 40 960 wall faces (4.8 %) tilted up to 8 degrees; `simpleFoam` 8.3 % off the uniform stream; the feature-edge surface gives 0 tilted faces and 5.7e-16 | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1756-L1775) |
| the corrected-mesh ladder diverges late | N = 64 and 96 at steps 1154 and 1308 (one sliver cell on a box edge); the 0.0195 mesh at step 1026, fake zero set at step 508 | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1801-L1839) |
| the transport alone fails | sigma = 0: the velocity stays uniform (L2 3.4e-15), the fake zero set appears at step 508, the volume error reaches 379 % | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1840-L1848) |
| the boundary-face stencil defect (survives the void) | normal transport fraction in the first layer: inlet 0.082 polyhedral and 0.198 hexahedral, outlet 0.240 and 0.491 | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1892-L1915) |
| `stencilBoundaryFaces inflowOnly` in 2D | outlet fraction 0.962 to 1.000; every droplet metric identical to `include` | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1942-L1948) |
| the clip removes the far-field failure | sigma = 0, 1563 steps: no fake zero set, volume error 1.08e-6 against 379 % | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2026-L2044) |
| the clip is not inert on the 2D interface | N = 64: volume error +30.4 %, centroid +12.5 %, shape +4.9 % | MEASURED, [`popinet2D_clipStencilGate.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/popinet2D_clipStencilGate.yaml#L4-L7), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L380) |

## Why it failed, or why we think so

The 2D benchmark did not fail. The polyhedral 3D benchmark failed three times, each time for a different reason:

1. The quadratic fit on cfMesh's boundary-layer cells was ill-conditioned ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1568-L1602)). The admissibility test fixed it, with the default `quadraticPivotTol 0.3` ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1698), [[models/sl-reconstruction]]).
2. The STL surface put four side walls into one solid, so cfMesh tilted wall faces into the flow. The mesh itself did not carry a uniform stream. VOID ([[retractions/polyhedral-popinet-3d-mesh-defect]]).
3. On the corrected meshes the SL update amplifies a checkerboard error in small, one-sided boundary cells near the outlet. The record reads the rate as about `(U dt/h)` times the stencil asymmetry per step, 0.5 to 1 % per step ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1848-L1852)). The amplification bound and the spectral radius of the fit measure it ([[concepts/polyhedral-fit-amplification]], [[decisions/mesh-family-hexahedral]]).

## Decisions

- The Popinet cases run in SI units at his groups (2026-09-08); the conversion is measured, not assumed ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2059-L2075)).
- L_inf stays in the Popinet table only because Popinet tabulates it; it carries no verdict ([`make_popinet_table.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_popinet_table.py#L15-L17)).
- `stencilBoundaryFaces` stays `include` by default; `inflowOnly` is selectable ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1942-L1948)).
- `SL_CLIP false`: the clip is not in the best configuration ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L380), [[decisions/sl-clip-and-value-bound-off]], [[retractions/clip-damage-is-the-narrow-band]]).

## Open questions

1. No polyhedral run of this benchmark completes the horizon with the production method. The discriminators that remain are a different polyhedral mesh generator and an update that cannot amplify a checkerboard ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1999-L2004)).
2. The subsection "RESULT: the Popinet-3D polyhedral ladder DIVERGES late" has no VOID marker, although its rungs ran on the tilted mesh ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1864-L1890)).
3. The SL article sends the reader to `sec:popinet` for the polyhedral divergence at step three, but `sec:popinet` covers only the 2D benchmark ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L438-L441)).
4. The first STATUS table gives an L2 order of 0.88, the SI table 0.91, and STATUS calls the orders "unchanged" ([STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L157), [L191](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L191)). DERIVED (2026-09-29): the rounded values 4.4e-3, 2.5e-3, 1.3e-3 give 0.88; the CSV values give 0.91. The 0.91 is the right number.
5. The full-horizon 3D hexahedral twin ran in Popinet's units only; the SI check of 3D covers 78 steps ([STATUS 4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2117-L2123)).

## Related

- Hubs: [[hubs/mass-flux]], [[hubs/verification]].
- Cases: [[cases/translating-droplet]], [[cases/benchmark-cases]].
- Retractions: [[retractions/polyhedral-popinet-3d-mesh-defect]], [[retractions/distance-cone-bound-as-transport-bound]], [[retractions/clip-damage-is-the-narrow-band]].
- Concepts: [[concepts/error-vector-and-read-out-instants]], [[concepts/polyhedral-fit-amplification]], [[concepts/density-ratio-amplifier]], [[concepts/wrong-setup-voids]].
- Models and decisions: [[models/sl-reconstruction]], [[models/sl-value-bound]], [[models/curvature-extension]], [[models/viscosity-face-model]], [[decisions/mesh-family-hexahedral]], [[decisions/sl-clip-and-value-bound-off]].
- Study: [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-29
Created from STATUS 0 and 4, the SL article, the curated Popinet CSV, the case files and the configs.
