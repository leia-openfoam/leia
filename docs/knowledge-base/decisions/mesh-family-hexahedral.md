---
title: "The mesh family: hexahedral production, polyhedral as the amplification rung"
description: "Hexahedral meshes are the production family, decided 2026-09-09 by the amplification probes: Lambda_max 1.0527 on blockMesh, 1.0549 on cfMesh's own hexahedral mesher and 1.0566 on snappyHexMesh with layers, against 1.2608 on pMesh polyhedra, whose growth exponent rho(B) - 1 is 2.3x larger; the defect is polyhedral cells, not cfMesh and not boundary layers"
aliases: [mesh family hexahedral]
kind: decision
status: settled
part: verification
tags: [decision, part/verification]
date: 2026-09-28
date_settled: 2026-09-09
decided_by: [workflow/scripts/fit_amplification_probe.sh, workflow/scripts/transport_spectrum_probe.sh, config/popinet3D_poly_sigma0_dtSweep.yaml, "author decision 2026-09-09"]
code: [workflow/scripts/fit_amplification_probe.sh, workflow/scripts/sl_fit_amplification_census.py, workflow/scripts/transport_spectrum_probe.sh, applications/test/leiaTestTransportSpectrum]
sources: ["METHOD 8.1 row mesh family (L399)", "METHOD 8.2 (L404-L488)", "STATUS 4 (L2121-L2231, L2233-L2310)", "sl_fit_amplification.csv", "CLAUDE regression set (L583-L645)"]
---
# The mesh family: hexahedral production, polyhedral as the amplification rung

> Hexahedral meshes are the production family; polyhedral meshes stay as the rung of the regression set that shows the amplification defect ([METHOD 8.1 L399](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L399), [CLAUDE L583-L645](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L583-L645)). Decided on 2026-09-09 by the one-mesh-pass amplification probe ([STATUS L2233-L2310](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2233-L2310)): the amplification bound of the fit `Lambda_c = abs(1 - sum g_j) + sum abs(g_j)` reaches 1.0527 on blockMesh, 1.0549 on cfMesh's own hexahedral mesher `cartesianMesh` with the same `meshDict` and surface, and 1.0566 on snappyHexMesh with three wall layers, against 1.2608 on pMesh polyhedra ([sl_fit_amplification.csv](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/sl_fit_amplification.csv)). The three hexahedral meshers cluster within 0.004 of each other, so the defect is polyhedral cells, not cfMesh and not boundary layers. The spectral radius of the frozen-velocity transport operator confirms the ordering: `rho(B) - 1` is 4.41e-03 on production hexahedra and 1.03e-02 on production pMesh, 2.3x in the exponent and 9.0e3x in the amplification over 1563 steps ([METHOD 8.2 L404-L488](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L404-L488)). A smaller time step cannot repair the polyhedral failure: the failure time is fixed to 6.5 % across a factor four in the step ([STATUS L2184-L2231](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2184-L2231)).

## The question

Every polyhedral conclusion of the campaign rested on one mesher, pMesh, so "is it only pMesh that is so bad?" was untestable ([STATUS L2233-L2238](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2233-L2238)). The sigma = 0 passive control on the polyhedral Popinet 3D mesh grows a false zero set at step 508 while the velocity is an exact uniform stream ([STATUS L1822-L1827](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1822-L1827)). The bound `Lambda` depends on the stencil geometry alone, so one mesh pass per mesher answers the question without a coupled run ([STATUS L2240-L2244](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2240-L2244)); `rho(B)` from the power iteration of `leiaTestTransportSpectrum` decides, because `Lambda` is a screen and not a proxy ([METHOD 8.2 L421-L435](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L421-L435), [[concepts/polyhedral-fit-amplification]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| blockMesh hex, 524 288 cells | `Lambda` median / p99.9 / max; cells > 1.10; demoted; smallest pivot | 1.0164 / 1.0430 / 1.0527; 0; 0; 0.637 | MEASURED, [STATUS L2245-L2252](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2245-L2252) |
| cartesianMesh (cfMesh hex), 281 216 cells, same `meshDict`, surface and cell size as pMesh | the same | 1.0131 / 1.0499 / 1.0549; 0; 0; 0.629 | MEASURED, same |
| snappyHexMesh + 3 wall layers, 617 984 cells | the same | 1.0164 / 1.0430 / 1.0566; 0; 0; 0.410 | MEASURED, same |
| pMesh (cfMesh polyhedral), 674 493 cells | the same | 1.0205 / 1.0815 / 1.2608; 588; 47 622; 1.4e-05 | MEASURED, same |
| the excess `Lambda_max - 1` | ratio | 0.053 / 0.055 / 0.057 / 0.261: pMesh 4.6 to 4.9x worse | MEASURED, [STATUS L2254-L2255](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2254-L2255) |
| cell-size histogram | cells below 0.4 h | cartesianMesh none below 0.6 h; snappy with layers none below 0.4 h; pMesh 27 191, and its five worst amplifiers at 0.32 to 0.33 h | MEASURED, [STATUS L2262-L2268](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2262-L2268) |
| power iteration, frozen uniform stream, Courant matched to 3 %: blockMesh n = 64 / snappy channel / pMesh production | `rho(B)`; `rho - 1`; factor over 1563 steps | 1.00441 / 1.00600 / 1.01028; 4.41e-03 / 6.00e-03 / 1.03e-02; 9.7e2 / 1.2e4 / 8.7e6 | MEASURED, [METHOD 8.2 L409-L417](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L409-L417) |
| `rho - 1` on hexahedra at n = 32 / 48 / 64 | invariance under refinement at fixed Courant | 4.39e-03 / 4.41e-03 / 4.41e-03; growth per physical time scales as 1/h | MEASURED, [METHOD 8.2 L466-L474](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L466-L474) |
| sigma = 0 polyhedral control at dt, dt/2, dt/4 | failure time; growth rate | 4.877e-03 / 4.660e-03 / 4.558e-03 s (6.5 %); 132.78 / 128.59 / 126.05 1/s | MEASURED, [STATUS L2204-L2214](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2204-L2214) |

Pre-registered read-out of the dt sweep ([STATUS L2171-L2175](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2171-L2175)): the control fails at the same physical time (steps 508 / 1016 / 2032) if the bound analysis holds; the measured steps are 528 / 1009 / 1974. The probe's own trap: the first snappyHexMesh arm was vacuous (`minVol 1e-13` removed every layer, exit 0) and reproduced the background blockMesh to eight digits; the artefact must be checked, never the exit code ([STATUS L2288-L2296](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2288-L2296)).

## What it does not cover

1. Polyhedral meshes stay supported: the fit uses the cell-face-cell stencil there, the transport orders match hexahedra (2.95 / 3.28 on the 3D shear case), and the polyhedral 3D shear rung is part of the standing regression set with its measured mesh-noise floor ([METHOD 2.2 L123-L125](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L123-L125), [METHOD 8 L343-L347](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347), [[concepts/advection-regression-set]]).
2. Hexahedra are not stable in the `rho <= 1` sense either; they carry a 2.5x smaller exponent ([METHOD 8.2 L461-L465](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L461-L465)). Snapping costs 36 % in `rho - 1`, which `Lambda` does not see; a microfluidic demonstration must measure `rho` on its own mesh ([METHOD 8.2 L481-L488](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L481-L488)).
3. The polyhedral Popinet 3D results are void for a different reason, tilted wall faces ([[retractions/polyhedral-popinet-3d-mesh-defect]]).
4. `Lambda` must not rank resolutions of one mesher or score a snapped mesh against an aligned one ([METHOD 8.2 L421-L435](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L421-L435)).
5. The route to a polyhedral production mesh is a rank reduction of the fit keyed on `Lambda <= 1`, untested ([STATUS L2280-L2286](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2280-L2286), [[decisions/sl-fit-normal-equations]]).

## Related

[[hubs/verification]] - [[hubs/advection]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/advection-regression-set]] - [[concepts/value-bounds-and-clips]] - [[models/sl-reconstruction]] - [[models/sl-value-bound]] - [[decisions/sl-clip-and-value-bound-off]] - [[decisions/sl-fit-normal-equations]] - [[retractions/polyhedral-popinet-3d-mesh-defect]] - [[cases/popinet-translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the probes of 2026-09-08 and 2026-09-09. Entered in [[decision-log#2026-09]].
