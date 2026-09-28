---
title: "volumeCorrection: noVolumeCorrection, newtonShift"
description: "The global uniform shift restores the total phase volume without touching the gradient, but it is a crossover: it helps where the volume error is large and lowers the order where the volume error is already small; production runs noVolumeCorrection."
aliases: [volumeCorrection, VOLUME_CORRECTION, newtonShift, globalShift]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-22
decided_by: [config/sdplsVolCorr3Dshear.yaml, config/sdplsVolCorr2DvortexRev.yaml]
code: [src/leiaLevelSet/volumeCorrection/volumeCorrection.H, src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.H]
sources: ["SDPLS article sec:volcorr (L2507-L2675)", "METHOD 8.1 row VOLUME_CORRECTION (L386)", "DP L196-L200", "STATUS 11.12 (L3566-L3569)"]
---
# volumeCorrection: noVolumeCorrection, newtonShift

> Verdict (2026-09-28). `newtonShift` adds one scalar `eps` to the whole level set, the root of `f(eps) = sum_c alpha_c(psi + eps) abs(Omega_c) - V_target`, found by Newton-Raphson with bisection ([SDPLS article L2510-L2520](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2507)). A constant shift leaves the discrete gradient unchanged bit for bit, so `abs(grad psi)` and both redistancing triggers are untouched ([SDPLS article L2521-L2532](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2507)). The measurement is a crossover: on 3D shear with the R source the shape order rises from +0.916 to +1.293, on the sourceless 3D shear the volume order falls from +3.928 to +1.050 and the shape order from +1.866 to +1.745 ([SDPLS article L2551-L2600](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2537)). Production runs `noVolumeCorrection` ([METHOD 8.1 L386](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L386)). It must never be used on the coupled droplet ([SDPLS article L2667-L2675](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2657)).

## What it is

The conserved functional is `V(psi) = sum_c alpha_c(psi) abs(Omega_c)`, evaluated through the runtime-selected phase indicator, never through an indicator of the corrector's own ([volumeCorrection.H L35-L54](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L35-L54)). The target is captured once from the alpha field at start-up so that it equals the `alphaV0` baseline of the error metric ([volumeCorrection.H L125-L130](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L125-L130)). The applied shift is kept so that a CSV column can report it: a volume error of 0 with a growing shift is a failing scheme with its symptom removed ([volumeCorrection.H L132-L138](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L132-L138)).

The key is `levelSet.volumeCorrection.type`, read with `subOrEmptyDict`, so a dictionary without the entry selects `noVolumeCorrection` and runs bit-identically ([volumeCorrection.H L62-L78](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L62-L78), [volumeCorrection.C L69](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.C#L69)). The token is `VOLUME_CORRECTION` ([DP L196-L200](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L196-L200)). Three solvers build it: `leiaLevelSetFoam`, `leiaLevelSetTwoPhaseFoam` and `leiaSemiLagrangianLevelSetTwoPhaseFoam` ([leiaLevelSetFoam.C L136](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetFoam/leiaLevelSetFoam.C#L136), [volumeCorrectionFields.H L32](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/volumeCorrectionFields.H#L32), [volumeCorrectionFieldsSL.H L47](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/volumeCorrectionFieldsSL.H#L47)).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| no correction | `noVolumeCorrection` | settled, production | The base class; `correct()` returns false and never touches psi. | [volumeCorrection.H L145](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L145), [METHOD 8.1 L386](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L386) |
| global uniform shift | `newtonShift` (alias `globalShift`) | settled, a crossover | Restores the total volume exactly; gradient untouched; ill-posed for several bodies; trades local shape for one global number. | [newtonShiftVolumeCorrection.H L30-L117](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.H#L30-L117), [alias L44-L57](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.C#L44-L57) |

The dictionary words are the `TypeName` strings ([noVolumeCorrection L145](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/volumeCorrection.H#L145), [newtonShift L298](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.H#L298)). The bracket of the bisection is `maxShiftCells` times `h` ([newtonShift H L58-L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.H#L58-L60)).

## Why it matters

Semi-Lagrangian point updates are non-conservative; the volume error is a bounded diagnostic and never telescopes to zero ([SL article L1182-L1186](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1164)). A correction that restores the volume number while the interface degrades hides a defect. The SDPLS article states the rule: any study that uses this model reports the shape error next to the volume error and the applied shift ([SDPLS article L2673-L2675](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2657)).

## Where in the code

- Family: `src/leiaLevelSet/volumeCorrection/`, library `libleiaVolumeCorrection`.
- The boundary check: `correct()` warns by name if any psi patch fixes its value, because a `fixedValue` patch does not shift with the interior ([newtonShift H L70-L74](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/volumeCorrection/newtonShiftVolumeCorrection.H#L70-L74)).

## Evidence

| claim | number | where |
|---|---|---|
| Gradient order untouched by the shift, all four pairs | 3D shear: -0.171 to -0.161 (noSource), +0.679 to +0.707 (R); 2D vortex: -0.116 to -0.107, +1.112 to +1.119 | MEASURED, [SDPLS article Table volcorr, L2551-L2567](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2537) |
| Where the volume error is large (3D shear, R, volume order +0.230) | shape order +0.916 to +1.293 | MEASURED, same table |
| Where the volume error is small (3D shear, noSource) | volume order +3.928 to +1.050; shape +1.866 to +1.745 | MEASURED, same table |
| 2D reversed vortex, R | shape order +1.883 to +1.671; the fitted volume value 0.835 is a floor, not a rate | MEASURED, same table and [L2642-L2655](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2642) |
| The floor: N = 32, t = 0.05 smoke | relative volume error 8.5917e-04 to 3.7357e-16 (simpleLinearImplicit), 2.6103e-03 to 2.4693e-13 (exponential), in two Newton iterations, `abs(eps)` 0.001 to 0.003 cell widths | MEASURED, [SDPLS article L2647-L2652](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2642) |
| The uncorrected R arm reproduces three independent studies digit for digit | +1.112, +1.883, +2.090 | MEASURED, [SDPLS article L2574-L2583](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2537) |
| The coupled translating study `volumeCorrectionTranslating2D` | VOID (closed box) | VOIDED, [STATUS L3566-L3569](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3566-L3569), [[retractions/closed-box-translating-droplet]] |

## Why it failed, or why we think so

The shift is one scalar for the whole field, so it displaces the zero contour by `eps/abs(grad psi)` everywhere at once. It can only remove a volume error that is itself a nearly uniform interface displacement. Where the uncorrected error is already small, the root is fitted to residual noise, and the shift displaces the interface by that noise ([SDPLS article L2614-L2630](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2602)). For several disconnected bodies one `eps` solves the sum of `m` equations; the remedy, a per-body `eps_i` from connected-component labelling, is not implemented ([SDPLS article L2658-L2666](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2657)).

## Decisions

- `VOLUME_CORRECTION noVolumeCorrection` in the best configuration ([METHOD 8.1 L386](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L386)); never on the coupled droplet, where the volume loss is the symptom of the mode-4 current amplified 260x ([SDPLS article L2667-L2672](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2657)).

## Open questions

1. A per-body correction (connected-component labelling of `{psi < 0}`) is documented and not implemented.
2. The fitted floor values (5.923, 6.432, 0.835) must be suppressed by the curation script rather than printed ([SDPLS article L2642-L2655](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2642)).

## Related

[[hubs/advection]] - [[models/phase-indicator]] - [[models/sdpls-source]] - [[models/redistancer]] - [[concepts/error-vector-and-read-out-instants]] - [[concepts/sdpls-source-eulerian]] - [[studies/sdpls-pre-print]] - [[retractions/closed-box-translating-droplet]]

## Log

### 2026-09-28
Created from the SDPLS article section on the global volume correction and the family headers.
