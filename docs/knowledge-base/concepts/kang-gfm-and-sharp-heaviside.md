---
title: "Kang GFM correction and the sharp Heaviside"
description: "The Kang face weights improve the static balance of the N = 64 droplet 6 times, 58 times with the sharp Heaviside, and both arms fail earlier than the default; the 58x belongs to the two options together (retracted, 2026-09-29)."
aliases: [Kang GFM, interfaceWeighted, sharpHeaviside, LS-SSF, Kang face curvature, the 58x attribution flag]
kind: concept
status: retracted
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-29
date_settled:
decided_by:
code: [src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.H, src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C, src/leiaLevelSet/surfaceTensionForce/correctionKang.C, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H]
sources: [SL article sec:droplet delivery study, SL negative deck 4/0 and 4/3, SL deck 4/28 and 4/44, PCS 3 item 3, T2 L207-L217, RM face gate, RM replay gate]
---
# Kang GFM correction and the sharp Heaviside

> Verdict (2026-09-29). The Kang/GFM face curvature interpolates the two cell curvatures to the `psi = 0` crossing, `kappa_f = (kappa_P |psi_N| + kappa_N |psi_P|)/(|psi_P| + |psi_N|)` ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1951-L1956), [doi:10.1023/A:1011178417620](https://doi.org/10.1023/A:1011178417620)). The sharp Heaviside puts the force only on the faces where `psi` changes sign ([`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L188-L210)). On the N = 64 stationary droplet of the Euler-era pipeline, Kang with the geometric alpha lowers the initial peak current 6 times, from 4.5e-3 to 7.0e-4. The current then grows to about 4e-3 and still rises, while the default decays to 1.3e-5 ([negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L211-L214)). Kang with the sharp Heaviside lowers it 58 times, to 7.7e-5. That arm diverges at t = 0.073 s at N = 64 and at 0.049 s at N = 128, against 0.44 s and 0.105 s for the default ([negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [negative deck 4/3](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/3)). The record gives only peak values (L_inf) for these arms. Both arms replicate the meta-law of [[concepts/integral-surface-tension-cst]]: better static balance, higher dynamic gain. Attribution flag: the plan and one deck slide give "58x, blows at 0.07" to Kang alone ([PCS 3](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-curvature-stabilization.md#L369-L370), [SL deck 4/28](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/28)). Thirteen case templates give "58x, blows at t=0.073" to `forceWeight sharpHeaviside` alone ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L215-L217)). The record supports neither single attribution: the number belongs to the arm with both options. The sharp Heaviside with the arithmetic face value was never measured, and the code cannot run it, because the arithmetic branch ignores `forceWeight` ([`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L218-L227)). Production keeps the code defaults `faceInterpolation arithmetic` and `forceWeight alpha` ([[decisions/surface-tension-reconstructed-curvature]]).

## What it is

The level-set curvature of a cell is the curvature of the contour through the cell centre, not of the interface. For parallel curves `kappa_d = kappa/(1 + d kappa)`, which is nearly linear in `d` over `+-h`. The inverse-distance weight evaluates that law at `d = 0`, so the Kang value is close to the interface curvature on a face that straddles `psi = 0` ([SL deck 4/44](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/44), [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L121-L133)). Kang, Fedkiw and Liu introduced the weight; Abadie et al. showed exact balance with it in the ghost-fluid setting ([doi:10.1016/j.jcp.2015.04.054](https://doi.org/10.1016/j.jcp.2015.04.054)).

The sharp Heaviside is `H_s = (1 - sign(psi))/2`. Its `snGrad` is non-zero only on the sign-change faces, so every force face straddles `psi = 0` and the Kang formula is a true interpolation on each of them. The jump integral `sigma kappa [H]` is unchanged; only the localisation of the force changes. Density, viscosity and the mass flux keep the geometric alpha ([`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L188-L197)). With the geometric alpha weight the force support has two to three faces. The faces on one side of the interface then receive a Kang value that is not an interpolation to a crossing ([`reconstructedCurvature.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.H#L91-L95), [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L190-L195)).

The two words are sub-options of `reconstructedCurvature`: `faceInterpolation arithmetic | interfaceWeighted` and `forceWeight alpha | sharpHeaviside`, with the defaults `arithmetic` and `alpha` ([`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L73-L80)). A separate model, `correctionKang`, applies the same weights to a finite-volume curvature ([[models/surface-tension-force]]).

## Why it matters

1. It is the first of the two published cures for the static imbalance that the record tested ([SL deck 4/44](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/44)). It shows that a better static balance is not evidence of a stable coupling. The deck counts the two Kang arms as replications 1 and 2 of that tension ([negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0)).
2. The attribution flag shows how the number of a combined arm moved to one option in copied comments. A result must name every option of its arm.

## Where in the code

1. The two words and their meaning: [`reconstructedCurvature.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.H#L82-L95).
2. The Kang branch, with coupled patches weighted from the neighbour-side `psi` and `kappa`: [`reconstructedCurvature.C`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L121-L186).
3. The sharp Heaviside inside the Kang branch ([L188-L210](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L188-L210)) and with a registered face curvature ([L100-L114](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L100-L114)).
4. The arithmetic branch always uses `alpha` ([L218-L227](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L218-L227)). The constructor prints the two words and does not check them ([L73-L84](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/surfaceTensionForce/reconstructedCurvature.C#L73-L84)), so `forceWeight sharpHeaviside` with `arithmetic`, or a misspelt word, runs the arithmetic CSF without a warning (DERIVED). The code is the same since the import 8de7503 (2026-07-30).
5. The case templates set neither word; they render only `type` and `faceCurvatureSource` ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L231-L275)).
6. The inverse-`|psi|` Kang weight survives in the per-face delivery `stabilizedFootPointFace`. There it combines the two sides of the foot distance and of the Gaussian curvature ([`stabilizedFootPointFaceCurvature.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H#L130-L148)). The plan records that it uses `|psi|` as a length ([PCS 11.3](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-curvature-stabilization.md#L880-L884)).

## Evidence

The coupled arms ran on the N = 64 stationary droplet of the negative-results deck: Euler momentum, `dt_sigma/4`, the cell-local symbolic curvature and CSF ([negative deck 0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/0)).

| claim | number | where |
|---|---|---|
| the default decays | initial 4.5e-3; settles at 1.3e-5; stable to 0.44 s | [negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [SL deck 4/44](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/44), MEASURED (L_inf only) |
| Kang with the geometric alpha grows | initial 7.0e-4; grows to about 4e-3 and rises (the template says 4.4e-3); no blow-up time recorded | [negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L211-L214), MEASURED (L_inf only) |
| Kang with the sharp Heaviside diverges earlier | initial 7.7e-5; divergence at 0.073 s (N = 64) and 0.049 s (N = 128); default 0.44 s and 0.105 s | [negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [negative deck 4/3](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/3), MEASURED |
| the two factors | 4.5e-3 / 7.7e-5 = 58.4 (the "58x"); 4.5e-3 / 7.0e-4 = 6.4 (the "6x") | DERIVED from [negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0) |
| the article states the pairing | "up to 58x (paired with a sharp Heaviside support)"; the sharp pairing diverges near 0.07 s | [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1956-L1959), MEASURED |
| Kang does not fix the static face accuracy | Kang h^1.14, 8.42 1/m at N = 512; arithmetic h^1.13, 11.35; stabilised foot point h^2.04, 0.105 | [RM face gate](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/capillary-level-set-research-roadmap.md#L1310-L1316), MEASURED (band L2) |
| the one-step replay on the translating droplet | quadratic + Kang 1.19e-3 against 7.11e-3 for the cell-centre value at t = 0; 1.08 against 2.03 m/s at t = 0.05 s | [RM replay](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/capillary-level-set-research-roadmap.md#L735-L743), MEASURED (L_inf only) |
| `correctionKang` on the N = 64 translating matrix | reaches 0.05 s with a peak disturbance of 0.167 m/s, above U0 = 0.05 m/s | [`transISTN64ForceFluxModelMatrix.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/method-comparison/method-comparison-article/data/tables/transISTN64ForceFluxModelMatrix.csv#L7), MEASURED (L_inf only) |

### What the record supports and what it does not

| place | option it names | what it states |
|---|---|---|
| [PCS 3 item 3](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-curvature-stabilization.md#L369-L370) | Kang alone | "statics 58x better, blows earlier (0.07 vs 0.44)" |
| [SL deck 4/28](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/28) | Kang alone | "statics 58x better, coupled loop worse (blows 0.07 vs 0.44)" |
| [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L215-L217) and twelve more case templates | `forceWeight sharpHeaviside` alone | "58x better initial balance, blows at t=0.073 (N=64)" |
| [negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [4/3](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/3), [SL deck 4/44](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/44), [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1956-L1959) | both options together | 7.7e-5 initially, divergence at 0.073 s; "58x better"; 0.073 / 0.049 s |

1. Supported: 58x and 0.073 s are the arm with `faceInterpolation interfaceWeighted` and `forceWeight sharpHeaviside` together (MEASURED; the factor is DERIVED above).
2. Supported: Kang with the alpha weight is 6x better statically and its current rises; the record has no blow-up time for it (MEASURED).
3. Not supported: any number for the sharp Heaviside with the arithmetic face value. No arm is recorded, and the code runs the alpha weight in that case (DERIVED).
4. Not supported: "Kang: 58x, blows at 0.07" and "sharpHeaviside: 58x, blows at 0.073". Each gives the number of the combined arm to one option.
5. Not recorded: the configs, the commits and the binaries of the four arms.

## Why it failed, or why we think so

1. The cell values of `kappa` still carry the aliasing response of the fit to the distortion of `psi`. The Kang weight cannot repair corrupted inputs, and a sharper delivery raises the gain from `psi` to the force. The arithmetic average and the smeared support act as a spatial low-pass filter ([negative deck 4/0](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set-negative-results.html#/4/0), [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1983-L1989)). This is the record's reading of the arms, not a separate measurement (HYPOTHESIS); the current mechanism is in [[concepts/parasitic-current-mechanism]].
2. The static accuracy did not change the order: with the alpha weight the Kang face value stays first order, 8.42 against 11.35 1/m at N = 512 ([RM face gate](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/capillary-level-set-research-roadmap.md#L1310-L1316)). The per-face re-referencing through the stabilised foot point gave second order later ([[concepts/face-curvature-deliveries]]).
3. The reference times 0.44 s and 0.105 s are Euler-era blow-up times. A blow-up time is not a proxy for the growth rate, and no result may be compared with an older blow-up table ([[retractions/t-blow-baseline]]).

## Decisions

- Production uses `faceInterpolation arithmetic` and `forceWeight alpha`, the code defaults ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L391-L392), [[decisions/surface-tension-reconstructed-curvature]]).
- The case templates say "do not use coupled" for the sharp Heaviside ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/fvSolution.template#L215-L217)).
- A face delivery is now scored on its noise gain and on the ellipse order, not on the static balance ([[concepts/face-curvature-deliveries]]).

## Open questions

1. The attribution in PCS 3 item 3, SL deck 4/28 and the thirteen case templates is not corrected at d1e3414.
2. No arm with the sharp Heaviside and the arithmetic face value exists. The code ignores `forceWeight` there without a warning; a check of the two words at construction would make the combination explicit.
3. None of the arms ran on the current pipeline (BDF2, `cellCentreInverse`, rhoLENT); the record has no Kang result after the Euler era except the replay and the `correctionKang` matrix.

## Related

[[hubs/surface-tension]], [[models/surface-tension-force]], [[concepts/integral-surface-tension-cst]], [[concepts/face-curvature-deliveries]], [[concepts/balanced-force-csf-flux]], [[concepts/parasitic-current-mechanism]], [[concepts/curvature-corrugation-and-the-fit]], [[decisions/surface-tension-reconstructed-curvature]], [[retractions/t-blow-baseline]], [[cases/stationary-droplet]], [[studies/sl-quadratic-pre-print]], [[studies/curvature-stabilization-campaign]].

## Log

### 2026-09-29
Created from the SL article delivery study, the two SL decks, PCS section 3, the case template comments, the roadmap gates and the `reconstructedCurvature` source, all read at d1e3414. The attribution flag (raw-material audit item 11) is resolved as stated in the Evidence section.
