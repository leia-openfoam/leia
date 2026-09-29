---
title: "A wrong setup voids its data"
description: "When a case is found to be set up wrong, every number it produced is void, including the arms that look unaffected; rename the study _VOID_<reason>_<date>, fix the setup, re-run; seven setups of this repository were wrong in this way"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-28
date_settled: 2026-09-02
decided_by: [author decision 2026-09-02]
code: [cases/translatingDroplet2D/system/blockMeshDict.template, cases/popinetTranslating3D_poly.parameter, cases/translatingDroplet3D/system/fvSolution.template, cases/default.parameter]
sources: [CLAUDE wrong setup section, STATUS 0, STATUS 4 tilted-wall section, STATUS 11.4, STATUS 11.14, CLUSTER cancel section]
---
# A wrong setup voids its data

> Verdict (2026-09-28). When a case is found to be set up wrong, every number it produced is void, including the ones that look unaffected and the arms that were merely adjacent to the defect. The study directory is renamed with a `_VOID_<reason>_<date>` suffix so it can never be curated by accident, the setup is fixed, and the study is re-run from scratch; exactness at the cost of CPU hours ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L661-L706)). Nobody reasons about which metrics "should still be valid" under the defect: that reasoning is as unreliable as the setup was, and a partially trusted table is worse than none. The rule was written on 2026-09-02, when `translatingDroplet2D` was found to be a closed slip box: its mesh had no inlet and no outlet while every field declared both, OpenFOAM silently ignored the field entries, and the pressure projection annihilated the uniform stream on step 1, so `max|U-U0|` measured the dead free stream at about `2 U0` on every arm of every study ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L54), [[retractions/closed-box-translating-droplet]]).

## What it is

The procedure:

1. Preserve first: rename the study directory with the suffix; solvers truncate their metrics CSV on restart ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L761-L765)).
2. Move the curated tables into a `VOID_<reason>_<date>/` folder with a README that says what was wrong and what replaces them ([VOID README](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/VOID_closedBox_20260902/README.md)).
3. Mark every conclusion that rests on the data as VOID or RETRACTED where it stands: STATUS.md, METHOD.md, the token comments, the configs ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L248-L269), [[concepts/bound-rho]]).
4. Stop a run in flight on the same case by job id, even when its own argument looks boundary-independent: a setup wrong in one way is not assumed wrong in only that way ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L691-L695)).
5. Gate the repair with a pre-registered read-out before any production arm is launched (`translatingFreeStreamGate2D`, 200 steps, N = 128, minutes on the laptop, [config](https://github.com/leia-openfoam/leia/blob/8867581/config/translatingFreeStreamGate2D.yaml#L1-L45)).
6. Re-run.

The check that belongs in every new case's gate: `constant/polyMesh/boundary` is the authority on what exists; the `boundaryField` entries of `0.org` are a wish list, and a field patch name not in the mesh does nothing, silently ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L696-L706)). No case in `cases/` commits a `boundary` file; the mesh is built per case, so the check runs after `blockMesh` or `pMesh`.

## The setups that were wrong

| setup | what was wrong | found | what was voided |
|---|---|---|---|
| `translatingDroplet2D`, closed box | all four sides in one `walls` patch; `inlet` and `outlet` existed only in the fields; first-step continuity error 1.00e-05 against 2.32e-20 after the fix, `max|U-U0|` at step 1 1.031e-01 against 2.159e-03 | 2026-09-02, fix 440107f | every translating result before the fix: `alphaFTest_*`, `bestConfigTranslating2D`, `rhoBoundGate2D`, `rhoDdtGate2D`, `rhoLENTGate2D`, `translatingMap2D`, `wellBalancedTranslating2D` and the studies of the same morning ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L34-L62)) |
| the polyhedral Popinet 3D mesh | the box STL put the four side walls in one solid, so the edges were not feature edges and cfMesh tilted 1 984 of 40 960 wall faces (4.8 %) up to 8 degrees into the flow; `simpleFoam` could not hold a uniform stream (8.3 % in the corner cells); with the feature-edge surface 0 tilted faces and `|U-U0| <= 5.7e-16` | 2026-09-05 | every polyhedral Popinet 3D result on the STL mesh: two ladder rungs, the smokes, the field dumps, the corrector sweep's translating arms, a census row; renamed `_VOID_tiltedWallFaces_20260905` ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1738-L1771), [[retractions/polyhedral-popinet-3d-mesh-defect]]) |
| stale arms in a resubmitted sweep | a driver cancelled seconds after submission had rendered two arms; the corrected driver reused them: two arms of a different case, dead at step 0, inside a sweep that looked complete | 2026-09-05 | the sweep, re-run with the stale arms removed ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/8867581/CLUSTER.md#L134-L146), [STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1952-L1964)) |
| a vacuous snappyHexMesh arm | the arm looked like a pass while it had not tested what it was meant to test | 2026-09-09 | the arm ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2289)) |
| `oscillatingDroplet2D`, algebraic psi | `implicitEllipsoid` gives `psi = sum (x_i - c_i)^2/a_i^2 - 1`, so `|grad psi|` is 1.8e3 to 2.2e3 at the interface; every gradient column and every band criterion in psi units had no meaning | 2026-09-26 | open: void the past oscillating studies or keep them with the caveat; the gates pin `signedDistanceEllipse` ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3271-L3280), [[cases/oscillating-droplet]]) |
| `translatingDroplet3D`, no reference velocity | the metrics writer used `(0 0 0)`, so the disturbance columns reported the translation itself (`meanMagUPrime = 0.0500 = U0`) and the zero-set and centroid errors the displacement | 2026-09-27 | those columns of `traceTranslating3Dhex` and `traceTranslating3Dpoly_r10p0/r12p7/r15p8`; no conclusion cites them ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3716-L3724)) |
| the Eulerian two-phase solver, frozen density | `rho`, `rhoPhi` and `muf` at their `t = 0` values for the whole run | 2026-09-27 | open: the five Eulerian coupled studies ([[concepts/eulerian-solver-mass-flux-port]]) |
| `2Dtranslation`, reversed flow | made from `2Dvortex` without the `oscillation` line (cd97e6da, 2026-09-01); the velocity model's default multiplied U by cos(pi t/T), so the circle came back at T; with no `psiEnd` the metrics compared T with the start, silently | 2026-09-29 | every study of the case: `kinematicTranslation2D`, `coneBoundMesh2Dtranslation` and its arms, `advConv2Dtranslation`, the translation rows of the regression gates ([STATUS 11.19](https://github.com/leia-openfoam/leia/blob/f47fc939/STATUS.md#L4177-L4189), [[retractions/reversed-2dtranslation]]) |
| `2Dtranslation`, fixed psi on a patch (the first two repairs) | the exact psi on the outflow patch made the SL update unstable at the outflow edge (4.5e5 at T at N = 256); on the inflow patch it failed at CFL 1, where the departure point leaves the domain | 2026-09-29 | the two re-runs `_VOID_outflowDirichlet_20260929` and `_VOID_inflowDirichlet_20260929` ([STATUS 11.19](https://github.com/leia-openfoam/leia/blob/f47fc939/STATUS.md#L4202-L4226)) |

## Why it matters

An entire `div(rhoPhi,U)` scheme comparison, a droplet-leaves-the-domain mechanism and a "scheme-independent kick" were read off the closed box before anyone checked the mesh; the rule "Is the interface still inside the domain?" was first written as "the leading edge reaches the outlet at t = 0.08" when there was no outlet ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L484-L501)). The tilted faces hid an "outlet checkerboard amplifier" that lived in the slab next to them and had to be re-measured ([STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1760-L1771)).

## Evidence

| claim | number | where |
|---|---|---|
| the closed box annihilated the stream | first-step continuity error 1.00e-05 (closed) against 2.32e-20 (repaired); `mean|U-U0|` 5.011e-02 against 1.160e-05; `maxMagU` 0.05987 against 0.05207 = U0 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L34-L44), MEASURED |
| the repaired case has a real parasitic current | `max|U-U0|` settles around 5e-03; `mean|U-U0|` drifts 1.2e-05 to 7.6e-05 over 200 steps, about 50x below the artefact | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L48-L54), MEASURED |
| the tilted faces | 1 984 of 40 960 wall faces, up to 8 degrees; 0 with the feature-edge surface | [STATUS](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1738-L1755), MEASURED |
| a mid-run row was published as final | the aggregator published the last row of a diverged run (`alphaFTest_averagedPlanes`, step 4634) as its final value until 2026-09-26 | [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3244-L3249), MEASURED |

## Decisions

- The rule and its corollary (author decision, 2026-09-02): [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L661-L706).

## Open questions

1. The algebraic-psi oscillating studies and the frozen-rho Eulerian studies: void or caveat (author decisions).
2. The parallel SL two-phase studies before 28d13f0 are a different class (a code defect, not a setup defect) and get their own decision ([[concepts/coupled-face-density-defect]]).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/bit-identity-and-inertness-gates]], [[concepts/cluster-provenance-and-binaries]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/method-gates]].
- Retractions: [[retractions/closed-box-translating-droplet]], [[retractions/polyhedral-popinet-3d-mesh-defect]], [[retractions/late-translating-instability-is-the-outlet]].
- Cases: [[cases/translating-droplet]], [[cases/popinet-translating-droplet]], [[cases/oscillating-droplet]].

## Log

### 2026-09-28
Created.

### 2026-09-29
Added the reversed `2Dtranslation` and the two fixed-psi repairs of the same day to the table of wrong setups.
