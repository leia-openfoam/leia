---
title: "The first 2D method-gate campaign (2026-09-26 to 2026-09-27)"
description: "Six gradient-control candidates against the production method in the fixed 2D gate on Lichtenberg: every candidate FAILS; the pre-fix gate, the fixed gate, the seam checks, the scoring corrections, and where the data live"
kind: study
status: settled
part: gradient-control
tags: [study, part/gradient-control]
date: 2026-09-28
sources: [STATUS 11.11 to 11.16, config/gates/methodGate2D.yaml, config/candidates/*.yaml]
---
# The first 2D method-gate campaign (2026-09-26 to 2026-09-27)

> Settled 2026-09-27 (measured). The gate [methodGate2D.yaml](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml)
> ran on Lichtenberg twice: with the pre-fix binaries (orchestrator 55044205, summaries
> `methodGate2D_summary_pre-20260927-020856`) and with the two parallel fixes of 2026-09-27
> (orchestrator 55048916, stamp `shared-method-config-2026-09-01-192-g1150e68`). Every candidate
> FAILS in both; the kinematic arms of the fixed gate reproduce the pre-fix ones byte for byte
> (177 CSV pairs) ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)).

## The candidates

| candidate | velocity extension (trace flux) | source law (strain weight none) | band |
|---|---|---|---|
| baseline | none | none | — |
| S1 | none | softWall, $C_\kappa = 1.25$, $\delta_s = 0.08$, $p = 5$, $\gamma = \mathrm{artanh}\,0.9$ | 3h |
| HL0 | haloLimited, $R = 1h$, $m = 2$, $\beta = 1$ | none | — |
| HL1q | haloLimited, $R = 1h$ | linearQ, $\mu = 1/T_\mathrm{ref}$ | 3h |
| HL1z | haloLimited, $R = 1h$ | linearZ, $\mu = 1/T_\mathrm{ref}$ | 3h |
| HL2 | haloLimited, $R = 1h$ | softWall (as S1) | 3h |
| FP0 | closestPoint | none | — |

The configs with their pre-registered read-outs: [config/candidates/](https://github.com/leia-openfoam/leia/blob/8867581/config/candidates).

## The verdicts (fixed gate, corrected scoring)

| candidate | verdict | target ratio | regressions | orders | seam | completion |
|---|---|---|---|---|---|---|
| FP0 | FAIL | 0.939 | 3 | 8 | 5 | 5 |
| HL0 | FAIL | 0.828 | 16 | 13 | 0 | 0 |
| HL1q | FAIL | 0.204 | 22 | 18 | 0 | 0 |
| HL1z | FAIL | 6.490 | 23 | 19 | 0 | 0 |
| HL2 | FAIL | 0.754 | 7 | 5 | 0 | 3 |
| S1 | FAIL | 0.707 | 12 | 5 | 1 | 1 |

Coupled seam check (`translatingSeamNp1`, tolerance 1e-5): baseline 3.5e-7, HL0 4.7e-6, HL1q
4.9e-7, HL1z 5.3e-7 PASS; S1 0.76 FAIL; FP0 1.11 FAIL; HL2 NOT_COMPARABLE. Per-arm numbers and the
mechanisms: [[concepts/why-the-candidates-failed]]. The 3D gate does not run: no 2D pass.

## What the campaign changed besides the verdicts

- Two parallel defects of the coupled SL solver found and fixed ([[concepts/coupled-face-density-defect]]);
  the Eulerian rhoLENT port ([[concepts/eulerian-solver-mass-flux-port]]).
- Two scoring corrections: `rhoClipFraction` reported, not scored; the oscillating arm scored by
  period and damping ([[concepts/method-gates]]).
- The translating horizon cut from the pre-registered 0.1 s to 0.05 s after the baseline diverged
  (a post-hoc change, stated as such); the late instability studied on the laptop
  ([[cases/translating-droplet]]).
- The baseline's oscillating arm found unstable at $N = 200$ ([[cases/oscillating-droplet]]).
- Two gate blind spots found on 2026-09-28 ([[concepts/method-gates]]).

## Where the data live

The archive `docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/`
(`gate/`, `prefix/`, `laptop/`; README and MANIFEST) ([[concepts/data-archive-per-version]]); the
tables `data/tables/methodGate2D_*`; the raw runs on Lichtenberg under `studies/methodGate2D_*`.
Documents: [[studies/gcls-pre-print]] and the technical report.

## Related

[[hubs/gradient-control]], [[concepts/method-gates]], [[models/sl-source]],
[[concepts/halo-limited-extension]], [[sessions/2026-09-27-gcls-first-campaign]].

## Log

### 2026-09-28
Created.
