---
title: "Session 2026-09-26 to 2026-09-28: gradient control, first campaign"
description: "What the session built and changed, by commit: the law family, the SL source step, the halo-limited extension, the method gates, two parallel defects fixed, the Eulerian rhoLENT port, the pre-print, the knowledge base and the technical report"
kind: session
status: settled
part: gradient-control
tags: [session, part/gradient-control]
date: 2026-09-28
---
# Session 2026-09-26 to 2026-09-28: gradient control, first campaign

> Branch `feature/gradient-controlled-level-set`; the plan
> [docs/plan-halo-limited-gradient-control.md](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md)
> (approved 2026-09-26); the record [STATUS 11](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3218-L4041).

## What the session built, by phase

| phase | what | where |
|---|---|---|
| A | the rules, the plan, the docs theme; the combined-source plan subsumed | [STATUS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3228-L3234) |
| B1 | workflow repairs; the 2D advection orders of METHOD 8.3.7 retracted (3/2 too high) | [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3260), [[retractions/advection-orders-3-2-factor]] |
| B3 to B6 | the band metric uses the unlimited scheme; the oscillating surface is a token; the 2D and 3D gates exist | [STATUS 11.3 to 11.5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3261-L3306) |
| C | the extension points of the SL solver, 53 cases bit-identical, 250 renders | [STATUS 11.7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3318-L3370) |
| D1 to D4 | the eight laws and the strain weights, the SL source step, the SDPLS gradientControl source (equals R bit for bit), the halo-limited extension | [STATUS 11.9 to 11.10](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3403-L3471), [[models/gradient-control-law]], [[models/sl-source]], [[concepts/halo-limited-extension]] |
| gate | the 2D gate on Lichtenberg with six candidates; every candidate FAILS | [STATUS 11.12 and 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3504-L3544), [[studies/method-gate-2d-campaign-2026-09]] |
| fixes | the coupled-face density (28d13f0) and the droplet metrics on processor faces (b1798c3); the Eulerian rhoLENT port (c094bd8); the 3D translating reference velocity | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3624-L3744), [[concepts/coupled-face-density-defect]], [[concepts/eulerian-solver-mass-flux-port]] |
| record | METHOD.md corrected (4.1, 4.3, 6, 8.1, 10); the translating arm re-run with `cellCentreInverse` to 0.05 s | [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3545-L3623) |
| late instability | the outlet triggers the fast growth; a slower interior growth shrinks with h | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023), [[cases/translating-droplet]] |
| pre-print | the numerical paper with its data archived per software version | [STATUS 11.16](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L4024-L4041), [[studies/gcls-pre-print]] |
| 2026-09-28 | the knowledge base, the paper split, the technical report, the CLAUDE.md rules, the record corrections | [[sessions/current]] |

## What is open, in order

See [[sessions/current]] and [[concepts/gradient-control-open-decisions]].

## Traps that cost time

- The spend limit of the API stopped eight parallel agents on 2026-09-28; relaunched after the reset.
- `lcluster5` stopped answering at about 02:00 on 2026-09-27; `lcluster4` was used.
- The pre-registered translating horizon (0.1 s) was cut to 0.05 s after the baseline diverged; stated as a post-hoc change in the record.
- The pre-registered HL0 read-out conflated the interface with the band ([[concepts/extension-strain-relocation]]).

## Where the numbers live

The archive `shared-method-config-2026-09-01-192-g1150e68` ([[concepts/data-archive-per-version]]).

## Log

### 2026-09-28
Frozen.
