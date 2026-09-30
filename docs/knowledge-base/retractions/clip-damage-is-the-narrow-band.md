---
title: "\"The clip's damage is the narrow band\" (2026-09-09)"
description: "RETRACTED 2026-09-09 - the +30 % volume error of the quasi-monotone clip was attributed to clipping the narrow band; a band-aware clip removes none of it, because the clip fires only in the six cells where psi is its own stencil extremum"
aliases: [band-aware clip retraction, clipRegion outsideBand]
kind: retraction
status: retracted
part: advection
tags: [retraction, part/advection]
date: 2026-09-28
code: [config/popinet2D_clipRegionGate.yaml, config/popinet2D_clipFiringProbe.yaml, config/popinet2D_clipStencilGate.yaml, config/popinet3D_poly_sigma0_clipGate.yaml, cases/default.parameter]
sources: [STATUS 4 (2026-09-08 and 2026-09-09), ROADMAP-poly3D 2 and 2a, METHOD 8.1 rows SL_CLIP, DP SL_CLIP_REGION]
---
# "The clip's damage is the narrow band" (2026-09-09)

> RETRACTED 2026-09-09. The claim was "the clip stops the polyhedral far-field failure but is not inert at the interface (+30 % volume error); the clean form is a band-aware clip that acts outside the narrow band only, so the interface metrics are untouched by construction" ([STATUS 4, 2026-09-08](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2033-L2039); [ROADMAP-poly3D 2](https://github.com/leia-openfoam/leia/blob/8867581/ROADMAP-poly3D.md#L35-L44), requirement 3). The measurement: `popinet2D_clipRegionGate` on Popinet's 2D hexahedral translating droplet at N = 64 over 1563 steps. `clipRegion outsideBand` withholds the clip from 144 band cells and moves every interface metric by exactly what the global clip moves it: volume +30.4 %, centroid +12.5 %, shape +4.9 %, L1|u'| +8.9 %, to every printed digit ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2398-L2407)). The firing probe located the damage: the mesh has exactly six cells where psi is the extremum of its own stencil, and every firing was in one of them ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2409-L2415)). Scope: the mechanism and the band-aware remedy. The +30 % cost stands. The replacement, the extremum exemption, was falsified by gate G4 the same day ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2312-L2338)). No data is void.

## The claim, and where it lived

- [STATUS 4, DECIDED 2026-09-08](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2033-L2039), `STATUS.md`: "The clean form is a BAND-AWARE clip ... The interface metrics are then untouched by construction."
- [ROADMAP-poly3D 2](https://github.com/leia-openfoam/leia/blob/8867581/ROADMAP-poly3D.md#L35-L61), `ROADMAP-poly3D.md`: requirement 3, now headed "requirement 3 is RETRACTED".
- [cases/default.parameter](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L64-L76), `cases/default.parameter`: the `SL_CLIP_REGION` comment still describes `outsideBand` as "the form for the polyhedral far-field defect".

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The band is where the clip damages the interface. | `outsideBand` gives the same +30.4 / +12.5 / +4.9 / +8.9 % as `all`, with 144 band cells withheld. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2398-L2407) |
| Where the clip fires. | Six cells: the four box corners (maxima of the distance field) and the two cells at the droplet centre (the apex of the distance cone). The `all` and `outsideBand` arms fire in the same cells. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2409-L2415) |
| Why a monotone bound cannot hold an extremum. | Where psi_c is the stencil minimum, `lo` is psi_c, so every reconstructed value below it is pulled back up. A smooth quadratic undershoots at the apex. | DERIVED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2417-L2421) |
| The classical exemption of a stencil extremum repairs it. | Damage 30x lower: volume +1.03 %, shape +0.2 %, centroid +0.4 %; activity from 44 per step to 3 cell-steps in total. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2431-L2438) |
| The residual of the exemption. | One apex cell. The symmetric tie breaks at (psi_c - lo)/(hi - lo) = 0.004; the first byte difference is at step 20, and three cell-steps of clipping grow to +1.03 % volume error by step 1562. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2452-L2462) |
| The outlet layer is the residual. | `inflowOnly` removes those eighteen firings and leaves the volume error at +1.027 % against +1.0 %: retracted too. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2464-L2473) |
| The exemption keeps the polyhedral fix (gate G4). | With the exemption the polyhedral far field fails at step 506 against 527 without any clip (4 %, inside the 5 to 38 % scatter). 59.2 % of the firings are in stencil extrema (2 680 917 of 6 578 187 cell-steps). | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2312-L2346) |

The two requirements are in direct conflict under this formulation: a growing checkerboard has extrema at its own peaks ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2340-L2353)). The apex detector proposed as the next step was retracted before it was built ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2355-L2360)).

## What survives

1. The global clip removes the polyhedral far-field failure. sigma = 0, polyhedral N = 64, 1563 steps: no fake zero set, zero-set error 3.98e-04 against 8.66e-01, volume error 1.08e-06 against 3.79 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2012-L2026)). In G4 the clipped arm has no failure in 650 steps; its zero-set error is 380x lower and its volume error 26x lower, but the zero-set error still grows from 4.46e-04 to 8.97e-04 ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2329-L2334), [L2383-L2387](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2383-L2387)).
2. The cost on hexahedra (+30.4 %) is real, so `SL_CLIP false` is the production setting ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L378-L380)).
3. The tokens `SL_CLIP_REGION` and `SL_CLIP_KEEP_EXTREMA` shipped inert: the four clip-off arms are bit-identical ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2445-L2450)).
4. Two scale-based discriminators of a genuine extremum are proposed and unmeasured: |grad psi| falls toward 0 at a genuine extremum, and a genuine extremum survives one more stencil layer ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2362-L2374)).
5. The process lesson: read the cumulative firing counter, never a write-time sample ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2475-L2481)).

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: [RETRACTED AND REPLACED](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2389-L2495) and the [G4 falsification](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2312-L2387).
- [x] `ROADMAP-poly3D.md`: [section 2](https://github.com/leia-openfoam/leia/blob/8867581/ROADMAP-poly3D.md#L35-L61) and [2a](https://github.com/leia-openfoam/leia/blob/8867581/ROADMAP-poly3D.md#L136).
- [x] `METHOD.md`: [8.1 rows SL_CLIP, SL_CLIP_REGION, SL_CLIP_KEEP_EXTREMA](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L378-L380) and [9 item 4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L777-L790).
- [x] `cases/default.parameter`: [SL_CLIP_KEEP_EXTREMA](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L77-L92) and [SL_VALUE_BOUND stencilBounds](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L112-L118).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] `cases/default.parameter` [L57-L63](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L57-L63) still says "true = REQUIRED on general polyhedra ... preserves the 2nd-order convergence". The measured record says the clip costs +30.4 % volume error on hexahedra and fails G4 with the exemption.
- [ ] Whether `SL_CLIP_REGION` earns its place on polyhedra, or should be dropped as a second knob ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2440-L2443)).
- [ ] The G4 clipped arm beyond 650 steps: its zero-set error still grows ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2383-L2387)).
- [ ] The scale-based extremum rule is unmeasured ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2362-L2374)).

## Related

Hubs: [[hubs/advection]]. Siblings: [[models/sl-value-bound]], [[concepts/value-bounds-and-clips]], [[concepts/polyhedral-fit-amplification]], [[decisions/sl-clip-and-value-bound-off]], [[decisions/mesh-family-hexahedral]], [[retractions/distance-cone-bound-as-transport-bound]], [[retractions/polyhedral-popinet-3d-mesh-defect]].

## Log

### 2026-09-28
Written from STATUS 4 (2026-09-08 and 2026-09-09) and ROADMAP-poly3D section 2. Retracted 2026-09-09. Entered in [[retraction-log#2026-09]].
