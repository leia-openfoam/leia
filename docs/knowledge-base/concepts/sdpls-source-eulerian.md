---
title: "The SDPLS source in the Eulerian transport: kinematic gain, coupled failure"
description: "The strain source R keeps |grad psi| near one in kinematic flows and amplifies the curvature-seeded spurious current 260 times under coupling; exact curvature removes the seed"
kind: concept
status: settled
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [SDPLS article, plan-combined-source-terms, STATUS 9]
---
# The SDPLS source in the Eulerian transport: kinematic gain, coupled failure

> Settled (measured 2026-08/09). Kinematic: `R` raises the band gradient order from −0.26 to
> +0.74 on the non-reversing vortex and from −0.094 to +0.668 on the 3D shear, at 3.7x worse shape
> and 18x worse volume ([SDPLS article, verification](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L232)).
> Coupled: with the production curvature `R` drains the droplet at $N = 32$ and diverges at 64 and 128,
> the mesh-locked mode-4 current is amplified 260 times, and with exact curvature the current is
> 0.000 in 12 of 12 arms ([SDPLS article, two-way coupling](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L842)).
> The seed is the curvature error; the source reads the strain of the spurious current and closes the loop.

## What it is

$\partial_t\psi + \mathbf{u}\cdot\nabla\psi = \psi\,a$ with $a = \mathbf{n}\cdot\nabla\mathbf{u}\cdot\mathbf{n}$:
along a characteristic $q$ obeys $\dot q = -a q$ for a signed-distance start, so the source
$\psi a$ cancels the drift at first order without moving the zero set in the continuum. The
continuum diagnosis and the noise-gain bound $1 + \lambda\Delta t\,\pi W$ with the design window
$\lambda\Delta t \lesssim 1/(\pi W)$ are in
[plan-combined-source-terms](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L97-L165).

## Evidence

| claim | number | where |
|---|---|---|
| order ceiling under second-order transport | about 1.2 | [SDPLS article, where the order is not](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L406) MEASURED |
| R amplifies solver-tolerance noise | 7 orders over 130 steps | [STATUS 9](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2836-L2856) MEASURED |
| Rdiv (material form) kinematic orders | −3.8, −1.3, −1.5; withdrawn | [SDPLS article](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L1263) MEASURED |
| exponentialImplicit discretisation | shape order +1.603 against +1.360 | [SDPLS article](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2162) MEASURED |
| the np 8 seam residual of the Eulerian coupled solver | 2.8e-3 | same article, coupled-patch defects MEASURED |

## Why it failed, or why we think so

Under coupling the source reads $a$ of the spurious current; that current is seeded by the
curvature error ([[concepts/parasitic-current-mechanism]]) and the source turns the seed into a
per-cell geometric correction of $\psi$. This is a loop through the flow, and it is the mechanism
the gcls pre-print wrongly transferred to the linear laws of the SL source step, which read no
velocity ([[retractions/gcls-coupled-loop-reading]]). For the soft wall, which reads $\sigma$, the
SDPLS reading does hold ([[concepts/source-discretisation-defects]]).

## Decisions

`noSource` is the default; the weighted-SDPLS route (a strain weight $w \le 1$) is abandoned for
now: the weight reduces the coupled gain by at most $w$, and $w = 1$ in the geometry cells by
design (technical report).

## Related

[[hubs/gradient-control]], [[models/sdpls-source]], [[models/gradient-control-law]],
[[concepts/why-the-candidates-failed]], [[studies/sdpls-pre-print]].

## Log

### 2026-09-28
Created.
