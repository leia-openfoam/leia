---
title: "sdplsSource: noSource, R, beta, Rdiv, RdivStrictSp, gradientControl"
description: "The source-term family of the Eulerian level-set transport (SDPLS): R converges kinematically and fails under coupling; beta's target is structurally wrong; Rdiv is withdrawn; gradientControl reproduces R bit for bit"
kind: model
status: settled
part: gradient-control
tags: [model, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/sdplsSource/]
sources: [SDPLS article, STATUS 9, STATUS 11.9, cases/default.parameter]
---
# sdplsSource: noSource, R, beta, Rdiv, RdivStrictSp, gradientControl

> Settled (measured, 2026-08 to 2026-09-26). The family adds $\psi F$ to the Eulerian level-set
> equation of `leiaLevelSetFoam` and `leiaLevelSetTwoPhaseFoam`. The strain source `R` keeps
> $\lvert\nabla\psi\rvert$ near one in kinematic flows (band gradient order +0.74 against −0.26
> without it), but under two-phase coupling it reads the strain of the spurious current and
> amplifies the mesh-locked mode-4 current 260 times; with exact curvature the loop has no seed
> ([[concepts/sdpls-source-eulerian]]). Details and every number: [[studies/sdpls-pre-print]].

## What it is

| member | dictionary word | status | one line | evidence |
|---|---|---|---|---|
| no source | `noSource` | settled | the default | — |
| strain source | `R` | settled (kinematic positive, coupled negative) | $F = a = \mathbf{n}\cdot\nabla\mathbf{u}\cdot\mathbf{n}$; band gradient order +0.74 against −0.26 on the non-reversing vortex (31x); on the 3D shear +0.668 against −0.094 but 3.7x worse shape and 18x worse volume | [SDPLS article, verification](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L232) |
| restoring source | `beta` | retracted | a linear restoring law in $\lvert\nabla\psi\rvert$ whose target is structurally wrong (fixed point $\beta - a$, band mean 1.516 to 1.481, flat); explicit form diverges at every CFL | [SDPLS article](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2862) |
| material form | `Rdiv` | retracted | $R + \nabla\cdot\mathbf{u}$; diverges kinematically (orders −3.8, −1.3, −1.5); withdrawn | [SDPLS article, Rdiv and a retraction](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L1263) |
| material form, strict Sp | `RdivStrictSp` | retracted | the implicit variant of Rdiv | same |
| gradient control | `gradientControl` | settled | the law family of [[models/gradient-control-law]]; with law `none` and weight `full` it equals `R` bit for bit (189 checks) | [STATUS 11.9](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3403-L3436) |

Strategies: `discretization` {none, explicit, simpleLinearImplicit, strictNegativeSpLinearImplicit,
exponential, exponentialImplicit} (exponentialImplicit repairs the shape order, +1.603 against
+1.360), `gradPsi` {fvc, narrowLS}, `mollifier` {none, m1, band}
([sdplsSource/](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/sdplsSource)).

## Why it matters

It is the first source-term line of leia, and its coupled failure is the reference for every
later source: a source that reads the flow closes a loop through the capillary force whose seed
is the curvature error ([[concepts/parasitic-current-mechanism]]). The SL source step of the first
gradient-control campaign has a different primary defect ([[concepts/source-discretisation-defects]]).

## Evidence

| claim | number | where |
|---|---|---|
| coupled: R drains the droplet at $N = 32$, diverges at 64 and 128; mode-4 current amplified | −1.000; 260x | [SDPLS article, two-way coupling](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L842) MEASURED |
| exact curvature: the current with R | 0.000 in 12 of 12 arms | same, MEASURED |
| R amplifies the solver-tolerance noise | 7 orders over 130 steps; the seam gate needs a psi tolerance of 1e-14 | [STATUS 9](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L2836-L2856) MEASURED |
| Rdiv at $N = 128$ | DIVERGED | [cases/default.parameter](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L186-L189) MEASURED |

## Decisions

`noSource` is the default; `R` is a research instrument. The combined-source plan is subsumed
([[concepts/combined-source-note]]).

## Open questions

Whether the strain reading can enter any coupled run before the curvature-error seed is smaller:
the technical report says no (the weighted-SDPLS route is abandoned for now).

## Related

[[hubs/gradient-control]], [[concepts/sdpls-source-eulerian]], [[models/gradient-control-law]],
[[models/sl-source]], [[studies/sdpls-pre-print]], [[models/level-set-advection]].

## Log

### 2026-09-28
Created.
