---
title: "slSource: the semi-Lagrangian source step"
description: "The source step psi <- psi exp(dt F) after the characteristic step, in a band of 3h, clamped at |dt F| = 30; candidate, and the primary suspect of the first campaign's failures"
kind: model
status: candidate
part: gradient-control
tags: [model, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/semiLagrangian/source/slSource.H, src/leiaLevelSet/semiLagrangian/source/slGradientControlSource.C, src/leiaLevelSet/semiLagrangian/slAdvection.C]
sources: [STATUS 11.9, gcls article sec:source-step, technical report sec 3]
---
# slSource: the semi-Lagrangian source step

> Candidate (2026-09-26); every gate arm with it FAILS (2026-09-27). Members: `none` (default,
> the identity) and `gradientControl` (a [[models/gradient-control-law]] in a band of
> `bandCells` $h$). The update is $\psi_c \leftarrow \psi_c \exp(\Delta t\,F_c)$ after the
> characteristic step of the same time step
> ([slGradientControlSource.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/source/slGradientControlSource.C#L140-L200)).
> The technical report names the discretisation of this step, not the laws, as the primary
> mechanism of the coupled failures ([[concepts/source-discretisation-defects]]).

## What it is

| member | dictionary word | status | one line |
|---|---|---|---|
| none | `none` | settled | the identity; the production method |
| gradient control | `gradientControl` | candidate | $q$ from `fvc::grad(psi, "gradPsiSource")` (a CENTRED `leastSquares` gradient in every case template), $\mathbf{D}$ from `fvc::grad(U, "gradUSource")` only when the law needs it, the band test $\lvert\psi_c\rvert \le \texttt{bandCells}\,h_c\,q_c$ (a hard switch), $F = $ law$(q^2, \sigma, a)$, the clamp $\lvert\Delta t F\rvert \le 30$, then the multiplication |

Tokens: `SL_SOURCE none|gradientControl`, `SL_SOURCE_BAND_CELLS 3`, and the law tokens. The step
runs inside `slAdvection::advect` after the characteristic update and before the narrow band, the
phase indicator and the curvature are rebuilt
([slAdvection.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/slAdvection.C#L95-L115)). In the coupled solver
with `psiOuterCorrectors yes` it runs once per outer pass on the restored $\psi^n$, so it does not
compound ([slAlphaEqn.H](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L145-L165)).

## Why it matters

In the continuum the update keeps the sign of $\psi$ and the zero set does not move. On the mesh
it moves the interpolated zero crossing wherever $F$ differs between the two cells that hold it
([[concepts/sl-source-step]]), and, more important, the centred gradient makes the source's own
map linearly unstable ([[concepts/source-discretisation-defects]]).

## Where in the code

`slSource` is the runtime-selectable base, `slGradientControlSource` the model; the extension point
was opened in Phase C with 53 cases bit-identical ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3318-L3370)); the
model has 49 unit checks, and its halo mutant fails on four ranks, which shows the check sees the
processor patches ([STATUS 11.9](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3403-L3436)).

## Evidence

| claim | number | where |
|---|---|---|
| S1 (soft wall alone), shear arm, shape error at $N = 136$ | 0.254 against 6.19e-4, order 0.06 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023) MEASURED |
| HL1q, shear arm, band gradient error ratio to the baseline / shape error ratio | 0.204 / 173 | same, MEASURED |
| HL1q and HL1z, stationary droplet at $N = 200$ | shape error 3.8 R and 3.7 R, volume error 76 % and 84 %; the runs COMPLETE | same, MEASURED |
| growth rates on the stationary droplet (HL1q): curvature error, spurious current | about 290 1/s and 320 to 350 1/s over 10 to 20 ms, against $\mu = 270$ 1/s | archive `gate/histories_stationary.csv`, MEASURED |
| the clamp is never reached | HL1q at $q = 190$: $\Delta t F = -0.2$ | technical report, MEASURED |
| coupled seam check (translatingSeamNp1) | S1 0.76 FAIL; HL1q 4.9e-7, HL1z 5.3e-7 PASS | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023) MEASURED |

## Why it failed, or why we think so

See [[concepts/source-discretisation-defects]] (the centred-q eigenmode at rate $\mu$, the hard band
edge under an accumulated O(1) factor, the soft-wall slope) and [[concepts/why-the-candidates-failed]].

## Decisions

Inert by default; no candidate is production. The repair order is in [[concepts/gradient-control-next-experiments]]:
a one-sided $q$, a smooth band taper, a source-CFL guard, then a zero-flow static gate before any
coupled run.

## Open questions

Whether $F$ evaluated at the algebraic foot and copied along the normal (E2.4) is still needed
once the discretisation is repaired.

## Related

[[hubs/gradient-control]], [[models/gradient-control-law]], [[concepts/sl-source-step]],
[[concepts/source-discretisation-defects]], [[concepts/sdpls-source-eulerian]],
[[studies/method-gate-2d-campaign-2026-09]].

## Log

### 2026-09-28
Created.
