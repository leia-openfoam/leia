---
title: "The combined source F = a + lambda (1 - q): re-opened"
description: "Re-opened 2026-10-13 by the author: the strain term a of the SDPLS source (D. Bothe) plus the restoring term lambda (1 - q) (K. Soga); executable as gradientControl linearQ with strain weight full; first the discretisation of q is repaired, then the kinematic ladders, then the coupled stationary droplet"
kind: concept
status: candidate
aliases: [combined source, Bothe and Soga source, strain plus restoring source]
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [docs/combined-source-terms/levelset_combined_source_note.tex, docs/plan-combined-source-terms.md]
---
# The combined source F = a + lambda (1 - q): re-opened

> RE-OPENED 2026-10-13 by the author's decision: the combined source is the first line of the
> gradient-control work. $F = a + \lambda(1 - q)$, with the strain term
> $a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n}$ of the SDPLS source (D. Bothe) and the restoring term
> $\lambda(1 - q)$ (K. Soga), $q = \lvert\nabla\psi\rvert$. Executable form: `SL_SOURCE
> gradientControl`, `GC_LAW linearQ` ($\lambda$ = `GC_MU`), `GC_STRAIN_WEIGHT full`
> ([[models/gradient-control-law]]). Step 1 repairs the discretisation of $q$: the centred $q$ fails at
> zero velocity ([[retractions/centred-q-eigenmode-growth-rate]]); `GC_Q_SCHEME upwindFit` is under
> test (gate E2.1, STATUS 11.26). Step 2: the kinematic ladders with the one-sided $q$ on the three
> rungs of the advection regression set. Step 3: the coupled stationary droplet at
> $M_\mu \in \{0.1, 1, 10\}$.

## Why the strain term carries the zero set and the restoring term the band

DERIVED. For $\psi_t + \mathbf{u}\cdot\nabla\psi = a\,\psi$ with $a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n}$,
$$D_t\lvert\nabla\psi\rvert = \psi\,\mathbf{n}\cdot\nabla a,$$
for any $\psi$: the stretching term $-\lvert\nabla\psi\rvert\,\mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n}$ of pure
transport cancels with $a\lvert\nabla\psi\rvert$. On the zero set ($\psi = 0$) the gradient is kept
exactly; off it, the drift is $\psi\,\mathbf{n}\cdot\nabla a$, of the order of the distance. The term
$\lambda(1 - q)$ removes that drift and the discretisation errors. Because the source multiplies
$\psi$, neither term moves the zero set in the continuum.

## The 2026-08-06 note and why it was subsumed

SUBSUMED 2026-09-26 (until 2026-10-13). The note
[levelset_combined_source_note.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/combined-source-terms/levelset_combined_source_note.tex)
(2026-08-06, analysis only) shows that $F = a + \lambda(1 - \lvert\nabla\psi\rvert)$ gives a logistic
law for $q$ on the interface and an invariant interval for $q$ in the band. Its plan
[plan-combined-source-terms.md](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md) (WP2, WP3) and the brief
`improvement-sdpls-combined.md` were subsumed: the combined source is the law `linearQ` with
strain weight `full` of [[models/gradient-control-law]] ([STATUS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3228-L3234)).

## What survives

The continuum diagnosis of the SDPLS source, the noise-gain bound and the design window
$\lambda\Delta t \lesssim 1/(\pi W)$ ([plan-combined-source-terms](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L117-L165)),
and the dead-end list ([same](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L594-L634)). The WP0
$\delta_h$ diagnostic of the plan never ran.

## Why it was closed (2026-09-28 to 2026-10-13)

No candidate of the first campaign used the strain weight (`none` everywhere), and the technical
report abandons the strain-reading route until the curvature-error seed is smaller
([[concepts/sdpls-source-eulerian]]).

## Related

[[hubs/gradient-control]], [[models/gradient-control-law]], [[concepts/gradient-control-overview]].

## Log

### 2026-09-28
Created.

### 2026-10-13
RE-OPENED by the author's decision: the combined source is the first line of the gradient-control
work. Status candidate. The derivation of $D_t\lvert\nabla\psi\rvert$ added.
Step 1 (the one-sided q, gate E2.1) PASSED. Step 2 (the kinematic arms of the 2D method gate)
FAILED for CS1 and both single terms on the shape order; the 2D cause is open (STATUS 11.26).
