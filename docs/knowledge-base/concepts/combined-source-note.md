---
title: "The combined source F = a + lambda (1 - q): subsumed"
description: "The 2026-08-06 analysis note of the combined strain-plus-restoring source (logistic interfacial law, invariant interval) and its plan; subsumed 2026-09-26 by gradientControl linearQ with strain weight full"
kind: concept
status: retracted
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [docs/combined-source-terms/levelset_combined_source_note.tex, docs/plan-combined-source-terms.md]
---
# The combined source F = a + lambda (1 - q): subsumed

> SUBSUMED 2026-09-26. The note
> [levelset_combined_source_note.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/combined-source-terms/levelset_combined_source_note.tex)
> (2026-08-06, analysis only) shows that $F = a + \lambda(1 - \lvert\nabla\psi\rvert)$ gives a logistic
> law for $q$ on the interface and an invariant interval for $q$ in the band. Its plan
> [plan-combined-source-terms.md](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md) (WP2, WP3) and the brief
> `improvement-sdpls-combined.md` are subsumed: the combined source is the law `linearQ` with
> strain weight `full` of [[models/gradient-control-law]] ([STATUS 11.1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3228-L3234)).

## What survives

The continuum diagnosis of the SDPLS source, the noise-gain bound and the design window
$\lambda\Delta t \lesssim 1/(\pi W)$ ([plan-combined-source-terms](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L117-L165)),
and the dead-end list ([same](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L594-L634)). The WP0
$\delta_h$ diagnostic of the plan never ran.

## Why it is closed

No candidate of the first campaign used the strain weight (`none` everywhere), and the technical
report abandons the strain-reading route until the curvature-error seed is smaller
([[concepts/sdpls-source-eulerian]]).

## Related

[[hubs/gradient-control]], [[models/gradient-control-law]], [[concepts/gradient-control-overview]].

## Log

### 2026-09-28
Created.
