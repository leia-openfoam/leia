---
title: "Gradient control: the direction, the dossiers, the plan, the pre-print, the report"
description: "Why leia tries to keep |grad psi| near one without redistancing, which documents define the direction (two dossiers outside git, the plan, the pre-print, the technical report), and the candidate hierarchy"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [docs/plan-halo-limited-gradient-control.md, gcls article, technical report]
---
# Gradient control: the direction, the dossiers, the plan, the pre-print, the report

> Open (2026-09-28). The level set of the production method is not reinitialised, so $q = \lvert\nabla\psi\rvert$
> drifts with the normal strain of the flow ([METHOD 9.2](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L757-L765)): in the band
> $q$ spreads to 0.84 to 1.37 under transport, and the oscillating droplet's drift turns
> unstable at $N = 200$ ([[cases/oscillating-droplet]]). Gradient control is the family of
> methods that keep $q$ near one without rebuilding the band from geometry: a source term
> $\psi F$ ([[models/sl-source]], [[models/sdpls-source]]) or a velocity extension
> ([[models/velocity-extension]]). The first campaign failed ([[studies/method-gate-2d-campaign-2026-09]]).

## What it is

The direction is defined by two technical dossiers of the authors (not in git, cited by title):
"Modified level set, technical dossier: halo-limited directional extension" (2026-09-26) and
"OpenFOAM FVM foot-point source dossier" (Bothe, Marić, Soga). The dossiers prescribe the
candidate hierarchy HL0 (extension alone), HL1 (extension plus a weak linear law), HL2
(extension plus the soft wall), FP0 (closest-point reference), ALG (the uncapped algebraic
sample with a strain taper) and S1 (the soft wall alone), then the routes S2 (a scalar memory
$Q$ that filters $q$) and V1 (a vector orientation memory) after the scalar routes. They
also prescribe Tier-I tests (uniform translation, rigid rotation, planar strain with
$q_0 = 0.92 / 1.00 / 1.08$, perturbed scaling) that were NOT run before the coupled gate; the
technical report makes them the first rung of the next campaign.

The executable plan is [plan-halo-limited-gradient-control.md](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md)
(approved 2026-09-26): the SL line first, one 2D and one 3D gate that test every method, the
dossiers stay out of git. The pre-print of 2026-09-27 ([[studies/gcls-pre-print]]) gives the
models, their discretisation and the gate results. The technical report of 2026-09-28
([gclsTechnicalReport.tex](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex))
assesses why each candidate failed and what to do next; the knowledge-base notes
[[concepts/why-the-candidates-failed]], [[concepts/source-discretisation-defects]],
[[concepts/extension-strain-relocation]] and [[concepts/gradient-control-next-experiments]] carry
its content.

## Why it matters

A working gradient control removes the last reason for a reinitialisation, which the earlier
lines showed to be injurious under transport ([[concepts/redistancing-geometric-grl]]). It must
satisfy the constraints of every leia method: unstructured FVM with compact stencils under MPI,
inert by default, no filter, judged on the moving-interface gates with the whole error vector
([CLAUDE.md, constraints](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L644-L659)).

## Where in the code

`libleiaGradientControl` (the laws), the SL source step, the SDPLS gradientControl source, the
halo-limited extension; the gate: `config/gates/methodGate2D.yaml`, `config/candidates/*.yaml`,
`workflow/Snakefile.gate`.

## Open questions

Everything listed in [[concepts/gradient-control-open-decisions]].

## Related

[[hubs/gradient-control]], [[concepts/sdpls-source-eulerian]], [[concepts/sl-source-step]],
[[concepts/halo-limited-extension]], [[concepts/combined-source-note]].

## Log

### 2026-09-28
Created.
