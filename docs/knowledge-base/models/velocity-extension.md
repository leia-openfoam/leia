---
title: "velocityExtension: the eight extensions"
description: "none, anisotropicDiffusion, pseudoTime, steadyUpwind, steadyUpwindLinear, closestPoint, meshWave and haloLimited: what each does, and the verdicts (legacy models dominated on pure advection; closestPoint decomposition-dependent; haloLimited relocates the strain)"
kind: model
status: settled
part: gradient-control
tags: [model, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/velocityExtension/]
sources: [VE deck, MC article, roadmap, STATUS 11.10, STATUS 11.15]
---
# velocityExtension: the eight extensions

> Settled 2026-09-28: `none` is production. Every extension measured on pure advection has the
> worst accuracy at 12 to 27 times the cost ([plan-combined-source-terms](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L623-L625)),
> `closestPoint` is decomposition-dependent (seam FAIL 0.45 kinematic, 1.11 coupled), and the
> halo-limited extension of the first gradient-control campaign relocates the normal strain into
> a shell instead of removing it ([[concepts/extension-strain-relocation]]).

## What it is

An extension replaces the trace velocity of the level set near the interface by a field that is
constant along the normals ($\mathbf{n}\cdot\nabla\mathbf{u}_\mathrm{ext} = 0$), so that the level
sets move as parallel surfaces and $q$ stays one. The dictionary words are the `TypeName` strings
under [velocityExtension/](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension); the token is
`VELOCITY_EXTENSION`, and in the SL solvers `SL_TRACE_FLUX extension` traces with the extended flux.

| member | dictionary word | status | one line | evidence |
|---|---|---|---|---|
| none | `none` | settled | the identity; production | — |
| anisotropic diffusion | `anisotropicDiffusion` | retracted | a tangential-diffusion PDE; cannot converge in the static test | [VE deck](https://leia-openfoam.github.io/leia/decks/velocity-extension.html) |
| pseudo time | `pseudoTime` | retracted | iterates the extension equation in pseudo time; one-sided errors | VE deck |
| steady upwind | `steadyUpwind` | retracted | the steady extension equation with upwind; `UEXT_DIV upwind` because a deferred correction diverges on the steady equation ([default.parameter](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L16-L18)); translating survival 0.0388 s against 0.0103 s for none at $N = 32$ ([roadmap](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L121-L134)) | VE deck, roadmap |
| steady upwind, linear | `steadyUpwindLinear` | retracted | a divergent cascade | VE deck |
| closest point | `closestPoint` | retracted as production | samples $\mathbf{u}$ at the closest interface point; the fallback steady solve where the search fails; FP0 in the gate: shear shape error 34x at order 1.38, seam FAIL, stationary and oscillating $N = 200$ timed out | [[concepts/closest-point-extension]] |
| mesh wave | `meshWave` | retracted | a fast-marching-type propagation; survival 0.0131 s at $N = 32$ | roadmap |
| halo limited | `haloLimited` | retracted at $R = h$ | `radiusCells <= 1`; travel `capped`, weight `fractionReached`, direction `levelSet`, sampler `stencilFit`; unit gates to 2e-16; HL0 FAILS the gate | [[concepts/halo-limited-extension]] |

The intermediate base `interfaceExtension` (nLayers, nDescent, nAnchorLayers, projectFlux,
fadeMode, maxScale) serves the propagation models; `haloLimited` derives from the base class
directly ([haloLimited.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension/haloLimited/haloLimited.C)).

## Why it matters

The extension is the geometric alternative to a source term: it removes the cause (the normal
strain at the interface) instead of correcting the effect. Its limits are the halo (a wide
extension costs a halo per layer under MPI) and the tangential shear that every normal-constant
velocity puts into the band (level sets at distance $d$ rotate at a different rate;
[[concepts/extension-strain-relocation]]).

## Evidence

| claim | number | where |
|---|---|---|
| the legacy models on the reversed vortex ladders | worst accuracy at 12 to 27x cost | [plan-combined-source-terms](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L623-L625) MEASURED |
| method comparison decision table: VE closestPoint | 1.44e-3 at 1.61e4 s against SL 1.49e-5 at 595 s | [MC article](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L124) MEASURED |
| "closestPoint 2185x worse" | DISPUTED; serial re-measure 0.172 / 0.267 | [[studies/sdpls-pre-print]] |
| HL0 gate | shear 14x at order 1.39; oscillating gradient band 2e7; translating 64x | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023) MEASURED |

## Decisions

`none` ([[decisions/sl-trace-velocity-projected-flux]] holds the trace velocity decision that
made the extension unnecessary for stability: the reconstruct operator carries the win).

## Open questions

A mesh-resolved extension ($R \ge 3h$ with the dossier's taper; the ALG route) needs a sampler
with a 3-layer halo; the technical report puts it after the source repair and only for droplets.

## Related

[[hubs/gradient-control]], [[concepts/halo-limited-extension]], [[concepts/closest-point-extension]],
[[concepts/extension-strain-relocation]], [[concepts/trace-velocity-projected-flux]],
[[studies/velocity-extension-pre-print]].

## Log

### 2026-09-28
Created.
