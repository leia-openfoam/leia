---
title: "The closest-point extension (FP0)"
description: "The reference extension that samples the velocity at the closest interface point: the same order loss as HL0 in the shear flow, decomposition-dependent (45 %), and a collapse at the shear tail"
kind: concept
status: retracted
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/velocityExtension/closestPoint.C, src/leiaLevelSet/velocityExtension/interfaceExtension.C]
sources: [STATUS 11.15, MC article, SDPLS article]
---
# The closest-point extension (FP0)

> Retracted as a production candidate (2026-09-27); kept as a reference arm. FP0 (`closestPoint`)
> extends $\mathbf{u}$ by its value at the closest interface point, with a steady solve as the
> fallback where the search fails, and traces with `interpolate(Uext) . Sf` (no shell). In the
> gate it loses the same order as HL0 in the shear flow (1.38; 34 times the shape error), its
> kinematic seam checks FAIL (0.45 and 0.48: the closest-point search depends on the decomposition
> at the processor boundaries), its coupled seam check fails at 1.11, and its stationary and
> oscillating $N = 200$ runs timed out at the 10 h limit
> ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)).

## Evidence

| claim | number | where |
|---|---|---|
| shear arm shape error at $N = 136$, order | 2.12e-2 against 6.19e-4; 1.38 | archive verdict, MEASURED |
| shear band gradient ratio (target 0.8) | 0.939 | same, MEASURED |
| exact-1D `qBandMean` | 0.96 to 1.09 (the closest-point extension keeps $q$ near one in 1D) | archive, MEASURED |
| the shear tail | band error 0.001 to 0.03 until $t = 0.1$ s, then 0.66 with volume error 3.6e-2 by 0.4 s | archive histories, MEASURED |
| method comparison: VE closestPoint on pure advection | 1.44e-3 at 1.61e4 s against SL 1.49e-5 at 595 s: dominated | [MC article](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L124) MEASURED |
| "closestPoint 2185x worse" | DISPUTED (serial re-measure 0.172 / 0.267) | [[studies/sdpls-pre-print]] |

## Why it failed, or why we think so

FP0 has no shell, so within-cell variation of a flux correction is not necessary for the order
loss; the loss comes from the tangential shear that every normal-constant extension puts into the
band and from the distance function's kinks entering the band at the shear tail (DERIVED;
[[concepts/extension-strain-relocation]]). The decomposition dependence is a defect of the search.

## Related

[[hubs/gradient-control]], [[models/velocity-extension]], [[concepts/halo-limited-extension]],
[[concepts/seam-checks-and-decomposition-invariance]].

## Log

### 2026-09-28
Created.
