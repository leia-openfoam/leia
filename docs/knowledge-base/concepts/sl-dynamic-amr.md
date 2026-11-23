---
title: "Semi-Lagrangian level set on cubic cells with dynamic AMR"
description: "The AMR line of the pre-print, version 1: blastAMR with two levels and a uniform band of finest cells, the unit-cube fit operator, a quadratic mesh transfer; the first 4-rank runs complete and do not depend on the decomposition or on a redistribution."
aliases: [AMR line, SLLS AMR, cubicLatticeQuadraticWeightedLeastSquares, slMeshTransfer, slRefinementBand, amrBand]
kind: concept
status: open
part: advection
tags: [concept, part/advection]
date: 2026-11-23
decided_by: [config/slAmr/amrSmokeTranslation3D.yaml, config/slAmr/amrDecompSerial3D.yaml, config/slAmr/amrDecompNp4_3D.yaml, config/slAmr/amrBalance3D.yaml]
code: [src/leiaLevelSet/semiLagrangian/cubicLatticeQuadraticWeightedLeastSquaresReconstruction.C, src/leiaLevelSet/semiLagrangian/meshTransfer/slQuadraticFitMeshTransfer.C, src/leiaLevelSet/semiLagrangian/meshTransfer/slRefinementBand.C, applications/solvers/leiaSemiLagrangeLevelSetFoam/leiaSemiLagrangeLevelSetFoam.C]
sources: ["STATUS 11.44", "docs/plan-sl-amr-preprint.md"]
---
# Semi-Lagrangian level set on cubic cells with dynamic AMR

> Verdict (2026-11-23, MEASURED, one small case). The kinematic semi-Lagrangian solver runs on a mesh that blastAMR adapts in every step, with two refinement levels and load balancing, on 4 ranks. The result does not depend on the decomposition or on a redistribution: the differences are at round-off (shape error 1.9e-12 relative, volume 6.8e-16 of the droplet volume). No convergence study exists yet, so no accuracy statement is made.

## What it is

The line builds on the ghost-cell fit ([[concepts/sl-ghost-cell-reconstruction]]). Four new run-time-selectable parts, all inert by default:

- **The unit-cube operator.** On a cube with all 26 lattice positions filled, the anchored weighted least-squares fit is one 9 x 26 matrix, computed once. Cells with a boundary face in the stencil keep the fit per cell.
- **The mesh transfer** `quadraticFit`. A child cell gets the value of the parent's quadratic at its centre. A merged cell gets the mean of its children minus (h_c^2/8) tr(H). Both are exact for a quadratic field. OpenFOAM's default copies the parent value, which is first order.
- **The refinement band.** A field marks the cells within 6 layers of the interface: 3 for the face curvature, 1 for the fit stencil, 1 guard against ghost values, 1 for the travel at CFL at most 1. The count follows from the stencils (no free parameter); the curvature study C2 of the plan must confirm it. A buffer of 2 cells separates the levels.
- **The rescaling** as a source of the `slSource` family ([[concepts/algebraic-rescaling-eulerian]]), with a deterministic trigger. Written, not yet measured.

## The numbers

All numbers are in STATUS 11.44. Base N = 16, finest h = 1/64, T = 0.1, CFL 0.5, 4 ranks: 56 steps, 20 mesh changes, 46104 to 50472 cells, relative volume error 2.03e-5, relative shape error 3.73e-3. The matrix symbol with three cell sizes ([[concepts/sl-amplification-symbol]]) gives the buffer of 2 cells a growth factor - 1 of at most 2.4e-3 in 3D, below the criterion 4.41e-3.

## Open

- The bit-identity run of the inert composition root on the advection regression set.
- The unit gates of the transfer and of the unit-cube operator.
- The advection ladders against the uniform control, without and with the rescaling.
- The restriction of the ghost values is the plain mean of the children; a restriction that is exact for a quadratic is a roadmap candidate.

## Log

- 2026-11-23 · first 4-rank smoke passes; a cell with hanging nodes has no hexahedral shape, so the lattice comes from its 8 anchor points (4768 cells without a lattice before, 0 after) · STATUS 11.44
- 2026-11-23 · serial against 4 ranks and the forced load balancing (15 redistributions against 5): equal cell counts at every mesh change, differences at round-off; both read-outs fail as written, because a relative difference of the volume error is not bounded where the error crosses zero · STATUS 11.44
