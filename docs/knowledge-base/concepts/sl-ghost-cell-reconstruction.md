---
title: "Ghost cells for the semi-Lagrangian fit at a 2:1 transition"
description: "The original quadratic fit is destroyed when the interface crosses a static hexRef8 2:1 transition (shape error 5.6e-3 to 5.1e-1 over N = 32 to 256); the ghost-cell class keeps the error of the uniform coarse mesh, shape order 1.97 over four rungs."
aliases: [C-Cp/L2+g, ghostCellQuadraticWeightedLeastSquares, ghost-cell reconstruction, refined-box gate]
kind: concept
status: open
part: advection
tags: [concept, part/advection]
date: 2026-11-13
decided_by: [config/gates/ghostCellRefinedBox3D.yaml, config/gates/ghostCellRefinedBoxUniform3D.yaml, config/gates/ghostCellRefinedBoxN256Ghost3D.yaml, applications/test/leiaTestGhostCellReconstruction]
code: [src/leiaLevelSet/semiLagrangian/ghostCellQuadraticWeightedLeastSquaresReconstruction.C, src/leiaLevelSet/semiLagrangian/ghostCellQuadraticWeightedLeastSquaresReconstruction.H]
sources: ["STATUS 11.32"]
---
# Ghost cells for the semi-Lagrangian fit at a 2:1 transition

> Verdict (2026-11-13, MEASURED). On a statically refined hexahedral mesh the original scheme ([[models/sl-reconstruction]], `uncachedQuadraticWeightedLeastSquares`) does not survive an interface that crosses a 2:1 transition: on the refined-box translation its shape error grows from 5.588e-3 at N = 32 to 5.059e-1 at N = 256. The class `ghostCellQuadraticWeightedLeastSquares` completes all four rungs with at most 0.89 times the shape error and 1.27 times the volume error of the uniform coarse mesh, and a shape order of 1.97. The pre-registered gate G2 passes; the class is a candidate, not a default.

## What it is

Every cell fits its quadratic on a uniform stencil at its own cell size. A fine cell next to a coarser cell gets ghost data at its own lattice positions, evaluated from the reconstruction of the coarse cell (prolongation). A coarse cell next to finer cells gets one datum per refined parent, the mean of the children (restriction). The fits run level by level from the coarsest, and the coefficients cross processor patches through the stencil exchange. The class reads `cellLevel` and `pointLevel` of hexRef8 and is inert on a mesh without them. The restriction is not exact for a quadratic (error `h_fine^2/8` times the Laplacian), which the exact-solution gate measures as 2.7e-6 at a transition.

## Why it matters

The goal is an adaptive solver. The matrix symbol of the original scheme predicts a growth factor of 1.096 per step at a planar hexahedral 2:1 transition at Co 0.5 and at most 1.0044 with ghost cells ([[concepts/sl-amplification-symbol]], analysis of 2026-11-12 in `docs/analysis-amplification-sl`). This gate is the first solver test of that prediction, and it agrees in sense: the original scheme is destroyed at the transition, ghost cells are not. This is a different defect from the polyhedral one of [[concepts/polyhedral-fit-amplification]].

## The numbers

All numbers and the limits of the gate are in [STATUS 11.32](https://github.com/leia-openfoam/leia/blob/feature/sl-ghost-cells/STATUS.md). Shape error (L1) at T = 1 over N = 32, 64, 128, 256: original 5.588e-3, 3.736e-2, 2.512e-1, 5.059e-1; ghost cells 5.695e-3, 2.287e-3, 3.488e-4, 1.122e-4; uniform coarse control 6.403e-3, 2.766e-3, 5.687e-4, 1.361e-4. At N = 64 and T = 1, 71 % of the summed level-set error of the original scheme lies within 3 coarse cells of a transition. Cost of the class: 1.4 to 1.7 times the solve wall clock, and more than 4.9 GB per rank at 10.2e6 cells on 8 ranks against 2.7 GB.

## Open

1. A serial against np 4 equivalence run of the solver arm (the unit gate ran both).
2. A restriction that is exact for a quadratic.
3. More than one refinement level, refinement that follows the interface, 2D, the coupled solvers.
4. Memory and cost of the class.
5. The advection regression set ([[concepts/advection-regression-set]]) with the class on uniform meshes, where it must equal the original scheme.

## Log

- 2026-11-13 · gate `ghostCellRefinedBox3D` on Lichtenberg: G1 holds, G2 passes, the falsification condition is not met (STATUS 11.32).

[[hubs/advection]] - [[models/sl-reconstruction]] - [[concepts/sl-amplification-symbol]] - [[concepts/polyhedral-fit-amplification]] - [[concepts/advection-regression-set]]
