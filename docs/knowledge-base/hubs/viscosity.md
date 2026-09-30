---
title: "Viscosity"
description: "Hub: the face viscosity of the two-phase momentum equation; alg_lin is the decided model in 2D, 3D is open"
kind: hub
status: settled
part: viscosity
tags: [hub, part/viscosity]
date: 2026-09-28
---
# Viscosity

[[index]] <- back

## The question this part answers

Which face viscosity $\mu_f$ of the two-phase momentum equation is consistent with the interface representation, so that the viscous term converges across the density and viscosity jump and does not feed the spurious current?

## Current verdict (2026-09-28)

`VISCOSITY_FACE_MODEL alg_lin` (the algebraic phase indicator interpolated linearly to the face) is the decision of 2026-09-03 on the 36-arm `mufGrid2D` ladder: the only model with a positive order in all three metrics at both viscosity jumps, $L^1$ error $1.400\times10^{-3}$ at order 1.10 at ratio 1000 ([cases/default.parameter](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L701-L763), [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L397)). The harmonic and blended members are worse everywhere ([[models/viscosity-face-model]]). A switch to `geo_lin` was RETRACTED: the algebraic arms had a frozen $\mu_f$ (bug 39e59b3). Popinet's benchmark uses the same face properties ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L139-L143)). The Eulerian two-phase solver froze $\mu_f$ at its $t = 0$ value until c094bd8 ([[concepts/eulerian-solver-mass-flux-port]]). The measured role of viscosity in the parasitic current: it spreads the current rather than damping it ([[concepts/parasitic-current-mechanism]]).

## Map

| kind | note | status | one line |
|---|---|---|---|
| model | [[models/viscosity-face-model]] | settled | the six-name grid alg/geo x lin/harm/blend, the ladder, the retracted switch |
| concept | [[concepts/viscosity-open-items]] | open | the 3D token, the banner default, the sharp geometric $\alpha_f$ |
| decision | [[decisions/viscosity-face-model-alg-lin]] | settled | 2D decided; 3D open |
| case | [[cases/stationary-droplet]] | settled | the case of the ladder |

## Open, in order

1. The 3D case templates carry no `VISCOSITY_FACE_MODEL` token ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3545-L3623)).
2. `mufGrid2D` ran on 4, 8 and 16 ranks (N = 128, 256, 512) before the coupled-face density fix of 2026-09-27, with `CURVATURE_EXTENSION none` to 0.02 s ([[concepts/coupled-face-density-defect]]); the decision holds until a re-run says otherwise.
3. A geometric sharp $\alpha_f$ consistent with the interface plane, untested.

## Log

### 2026-09-28
Created.

### 2026-09-29
CORRECTED the rank counts of the mufGrid2D ladder (4, 8 and 16, not 8 only; concepts/viscosity-open-items).
