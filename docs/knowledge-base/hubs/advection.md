---
title: "Interface advection and the phase indicator"
description: "Hub: how the level set is transported (quadratic semi-Lagrangian, production) and how the phase indicator is built from it; what is settled, what failed, what is open"
kind: hub
status: settled
part: advection
tags: [hub, part/advection]
date: 2026-09-28
---
# Interface advection and the phase indicator

[[index]] <- back

## The question this part answers

How is the level set $\psi$ transported without a reinitialisation, at which order, on which meshes, and how is the phase indicator $\alpha$ built from $\psi$ so that the two-phase solver sees a sharp interface?

## Current verdict (2026-09-28)

The production transport is the quadratic semi-Lagrangian scheme: a constant-free quadratic weighted-least-squares VALUE fit ([[models/sl-reconstruction]], `uncachedQuadraticWeightedLeastSquares`, normal equations, pivot tolerance 0.3), a departure-centred AB2 foot ([[concepts/departure-foot-ab2-centring]]), the trace velocity from the reconstructed projected flux ([[concepts/trace-velocity-projected-flux]]), no clip and no value bound ([[models/sl-value-bound]], [[decisions/sl-clip-and-value-bound-off]]), no filter ([[decisions/psi-filter-none]]) and the Detrixhe-Aslam indicator ([[models/phase-indicator]]). Measured orders on the kinematic cases: shape 2.97 and volume 2.59 on the 2D vortex, 2.95 and 3.28 on the 3D shear, 1.36 and 1.46 on the 3D deformation ([METHOD 8](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347)); the fixed 2D gate's shear arm converges at order 3.07 in the shape and 3.60 in the volume ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)). Every alternative transport is measured and closed or ranked: the Eulerian FV transport is a robust second at 20 times the error ([[concepts/eulerian-fv-transport]]), the linear SL fits reach order about 1.1 ([[concepts/linear-semi-lagrangian]]), the normal-projected SL is closed ([[concepts/normal-projected-sl]]), geometric redistancing is closed ([[concepts/redistancing-geometric-grl]], [[models/redistancer]]), and every clip or value bound is falsified as a fix ([[concepts/value-bounds-and-clips]]). The scheme has no maximum principle: the reconstruct-and-evaluate operator amplifies the grid-scale mode on every mesh, $\rho(B) = 1.0044$ on hexahedra and $1.0111$ on pMesh polyhedra ([METHOD 9.3](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L771-L781)), and the polyhedral fit carries an amplification bound $\Lambda = 1.26$ against $1.05$ on hexahedra ([[concepts/polyhedral-fit-amplification]]).

## Map

| kind | note | status | one line |
|---|---|---|---|
| model | [[models/sl-reconstruction]] | settled | the value fit: why a value fit, why degree two, the members and their orders |
| model | [[models/sl-scheme]] | settled | pointValue (production), fluxForm, normalProjected (closed) |
| model | [[models/sl-value-bound]] | settled | none in production; stencilBounds and lipschitzCone falsified |
| model | [[models/level-set-advection]] | settled | eulerian against semiLagrangian in the kinematic solver |
| model | [[models/phase-indicator]] | settled | detrixheAslam in production; the first-order offset is open |
| model | [[models/narrow-band]] | settled | which band is used where; the dilation incident |
| model | [[models/redistancer]] | retracted | the line is closed: PDE reinit divergent, planeFootWave second order static only |
| model | [[models/volume-correction]] | settled | the global correction is a crossover, not a fix |
| concept | [[concepts/departure-foot-ab2-centring]] | settled | the arrival form leaves $+\Delta t^2 \partial_t u$; taylor equals rk2 |
| concept | [[concepts/polyhedral-fit-amplification]] | settled | $\Lambda$ 1.26 on pMesh; not a time-step effect |
| concept | [[concepts/idec-defect-correction-failure]] | retracted | iDEC diverges, $\rho > 1$ |
| concept | [[concepts/trace-velocity-projected-flux]] | settled | the reconstruct operator carries 70 % of the projectedFlux win |
| concept | [[concepts/psi-outer-correctors]] | settled | re-advect in every outer iteration; default yes |
| concept | [[concepts/static-local-refinement]] | settled | 3.4 to 6.9 times fewer core-hours at equal orders |
| concept | [[concepts/linear-semi-lagrangian]] | settled | the linear line at order about 1.1 |
| concept | [[concepts/normal-projected-sl]] | retracted | trace clean, write-back diverges; closed |
| concept | [[concepts/eulerian-fv-transport]] | settled | the limiter drops the order 3.0 to 0.9; flux form loses 17 % volume |
| concept | [[concepts/redistancing-geometric-grl]] | retracted | closed line |
| concept | [[concepts/value-bounds-and-clips]] | retracted | what the clips and bounds did, and the mesh-noise floor |
| decision | [[decisions/sl-reconstruction-uncached-qwls]] | settled | SL_RECONSTRUCTION |
| decision | [[decisions/sl-fit-normal-equations]] | settled | SL_FIT: QR is bit-identical and blows up identically |
| decision | [[decisions/sl-trace-velocity-projected-flux]] | settled | SL_TRACE_VELOCITY projectedFlux |
| decision | [[decisions/sl-clip-and-value-bound-off]] | settled | SL_CLIP false, SL_VALUE_BOUND none |
| decision | [[decisions/psi-filter-none]] | settled | no filtering in production |
| decision | [[decisions/phase-indicator-detrixhe-aslam]] | settled | PHASE_INDICATOR detrixheAslam |
| retraction | [[retractions/distance-cone-bound-as-transport-bound]] | retracted | the lipschitzCone gain was a coarse-mesh artefact |
| retraction | [[retractions/clip-damage-is-the-narrow-band]] | retracted | the clip fires in stencil-extremum cells |
| retraction | [[retractions/gradu-coupled-patch-contamination]] | voided | 31 parallel kinematic studies before 2026-08-26 |
| retraction | [[retractions/advection-orders-3-2-factor]] | retracted | the 2D orders of METHOD 8.3.7 were 3/2 too high |
| case | [[cases/kinematic-advection-cases]] | settled | the vortex, translation, shear, deformation, rotation cases |
| study | [[studies/sl-quadratic-pre-print]], [[studies/sl-linear-pre-print]], [[studies/grl-pre-print]], [[studies/npsl-design]], [[studies/method-comparison]] | settled | the documents |
| concept | [[concepts/advection-regression-set]] | settled | the standing hex 2D, hex 3D, poly 3D set |

## Open, in order

1. The amplification $\rho(B) > 1$ of the reconstruct-and-evaluate operator on every mesh: the acceptance criterion $\rho \le 1$ is met nowhere ([METHOD 9.3 and 9.4](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L771-L790)); a quadratic-exact rule with $\Lambda = 1$ is impossible at an off-node point.
2. The gradient drift without redistancing: the oscillating droplet's band gradient error grows exponentially at $N = 200$ after $t = 0.05$ s ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)); [[hubs/gradient-control]] holds the candidates.
3. The translation `none` arm saturates at $N = 256$ ([METHOD 8](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L648-L653)).
4. The signed-distance assumption survives in the phase indicator's first-order offset, unmeasured ([[models/phase-indicator]]).
5. The curated `advConv2D*_convergence.csv` files still carry the 3/2 orders; regenerate on Lichtenberg ([STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3260)).

## Log

### 2026-09-28
Created from STATUS sections 0, 1 and 11, METHOD sections 2, 8 and 9, and the SL article.
