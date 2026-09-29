---
title: "Surface tension"
description: "Hub: the balanced-force capillary flux from the reconstructed curvature, the parasitic-current mechanism (source: the curvature estimator; amplifiers: translation and density ratio), and what is open"
kind: hub
status: settled
part: surface-tension
tags: [hub, part/surface-tension]
date: 2026-09-28
---
# Surface tension

[[index]] <- back

## The question this part answers

How is the capillary force built so that a droplet at rest stays at rest to round-off, and why does the stationary, translating and oscillating droplet develop a spurious current that grows?

## Current verdict (2026-09-28)

The force is a balanced-force CSF flux in the pressure-flux space ([[concepts/balanced-force-csf-flux]]), with the curvature from the quadratic fit of $\psi$ ([[concepts/curvature-from-the-fit]]) delivered by `cellCentreInverse` with the Gaussian-curvature-aware inverse ([[concepts/cell-centre-inverse-curvature]], [[decisions/curvature-extension-cell-centre-inverse]], [[decisions/curvature-inverse-gaussian]]); the Popinet translating family uses `none` (case-dependent). The pressure projection absorbs 99.85 to 99.98 % of the delivered force ([SL article sec. 2200](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2200)); with exact curvature the velocity is zero to $3.8\times10^{-10}$ in every arm ([[concepts/well-balanced-exact-curvature-gate]]). The parasitic current has one source and two amplifiers ([[concepts/parasitic-current-mechanism]]): the source is the curvature estimator (the step-1 kick is independent of $U_0$ to 0.6 %; exact curvature removes it $2\times10^4$ to $7\times10^4$ times), the amplifiers are translation and the density ratio, $\max|U|(T) = u_0(h)\exp G(h)$ ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L74-L129)). No filter, smoothing or relaxation is allowed in production ([[decisions/psi-filter-none]]); the capillary time step is 0.2323 of the Brackbill limit ([[concepts/capillary-time-step]]); the force is centred at the end of the step with BDF2 momentum ([[concepts/force-time-centring]], [[decisions/momentum-schemes-bdf2-upwind]]). Two limits of the production method stand: the curvature of the MOVING interface does not converge (2.4 to 3.6 % of $1/R$ at every rung of the translating arm, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)), and the oscillating droplet's gradient drift turns unstable at $N = 200$ within ten periods.

## Map

| kind | note | status | one line |
|---|---|---|---|
| model | [[models/surface-tension-force]] | settled | the twelve force models; reconstructedCurvature in production |
| model | [[models/curvature-extension]] | settled | the delivery words; cellCentreInverse in production |
| model | [[models/semi-implicit-capillary-force]] | open | the fvOption; the collective deadlock; not needed with projectedFlux |
| concept | [[concepts/balanced-force-csf-flux]] | settled | forces in flux space; what the projection absorbs |
| concept | [[concepts/curvature-from-the-fit]] | settled | symbolic curvature, the parallel-curve offset correction |
| concept | [[concepts/cell-centre-inverse-curvature]] | settled | +2.01 non-gradient content; second order on constant curvature only |
| concept | [[concepts/face-curvature-deliveries]] | retracted | gain against accuracy; the ellipse gate |
| concept | [[concepts/parasitic-current-mechanism]] | settled | source and amplifiers; the two-factor law |
| concept | [[concepts/curvature-corrugation-and-the-fit]] | settled | grid-scale modes; the filter as instrument |
| concept | [[concepts/integral-surface-tension-cst]] | retracted | better static balance, higher dynamic gain |
| concept | [[concepts/kang-gfm-and-sharp-heaviside]] | retracted | 58x better statics, earlier blow-up |
| concept | [[concepts/force-time-centring]] | settled | endStep; midpoint diverges 32 % earlier |
| concept | [[concepts/capillary-time-step]] | settled | 0.2323 Brackbill; the dt sweep |
| concept | [[concepts/pressure-projection-and-linear-solvers]] | settled | the projection gates; the non-orthogonal caveat |
| concept | [[concepts/variational-capillary-force]] | candidate | the proposal, no code |
| concept | [[concepts/well-balanced-exact-curvature-gate]] | settled | zero velocity with exact curvature |
| decision | [[decisions/surface-tension-reconstructed-curvature]] | settled | SURFACE_TENSION |
| decision | [[decisions/curvature-extension-cell-centre-inverse]] | settled | CURVATURE_EXTENSION, case-dependent |
| decision | [[decisions/curvature-inverse-gaussian]] | settled | the K-aware inverse |
| decision | [[decisions/momentum-schemes-bdf2-upwind]] | settled | BDF2, upwind |
| retraction | [[retractions/cell-mean-delivery-adoption]] | retracted | the ellipse gate collapsed it to first order |
| retraction | [[retractions/force-at-n-not-n-plus-1]] | retracted | the force is built from $\psi^{n+1}$ |
| retraction | [[retractions/t-blow-baseline]] | retracted | $t_\mathrm{blow}$ is not a proxy |
| retraction | [[retractions/psi-filter-seam-bug]] | voided | filtered results before 2026-08-19 |
| case | [[cases/stationary-droplet]], [[cases/oscillating-droplet]], [[cases/curvature-static-gates]] | settled | the cases |
| study | [[studies/curvature-stabilization-campaign]], [[studies/shannon-parasitic-currents-campaign]], [[studies/poly3d-roadmap]] | settled | the campaigns |

## Open, in order

0. The production curvature `cellCentreInverse` was never scored on the varying-curvature ellipse gate, so the acceptance criterion (`G h^2 <= 0.65`, order `>= 1.9` on the ellipse) is not demonstrated for it ([METHOD 8.1 L393](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393), [[cases/curvature-static-gates]]). The cheapest open item: a static gate of minutes.
1. Anti-convergence of the settled current under refinement ([METHOD 9.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L754-L756)); the two-factor law names the amplifier $G(h) \sim h^{-3.27}$ in 3D.
2. The curvature of the moving interface does not converge on the translating arm ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)).
3. The oscillating droplet's gradient drift at $N = 200$: the case for the gradient-control candidates ([[hubs/gradient-control]]).
4. The skew-mesh residual with exact curvature, insensitive to every lever but the `fvc::reconstruct` of the velocity correction ([METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L795-L801)).
5. The semi-implicit capillary force, untested on the translating droplet since its collective fix; the 3D offset correction; a torus gate for the K-aware inverse; the variational force.

## Log

### 2026-09-28
Created.

### 2026-09-29
Added open item 0 (the ellipse gate was never run on the production curvature).
