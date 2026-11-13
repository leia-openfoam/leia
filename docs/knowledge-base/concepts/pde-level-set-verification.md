---
title: "The PDE level-set direction (Bothe and Soga): verify the parts one at a time"
description: "Open 2026-11-12 - the model of D. Bothe and K. Soga (velocity extension plus a gradient source) on the Eulerian FVM line, split into sub-problems with exact references; an audit of the same day finds the earlier velocity-extension verdicts weak or void (reversed flows read at T, two changes at once, a far-field flux replacement in the code, a band-edge jump) and re-opens them"
aliases: [PDE level set, Bothe-Soga model, velocity extension verification, VE verification]
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-11-12
decided_by: author decision 2026-11-12
code: [src/leiaLevelSet/velocityExtension/, src/leiaLevelSet/advection/eulerianAdvection.C, applications/test/leiaTestVelocityExtension/]
sources: [docs/plan-velocity-extension-verification.md, docs/pde-level-set/pde-level-set-presentation/pde-level-set.template.html, STATUS 11.29, D. Bothe's meeting notes of 2026-10-29]
---
# The PDE level-set direction (Bothe and Soga): verify the parts one at a time

> Open, 2026-11-12, by the author's decision. The goal is a PDE level-set method on the Eulerian
> FVM line; since the third round of 2026-11-12 the goal is an accurate, resolved, adaptive solver,
> judged against exact solutions and observed orders (interFoam is no reference). The model is
> $\partial_t\phi + \mathbf u_e\cdot\nabla\phi = \beta\,\phi\,g(\lvert\nabla\phi\rvert)$ with the extended velocity
> $\mathbf u_e = \mathbf u(\Pi(x))$ ([plan, section 1](https://github.com/leia-openfoam/leia/blob/development/docs/plan-velocity-extension-verification.md)).
> The work verifies five parts one at a time, each against an exact reference: the transport
> (step 0), the extension operator (A0 repairs, A static), the transport with the exact and with the
> computed extension (B, C), the source (D), the assembly and the droplets against exact solutions, then on refined meshes (E).
> An audit of the same day re-opens the earlier extension verdicts (below).
> The deck for the supervisor: `docs/pde-level-set/pde-level-set-presentation/`.

## What it is

The direction of the dossiers of D. Bothe (2026-09-26) and of D. Bothe, T. Marić and K. Soga
(2026-09), and of D. Bothe's notes of 2026-10-29. The notes rank three routes: (0) an algebraic
improvement of $\phi$ ($\phi/\lvert\nabla\phi\rvert$, the second-order distance $\hat d$), (1) the velocity
extension with the Newton point $P_x = x - \phi\nabla\phi/\lvert\nabla\phi\rvert^2$, a cap and a blend
$\delta\mathbf v(x) + (1-\delta)\mathbf v(P_x)$, (2) the source $\beta\phi(1 - \lvert\nabla\phi\rvert)$. The
semi-Lagrangian line is out of scope for this direction.

## Why it matters

Pure transport changes $q = \lvert\nabla\phi\rvert$ by the normal strain, $D_t q = -a q$. The extension
removes the cause; the source removes the rest. Every earlier attempt measured several parts at once
([[hubs/gradient-control]]), so no failure could be located.

## Evidence

| claim | number | where |
|---|---|---|
| no extension was ever compared with the exact extension $\mathbf u(\Pi(x))$; the static test measures only $\lvert\hat{\mathbf n}\cdot\nabla\mathbf u_e\rvert$ | — | [plan, section 12](https://github.com/leia-openfoam/leia/blob/development/docs/plan-velocity-extension-verification.md), CONFIRMED in the code |
| every model except `none` and `haloLimited` transports $\phi$ with `interpolate(Uext) & Sf` on every face, so the far-field flux is replaced too | — | [`interfaceExtension.C`](https://github.com/leia-openfoam/leia/blob/456cae4f/src/leiaLevelSet/velocityExtension/interfaceExtension.C#L435-L468), [`projectionBlend.C`](https://github.com/leia-openfoam/leia/blob/456cae4f/src/leiaLevelSet/velocityExtension/projectionBlend.C#L172), CONFIRMED |
| `closestPoint` and `meshWaveExt` write the extension only inside the band; the fade blends U with U, so the velocity jumps at the band edge | — | same, CONFIRMED |
| at T/2 the extensions win in every data set found; the losses are read at T of reversed flows | bench N = 512: band gradient 0.87 against 1.67 for `none` | [plan, 12.1](https://github.com/leia-openfoam/leia/blob/development/docs/plan-velocity-extension-verification.md), MEASURED (old data) |
| the Eulerian baseline is not verified: its band gradient error does not converge; serial and np 4 differ at step 1 | orders -0.3 to 0.0; 1.6e-6 | STATUS 11.28 items 5 and 9, MEASURED |

## Decisions

The advected verdicts on velocity extension are re-opened: they are kept as verdicts on a
configuration and measured again after the A0 repairs ([[decision-log]]). Step 0 (the transport) and
A0 (the repairs) come first.

## Open questions

The four extension variants a) to d) of D. Bothe's notes; what sets $\beta$; whether the theorem of
the velocity-extension paper holds with the blend; which range of $q$ is enough for multiphase flow
(plan, section 9 idea 4).

## Related

[[hubs/gradient-control]], [[models/velocity-extension]], [[concepts/projection-blend-extension]],
[[concepts/closest-point-extension]], [[concepts/halo-limited-extension]],
[[concepts/algebraic-rescaling-eulerian]], [[concepts/static-source-gate-eulerian]],
[[concepts/combined-source-note]].

## Log

### 2026-11-12
Created with the plan, the deck and the audit (STATUS 11.29).
The supervisor's second round: the reading of D. Bothe's notes is right; the q range waits for the static gate of plan step 0.5 (the function-reading operators are invariant under a uniform scaling of psi, the distance-reading ones are not); step 0 starts on Lichtenberg.
The supervisor's third round: forget interFoam; an accurate, resolved, adaptive solver first. Plan rows 1 and 10, step E and the deck are updated.

### 2026-11-12 (results)
Step 0.1 PASS (the transport keeps a distance to second order and better); the flux confound is real and small (volume 1 to 4 %); closestPoint is the only second-order extension against u(Pi(x)), but its face flux is first order; the seam difference of the baseline goes to a tolerance test (STATUS 11.30).

### 2026-11-13
The Eulerian baseline's serial/np 4 difference is solver noise (it falls by 2e4 from psi tolerance 1e-10 to 1e-14); closestPoint's first-order face flux is the band-edge jump, the flux is second order inside the band (STATUS 11.30 items 8 and 9).
The supervisor's go for A0 repair 2 (the band edge) and then step A on q != 1, the ellipse and 3D; jobs on project special00005. `fadePlacement` (`legacy` default, `extendThroughFade`, `insideBand`) and the pre-registered gates `aStatic2Dedge`, `a2Seam2Dvortex` (STATUS 11.31).
A0.2 results (STATUS 11.31 items 4 to 8): the default is inert (36 of 36 identical); extendThroughFade removes the jump inside L (4h-shell flux error 24 times smaller, order 2 instead of 1) and makes closestPoint 330 times less decomposition-dependent; against the reference with the fade the order stays 1 (falsified as pre-registered; derived cause: a fade whose width scales with h is first order by itself); on np 4, 7 % of the wider band falls back to the steadyUpwind fill (the one-layer halo), which repair 3 must cover.
