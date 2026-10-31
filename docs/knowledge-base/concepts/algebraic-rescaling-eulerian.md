---
title: "The algebraic rescaling of psi on the UFVM line (D. Bothe, use A)"
description: "Settled 2026-10-31 - psi -> c psi with c > 0, first order psi/|grad psi| and second order (the quadratic root): applied once with the centred gradient it gives a second-order distance and keeps the zero set and the DA volume; applied every step it diverges (centred) or settles 1.1 h from d (one-sided); not an alternative to a source term on the UFVM line"
aliases: [rescaling, rescalingRedistancer, REDISTANCER rescaling, use A, (0a), (0b), dhat]
kind: concept
status: settled
outcome: negative
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-10-29
date_settled: 2026-10-31
decided_by: the static rescaling studies on Lichtenberg (STATUS 11.28 item 8)
code: [src/leiaLevelSet/redistancer/rescalingRedistancer.H, src/leiaLevelSet/profile/rescaledDistance.H, config/staticSource/staticSourceRescaleOnce2D.yaml, config/staticSource/staticSourceRescaleEvery2D.yaml, config/staticSource/staticSourceRescaleGauss2D.yaml]
sources: [STATUS 11.28, gcls-static-source-report sections 5 and 6, D. Bothe's notes of 2026-10]
---
# The algebraic rescaling of psi on the UFVM line (D. Bothe, use A)

> Settled 2026-10-31 (STATUS 11.28 item 8; report
> `docs/gradient-controlled-level-set/gcls-static-source-report/`, sections 5 and 6). MEASURED on
> the static circle of [[concepts/static-source-gate-eulerian]], uniform hexahedra, N = 64 to 512.
> Once, centred gradient, second order: exact for r^2 - R^2 (round-off), order 2 otherwise, the
> zero crossings and the DA volume of d kept. Every step: every variant fails. Centred (first or
> second order, leastSquares or Gauss Hessian) diverges from the first events; one-sided
> (Rouy-Tourin) stays bounded but settles at dist3_2 = 1.13 h on every rung and moves the zero
> crossings 0.24 to 0.41 h. The FVM Hessian is not what fails first: the centred gradient is.

## What it is

Two pointwise maps that multiply psi by a positive factor, so the sign of psi and the zero set stay
as they are (D. Bothe's notes, the author's decision of 2026-10-29):

- (0a) $r = \psi/|\nabla\psi|$, a first-order distance: $r = d - \frac{B}{2A}d^2 + O(d^3)$.
- (0b) $\hat d = 2\psi/(|\nabla\psi| + \sqrt{(|\nabla\psi|^2 - 2\psi\,\partial_{nn}\psi)^+})$,
  the stable root of the second-order Taylor expansion of $\psi(x - dn) = 0$: $\hat d = d + O(d^3)$.

"Use A" applies the map to psi as a redistancer event (once, or every step). "Use B" uses the
distance only inside the velocity extension ([[concepts/projection-blend-extension]]). The
semi-Lagrangian version with the quadratic fit is task T6 of [[sessions/sl-session-handover]].

## Why it matters

The static source test showed that every Eulerian source law grows grid-scale extrema at its G(0)
with the centred gradient. The rescaling was the candidate alternative: a map with no rate.

## Where in the code

`REDISTANCER rescaling` (`rescalingRedistancer`, the redistancer family of `leiaLevelSetFoam`),
the shared formula `src/leiaLevelSet/profile/rescaledDistance.H`, tokens `RESCALE_*` in
`cases/default.parameter`. Euler in pseudo-time (`PSI_DDT_SCHEME PsiDdtEuler`): Crank-Nicolson
re-applies a change made outside the solve.

## Evidence

| claim | number | where |
|---|---|---|
| once, centred, second order, from r^2 - R^2 | dist3_2 4.9e-15 to 4.5e-14 h (round-off) | STATUS 11.28 item 8; `rescale_verdicts.csv` R1 |
| once, centred, second order, other starts | order 1.93 to 2.01 | same |
| zero crossings and DA volume after one event | equal to those of d (difference 0.000) | R2, R3 |
| d a fixed point, once, centred | 1.6e-5 h at N = 512 | R4 |
| once, one-sided | order 1; moves d by 3.3e-3 h at N = 512 | R1, R4 |
| every step, centred | band mean q 37.6 after 13 events (N = 64); dist3_2 3e4 to 6e6 h at tau = 10 | R5; `rescale_survival.csv` |
| every step, one-sided first order | dist3_2 1.13 h on every rung, crossings moved 0.24 to 0.41 h | R5 |
| Gauss against leastSquares Hessian | ratio 0.89 to 1.02 at tau = 10, both diverged | R6 |

## Why it failed, or why we think so

DERIVED (report section 5, item 2): linearised about d, the first-order map with the centred
gradient is one step of the `linearQ` source with mu dt = 1, so it doubles the mode
$(-1)^{i+j}d$ per application. The discretisation error of one event (2.6e-4 h at N = 128) seeds
it. The one-sided gradient does not vanish at a grid-scale extremum, so it does not amplify the
mode, but it is a first-order estimate, and its fixed point is not d.

## Decisions

- Every-step rescaling with these gradients is a negative route on the UFVM line
  ([[negative-log]] 2026-10-31).
- The once-map stays available (inert by default) for initialisation.

## Open questions

- An event every k steps (k > 1) with transport between events.
- The polyhedral prediction: the FVM Hessian is inconsistent there (report section 5, item 4).
- The semi-Lagrangian version with the quadratic fit (T6), whose fit is not blind to the odd-even mode.

## Related

[[concepts/static-source-gate-eulerian]], [[concepts/projection-blend-extension]],
[[hubs/gradient-control]], [[sessions/sl-session-handover]], [[retractions/centred-q-eigenmode-growth-rate]]

## Log

### 2026-10-29
Created with the models and the pre-registered read-out R1 to R6.

### 2026-10-31
Settled on the Lichtenberg ladder (STATUS 11.28 items 7 and 8). R1 for r^2 - R^2 (second order)
reads FAIL as written because the order of a round-off error is meaningless and the read-out had
no clause for it; the value is round-off. The report's earlier phrase "a pointwise map with no
growth mechanism" was wrong and is corrected there.
