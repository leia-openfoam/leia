---
title: "The static source test: the Eulerian source terms alone, at zero velocity"
description: "Settled 2026-10-29 - every Eulerian source of leiaLevelSetFoam alone at zero velocity on a circle, eight starts, N = 64 to 512: with the centred gradient every law grows grid-scale extrema at its G(0), no law restores d without moving the zero set, and the DA volume depends on the profile at O(h^2); next: the algebraic rescaling of D. Bothe"
aliases: [static source test, staticSource, 2DstaticCircle, source terms at zero velocity]
kind: concept
status: settled
outcome: negative
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-10-14
code: [cases/2DstaticCircle, cases/2DstaticCircle.parameter, config/staticSource, workflow/scripts/static_source_metrics.py, workflow/scripts/static_source_smoke.sh]
sources: [STATUS 11.27, the plan of 2026-10-14, the config headers of config/staticSource]
---
# The static source test: the Eulerian source terms alone, at zero velocity

> Settled 2026-10-29 (STATUS 11.27 part 3; report `docs/gradient-controlled-level-set/
> gcls-static-source-report/`). MEASURED, 608 runs, N = 64 to 512: with the centred gradient every
> law grows grid-scale extrema of psi at its G(0) (linearQ 0.98 to 1.00, linearZ 0.46 to 0.50,
> saturatedLinearZ 0.331 to 0.337, twoThirdsZReg 0.59 to 0.62, softWall 1.00 to 1.02); no law
> restores d from a perturbed profile without moving the zero set; the cubic laws and softWall keep
> the zero set of d only because they hardly act near q = 1; on the whole domain every law fails by
> tau = 10 to 50 from the kink of d. From d in the band the zero crossings move 0.46 to 0.032 h at
> tau = 80, order 2.35 (the slow drift converges). The Detrixhe-Aslam volume of
> (r^2 - R^2)/(2R) differs from that of d by O(h^2) (order 1.97). The first run was void
> ([[retractions/static-source-solver-dead-zone]]). Next: the algebraic rescaling of D. Bothe
> (STATUS 11.28).

## What it is

The author's request of 2026-10-13/14, before the work returns to advection: measure each
level-set source term on its own. At zero velocity the level-set equation of `leiaLevelSetFoam`
reduces to $\partial_\tau\psi = S(\psi)$, so a run integrates the source alone in pseudo-time
$\tau = t$ (every rate is 1 1/s).

- Case `cases/2DstaticCircle`: the unit square, $N \times N \times 1$ cells, $N$ = 64, 128, 256,
  512; a circle of radius $R = 0.25$ with the centre (0.5013, 0.4987) off the cell lattice;
  $\Delta\tau = 0.5h$; $\tau$ to 80.
- Eight initial fields, all with the exact circle as zero set (OpenFOAM `setExprFields`):
  the signed distance $d$; the blends $d + \varepsilon d^2/(2R)$ with $\varepsilon = (h/R)^2$,
  $h/R$, 0.1; $(r^2 - R^2)/(2R)$; $r^2 - R^2$ ($q = 2R$ on the circle); $d(1 + 0.1\cos m\theta)$,
  $m$ = 2, 4.
- Sources: the gradient-control laws `linearQ` (Soga's restoring term, the Bothe-Fricke-Soga
  bounded-gradient equation), `linearZ`, `cubicQ`, `cubicZ`, `twoThirdsZReg`, `saturatedLinearZ`
  (the Bothe-Marić-Soga dossier laws), `softWall` (Bothe S1a, `SW_EPS_D 0.8`); `beta`; and the
  strain sources `R`, `Rdiv`, `RdivStrictSp` (Fricke, Marić, Vučković, Roisman, Bothe), which are
  zero at $u = 0$ and serve as null checks. Source regions: the whole domain (`MOLLIFIER none`)
  and a 3-layer band.
- At the zero set the source obeys the ODE $dq/d\tau = q\,G(q)$ exactly, so each law has a
  closed-form prediction there.

## Why it matters

It separates what a source does to the profile from what the transport does, and it measures
whether the Detrixhe-Aslam volume changes when only the profile changes: the DA phase indicator
depends on $\psi/|\nabla\psi|$, not only on the zero set.

## Where in the code

`cases/2DstaticCircle`, the eight configs of `config/staticSource/` (stage A: Control, NullR,
NullStrain; stage B: GcNone, GcBand, Beta; stage C: the $N = 512$ rungs), the metrics script
`workflow/scripts/static_source_metrics.py` (the exact cell areas of the circle, the zero
crossings against the exact circle, L1 and L2 norms only, a verdict row for every pre-registered
line), and the laptop smoke `workflow/scripts/static_source_smoke.sh`.

## Evidence

| claim | number | where |
|---|---|---|
| The psi solver must iterate in this test. | First run: 68 % of the psi solves with a nonzero source did no iteration; 39 cases never changed psi. | MEASURED, [[retractions/static-source-solver-dead-zone]] |
| `minIter 1` fixes it. | The case that froze iterates in all 1536 solves and moves; S0 PASS in 301 of 301 smoke cases. | MEASURED, STATUS 11.27 part 2 |
| `minIter 1` is not bit-inert where the solver already iterated. | 149 of 190 such smoke cases differ, at round-off (at most 4.3e-11 relative in psi), except softWall on the whole domain (0.27 to 0.59). | MEASURED, STATUS 11.27 part 2 (G3) |
| softWall on the whole domain amplifies round-off. | Serial against np 4: 7.5e-4 in the volume columns, 20 % in the gradient columns at $N = 64$, $\tau = 2$; the four other non-cubic laws agree to round-off. | MEASURED, STATUS 11.27 part 2 |

## Open questions

The pre-registered read-out (P0 to P9 in the config headers): the representation orders of the
eight fields, the zero-set ODE, the DA volume of the normalised algebraic field against $d$, the
centred-$q$ damage of E1 transferred to the Eulerian line, the whole-domain arm.

## Related

[[hubs/gradient-control]] - [[concepts/sdpls-source-eulerian]] - [[models/gradient-control-law]] - [[models/sdpls-source]] - [[concepts/gradient-control-next-experiments]] - [[concepts/source-discretisation-defects]] - [[models/phase-indicator]]

## Log

### 2026-10-14
Created. The first cluster run (010b553a) voided; the fix edb1d5bb gated on the laptop; the
re-run at 3727b9e9 submitted on Lichtenberg: build 55441730 (`CI: PASS`, stamp
`shared-method-config-2026-09-01-338-g3727b9e9`), stage A and B drivers 55441743 to 55441748
(`leia-s2-*`); stage C after stage A passes.

### 2026-10-29
SETTLED with the final metrics (job 55753419, 608 runs): growth at G(0) for every law, no recovery
of d without a moved zero set, the DA volume depends on the profile; P4.2 fails as written (the
shift of d converges, 0.032 h at N = 512). Negative-log lines added. Next: the rescaling (STATUS 11.28).
