---
title: "The first static source test: the psi solver did nothing below its tolerance (2026-10-14)"
description: "VOIDED 2026-10-14 - every result of the static source studies at commit 010b553a (stages A, B and the cancelled C on Lichtenberg, and the laptop smoke): at u = 0 the psi matrix is diagonal, the residual is normalised by the whole field, and a small source update started below tolerance 1e-10, so PBiCGStab returned psi unchanged; 68 % of the psi solves with a nonzero source did no iteration"
aliases: [psi solver dead zone, static source void, staticSource_VOID_psiSolverDeadZone_2026-10-14, minIter 1]
kind: retraction
status: voided
part: gradient-control
tags: [retraction, part/gradient-control]
date: 2026-10-14
code: [cases/2DstaticCircle/system/fvSolution.template, cases/2DstaticCircle.parameter, workflow/scripts/static_source_metrics.py, workflow/scripts/static_source_smoke.sh]
sources: [STATUS 11.27 part 2, "CLAUDE.md, Solver convergence in a cancellation-dominated problem", "CLAUDE.md, A wrong setup voids its data", commit edb1d5bb]
---
# The first static source test: the psi solver did nothing below its tolerance (2026-10-14)

> VOIDED 2026-10-14. The static source test (STATUS 11.27) integrates `ddt(psi) == S(psi)` at zero
> velocity in `leiaLevelSetFoam`, so that each Eulerian source term is measured alone. The first
> run, at commit 010b553a, did not integrate the source: wherever the source update was small, the
> linear solver of the psi equation returned psi unchanged. Every number of that run is void,
> including the null checks that look unaffected. Fixed in commit edb1d5bb (`minIter 1` for the
> psi solver of `cases/2DstaticCircle`, token `PSI_MIN_ITER`).

## The claim, and where it lived

No result of the run was published. The void covers what the run produced and what was read
from it before the defect was found:

1. The Lichtenberg studies `studies/staticSource/*` of the worktree `leia-static` (stage A:
   Control, NullR, NullStrain; stage B: GcNone, GcBand, Beta; 496 cases, all COMPLETED), now
   `studies/staticSource_VOID_psiSolverDeadZone_2026-10-14` with a README and the iteration counts.
2. Stage C (the N = 512 rungs of GcNone and GcBand), cancelled at 20:28 with its 112 running jobs
   and both drivers (55440959, 55440960), by id.
3. The stage A null check (job 55440912, 114 bit-identical pairs) and the first metrics job
   (55441345, metrics commit 4025143f).
4. The laptop smoke of STATUS 11.27 part 1. Its reading "the whole-domain sources amplify
   round-off far from the interface; serial and np 4 differ by 17 % in the gradient columns at
   tau = 2, N = 64" was re-examined with the fix and CORRECTED (see the log below).

## Why it was wrong

| claim | number | where |
|---|---|---|
| The run integrates the source. | linearQ in the band on the signed distance at N = 256: zero iterations in all 122880 psi solves; psi at tau = 80 bit-identical to psi0. | MEASURED, STATUS 11.27 part 2 |
| The defect is rare. | 415 of the 416 cases with a nonzero source (gradient control and beta) have psi solves without an iteration; 19.8 of 29.0 million such solves did no iteration (68 %; median per case 0.70). | MEASURED, the solver-performance lines of every log |
| The fine rungs show the converged behaviour. | 39 cases never changed psi: every cubicQ, cubicZ and softWall band case on d and on the h^2 blend at every rung (G'(1) = 0 makes F tiny near q = 1), and linearQ, linearZ, saturatedLinearZ and beta in the band on d and on the h^2 blend at N = 256. The band arm of linearQ on d gave a zero-crossing shift of 0.44 h, 0.30 h and exactly 0 at N = 64, 128 and 256: a false convergence. | MEASURED |
| The smoke checked the setup. | 60 of its 176 np-4 cases already had zero-iteration solves; the smoke checked COMPLETED, check-init, null and seam, not whether the solver did anything. | MEASURED |

The mechanism. At u = 0 the psi matrix holds only `ddt` and the implicit source `Sp(F, psi)`,
so it is diagonal. OpenFOAM normalises the residual by sum(|A psi - A psibar| + |b - A psibar|),
here of the size of the whole field over dt. A source that is nonzero only in a band, or of order
h^2 on a signed distance, starts below the tolerance 1e-10, and PBiCGStab stops before its first
iteration. This is the cancellation rule of `CLAUDE.md` seen from the other side: below the
tolerance the solver does not under-converge, it does nothing.

## What survives

Nothing of the run is used. What the run taught is kept: the mechanism, the iteration counts as
the measurement of the defect's size, and the check that now runs with every analysis.

## The fix

1. `cases/2DstaticCircle/system/fvSolution.template`: `minIter @!PSI_MIN_ITER!@;` in the psi
   solver; `PSI_MIN_ITER 1` in `cases/2DstaticCircle.parameter`. One DILU-PBiCGStab iteration
   solves a diagonal system exactly; with F = 0 the residual is exactly zero and PBiCGStab stops
   before any update, so the null arms keep psi.
2. `static_source_metrics.py iterations` counts, per case, the psi solves and those without an
   iteration (the solver-performance lines only; not a state classifier). `analyse` adds the row
   S0: with F != 0 no zero-iteration solve, with F = 0 no iteration.
3. `static_source_smoke.sh` runs the case that froze (linearQ, band, d, N = 256, tau = 1, np 4).

## Propagation (checklist)

- [x] `STATUS.md` 11.27 part 2.
- [x] [[retraction-log]], [[decisions/process-operating-rules]] (a new row), [[decision-log]].
- [x] Lichtenberg: the studies renamed `_VOID_psiSolverDeadZone_2026-10-14`, with a README.
- [x] The session memory of the laptop (a pointer only).

## Related

[[concepts/pressure-projection-and-linear-solvers]] - [[concepts/bit-identity-and-inertness-gates]] - [[concepts/sdpls-source-eulerian]] - [[hubs/gradient-control]] - [[concepts/wrong-setup-voids]]

## Log

### 2026-10-14
VOIDED: found in the first metrics of the run; the mechanism confirmed from the solver-performance
lines of all 608 logs; stage C cancelled; the fix edb1d5bb with the gate G1 to G4 pre-registered
in its commit message.
GATE (laptop, 301 cases): G1, G2 PASS (the dead-zone case iterates and moves; S0, check-init and
null all PASS); G3 FAIL as pre-registered: the fix is not bit-inert where the solver already
iterated, because PBiCGStab also tests minIter at its early exit after the first half-step
(`PBiCGStab.C` line 223, v2606); the change is round-off (at most 4.3e-11 relative in psi at
tau = 2) except for softWall on the whole domain (0.27 to 0.59). G4 as in part 1.
CORRECTED: part 1's reading "the whole-domain sources amplify round-off" holds for softWall (serial
against np 4: 7.5e-4 in the volume columns, 20 % in the gradient columns) and less for cubicQ
(4.6e-4); linearQ, linearZ, saturatedLinearZ and twoThirdsZReg agree to round-off. The same per-law
pattern is in the void smoke, so the dead zone did not cause it.
