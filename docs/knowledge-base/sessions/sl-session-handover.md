---
title: "Handover to the semi-Lagrangian session"
description: "What the separate semi-Lagrangian session must know from the gradient-control session: the two parallel fixes, the metric fixes, the Eulerian port, the translating-droplet box study, the oscillating drift at N = 200, the gate infrastructure, the scoring corrections, the record corrections, and the open decisions"
kind: session
status: open
part: advection
tags: [session, part/advection]
date: 2026-09-28
---
# Handover to the semi-Lagrangian session

> Written 2026-09-28 by the gradient-control session. The four SL-baseline findings below moved
> from the gcls pre-print into the SL article on this date (the subsections "Consistency under
> domain decomposition", "Mass-momentum-consistent flux (rhoLENT)", the late-instability paragraphs
> of the translating droplet, and the oscillating drift in the Limitations). STATUS 11.13 to 11.15
> holds the lab record.

## What the SL session must know

1. **Two parallel defects of the coupled SL solver, fixed 2026-09-27** ([[concepts/coupled-face-density-defect]]). The face density on a processor face was rank-local (28d13f0): $\rho_f$ up to 90 % apart across a seam, np 4 against serial 1e-5 to 5e-4 before the fix and 1e-8 to 1e-10 after. The droplet metrics used internal faces only (b1798c3): the shape error differed 3.42e-2 between serial and np 4 before, 1.2e-7 after. Every SL two-phase result on more than one rank before that date carries both; the author decides which studies to re-run ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3624-L3744)).
2. **The Eulerian two-phase solver now uses the SL mass flux** (c094bd8, [[concepts/eulerian-solver-mass-flux-port]]): it had frozen $\rho$, $\rho\phi$ and $\mu_f$ and moved the heavy droplet at 42 % of the stream. The mass-flux code is shared through three headers; the SL files reconstruct byte for byte.
3. **The translating droplet at long times** ([[cases/translating-droplet]]): an outlet-triggered fast growth (10 mm box, onset at 0.063 s with the leading edge 35 cells from the outlet), and a slower interior growth (20 mm box: 30 1/s from 0.10 s; 40 mm box: completes 0.3 s, degradation 2.3x smaller in the current, 5.4x in the volume and 2.1x in the lead at $N = 142$ than at 100). Two rungs only; the third belongs on the cluster. The curvature of the moving interface does not converge (2.4 to 3.6 % at every rung). The gate's horizon is 0.05 s with `cellCentreInverse`.
4. **The oscillating droplet's gradient drift is unstable at $N = 200$** within ten periods ([[cases/oscillating-droplet]]): the band gradient error grows exponentially after 0.05 s and reaches 5.1 at T; the arm's finest rung cannot carry a verdict. The oscillating surface is now a token (`DROPLET_SURFACE`, default still the algebraic `implicitEllipsoid`, so a new oscillating study that does not set it runs the algebraic psi); whether the earlier algebraic-psi studies are void is an OPEN author decision ([STATUS 11.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3289-L3297)). CORRECTED 2026-09-29: this note first said they are void.
5. **The gate infrastructure** ([[concepts/method-gates]]): `config/gates/methodGate2D.yaml` and `3D`, the seam arms including the coupled `translatingSeamNp1` (baseline 3.5e-7 at tolerance 1e-5), NOT_COMPARABLE for a diverged reference, and two scoring corrections: `rhoClipFraction` is reported, not scored (it counts round-off clips), and the oscillating arm scores period and damping, not the velocity norms.
6. **The record corrections of METHOD.md** (4.1, 4.3, 6, 8.1, 10 carry CORRECTED notes; [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3545-L3623)); the 2D advection orders of 8.3.7 were 3/2 too high ([[retractions/advection-orders-3-2-factor]]); the curated `advConv2D*` CSVs must be regenerated on Lichtenberg.
7. **The 3D translating case** had no `dropletReferenceVelocity`: the disturbance columns of `traceTranslating3D*` are void ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3624-L3744)).

## Record inconsistencies found on 2026-09-29 (for the SL session to correct at the source)

The note writers of 2026-09-29 read the record against itself and found these; none is corrected yet.

1. **The production curvature was never scored on the ellipse gate.** `cellCentreInverse` is production
   since c935883 (2026-09-01), but the varying-curvature ellipse gate never scored it
   ([METHOD 8.1 L393](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393)); its second order is measured on constant curvature only.
   The acceptance criterion of STATUS section 7 (`G h^2 <= 0.65` and order `>= 1.9` on the ellipse) is
   therefore not demonstrated for the current method ([[cases/curvature-static-gates]]).
2. **The 6R box certification predates the seam fixes.** The 2D domain-size control (6R against 10R, 4R
   fails) ran with `PSI_FILTER biharmonicBand` on 8 ranks on 2026-08-18, before f83a1ab, 28d13f0 and
   b1798c3; no filter-off re-run is recorded, and the 3D gate and ladders use this box
   ([[cases/stationary-droplet]], [[retractions/psi-filter-seam-bug]]). Open, not void.
3. **The SL article's circle text disagrees with its own table:** the text says `O(h^1.2)` and 35 % to
   1 %; the table, refreshed on 2026-08-06 after the plane-fit guard fix, gives order 1.07 and 1.7 %.
4. **METHOD 4.2, the N = 512 row:** 0.890 % uncorrected against 1.51 % in the curated CSV (and 1.46 % from
   `3/4 h/R`); the rows N = 32 to 256 agree. HYPOTHESIS: the row and the article text predate the guard fix
   (the roadmap records that an N = 512 circle row read low before it).
5. **A stale jump probe:** the SL article's pressure-jump error "4.2 % to 1.5 %" is the old alpha = 1/2
   probe that METHOD section 7 calls an artefact; the comment of `createDropletMetricsFile.H` lines 25-26
   still describes it, while the code uses the pure phases.
6. **Stale header comments:** the first lines of `cases/oscillatingDroplet2D.parameter` and lines 17-20
   of the 3D stationary `blockMeshDict.template`.

## What is open for the author

`boundRho` ([[decisions/mass-flux-bound-rho]]); the re-run list; the void of the frozen-density Eulerian studies; a box-length token; the oscillating horizon; a baseline check of the gradient drift over time.

## Log

### 2026-09-28
Created.

### 2026-09-29
CORRECTED the void statement of item 4 (the void is OPEN); added the record inconsistencies found by the note writers.
