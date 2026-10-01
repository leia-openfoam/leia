---
title: "Handover to the semi-Lagrangian session"
description: "What the separate semi-Lagrangian session must do and know from the gradient-control session: three tasks (the frozen-viscosity table row, the clip-on polyhedral orders, fixed psi patch values), the results of 2026-09-29, the two parallel fixes, the metric fixes, the Eulerian port, the translating-droplet box study, the oscillating drift at N = 200, the gate infrastructure, the scoring corrections, the record corrections, and the open decisions"
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

## Tasks for the SL session (2026-09-29)

Three findings need a decision of the author in the SL line. Each task names the evidence, the
decision and the cheapest discriminator.

**T1. A row of the viscous-term table ran with a frozen face viscosity** (item 7 below).

- Evidence: the `mu_f interpolated` row of the viscous-term table (`sec:viscous`,
  [SL article L801-L829](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L801-L829),
  `docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`) comes from the
  run at 4267d7b. At that commit the viscous term read `muf` for every model, and `muf` was rebuilt
  for the geometric models only. The row entered the article at a2cfb2a, before the fix 39e59b3
  ([[concepts/viscosity-open-items]]). The claims that rest on it: the 36 % separation after 2667
  steps, and the gains 0.64x (L1), 0.81x (L2) and 3.3x (shape).
- Decision: re-run the arm on a post-fix binary and replace the row, or retract the row and its
  claims by the retraction rule.
- Cheapest discriminator: re-run the table's matrix (the `mu_f interpolated` arm and, because a
  setup wrong in one way is not assumed wrong in one way only, the geometric and harmonic arms of
  the same run) with the current binaries at the table's N and horizon, on the laptop at np 4, and
  compare with the published rows.

**T2. The published polyhedral transport orders ran with the clip on** (item 17 below).

- Evidence: the orders 3.28 (3D shear) and 1.46 (3D deformation) come from configs with
  `SL_CLIP true` and the face stencil
  ([`uncachedConv3DshearPoly.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/uncachedConv3DshearPoly.yaml),
  [[concepts/value-bounds-and-clips]]). The production default is `SL_CLIP false`
  ([METHOD 8.1 row SL_CLIP](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L410)), and with it the coarsest polyhedral 3D
  shear rung diverges at step 198 ([METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L608-L668), finding 3).
  The article's "Boundedness" subsection
  ([SL article L543-L549](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L543-L549))
  still calls the clip order-preserving.
- Decision: state the clip setting next to the orders and correct "Boundedness", or re-run the
  polyhedral ladders with the production default and report the divergence as the result.
- Cheapest discriminator: the first option needs no run (a text change); for the second, the
  coarsest rung with `SL_CLIP false` already has a recorded answer (divergence at step 198).

**T3 (new). A fixed psi value on a patch makes the SL update unstable.**

- Evidence ([STATUS 11.19 item 3](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4201-L4226), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/CLAUDE.md#L528-L532)):
  the SL fit reads every physical patch value as a stencil datum (`SL_STENCIL_BOUNDARY_FACES
  include`, the default). On the one-way `2Dtranslation`, the exact psi on the outflow patch made
  the outflow-edge error grow by a factor of about 1.1 per step (4.5e5 at T at N = 256), and the
  exact psi on the inflow patch failed at CFL 1, where the departure point of the first cell column
  lies outside the domain (`E_GEOM_ALPHA_REL` 1.79 at N = 128). zeroGradient everywhere is stable.
  No production case gives psi a fixed value.
- Decision: whether the SL scheme must accept fixed psi data (an inflow of the second phase, a jet,
  a filling), and then how it treats a departure point outside the domain.
- Cheapest discriminator: `workflow/scripts/translation_bc_probe.sh` (with `_scan.py`) with
  `stencilBoundaryFaces inflowOnly` on the `exactAll` variant at N = 128, CFL 0.5 (does the outflow
  instability go when the outflow faces leave the stencil?), and on the `exactInflowOnly` variant at
  CFL 1 (does the inflow need the boundary value at the departure point?). Seconds per run, serial.

**T4 (2026-10-01). The CFL 1.0 row of the 2D convergence table is retracted** ([[retractions/sl-stale-band-alpha-metrics]]).

- Evidence: `leiaSemiLagrangeLevelSetFoam` computed alpha with the narrow band of psi^n (fixed in
  6f63418a). On the published ladder `config/uncachedConv2Dvortex.yaml` the fixed solver gives the
  CFL 1.0 shape order 2.465 (published 2.378) and volume order 3.193 (published 3.542); the shape
  error at T is 10.5 to 23.3 % lower at five of seven rungs. The CFL 0.5 row stands (orders within
  0.005). psi and every gradient metric are unchanged ([STATUS 11.23](https://github.com/leia-openfoam/leia/blob/d2984c5e/STATUS.md#L4947-L5055)).
- Decision (author): regenerate the curated tables of this solver. 19 tables carry its rows; 8 have
  arms at CFL 0.8 or 1.0 (`uncachedConv2Dvortex`, `npslConv2Dvortex`, `nslConv2Dvortex`,
  `sdCompare2D`, `linearConv2Dvortex`, `linearConv2DvortexClip`, `linearConv3Dshear`,
  `kinematicTranslation2D`). Then the generated `convergence_orders*.tex`, the deck and the prose.
- Cheapest discriminator: the 2D studies with the current binaries on the laptop (the published 2D
  ladder runs in about one minute per arm pair at np 4); the 3D studies on Lichtenberg.

## Results of 2026-09-29 (information for the SL session)

1. **`2Dtranslation` translates one way, and every earlier number of the case is void**
   ([[retractions/reversed-2dtranslation]], [STATUS 11.19 items 1 and 2](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4176-L4200)).
   The author decided on 2026-09-29 that the case translates left to right (item 16 below). It now
   has `OSCILLATION off`, the exact end reference `psiEnd`/`alphaEnd` from
   `workflow/scripts/write_end_reference.py` (token `END_REFERENCE`, inert default `none`) and psi
   zeroGradient on all four patches. The three kinematic solvers print which error reference they use.
2. **The new baseline of the hex 2D regression rung** is `studies/advConv2Dtranslation` of 2026-09-29
   (laptop, np 8; [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L685-L704)): `none` 5.364e-02 / 7.210e-03 /
   1.609e-03 / 4.442e-04 at N = 32 to 256, orders 2.90, 2.16, 1.86. The "saturation at N = 256" is
   retracted. The curated table is `advConv2Dtranslation_convergence.csv` (method-comparison data).
   The finalize rule also writes `advConv2Dtranslation_errors.csv` into the SL article's and deck's
   `data/tables/`; those copies are not committed, and the SL session decides whether they belong there.
3. **The cone bound on the one-way translation** is 14.1x (unity) and 9.5x (stencil) worse than `none`
   at N = 64 and does not converge over N = 32 to 256 ([METHOD 8.3.4](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L608-L668)).
   The falsification stands and is stronger.
4. **`kinematicTranslation2D`** re-ran ([STATUS 11.19 item 5](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4234-L4268)):
   second order at fixed CFL 0.25 and 0.5 on the geometric alpha error; the prediction "the error
   collapses at CFL 1" is falsified again; the null control passes (5.0e-12 column-scaled).
5. **The production curvature `cellCentreInverse` is scored on the signed-distance ellipse and
   ellipsoid** (item 1 below; [METHOD 4.1](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L169-L193),
   [STATUS 11.19 item 6](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4269-L4302)): second order on both (ellipse 1.98 and
   2.00 over N = 128 to 512; ellipsoid 2.10 over N = 50 to 128, 1.00 without the Gaussian term), 14 to
   26 % below the per-face inverse. Its gain equals the per-face inverse's: G h^2 0.647 at N = 512, so
   the STATUS 7 criterion as the curvature plan applied it (the finest rung) is met. The stricter
   per-rung form fails at N = 128 (0.651) and 256 (0.673); at N = 256 every second-order delivery fails
   it, the per-face inverse included; only the first-order cell-mean and symmetric face-mean deliveries
   stay below 0.65.
   On the implicit psi every delivery is first order.
6. **Traps** ([STATUS 11.19 item 7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/STATUS.md#L4303-L4312)): `leiaSetFields` is not
   idempotent on a non-pristine `0/`; a verification study through the full workflow runs the finalize
   rule, which overwrites curated outputs (use `--until solve`); the per-value relative test of
   `compare_metrics_csv.py` fails on round-off columns (read the column-scaled difference).

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

1. **RESOLVED 2026-09-29 (Results 5): the production curvature was never scored on the ellipse gate.** `cellCentreInverse` is production
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
7. **Task T1. A table of the SL article carries the frozen face-viscosity bug.** The `mu_f interpolated` row of
   the viscous-term table (section `sec:viscous`, the table near L805-L829 at d1e3414) comes from the
   4267d7b run, when `muf` was rebuilt for the geometric models only; the row entered the article at
   a2cfb2a, before the fix 39e59b3. The claims built on it (36 % after 2667 steps; 0.64x, 0.81x, 3.3x)
   need a re-run or a retraction: an author decision ([[concepts/viscosity-open-items]]).
8. **The same bug is still reachable.** `updateFaceDensity.H` line 10 skips `faceViscosity.H` when the
   mass flux is `interpolatedDensity` or `SL_FREEZE_RHOPHI` is set, and `interpolatedDensity` is the code
   default of `createMassFluxFields.H`; no config selects it today.
9. **The "58x" attribution.** 58x and 0.073 s belong to the arm with Kang AND the sharp Heaviside
   together (4.5e-3 / 7.7e-5 = 58.4); `docs/plan-curvature-stabilization.md` and the SL deck give it to
   Kang alone, thirteen case templates to `sharpHeaviside` alone, and `sharpHeaviside` with the
   arithmetic face value cannot run ([[concepts/kang-gfm-and-sharp-heaviside]]).
10. **Stale claims in comments and headers:** thirteen templates call CST "MEASURED SUPERIOR on the static
    droplet"; `UEqn.H` lines 25-35 still name "interpolated" as the default and keep a retracted "2.08x
    worse"; the `faceViscosity.H` header names keys the solver no longer reads; the header of
    `config/variationalForce2D.yaml` calls interFoam's curvature an exact variation (plan-shannon lines
    616-619 corrects it); `cases/translatingDroplet2D/system/fvSolution.template` lines 405-407 keep a
    closed-box "MEASURED 2026-09-01 ... outlet at t = 0.08" comment without a VOID marker; nine
    `.parameter` files start with the stationary-droplet header.
11. **A closed-box trap in a committed mesh file:** `cases/transISTDroplet2D/system/blockMeshDict` puts all
    four sides in one `walls` patch; the `.template` is right and the workflow overwrites the file, but a
    run from the unrendered case gets a closed box.
12. **Two numbers that disagree inside STATUS 11.15:** the 40 mm box at N = 100 and t = 0.3 s has a volume
    change of 1e-1 in the first table and its sentence, and 1.5e-1 in the second table (the archived
    history: 0.147; 15 % is right); the Popinet L2 order is 0.88 at one place and 0.91 at another (0.88
    comes from the rounded table values).
13. **A wrong cross-reference in the SL article:** lines 438-441 point to `sec:popinet` for the 3D
    polyhedral divergence at step 3, but that section covers the 2D benchmark only.
14. **Case bookkeeping:** `3Dcontactline` is a 2D mesh (one cell thick); four cases have no study config
    (`3Dtranslation`, `2Dcontactline-periodic`, `2Dcontactline-vortex`, `3Dcontactline`)
    ([[cases/benchmark-cases]]).
15. **CLAUDE.md said a second clone "ran" a library about 200 commits ahead** (2026-09-09); STATUS 9.5
    records that no job ran with it. Corrected in CLAUDE.md and AGENTS.md on 2026-09-29; the comment in
    `src/leiaLevelSet/leiaVersionRegistry.H` line 11 still says "ran" (a code comment; left for the next
    library change, since an edit changes the library stamp).
16. **RESOLVED 2026-09-29 (Results 1 to 4; the author decided one-way): `2Dtranslation` is a reversed translation.** Its `velocityModel` block sets no `oscillation`
    entry, the default is on ([velocityModel.C L49-L50](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/velocityModel/velocityModel.C#L49-L50)),
    and `tau` defaults to `endTime` = 0.5 s, so U = U0 cos(pi t/tau): the circle moves at most
    tau/pi = 0.16 and returns at T, while the case comment describes a one-way translation with the
    exact solution psi0(x - U t). The ladder `advConv2Dtranslation` (METHOD 8.3.7) read its errors at T,
    where the reversal cancels errors. DERIVED from the code; a solver log confirms it. OPEN, author
    decision: a wrong setup (re-run with `oscillation off`) or a reversed-flow gate read at the right
    instants ([[cases/kinematic-advection-cases]], [[hubs/advection]]).
17. **Task T2. The published polyhedral orders ran with the clip on:** 3.28 (3D shear) and 1.46 (3D deformation)
    come from configs with `SL_CLIP true` and the face stencil; with the production default the
    coarsest polyhedral 3D shear rung diverges at step 198; the SL article's "Boundedness" subsection
    (L543-L549) still calls the clip order-preserving ([[concepts/value-bounds-and-clips]]).
18. **Small items:** `workflow/README.md` line 354 still gives 128 x 128 x 256 cells for the 3D shear (the
    case is N^3 now); `cases/3Dtranslation_poly.parameter` sets END_TIME 3.0, which moves the sphere
    out of the box at t = 1.65 (no study runs `3Dtranslation`); the GRL article's text and table give
    different one-step volume changes (1.05e-4 against 1.11e-5 at the coarsest rung).

## What is open for the author

`boundRho` ([[decisions/mass-flux-bound-rho]]); the re-run list; the void of the frozen-density Eulerian studies; a box-length token; the oscillating horizon; a baseline check of the gradient drift over time.

## Log

### 2026-09-28
Created.

### 2026-09-29
CORRECTED the void statement of item 4 (the void is OPEN); added the record inconsistencies found by the note writers.
Added items 7 to 15 from the reports of the note writers.
Added items 16 to 18 (the reversed 2Dtranslation, the clip-on polyhedral orders, small items).
Added the tasks T1 to T3 and the results of 2026-09-29; items 1 and 16 are resolved.

### 2026-10-01
Added task T4 (the stale narrow band; the CFL 1.0 row of the 2D convergence table retracted).
