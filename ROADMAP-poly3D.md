# ROADMAP — 3D polyhedral verification: stationary droplet and translating droplet

Hand-off for another session. Written 2026-09-09 at commit `f4324c2`, branch `development`.
Read [CLAUDE.md](CLAUDE.md) first (rules), then [CLUSTER.md](CLUSTER.md) (cluster), then
[STATUS.md](STATUS.md) (full state). This file is the polyhedral 3D plan only.

All communication uses ASD-STE100 Simplified Technical English.

---

## 1. The one open defect

The semi-Lagrangian quadratic fit amplifies the level set in cfMesh's small one-sided cells.
The far field grows a false zero set, and the run then fails. Three measurements fix the
mechanism:

1. **The static amplification bound** (`fitProbeDisplacement`, library diagnostic, writes
   `slFitAmplification` into `0/`). At `d = -u dt` the bound is `Lambda_max` 1.0527 on hex and
   1.2608 on poly. The whole polyhedral excess sits in cells of 0.32-0.33 h. `Lambda` falls
   monotonically with cell size.
2. **`Lambda - 1` is proportional to the displacement**, hence to dt. The accumulated bound
   over a fixed physical time is `exp(c U T)` and does not contain dt.
3. **The dt sweep confirms 2 experimentally** (job 54503136, `config/popinet3D_poly_sigma0_
   dtSweep.yaml`). The sigma = 0 control failed at the steps 528 / 1009 / 1974 for dt, dt/2,
   dt/4, at the times 4.877e-03 / 4.660e-03 / 4.558e-03 s. The growth rate per unit physical
   time is constant to 5 %. The pre-failure trajectories agree to 0.3 % at matched times.

**Conclusion: a smaller time step cannot repair this. The fix must change the update operator.**

The Courant number is NOT the cause and is evaluated correctly: hex 0.0164 everywhere, poly
bulk 0.0230, and the small far cells reach 0.1507. That is 6.5 times the bulk value and far
below any classical limit. The step comes from the capillary law, not from a Courant limit
(`adjustTimeStep no`).

## 2. The fix — built, and requirement 3 is RETRACTED

**RETRACTION, 2026-09-09.** Requirement 3 below said the clip must act in **the far field
only**, because a global clip costs +30 % volume error at the interface. The +30 % is real and
reproduces exactly. **Its mechanism was wrong.** The damage does not come from clipping the
narrow band, and withholding the clip from the band removes NONE of it: measured on Popinet's
2D hexahedral translating droplet at N = 64, `clipRegion outsideBand` moved every interface
metric by exactly what the global clip moved it — volume +30.4 %, centroid +12.5 %, shape
+4.9 %, L1|u'| +8.9 %, to every printed digit — while correctly withholding the clip from 144
band cells.

`config/popinet2D_clipFiringProbe.yaml` located the firings. That mesh has exactly **six**
cells where psi is the extremum of its own stencil — the four box corners (local maxima of the
distance field) and the two cells at the droplet centre (the apex of the distance cone, a
local minimum) — and **every** cell the clip fired in was one of them. Nothing else fired, at
any step.

A quasi-monotone bound cannot represent an extremum. Where psi_c is already the stencil
minimum, `lo` IS psi_c, so every reconstructed value below it is pulled back up and the
extremum is flattened at every step. The fit undershoots at the apex because a smooth
quadratic cannot follow a non-differentiable minimum. So the clip fires on a **uniform
hexahedral mesh**, which requirement 1 assumed it would not, and the band is irrelevant to it.

**The second rule, and it is classical**: exempt a cell that is ALREADY its stencil's
extremum — the standard exemption of a monotonicity-preserving limiter at a smooth extremum.
It carries no coefficient: the test is `psi_c == lo` or `psi_c == hi`, exact in floating point
because lo and hi are taken over a stencil that contains psi_c.

### What is implemented (all runtime-selectable, all inert by default)

| entry (`fvSolution.levelSet.semiLagrangian`) | token | default | what it does |
|---|---|---|---|
| `clipToStencilBounds` | `SL_CLIP` | false | the bound itself (pre-existing) |
| `clipRegion` | `SL_CLIP_REGION` | `all` | `outsideBand` withholds the clip from the sign-change narrow band |
| `clipKeepExtrema` | `SL_CLIP_KEEP_EXTREMA` | false | exempt a cell that is already its stencil's extremum |
| `valueBound` | `SL_VALUE_BOUND` | `fromClipSwitch` | selects the bound: `none` \| `stencilBounds` \| `lipschitzCone`. The sentinel follows `clipToStencilBounds`, so every pre-2026-09-10 case is unchanged |
| `lipschitzMode` | `SL_CONE_L_MODE` | `unity` | where L comes from: `unity` (L = 1, coefficient-free) \| `stencil` (measured, a diagnostic) |
| `lipschitzConstant` | `SL_CONE_L` | 1 | the eikonal value; anything else is a tuned coefficient |
| `onInadmissible` | `SL_CONE_INADMISSIBLE` | `cellOnly` | recovery when the cone interval is EMPTY: owner-only interval \| no bound |

Code: `slCorrector::robustEvaluate` (the bound and the exemption),
`slCorrector::buildClipMask` (the region), `slReconstruction` (the entries),
`slCorrector::reportClipActivity` (the counters and `slClipEligible` /
`slClipFired` / `slClipFiredEver`).

### Measured, 2026-09-09, Popinet 2D hex N = 64, full horizon 1563 steps

| | L1\|u'\| | L2\|u'\| | volume err | shape L2/R | centroid/R | kappa L2 band | clip activity |
|---|---|---|---|---|---|---|---|
| baseline, clip off | 9.0723e-05 | 2.0419e-04 | 1.2123e-03 | 4.0755e-03 | 4.1908e-03 | 1.7498e+02 | — |
| clip, no exemption | +8.9 % | +3.3 % | **+30.4 %** | +4.9 % | +12.5 % | -3.1 % | 44/step |
| clip + exemption | +0.3 % | +0.2 % | **+1.03 %** | +0.2 % | +0.4 % | +0.1 % | **3 cell-steps, total** |

The exemption cuts the damage by a factor 30 and the clip's activity from over 600 cell-steps
to THREE in a 1563-step run. The four clip-off arms are bit-identical to each other, so both
new tokens are inert (G0 and G1 also pass: the rendered case differs by the token lines only,
and the metric CSV is bit-identical).

### The residual, characterised exactly

The candidate is NOT yet bit-identical on hex, so requirement 1 is not met. The residual is
**one cell**: 4000, the droplet-centre apex. The apex minimum is shared by a symmetric PAIR of
cells (4000 and 4128, equal to the last digit at t = 0). The exact equality test holds while
that tie is exact and fails the moment a rounding-level perturbation breaks it — the survivor
sits at `(psi_c - lo)/(hi - lo) = 0.004`, a hair off the minimum, so it is clipped, which
pushes it further off. Three cell-steps of clipping in the first steps then amplify to +1.03 %
volume error by step 1562: identical to eleven digits through step 100, first byte difference
at step 20, then growth.

G4 has now settled what to do about this residual: nothing, until the conflict in section 2a
is resolved. The apex detector proposed here is retracted as stated — see 2a.

### One diagnosis that was tested and did NOT hold

The full-horizon probe found the pre-exemption survivors in a strip at `i = 126`, the single
cell layer just inside the OUTLET, all near-minima at `pos ≈ 0.006`. `slReconstruction.H`
records that `stencilBoundaryFaces include` makes the outlet carry a stationary alternating
error, and `inflowOnly` is the gated remedy. It works on the firings and **not** on the damage:
`config/popinet2D_clipStencilGate.yaml` removed those eighteen firings and left the volume
error at **+1.027 %**, against +1.0 % with `include`. The `inflowOnly` baseline's own final
metrics match the `include` baseline to nine digits. So the outlet layer is not the residual.

Requirements 1, 2, 4 and 5 stand unchanged:

1. **It must not act on hexahedral meshes at all.** Hexahedral meshes run correctly and must
   stay bit-identical. NOT YET MET — see "The residual" above.
2. **No mesh-type switch.** MET: both triggers are local algebraic tests, and neither names a
   mesh type.
3. ~~**The far field only.**~~ RETRACTED, see above. The band rule is implemented and measured
   to remove zero firings on hex; whether it earns its place must be argued on the polyhedral
   mesh or the entry should be dropped rather than kept as a second knob.
4. **Runtime-selectable and inert by default.** MET and gated.
5. **Unstructured FVM only**: per cell, compact stencil, MPI-decomposable, one implementation
   for 2D and 3D. MET by construction; the 4-rank gate passed and the serial-versus-np4 check
   is `config/polyDroplet3Drefined_clipRegionGate.yaml`.

The theory behind the clip is unchanged: a reconstruction with `Lambda = 1` in a cell cannot
create a new extremum, because the update is then a convex combination of the stencil values.
The clip enforces that property where the fit does not have it — and the exemption keeps it
from destroying the extrema the level set genuinely has.

## 2a. G4 FAILED: the bound and the exemption are in direct conflict (2026-09-09)

`popinet3D_poly_sigma0_clipGate`, four arms, np 32, all COMPLETED at 650 steps, sigma = 0 so
only the transport is under test. Every arm ran the `curvature-v2512` binary.

| arm | clip / keepExtrema | failure step | zeroSet/R at T | volume err at T |
|---|---|---|---|---|
| 00000 | false / false | **527** | 3.4268e-01 | 1.0112e-04 |
| 00001 | false / true | **527** | 3.4268e-01 | 1.0112e-04 |
| 00002 | true / false | **NONE** | 8.9732e-04 | 3.9322e-06 |
| 00003 | true / true | **506** | 1.2678e-01 | 2.2880e-05 |

The controls hold: arm 00000 reproduces the recorded failure at step 528 to one step (so the
clip code is bit-inert on a polyhedral mesh too), arms 00000 and 00001 are bit-identical, and
arm 00002 re-establishes in SI that the global clip removes the defect.

**The candidate fails at essentially the same step as no clip at all** — 506 against 527 is
4 %, inside the 5-38 % scatter this campaign documents for genuine instabilities.

**Why, and it is not repairable by a better extremum test.** The exemption withheld 59.2 % of
the clip's firings on this mesh (2 680 917 cell-steps against 6 578 187 over 1950 corrector
calls). So **59 % of the cells the clip must bound on a polyhedral mesh are themselves stencil
extrema of psi.** A growing checkerboard has extrema at its own peaks; any rule that exempts
"a cell that is its stencil's extremum" hands the defect exactly the cells it needs.

| | hexahedral interface | polyhedral far field |
|---|---|---|
| clip, no exemption | +30.4 % volume error | defect REMOVED |
| clip + exemption | +1.03 % volume error | defect RETURNS at step 506 |

**RETRACTED: the apex detector as proposed.** "Exempt a cell whose own quadratic has a
stationary point inside the cell" would very likely withhold the same 59 % and fail G4 the
same way, because a checkerboard peak's fit also has a stationary point inside its cell. It is
a better DETECTOR of an extremum, and the measurement says detecting extrema is not the problem.

**What separates them is SCALE, and that is the next measurement.** The genuine extrema are the
medial axis of the signed distance — the apex of the distance cone, the far box corner. They
are WIDE: the apex sits 12.8 cells from the interface and the curvature of psi there is of
order 1/R. A spurious extremum is one cell wide, wavelength 2h. Two threshold-free properties
follow, neither tested:

1. **`|grad psi|` falls toward 0 at a genuine extremum** of a signed distance, because the
   field is smooth and stationary there; at a checkerboard peak it is of order 1 or larger,
   because the oscillation steepens the gradient. `|grad psi| = 1` is the normalisation the
   signed distance already carries, so this is a dimensionless statement, not a tuned length.
2. **A genuine extremum survives widening the stencil by one layer**; a mesh-scale oscillation
   does not. A layer count is discrete geometry, not a fitted coefficient.

**Also worth naming**: the hexahedral firing probe found the clip's WHOLE cost on hex in SIX
cells. A rule that names those six correctly satisfies both requirements at once, and the two
properties above are candidates for naming them.

**A caveat on arm 00002 that must not be lost.** Its zero-set error is still GROWING slowly —
4.46e-04 at step 1, 7.38e-04 at 500, 8.97e-04 at 649 — and 8.97e-04 sits just below the 1e-3
threshold. The global clip suppresses the growth by a factor 380; it does not stop it. Extend
END_TIME on that arm alone before the clip is called a fix. That is cheap.

## 3. The gate ladder for the fix, cheapest first

| # | gate | command / config | pass criterion | state 2026-09-09 |
|---|---|---|---|---|
| G0 | inertness, render | `--until generate_case` at HEAD~ and HEAD | only the new token line changes | **PASS** — only `clipRegion` / `clipKeepExtrema` appear |
| G1 | inertness, 2D hex | `config/popinet2D_La12000_N64.yaml` | metric CSV **bit-identical** to the pre-change run | **PASS** — 1563 steps, bit-identical |
| G2 | inertness, 3D hex | `config/popinet3D_La12000_hex_smoke4.yaml` (R/h 12.8) | metric CSV bit-identical | **PASS** — 78 steps, bit-identical |
| G3 | 4-rank parallel gate | `config/polyDroplet3Drefined_clipRegionGate.yaml`, `mpirun -np 4` | serial and np4 agree; see CLAUDE.md | code path exercised on 4 ranks; serial companion open |
| G4 | the defect is removed | `config/popinet3D_poly_sigma0_clipGate.yaml`, sigma = 0 | no failure to END_TIME; `zeroSetRadialL2/R` stays below 1e-3 | **FAILED — the candidate fails at step 506, the control at 527** |
| G5 | inertness at a MOVING interface | `config/popinet2D_clipRegionGate.yaml` (hex), `config/polyDroplet3Drefined_clipRegionGate.yaml` (poly) | volume error unchanged to 1 % | hex **+1.03 %** at the threshold; poly **PASS**, 0 cells bounded in 458 steps |
| G6 | the coupled polyhedral case runs | `config/popinet3D_La12000_poly_r12p8_mcs0195.yaml` in SI | COMPLETED at 1563 steps | not run |
| G7 | order is preserved | the poly ladder, section 5 | orders within +-0.3 of the hexahedral ladder | not run |
| **B0** | render inertness, the bound family | `--until generate_case` at HEAD~ and HEAD | only the new token lines change | **PASS** — only the four new entries appear |
| **B1** | bit-identity, 2D hex, EVERY bound sub-configuration | `config/popinet2D_clipRegionGate.yaml`, 8 arms, 1563 steps, np = 4 | every arm byte-identical to the pre-change study | **PASS** — 8/8 identical to `studies/popinet2D_clipRegionGate_preBound_20260910`; covers `none` AND `stencilBounds` in all four (region, keepExtrema) forms |
| **B3** | exact-field unit gate, cone bound | `leiaTestSLReconstruction`, unit square N = 64 and the SI Popinet mesh | 0 empty intervals, truth inside, no error increase, apex cells reached | **PASS** — 0 empty, truth inside to 1.1e-16, ZERO error increase; at 96 apex cells the monotone clip errs ONE CELL WIDTH and the cone errs 0 |
| **B5** | the NONLINEAR map's amplification | `leiaTestTransportSpectrum -mode growth`, 2D hex, amp 0.1h/1h/10h | per-step growth no worse than `none` | **FAILED for the coefficient-free arm** — L = 1 is 1.5-2.4x WORSE than no bound at every amplitude. Only the `stencil` L mode damps, and only at amp >= h |
| **B7a** | inertness at a MOVING interface, 2D hex, N = 64 | `config/popinet2D_coneBoundGate.yaml`, 12 arms, 1563 steps | interface metrics within 1 % of `none` | every interface metric IMPROVES 25-68 % at THIS rung. Retracted as a general result by B6 and B7b below |
| **B6** | ADVECTION ladder, hex and poly, one mesh per case | `2Dtranslation`, `2Dvortex`, `3Dshear` hex and poly, all arms on the IDENTICAL mesh | the bound must not degrade pure transport | **FALSIFIED** — 2.3-3.6x worse on uniform translation where L = 1 is exactly valid, 9-19x worse in strained flow, and the `stencil` mode LOSES THE PHASE on 3D shear hex (E_VOL_REL = 1.0000). The unbounded POLY arm diverges at step 198 and every bound prevents it, but `stencilBounds` does so with the best volume error. The falsified monotone clip beats the cone bound on every advection gate |
| **B7b** | RESOLUTION ladder, coupled 2D hex | `config/popinet2D_coneBoundLadderN128.yaml`, N = 64 and 128, 4420 steps at the fine rung | the advantage persists and the order is not degraded | **FALSIFIED, both criteria** — the eikonal error is 7.119e-03 at N = 64 and 7.066e-03 at N = 128 (order 0.01: a FLOOR, not a trend), the spurious-current order falls 0.24 to 0.03, the shape order 1.96 to 1.30, and the centroid error REVERSES to +119 % worse than `none`. The unbounded run reaches the bound's floor at about N = 256 |

G1 and G2 gate everything else, and they PASS. **G4 RAN AND FAILED**, exactly as its
pre-registered falsifier said it might. See section 2a.

The baselines for G1 and G2 were re-run at HEAD before any clip code existed, and both came
back bit-identical to the studies of 2026-09-08 — which also proves commit `7bf52eb`, the
amplification diagnostic, bit-inert at its default `fitProbeDisplacement (0 0 0)`. The
pre-change studies are preserved as `studies/*_preDiag_20260908` and `studies/*_B0`.

## 4. Stationary droplet 3D polyhedral — state and work

The case `stationaryDroplet3D` is ALREADY in SI (DOMAIN_LENGTH 0.006 m, DROPLET_RADIUS 1e-3 m).
The SI conversion of 2026-09-08 touched the Popinet cases only.

**Complete, on the cluster in `/work/scratch/tm83tomy/leia-curvature/studies/`:**

| study | config | steps | state |
|---|---|---|---|
| `polyDroplet3D_r13p8` | `polyDroplet3D_r13p8.yaml` | 3813 | COMPLETED |
| `polyDroplet3D_r18p9` | `polyDroplet3D_r18p9.yaml` | 6028 | COMPLETED |
| `polyDroplet3D_r25p6` | `polyDroplet3D_r25p6.yaml` | 9650 | COMPLETED |
| `polyDroplet3Drefined_r13p8` | `polyDroplet3Drefined_r13p8.yaml` | 3813 | COMPLETED |
| `polyDroplet3D_r13p8_clip` | `polyDroplet3D_r13p8_clip.yaml` | 200 | COMPLETED, short |

The uniform ladder is complete at R/h = 12.6 / 18.0 / 25.2 measured at the interface.

**Open work:**

1. **Re-pin `polyDroplet3D_r13p8`.** `N_CELLS = 84` describes h = 7.143e-05 m, and cfMesh built
   7.937e-05 m. `leia_refine.band_check` reports `pinRelError` 0.100 and suggests 76. The rung
   therefore runs 16 % below its own capillary law and is NOT dt-matched to its hexahedral
   twin. Re-pin to 76 and re-run before any dt-matched hex-poly comparison. Treat the present
   arm as valid for its own ladder only.
2. **The clip gate here is weak.** `polyDroplet3D_r13p8_clip` ran 200 steps with all nine
   metrics identical to seven digits. The interface does not move in that case, so it proves
   inertness for a QUIET interface only. G5 above must add a moving interface.
3. **The polyhedral stationary ladder is the well-balanced reference** for the fix: the clip
   must not change these numbers.

## 5. Translating droplet 3D polyhedral (Popinet) — state and work

**Every completed polyhedral translating arm predates the SI conversion, and every unclipped
one DIVERGED.**

| study | dt in its controlDict | steps | state |
|---|---|---|---|
| `popinet3D_La12000_poly_r12p8` | 2.55996e-04 (Popinet units) | 1154 | DIVERGED |
| `popinet3D_La12000_poly_r12p8_mcs0195` | 2.55996e-04 | 1026 | DIVERGED |
| `popinet3D_La12000_poly_r19p2` | 1.39347e-04 | 1308 | DIVERGED |
| `popinet3D_La12000_poly_r12p8_mcs0195_sigma0` | 2.55996e-04 | 1563 | ran, FAILED at step 509 |
| `popinet3D_La12000_poly_r12p8_mcs0195_sigma0_clip` | 2.55996e-04 | 1563 | COMPLETED, clean |
| `popinet3D_poly_sigma0_dtSweep` (SI) | 9.23685e-06 / 2 / 4 | 650/1299/2598 | COMPLETED, all failed |

The configs in `config/` now render SI. The runs above do not match them. That is not a void
setup: the SI equivalence is proven (2D ladder to 0.05 %, 3D smoke to 0.33 %). It means the
arms are not reproducible from the current configs, so **re-run the polyhedral translating
ladder in SI once the fix passes G1-G5**.

**Open work, in order:**

1. Fix first (section 2), gates G1-G5.
2. `popinet3D_La12000_poly_r12p8_mcs0195.yaml` in SI with the fix ON — G6.
3. The ladder `r12p8 / r19p2 / r25p6` in SI with the fix ON — G7. Report L1 and L2 only.
4. The hexahedral twins are already SI and correct. Do not re-run them.

## 6. Traps that have already cost time

1. **`zeroSetRadialL2`, `zeroSetRadialLinf`, `centroidError` and the `m2*` columns are ABSOLUTE
   LENGTHS in metres.** In SI they are 400 times smaller than in Popinet's units. A fixed
   threshold reads every SI run as clean. Divide by `DROPLET_RADIUS` first. This error misread
   the whole dt sweep on the first pass.
2. **Never report L_inf** for a verdict. L1 and L2 only.
3. **`MAX_CELL_SIZE` is a pin, not a mesh size.** On a polyhedral mesh `N_CELLS` only sets dt.
   Measure the built cell size with `leia_refine.band_check` and pin `N_CELLS` from it.
4. **cfMesh is not bit-reproducible under a change of scale.** The SI and the dimensionless
   polyhedral meshes both have 674 493 cells and demote the same 47622 cells to a linear fit,
   but they are not identical. Expect 0.3 % differences, not round-off.
5. **`quadraticPivotTol 0.3` demotes 47622 cells** of this mesh to a linear fit. That is the
   existing protection and it is not enough on its own.
6. **A finished case's `0/` is not its initial state.** Regenerate from `0.org` plus
   `leiaSetFields` for a bit-identity re-run.
7. **Rebuild the library, then relink every application**, and verify by timestamp and symbol.
8. **`scancel -u $USER` is forbidden.** The account is shared. Cancel by id from `.my_jobs`.

## 7. Where things are

- Cluster clone for this work: `/work/scratch/tm83tomy/leia-curvature`, with
  `WM_PROJECT_USER_DIR=$HOME/OpenFOAM/curvature-v2512`. Never touch `/work/scratch/tm83tomy/leia`.
- `profiles/slurm/config.yaml` and `run-studies.sbatch` are modified on the cluster on purpose.
  Use `git stash push -- profiles/slurm/config.yaml run-studies.sbatch; git pull --ff-only;
  git stash pop`.
- The dt sweep CSVs are local in `studies/popinet3D_poly_sigma0_dtSweep/` (24 MB). The full
  6 GB study stays on the cluster.
- The amplification diagnostic: `src/leiaLevelSet/semiLagrangian/
  uncachedQuadraticWeightedLeastSquaresReconstruction.C`, entry `fitProbeDisplacement`,
  default `(0 0 0)`. Documented in `workflow/README.md`.
