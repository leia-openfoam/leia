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

## 2. The fix to build (not started)

**A band-aware quasi-monotone clip of the semi-Lagrangian update.**

Requirements, all mandatory:

1. **It must not act on hexahedral meshes at all.** The user stated this directly. Hexahedral
   meshes run correctly and must stay bit-identical.
2. **No mesh-type switch.** The trigger must be a local geometric or algebraic property (for
   example the per-cell `Lambda`, or a one-sided stencil measure), never `mesh == poly`.
3. **The far field only.** A global clip is NOT inert at the interface: it costs +30 % volume
   error on the 2D hexahedral translating case (measured).
4. **Runtime-selectable and inert by default** (a dictionary entry in
   `fvSolution.levelSet.semiLagrangian` plus a case token), so every existing study is
   unaffected and the clip can be a study axis.
5. **Unstructured FVM only**: per cell, compact stencil, MPI-decomposable, one implementation
   for 2D and 3D.

The theory behind the clip: a reconstruction with `Lambda = 1` in a cell cannot create a new
extremum, because the update is then a convex combination of the stencil values. The clip
enforces that property where the fit does not have it.

## 3. The gate ladder for the fix, cheapest first

| # | gate | command / config | pass criterion |
|---|---|---|---|
| G0 | inertness, render | `--until generate_case` at HEAD~ and HEAD | only the new token line changes |
| G1 | inertness, 2D hex | `config/popinet2D_La12000_N64.yaml` | metric CSV **bit-identical** to the pre-change run |
| G2 | inertness, 3D hex | `config/popinet3D_La12000_hex_smoke4.yaml` (R/h 12.8) | metric CSV bit-identical |
| G3 | 4-rank parallel gate | any 3D arm, `mpirun -np 4` | serial and np4 agree; see CLAUDE.md |
| G4 | the defect is removed | `config/popinet3D_poly_sigma0_dtSweep.yaml` first arm, sigma = 0 | no failure to t = T_U; `zeroSetRadialL2/R` stays below 1e-3 |
| G5 | inertness at a MOVING interface | 2D hex translating, and the poly stationary rung | volume error unchanged to 1 % |
| G6 | the coupled polyhedral case runs | `config/popinet3D_La12000_poly_r12p8_mcs0195.yaml` in SI | COMPLETED at 1563 steps |
| G7 | order is preserved | the poly ladder, section 5 | orders within +-0.3 of the hexahedral ladder |

G1 and G2 gate everything else. Do not run G4 before they pass.

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
