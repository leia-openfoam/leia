# Plan: SLLS on cubic meshes with dynamic AMR — pre-print, version 1 (7 days)

## Context

The supervisor wants a first arXiv pre-print by 2026-11-30. Scope, as decided on 2026-11-23:

1. Semi-Lagrangian level set (SLLS) on **cubic cells with dynamic AMR** (blastAMR with load
   balancing, **two refinement levels**, 3D only). The refinement keeps a **uniform band of
   finest cells around the interface**, wide enough for the curvature with the quadratic
   reconstruction. Advection errors are reported.
2. D. Bothe's **algebraic rescaling** of psi with the SL fit, applied on a **deterministic
   trigger** (not every step, not case-specific), to keep |grad psi| away from 0 and bounded.
3. **Static curvature** for a sphere and an ellipsoid: uniform meshes, and AMR meshes with
   the uniform band.
4. **Stationary and oscillating droplet** (3D, coupled). The translating droplet is out.

What exists (MEASURED, STATUS 11.35 to 11.43): the ghost-cell reconstruction passes a static
2:1 refined-box gate (shape order 1.97); the 3D stationary droplet is accurate (jump error
0.03 %); the oscillating period is 1.4 to 2.7 % short. What does not exist: any dynamic-mesh
handling in the SL path, the rescaling on the SL line, curvature on refined meshes, the
article sections. Three branches are unmerged.

Supervisor's decisions already taken: own clone of blastAMR (plain fork `9c4f667`) built for
v2606; two levels from the start; 3D only; a new focused article; droplets on uniform meshes;
static refinement along the interface path as the fallback if dynamic AMR fails on day 4.

## 0. Rules that shape the work (CLAUDE.md)

- Every run is on Lichtenberg (account `special00005`, `--mem-per-cpu`, ledger `.my_jobs`).
- Every study has its prediction and falsifier in the config header before the submit.
- New behaviour is a new run-time-selectable (RTS) class, inert by default. A solver gets a
  new composition root as its own commit with a bit-identity run.
- Advection: at least 3 rungs, orders for L1 and L2, never L_inf; the standing regression
  set (2D hex, 3D shear hex, 3D shear poly) must stay byte-identical for inert changes.
- 4-rank check before any cluster study; serial against 4 ranks must agree.
- No free parameter: a width is a cell count derived from stencils and checked by an
  exact-solution gate; a threshold has one global dimensionless value, shown over a range.
- "Score every candidate with reinitialisation OFF": the pre-print reports every advection
  result **without and with** the rescaling. The supervisor's decision to include the
  rescaling in the version-1 method goes into the decision log.
- Lean output on every new study (`keep_fields: last`, `PURGE_WRITE 2`); check `cquota`
  before every campaign.

## 1. Work packages

### WP0 — Base (day 1, morning)

1. Merge into `development`, in this order, each fast-forward or with STATUS renumbering:
   `fix/sl-remote-patch-faces`, `work/droplets3D` (corrected Lamb period, 3D tables),
   `feature/lean-output`. Run `check_kb.py`, `diff CLAUDE.md AGENTS.md`.
2. Work branch `feature/slls-amr` from the merged `development`.
3. Lichtenberg: clone `tmaric/blastAMR` at `9c4f667` into `$HOME/OpenFOAM/blastAMR-v2606`
   (not the `peng` directory), build against OpenFOAM-v2606. Known risk: the private-member
   macro in `fvMeshRefiner.C:32-34` and the `Allwmake` branch check. Verify the seven
   `libamr*` libraries and run the tutorial `tests/testCases/hex3D` on 4 ranks with balance.
4. New cluster clone `/work/scratch/tm83tomy/leia-amr` (branch `feature/slls-amr`), built
   with `workflow/scripts/integration/build.sbatch`. `etc/leia-env.sh` gets the blastAMR
   library directory as a sibling build (same pattern as cfMesh and TwoPhaseFlow).

### WP1 — Dynamic mesh in the kinematic SL solver (days 1 to 3) — the critical path

Files: `applications/solvers/leiaSemiLagrangeLevelSetFoam/{leiaSemiLagrangeLevelSetFoam.C,createFields.H,Make/options}`,
new files in `src/leiaLevelSet/semiLagrangian/`, case templates of `3Dtranslation`,
`3Dshear`, `3Ddeformation`.

1. **Composition root (own commit, inert).** `createMesh.H` → `createDynamicFvMesh.H`;
   `mesh.update()` at the start of each step. Default `staticFvMesh` (a
   `constant/dynamicMeshDict` token `DYNAMIC_MESH staticFvMesh` in the three templates).
   Gate: the regression set byte-identical (`compare_metrics_csv.py --tol 0`).
2. **Rebuild after a topology change.** When `mesh.topoChanging()`: reset and construct
   again `slAdvection` (with its reconstruction, stencil map, centre tails, limiter field),
   the narrow band, the velocity model (so that U and phi are the exact prescribed field on
   the new mesh, not mapped values), and the phase indicator. No class is modified; the
   solver owns the `autoPtr`s. The stencil `MeshObject` is deleted by OpenFOAM on a topology
   change, so a rebuild is required, not optional.
3. **Second-order transfer of psi (new RTS family `slMeshTransfer`, base = `none` =
   OpenFOAM's mapping).** Model `quadraticFit`:
   - before `mesh.update()`: store the fit gradient `g`, Hessian `H` (`fitDerivatives`), the
     old cell centre and the old level as registered fields; blastAMR copies them to
     children and volume-averages them on unrefinement;
   - after: a child cell gets `psi = psi_p + g.d + 0.5 d.H.d` with `d` = new centre minus
     the mapped old centre; a merged parent gets `mean(psi_children) - (h_c^2/8) tr(H)`;
     an unchanged cell has `d = 0` and keeps its value bit for bit.
   - Exact for a quadratic field (unit gate below). OpenFOAM's default gives children the
     parent value, which is first order.
4. **Refinement criterion: the uniform band (derived, not tuned).** A registered indicator
   field marks the cells within `N_band` cell layers of the sign change (existing
   `redistancer::topologicalMask` / `narrowBand` layer count). blastAMR `errorEstimator
   fieldValue` refines where the indicator is 1 and unrefines where it is 0;
   `maxRefinement 2`; `refineInterval 1`; the decomposition constraint `refinementHistory`;
   `balance yes`, `allowableImbalance 0.2` (blastAMR default).
   - `N_band` (finest cells on each side) = 3 (face curvature reads layers l <= 3) + 1 (the
     quadratic fit stencil of the last of those) + 1 (no band cell reads a ghost value)
     + ceil(CFL x refineInterval) (interface travel between adaptations) = **6** for
     CFL <= 1. This equals the default `REFINE_BAND_CELLS 6` of the static script.
   - Level-1 buffer between level 2 and level 0: **2 cells** of level 1, so that no cell
     has a finer and a coarser neighbour in one stencil.
   - Both counts are confirmed by WP4 gate C2 (sweep 4, 5, 6, 8) before they are used.
5. **Levels from the registry.** The existing ghost class reads `cellLevel`/`pointLevel`
   from disk once; after an adaptation the files do not exist. The new class of WP2 reads
   the registered `labelIOList`s that blastAMR keeps under `polyMesh`.

### WP2 — Cached unit-cube reconstruction (days 1 to 2, parallel to WP1)

New class `cubicLatticeQuadraticWeightedLeastSquares` in
`src/leiaLevelSet/semiLagrangian/` (one line in `Make/files`, selected by the existing token
`SL_RECONSTRUCTION`). The existing ghost class stays unchanged as the reference.

- For a cube, the anchored weighted least-squares operator on the regular 26-point lattice
  is one **9 x 26 matrix** (5 x 8 in 2D), computed once for the unit cube. Weights are
  1/|d|^2, so they are the same on every lattice. Scaling: gradient rows by 1/h, second
  derivative rows by 1/h^2. Each cell stores only a permutation from stencil slot to
  lattice offset (26 labels) and its level.
- The lattice values come from the same three sources as in the ghost class: same-level
  centres, means of children at the parent point, ghost values from the coarser cell's
  polynomial (level-by-level fit, coarsest first, coefficient exchange between levels).
- **Not regular, solved per cell as today:** cells with a physical patch face in the stencil
  (data at h/2), and cells where a ghost position has no containing coarser cell. The
  supervisor's statement holds for all interior cells, also at transitions; it does not
  hold at the domain boundary.
- Adds `fitDerivatives` (needed by WP1.3 and WP3) and the registry levels (WP1.5), and
  supports three levels (0, 1, 2).
- Gates: (a) exact for a quadratic field on a 0/1/2-level mesh, serial and 4 ranks
  (extend `applications/test/leiaTestGhostCellReconstruction`); (b) the N = 64 refined-box
  run equals the ghost class to 1e-10 relative in every metric; (c) wall-clock per step
  against the ghost class on N = 128 (the profile of STATUS 11.36 gave 60 % of the time to
  the fit).

### WP3 — Algebraic rescaling with the SL fit (days 2 to 3)

New model `rescaling` in the existing `slSource` family
(`src/leiaLevelSet/semiLagrangian/source/`, one line in `Make/files`, token `SL_SOURCE`).
The family already runs after each advection step and receives the fit, so the solver is not
edited.

- Map, in the band cells only, psi_c multiplied by a positive factor (zero set unchanged):
  first order `psi/|g|`; second order `2 psi/(|g| + sqrt((|g|^2 - 2 psi n.H.n)^+))`, with
  `g`, `H` from the fit at the cell centre (the formula of
  `src/leiaLevelSet/profile/rescaledDistance.H`; the fit already has it as `signedOffset`).
- **Trigger (deterministic, global, one constant):** act in a step only if, over the band,
  `max | |g| - 1 | > theta`. `theta` is one dimensionless number for every case, default
  0.1 (token `SL_RESCALE_THETA`); the studies run theta = 0.05, 0.1, 0.2 and must show the
  same order. Tokens `SL_RESCALE_ORDER first|second`, default arm second.
- Record: the Eulerian form failed when applied every step and worked once (negative log,
  2026-10-31; its re-open clause "an event interval k > 1" is met by the trigger). A
  negative result here is possible and is reported as such.
- Gates: static circle/sphere from the eight initial fields of the static source gate
  (once: distance error order 2, zero set and volume unchanged); then the advection ladders
  of WP4 with the event count per run in the log.

### WP4 — Studies (pre-registered; configs under `config/slAmr/`)

| id | study | arms | rungs | read-out |
|---|---|---|---|---|
| U1 | transfer unit gate | linear, quadratic field; refine then unrefine; serial, 4 ranks with balance | one mesh, levels 0/1/2 | error at round-off for the quadratic field |
| U2 | cached operator | WP2 gates (a) to (c) | — | as in WP2 |
| D1 | serial against 4 ranks with balancing | 3D translation, AMR | effective N = 64 | every metric within 1e-9 relative |
| A1 | 3D translation, one way, T = 1 | uniform fine (control), dynamic AMR; rescaling off / on | effective N = 64, 128, 256 (512 if the cost allows) | shape, volume, band gradient L1 and L2 with orders; AMR error <= 1.5 x uniform at equal fine h; cell count and wall clock |
| A2 | 3D shear, reversed, T = 3 | same | same | shape and volume at T, gradient at T/2 and T |
| A3 | 3D deformation, reversed, T = 3 | same | same | same |
| R1 | rescaling threshold | theta = 0.05, 0.1, 0.2; order first / second | A2 ladder | same order for every theta; event counts |
| C1 | static curvature, uniform | sphere, ellipsoid (signed distance and implicit quadric); rescaled once / not | N = 32, 50, 80, 128 | L1 and L2 curvature error, order (existing configs `faceCurvatureSphere3D`, `faceCurvatureEllipsoid3D`) |
| C2 | static curvature, banded AMR | N_band = 4, 5, 6, 8; two levels | same fine h as C1 | error within 5 % of the uniform mesh at equal fine h for N_band >= the derived value |
| P1 | stationary droplet 3D | uniform; MEASURED (STATUS 11.42) | N = 60, 78, 102 | tables from the committed CSVs |
| P2 | oscillating droplet 3D | uniform; MEASURED; optional horizon 0.05 s | N = 60, 78, 102 | period against Lamb 8.232 ms, damping against 5.0 1/s |

Enablers for C1/C2: `N_CELLS_BASE` in the blockMesh templates of `3Ddeformation` and
`ellipsoidDroplet3D`; the fine h (not `cbrt(V/nCells)`) in `leiaTestMeanCurvature`'s
`DELTA_X`; the band-refined meshes from `workflow/scripts/leiaRefineHexMesh.py`
(`REFINE_SOURCE interface`, `REFINE_BAND_CELLS`, `REFINE_LEVELS 2`). Uniform results must
stay identical after these changes.

Every advection table also carries the row "original scheme without ghost cells on the
refined mesh" from STATUS 11.35 as the reason for the ghost values.

### WP5 — The pre-print (days 5 to 7; skeleton on day 1)

- New theme `sl-amr`: `docs/sl-amr/sl-amr-article/` (elsarticle preprint), data in
  `docs/sl-amr/sl-amr-article/data/{tables,figures}`; one row in `workflow/scripts/paths.py`.
- Sections: introduction; method (quadratic fit, ghost values at 2:1 transitions, unit-cube
  operator, mesh transfer, band criterion and its derivation, rescaling and trigger);
  verification (U1/U2, A1 to A3, R1, C1/C2, P1/P2); limitations (translating droplet,
  boundary cells, period drift, cubic cells only); conclusions.
- Reuse the still-valid method text of
  `docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`
  (sections at lines 214 to 716); do not copy the five superseded statements listed in the
  knowledge-base note `studies/sl-quadratic-pre-print`.
- Tables and figures by committed scripts (`workflow/scripts/make_sl_amr_tables.py`, one
  script, reads the study CSVs); `make article-sl-amr` (latexmk) and `make arxiv-sl-amr`
  (flattened sources and a tar file).

## 2. Schedule

| day | date | work | exit condition |
|---|---|---|---|
| 1 | Mon 23 | WP0 complete; WP1.1 (inert root) and its identity run; WP2 class written; article skeleton; C1 submitted | blastAMR v2606 runs its tutorial on 4 ranks; regression set identical |
| 2 | Tue 24 | WP1.2 to 1.5; WP2 gates; WP3 class written; C2 enablers and submit | U2 passes; first AMR run completes on 4 ranks |
| 3 | Wed 25 | U1, D1; WP3 static gate; A1 first rungs | U1 exact, D1 within 1e-9 |
| 4 | Thu 26 | **Decision point.** A1 to A3 ladders, R1 | pass: continue. fail: fallback = static two-level band along the interface path (existing script), same studies |
| 5 | Fri 27 | remaining rungs, optional P2 horizon, all tables and figures | every table generated by script |
| 6 | Sat 28 | text: method and verification | full draft |
| 7 | Sun 29 | text: limitations, checks, arXiv build; STATUS, METHOD, knowledge base | PDF and arXiv tar; `check_kb.py` passes |
| buffer | Mon 30 | supervisor's read, corrections | — |

Parallel streams: WP1 (critical path, main session), WP2, WP3 and the C1/C2 enablers as
agents with narrow briefs and incremental commits, each on the work branch.

## 3. Risks and the answer to each

1. **blastAMR does not build or run on v2606.** Day-1 exit condition. Answer: fix the
   build in the own clone (a branch of the fork); if not solved by day 2, fallback mesh.
2. **Transfer or rebuild is wrong in parallel after load balancing.** U1 and D1 run before
   any ladder. Answer: balancing off as an intermediate state; it is then stated as a limit.
3. **Two levels break the ghost construction.** U2 gate (a) on a 0/1/2 mesh before any
   study. Answer: one level, stated.
4. **Rescaling does not help or harms.** It is an arm, not the base: every result has the
   arm without it. A negative result is reported.
5. **File quota.** Lean output; AMR runs write two times only; `cquota` before each batch.
6. **Oscillating period.** Reported as measured, with the corrected reference.

## 4. Verification of the plan's result

```bash
# inert composition root and every later inert change
workflow/scripts/compare_metrics_csv.py A.csv B.csv --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME
# unit gates (serial and mpirun -np 4, on Lichtenberg)
leiaTestGhostCellReconstruction          # extended: 0/1/2 levels, cached class, transfer
# state of every case
workflow/scripts/foam_log_state.sh <log>
# studies
sbatch -J leia-amr --export=ALL,STUDY=slAmr/<study> run-studies.sbatch
# record and article
python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base
diff CLAUDE.md AGENTS.md
make article-sl-amr && make arxiv-sl-amr
```

End-to-end check: the AMR ladders A1 to A3 complete on 32 ranks with load balancing; the
orders are printed next to the errors; the uniform control at equal fine h is in the same
table; the curvature on the banded mesh equals the uniform one within the stated tolerance;
the PDF builds from committed scripts and configs only.

## 5. Record

Each result goes into `STATUS.md` (sections from 11.44), `METHOD.md` (the band width, the
buffer, theta, the reconstruction of the AMR line), and the knowledge base (a new note
`concepts/sl-dynamic-amr`, log lines in `concepts/sl-ghost-cell-reconstruction` and
`concepts/algebraic-rescaling-eulerian`, decision and negative logs) in the same commit as
the result.

## 6. Decisions after the approval (2026-11-23, supervisor)

1. **The scheme of version 1 is the least-squares fit** (`C-Cp/L2+g`, class
   `cubicLatticeQuadraticWeightedLeastSquares`). It is also the route to general
   hexahedral meshes: on perturbed hexahedra at Co 1 it amplifies about ten times less
   than the interpolating scheme (matrix symbol, MEASURED).
2. **Roadmap, not version 1: the interpolating scheme on cubes** (`C-Cp/I+g`, the
   tensor-product quadratic on the 27 lattice points). The matrix symbol gives it a
   growth factor of 1 at Co >= 0.3 on a cubic 2:1 transition. It has never run in the
   solver.
3. **The Courant number stays at most 1.** The solver evaluates the polynomial of the
   arrival cell. The literature survey of 2026-11-23 states that this is unstable above
   CFL 1 and that larger steps need a stencil centred on the departure point.
4. **The stability analysis is numerical** (matrix symbol on our meshes), not symbolic.
   It reports the growth factor and the largest norm of the powers of the update matrix.
5. **Roadmap candidate (proposed, not decided): a restriction that is exact for a
   quadratic.** The mean of the eight children is the centre value plus (h_c^2/8) tr(H).
   The mesh transfer corrects this term; the ghost restriction does not. A corrected
   restriction is a method change with its own gate.
6. **The advection ladders run at CFL 0.5**, and every table states the time-step law:
   the semi-Lagrangian error has the form h^(p+1)/dt, so a smaller time step at fixed h
   increases the error.
7. **References.** The survey marks several of its references as not verified. Each
   reference is checked against its source before it enters the article.
