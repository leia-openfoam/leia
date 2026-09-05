# leia test-suite — Snakemake workflow

KISS Snakemake replacement for the old `pyFoamStudy` parameter-study flow. It
generates parameter variations from the existing `cases/<Case>.parameter` files
+ `@!TOKEN!@` templates, runs each variation through OpenFOAM, and aggregates
the per-case CSVs into one study database. The **same** workflow runs locally
(`mpirun`) or on SLURM (one `sbatch` per case) — you only switch the profile.

## Layout

```
workflow/Snakefile          generate_case -> mesh -> [decompose] -> solve -> aggregate
workflow/scripts/
  foam_param.py             parse <Case>.parameter + default.parameter (+ which @!TOKEN!@ each case uses)
  decompose.py              np -> (nx ny nz) decomposition (2D/1D aware)
  materialize.py            render templates, write np-driven decomposeParDict + case_params.json
  aggregate.py              join per-case CSVs + parameter vector -> <study>_database.csv
config/config.yaml          which case/mesh/mode/np + scope (smoke subset by default)
config/bulkVortex.yaml      FAST reversed-vortex suite (N<=128, 6 models, ~5 min; export_slides: false)
config/bulkVortexHighRes.yaml  opt-in deep convergence (+N=256, 6 models, ~46 min; regenerates the deck)
config/phaseIndicatorConvergence.yaml  geometric vs detrixheAslam (N<=128 by design)
config/faceCurvatureDroplet2D.yaml  FACE-CENTERED curvature convergence: kappa_f as the
                             CSF force applies it (active snGrad(alpha) faces, static
                             circle, exact SDF), every curvature model with/without the
                             stabilized foot point (leiaTestMeanCurvature, serial,
                             seconds; figure + orders into the method-comparison theme)
config/stationaryDropletStableFoot.yaml  parasitic currents with the SECOND-ORDER face
                             curvature (curvatureExtension stabilizedFootPointFace +
                             faceCurvatureSource registered): N=256,512 (both >20
                             cells/radius), np 4; observable = per-step maxMagU
                             evolution vs the arithmetic-delivery baseline.
                             MEASURED: N=256 blows up at t=0.0348 (arithmetic twin
                             t=0.0526) -- the reinit-free psi-profile feedback, not
                             the delivery, is the active channel.
config/stationaryDropletStableFootFiltered.yaml  the profile-control gate: same
                             delivery + psiFilter biharmonicBand theta 0.05,
                             N=64..512 -- does minimal per-step band damping close
                             the parasitic feedback loop? (evolution CSVs;
                             make_droplet_evolution_fig.py)
                             MEASURED: N=64/128/256 all reach the full t=0.1 with
                             the band |grad psi| pinned and the curvature error
                             flat -- the loop is closed; residual slow m=1-like
                             drift, weakening under refinement.
config/stationaryDropletStableFootLong.yaml  horizon extension of the filtered
                             gate: N=128,256 to t=0.3 (np 8) -- does the residual
                             drift saturate, oscillate, or grow?
config/staticExtension.yaml  static t=0 extension verification: e=|n.grad(Uext)| vs h
                             (solver: leiaTestVelocityExtension; 7 models x 2 div schemes)
config/steadyVortex2D.yaml   NON-REVERSING stress test: steady vortex (oscillation off,
                             T=3), 6 models; error-vs-time crossover t*; shape error vs a
                             marker-traced reference (scripts/marker_ref.py) since there
                             is no analytic final interface without reversal
config/steadyVortex2DHighRes.yaml  + N=256 tier (opt-in; refreshes the deck's steady_* figures)
config/bulkVortexSL.yaml     SEMI-LAGRANGIAN solver (solver: leiaSemiLagrangeLevelSetFoam)
                             on the reversed vortex; sweeps SL_RECONSTRUCTION
                             (linearTaylor/linearWeightedLeastSquares/quadraticTaylor/quadraticWeightedLeastSquares) x CFL{0.5,1.0} -- the SL
                             analog of the VELOCITY_EXTENSION model sweep. plots.py emits the
                             reconstruction-convergence figure + an sl_vs_extension cross-study
                             overlay (reads the bulkVortexHighRes velocity-extension study).
config/bulkVortexSLHighRes.yaml  + N=256 tier (opt-in; refreshes the deck's sl_* figures)
profiles/local/config.yaml  executor: local   (mpirun; jobs x np = 24 ranks, no oversubscription)
profiles/local20/config.yaml executor: local  (20 of 24 cores; leaves 4 free for interactive use)
profiles/local8/config.yaml executor: local   (8 of 24 cores, ONE arm at a time: the tier to use
                             while another session holds the rest of the machine; memory, not
                             cores, is what bounds 3D here -- see the profile's header)
profiles/slurm/config.yaml  executor: slurm    (one sbatch per case; srun + module env)
studies/<study>/            generated cases + <study>_database.csv  (git-ignored)
```

Named studies write to their own `studies/<study_name>` — never point configs with
different axes at one study dir (cartesian-product indices would remap and stale
cases would be silently reused). Cost scaling: +1 resolution level in 2D = 8x
(4x cells x 2x CFL steps) — that is why N=256 is a separate opt-in config.

## Install

```bash
python3 -m pip install --user --break-system-packages "snakemake>=8" snakemake-executor-plugin-slurm
export PATH="$HOME/.local/bin:$PATH"
```

OpenFOAM must be sourced for the run steps. Locally this is done by
`env_preamble` in `config/config.yaml` (defaults to OpenFOAM-v2512 — adjust to
your install, or set it to `""` if you already source OpenFOAM in your shell).

## Run

Edit `config/config.yaml` (`case`, `mesh`, `mode`, `np`, scope), then:

```bash
cd <repo>
snakemake --workflow-profile profiles/local                 # this machine, mpirun
snakemake --workflow-profile profiles/slurm                  # cluster, one sbatch/case
snakemake --workflow-profile profiles/local -n               # dry-run: preview job count
```

Override anything on the command line, e.g. a different case in parallel:

```bash
snakemake --workflow-profile profiles/local \
  --config case=3Dshear mesh=hex mode=parallel np=8 study_name=shear-par
```

### Scope (smoke vs full)

`config/config.yaml` defaults to a small **smoke** grid via `axes_override` +
`collapse_other_axes: true`. For the full cartesian product from
`cases/<Case>.parameter`:

```yaml
axes_override: {}
collapse_other_axes: false
```

Always preview the size first with `-n` — full grids are hundreds–thousands of
cases (e.g. 3Drotation = 1296, 3Ddeformation = 31104). `foam_param.py` only
sweeps tokens the case's templates actually use, so unused axes
(`NEIGHBOURS`/`NARROWBAND` for 3Drotation) don't create duplicate cases.

### Comparing phase-indicator methods

Every case sweeps `PHASE_INDICATOR ( geometric detrixheAslam )` (the `@!PHASE_INDICATOR!@`
token in `fvSolution.template`), so a study compares the polygon-clipping
`geometric` indicator against the analytic tetrahedral-fill `detrixheAslam`
indicator (and `heaviside` if added to a `.parameter`). The resulting
`PHASE_INDICATOR` column in the database lets you diff their error metrics. The
default smoke config sweeps both; drop it from `axes_override` to fix one.

## Notes

- Advection studies also aggregate the per-step CSV row nearest **t = T/2**
  (`half.*` database columns; `gradientErrorBandHalf` etc. in `*_errors.csv`):
  the state at maximal deformation, before any reversal cancellation. The
  `maxdef_convergence.png` figure contrasts it with the final-time reading —
  the reversed benchmark's final row credits `none` with error cancellation
  that no extension model receives.
- `mesh` ∈ {hex, perturbed, poly, hexRefined, polyRefined}; `poly` and
  `polyRefined` are 3D-only and use `cases/<Case>_poly.parameter`. `perturbed`
  adds `-fluxCorrection`. `hexRefined` / `polyRefined` are statically refined
  around the interface by `workflow/scripts/leiaRefineHexMesh.py` /
  `leiaRefinePolyMesh.py` (need `REFINE_LEVELS >= 1`; see the static-refinement
  entries below). If the
  case exposes `N_NON_ORTHOGONAL_CORRECTORS`, materialisation also enforces a
  minimum of 8: the frozen-circle velocity sweep is converged at 8--64 on
  deterministic 10%-perturbed N=32,64,128 meshes, whereas 1 is insufficient.
- `np` is the single source of truth: it regenerates `system/decomposeParDict`,
  drives `mpirun -np {np}` / SLURM `--ntasks`, so the rank count and the
  decomposition can never disagree.
- `3Dtranslation` runs its `system/init_End/init_End.<mode>` automatically.
- Non-parametric utility tests (1DredistanceTest, 2DgradTest):
  `snakemake --workflow-profile profiles/local utilities`
  (requires the `leiaTestRedistance` / `leiaTestGradScheme` apps to be built).
- On SLURM, edit the account and the module/`source` lines in
  `profiles/slurm/config.yaml` for your cluster.

### Capillary pressure-compatibility workflow

The canonical entry point for the frozen-circle pressure studies is Snakemake:

    snakemake -s workflow/Snakefile.pressure-compatibility \
      --workflow-profile profiles/local

Its DAG runs the non-orthogonal-correction convergence sweep, the constant-
curvature CSF/pressure-potential comparison, the paired corrected versus
uncorrected `snGrad`/pressure-Laplacian gate, the physical/constant `rAUf`
oracle, the pressure-algebra tolerance sweep, and the GAMG/PCG solver gate.
The similarly named Make targets are thin aliases only; they do not own study
logic or freshness.

## 3D semi-Lagrangian convergence (quadraticWeightedLeastSquares)

Plain snakemake, one config per case (both sweep N=32/64/128, CFL 0.5 & 1.0):

    snakemake --workflow-profile profiles/local --configfile config/3DdeformationSL.yaml
    snakemake --workflow-profile profiles/local --configfile config/3DshearSL.yaml

The local profile caps the global `tasks` budget to the core count, so concurrent
np=8 solves are bounded to 3 (no oversubscription) — which is why plain snakemake is
enough. 32/64 are fast; each 128^3 case is ~1.5–2 h (the 3D departure-foot search
dominates).

Memory (this box is 15 GB; the `quadraticWeightedLeastSquares` pseudo-inverse cache is single precision):
  * deformation 128^3 (2.1M cells) fits at np=8 — give a 128^3 case the whole box by
    running one solve at a time:

        snakemake ... --configfile config/3DdeformationSL.yaml --resources tasks=8

  * shear 128^3 (128×128×256 = 4.2M cells) OOMs at np=8 on 15 GB; use fewer ranks:

        snakemake ... --configfile config/3DshearSL.yaml --config np=4 --resources tasks=4

    or run it on a >~24 GB node.

Resuming a long 128^3 run (snakemake restarts an interrupted case from t=0) — resume the
solver directly instead:

    cd studies/<study>/<case_dir>
    sed -i 's/startFrom       startTime/startFrom       latestTime/' system/controlDict
    mpirun -np <np> leiaSemiLagrangeLevelSetFoam -parallel   # continues from latestTime
    reconstructPar -withZero                                 # then re-run snakemake to aggregate

Regenerate the deck figures: `python3 workflow/scripts/make_sl_3d_fig.py`

## Geometrically redistanced level set (leiaRedistancedLevelSetFoam)

The third research line: Eulerian psi advection + criterion-gated geometric
redistancing from the phase indicator's own least-squares planes. Theme
`geometrically-redistanced-levelset`; results agglomerate into
`docs/geometrically-redistanced-levelset/grl-level-set-article/data/` and both
decks (`geometrically-redistanced-level-set.template.html`, `geometrically-redistanced-level-set-negative-results.template.html`).

Studies (each config header documents its axes and purpose):

    # static redistancing gate (circle vs analytic SDF; the acceptance test
    # every fill change must pass: post-event band L1 <= pre-event)
    snakemake --workflow-profile profiles/local --configfile config/redistanceCircle2D.yaml

    # trigger ablation: every-step (interval) vs gradPsiThreshold criterion
    snakemake --workflow-profile profiles/local --configfile config/vortexTriggerGRL.yaml

    # advection: reversed 2D vortex / 3D shear / 3D deformation,
    # REDISTANCER axis = [noRedistancing, PDE, planeFootWave, anchoredEikonal]
    snakemake --workflow-profile profiles/local --configfile config/bulkVortexGRL.yaml
    snakemake --workflow-profile profiles/local --configfile config/3DshearGRL.yaml
    snakemake --workflow-profile profiles/local --configfile config/3DdeformationGRL.yaml

or all of them: `make studies-grl` (repo root). Sweep tokens live in
`cases/default.parameter` (`REDISTANCER`, `REDIST_TRIGGER`, `REDIST_THRESHOLD`,
`REDIST_INTERVAL`) and render into the `levelSet.redistancer` subdict of the
case `fvSolution.template`s. Report scripts: `make_redistance_table.py`
(static gate table + figure), `make_grl_fig.py` (advection convergence).

Unit/acceptance tests outside snakemake: `leiaTestRedistance` (single event,
band errors vs the analytic SDF from `levelSet.implicitSurface` -> CSV),
`leiaTestLevelSet` (planar invariance per model contract),
`cases/1DredistanceTest/Allrun_variants.sh` (all models on the 1D plane).
The per-case `Allrun.sh`/`Allclean` scripts are DEBUG conveniences for a
single variant; snakemake is the canonical, reproducible path.

CALIBRATION NOTE (measured, 2026-07-17): the automatic gradPsiThreshold
default (h/L)^2 lies far BELOW the geometric fill's achievable post-event
band floor (~5e-3 at N=64 on the static gate) -> the criterion fires every
step and per-event interface displacements compound over long runs
(bulkVortexGRL T=8: volume drift for every redistancer). Set REDIST_THRESHOLD
explicitly above the measured floor for long advection studies.

## Static local refinement around the interface (mesh: hexRefined | polyRefined)

Refine ONLY a band around the interface, sized from the stencils the method
uses, and leave the far field coarse. Pre-processing only -- the solver is not
changed and `fvSolution.levelSet` is untouched -- realised by two stdlib-Python
drivers over existing OpenFOAM/leia/cfMesh apps, called from the `mesh` rule
(under `profiles/slurm` that rule is already a serial `case_pre` job):

    workflow/scripts/leiaRefineHexMesh.py    blockMesh, then REFINE_LEVELS passes of
                                             [0/ := 0.org; leiaSetFields; topoSet (seed
                                             0 < alpha < 1 = the snGrad(alpha) support,
                                             face dilations added until the psi on the
                                             current mesh proves REFINE_BAND_CELLS
                                             complete fine layers at the worst point);
                                             refineHexMesh (hexRef8, cellLevel/pointLevel
                                             PERSISTED -> 2:1 across passes)]
    workflow/scripts/leiaRefinePolyMesh.py   pMesh, then REFINE_LEVELS passes of
                                             [0/ := 0.org; leiaSetFields; psi = 0
                                             iso-surface as STL; pMesh with
                                             surfaceMeshRefinement { additionalRefinementLevels i;
                                             refinementThickness REFINE_BAND_CELLS*h }]
                                             (cfMesh has no in-place refiner: re-meshed)
    workflow/scripts/leia_refine.py          shared: runner, 0/ reset, ascii readers,
                                             the band check, refinement.csv/refinedBand.csv
    workflow/scripts/check_refined_band.py   re-check an existing case

Both end with `0/ := 0.org; leiaSetFields` on the FINAL mesh: every field in `0/`
is either its `0.org` value or freshly computed there -- nothing mapped through a
refinement survives (a mapped alpha is a smeared alpha). The band check then
FAILS (exit 2, the mesh rule fails) when an interface cell lies outside the fine
region, when the first coarse cell centre is < 4 fine cells from the interface,
or when `N_CELLS` -- the capillary-dt handle, `adjustTimeStep no` -- does not
encode the FINE spacing (hex: `DOMAIN_LENGTH/N_CELLS` vs the built spacing;
poly: the pin the driver prints as `N_CELLS_suggested`).

Tokens (`cases/default.parameter`, inert defaults): `REFINE_LEVELS 0`
(halvings of the near-interface size; >= 1 only with the refined mesh kinds,
materialize asserts the pairing), `REFINE_BAND_CELLS 6` (complete fine layers
each side at the worst point -- measured, since face dilation is Manhattan growth
and buys only ~1.15 fine cells per dilation along a sphere's diagonals; the
stencil minimum is 4), `REFINE_SOURCE interface | ball` (ball = the control
that refines the whole droplet interior). Derived: `N_CELLS_BASE =
N_CELLS/2^REFINE_LEVELS` is what blockMesh renders.

Studies (each config header carries the pre-registered prediction and gate):

    # G-1: uniform N=30 arm run before/after the template change; CSVs cmp-identical
    snakemake --workflow-profile profiles/local8 --configfile config/stationaryDroplet3DbitIdentity.yaml
    # G0/G1/G2/GC: constant-curvature (2/R) well-balanced gate on the refined mesh,
    # interface band + ball control, np 4; its serial twin; its uniform control
    snakemake --workflow-profile profiles/local8 --configfile config/stationaryDroplet3DrefinedWB.yaml
    snakemake --workflow-profile profiles/local8 --configfile config/stationaryDroplet3DrefinedWBserial.yaml
    snakemake --workflow-profile profiles/local8 --configfile config/stationaryDroplet3DuniformWB.yaml
    # G3/G4 (cluster): refined ladder N = 60/76/96/120, its two controls, uniform twins
    snakemake --workflow-profile profiles/slurm --configfile config/stationaryDroplet3Drefined.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/stationaryDroplet3DrefinedL2.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/stationaryDroplet3DrefinedBall.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/stationaryDroplet3Duniform.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/stationaryDroplet3Duniform120.yaml
    # P0-P2: polyhedral twin of polyDroplet3D_r13p8 (pin N_CELLS from the driver first)
    snakemake --workflow-profile profiles/slurm --configfile config/polyDroplet3Drefined_r13p8.yaml

Laptop poly meshes: cfMesh needs jemalloc preloaded on glibc 2.39, scoped to pMesh only
(`LEIA_PMESH_PRELOAD=/path/libjemalloc.so.2` in `env_preamble`, a bare path because
`--config` strips quotes; a global `LD_PRELOAD`
segfaults the MPI solver at startup with an empty log -- measured).

Curation: `make_refined_mesh_table.py` (refinement.csv + refinedBand.csv +
checkMesh -> `refined_mesh_stats.csv`), `make_refined_equivalence_table.py`
(refined vs uniform at matched N_CELLS, matched t, equal steps; L1 and L2 only,
never L_inf), `make_refined_mesh_fig.py` (mid-plane slice coloured by cellLevel).

## Popinet's translating droplet in 3D on polyhedral meshes (mesh: poly)

`cases/popinetTranslating3D` is the 3D twin of `cases/popinetTranslating2D` (Popinet 2009,
Sec. 6.2.2): D = 0.4 in a 2 x 1 x 1 box, inflow U = 1 at x = 0, pressure outlet at x = 2,
free-slip sides, equal densities and viscosities (La = 12000, We = 0.4, sigma = 1), horizon
T = D/U = 0.4. The box STL is generated with named solids so that cfMesh produces exactly
the patches the fields declare (`workflow/scripts/make_box_stl.py --xlen 2 --ylen 1
--zlen 1`: solids `inlet`, `outlet`, `walls`); `meshDict.template` only sets their types.
A `blockMeshDict.template` for a hexahedral twin is included.

Resolution is set at the INTERFACE: cfMesh's dual cells there measure 2^(-1/3) x maxCellSize
(measured on every polyhedral rung of the stationary droplet), so `MAX_CELL_SIZE =
h/0.7937` and `N_CELLS` (the capillary-dt handle) is pinned to `1/h`; the pin is verified on
the built mesh with `check_refined_band.py --mode poly`.

    # 4-rank smoke (78 steps) on the coarsest mesh -- the gate before the cluster
    snakemake --workflow-profile profiles/local8 --configfile config/popinet3D_La12000_poly_smoke4.yaml \
        --config env_preamble="...; export LEIA_PMESH_PRELOAD=$HOME/miniconda3/lib/libjemalloc.so.2"
    # the ladder: R/h = 12.8 / 19.2 / 25.6 at the interface (N = 64 / 96 / 128; ~0.7M / 2.4M / 5.6M cells)
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r12p8.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r19p2.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r25p6.yaml

Read-out as for the 2D reproduction (`make_popinet_table.py`): maximum over time of the L1
and L2 norms of |u - U0|/U0, plus volume, shape, Laplace jump and band curvature error at T.

### Geometric admissibility of the quadratic fit (`SL_QUAD_PIVOT_TOL`, 2026-09-05)

The first polyhedral Popinet-3D run diverged at step 3 while its hexahedral twin
(`config/popinet3D_La12000_hex_smoke4.yaml`) was clean. Root cause, in the library: the
cell-point-cell quadratic stencil of cfMesh's wall boundary-layer and size-transition
polyhedra offers as few as 10 point-neighbours for 9 coefficients, so the weighted normal
matrix has a scaled condition number of 1e7-1e12 and the fit carries curvature
coefficients of 1e3-1e8: exact at the cell centre, wrong by coefficient x (U dt)^2 at the
departure foot. The stationary droplet (|u| ~ 1e-5) and the kinematic gates (velocity zero
on the walls, where these cells sit) cannot see it; a droplet translating at U = 1 through
those cells can. `uncachedQuadraticWeightedLeastSquaresReconstruction` now decides once
per mesh, from stencil POSITIONS alone, whether a cell's quadratic fit is admissible: the
smallest Cholesky pivot of the Jacobi-scaled weighted normal matrix must exceed
`quadraticPivotTol`; below it the cell uses the linear fit, and a linear stencil that fails
the same test keeps the cell value.

    SL_QUAD_PIVOT_TOL   0.3     # default; 0 = never fall back (the pre-fix behaviour, the control arm)
    SL_WRITE_FIT_ORDER  false   # true: write slFitOrder (2/1/0) and slFitPivot once, for the census

The default was read off the pivot census of five meshes (`sl_fit_pivot_census.py`, table
`docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/sl_fit_pivot_census.csv`):
band cells >= 0.74 on every mesh, hex >= 0.64 everywhere (hanging-node cells included),
polyhedral interior >= 0.50, while the degenerate stencils form a tail from 0 to ~0.3 with
an almost empty 0.3-0.4 bin. Stability needs the admitted condition c to satisfy
c (U dt / h)^2 < 1 (1e-2, c ~ 1e4, still diverged at step 16 at U dt / h = 0.016); 0.3 admits
c <~ 100, i.e. U dt / h < 0.1, and demotes ~7 % of a uniform cfMesh box, all in the wall
boundary layer, none within 12h of an interface. Gates on the 0.3 binary: hex bit-identity
(`config/stationaryDroplet3DbitIdentity.yaml`) cmp-identical to the pre-fix CSV; the
polyhedral refined smoke re-run on its own mesh and decomposition; the 78-step polyhedral
Popinet-3D smoke (`popinet3D_La12000_poly_smoke4`). Diagnostic configs kept as the record
of the isolation: `popinet3D_La12000_poly_dump4*.yaml` (4-step field dumps; `_oc1` one
outer corrector + `SL_WRITE_FIT_ORDER true`, `_qr` Householder QR, `_serial` one rank) and
`popinet3D_La12000_poly_smoke4_ccTrace.yaml` (cell-centred trace) all reproduced the
divergence, so neither the solver, the trace nor the rank count was the cause.

Census recipe (serial, in a COPY of a rendered case; the script prints the histogram by
class and appends one CSV row per mesh):

    cp -r <case>/{0,0.org,constant,system} <copy> && cd <copy>
    foamDictionary -entry levelSet/semiLagrangian/writeFitOrder -set true system/fvSolution
    foamDictionary -entry endTime -set <deltaT> system/controlDict
    leiaSemiLagrangianLevelSetTwoPhaseFoam > log.census
    python3 workflow/scripts/sl_fit_pivot_census.py . --label "<mesh>" --csv <table.csv>
