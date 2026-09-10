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

## Popinet's translating droplet: the SI parameter set

Popinet (JCP 228, 2009, Sec. 6.2.2) writes the benchmark in dimensionless form
(rho = sigma = U = 1, D = 0.4 of the box height). This repository ran it that way until
2026-09-08 and now runs the DIMENSIONAL twin at the same dimensionless groups, so that
every dictionary holds a real droplet and a real speed:

| quantity | SI value | fixed by |
|---|---|---|
| D (droplet diameter) | 1.0e-3 m (R = 5e-4) | chosen |
| rho, both phases | 1000 kg/m^3 | chosen (density ratio 1, as in Popinet) |
| nu, both phases | 1.0e-6 m^2/s | chosen (viscosity ratio 1); mu = 1.0e-3 Pa s |
| sigma | 0.012 N/m | La = sigma D/(rho nu^2) = 12000 |
| U | 0.0692820323 m/s | We = rho U^2 D/sigma = 0.4 |
| box | 5 x 2.5 x 2.5 mm | height = D/0.4, length = POPINET_XLEN heights |
| T_U = D/U | 0.01443375673 s | one diameter of travel |
| CAPILLARY_DT_COEFF | 4.729265146e-3 | 0.2323 sqrt((rho1+rho2) H^3/(2 pi sigma)) |

Five dimensional quantities carry two groups, so three are free: D, rho and nu are the
choice (water at 20 C). Re = sqrt(La We) = 69.28 and Oh = 1/sqrt(La) = 9.13e-3 follow, so
every number Popinet reports still applies. Scale factors from his units: length x 2.5e-3 m,
velocity x 0.0692820323 m/s, time x 0.03608439182 s.

`workflow/scripts/popinet_si.py` owns the set:

    python3 workflow/scripts/popinet_si.py print                  # the set, dt per rung
    python3 workflow/scripts/popinet_si.py yaml                   # the axes_override lines
    python3 workflow/scripts/popinet_si.py check config/popinet*.yaml   # exits 1 on a mismatch
    python3 workflow/scripts/popinet_si.py rewrite config/x.yaml  # convert a dimensionless config

`check` verifies the density, viscosity, surface tension, speed, radius, box, horizon and
the polyhedral cell-size pin of a config, and it is the gate to run after any edit.
`POPINET_XLEN` is an ASPECT RATIO (box length in box heights), not a length: both
`blockMeshDict.template`s compute `xlen` with `#eval{ POPINET_XLEN * ylen }`.

## Popinet's translating droplet in 3D on polyhedral meshes (mesh: poly)

`cases/popinetTranslating3D` is the 3D twin of `cases/popinetTranslating2D`: D = 1 mm in a
5 x 2.5 x 2.5 mm box, inflow U = 0.0692820323 m/s at x = 0, pressure outlet at x = 5 mm,
free-slip sides, equal densities and viscosities (La = 12000, We = 0.4, sigma = 0.012 N/m),
horizon T_U = D/U = 0.01443375673 s. The box STL is generated with named solids so that
cfMesh produces exactly the patches the fields declare
(`workflow/scripts/make_box_stl.py --xlen 5e-3 --ylen 2.5e-3 --zlen 2.5e-3 --out
cases/popinetTranslating3D/box5x2p5x2p5mm.stl`: solids `inlet`, `outlet`, `walls`);
`meshDict.template` only sets their types.
**The mesher must be fed the FEATURE-EDGE surface `box5x2p5x2p5mm.fms`**
(`surfaceFeatureEdges box5x2p5x2p5mm.stl box5x2p5x2p5mm.fms -angle 45`, committed with the
case): the plain STL keeps the
four side walls in one solid, cfMesh's Voronoi dual then wraps faces around the four edges
between them, and 4.8 % of the wall faces end up tilted into the flow (up to 8 deg) -- a
uniform stream is no longer a discrete solution (`simpleFoam`: 8 % velocity error at the
corners, pressure +-0.2). With the `.fms`, 0 tilted faces and the stream is held to 1e-16
(measured 2026-09-05). Check any new polyhedral box the same way (face normals against the
patch plane) before trusting a translating case on it.
A `blockMeshDict.template` for a hexahedral twin is included.

Resolution is set at the INTERFACE: cfMesh's dual cells there measure 2^(-1/3) x maxCellSize
(measured on every polyhedral rung of the stationary droplet), so `MAX_CELL_SIZE =
h/0.7937` and `N_CELLS` (the capillary-dt handle) is pinned to `DOMAIN_LENGTH/h`; the pin is
verified on the built mesh with `check_refined_band.py --mode poly`. At N = 64 the interface
cell is h = 2.5e-3/64 = 3.90625e-05 m and `MAX_CELL_SIZE` is 4.921566601e-05 m.

**The capillary time step is pinned to the INTERFACE cell, not to the smallest cell.** The law
`dt = CAPILLARY_DT_COEFF/N_CELLS^1.5` reproduces 0.2323 of the Brackbill limit
`sqrt((rho1+rho2) h^3/(2 pi sigma))` at the interface spacing, and the band check confirms the
pin (measured: interface cells 0.991 h at maxCellSize 0.0195, dt/limit 0.236). cfMesh's dual
also makes much smaller cells along the box feature edges (0.29 h) and in the boundary slab
(0.25 h on the stationary box), where the same step is 1.5 to 1.9 times ABOVE the local limit.
That is harmless while the interface stays away from those cells -- the capillary limit is a
condition on cells that carry surface tension -- but any case whose interface can reach a wall
or an edge must pin the step to the smallest BAND cell instead.

    # 4-rank smoke (78 steps) on the coarsest mesh -- the gate before the cluster
    snakemake --workflow-profile profiles/local8 --configfile config/popinet3D_La12000_poly_smoke4.yaml \
        --config env_preamble="...; export LEIA_PMESH_PRELOAD=$HOME/miniconda3/lib/libjemalloc.so.2"
    # the ladder: R/h = 12.8 / 19.2 / 25.6 at the interface (N = 64 / 96 / 128; ~0.7M / 2.4M / 5.6M cells)
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r12p8.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r19p2.yaml
    snakemake --workflow-profile profiles/slurm --configfile config/popinet3D_La12000_poly_r25p6.yaml

Read-out as for the 2D reproduction (`make_popinet_table.py`): maximum over time of the L1
and L2 norms of |u - U0|/U0, plus volume, shape, Laplace jump and band curvature error at T.

### The amplification bound of the fit (`fitProbeDisplacement`, 2026-09-08)

The semi-Lagrangian update is LINEAR in the stencil values. With `g_j = b(d)^T M^-1 w_j^2
b(d_j)` for a departure displacement `d`,

    psi^{n+1}_c = fit_c(x_c + d) = (1 - sum_j g_j) psi_c + sum_j g_j psi_j,

so one step cannot amplify anything by more than the Lebesgue constant of the fit,

    Lambda_c = |1 - sum_j g_j| + sum_j |g_j|.

`Lambda_c = 1` at `d = 0`, and `Lambda_c = 1` whenever every weight is non-negative: the
update is then a convex combination of the stencil values and creates no new extremum, so no
mode can grow. `Lambda_c > 1` is the necessary condition for the checkerboard growth that
destroys the far field of a translating polyhedral case. Lambda depends on the stencil
GEOMETRY alone, so it costs one mesh pass and no time steps.

    // system/fvSolution, levelSet/semiLagrangian
    fitProbeDisplacement (-6.3995e-07 0 0);   // = -u dt; default (0 0 0) = diagnostic off

The reconstruction then writes the field `slFitAmplification` once, at construction, and
prints `max Lambda`. Recipe, in a COPY of a rendered case, serial, one step:

    foamDictionary -entry levelSet/semiLagrangian/writeFitOrder -set true system/fvSolution
    foamDictionary -entry levelSet/semiLagrangian/fitProbeDisplacement -set "(-6.3995e-07 0 0)" \
        system/fvSolution
    foamDictionary -entry endTime -set <deltaT> system/controlDict
    leiaSemiLagrangianLevelSetTwoPhaseFoam          # slFitAmplification lands in 0/

MEASURED on the SI Popinet meshes at `d = -U dt` (both N = 64, same dt, same h):

| mesh | class | cells | median | p99.9 | max | share > 1.10 |
|---|---|---|---|---|---|---|
| hex | band | 26 504 | 1.0164 | 1.0164 | 1.0164 | 0 |
| hex | far-interior | 497 784 | 1.0164 | 1.0430 | 1.0527 | 0 |
| poly | band | 27 027 | 1.0205 | 1.0205 | 1.0205 | 0 |
| poly | far-interior | 495 118 | 1.0205 | 1.0687 | 1.0963 | 0 |
| poly | far-small (cfMesh slab and edges) | 152 348 | 1.0329 | 1.2348 | 1.2608 | 0.39 % |

The bound is 5 times larger on the polyhedral mesh, and its whole excess sits in the small
one-sided cells: the five worst cells are all 0.32-0.33 h in size, and Lambda falls
monotonically with cell size (median excess 0.039 below 0.4 h against 0.021 at h). Those are
the cells where the sigma = 0 control produced its fake zero set.

`Lambda - 1` is PROPORTIONAL to the displacement, hence to the local Courant number: at
`|d| / U dt = 0.25 / 0.5 / 1 / 2 / 4` the maximum reads 1.0133 / 1.0265 / 1.0527 / 1.1041 /
1.2030 (hex) and 1.0684 / 1.1347 / 1.2608 / 1.4871 / 1.8413 (poly). Over a fixed physical
time the accumulated bound is `Lambda^(T/dt) = exp((Lambda - 1) T/dt) = exp(c U T)`, which
does not contain dt -- a smaller time step cannot remove the growth. `config/popinet3D_poly_
sigma0_dtSweep.yaml` tests that prediction directly.

**CONFIRMED 2026-09-09** (job 54503136, three arms, all COMPLETED). The sigma = 0 control ran
at dt, dt/2 and dt/4 on the SI polyhedral mesh at N = 64. It failed in every arm, at the steps
528 / 1009 / 1974 (ratio 1 : 1.91 : 3.74) and at the times 4.877e-03 / 4.660e-03 /
4.558e-03 s. The growth rate of `zeroSetRadialL2/R` per unit physical time is 132.78 / 128.59
/ 126.05 per second; the same rate per STEP is 0.1227 / 0.0594 / 0.0291 %, which halves and
quarters with the step. Before the failure the three trajectories are ONE function of time:
`zeroSetRadialL2/R` at matched physical times agrees to 0.3 % from t = 1e-4 s to t = 4.5e-3 s.
A smaller time step therefore cannot repair the polyhedral failure, and the near-wall cells do
not need their own step limit.

**Normalise before you threshold.** `zeroSetRadialL2`, `zeroSetRadialLinf`, `centroidError`
and the `m2*` columns are absolute lengths in metres. In SI they are 400 times smaller than in
Popinet's units, so a fixed threshold reads every SI run as clean. Divide by `DROPLET_RADIUS`
first, in every gate, curated table and comparison across unit systems.

Lambda > 1 does not PROVE instability: it is an upper bound, and the hexahedral mesh reaches
1.05 and runs the horizon. It proves the opposite, though: a mesh whose Lambda is 1
everywhere cannot grow a new extremum at all.

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

### Physical boundary faces in the reconstruction stencil (`SL_STENCIL_BOUNDARY_FACES`, 2026-09-05)

The two Popinet-3D polyhedral ladder rungs diverged late (steps 1151 / 1301) from a fake zero
level set that appeared near the OUTLET while the velocity was still quiet. Underneath: the
cell-to-cell stencils count every non-coupled, non-empty boundary FACE as a data point, and
with `psi zeroGradient` each face brings the cell's own value at h/2 (0.19 h in cfMesh's
boundary cells) at a 1/d^2 weight four to thirty times a neighbour's -- the fit's
boundary-normal gradient collapses toward zero. Measured on the two 78-step smokes at
t = 0.02 (measured / exact change of psi in the first cell layer): side walls 1.002 on both
meshes (tangential transport exact); inlet 0.082 (poly) / 0.198 (hex); outlet 0.240 / 0.491.
The polyhedral mesh then amplifies the outlet error in the size-transition cells behind its
boundary layer (e-fold ~100 steps) until psi changes sign there. Invisible on every
closed-box case and every kinematic gate (u = 0 on the walls).

    SL_STENCIL_BOUNDARY_FACES include   # the former behaviour (default until the exclude gates are recorded)
    SL_STENCIL_BOUNDARY_FACES exclude   # the fit uses CELL data only

`exclude` drops the boundary-face entries in `slReconstruction` for every reconstruction model
(the base-class accessors are the only path to the stencil); processor patches are remote
cells, never boundary slots, so MPI is untouched, and empty patches were never included. The
measurement is the transport fraction above (`psi_check`: PASS = ~1.00 at inlet and outlet on
both smokes with the droplet metrics unchanged to round-off), then the 2D Popinet horizon run,
then a long polyhedral run past the onset step before any rung is resubmitted.


### The value bound is a runtime-selectable family (`SL_VALUE_BOUND`, 2026-09-10)

The bound on the reconstructed value is a runtime-selected model, `slValueBound`, so `none`,
the falsified quasi-monotone clip and the distance-cone guard are ONE study axis instead of
three branches inside `slCorrector::robustEvaluate`. Types, which are also the dictionary
values:

| `valueBound` | what it does |
|---|---|
| `none` | no bound, only the runaway cap (stencil mid +- 10*stencil range) that predates the family. Reproduces `clipToStencilBounds false`. |
| `stencilBounds` | the quasi-monotone clip to the stencil [min, max], with `clipRegion` and `clipKeepExtrema`. Reproduces `clipToStencilBounds true`. FALSIFIED as a fix; kept so its arms stay reproducible. |
| `lipschitzCone` | the distance-cone bound, `psi_c^{n+1} = clip(H_c, l_c, u_c)` with `l_c = max_j(psi_j - L|x_d - x_j|)` and `u_c = min_j(psi_j + L|x_d - x_j|)`. |

**The default is a sentinel, not a type.** `valueBound fromClipSwitch` resolves to
`stencilBounds` when `clipToStencilBounds` is true and to `none` otherwise, so a case that
predates the family keeps its behaviour whether or not its template carries the token. The
pattern is the one `slopeLimiter` already used for `limitSlope`. MEASURED: all 8 arms of
`config/popinet2D_clipRegionGate.yaml` byte-identical over 1563 steps at np = 4 against
`studies/popinet2D_clipRegionGate_preBound_20260910`.

**Scope.** The bound acts in `robustEvaluate`, so it reaches the `pointValue` scheme only --
which is the production default. `fluxFormScheme` calls `slReconstruction::evaluate()`
directly and `normalProjectedScheme` never uses the corrector, so both stay unbounded.

**The cone-specific entries** (`lipschitzMode`, `lipschitzConstant`, `onInadmissible`) are
read by `lipschitzCone` only. That is a self-validating control: the four arms of each
non-cone bound in a sweep of those axes MUST be byte-identical to each other.

**Diagnostics**, written at write time next to `slClipEligible`/`slClipFired`/
`slClipFiredEver`: `slBoundDelta` (the change the bound made, in the units of psi, so the
interface damage is measurable), `slBoundSlack` (signed headroom over the displacement;
negative means the bound fired) and `slBoundInadmissible` (the empty-interval flag, which is
a direct measurement of how far psi has drifted from a distance function).

**Recipes.**

```bash
# exact-field unit gate: the review's 1D case, plane, sphere, and the cone apex
cd <a rendered case> && blockMesh && leiaSetFields && leiaTestSLReconstruction

# the NONLINEAR map's amplification. Lambda may NOT be used to certify a clip, and
# the power iteration needs linearity, so a bound is scored with -mode growth.
leiaTestTransportSpectrum -mode growth -seed checkerboard -amp 1 -nIter 2000

# the coupled 2D matrix
make studies-one STUDY=popinet2D_coneBoundGate
```

**MEASURED, and `lipschitzCone` is FALSIFIED as a transport bound** -- see METHOD.md 8.3 for
every table. It is exact at a distance cusp, and that result stands; it does not follow that
it transports well. On PURE ADVECTION, with every arm on the identical mesh, it is 2.3 to 3.6
times WORSE than no bound on uniform translation -- the one flow where its `L = 1` is exactly
valid -- 9 to 19 times worse in strained flow, and with `lipschitzMode stencil` it loses the
phase completely on 3D shear (`E_VOL_ALPHA_REL` = 1.0000). The mechanism is the review's own
bound `L(t) <= L(0) exp(int ||grad u||_inf dt)`, which this implementation does not carry: in
a strained flow the true Lipschitz constant grows and clamping to 1 destroys the field.

The coupled 2D hexahedral gain (every interface metric better by 25 to 68 % at N = 64) does
not survive refinement either: the eikonal error is 7.119e-03 at N = 64 and 7.066e-03 at
N = 128, an order of 0.01 -- a FLOOR, not a converging error -- and the centroid error
reverses to 119 % worse than no bound.

The unbounded POLYHEDRAL advection arm does diverge (floating-point exception at step 198)
and every bound prevents that, but `stencilBounds` prevents it with the best volume error of
the three. **The falsified monotone clip beats the cone bound on every advection gate.**

So `SL_VALUE_BOUND` stays at the sentinel that resolves to `none`. What survives is the
FAMILY -- the seat for the review's Rank 1 -- and the `-mode growth` instrument.
