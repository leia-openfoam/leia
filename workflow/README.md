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
config/staticExtension.yaml  static t=0 extension verification: e=|n.grad(Uext)| vs h
                             (solver: leiaTestVelocityExtension; 7 models x 2 div schemes)
config/steadyVortex2D.yaml   NON-REVERSING stress test: steady vortex (oscillation off,
                             T=3), 6 models; error-vs-time crossover t*; shape error vs a
                             marker-traced reference (scripts/marker_ref.py) since there
                             is no analytic final interface without reversal
config/steadyVortex2DHighRes.yaml  + N=256 tier (opt-in; refreshes the deck's steady_* figures)
config/bulkVortexSL.yaml     SEMI-LAGRANGIAN solver (solver: leiaSemiLagrangeLevelSetFoam)
                             on the reversed vortex; sweeps SL_RECONSTRUCTION
                             (linearTaylor/nestedLSQ/quadraticWLSQ) x CFL{0.5,1.0} -- the SL
                             analog of the VELOCITY_EXTENSION model sweep. plots.py emits the
                             reconstruction-convergence figure + an sl_vs_extension cross-study
                             overlay (reads the bulkVortexHighRes velocity-extension study).
config/bulkVortexSLHighRes.yaml  + N=256 tier (opt-in; refreshes the deck's sl_* figures)
profiles/local/config.yaml  executor: local   (mpirun; jobs x np = 24 ranks, no oversubscription)
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
`env_preamble` in `config/config.yaml` (defaults to OpenFOAM-v2506 — adjust to
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
- `mesh` ∈ {hex, perturbed, poly}; `poly` is 3D-only and uses
  `cases/<Case>_poly.parameter`. `perturbed` adds `-fluxCorrection`.
- `np` is the single source of truth: it regenerates `system/decomposeParDict`,
  drives `mpirun -np {np}` / SLURM `--ntasks`, so the rank count and the
  decomposition can never disagree.
- `3Dtranslation` runs its `system/init_End/init_End.<mode>` automatically.
- Non-parametric utility tests (1DredistanceTest, 2DgradTest):
  `snakemake --workflow-profile profiles/local utilities`
  (requires the `leiaTestRedistance` / `leiaTestGradScheme` apps to be built).
- On SLURM, edit the account and the module/`source` lines in
  `profiles/slurm/config.yaml` for your cluster.
