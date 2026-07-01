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
profiles/local/config.yaml  executor: local   (mpirun)
profiles/slurm/config.yaml  executor: slurm    (one sbatch per case; srun + module env)
studies/<study>/            generated cases + <study>_database.csv  (git-ignored)
```

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
