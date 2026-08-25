# leia

The **leia** project is an [OpenFOAM](https://develop.openfoam.com/Development/openfoam)
module implementing unstructured Lagrangian / Eulerian Interface Approximation (LEIA)
methods for multiphase flow in complex geometries — currently an **unstructured
Finite-Volume Level Set Method** for interface advection and two-phase flow.

[![Build Tests](https://github.com/leia-openfoam/leia/actions/workflows/build.yml/badge.svg)](https://github.com/leia-openfoam/leia/actions/workflows/build.yml)
[![Documentation](https://github.com/leia-openfoam/leia/actions/workflows/docs.yml/badge.svg)](https://leia-openfoam.github.io/leia/)

## Method

The transported field is the level set `psi`, advected by

```
d(psi)/dt + div(phi, psi) = source
```

The Heaviside / volume-fraction field `alpha` is recovered from `psi` **geometrically**:
in each interface (narrow-band) cell a linear least-squares (LLS) plane is reconstructed
from the cell-centre `psi` values, and the cell volume fraction on the negative side of
that plane is computed exactly. This avoids a first-order algebraic Heaviside and does not
require `psi` to be a signed-distance field; if `psi` is second-order accurate, so is the
recovered interface.

The level-set components are **runtime-selectable** (chosen in `system/fvSolution`):

- **phaseIndicator** (`src/leiaLevelSet/phaseIndicator/`) — `alpha` from `psi`:
  - `geometric` — LLS plane + polygon-clipping cell intersection;
  - `detrixheAslam` — same LLS plane, but the per-cell volume fraction is computed
    analytically by decomposing the cell into tetrahedra (cell centroid + face centroids +
    corner points) and applying the closed-form simplex formulas of Detrixhe & Aslam
    (2016). Tolerance-free and robust in 2D. Agrees with `geometric` to machine precision
    on the same reconstruction;
  - `heaviside` — first-order algebraic (smoothed) indicator.
- **narrowBand**, **redistancer**, **sdplsSource**, **profile** — interface-band marking,
  reinitialization, signed-distance-preserving source, and initial-field profiles.
- **velocityModel** (`src/leiaLevelSet/velocityModel/`) — prescribed *verification*
  velocity (`rotation`, `deformation3D`, `shear3D`, `translation`, `vortex2D`,
  `periodic2D`, `uniaxialStrain`); sets `U` and the flux `phi`. All are solenoidal
  except `uniaxialStrain`, `v(x) = a (x - x0).e_x e_x`, for which `div v = a != 0`
  intrinsically (a 1D flow has no tangential direction in which to compensate the
  normal strain); the psi equation is assembled in advective form
  (`- fvm::Sp(fvc::div(phi), psi)`), so this is handled correctly, but such a case
  must not be run with the solver option `-fluxCorrection`, which projects `phi`
  onto `div(phi) = 0`.
- **velocityExtension** (`src/leiaLevelSet/velocityExtension/`) — optional velocity
  *correction* that, in a narrow band, replaces the advecting velocity with one constant
  along the interface normal (an extension velocity). Non-invasive: it reads `U`, emits a
  corrected field `Uext` and an advection flux, and leaves `U`/`phi` untouched.
  Types: `none` (default, identity), `anisotropicDiffusion`, `pseudoTime`,
  `steadyUpwind`, `steadyUpwindLinear`, `closestPoint` (the statically
  ~O(h^2)-convergent geometric reference), `meshWave` (parallel-robust wave). See
  [Velocity extension](#velocity-extension).

Solvers (`applications/solvers/`):

- **`leiaLevelSetFoam`** — Eulerian advection verification;
- **`leiaLevelSetTwoPhaseFoam`** — Eulerian two-phase flow with surface tension;
- **`leiaSemiLagrangeLevelSetFoam`** — semi-Lagrangian advection verification (transports
  `psi` along characteristics, `psi^{n+1}(x_c) = psi^n(x_d)`, with a value least-squares
  reconstruction at the departure foot — no flux, no linear solve);
- **`leiaSemiLagrangianLevelSetTwoPhaseFoam`** — semi-Lagrangian two-phase flow: the same
  reinitialisation-free transport coupled to one-field incompressible Navier–Stokes, with
  the interface curvature evaluated **symbolically** from the quadratic reconstruction
  (gradient + Hessian) and extended constant along the normal for the continuum
  surface-tension force.

## Authors

* **Tomislav Maric** [ORCID](https://orcid.org/0000-0001-8970-1185)
* **Julian Reitzel** [ORCID](https://orcid.org/0000-0002-3787-0283)

## Installation

### Dependencies

- A C++17 compiler (tested with GCC 10.x / 11.x).
- **OpenFOAM-v2512** (current standard version; also builds against v2206/v2506).
  Source it before building (`source .../OpenFOAM-v2512/etc/bashrc`). On the
  laptop (WSL) it lives in `$HOME/OpenFOAM/OpenFOAM-v2512`; on Lichtenberg it is
  the source build in `$HOME/OpenFOAM/OpenFOAM-v2512` — see [CLUSTER.md](CLUSTER.md).
- [cfmesh](https://cfmesh.com/cfmesh/) (OpenFOAM sub-module) — polyhedral (`pMesh`) cases.

### Build

```
leia> git submodule update --init        # pyFoamStudy (legacy), cfMesh
leia> ./Allwmake
```

This builds the libraries (`libleiaLevelSet`, …) and the solvers/utilities
(`leiaLevelSetFoam`, `leiaLevelSetTwoPhaseFoam`, `leiaSetFields`, `leiaPerturbMesh`, …)
into `$FOAM_USER_*BIN`. Doxygen docs: <https://leia-openfoam.github.io/leia/>.

## Running the verification test suite (Snakemake)

**One-command reproduction.** The repo-root `Makefile` rebuilds everything from source —
the library and solvers, every convergence/verification study, and the figures, tables,
slide decks and article PDFs they feed (source OpenFOAM first, and have `snakemake`,
`foamlib`, `vtk`, `numpy`, `matplotlib` and `latexmk` on `PATH`):

```
leia> make build        # ./Allwmake: library + all solvers
leia> make studies      # run every study (SL advection, two-phase droplet, velocity-ext)
leia> make docs         # rebuild the reveal decks + compile the articles from the data
leia> make all          # build + studies + docs
```

`make studies` runs three groups: `studies-sl` (semi-Lagrangian convergence),
`studies-droplet` (the two-phase stationary droplet), and `studies-ve` (velocity
extension). Each study's results agglomerate into its theme's single data source
(`docs/<theme>/<slug>-article/data/{figures,tables}`), which both the deck and the
article consume, so slides and paper always reflect the newest numbers.

The parameter-study test suite is driven by **Snakemake** (replacing the old
`pyFoamStudy` flow). The same workflow runs locally (`mpirun`) or on SLURM (one job per
case) — you switch only the *profile*. See [workflow/README.md](workflow/README.md) for
details.

```
# one-time
pip install --user --break-system-packages "snakemake>=8" snakemake-executor-plugin-slurm

# choose the study in config/config.yaml (case, mesh, mode, np, scope), then:
leia> snakemake --workflow-profile profiles/local      # this machine (mpirun)
leia> snakemake --workflow-profile profiles/slurm       # cluster (one sbatch per case)
leia> snakemake --workflow-profile profiles/local -n    # dry-run: preview job count
```

Each study materializes one case per parameter combination from the `cases/<case>`
template + `cases/<case>.parameter` grid (`@!TOKEN!@` substitution), runs
mesh → (decompose) → solve, and aggregates the per-case CSVs into
`studies/<study>/<study>_database.csv`. Every case sweeps
`PHASE_INDICATOR ( geometric detrixheAslam )`, so a study directly compares the two
geometric indicators (a `PHASE_INDICATOR` column appears in the database). The default
`config/config.yaml` is a small smoke study; set `axes_override: {}` and
`collapse_other_axes: false` for the full grid (preview the size with `-n` first).

### Fast vs. high-resolution convergence studies

Adding one resolution level to a time-reversal benchmark costs **8× in 2D**
(4× cells × 2× CFL time steps), so the reversed-vortex convergence study is split
into a fast daily suite and a deliberate deep-convergence run:

```
# FAST default (N = 32/64/128, 6 models, 72 cases): ~5 min measured.
# Smoke/regression; keeps figures in studies/, never touches the slide deck
# (export_slides: false).
leia> snakemake --workflow-profile profiles/local --configfile config/bulkVortex.yaml

# HIGH-RES opt-in (adds N = 256; 6 models, 96 cases): ~46 min measured, floored by the
# heaviest cases (T=8, N=256 ~ 12k steps). THE deck-regeneration run: rebuilds the
# triptychs, the alpha / |grad psi|-1 field atlas and BOTH standalone deck variants:
# docs/velocity-extension/velocity-extension-presentation/velocity-extension.html
# (vertical section stacks) and its -linear variant (flat front-to-back flow).
leia> snakemake --workflow-profile profiles/local --configfile config/bulkVortexHighRes.yaml

# Phase-indicator comparison (geometric vs detrixheAslam), N <= 128 by design:
# the indicators coincide to ~1e-11 and psi dynamics duplicate the bulk `none` cases.
leia> snakemake --workflow-profile profiles/local --configfile config/phaseIndicatorConvergence.yaml

# STATIC (t=0) extension verification: apply each velocityExtension model ONCE to the
# initial configuration (no advection) and measure e = |n.grad(Uext)| convergence with h
# (leiaTestVelocityExtension; 7 models x 4 resolutions x 2 divergence schemes = 56 cases,
# ~1 min). The gatekeeper: what does not converge statically cannot converge advected.
leia> snakemake --workflow-profile profiles/local --configfile config/staticExtension.yaml

# NON-REVERSING stress test: STEADY vortex (oscillation off, T=3, 6 models). The
# reversed benchmark's final-time reading credits `none` with error cancellation
# (measured ~1300x at T=2/N=256); here the strain never reverses, the headline is
# error-vs-TIME per model (crossover t*), and the shape error is measured against a
# marker-traced quasi-exact reference (workflow/scripts/marker_ref.py) because no
# analytic final interface exists. Fast tier N <= 128 (minutes); the HighRes config
# adds N=256 and refreshes the deck's steady_* figures.
leia> snakemake --workflow-profile profiles/local --configfile config/steadyVortex2D.yaml
leia> snakemake --workflow-profile profiles/local --configfile config/steadyVortex2DHighRes.yaml

# SEMI-LAGRANGIAN advection (leiaSemiLagrangeLevelSetFoam): the parallel research
# line to velocity extension. Instead of correcting the velocity it solves
# Dpsi/Dt = 0 along characteristics (psi^{n+1}(x_c) = psi^n(x_d)); the reconstruction
# of psi^n at the departure foot is the swept axis (linearTaylor O(h) baseline,
# quadraticTaylor / quadraticWeightedLeastSquares O(h^2)). Same reversed 2Dvortex study; quadraticWeightedLeastSquares
# beats every velocity-extension model on shape AND |grad psi| with no linear solve.
leia> snakemake --workflow-profile profiles/local --configfile config/bulkVortexSL.yaml
leia> snakemake --workflow-profile profiles/local --configfile config/bulkVortexSLHighRes.yaml

# TWO-PHASE stationary droplet (leiaSemiLagrangianLevelSetTwoPhaseFoam): the
# parasitic-current benchmark of Lippert et al. (2022). A water/air drop (R=1 mm,
# sigma/R = 72.74 Pa, g=0) is at rest; the metrics are the spurious currents max|U|,
# mean|U| and the measured Laplace jump, written every step to the solver CSV. Mesh
# ladder N = 32/64/128/256; the report step fits the log-log order of max|U| vs h and
# regenerates data/tables/droplet_parasitic.tex + data/figures/sl_droplet_parasitic.png.
# Surface tension is explicit, so the step is capillary-limited: maxDeltaT is a DERIVED
# token, CAPILLARY_DT_COEFF / N^1.5 (~dt_sigma/4). 256^2 is the long pole (~33k serial
# steps, hours); use the slurm profile or reduce the ladder for a quick check.
leia> snakemake --workflow-profile profiles/local --configfile config/stationaryDroplet2D.yaml
# or, standalone at the base resolution:  cases/stationaryDroplet2D> ./Allrun.sh

# Fine-mesh DIVERGENCE MECHANISM probe (why the currents grow at N>=128): reruns the
# droplet at N=64 (stable) and N=128 (diverging) with the solver logging, per step,
# min|grad psi| and the band curvature error alongside max|U|, plus an N=128 control
# with the interface FROZEN (SL_FREEZE_INTERFACE=1, advection skipped). Shows the
# signed-distance loss LEADS the current runaway, and that freezing psi removes it --
# a transport-coupled feedback, not a force imbalance. Writes the compact CSVs under
# docs/.../data/mechanism/ and (via make_droplet_mechanism.py) sl_droplet_mechanism.png.
leia> bash workflow/scripts/droplet_mechanism_probe.sh   # WSL, OpenFOAM sourced
```

The two vortex configs write to **separate study directories**
(`studies/bulkVortexConvergence` vs `studies/bulkVortexHighRes`) — never point
configs with different axes at the same study dir: the cartesian-product case
indices would remap and stale cases would be silently reused. Per-model linear
solvers for the extension solves live in
`cases/2Dvortex/system/fvSolution(.template)` (`solvers → "Uext.*"`, selected via
the `VELOCITY_EXTENSION` token; discretization schemes in `fvSchemes(.template)`
via the `UEXT_DIV` aliasDict token — keep it `UpwindUext`: the linearUpwind
deferred correction diverges on the steady extension equation). Need more levels
than N=256? Use `--workflow-profile profiles/slurm` on a cluster.

### Legacy single-case runs

The committed per-case scripts still work for manual runs:

```
case> ./Allrun_hex_serial.sh        # also _hex_parallel / _perturbed_* / _poly_*
```

## Velocity extension

To advect `psi` with an interface-normal-constant velocity near the interface, select a
`velocityExtension` model in `system/fvSolution` (defaults to `none`):

```
levelSet
{
    velocityExtension
    {
        // none | anisotropicDiffusion | pseudoTime | steadyUpwind
        // | steadyUpwindLinear | closestPoint | meshWave
        // Statically verified (staticExtension study): closestPoint converges at
        // ~O(h^2) (the geometric reference); meshWave is the parallel-robust O(h)
        // wave fallback; steadyUpwind is the one-solve PDE transport (keep the
        // upwind scheme: the linearUpwind deferred correction diverges on the
        // steady equation); the legacy PDE models do not converge.
        type        closestPoint;
        narrowBand  NarrowBand;  nLayers 3; // |psi|-based band + fade length scale
        epsilon     1.0;  epsilonModel cellSize; // anisotropicDiffusion: eps_c = c*h
        nIterations 10;  deltaTau 0.3;       // pseudoTime
        nDefCorrExt 5;                       // steadyUpwind family
        projectFlux false;                   // keep false: projection destroys (n.grad)Uext = 0
        cpHaloReach 1.5;                     // closestPoint MPI halo: Taylor validity radius (x cellSize)
    }
}
```

The model writes a corrected velocity `Uext` (visualizable, written with `U`) and advects
the level set with the corresponding flux, without modifying `U`. `none` reproduces the
unmodified solver exactly.

**closestPoint in MPI-parallel runs:** the Newton foot-point walk is processor-local; a
`zoneDistribute` cell-point-cell halo (as in geometricVoF/plicRDF) carries one layer of
remote `(C, psi, grad psi, U, grad U)` around the band, and a walk that steps off the
local mesh continues on first-order Taylor data anchored at the nearest halo cell
(sub-stepped at ~one cell per step so the halo anchor tracks the true exit cell). Foot
points deeper than the one-layer halo fall back to the pinned `steadyUpwind` fill.
Measured on the steady vortex (N=64, np=4, scotch): fallback share 5.3% -> **1.3%**
(serial floor 0.15%); serial results are bit-identical (the halo branch never runs);
the parallel solution moves 2-3.5x closer to the serial one on the field metrics. The
per-step counters print as `closestPoint: N foot-pointed band cells (M via halo), K
fallback cells`.

**Coupling design (validated on the reversed vortex):** the psi equation is solved in its
**advective** form, `ddt(psi) + div(phiExt, psi) - psi*div(phiExt)` — the level set obeys
`Dpsi/Dt = 0`, and a normal-constant extension velocity is generically *not* solenoidal
(`div(v0*n) = v0*kappa`), so the conservative form alone hides a spurious compression
source. The extension flux is **seamless**: `phiExt = interpolate(Uext) & Sf` on all
internal faces, with `Uext` fading into `U` by a C1 cosine ramp over `|psi| in [L, 2L]`,
`L = nLayers*h` (a hard band-edge flux switch is a velocity discontinuity whose
`|grad psi|` kinks return onto the interface under flow reversal). Do **not** re-enable
`projectFlux`: Helmholtz-projecting `phiExt` to divergence-free destroys the
`(n.grad)Uext = 0` property that preserves `|grad psi| = 1` near the interface — with
projection the `|grad psi| = 1` channel at the interface measurably disappears.

## License

GPL-3.0 — see [LICENSE.md](LICENSE.md).

## Contributing

Fork the project and submit a merge request.
