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
  `periodic2D`); sets `U` and the flux `phi`.
- **velocityExtension** (`src/leiaLevelSet/velocityExtension/`) — optional velocity
  *correction* that, in a narrow band, replaces the advecting velocity with one constant
  along the interface normal (an extension velocity). Non-invasive: it reads `U`, emits a
  corrected field `Uext` and an advection flux, and leaves `U`/`phi` untouched.
  Types: `none` (default, identity), `anisotropicDiffusion`, `pseudoTime`. See
  [Velocity extension](#velocity-extension).

Solvers (`applications/solvers/`): **`leiaLevelSetFoam`** (advection verification) and
**`leiaLevelSetTwoPhaseFoam`** (two-phase flow with surface tension).

## Authors

* **Tomislav Maric** [ORCID](https://orcid.org/0000-0001-8970-1185)
* **Julian Reitzel** [ORCID](https://orcid.org/0000-0002-3787-0283)

## Installation

### Dependencies

- A C++17 compiler (tested with GCC 10.x / 11.x).
- **OpenFOAM** — tested with OpenFOAM-v2206 and OpenFOAM-v2506/v2512. Source it before
  building (`source .../OpenFOAM-vXXXX/etc/bashrc`).
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

### Legacy single-case runs

The committed per-case scripts still work for manual runs:

```
case> ./Allrun_hex_serial.sh        # also _hex_parallel / _perturbed_* / _poly_*
```

## Velocity extension

To advect `psi` with an interface-normal-constant velocity in a narrow band, select a
`velocityExtension` model in `system/fvSolution` (defaults to `none`):

```
levelSet
{
    velocityExtension
    {
        type        anisotropicDiffusion;   // none | anisotropicDiffusion | pseudoTime
        narrowBand  NarrowBand;  nLayers 3; // seed band + dilation
        epsilon     1e-3;                    // anisotropicDiffusion regularization
        nIterations 10;  deltaTau 0.3;       // pseudoTime
        projectFlux false;                   // true => correctFlux(phiExt) (volume cons.)
    }
}
```

The model writes a corrected velocity `Uext` (visualizable, written with `U`) and advects
the level set with the corresponding flux, without modifying `U`. `none` reproduces the
unmodified solver exactly.

## License

GPL-3.0 — see [LICENSE.md](LICENSE.md).

## Contributing

Fork the project and submit a merge request.
