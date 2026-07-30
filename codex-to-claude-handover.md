# Codex to Claude Code handover

Date: 2026-07-26

Repository: `/home/tmaric/OpenFOAM/repos/leia`

Branch: `feature/velocity-extension`

## Read this first

This is a heavily dirty worktree containing substantial pre-existing user work.
Do **not** reset, clean, discard, stage, or commit broad groups of files. Several
files touched in this investigation were already added or modified before the
rhoLENT work. Inspect targeted diffs and preserve unrelated changes.

There is currently no running `leiaSemiLagrangianLevelSetTwoPhaseFoam` process.
No commit was created.

## Research objective

The project is developing a cell-collocated, unstructured finite-volume level-set
method for surface-tension-driven two-phase flow in OpenFOAM. The intended
pipeline is:

1. semi-Lagrangian level-set transport on arbitrary polyhedral meshes;
2. interface-normal-constant velocity extension for transporting `psi`;
3. the Detrixhe--Aslam geometric cell phase indicator;
4. a geometrical face phase fraction/density;
5. mass/momentum-consistent large-density-ratio coupling; and
6. a conservative, pressure-compatible surface-tension force.

The unresolved scientific target is stable stationary, translating, and
oscillating capillary droplets without hiding instability behind short end times.

The consolidated roadmap is in
`docs/capillary-level-set-research-roadmap.md`. Read it before changing the force
formulation or declaring a gate passed.

## Attached rhoLENT paper and Gate 2 conclusion

The supplied paper was:

`C:\Users\TomislavMaric\Downloads\1-s2.0-S0021999123005211-main.pdf`

It is Liu et al., *An unstructured finite-volume level set / front tracking
method for two-phase flows with large density-ratios*, Journal of Computational
Physics 493 (2023) 112426.

The user's proposed algorithm--solve an additional conservative density equation
during pressure/velocity coupling, then reset density from the geometric
interface--is exactly rhoLENT Eq. 40 and Algorithm 1:

\[
\frac{\rho_c^{n+1}-\rho_c^n}{\Delta t}
+\frac{1}{V_c}\sum_f \rho_f^{n+1}F_f^o=0,
\qquad
\rho_f^{n+1}=\alpha_f^{n+1}\rho_1
+(1-\alpha_f^{n+1})\rho_2.
\]

The load-bearing part is not merely solving `rhoEqn`. On every PIMPLE outer
iteration, the **same** `rhoPhi = rhof*phi` must be used in both

```cpp
fvm::ddt(rho) + fvc::div(rhoPhi)
```

and

```cpp
fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
```

After the outer loop converges, current-time cell density is reset to

```cpp
rho = alpha1*rho1 + (1 - alpha1)*rho2;
```

so the next step's `rho.oldTime()` is geometrically consistent with the
Detrixhe--Aslam interface. Do not clip auxiliary density: clipping destroys the
discrete free-stream identity.

## Source implementation

### New transport setup

`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createTransportFields.H`

- Constructs the runtime-selectable `velocityExtension` after `U`, `phi`, `psi`,
  `alpha1`, and `rho` exist.
- Creates restartable `UextOld` with `READ_IF_PRESENT` and `AUTO_WRITE`.
- Supports mass-flux modes:
  - `interpolatedDensity`
  - `geometricFaceDensity`
  - `rhoLENT`
- Creates `alphaf` and `rhof`.
- Uses central owner/neighbour face-plane averaging by default for rhoLENT.
- Stores final-outer-iteration mass/reset diagnostics for the time-step CSV.

The mass-flux setup was removed from `createSLFields.H` because velocity
extension must be constructed only after the physical transport fields exist.

### Level-set and mass-flux time ordering

`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H`

On the first PIMPLE outer iteration only:

1. recompute `Uext` from current physical `U` and pre-advection `psi`;
2. advect with `slAdv->advect(psi, slVelExt->Uext(), UextOld)`;
3. save the current pre-momentum extension into `UextOld`;
4. update the narrow band, optional filter/redistancing, DA `alpha1`, and
   curvature;
5. construct geometrical `alphaf` and `rhof` at the new interface location.

On **every** PIMPLE outer iteration:

1. rebuild `rhoPhi` with the current pressure-corrected `phi`;
2. if `rhoLENT`, solve the auxiliary density equation;
3. use that same `rhoPhi` in the existing momentum equation.

The interface and `rhof` remain fixed within a time step; `phi^o` and therefore
`rhoPhi^o` change across outer iterations.

### Auxiliary density equation

`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H`

- Solves `fvm::ddt(rho) + fvc::div(rhoPhi) = 0` with the configured diagonal
  solver.
- Does not bound or clip `rho`.
- Logs absolute L1, Linf, and scale-relative L1 mass residuals.

### Post-coupling reset

`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/resetRhoLENT.H`

- Runs after the complete PIMPLE outer loop.
- Measures `mean |rhoAux-rhoGeom|`.
- Resets `rho` from DA `alpha1`.
- Reconstructs the diagnostic pressure `p = p_rgh + rho*gh` and maintains the
  reference-pressure correction.

### Solver wiring

`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/leiaSemiLagrangianLevelSetTwoPhaseFoam.C`

- Includes `velocityExtension.H`.
- Includes `createTransportFields.H` after the sibling `createFields.H`.
- Includes `resetRhoLENT.H` after PIMPLE convergence.
- Writes the expanded droplet diagnostics after the density reset.

### Direct normal trajectory velocities

`src/leiaLevelSet/semiLagrangian/pointValueScheme.H` and
`pointValueScheme.C` now support
`levelSet.semiLagrangian.trajectoryVelocity` values:

- `input`: the legacy supplied `U`/`Uext` path, unchanged by default;
- `normalProjection`: insert `(U dot n)n` directly into the Taylor displacement;
- `normalClosestPoint`: walk to the interpolated `psi=0` closest point, sample
  physical `U`, retain its normal component there, and use the local projection
  as a guarded fallback.

For both direct modes, the scheme projects the physical velocity field named by
`trajectoryU` (default `U`) rather than accidentally re-projecting an Eulerian
`Uext`. It constructs the current and old projected velocities from the
corresponding `psi` and `psi.oldTime()` normals and computes
`grad(projectedVelocity)` for the Taylor convective acceleration. Projecting
only the leading displacement while retaining `grad(U)` would be inconsistent.

The regularised normal is

```text
n = grad(psi)/sqrt(magSqr(grad(psi)) + epsilon^2).
```

The current closest-point implementation uses `interpolationCellPoint` for
`psi`, `grad(psi)`, and `U`; it is not the Detrixhe--Aslam plane geometry. A
processor-crossing walk currently falls back to the local projection. Treat
this mode as a diagnostic prototype until shared DA geometry and distributed
point sampling are implemented.

## Diagnostics added

`createDropletMetricsFile.H` and `writeDropletMetrics.H` now write:

- `maxMagU`, `meanMagU`;
- `pLaplace`;
- curvature L2/Linf error in the force-support band;
- minimum and maximum `|grad(psi)|` in that band;
- band L2 error of `|grad(psi)|-1`;
- band maximum/minimum gradient ratio;
- rhoLENT absolute L1 and Linf mass residuals;
- rhoLENT relative L1 residual;
- the auxiliary-to-geometrical density reset mismatch.

The standalone checker is
`workflow/scripts/check_capillary_gate.py`. It rejects incomplete/non-finite
CSVs and accepts explicit thresholds for end time, velocity, gradient quality,
mass residual, and pressure jump.

Example:

```bash
python3 workflow/scripts/check_capillary_gate.py CASE/leiaSemiLagrangianLevelSetTwoPhaseFoam.csv \
  --target-time 0.03 \
  --max-velocity 2 \
  --min-grad 0.5 \
  --max-grad 2 \
  --expected-pressure-jump 72.74 \
  --pressure-jump-tolerance 5
```

Old CSVs written before these columns were added will intentionally fail the
checker as incompatible evidence.

## Cases and workflow configuration

Updated translating and oscillating IST cases/templates:

- added a diagonal solver for `rho`;
- added a `Uext` linear solver;
- added `div(velExtPhiW,Uext) Gauss upwind`;
- added `velocityExtension` and `massFlux` dictionaries;
- set the current capillary default to `steadyUpwind` plus `rhoLENT`;
- set `nDefCorrExt 5`, `defCorrRelax 0.7`, and
  `zeroInflowGuard 1e-6`.

Relevant paths:

- `cases/transISTDroplet2D*`
- `cases/oscISTDroplet2D*`
- `config/transISTDroplet2D.yaml`
- `config/oscISTDroplet2D.yaml`

The controlled Gate-1/Gate-2 matrix is
`config/transISTTransportGate.yaml`. It expands to eight cases:

```text
velocity extension = none, meshWave, closestPoint, steadyUpwind
mass coupling       = geometricFaceDensity, rhoLENT
N                   = 32
end time            = 0.03 s
```

Dry-run command:

```bash
snakemake -s workflow/Snakefile \
  --configfile config/transISTTransportGate.yaml \
  -n --cores 1
```

The dry run succeeded and produced 8 `generate_case`, 8 `mesh`, and 8 `solve`
jobs. Failed capillary cases can collapse to extremely small time steps rather
than terminate promptly, so run the real matrix with per-case wall-clock limits
or a monitor.

The direct-trajectory matrix is `config/transISTTrajectoryGate.yaml`. It holds
Eulerian velocity extension at `none`, holds rhoLENT fixed, and compares
`input`, `normalProjection`, and `normalClosestPoint`. Its Snakemake dry run
successfully expands to three cases.

## Verification performed

### Build

```bash
source /home/tmaric/OpenFOAM/OpenFOAM-v2506/etc/bashrc
wmake applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam
```

This completed successfully with OpenFOAM v2506 and no compiler errors.

### Fresh smoke run

Temporary case:

`tmp/gates/transISTDiagnosticsSmoke`

The fresh `steadyUpwind + rhoLENT` run completed. During its stable part:

```text
absolute mean rho mass residual  ~ 1e-10
relative mean rho mass residual  ~ 2e-14 to 5e-14
min |grad(psi)| in band          ~ 0.94
max |grad(psi)| in band          ~ 1.00
```

The gate checker passed with explicit velocity, gradient, pressure-jump, and
mass-residual thresholds.

### Restart

The smoke case was restarted from its latest written time. `UextOld` was read
from disk, the solver completed, and the new CSV row remained finite. This
verifies the new restart path; it is not yet a bitwise uninterrupted/restarted
equivalence test.

## Controlled translating-droplet result

The completed A/B work so far covers the **translating droplet only**. The runs
used `N=32`, rhoLENT, the central geometric face fraction, quadratic SL
reconstruction, the structured integral surface-tension model, and a
capillary-limited step. The approximate onset/final runaway times were:

| Transport velocity | Approximate runaway time | Outcome |
|---|---:|---|
| physical velocity (`none`) | 0.0103 s | failed |
| `meshWave` | 0.0131 s | failed after a small delay |
| `closestPoint` | 0.0258 s | failed after a substantial delay |
| `steadyUpwind` | 0.0388 s | longest lived, then failed |

The exact onset definition is not automated yet; killed/stalled CSVs may end
with partial final rows. Treat these times as approximate nonlinear-runaway
windows, not high-precision benchmark values.

The scientific conclusion is:

> Velocity extension substantially delays the translating-droplet instability,
> but does not remove it. rhoLENT closes the mass/momentum consistency gap, but
> does not remove the capillary force/profile feedback.

Do **not** generalise this table to stationary or oscillating droplets. A
controlled extension-by-rhoLENT comparison has not yet been run for those two
cases. The statement currently supported by evidence is only that translating
droplet survival improved by about `1.3x`, `2.5x`, and `3.8x` for meshWave,
closestPoint, and steadyUpwind respectively.

### Direct Taylor-projection result

Both new modes passed fresh smoke runs. In the initial band,
`normalClosestPoint` sampled 36/36 cells successfully at both time levels with
zero fallback, and rhoLENT remained near `5e-14` relative residual.

Longer translating-droplet results with Eulerian extension disabled were:

| Taylor velocity | Outcome |
|---|---|
| `normalProjection` | computationally completed `t=0.03`, but was physically invalid there; failed later near `t=0.04365` |
| `normalClosestPoint` | failed near `t=0.01309` after fallback grew above 560 cells |

At `t=0.03`, `normalProjection` gave:

```text
maxMagU                 0.2478 m/s
pLaplace                68.10 Pa
kErrL2Band              659 1/m
min/max |grad(psi)|     0.690 / 1.819
rho relative residual   1.2e-14
```

Because `U0=0.05 m/s`, the total maximum implies
`max|U-U0| >= 0.198 m/s`; do not describe this as a clean or stable result.
It is less violent than the `steadyUpwind` field extension at the same time
(`maxMagU` about `1.15`, `pLaplace` about `51.2 Pa`, curvature error about
`1142 1/m`) and gives the longest translating numerical-survival window
measured so far.
It is still not a cure: by `t~0.04365`, velocity is `O(10^3) m/s` and
`max|grad(psi)|` is `O(10^2)`.

### Co-moving diagnostics and force-mechanism gate (2026-07-26)

`createDropletMetricsFile.H` and `writeDropletMetrics.H` now read
`levelSet.dropletReferenceVelocity` (the translating case sets
`(0.05 0 0)`) and write:

- `maxMagUPrime`, `meanMagUPrime`, and `disturbanceKineticEnergy` for
  `U'=U-Uref`;
- `phaseVolume`, `phaseVolumeRelError`;
- `centroidX/Y/Z` and `centroidError` relative to `center0+Uref*t`.

`workflow/scripts/check_capillary_gate.py` now supports
`--max-disturbance`, `--max-phase-volume-relative-error`, and
`--max-centroid-error`.

Three reproducible one-case configs hold `normalProjection + rhoLENT` fixed:

| Config | Force control | Result |
|---|---|---|
| `config/transISTFreeStreamGate.yaml` | `sigma=0` | reaches `t=0.05`; max `|U'|=2.33e-15 m/s` |
| `config/transISTConstantCurvatureGate.yaml` | exact `kappa=1000 1/m`, balanced `sigma*kappa*snGrad(alpha)` | reaches `t=0.05`; max `|U'|=1.49e-7 m/s` |
| `config/transISTIntegralForceGate.yaml` | current `integralSurfaceTension` | FPE near `t=0.0375`; `|U'|=5.65e-2 m/s` after the first step and max complete-row value `0.523 m/s` |

The first two controls have essentially identical interface histories. Both
reach a maximum DA phase-volume change of `6.04%` and centroid error of
`0.0853 mm` by `t=0.05`, while the band curvature L2 error grows to about
`3003 1/m`. Thus:

1. uniform-velocity preservation is excellent on this regular mesh;
2. exact-curvature CSF/pressure compatibility is excellent on this mesh;
3. the current integral curvature/tangent force injects the immediate
   spurious current and causes the capillary runaway;
4. independently, the SL/DA kinematics still fail phase-volume and centroid
   fidelity, even when `U` stays exact.

The integral case's auxiliary rho becomes negative during the final runaway.
Treat that as a downstream consequence/accelerant, not evidence that rhoLENT
initiated the instability: the co-moving current is already `O(1e-2)` in the
first step, while the zero/exact-force rhoLENT controls remain stable.

### Curvature-source and curvature-free integral-force split (2026-07-26)

The current structured `integralSurfaceTension` is not fundamentally a pure
curvature-free line-integral operator. Its crossing-tangent contribution is
geometrical, but the paper transcription also contains `sigma*kappa*s`
pressure-jump corrections. I added `curvatureSource field|constant|none` so
that term can be isolated without changing the remaining stencil.

With physical `input` trajectories, rhoLENT, `N=32`, and the capillary
timestep, the controlled results are:

| Force | First-step max `|U'|` | Result |
|---|---:|---|
| structured CST, transported `kappa` | `0.0564534 m/s` | attempted failure at `t=0.00293996`; max completed-row `3.542 m/s` |
| structured CST, exact `kappa=1000 1/m` | `0.0370659 m/s` | attempted failure at `t=0.00359995`; max completed-row `5.041 m/s` |
| structured CST, correction disabled | `0.336614 m/s` | worse; attempted failure at `t=0.00263996` |
| legacy balanced CSF, exact `kappa=1000 1/m` | roundoff scale | completes `t=0.05`; max `|U'|=1.43545e-7 m/s` |

Interpretation: exact curvature only modestly reduces the initial CST defect,
so transported curvature is not its root cause. Conversely,
`curvatureSource none` merely deletes a required correction from that
structured paper stencil; it is not the desired curvature-free integral
formulation and predictably performs worse.

I therefore added a separate runtime model,
`integralConormalSurfaceTension`, plus a default-zero
`cellSurfaceTensionForce()` interface on the force base class. It reconstructs
one owner-oriented conormal traction on every cut internal face and uses
`fvc::surfaceIntegrate` for exactly equal-and-opposite owner/neighbour cell
forces. `UEqn.H` accepts this optional vector cell force; legacy models keep
their existing face-normal path.

An important orientation bug was found during this experiment. Conormal
orientation used `(crossingPoint-ownerCellCentre)`, which is ambiguous at a
corner cut and produced an artificial `2*sigma/h^2` endpoint-force scale. All
three conormal helpers in `isoGeometricCurvature.H` now use the actual
owner-oriented `mesh.Sf()[face]`.

After that fix, the static `N=32` circle gives:

```text
max cell-force density                  2.22727e5 N/m^3
net force                               O(1e-20) N
conormal kappa L1/L2/Linf errors         364.0 / 390.7 / 640.7 1/m
raw identity-weighted residual L2       1.40518e4 N/m^3
raw relative residual L2                0.494514
raw residual Linf                       9.49977e4 N/m^3
```

The pressure-range diagnostic is in `YoungLaplaceEqn.H`. It compares the
initial cell force against the reconstructed gradient of an identity-weighted
Young--Laplace pressure solution, using a cell-volume-weighted L2 norm. The
large `49.45%` residual is evidence of a non-gradient component, but it is not
the definitive momentum-projection residual. The actual `pEqn.H` path first
multiplies the cell predictor by `rAU`, uses `rAUf` in the pressure Laplacian,
and applies OpenFOAM's face interpolation/reconstruction. At water/air density
ratio this weighting is material. Replace this initializer-only diagnostic with
the correct weighted projection before quoting a final pressure-compatibility
percentage. Neither projection residual nor `max|U|` replaces a direct
surface-tension-force accuracy test. Shared-face conservation alone still does
not establish a mimetic coupling.

The primary next metric must be independent of the CFD solver. For an analytic
circle in pseudo-2D, compute the exact cell-integrated force from the oriented
interface-arc endpoint tangents,

```text
F_exact,c = sigma * emptyThickness * (t_exit - t_entry)
F_num,c   = V_c * fSigma,c
```

and report vector L1/L2/Linf error and convergence. Also compare each shared
face traction with its analytic value before summing it into cells. This
separates interface reconstruction, conormal quadrature/orientation, and
deposition errors. Only after these errors converge should pressure projection
be tested; maximum cell-centred velocity is the final end-to-end acceptance
observable, not the force-error measurement.

As a diagnostic bridge, `isoCurvature` method `conormalCSF` projects the same
cell force onto the interface normal, infers a scalar curvature, and sends it
through the existing pressure-compatible
`sigma*interpolate(kappa)*snGrad(alpha)` path. With corrected sign it starts
with `kappa=[1021.43,1538.47] 1/m`, gives first-step
`max|U'|=0.0285425 m/s`, reaches a maximum completed-row disturbance of
`1.4188 m/s`, and fails near attempted `t=0.00455993`. The curvature later
becomes noisy and changes sign. Balanced delivery improves the initial impulse,
but this reconstructed curvature is not accurate or dynamically stable enough.

The exact-curvature balanced-CSF control was rerun after the force API and
solver changes. It still completes `t=0.05` with max
`|U'|=1.43545e-7 m/s`, final `|U'|=4.973e-8 m/s`, and maximum DA-volume change
`6.394%`; this is the regression guard for future operator work.

Implementation/configuration files added or changed for this split:

- `src/leiaLevelSet/surfaceTensionForce/integralSurfaceTension.{H,C}`;
- `src/leiaLevelSet/surfaceTensionForce/integralConormalSurfaceTension.{H,C}`;
- `src/leiaLevelSet/surfaceTensionForce/isoGeometricCurvature.H`;
- `src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.{H,C}`;
- `applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H`;
- `applications/solvers/leiaLevelSetTwoPhaseFoam/YoungLaplaceEqn.H`;
- `config/transISTIntegralCurvatureSourceGate.yaml`;
- `config/transISTConstantCurvatureInputGate.yaml`;
- `config/transISTIntegralNoCurvatureGate.yaml`;
- `config/transISTIntegralConormalGate.yaml`;
- `config/transISTIntegralConormalCSFGate.yaml`.

The library and two-phase solver build successfully. The generated study logs
contain the experimental runs, but failed Snakemake jobs can remove buffered
CSV files; use the values above and rerun the named configurations for fresh
tables.

### Zero-force `input` versus `normalProjection` ablation

`config/transISTKinematicGate.yaml` is a two-case Gate-1 experiment. It holds
`sigma=0`, rhoLENT, `VELOCITY_EXTENSION=none`, the mesh, reconstruction, and
time step fixed and changes only `TRAJECTORY_VELOCITY`:

| Taylor velocity | max `|U'|` | max DA volume change | max centroid error | final curvature L2 error |
|---|---:|---:|---:|---:|
| `input` | `2.24e-15 m/s` | `6.394%` | `0.08512 mm` | `3163 1/m` |
| `normalProjection` | `2.33e-15 m/s` | `6.039%` | `0.08529 mm` | `3003 1/m` |

Both reach `t=0.05`. At `N=32`, `h=0.3125 mm`, the droplet diameter is only
`6.4 h`, and it travels `8 h`; the approximately `0.085 mm` centroid error is
`0.27 h`, `8.5%` of the radius, or `3.4%` of the travelled distance. The 1%
phase-volume and `0.1 h` centroid limits initially applied to these runs were
too strict and have been withdrawn. Do not label either coarse run a Gate-1
failure from those cutoffs.

That near-equivalence is specific to the original capillary-limited
`Co=0.0096`. The completed Courant and refinement studies below supersede any
broader interpretation.

### Courant sweep and fixed-Courant refinement

The workflow now supports `TIME_STEP_CONTROL=advective`, deriving
`dt=TARGET_ADVECTIVE_CO*DOMAIN_LENGTH/(N_CELLS*REFERENCE_SPEED)`. Capillary mode
remains the default and still derives `dt~h^1.5`.

At `N=32`, physical `input` improves monotonically as `Co` increases from 0.01
to 0.5: signed volume change falls from `+6.36%` to `+3.24%`, and moving-circle
zero-set L2 error falls from `123.5` to `61.2 um`. This is accumulated
point-value-remap error.

Direct `normalProjection` has a different tradeoff. Its volume change grows
from `+5.97%` at `Co=0.01` to `+24.22%` at `Co=0.5`; zero-set L2 is best near
`Co=0.1` (`91.0 um`) and worsens to `143.6 um` at `Co=0.5`. Fewer remaps improve
its profile/curvature, but large projected-characteristic steps distort the
interface.

The fixed-`Co=0.5`, fixed-physical-time refinement is:

| trajectory | N=32 volume / zero L2 | N=64 volume / zero L2 | N=128 volume / zero L2 |
|---|---:|---:|---:|
| input | `3.236% / 61.2 um` | `0.582% / 10.4 um` | `0.099% / 2.81 um` |
| normalProjection | `24.22% / 143.6 um` | `8.890% / 63.9 um` | `5.056% / 29.6 um` |

For input, max-volume orders are `2.47,2.55` and zero-set L2 orders are
`2.55,1.89`: approximately second order. Projection gives volume orders
`1.45,0.81` and zero-set orders `1.17,1.11`: roughly first order and much less
accurate. Input curvature L2 only changes `367 -> 263 -> 238 1/m`; projected
curvature worsens `457 -> 550 -> 591 1/m`, so curvature remains a separate open
problem.

New reproducible configs:

- `config/transISTKinematicTimeStepGate.yaml` (`input` Courant sweep);
- `config/transISTKinematicProjectionTimeStepGate.yaml`;
- `config/transISTKinematicRefinementGate.yaml`.

The droplet CSV now also writes signed `phaseVolumeChange` and
`zeroSetRadialL2/Linf` from linearly interpolated sign crossings relative to the
analytically translated circle.

Directly removing tangential displacement delays the integral-force capillary
runaway, but the new kinematic gates show that this is not a generally superior
transport method: it is timestep-sensitive and only roughly first order in the
fixed-Courant translating-circle test. Keep it as a diagnostic. The on-demand
closest-point sampler is also not robust once the discrete profile develops
extra/ambiguous zero-set structure.

## Interpretation of the remaining instability

The failure loop is still

```text
curvature/tangent error
    -> non-gradient capillary residual
    -> spurious velocity
    -> changed zero set and level-set profile
    -> larger or more coherent curvature/tangent error
```

rhoLENT prevents a separate density-flux inconsistency from injecting momentum.
It is necessary for large-density-ratio credibility but cannot turn a vortical
surface-force residual into a pressure gradient.

Velocity extension reduces the rate at which off-interface velocity degrades
`psi`. In the translating case this postpones the feedback. It does not preserve
a sufficiently well-conditioned profile indefinitely. With `closestPoint`, the
number of fallback cells increased as the profile deteriorated. With
`steadyUpwind`, the run survived longest but ultimately reached severe curvature
error, gradient collapse, large velocity, and timestep stagnation.

Do not clip auxiliary `rho` when this occurs. Negative/overshooting auxiliary
density during full runaway is an important diagnostic and a downstream
accelerant; clipping would conceal the mass/momentum identity being tested.

## Major remaining technical issues

### 1. Stationary and oscillating controlled tests are missing

Run the same `none/meshWave/closestPoint/steadyUpwind` comparison for stationary
and oscillating droplets before making a general statement about velocity
extension. Use identical force, density coupling, resolution, time-step law, and
end-time criteria.

### 2. The surface-tension model is not unstructured

`src/leiaLevelSet/surfaceTensionForce/integralSurfaceTension.C` explicitly builds
structured `(i,j)` addressing and requires a uniform, single-block 2-D hex mesh.
It is not the final collocated unstructured formulation.

The existing `isoCurvatureSurfaceTension`/conormal variants are negative-result
prototypes. Do not present them as the solution. A valid arbitrary-polyhedral
operator needs one shared, pairwise-conservative internal-face traction whose
normal component is compatible with the same discrete pressure gradient used by
the collocated pressure equation.

### 3. Geometrical face fractions need a decomposed-mesh gate

`faceAreaFraction.H` produces one central owner/neighbour average on serial
internal faces, but boundary faces currently use the owner reconstruction. A
proper processor/cyclic implementation must exchange reconstructed face
geometry and produce the same `alphaf` on both sides. Add serial/decomposed
parity tests before claiming parallel rhoLENT support.

### 4. Central alpha-f is a symmetric approximation

The rhoLENT paper uses one geometric face fraction. The current central mode
averages face fractions obtained from the owner and neighbour cell planes. This
does create one numerical face value used identically by mass and momentum, but
it is an approximation to a unique shared interface reconstruction. Preserve
this distinction in publications and tests.

### 5. Dynamic meshes are not covered

The SL characteristic code assumes the static cached-stencil path. The current
implementation should not be advertised as a dynamic-mesh rhoLENT method.

## Recommended next sequence

1. Treat the three force controls above as the new diagnostic baseline; plot
   `maxMagUPrime`, disturbance energy, curvature/profile, phase volume, and
   centroid on the same time axis.
2. Treat physical `input` as the kinematic reference: its uniform-mesh
   fixed-Courant gate is approximately second order. Keep `normalProjection`
   diagnostic-only; its capillary survival benefit does not compensate for its
   first-order, timestep-sensitive translation error.
3. The curvature-source split is complete: exact `kappa` does not cure the
   structured CST, deleting its correction is worse, and a direct conservative
   conormal cell source has a `49.45%` *raw identity-weighted* initial residual.
   Preserve these as negative-result controls, but do not quote that number as
   the PIMPLE-weighted compatibility result.
4. Build the unstructured operator as a scalar owner-oriented face-normal force
   density/flux, not as a cell vector. It must enter `phig`, be weighted by the
   same `rAUf` as the pressure Laplacian, and be differenced against
   `p_rghEqn.flux()` before the shared `fvc::reconstruct` call.
5. Keep analytic traction/cell-force comparisons only as geometry diagnostics.
   Do not treat cell-source insertion plus `HbyA` interpolation as the final
   coupling; the decisive acceptance observable is maximum cell-centred
   velocity after the complete pressure solve.
6. Keep the exact-curvature balanced-CSF case as a mandatory regression and
   after direct force convergence require the new operator's weighted
   pressure-range residual to be small and convergent before running stationary
   or oscillating droplets. Use maximum cell-centred velocity only as the final
   coupled acceptance measurement.
7. Repeat free-stream/rhoLENT on skewed and polyhedral meshes; fix coupled-patch
   `alphaf` and verify serial/decomposed parity.
8. Only after the static operator gate passes, run stationary, translating, and
   oscillating refinement/density-ratio ladders. Keep `normalProjection`
   diagnostic-only.

## N=64 translating capillary gate and `rAU/rAUf` diagnostic (2026-07-26)

The new controlled matrix is
`config/transISTN64QuadraticCapillaryGate.yaml`. It uses `N=64`, quadratic SL,
rhoLENT, `Utrans=(0.05 0 0) m/s`, and physical `input` trajectories. This is a
laboratory-frame advection test: the solver advances a translating droplet and
the metric is computed cellwise as `max|U-Utrans|`. Do not describe it as a
static test or a moving-frame calculation.

| model | first `max|U-Utrans|` | largest completed value | outcome | first final-corrector `rAU/rAUf` residual |
|---|---:|---:|---|---:|
| exact constant-curvature balanced CSF | `1.09e-8 m/s` | `2.17e-7 m/s` | reaches `t=0.05` | `5.77e-9 m/s` |
| structured integral CST | `1.003e-2 m/s` | `4.63 m/s` | FPE attempting `t=0.001379` | `1.065e-2 m/s` |
| direct conormal cell source | `5.63e-1 m/s` | `3.99 m/s` | FPE attempting `t=0.001421` | `5.83e-1 m/s` |
| conormal-CSF bridge | `3.61e-2 m/s` | `9.03e-1 m/s` | FPE attempting `t=0.001506` | `3.91e-2 m/s` |

First-step volume-weighted relative-L2 `rAU/rAUf` residuals in the same order
are `3.70e-9`, `6.78e-3`, `6.18e-1`, and `1.74e-2`.

The exact control's maximum DA phase-volume relative error over `t=0.05` is
`1.91%` (`1.05%` at the end); its maximum translation-aware centroid error is
`11.3e-6 m` (`4.32e-6 m` at the end). Refining from `N=32` to
`N=64` improves the structured CST first impulse from approximately `5.65e-2`
to `1.00e-2 m/s`, so the old `N=32` result was too coarse quantitatively. It
does not change the scientific verdict: the current force still runs away in
about `1e-3 s`.

`applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H` now reports, once per
timestep on the final PIMPLE outer/inner correction, the velocity increment
left after applying the real cell `rAU`, face `rAUf`, pressure-equation flux,
and reconstruction. It is a coupled secondary diagnostic, not a replacement
for the analytic cell-integrated force error. The exact control being at the
pressure-solver tolerance scale while all three geometric models leave large
first-step residuals is strong evidence that their discrete forces are not in
the pressure-gradient range.

### Cell force averaged to the PIMPLE face-flux space (2026-07-27)

Implemented the requested bridge in
`integralConormalSurfaceTension.{H,C}`. Runtime entry
`surfaceTensionForce.delivery` accepts:

- `directCell`: original conservative cell source, retained as a negative
  control;
- `faceFluxFromCellForce`: compute the identical conservative conormal cell
  force, form `Gsigma = linearInterpolate(fCell) & Sf`, return
  `Gsigma/magSf` through the face-force API, and return a zero cell source.

There is no double counting. In the new mode the force is present in the
momentum predictor, `phig`, pressure equation and velocity recovery only through
the scalar face-flux path. The two-case reproducible matrix is
`config/transISTN64ConormalFaceFluxGate.yaml`.

| delivery | first `max|U-Utrans|` | largest completed value | first `rAU/rAUf` max / relative-L2 | outcome |
|---|---:|---:|---:|---|
| `directCell` | `0.5634097 m/s` | `3.98906 m/s` | `0.583247 / 0.618467` | FPE attempting `t=0.0014213` |
| `faceFluxFromCellForce` | `0.3204034 m/s` | `3.53849 m/s` | `0.317543 / 0.178854` | FPE attempting `t=0.0011243` |

The correct flux delivery reduces the first velocity impulse by about `43%`
and the first weighted relative-L2 residual by about `71%`, confirming that the
OpenFOAM force-flux path matters. It is not a cure and must not be promoted:
both runs have order-one cell velocities before failure, and the bridge fails
slightly earlier. Averaging a conservative cell force does not make its face
flux the discrete gradient of a scalar pressure jump. The next derivation must
construct that scalar face flux directly from the shared interface geometry.

## Herrmann JCP 2008 and unstructured PLIC conclusion

The detailed Stanford technical precursor to
[Herrmann's 2008 JCP paper](https://doi.org/10.1016/j.jcp.2007.11.002) was
reviewed together with the paper metadata/abstract. Preserve these lessons:

- Herrmann overlays a separately refinable equidistant Cartesian level-set
  geometry grid on either structured or unstructured flow meshes.
- The flow-cell volume fraction is integrated from the geometry-grid analytic
  phase fraction; curvature, treated as a surface quantity, is transferred
  separately.
- The face capillary force uses the same phase-fraction face difference as the
  pressure jump and participates explicitly in the Rhie--Chow-like
  face-velocity correction. In OpenFOAM terms, the crucial invariant is that
  the scalar capillary face flux and `p_rghEqn.flux()` are differenced before
  their common `fvc::reconstruct`; exact balance is facewise cancellation.
- Exact prescribed curvature gives machine-level equilibrium on unstructured
  prism/tetrahedron flow grids. Interface-projected and normally extended
  computed curvature gives approximately second-order convergence as the
  geometry grid is refined.
- The remaining volume-error floor is traced to non-divergence-preserving
  velocity interpolation from the flow grid to the geometry grid.

This says our closest-point concept is most defensible for evaluating curvature
at the interface and extending it normally, not for projecting away tangential
physical advection. It also says that a cell source plus generic `HbyA`
interpolation cannot be the final capillary coupling. We need a scalar
face-normal force flux in the same space and with the same `rAUf` weighting as
the pressure-Laplacian flux; only their difference is reconstructed.

The user's PLIC assessment is essentially right for **computed curvature on
general unstructured meshes**, but not as an impossibility result. Modern
unstructured PLIC/RDF methods can reconstruct the interface accurately, and
unstructured balanced VOF methods achieve solver-tolerance equilibrium with
prescribed exact curvature. What is still not broadly competitive with
structured height functions is the complete arbitrary-polyhedral chain from
PLIC through computed curvature to capillary equilibrium. Recent specialised
2-D triangular height-function work is an exception to watch, not yet a general
solution.

Acceptance must distinguish:

- completing a requested end time;
- retaining a conditioned level-set profile;
- converging pressure jump and parasitic currents under refinement; and
- genuine late-time stability.

## Repository documentation and Claude memory

The main local synthesis is:

`docs/capillary-level-set-research-roadmap.md`

It was assembled after reviewing the standalone articles/decks under `docs/`,
the solver/library implementation, and negative-result records.

There is no project-specific Claude Code memory for `leia`. The repository has
only `.claude/settings.local.json`. The directories under
`/home/tmaric/.claude/projects/*/memory` belong to unrelated projects and must
not be treated as leia research notes.

## Final state

- Solver build: passing.
- Fresh rhoLENT smoke run: passing.
- Restart smoke: passing.
- Gate checker: passing on the fresh smoke output.
- Eight-case Gate-1/Gate-2 workflow: dry-run passing, full matrix not run as one
  workflow campaign.
- Translating droplet: velocity extension delays but does not cure runaway.
- Direct normal projection has the latest translating numerical runaway so far,
  but is already physically invalid at `t=0.03` and fails near `t=0.04365`.
- The zero-force and exact-curvature controls reach `t=0.05` with co-moving
  disturbances below `2.4e-15` and `1.5e-7 m/s`, respectively.
- The integral-force control injects a first-step disturbance of `0.0565 m/s`
  and crashes near `t=0.0375`; it is the immediate capillary trigger.
- The structured CST curvature-source split is complete: exact curvature still
  fails, while deleting the paper's jump correction is substantially worse.
- A separate curvature-free shared-face conormal model now conserves global
  force to roundoff. Its raw identity-weighted projection gives a `49.45%`
  residual. The actual first-step final-corrector `rAU/rAUf` residual is
  `0.583 m/s`, versus `5.77e-9 m/s` on the first step of the exact balanced-CSF
  control (whose maximum over `t=0.05` is `3.94e-8 m/s`).
  Pressure-compatible delivery is the current Gate-3 blocker.
- The conormal-to-balanced-CSF bridge reduces the first disturbance to
  `0.0285 m/s` but remains dynamically unstable because its inferred curvature
  is noisy.
- Averaging the conormal cell force to faces and delivering
  `linearInterpolate(fCell)&Sf` through `phig` lowers the first N=64 disturbance
  from `0.563` to `0.320 m/s`, but still fails near `t=0.00112`; it is a useful
  bridge/negative result, not the final face-force construction.
- The exact-curvature balanced-CSF control still completes `t=0.05` after the
  solver/API changes, with max disturbance `1.435e-7 m/s`.
- At the original very small capillary Courant number, stable `N=32` controls
  show about 6% DA phase-volume change and 0.08 mm centroid lag. Fixed-`Co`
  refinement subsequently shows approximately second-order convergence for
  physical `input`; broader Gate-1 mesh/timestep coverage remains open.
- Fixed-Courant refinement is complete: physical `input` is approximately
  second order for volume and zero-set position; `normalProjection` is roughly
  first order and substantially less accurate.
- Direct closest-point normal sampling fails near `t=0.01309` as fallbacks
  proliferate.
- Stationary/oscillating velocity-extension verdict: intentionally deferred
  until the frozen-interface pressure-range gate passes.
- Arbitrary-polyhedral conservative surface tension: shared traction geometry
  and a cell-to-face-flux bridge now exist. The open derivation is a scalar
  capillary face flux constructed directly from that geometry and integrable by
  the `rAUf`-weighted pressure Laplacian.
- `N=64` translating gate: exact control completes; all three current
  geometrical-force models fail near `t=0.001`--`0.0015` despite the structured
  model's substantially smaller initial impulse than at `N=32`.

## Curvature-delivery replay completed (2026-07-27)

Implemented and ran the next controlled gate requested after the strict
`Gsigma` model matrix.

### Implementation

- `interfaceMeanCurvature.H` adds a diagnostic curvature extension named
  `interfaceMean`.  It computes the current reconstructed curvature, forms its
  `|snGrad(alpha)| |Sf|`-weighted face-support mean, and replaces the cell field
  by that one spatial value.  It still enters through
  `reconstructedCurvature::faceSurfaceTensionForceFlux()` and the complete
  PIMPLE path; it is not a solver bypass or a production curvature model.
- The translating-droplet template now parameterises `WRITE_INTERVAL`,
  `CURVATURE_EXTENSION`, `CURVATURE_EXTENSION_RELAX`,
  `CURVATURE_FACE_INTERPOLATION`, and `CURVATURE_FORCE_WEIGHT`.
- Reproducer: `config/transISTN64CurvatureDeliveryReplay.yaml` plus
  `workflow/scripts/run_curvature_delivery_replay.py`.
- Study output: `studies/transISTN64CurvatureDeliveryReplay`.
- Solver rebuilt successfully against OpenFOAM-v2506.  The exact baseline and
  all 28 one-step replay cases returned status 0.

### Experiment

The baseline uses N=64, uncached quadratic SL, no VE, physical `input`
trajectories, Detrixhe--Aslam, rhoLENT/central face density, and exact constant
curvature.  Snapshots near `t={0,0.01,0.025,0.05}` are therefore transported by
the real method without capillary-flow contamination.  Each replay copies only
the snapshot `psi` and `alpha.water`, resets U with the case's uniform
`Utrans=(0.05,0,0) m/s` field, freezes interface advection, and advances one
complete PIMPLE timestep.  The decision metric is only maximum cell-centred
`|U-Utrans|`.

| variant | t=0 | t~0.01 | t~0.025 | t~0.05 |
|---|---:|---:|---:|---:|
| exact constant | `7.60e-9` | `5.07e-9` | `9.50e-9` | `7.37e-9` |
| interface mean | `7.66e-9` | `5.11e-9` | `9.60e-9` | `3.84e-9` |
| quadratic cell centre | `7.11e-3` | `8.47e-3` | `1.08e-1` | `2.03` |
| closest-point Newton | `7.07e-3` | `8.43e-3` | `1.14e-1` | `1.25` |
| foot-point height function | `1.49e-4` | `2.98e-3` | `1.18e-1` | `1.87e-1` |
| quadratic Kang | `1.19e-3` | `2.06e-3` | `1.01e-1` | `1.08` |
| FVM Kang | `1.91e-3` | `2.89e-3` | `1.64e-2` | `2.50e-2` |

### Verdict for Claude

Spatial variation of computed curvature on the CSF support is now isolated as
the velocity-producing defect.  The mean-only control remains at pressure
solver tolerance at every snapshot even though the computed mean itself is not
physically accurate at late time.  Thus:

1. do not spend the next iteration changing PIMPLE, `rAUf`, or
   `fvc::reconstruct`;
2. do not promote `interfaceMean`--it suppresses real variable curvature and
   cannot predict oscillation dynamics or the correct pressure jump;
3. use the replay as the cheap mandatory development gate;
4. develop a shared zero-surface curvature representation with tangentially
   regularised local fits and normal-constant delivery;
5. only after the replay remains small across all four snapshots, return to the
   fully coupled translating and oscillating droplets.

## Connected zero-surface curvature candidate (2026-07-27, Codex)

Implemented the requested next step and deliberately stopped before a coupled
droplet run.

### Code path

- `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/connectedInterfaceCurvature.H`
  implements the current 2-D serial prototype.
- A point-interpolated `psi=0` contour supplies one globally unique crossing per
  shared internal face.  The two crossings of each cut cell form one interface
  element.  Do not revert to independent per-cell LSQ-plane exits: that first
  attempt fragmented a clean circle into 17 components.
- Every snapshot currently reconstructs 52 elements, one closed component,
  zero open endpoints, 52 successful ordered fits, and zero fallbacks.
- One curvature value is attached to each element after an ordered tangential
  quadratic fit, Helmholtz or biharmonic chain regularisation, and a planar
  Gauss--Bonnet additive-mode constraint.
- `createSLFields.H` registers `kappaInterfaceFace` before the force model is
  constructed.  Normal-ray delivery fills this surface field directly.
- `reconstructedCurvature` accepts
  `faceInterpolation connectedInterface` and returns
  `sigma*kappaInterfaceFace*snGrad(alpha)*magSf` as integrated `Gsigma`.
  There is no cell-to-face curvature interpolation on this path.
- The solver fills the connected field both initially and after every phase
  update.  `kappa` is also filled by the same normal-ray parent assignment for
  diagnostics, but is not the production force input.
- Current template defaults for this option are half-width 5, Helmholtz
  parameter 16, and 120 iterations.  They have no effect unless
  `curvatureExtension connectedInterface` is selected.

### Replay hardening and results

`workflow/scripts/run_curvature_delivery_replay.py` now parses the topology
audit line.  A connected candidate fails the runner if the reconstruction is
empty, has more than one component, has open ends, or uses a fallback.  This is
load-bearing: an early bug returned machine-zero velocity only because
`kappaInterfaceFace` was empty.

All five real candidates passed the topology/solver audit on all four N=64
snapshots.  Maximum cell-centred `|U-Utrans|` for the selected
`connectedInterface_h5_l16` candidate is:

| snapshot | max `|U-Utrans|` [m/s] |
|---|---:|
| `t=0` | `6.0116e-4` |
| `t~0.01` | `4.3806e-4` |
| `t~0.025` | `1.7906e-3` |
| `t~0.05` | `4.5485e-3` |

At `t~0.05`, compare FVM--Kang `2.4956e-2`, foot-point height function
`1.8697e-1`, quadratic Kang `1.0833`, and quadratic cell centre `2.0313` m/s.
Full provenance and all candidates are in
`studies/transISTN64CurvatureDeliveryReplay/curvature_delivery_replay.csv` and
the mirrored article table.

### Important limitation / exact next step

Do not launch the oscillating droplet yet.  Helmholtz lambda 16 is quiet but its
estimated transfer factor for the circle's physical `m=2` mode is only about
0.52, before attenuation by the half-width-5 fit.  The tested biharmonic
filters preserve low modes better but leave `3.50e-2` to `5.21e-2` m/s at the
last replay snapshot.

Next implement a manufactured perturbed-circle gate with known modes
`m=2,3,4`.  For each regularisation candidate, report the recovered curvature
mode amplitude/phase and keep the same four-snapshot maximum-velocity replay.
Admit a candidate to coupled translating/oscillating runs only if it satisfies
both gates.  Parallel and 3-D connected reconstruction remain unimplemented.

## 2026-07-28: shared curvature service, executed mode gate, constrained filter

The requested follow-up is implemented and executed.
This section supersedes the preceding provisional recommendation and its
half-width-5 replay selection; those entries are retained as experiment history.

- `surfaceTensionForce` now owns the registered face-curvature lookup and the
  integrated CSF assembly.  Set `faceCurvatureSource registered` and
  `faceCurvatureField kappaInterfaceFace` to make any scalar-curvature CSF model
  consume the shared face field.  Native estimators remain the default with
  `faceCurvatureSource model`.
- The registered path is wired through reconstructed, div-grad-psi, both trace
  variants, Kang, alpha-CSF, and iso-CSF.  Kang explicitly bypasses its cell
  interpolation on this path.  Integral conormal and exact constant curvature
  deliberately remain separate.
- The replay proves that reconstructed/div/trace/Kang/iso models are identical
  when curvature and geometric-alpha localisation are identical.  Alpha-CSF
  differs because it retains smoothed-alpha localisation and is worse at the
  late snapshot (`0.11393` versus `0.05209 m/s`).
- `applications/test/leiaTestConnectedCurvatureModes` calls the production
  connected implementation and exports its real elements, raw curvature and
  regularised curvature.  The Python driver runs 8 candidates, 3 modes, 3
  resolutions and 2 mesh families and publishes detailed/summary CSV files.
- Fixed gate for `N>=64`: transfer `[0.8,1.2]`, phase error <=10 degrees,
  relative L2 <=0.25, valid closed one-component topology.  N=32 is reported as
  a topology/under-resolution point.
- Unconstrained winner: half-width-3 biharmonic, 20 passes.  It preserves modes
  but its late replay is `0.05209 m/s`.
- New operator `helmholtzPreserveModes` strongly filters the graph and restores
  the arclength-Fourier projection of modes 2--4 before the Gauss--Bonnet mean
  correction.  Half-width 4, lambda 16, 120 iterations passes the mode gate
  with transfer `0.8190..1.1342`, max phase error `1.426 deg`, and max relative
  L2 `0.0390`.  Its late replay is `0.0128448 m/s`.
- Full N=64 translating config:
  `config/transISTN64ConnectedModePreserving.yaml`.  It completed `t=0.05` but
  failed dynamically: max/final `|U-Utrans|=0.08933 m/s`, maximum transient
  volume error `2.376%`, final radial Linf `0.2038 mm`.  The concise result is
  `docs/method-comparison/method-comparison-article/data/tables/transISTN64ConnectedModePreserving.csv`.
- The requested oscillating follow-up was run at N=64.  Before accepting its
  result, the legacy `implicitEllipsoid` initialiser was found to be algebraic,
  with initial band `|grad(psi)|~10^3`.  A separate
  `signedDistanceEllipse` runtime surface now computes the exact Euclidean
  closest-point distance without changing the legacy type.  The corrected
  case starts with band `|grad(psi)|=0.9877..0.99997`.
- The corrected oscillation still fails: `SIGFPE` at `t=0.020895 s` (requested
  `0.1 s`), maximum velocity `15.46 m/s`, mode-2 amplitude growth from
  `0.0949` to `0.5121 mm`, and `4.965%` maximum volume error.  Its first four
  crossings imply `T_2=8.801 ms` versus inviscid `9.508 ms`, but this is not a
  valid frequency result because the amplitude then grows catastrophically.
  Summary and trace are in
  `oscISTN64ConnectedModePreserving.csv` and
  `oscISTN64ConnectedModePreservingTrace.csv` in the method-comparison tables.

Next: add an exact variable-curvature velocity oracle.  Feed analytic ellipse
or perturbed-circle curvature to the same CSF faces, run the identical frozen
one-step PIMPLE solve, and compare cell-centred `max|U_candidate-U_oracle|`.
This preserves the real mode-2 restoring response while isolating curvature
delivery/localisation/coupling error.  Repeat on analytically translated
manufactured shapes before tuning the filter or attempting another long
coupled oscillation.

## 2026-07-28: exact variable-curvature velocity oracle executed

Provenance: after a repository-wide `./Allwclean`, the entire project rebuilt
successfully against OpenFOAM v2512.  The velocity-oracle matrix and the
manufactured curvature replay were then rerun with the v2512 executables and
reproduced the reported values bit-for-bit.

The follow-up oracle is implemented in the production connected-interface
path.  `curvatureExtension.estimator analyticImplicitSurface` evaluates the
configured surface's exact curvature at the same reconstructed elements and
uses the same normal-ray face assignment.  It deliberately bypasses the fit,
regulariser and Gauss--Bonnet correction; everything downstream is identical.

`workflow/scripts/run_variable_curvature_velocity_oracle.py` runs analytic and
candidate cases at N=32,64,128 on uniform and deterministic 10%-perturbed
meshes.  It parses the written internal U fields and reports the direct
cell-centred maximum velocity difference.  Do not substitute a force norm.
The perturbed rows were subsequently rerun with the velocity-converged eight
non-orthogonal pressure corrections and reproduced to the printed precision.

Uniform results `(oracle max U, candidate max U, max difference)` [m/s]:

- N32: `(0.02895, 0.04268, 0.01374)`;
- N64: `(0.03699, 0.04008, 0.004023)`;
- N128: `(0.008376, 0.009743, 0.001654)`.

Perturbed results:

- N32: `(2.364, 3.341, 2.230)`;
- N64: `(3.673, 3.648, 0.4643)`;
- N128: `(0.7119, 0.9000, 0.7118)`.

The analytic oracle itself is 82--99 times larger on perturbed than uniform
meshes.  Exact curvature is therefore not the controlling remaining error.
The true signed-distance input nevertheless reconstructs band minima
`|grad(psi)|=0.569,0.510,0.463` on the perturbed meshes, versus values tending
to one on uniform meshes.  The remaining shared suspects are interpolated
vertex psi, least-squares normals, Detrixhe--Aslam tangent planes,
`snGrad(alpha)`, density weighting and `rAUf`.

Next implement an analytic-geometry/localisation oracle: evaluate the analytic
surface directly at vertices, use analytic normals in the connected elements
and Detrixhe--Aslam planes, and rerun exactly this velocity matrix.  A collapse
of perturbed velocities identifies geometry/localisation; persistence isolates
pressure/`rAUf` coupling.

## 2026-07-28: analytic geometry/localisation oracle executed

The next separation is implemented.  Two new dictionary controls default to
the unchanged production path `levelSetField`:

- `curvatureExtension.geometrySource analyticImplicitSurface` evaluates the
  configured surface directly at vertices and uses analytic element normals;
- `phaseIndicator.geometrySource analyticImplicitSurface` supplies analytic
  signed-distance tangent planes to both Detrixhe--Aslam cell alpha and the
  rhoLENT face-area/density calculation.

The velocity runner now compares three paths on every mesh: connected fitted
curvature, exact curvature with numerical geometry, and exact curvature with
analytic geometry/localisation.  All 18 frozen one-step v2512 solves completed
with valid one-component closed topology and no fallback.

Maximum cell velocities `(exact curvature, + analytic geometry)` [m/s]:

- uniform N32/N64/N128: `(0.02895,0.03858)`, `(0.03699,0.02847)`,
  `(0.008376,0.006917)`;
- perturbed N32/N64/N128: `(2.364,3.921)`, `(3.673,3.263)`,
  `(0.7119,1.129)`.

The analytic-geometry perturbed/uniform ratios are `101.6`, `114.6` and
`163.2`.  There is no collapse.  Analytic vertices, normals and tangent-plane
offsets are therefore not the controlling cure.  The remaining numerical path
is `snGrad(alpha)` -> density/`rAUf` weighting -> pressure Laplacian -> common
`fvc::reconstruct`.

Next implement a constant-curvature pressure-potential gate on uniform and
perturbed N=32,64,128 meshes.  Compare the current integrated CSF flux
`sigma*kappa*snGrad(alpha)*magSf` with
`snGrad(sigma*kappa*alpha)*magSf`, route both through the same PIMPLE pressure
and velocity recovery, and measure only maximum cell-centred velocity.  Do not
return to translating or oscillating droplets until this constant-curvature
skew-mesh gate is quiet.

## 2026-07-28: frozen-circle pressure-compatibility gate executed

The requested gate is implemented in
`workflow/scripts/run_pressure_compatibility_gate.py`, configured by
`config/staticISTPressureCompatibilityGate.yaml`, and reproducible with
`make pressure-compatibility-gate`.  The oscillating-droplet template now has
parameterised in-plane axes; the gate sets both to 1 mm, keeps zero initial
velocity, freezes the interface, and runs one capillary step at N=32,64,128 on
uniform and deterministic 10%-perturbed meshes.

The new run-time model `constantCurvaturePressurePotential` constructs
`pSigma=sigma*kappa*alpha` as a cell field and returns the strict integrated
oriented scalar flux `snGrad(pSigma)*magSf`.  Its comparison is the existing
`constantCurvatureSurfaceTension` flux
`sigma*kappa*snGrad(alpha)*magSf`.  Both use identical analytic
Detrixhe--Aslam localisation, rhoLENT-central, pressure controls, `rAUf`,
PIMPLE and `fvc::reconstruct`.  The runner reads written U fields; no force
norm is a decision metric.

Maximum cell-centred velocities `(CSF product, pressure potential)` [m/s]:

- uniform N32/N64/N128: `(3.770e-9,3.770e-9)`,
  `(7.635e-9,7.635e-9)`, `(1.039e-8,1.039e-8)`;
- perturbed N32/N64/N128: `(8.584e-6,8.584e-6)`,
  `(7.328e-4,7.328e-4)`, `(1.392e-3,1.392e-3)`.

The maximum direct differences between the two written velocity fields are
`2.02e-14`, `4.25e-14` and `1.50e-13 m/s` on the perturbed meshes (and below
`8.4e-16 m/s` on uniform meshes).  Thus the constant-kappa product and
pressure-potential assemblies are equivalent to roundoff.  Localisation
product assembly is not the controlling error.  The potential path itself is
not quiet on skew meshes, and the error grows under refinement.

## 2026-07-28: perturbed-mesh pressure correction audited and enforced

Important provenance correction: the coupled droplet tests did **not** use
perturbed meshes all along. `stationaryDroplet2D`,
`transISTN64ConnectedModePreserving`, `oscISTN64ConnectedModePreserving` and
the legacy oscillating-droplet config all specify `mesh: hex`.  The
oscillating interface was a perturbed ellipse on a uniform mesh.  Explicitly
perturbed meshes entered in the curvature-mode gate (no pressure solve) and
the later variable-curvature and pressure-compatibility oracle matrices.
Those CFD oracle rows initially inherited `nNonOrthogonalCorrectors 1`.

The reproducible sweep is
`workflow/scripts/run_nonorthogonal_correction_sweep.py` / make target
`pressure-nonorthogonal-sweep`.  On perturbed exact-circle pressure-potential
cases, maximum cell velocities [m/s] for corrections `(0,1,2,4,8,16)` are:

- N32: `(1.358e-5,1.006e-5,8.536e-6,8.583e-6,8.584e-6,8.584e-6)`;
- N64: `(1.370e-3,7.519e-4,7.378e-4,7.328e-4,7.328e-4,7.328e-4)`;
- N128: `(1.183e-2,1.146e-3,1.341e-3,1.388e-3,1.392e-3,1.392e-3)`.

Values at 8,16,32,64 are converged to the printed precision.  Therefore one
correction was inadequate and eight is now the minimum for perturbed CFD
gates.  This correction changes the pressure-gate numbers, but does not make
them quiet; the converged N128 velocity is actually larger than the
one-correction value.

Enforcement is now explicit:

- `oscISTDroplet2D` and `transISTDroplet2D` parameterise
  `N_NON_ORTHOGONAL_CORRECTORS`;
- `materialize.py` raises it to at least 8 whenever `mesh == perturbed`;
- both custom perturbed CFD oracle runners set 8 and abort below 8.

That planned operator and `rAUf` work is now complete.  The pressure studies
have also moved to the canonical dedicated DAG:

```bash
snakemake -s workflow/Snakefile.pressure-compatibility \
  --workflow-profile profiles/local
```

Make targets with `pressure-*` names are thin aliases only.

Paired corrected/uncorrected maximum velocities [m/s] at N=32,64,128 are
`(8.584e-6,1.044e-6)`, `(7.328e-4,7.537e-4)` and
`(1.392e-3,1.815e-3)`.  The uncorrected pair is quieter only at N=32 and is
30% worse at N=128.  It does not cure the gate.

Replacing face-varying `rAUf` by its volume-weighted mean also fails as a cure.
For the corrected operator it changes N=32,64,128 from
`(8.584e-6,7.328e-4,1.392e-3)` to
`(3.163e-6,8.582e-4,1.603e-3)`.  It helps only the coarse grid and worsens both
resolved grids.  The cell `rAU` used in velocity recovery stays physical in
this diagnostic.

The decisive new result is the pressure linear solver.  PCG/DIC with
`tolerance 1e-11`, `relTol 0` gives
`(9.050e-6,4.029e-5,4.839e-5) m/s`, versus production GAMG
`(8.584e-6,7.328e-4,1.392e-3) m/s`.  This removes the refinement-growing GAMG
artefact by 18.2x at N=64 and 28.8x at N=128.  Tightening PCG to `1e-13` does
not reduce the remaining plateau.  Strict GAMG is not reliable here: it makes
N=64 worse and raises SIGFPE at N=32 and N=128.

Exact next step: use PCG/DIC as the perturbed-mesh reference, recheck its
non-orthogonal convergence, and add a manufactured pressure-initialization
variant that sets `p_rgh = sigma*kappa*alpha` (modulo the reference constant)
instead of solving `YoungLaplaceEqn.H`.  Continue using only maximum written
cell velocity.  Collapse means startup pressure initialization was the
remaining source; persistence means isolate the common `fvc::reconstruct`
projection.  Do not resume variable-curvature, translating or oscillating
droplets before this exact-circle gate is quiet.
