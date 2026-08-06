# Unstructured capillary level-set research roadmap

## Objective

Develop a collocated, unstructured finite-volume level-set method in OpenFOAM
that remains stable for surface-tension-driven flow at large density ratio. The
target method has four coupled parts:

1. cell-centred level-set transport on arbitrary polyhedral meshes;
2. a geometric Detrixhe--Aslam phase indicator from one reconstructed interface;
3. one geometrical face mass flux used identically by mass and momentum; and
4. a conservative, pressure-compatible surface-tension force that does not
   amplify harmless level-set-profile errors.

This document records the dependency order and falsifiable gates. Passing an
earlier gate is necessary but not evidence that a later one passes.

## Evidence hierarchy

The repository contains six method lines rather than one finished method:

- the quadratic semi-Lagrangian (SL) reconstruction is the best measured pure
  advection method once the interface is resolved;
- the `quadraticTaylor` SL reconstruction (formerly named `nestedLSQ`) builds the
  same quadratic expansion from twice-differentiated `psi` instead of from stencil
  values, and needs the stencil clip to stay bounded; no genuinely linear
  reconstruction measured here is both stable and convergent for
  reinitialisation-free SL advection;
- SDPLS controls steepening without moving the zero set, but does not repair
  under-resolution flattening;
- geometric redistancing has useful local constructions, but the full method
  can turn remote transport sign errors into false interfaces;
- interface-normal velocity extension is a transport operator, not a physical
  momentum velocity;
- the 2-D structured integral surface-tension prototype reaches a much better
  static equilibrium, but still has a dynamic profile-feedback channel and is
  not an arbitrary-polyhedral formulation.

The article and negative-result decks are the authority for measured claims.
The generated HTML files duplicate their templates and are not independent
evidence. There is no project-specific Claude Code memory directory for this
repository: the repository has only `.claude/settings.local.json`, while the
memory directories under `~/.claude/projects` belong to unrelated projects.

## Current discrete pipeline

```text
physical U^n, psi^n
        |
        +--> normal-constant velocity extension --> Uext^n
        |                                            |
        |                         semi-Lagrangian characteristics
        |                                            |
        |                                         psi^(n+1)
        |                                            |
        |                shared least-squares interface reconstruction
        |                         /                  \
        |             Detrixhe--Aslam alpha_c      geometric alpha_f
        |                         |                  |
        |                  geometric rho_c     geometric rho_f
        |                                            |
        +-------------------- PIMPLE outer iteration o -------------------+
                              rhoPhi^o = rho_f phi^o
                              ddt(rho_aux) + div(rhoPhi^o) = 0
                              momentum convection = div(rhoPhi^o,U)
                              pressure correction updates phi^o
        +-----------------------------------------------------------------+
                              |
                 rho <- alpha_c rho_1 + (1-alpha_c) rho_2
                    for output and the next time level
```

The physical `U` remains the momentum unknown. `Uext` is used only by the
characteristic departure-foot calculation.

## Gate 2 decision: auxiliary density plus geometric reset

Yes. This is the algorithm in Liu et al., *Journal of Computational Physics*
493 (2023) 112426, Eq. 40 and Algorithm 1:

\[
 \frac{\rho_c^{n+1}-\rho_c^n}{\Delta t}
 + \frac{1}{V_c}\sum_f \rho_f^{n+1} F_f^o = 0,
 \qquad
 \rho_f^{n+1}=\alpha_f^{n+1}\rho_1
 +(1-\alpha_f^{n+1})\rho_2.
\]

The auxiliary `rho` is solved on every outer pressure--velocity iteration with
the same `rhoPhi = rhof*phi` used in `div(rhoPhi,U)`. After PIMPLE convergence,
cell density is reset from the Detrixhe--Aslam cell fraction. Consequently,
`rho.oldTime()` at the next step is geometrically consistent with the interface,
while the current auxiliary density supplies the discrete mass--momentum
identity during coupling.

Load-bearing requirements are:

- `rhof` and `alphaf` are fixed at the newly advected interface during a time
  step, but `rhoPhi` is rebuilt when the PIMPLE flux changes;
- mass and momentum use the identical face field and volumetric flux;
- auxiliary density is not clipped, because clipping breaks the free-stream
  cancellation;
- the geometric reset occurs after all outer correctors, never before them;
- an SL restart reads the stored previous-time extension `UextOld`.

The implementation is in `createTransportFields.H`, `slAlphaEqn.H`,
`rhoLENTEqn.H`, and `resetRhoLENT.H` under
`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam`.

## Gate status

| Gate | Question | Acceptance criterion | Status (2026-07-26) |
|---|---|---|---|
| 0 | Is the experiment reproducible? | Template/config parity, clean build, finite CSV, restart equivalence | Infrastructure passes; repository worktree is intentionally not normalised by this work |
| 1 | Does interface transport remain usable? | bounded `psi`; convergent zero-set error; band `min/max(|grad psi|)` remain conditioned through the physical test window | Uniform zero-force `input` provisionally passes at fixed Co with about second-order volume/zero-set convergence; capillary-timestep, skewed, and polyhedral gates remain open; direct projection is inferior |
| 2 | Are mass and momentum discretely consistent? | same `rhof*phi` in auxiliary mass and momentum; relative mass residual near roundoff; post-coupling DA reset | Implemented and passing before nonlinear runaway |
| 3 | Is capillarity pressure compatible on arbitrary polyhedra? | constant-curvature force lies in the range of the discrete pressure gradient; pairwise conservative internal-face traction; no structured addressing | Structured constant-curvature control passes; arbitrary-polyhedral shared-traction operator remains open |
| 4 | Does the coupled static droplet converge? | bounded/decaying parasitic current, convergent Laplace jump, no late-time regrowth over a capillary-time-scaled window | Open |
| 5 | Does the method survive dynamic capillary flow? | translating and oscillating droplets on hex, perturbed, and polyhedral meshes; mesh refinement and density-ratio ladders | Open |

## New controlled Gate-1/Gate-2 result

The `N=32` translating-droplet experiments used the same integral force,
quadratic SL reconstruction, central geometric face fraction, rhoLENT coupling,
and capillary-limited step. Only the velocity extension changed. The time below
is the approximate onset of timestep collapse or runaway, not a passing end
time:

| SL transport velocity | Approximate runaway time | Interpretation |
|---|---:|---|
| physical velocity (`none`) | 0.0103 s | rhoLENT alone does not remove the trigger |
| `meshWave` | 0.0131 s | small delay |
| `closestPoint` | 0.0258 s | substantial delay; fallback cells proliferate as the profile degrades |
| `steadyUpwind` | 0.0388 s | longest-lived case; still fails with severe profile/curvature degradation |

A subsequent direct-trajectory experiment removed the Eulerian `Uext` field
and inserted a projected velocity consistently into both Taylor time levels and
`grad(Utrajectory)`:

| Taylor trajectory velocity | Approximate runaway time | Interpretation |
|---|---:|---|
| `(U dot n)n` (`normalProjection`) | 0.04365 s | computationally reaches 0.03, but is already physically invalid there; still gives the latest numerical runaway measured so far |
| closest-point normal speed (`normalClosestPoint`) | 0.01309 s | on-demand projection loses its valid support as the profile fragments; no improvement |

At `t=0.03`, `normalProjection` had `max|U|=0.248`, so
`max|U-U0| >= 0.198 m/s` for `U0=0.05 m/s`; this is not a stable translating
drop. The pressure jump was
`68.10 Pa`, curvature-band L2 error `659 1/m`, and band
`|grad(psi)|` in `[0.690,1.819]`. The corresponding `steadyUpwind` field
extension had already reached about `max|U|=1.15`, pressure jump `51.2 Pa`,
curvature error `1142 1/m`, and subsequently failed near `t=0.0388`.
The local projection therefore postpones numerical runaway, but does not extend
the physically credible interval nearly that far. The longer run becomes
unstable near `t=0.04365`, with `max|U|=O(10^3)` and
`max|grad(psi)|=O(10^2)`.

## Translating-droplet mechanism gate

The decisive follow-up held `N=32`, `normalProjection`, rhoLENT, the quadratic
SL reconstruction, and the time step fixed, while changing only capillary-force
delivery. Diagnostics use the co-moving velocity `U' = U-(0.05,0,0)` rather
than total `|U|`.

| Control at `t=0.05` | Outcome | `max |U'|` | DA volume / centroid result |
|---|---|---:|---|
| `sigma=0` | completes; uniform velocity preserved | `2.33e-15 m/s` | max relative volume change `6.04%`; max centroid error `0.0853 mm` |
| exact `kappa=1/R`, `sigma*kappa*snGrad(alpha)` | completes; discrete pressure-force balance stable | `1.49e-7 m/s` | identical kinematic drift within roundoff |
| current `integralSurfaceTension` | FPE during pressure solve near `t=0.0375` | already `5.65e-2 m/s` after the first step; max complete-row value `0.523 m/s` | auxiliary rho becomes negative only after the capillary/current runaway is established |

This separates two defects. First, the direct SL trajectory plus rhoLENT
preserves the uniform translating velocity, and the constant-curvature CSF lies
in the range of the discrete pressure gradient on this mesh. The immediate
spurious-current source is therefore the current integral force fed by the
transported curvature/tangent geometry, not the Taylor trajectory or PIMPLE
coupling. Second, even with no disturbance velocity, the transported level set
causes the Detrixhe--Aslam phase volume to grow by about 6% and the centroid to
lag by about `0.26 h`. Gate 1 remains open independently of capillarity.

The reproducible cases are:

- `config/transISTFreeStreamGate.yaml`;
- `config/transISTConstantCurvatureGate.yaml`;
- `config/transISTIntegralForceGate.yaml`.

The CSV now records `maxMagUPrime`, `meanMagUPrime`, disturbance kinetic energy,
phase volume/change, and centroid/error. The checker accepts
`--max-disturbance`, `--max-phase-volume-relative-error`, and
`--max-centroid-error`.

### Curvature-source split and curvature-free force prototype

The word "integral" must not be conflated with "curvature-free" in the current
structured CST implementation. Its tangent/crossing stress is geometrical, but
the transcription of Abu-Al-Saud, Popinet & Tchelepi also contains
`sigma*kappa*s` pressure-jump corrections at interface crossings. A runtime
`curvatureSource field|constant|none` split now isolates that term without
changing the rest of the stencil.

The translating-droplet ablation uses physical `input` trajectories, rhoLENT,
`N=32`, and the capillary timestep:

| Force construction | First-step max `|U'|` | Outcome |
|---|---:|---|
| structured CST, transported `kappa` | `5.645e-2 m/s` | fails near attempted `t=0.002940`; max completed-row `3.54 m/s` |
| structured CST, exact `kappa=1000 1/m` | `3.707e-2 m/s` | fails near attempted `t=0.003600`; max completed-row `5.04 m/s` |
| structured CST, jump correction removed | `3.366e-1 m/s` | worse; fails near attempted `t=0.002640` |
| balanced CSF, exact `kappa=1000 1/m` | at roundoff scale | completes `t=0.05`; max `|U'|=1.435e-7 m/s` |

Exact curvature reduces the first CST impulse only modestly and does not cure
the instability. Removing the jump correction is not a valid curvature-free
integral method: it destroys part of the pressure-jump balance while retaining
the structured stencil. Curvature estimation is therefore not the sole defect.

A separate `integralConormalSurfaceTension` prototype now constructs one
owner-oriented vector traction per cut internal face from the reconstructed
interface segment and deposits equal-and-opposite cell forces with
`fvc::surfaceIntegrate`. It contains no curvature. The face orientation must use
the actual owner-oriented `mesh.Sf()[face]`; the former crossing-to-cell-centre
test was ambiguous for corner cuts and produced an artificial
`2*sigma/h^2` force scale. With that correction, the `N=32` circle gives:

- maximum cell-force density `2.227e5 N/m^3`, of the expected
  `sigma*kappa/h` scale;
- global resultant force at roundoff (`O(1e-20) N`);
- conormal-inferred curvature errors
  `L1/L2/Linf = 364/391/641 1/m`;
- raw identity-weighted pressure-range residual
  `||f_sigma-Gp||_2/||f_sigma||_2 = 0.4945`, with
  `Linf=9.50e4 N/m^3`.

Thus pairwise and global force conservation are necessary but insufficient.
The raw projection leaves a large non-gradient component, so inserting this
cell source in `UEqn` is not automatically a mimetic capillary operator.
However, `0.4945` is not yet the definitive PIMPLE compatibility ratio: the
initializer solves an identity-weighted Laplacian and uses a volume-weighted
cell norm, whereas the actual momentum projection applies cell `rAU`, face
`rAUf`, and the corresponding interpolation/reconstruction operators. This is
especially material at water/air density ratio. The solver-weighted residual
must be measured inside `pEqn.H` from the corrected velocity (or its exactly
equivalent discrete force residual).

As a diagnostic bridge, projecting the same conormal cell force onto the level-
set normal to infer a scalar curvature and then using the existing balanced
`sigma*kappa_f*snGrad(alpha)` delivery reduces the first disturbance to
`2.854e-2 m/s`. It still fails near attempted `t=0.004560`: the inferred
curvature starts in `[1021,1538] 1/m`, then becomes noisy and eventually changes
sign. This confirms that pressure-compatible delivery helps, while the present
conormal geometry is not yet accurate or stable enough.

Reproducible configurations are:

- `config/transISTIntegralCurvatureSourceGate.yaml`;
- `config/transISTConstantCurvatureInputGate.yaml`;
- `config/transISTIntegralNoCurvatureGate.yaml`;
- `config/transISTIntegralConormalGate.yaml`;
- `config/transISTIntegralConormalCSFGate.yaml`.

### Zero-force trajectory ablation

`config/transISTKinematicGate.yaml` repeats the zero-force case with only the
Taylor velocity changed. Both cases preserve `U0` to roundoff and reach
`t=0.05`:

| Taylor velocity | max `|U'|` | max DA volume change | max centroid error | final curvature L2 error |
|---|---:|---:|---:|---:|
| full physical `U` (`input`) | `2.24e-15 m/s` | `6.394%` | `0.08512 mm` | `3163 1/m` |
| `(U dot n)n` (`normalProjection`) | `2.33e-15 m/s` | `6.039%` | `0.08529 mm` | `3003 1/m` |

This is a deliberately coarse diagnostic: `h=0.3125 mm`, so the `2 mm`
diameter spans only `6.4` cells and the drop travels `2.5 mm = 8 h`. The maximum
centroid error is about `0.27 h` (`8.5%` of the radius and `3.4%` of the travel
distance). The previously applied 1% volume and `0.1 h` centroid limits are
withdrawn as unjustified at this resolution. These `N=32` results are baselines,
not a Gate-1 rejection.

That comparison applies only at its capillary-limited `Co=0.0096`. Subsequent
Courant and refinement gates show that the apparent equivalence does not persist
at practical advective steps.

### Zero-force Courant sweep

`config/transISTKinematicTimeStepGate.yaml` and
`config/transISTKinematicProjectionTimeStepGate.yaml` use
`dt=Co*h/U0` at `N=32`. This removes the irrelevant capillary timestep from the
`sigma=0` problem. All cases preserve `U0` to roundoff.

| `Co` | input: signed volume change | projection: signed volume change | input: zero-set L2 | projection: zero-set L2 |
|---:|---:|---:|---:|---:|
| 0.01 | `+6.36%` | `+5.97%` | `123.5 um` | `119.5 um` |
| 0.05 | `+5.53%` | `+7.02%` | `95.7 um` | `102.2 um` |
| 0.10 | `+4.84%` | `+8.64%` | `79.0 um` | `91.0 um` |
| 0.25 | `+3.95%` | `+14.47%` | `70.7 um` | `99.4 um` |
| 0.50 | `+3.24%` | `+24.22%` | `61.2 um` | `143.6 um` |

For physical `input`, fewer point-value remaps improve every reported interface
metric. For `normalProjection`, fewer remaps improve profile and curvature at
first, but the projected characteristic produces an increasingly large
step-size-dependent volume/centroid error; its zero-set error is best near
`Co=0.1` and then rises. Direct projection is therefore not simply an equivalent
normal-speed rewrite at the discrete characteristic level.

### Fixed-Courant refinement result

`config/transISTKinematicRefinementGate.yaml` holds `Co=0.5`, displacement, and
physical end time fixed:

| trajectory | `N` | max `|Delta V|/V0` | max centroid error | final zero-set L2 | final curvature L2 |
|---|---:|---:|---:|---:|---:|
| input | 32 | `3.236%` | `73.2 um` | `61.2 um` | `367 1/m` |
| input | 64 | `0.582%` | `6.33 um` | `10.4 um` | `263 1/m` |
| input | 128 | `0.099%` | `3.23 um` | `2.81 um` | `238 1/m` |
| normal projection | 32 | `24.22%` | `105.8 um` | `143.6 um` | `457 1/m` |
| normal projection | 64 | `8.890%` | `8.69 um` | `63.9 um` | `550 1/m` |
| normal projection | 128 | `5.056%` | `5.88 um` | `29.6 um` | `591 1/m` |

Physical `input` has max-volume orders `2.47, 2.55` and zero-set L2 orders
`2.55, 1.89`: the zero-force kinematics converge approximately second order on
this uniform mesh when `Co` is held fixed. `normalProjection` has volume orders
`1.45, 0.81` and zero-set orders `1.17, 1.11`; it is roughly first order here
and substantially less accurate. Centroid order is not monotone enough to quote
asymptotic accuracy. Neither curvature sequence is satisfactory: input plateaus
near `238 1/m`, while projected-trajectory curvature error grows with refinement.

Gate 1 is therefore provisionally passed only for the physical-velocity
point-value SL transport on this uniform translating-circle test. Direct normal
projection remains a useful capillary-runaway diagnostic, not the production
trajectory. Skewed/polyhedral and capillary-timestep refinement remain open.

During the stable part of a fresh `steadyUpwind` smoke run, the auxiliary
equation gives an absolute mean residual of about `1e-10` and a relative mean
residual of about `1e-14`. These values establish the Gate-2 algebra. They do
not establish Gate 4: once the nonlinear capillary runaway reaches enormous
fluxes, absolute cancellation error and auxiliary-density undershoots also
grow. Those are consequences of the runaway and can then accelerate it.

The controlled matrix is `config/transISTTransportGate.yaml`. Solver output now
records the mass residual, reset mismatch, and bandwise level-set conditioning
in `leiaSemiLagrangianLevelSetTwoPhaseFoam.csv`. Use
`workflow/scripts/check_capillary_gate.py` to turn chosen physical thresholds
into an automated pass/fail decision.

## Why Gate 2 is necessary but insufficient

The older `geometricFaceDensity` option improves the violence of a density-ratio
failure but does not enforce the cell mass balance. rhoLENT closes that algebraic
gap. It cannot remove a non-gradient capillary force. A small curvature/tangent
error generates a vortical residual, the velocity moves `psi`, and the modified
profile changes the next force. The loop is

\[
 e_\kappa\;\longrightarrow\;u_{spurious}\;\longrightarrow\;
 e_\psi,\;e_{|\nabla\psi|}\;\longrightarrow\;e_\kappa^{new}.
\]

The structured integral force reduces the initial static imbalance strongly,
but its tangents and interpolated crossing curvature still respond to the
level-set profile. Better static balance can therefore coexist with larger
dynamic feedback gain.

## N=64 lab-frame translating capillary gate

The capillary test must actually advect the interface. The solver therefore
remains in the laboratory frame with `Utrans=(0.05 0 0) m/s`; it is not solved
in a translating frame. The reported disturbance is calculated afterwards as

\[
 U_{err,\infty}=\max_c|\mathbf U_c-\mathbf U_{trans}|.
\]

`config/transISTN64QuadraticCapillaryGate.yaml` repeats the quadratic-SL,
rhoLENT test at `N=64` (12.8 cells across the initial diameter). The exact
balanced-CSF control reaches `t=0.05`; all three geometrical-force candidates
run away and encounter a floating-point exception near `t=0.001`--`0.0015`.

| force model | first-step `max|U-Utrans|` | largest completed value | last completed time | final-outer `rAU/rAUf` residual on first step |
|---|---:|---:|---:|---:|
| exact `constantCurvatureSurfaceTension` | `1.09e-8 m/s` | `2.17e-7 m/s` | `0.049999 s` | `5.77e-9 m/s` |
| structured `integralSurfaceTension` | `1.003e-2 m/s` | `4.63 m/s` | `9.33e-4 s` | `1.065e-2 m/s` |
| direct shared-face conormal cell force | `5.63e-1 m/s` | `3.99 m/s` | `9.55e-4 s` | `5.83e-1 m/s` |
| conormal-to-balanced-CSF bridge | `3.61e-2 m/s` | `9.03e-1 m/s` | `1.421e-3 s` | `3.91e-2 m/s` |

The corresponding first-step volume-weighted relative-L2 projection residuals
are `3.70e-9`, `6.78e-3`, `6.18e-1`, and `1.74e-2` in table order. Across the
full exact-control run, the maximum cell residual is `3.94e-8 m/s`.

The exact control ends with a `1.05%` DA phase-volume error and a `4.32 um`
translation-aware centroid error; its worst intermediate values are `1.91%`
and `11.3 um`, respectively. This is a real advection result over 2.5 mm, not a
static-droplet surrogate. The stable velocity together with the measurable
geometry drift also separates capillary balance from the remaining SL/DA
transport error.

Increasing from `N=32` to `N=64` reduces the structured integral model's first
impulse from about `5.65e-2` to `1.00e-2 m/s`; the coarse grid was indeed
insufficient for a quantitative judgment. Resolution does not, however, make
the present model dynamically viable: it still becomes nonlinear within the
first `1.4e-3 s`. The conormal-CSF first impulse does not improve monotonically,
and the direct conormal cell source is immediately much worse.

The new diagnostic in `applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H`
uses the actual cell `rAU`, face `rAUf`, pressure-equation flux, and OpenFOAM
reconstruction on the final PIMPLE outer/inner corrector. Its exact-control
result shows that the collocated pressure machinery can suppress a compatible
capillary force to the pressure-solver tolerance scale. The many-orders-larger
residuals of the geometric models are therefore evidence of force/gradient
incompatibility, not evidence that PIMPLE is intrinsically incapable of
capillary equilibrium. Maximum cell-centred velocity remains the coupled
acceptance observable; analytic force comparisons are optional geometry
diagnostics as described under P3.

## Cell-force-to-face-flux bridge result

`integralConormalSurfaceTension` now has two runtime deliveries of the same
conservative conormal cell force:

- `directCell`: the original negative-control source in `UEqn`;
- `faceFluxFromCellForce`: form
  `Gsigma = linearInterpolate(fSigmaCell) & Sf`, return
  `Gsigma/magSf` through `faceSurfaceTensionForce()`, and return zero through
  `cellSurfaceTensionForce()` so the force is not counted twice.

The second path enters `phig`, receives the same `rAUf` weighting as the
pressure Laplacian, and is differenced against `p_rghEqn.flux()` before the
common `fvc::reconstruct`. The reproducible comparison is
`config/transISTN64ConormalFaceFluxGate.yaml`.

At `N=64`, with physical lab-frame translation and identical conormal geometry:

| delivery | first `max|U-Utrans|` | largest completed value | first weighted residual | first relative-L2 residual | failure attempt |
|---|---:|---:|---:|---:|---:|
| `directCell` | `0.5634 m/s` | `3.989 m/s` | `0.5832 m/s` | `0.6185` | `t=0.001421` |
| `faceFluxFromCellForce` | `0.3204 m/s` | `3.538 m/s` | `0.3175 m/s` | `0.1789` | `t=0.001124` |

The correct face-flux/reconstruction chain reduces the first velocity impulse
by `43%`, the first maximum weighted residual by `46%`, and its relative-L2
value by `71%`. This validates the importance of the delivery path. It does not
validate interpolation as the final operator: both cases reach order-one
velocities, and the interpolated bridge fails slightly earlier. The averaged
flux is still not the facewise gradient of one pressure jump. The next operator
must construct the scalar capillary face flux directly from shared interface
geometry, rather than constructing a cell vector first and attempting to
recover a balanced flux afterwards.

## What [Herrmann (JCP 2008)](https://doi.org/10.1016/j.jcp.2007.11.002) changes in the roadmap

Herrmann's balanced RLSG method is especially relevant because the flow mesh
may be structured or unstructured, while all level-set geometry is evaluated
on a separate, equidistant Cartesian geometry grid. Its successful balance is
not obtained by depositing a generic cell force. The method:

1. integrates an analytic geometry-grid phase fraction into each flow control
   volume;
2. forms the capillary force as a face quantity proportional to the same
   face-centred phase-fraction difference used by the pressure jump;
3. represents capillary and pressure terms in one compatible face-normal
   space before their common cell reconstruction;
4. includes the force in the Rhie--Chow-like face-velocity interpolation and
   uses the matching variable-coefficient pressure equation--in OpenFOAM this
   means differencing the capillary face flux and `p_rghEqn.flux()` before the
   shared `fvc::reconstruct`; and
5. evaluates curvature at a normal projection to the interface and extends it
   normally, rather than using curvature of an offset level-set isosurface.

With prescribed exact curvature, that construction reports machine-level
static balance even on unstructured prism/tetrahedron flow meshes and at very
large density ratios. With computed curvature, its error and parasitic-current
magnitude converge approximately second order as the Cartesian geometry grid
is refined. The remaining non-convergent volume-error floor is attributed to
interpolating a discretely divergence-free flow velocity onto a geometry grid
where it is no longer discretely divergence-free.

This sharpens P3. The closest-point idea belongs primarily in curvature
evaluation/normal extension; it is not a justification for deleting the
tangential part of the physical trajectory. More importantly, the final force
must share the pressure operator all the way from face difference through
cell reconstruction and face-velocity correction. A separately refined
Cartesian geometry grid is also a credible research branch: it may provide
high-order geometry without forcing the entire unstructured flow mesh to carry
that resolution.

## What modern unstructured PLIC does and does not establish

The practical claim is mostly correct but needs one qualification. Modern
[RDF-PLIC reconstruction](https://arxiv.org/abs/1801.05382) can be second-order
accurate on arbitrary meshes, yet PLIC
reconstruction by itself supplies neither a second-order curvature nor a
pressure-compatible capillary operator. Published comparisons still show the
structured height-function method converging where RDF/paraboloid curvature on
triangular/tetrahedral or general unstructured meshes is less accurate and may
be zero order; see the
[TwoPhaseFlow comparison](https://journal.openfoam.com/index.php/ofj/article/download/80/107/1780).

On the other hand, unstructured collocated VOF methods with *prescribed exact
curvature* have demonstrated velocity and pressure errors at solver tolerance,
including for moving interfaces. Thus unstructured pressure/capillary balance
is possible; the unresolved general problem is the combined
PLIC--computed-curvature--balanced-force chain. Structured height functions
remain the mature benchmark for computed-curvature capillary equilibrium.
Recent specialised two-dimensional triangular height-function variants are
promising exceptions, not yet a general arbitrary-polyhedral solution.

## Priority research sequence

### P0. Keep Gate 2 fixed while researching later gates

Treat rhoLENT as the reference density coupling. Keep
`geometricFaceDensity` only as an ablation. The regular structured free-stream
control now passes velocity preservation; repeat it on skewed/polyhedral meshes
and add a serial/parallel face-flux parity test. In parallel,
replace the processor-patch owner-only face fraction with an exchanged,
single-valued coupled-face construction.

### P1. Separate zero-set accuracy from profile conditioning

Use `transISTKinematicRefinementGate` with physical `input` as the uniform-mesh
baseline and report, at identical physical times:

- geometric interface and phase-volume errors;
- `min`, `max`, ratio, and L2 error of `|grad(psi)|` in the force-support band;
- curvature L2/Linf and angular spectrum;
- velocity-extension normal defect `(n dot grad)Uext`;
- fallback counts and their distance to the interface;
- rhoLENT relative residual and geometric-reset mismatch.

`steadyUpwind` remains the materialised capillary-case default because changing
the production choice should wait for stationary and oscillating comparisons.
The new `normalProjection` trajectory has nevertheless maximised the translating
numerical-survival window under the integral force. It is not physically stable
to `t=0.03`, and the fixed-Courant refinement shows only approximately
first-order zero-set convergence with much larger volume error. Do not promote
it to the production trajectory. Physical `input` is approximately second order
on the uniform zero-force gate. The pure-advection
method-comparison result still favours quadratic SL without velocity extension
once the interface is sufficiently resolved.

### P2. Profile maintenance must preserve the discrete interface

Do not enable every-step redistancing. Existing experiments show that geometric
bulk fills can promote remote sign errors to false interfaces. Evaluate only
triggered methods whose support and invariants are explicit:

1. a zero-set-vanishing SDPLS reaction composed with the SL characteristic
   update;
2. the band biharmonic filter as a grid-mode diagnostic, with measured zero-set
   displacement;
3. threshold redistancing only when the band profile metric exceeds a scale-aware
   bound, with before/after DA phase volume and interface position recorded.

The acceptance condition is not merely `|grad psi|` closer to one: zero-set and
DA phase-volume errors must not worsen at the same resolution.

### P3. Derive the unstructured conservative capillary operator

Do not generalise the structured `(i,j)` stencil by indexing arbitrary cells,
and do not deposit independent cell conormal forces. In OpenFOAM the decisive
object is the scalar face-normal force density returned by
`faceSurfaceTensionForce()`. With gravity omitted for clarity, define

\[
 G_{\sigma,f}=f_{\sigma,f}|S_f|,
 \qquad \phi_{\sigma,f}=rAU_f G_{\sigma,f}.
\]

PIMPLE adds `phiSigma` to `phiHbyA`, the pressure Laplacian produces
`p_rghEqn.flux()`, and the velocity correction is

\[
 \mathbf U=\mathbf H/A+rAU\;\operatorname{reconstruct}
 \left(G_{\sigma,f}-\frac{\phi_{p,f}}{rAU_f}\right).
\]

Exact discrete equilibrium is therefore facewise
`p_rghEqn.flux() == rAUf*Gsigma`. The cancellation occurs in the scalar
surface-normal flux space *before* `fvc::reconstruct`. This is more precise
than merely asking a cell force to resemble a pressure gradient.

The derivation should start from a cell-integrated surface-stress identity on
the reconstructed interface polygon:

\[
 \int_{V_c}\sigma\kappa\mathbf n\,\delta_\Gamma\,dV
 = -\int_{\partial V_c}\sigma
   (\mathbf I-\mathbf n\mathbf n)\,\delta_\Gamma\cdot\mathbf n_c\,dA,
\]

then enforce shared intersection geometry and opposite tractions on every
internal face. That geometrical traction must ultimately produce one
owner-oriented scalar `Gsigma` on every face, with the dimensions and
orientation expected by the pressure correction. The direct conormal
prototype already satisfies equal-and-opposite deposition and roundoff global
force balance, while its raw identity-weighted projection leaves a `49.45%`
residual. This is strong evidence that surface-integrating the traction to a
cell source is insufficient. Its cell-source branch is a negative control, not
a candidate balanced-force delivery.

The decisive test is the coupled static/translated-droplet maximum
cell-centred velocity. The secondary `rAU/rAUf` diagnostic measures the
face-flux difference after the real pressure solve and its reconstructed
velocity increment. Analytic geometry tests may additionally compare the
numerical **cell-integrated vector force**

\[
 \mathbf F_c^h = V_c\mathbf f_{\sigma,c}^h
\]

against the exact surface-stress integral over the portion of the analytic
interface inside that cell. In pseudo-2D with constant `sigma`, if an oriented
interface arc enters and leaves the cell with exact unit tangents
`t_a` and `t_b`,

\[
 \mathbf F_c^{exact}=\sigma b\,(\mathbf t_b-\mathbf t_a),
\]

where `b` is the empty-direction thickness. Such errors can isolate geometry,
quadrature, orientation, and deposition defects, but they cannot certify
OpenFOAM force balance: only the scalar face-flux representation and the
resulting maximum cell velocity do that.

The derivation must therefore construct `Gsigma` directly, or a rigorously
equivalent traction-to-face-flux mapping, before dynamic coupling. Sending a
cell vector through `HbyA` interpolation is not that mapping.

The geometrical integral identity is fundamentally curvature-free. Curvature
is needed only by comparison models and by the paper-specific pressure-jump
correction retained in `integralSurfaceTension`; it must not be assumed to be a
required input to the final unstructured operator. If a jump correction is
retained, its samples must come from the same reconstructed geometry as the
face traction. Mixing a cell-centred curvature stencil with a different face
intersection recreates the imbalance.

Before dynamic flow, the operator must pass:

1. planar interface: exactly zero capillary force;
2. constant-curvature circle/sphere: `Gsigma` is a single-valued,
   owner-oriented surface scalar field and the pressure solve cancels it in the
   same `rAUf`-weighted face-flux space;
3. maximum cell-centred velocity: bounded and convergent to solver tolerance
   for exact geometry, together with the correct pressure jump;
4. optional geometry diagnostics: face traction and cell-integrated force
   errors distinguish geometry from delivery errors;
5. global force and torque: roundoff zero on closed interfaces;
6. orientation reversal: invariant under `psi -> -psi` with phase labels swapped;
7. mesh tests: uniform hex, perturbed hex, triangles/tetrahedra, and general
   polyhedra, including decomposed runs.

The existing `isoCurvatureSurfaceTension` conormal variants are negative-result
prototypes, not substitutes for this derivation.

### P4. Coupled acceptance ladder

Advance only after the lower gate passes:

1. frozen-interface static droplet (force/pressure only);
2. freely evolving static droplet for a time window scaled by the capillary
   period, not only an early plateau;
3. translating droplet with matched inlet/outlet velocity;
4. small-amplitude oscillating droplet with frequency and damping errors;
5. density ratios `1`, `10`, `100`, and water/air;
6. uniform, perturbed, and genuinely polyhedral meshes;
7. serial/decomposed equivalence and mesh-refinement convergence.

Every report must distinguish “completed the requested end time” from
“converged under refinement” and “late-time stable.”

## Immediate implementation record

This roadmap iteration implemented:

- normal-constant `Uext`/`UextOld` use in the two-phase SL characteristics;
- persistent `UextOld` for restart equivalence;
- `rhoLENT` as a runtime mass-flux mode;
- geometric `alphaf`/`rhof`, rebuilt once per new interface;
- the auxiliary density equation on every PIMPLE outer iteration;
- the post-PIMPLE Detrixhe--Aslam density reset;
- persistent mass-balance, reset, signed-distance, curvature, pressure-jump,
  and parasitic-current diagnostics;
- co-moving disturbance energy, phase-volume, and centroid diagnostics;
- `steadyUpwind` capillary defaults and the controlled extension/rhoLENT A/B
  configuration;
- direct `normalProjection` and `normalClosestPoint` Taylor-trajectory modes,
  with both velocity time levels and the projected-velocity gradient treated
  consistently;
- a three-case direct-trajectory configuration in
  `config/transISTTrajectoryGate.yaml`;
- zero-force, exact-curvature balanced-CSF, and integral-force translating
  mechanism gates, with the first two completing `t=0.05` and the third failing;
- a zero-force `input`/`normalProjection` kinematic ablation showing that both
  look similar only at the very small capillary-limited Courant number;
- separate Courant sweeps plus `N=32,64,128` fixed-Courant refinement, showing
  approximately second-order zero-set convergence for physical `input` and
  roughly first-order, larger errors for direct normal projection;
- signed phase-volume change and translation-aware zero-set L2/Linf metrics;
- runtime `field|constant|none` curvature-source isolation in the structured
  CST, demonstrating that exact curvature does not cure it and deleting its
  jump correction makes it worse;
- a curvature-free, shared-face conormal prototype with corrected owner-face
  orientation and roundoff global force conservation;
- a raw identity-weighted pressure-range diagnostic giving a `49.45%` residual
  for the direct conormal cell force;
- the actual final-corrector `rAU/rAUf`-weighted PIMPLE diagnostic, giving
  `5.77e-9 m/s` on the first step of the exact control but `1.06e-2`, `5.83e-1`, and
  `3.91e-2 m/s` on the first step of the three geometrical-force candidates;
- the `N=64` lab-frame translating gate, which confirms that refinement lowers
  the structured force's first impulse substantially but does not prevent its
  early nonlinear runaway;
- a conormal-to-balanced-CSF bridge that halves the first-step disturbance but
  still fails as its inferred curvature becomes noisy;
- a selectable cell-force-to-face-flux bridge for the direct conormal force;
  at `N=64` it reduces the first disturbance from `0.563` to `0.320 m/s` and
  the first weighted residual from `0.583` to `0.318 m/s`, but still runs away
  near `t=0.0011`;
- a passing regression of the legacy exact-curvature balanced-CSF control after
  adding the optional cell-force path;
- an executable CSV gate checker.

It deliberately does **not** label either the existing structured integral
force or the new conservative conormal prototype as an unstructured solution.
Gate 3 is now narrowed to constructing an accurate scalar capillary face flux
that cancels the `rAUf`-weighted pressure-Laplacian flux before the common
`fvc::reconstruct`; Gates 4 and 5 cannot honestly be claimed before it and the
profile-maintenance gate pass.

## Curvature-delivery replay gate (2026-07-27)

The strict face-flux matrix was followed by a snapshot replay that separates
curvature delivery from nonlinear flow feedback.  A quiet N=64 translating
droplet was first advected to `t = 0, 0.01, 0.025, 0.05` with exact constant
curvature, uncached quadratic SL, Detrixhe--Aslam, rhoLENT, physical `input`
trajectories, and no velocity extension.  At each saved interface, fresh cases
reset `U=Utrans`, froze only the interface advection, and ran one complete
PIMPLE timestep.  The decision observable was exclusively the maximum
cell-centred `|U-Utrans|`.

| curvature delivery | t=0 | t~0.01 | t~0.025 | t~0.05 |
|---|---:|---:|---:|---:|
| exact constant | `7.60e-9` | `5.07e-9` | `9.50e-9` | `7.37e-9` |
| computed interface mean | `7.66e-9` | `5.11e-9` | `9.60e-9` | `3.84e-9` |
| quadratic cell centre | `7.11e-3` | `8.47e-3` | `1.08e-1` | `2.03` |
| closest-point quadratic | `7.07e-3` | `8.43e-3` | `1.14e-1` | `1.25` |
| foot-point height function | `1.49e-4` | `2.98e-3` | `1.18e-1` | `1.87e-1` |
| quadratic + Kang face value | `1.19e-3` | `2.06e-3` | `1.01e-1` | `1.08` |
| FVM + Kang | `1.91e-3` | `2.89e-3` | `1.64e-2` | `2.50e-2` |

This is the cleanest mechanism result so far.  Replacing the reconstructed
curvature by its `|snGrad(alpha)| |Sf|`-weighted spatial mean restores the
pressure-solver tolerance at every transported snapshot.  Therefore spatial
curvature variation on the force support, not mean-curvature bias, is the part
that produces parasitic velocity.  Mean accuracy still matters for the
physical Laplace jump: the diagnostic must never be presented as a production
model.

The foot-point height function reduces the clean-interface impulse by about
48x relative to cell-centred quadratic curvature, with no fit fallbacks.  Its
loss of accuracy on later quiet snapshots therefore comes from sensitivity of
the reconstructed feet/fits to the transported level-set profile, rather than
from flow feedback or missing fits.  FVM--Kang is the least bad late-snapshot
variable-curvature method, but its `0.025 m/s` disturbance is still half the
imposed translation speed.

This result specified the next CSF research step: construct curvature as a
surface quantity on one connected, shared reconstruction of the zero set;
regularise the local quadratic fit in tangential coordinates; attach one
curvature degree of freedom to each interface element; and extend that value
normally to every face in the CSF support.  That step is implemented and
replayed in the following section.  Every further candidate must keep passing
the four-snapshot one-step replay before a coupled translating or oscillating
run.
The curvature-free conormal branch remains separate and still requires a
genuinely face-native traction-to-`Gsigma` derivation instead of a cell-force
interpolation.

Reproduction:

- `config/transISTN64CurvatureDeliveryReplay.yaml`;
- `workflow/scripts/run_curvature_delivery_replay.py`;
- `studies/transISTN64CurvatureDeliveryReplay/curvature_delivery_replay.csv`.

## Connected-interface curvature implementation and replay (2026-07-27)

The requested surface representation is now implemented as the 2-D serial
`connectedInterface` curvature extension.  Its zero contour is globally
continuous: `psi` is interpolated once to mesh points, every zero crossing on
an internal face is computed once, and the identical crossing node is supplied
to both adjacent cut-cell elements.  The resulting interface has 52 elements
and, at every replay snapshot, exactly one closed component, zero open ends,
52 successful fits, and zero fallbacks.

Curvature is no longer a cell degree of freedom.  Each interface element gets
one value from an ordered tangential quadratic fit.  A selectable 1-D
Helmholtz or biharmonic operator regularises those element values along the
chain, and the additive mode on a closed component is fixed with
`integral(kappa ds) = +/-2 pi`.  The element value is then extended along its
normal ray directly to the registered `kappaInterfaceFace` surface field.
`reconstructedCurvature` consumes that field when
`faceInterpolation connectedInterface`; it does not interpolate a cell field
back to the CSF faces.  The cell `kappa` field is retained only for output and
diagnostics.

The replay runner now rejects false quiet results: connected candidates must
report a non-empty, single closed component with complete fit coverage.  This
caught and prevented promotion of an early implementation whose face field was
zero because no elements survived reconstruction.

| connected candidate | t=0 | t~0.01 | t~0.025 | t~0.05 |
|---|---:|---:|---:|---:|
| Helmholtz, half-width 3, lambda 1 | `1.67e-3` | `1.66e-3` | `1.79e-2` | `5.01e-2` |
| Helmholtz, half-width 4, lambda 4 | `9.98e-4` | `7.05e-4` | `7.88e-3` | `1.60e-2` |
| **Helmholtz, half-width 5, lambda 16** | **`6.01e-4`** | **`4.38e-4`** | **`1.79e-3`** | **`4.55e-3`** |
| biharmonic, half-width 3, 20 passes | `1.74e-3` | `9.44e-4` | `1.62e-2` | `5.21e-2` |
| biharmonic, half-width 3, 50 passes | `1.31e-3` | `1.18e-3` | `8.47e-3` | `3.50e-2` |

The strongest Helmholtz candidate is the replay winner.  At the last snapshot
it is 5.5 times quieter than FVM--Kang, 41 times quieter than the independent
foot-point fit, and 447 times quieter than cell-centred quadratic curvature.
It also keeps the late transported-contour element curvature bounded to about
`[978,1034] 1/m`, compared with sign-indefinite cell-local excursions.

This is not yet permission to claim an oscillating-droplet method.  On a
52-element circle, lambda 16 attenuates the physical azimuthal mode `m=2` by
roughly one half before the local-fit transfer function is counted.  The
biharmonic candidates preserve low modes better but did not match the late
replay velocity.  Therefore no coupled translating or oscillating run was
started in this iteration.  The next mandatory gate is a manufactured
perturbed-circle spectrum: measure the transfer function and curvature error
for modes `m=2,3,4` while retaining the four-snapshot velocity replay.  Only a
regularisation choice that is both replay-quiet and low-mode-accurate should
enter the coupled acceptance ladder.

## Shared face-curvature service and manufactured mode gate (2026-07-28)

The face-native curvature path is now shared by every scalar-curvature CSF
family.  `surfaceTensionForce` provides a common
`faceCurvatureSource model|registered` lookup and the common integrated
assembly

`G_sigma,f = sigma kappa_f snGrad(weight) magSf`.

`reconstructedCurvature`, `divGradPsiSnGradAlpha`, both trace models,
`correctionKang`, `divGradAlphaSnGradAlpha`, and the CSF branches of
`isoCurvature` can all consume `kappaInterfaceFace`.  Kang interpolation is
bypassed when the registered face field is selected.  The integral-conormal
model remains a separate traction-to-flux problem, and exact constant
curvature remains the pressure-coupling oracle.

The four-snapshot replay verifies the abstraction: reconstructed, divergence,
both trace variants, Kang, and iso-CSF return identical maximum cell-centred
`|U-Utrans|` to the printed digits when they use the same connected curvature
and `snGrad(alpha)`.  The alpha-CSF result differs because that model retains
its smoothed-alpha force localisation; at the last snapshot it gives
`0.1139276 m/s` versus `0.0520880 m/s` for the common geometric-alpha group.

The manufactured gate is implemented by
`leiaTestConnectedCurvatureModes` and
`workflow/scripts/run_curvature_mode_transfer_gate.py`.  It evaluates the
production connected reconstruction on

`r(theta) = R [1 + 0.02 cos(m theta)]`, `m=2,3,4`,

at `N=32,64,128` on uniform and deterministic 10%-perturbed meshes.  For
`N>=64`, the preregistered limits are curvature-amplitude transfer in
`[0.8,1.2]`, phase error below 10 degrees, relative L2 below 0.25, and valid
single-component closed topology.  Detailed and candidate-summary tables are
written into the method-comparison data directory.

The scalar filters expose a strict tradeoff: the replay-quiet half-width-5,
lambda-16 Helmholtz filter has transfer as low as `0.1556`; the 20-pass
half-width-3 biharmonic filter passes the mode gate but gives `0.05209 m/s` on
the late replay.  Widening or increasing the biharmonic filter crosses the
fixed `m=4` transfer limit before becoming replay-quiet.

A separate constrained operator, `helmholtzPreserveModes`, was therefore
added.  It applies strong graph Helmholtz damping, then restores the weighted
arclength Fourier projection of modes `m=2..4` on every closed component;
Gauss--Bonnet still fixes the mean.  The half-width-4, lambda-16, 120-iteration
candidate passes the manufactured gate: over all `N>=64` tests its amplitude
transfer is `0.8190..1.1342`, maximum phase error is `1.426 degrees`, and
maximum relative L2 is `0.0390`.  Its late replay is `0.0128448 m/s`, the best
among the mode-preserving candidates but still not pressure-tolerance quiet.

The candidate was advanced to the full N=64 translating run.  It reached
`t=0.05`, but `max|U-Utrans|` grew to `0.08933 m/s`, transient absolute volume
error reached `2.376%`, and the final zero-set radial Linf error was
`0.2038 mm`.  This fails the coupled translating gate.

The requested N=64 oscillating follow-up was nevertheless executed.  This
first exposed a benchmark defect: the legacy `implicitEllipsoid` returns the
algebraic function `x^2/a^2+y^2/b^2-1`, not a metric signed distance, and the
initial interface band consequently had `|grad(psi)|` of order `10^3`.
The new backward-compatible `signedDistanceEllipse` evaluates the Euclidean
distance to the closest ellipse point.  With it, the first-step band range is
`0.9877..0.99997`, so Detrixhe--Aslam localisation and SL transport start from
the intended field.

Even with that correction, the candidate fails the oscillating test.  The run
ends with `SIGFPE` at `t=0.0209 s` instead of `0.1 s`.  Before failure,
`max|U|=15.46 m/s`, the mode-2 amplitude grows from `0.0949 mm` to a maximum
`0.5121 mm`, volume error reaches `4.965%`, and the level-set band degrades to
`|grad(psi)|=0.0858..7.087`.  Four zero crossings give a mode-2 period of
`8.801 ms`, versus the inviscid cylindrical-drop value `9.508 ms` (`-7.44%`),
but the subsequent amplitude growth makes this an unstable, nonphysical
oscillation rather than a usable frequency result.

The next research issue is no longer the force-model class or Kang
interpolation.  It is the incompatibility between (a) preserving real low
curvature modes and (b) treating every low-mode distortion in the transported
interface as a parasitic-force defect.  The next controlled test should be a
variable-curvature velocity oracle: deliver exact analytic ellipse or
perturbed-circle curvature to the same faces, run the same frozen one-step
PIMPLE solve, and compare the candidate's cell-centred velocity with the
oracle velocity.  This retains the physical mode-2 response while measuring
only the velocity error created by curvature delivery/localisation/coupling.
Repeat on analytically prescribed translated shapes before further filter
tuning or another long coupled oscillation.

## Exact variable-curvature velocity oracle (2026-07-28)

Build and execution provenance: the repository was cleaned globally with
`./Allwclean`, rebuilt completely with OpenFOAM v2512, and both the velocity
oracle and manufactured mode-transfer gate were rerun under that v2512 build.

The next gate is implemented and executed.  The connected-interface service
now accepts `curvatureExtension.estimator analyticImplicitSurface`.  It builds
the same shared nodes, ordered elements and normal-ray face ownership as the
candidate, but attaches the exact curvature of the configured analytic
surface to each element.  Those oracle values bypass fitting, tangential
regularisation and the discrete Gauss--Bonnet shift.  The production force,
geometric alpha, rhoLENT, pressure tolerances and PIMPLE loop remain identical.

`run_variable_curvature_velocity_oracle.py` advances the analytic ellipse
oracle and `connectedFit+helmholtzPreserveModes` candidate for one frozen
PIMPLE step, reads the two written cell-centred velocity fields, and reports
only `max|U_candidate-U_oracle|`.  The matrix covers N=32,64,128 on uniform and
deterministic 10%-perturbed meshes.  After the pressure-correction audit, all
perturbed rows were rerun with eight non-orthogonal correctors; their reported
values reproduced to the printed precision.

On uniform meshes, the maximum velocity difference decreases from
`0.013739` through `0.004023` to `0.001654 m/s`.  At N=64 the candidate and
oracle maxima are `0.040079` and `0.036987 m/s`, respectively, and their direct
cell-field discrepancy is 10.88% of the oracle maximum.  The relative error is
not monotone (`47.5%, 10.9%, 19.7%`), so this is improvement rather than a
verified asymptotic rate.

Perturbed meshes expose the controlling defect.  Even the analytic-curvature
oracle gives maximum velocities `2.364, 3.673, 0.7119 m/s` at N=32,64,128,
which are approximately `81.7, 99.3, 85.0` times the corresponding uniform
oracle responses.  Candidate--oracle differences are `2.230, 0.4643,
0.7118 m/s`.  The oracle and candidate mesh point files are byte-identical for
each pair, so this is not a mesh-realisation mismatch.

Exact curvature therefore does not solve the skew-mesh problem.  The current
oracle still reconstructs point values from cell-centred psi, obtains
least-squares normals, constructs Detrixhe--Aslam tangent planes, and uses
`snGrad(alpha)`, density and `rAUf`.  Even though the analytic input is a true
signed distance, the minimum reconstructed band `|grad(psi)|` on perturbed
meshes is `0.569, 0.510, 0.463` with refinement, whereas the uniform values
approach one.

The next gate must keep the same velocity-only observable while replacing the
remaining geometric approximations: evaluate analytic psi at mesh vertices
and analytic normals at interface/cell points, use them in both connected
geometry and Detrixhe--Aslam planes, and rerun the matrix.  If the perturbed
oracle velocities collapse toward the uniform values, geometry/localisation
is the cause.  If they remain large, the pressure-gradient/`rAUf` coupling is
the next isolated subsystem.

## Analytic geometry/localisation oracle (2026-07-28)

The requested geometric separation is implemented and executed under
OpenFOAM v2512.  `curvatureExtension.geometrySource` and
`phaseIndicator.geometrySource` now accept `analyticImplicitSurface` as an
oracle-only alternative to the default `levelSetField` path.

The analytic case removes the following numerical geometry inputs before the
Navier--Stokes solve:

- the configured signed-distance ellipse is evaluated directly at mesh
  vertices, so connected face crossings do not use cell-to-point interpolation;
- connected-element normals are evaluated from the analytic surface gradient
  at each element centre instead of a cell-centred least-squares plane;
- Detrixhe--Aslam cell fractions use the analytic signed distance and normal at
  the cell centre to form their local tangent plane; and
- rhoLENT face area fractions use the same analytic tangent planes, removing
  the least-squares plane from the face-density path as well.

The experiment deliberately retains `snGrad(alpha)`, density weighting,
`rAUf`, the pressure Laplacian and the common `fvc::reconstruct`.  The runner
now executes three variants on each mesh: the connected mode-preserving
candidate, exact curvature on the numerical geometry, and exact curvature
plus analytic geometry/localisation.  All 18 frozen one-step solves completed;
all reconstructions retained one closed component with zero fallback.

Maximum cell-centred velocities `(exact curvature on numerical geometry,
exact curvature plus analytic geometry)` [m/s] are:

| N | uniform | 10%-perturbed |
|---:|---:|---:|
| 32 | `(0.02895, 0.03858)` | `(2.364, 3.921)` |
| 64 | `(0.03699, 0.02847)` | `(3.673, 3.263)` |
| 128 | `(0.008376, 0.006917)` | `(0.7119, 1.129)` |

The direct cell-field changes caused by analytic geometry on perturbed meshes
are `2.245, 1.281, 0.5924 m/s`, respectively.  More importantly, the fully
analytic-geometry oracle remains approximately `101.6, 114.6, 163.2` times
larger on the perturbed meshes than on the corresponding uniform meshes.
Analytic vertices and normals therefore do not cure the skew-mesh response;
at two resolutions they increase it.

This rejects vertex interpolation and least-squares normals as the controlling
defect.  It does not prove that the Detrixhe--Aslam tangent-plane approximation
is exact; it proves that replacing its numerical normal and offset by analytic
geometry is insufficient.  The remaining load-bearing subsystem begins with
the discrete localisation `snGrad(alpha)` and its coupling to the
`rAUf`-weighted pressure equation.

The next gate should use a frozen analytic circle with exact constant
curvature and zero initial velocity on N=32,64,128 uniform and perturbed
meshes.  Compare the current integrated CSF flux
`sigma*kappa*snGrad(alpha)*magSf` with a pressure-potential oracle assembled as
`snGrad(sigma*kappa*alpha)*magSf`; pass both through the identical pressure
equation and `fvc::reconstruct`, and continue to measure only maximum
cell-centred velocity.  For constant `kappa` the two fluxes should be the same
discrete gradient.  A quiet potential oracle with a noisy CSF path identifies
localisation assembly; both noisy on perturbed meshes identifies the
pressure/`rAUf`/reconstruction chain.

## Frozen-circle pressure-compatibility gate (2026-07-28)

The constant-curvature gate is implemented and executed under OpenFOAM v2512.
The new production model `constantCurvaturePressurePotential` constructs the
cell-centred capillary pressure `pSigma = sigma*kappa*alpha` and returns only
the strict integrated oriented scalar face flux
`snGrad(pSigma)*magSf`.  The comparison model remains
`constantCurvatureSurfaceTension`, which returns
`sigma*kappa*snGrad(alpha)*magSf`.  Both enter the unchanged
`rAUf`-weighted pressure equation and the same `fvc::reconstruct` velocity
recovery.

`workflow/scripts/run_pressure_compatibility_gate.py` freezes a 1 mm analytic
circle at zero initial velocity, uses analytic Detrixhe--Aslam geometry and
rhoLENT-central mass flux, and advances one capillary step at N=32,64,128 on
uniform and deterministic 10%-perturbed meshes.  The only decision observable
is maximum written cell-centred velocity.

| N | uniform CSF | uniform potential | perturbed CSF | perturbed potential | max direct field difference |
|---:|---:|---:|---:|---:|---:|
| 32 | `3.770e-9` | `3.770e-9` | `8.584e-6` | `8.584e-6` | `2.02e-14` |
| 64 | `7.635e-9` | `7.635e-9` | `7.328e-4` | `7.328e-4` | `4.25e-14` |
| 128 | `1.039e-8` | `1.039e-8` | `1.392e-3` | `1.392e-3` | `1.50e-13` |

For constant curvature the two complete velocity fields agree to roundoff.
The product form of the CSF localisation is therefore not a distinct source
of error in this gate.  Uniform meshes are quiet at approximately `1e-8 m/s`,
whereas the perturbed-mesh velocity grows with refinement and reaches
`1.392e-3 m/s`.  Exact curvature, analytic geometry and a capillary pressure
potential are all insufficient on the corrected skew-mesh pressure path.

## Perturbed-mesh pressure-correction audit (2026-07-28)

The coupled test history did not use perturbed meshes throughout.  The
stationary, translating and oscillating droplet configurations all specified
`mesh: hex`.  In the oscillating case, "perturbed" described the mode-2
ellipse, not mesh non-orthogonality.  Explicit uniform/perturbed mesh pairs
entered later in the curvature-mode gate (which has no pressure solve), the
variable-curvature velocity oracle and the frozen-circle pressure gate.  The
two CFD oracle runners initially inherited
`nNonOrthogonalCorrectors 1` from the case template.

`workflow/scripts/run_nonorthogonal_correction_sweep.py` now sweeps
`0,1,2,4,8,16,32,64` corrections on the perturbed exact-circle
pressure-potential cases at N=32,64,128.  Maximum written cell velocity [m/s]
is:

| N | 0 | 1 | 2 | 4 | 8 | 16 |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | `1.358e-5` | `1.006e-5` | `8.536e-6` | `8.583e-6` | `8.584e-6` | `8.584e-6` |
| 64 | `1.370e-3` | `7.519e-4` | `7.378e-4` | `7.328e-4` | `7.328e-4` | `7.328e-4` |
| 128 | `1.183e-2` | `1.146e-3` | `1.341e-3` | `1.388e-3` | `1.392e-3` | `1.392e-3` |

The values at 8,16,32,64 are identical to the relevant reported precision.
One correction was under-converged; eight is sufficient for these meshes.
However, exact convergence of the non-orthogonal correction does not remove
the parasitic current.  At N=128 it converges to a larger value than the
one-correction result.  Therefore an under-converged correction contaminated
the earlier pressure-gate numbers but was not the controlling mechanism.

The repository now enforces this result in three places:

- `oscISTDroplet2D` and `transISTDroplet2D` expose
  `N_NON_ORTHOGONAL_CORRECTORS` as a materialised parameter;
- general `mesh: perturbed` materialisation raises that value to at least 8;
- the pressure and variable-curvature oracle runners explicitly use 8 on
  perturbed meshes and abort if configured below it.

The next controlled step is the paired operator gate.  Compare corrected
`snGrad`/Laplacian at eight corrections with a paired uncorrected
`snGrad`/Laplacian construction.  If the orthogonal pair is quiet, the
corrected face-flux construction is not exactly compatible.  If it is not,
replace physical `rAUf` by a constant oracle before testing
`fvc::reconstruct` separately.

No variable-curvature model should return to a translating or oscillating
droplet until this exact-circle perturbed-mesh velocity is reduced to the
pressure-solver tolerance scale.

## Pressure operator, `rAUf` and linear-solver gates (2026-07-28)

The pressure programme is now a dedicated Snakemake DAG:

```bash
snakemake -s workflow/Snakefile.pressure-compatibility \
  --workflow-profile profiles/local
```

The similarly named Make targets are thin aliases; Snakemake owns inputs,
freshness, rule dependencies and generated tables.

The paired corrected/uncorrected operator gate used eight non-orthogonal loops
for both internally consistent pairs.  Maximum cell velocities `(corrected,
uncorrected)` [m/s] were `(8.584e-6,1.044e-6)`,
`(7.328e-4,7.537e-4)` and `(1.392e-3,1.815e-3)` at N=32,64,128.  The
uncorrected pair helps only at N=32 and is 30% worse at N=128.  Corrected
face-normal interpolation is not the sole cause.

The 2x2 physical/constant-`rAUf` gate also rejected coefficient variability as
the controlling mechanism.  Constant `rAUf` improves the N=32 corrected case
from `8.584e-6` to `3.163e-6 m/s`, but worsens N=64 from `7.328e-4` to
`8.582e-4 m/s` and N=128 from `1.392e-3` to `1.603e-3 m/s`.  The same trend
holds for the uncorrected pair.

Pressure algebra is load-bearing.  Production GAMG gives
`(8.584e-6,7.328e-4,1.392e-3) m/s`.  Strict PCG/DIC at tolerance `1e-11`,
`relTol 0` gives `(9.050e-6,4.029e-5,4.839e-5) m/s`: essentially unchanged
at N=32, but 18.2x and 28.8x smaller at N=64 and N=128.  Tightening PCG from
`1e-11` to `1e-13` leaves the velocity unchanged, so the remaining
approximately `4.8e-5 m/s` response is not an algebraic tolerance floor.
Strict GAMG is unusable as the oracle: it worsens N=64 and terminates with a
floating-point exception at N=32 and N=128.

The next gate must use PCG/DIC, verify its non-orthogonal-correction convergence,
then bypass `YoungLaplaceEqn.H` by directly initializing the manufactured
pressure `p_rgh = sigma*kappa*alpha` (up to the reference constant).  If the
frozen one-step maximum cell velocity collapses, the startup pressure solve is
responsible.  If the PCG plateau remains, isolate the common
`fvc::reconstruct` projection.  Variable-curvature, translating and
oscillating cases remain blocked until this manufactured gate is quiet.

## The momentum predictor, and why it does not compose with the linear solver (2026-07-30)

### The manufactured-pressure-initialisation gate declared above is null

Before running it: the solved `p_rgh` already *is* the manufactured field.  On the
perturbed exact-circle cases the written `p_rgh` differs from
`72.74*alpha.water` plus a reference constant by at most `5.35e-6` Pa at
N=128, and `pEqn.H` already prints the exact observable the gate proposes to
measure -- `rAU/rAUf capillary-pressure velocity residual: max=0.00139183637045
m/s` -- matching the written `max|U|` to seven digits.  Because the one-step run
re-solves `fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)` from scratch, a
manufactured initial field acts only as a Krylov initial guess plus, through
`momentumPredictor yes`, the first outer iteration's `HbyA`.  Whoever runs it will
spend the build confirming a null.  It is superseded by the facewise gate at the
end of this section, and the standing block on droplet runs should be re-scoped
accordingly rather than waiting on it.

### The predictor was hardcoded and never varied

`momentumPredictor yes` was fixed in all sixteen droplet dictionaries and appears
nowhere in the negative-result record.  With the predictor on, the capillary force
reaches the velocity twice per outer iteration: once as a cell source obtained by
`fvc::reconstruct` of the face-normal force inside the momentum predictor, where
no pressure gradient yet exists to cancel against, and once through the pressure
correction, where it enters as a facewise *difference* that nearly vanishes for a
compatible force.  The balanced-force cancellation is exact in the scalar
face-normal flux space, so the predictor is the one place in the loop that pushes
the un-cancelled force, of magnitude `sigma*kappa/h ~ 4.6e5` N/m^3, through a
reconstruction at full strength.  It is now the materialised token
`MOMENTUM_PREDICTOR` (default `yes`, so every existing study is unchanged), with a
`LEIA_MOMENTUM_PREDICTOR=yes|no` override honoured by the gate runners.

Frozen exact-circle pressure-compatibility gate, exact constant curvature,
analytic geometry, production GAMG, eight non-orthogonal correctors.  Maximum
written cell velocity [m/s], predictor on -> off:

| N | uniform | perturbed |
|---:|---|---|
| 32 | `3.770e-9` -> `6.766e-9` | `8.584e-6` -> `2.725e-6` |
| 64 | `7.635e-9` -> `6.770e-9` | `7.328e-4` -> `3.137e-5` |
| 128 | `1.039e-8` -> `4.081e-9` | `1.392e-3` -> `6.121e-5` |

Two distinct effects.  The perturbed-mesh response falls by `23x` on both
resolved meshes, comparable to the linear-solver gain reported above and reached
by a one-word dictionary change.  And the uniform sequence stops growing under
refinement and starts falling, so on orthogonal meshes the predictor was the
entire source of the (small) refinement growth.  The product and potential force
forms still agree to roundoff, `4e-16` to `1e-13`.

### The predictor is irrelevant to curvature-variation-driven currents

The four-snapshot replay was rerun with the predictor off, reusing the existing
transported baseline so the `psi` snapshots are byte-identical and only the
one-step solve differs.  At `t=0.05`, predictor on -> off:

| curvature delivery | on | off |
|---|---:|---:|
| exact constant | `7.37e-9` | `1.33e-9` |
| interface mean | `3.84e-9` | `6.94e-10` |
| quadratic cell centre | `2.03` | `2.031` |
| closest-point Newton | `1.25` | `1.252` |
| foot-point height function | `1.87e-1` | `1.870e-1` |
| quadratic + Kang | `1.08` | `1.083` |
| FVM + Kang | `2.50e-2` | `2.496e-2` |
| mode-preserving connected fit | `1.28448e-2` | `1.28448e-2` |

Every variable-curvature arm is unchanged to four significant figures; only the
two controls move, and only within solver tolerance.  This dissociates the two
incompatibilities.  The predictor's extra reconstruction amplifies
**mesh-induced** force/gradient incompatibility, which is the entire residual when
the curvature is exact and constant.  It does nothing for
**curvature-variation-induced** incompatibility, because there the force genuinely
lies outside the range of the discrete gradient and removing a reconstruction pass
cannot change the range of an operator.  Tangential curvature variation remains
the velocity-producing defect of the coupled problem, and it sits upstream of the
predictor.

### The two levers do not compose: a common floor

Registered before running: if the predictor and the linear solver attack
independent channels the combination should reach about `2e-6` at N=128; if they
attack the same residual both land on a common floor near `5e-5`.  Perturbed mesh,
exact constant curvature, eight correctors, `strictPCG1e11` (PCG/DIC, tolerance
`1e-11`, `relTol 0`).  Maximum cell velocity [m/s]:

| N | GAMG, pred. on | PCG, pred. on | GAMG, pred. off | PCG, pred. off |
|---:|---:|---:|---:|---:|
| 32 | `8.584e-6` | `9.050e-6` | `2.725e-6` | `1.914e-6` |
| 64 | `7.328e-4` | `4.029e-5` | `3.137e-5` | `3.509e-5` |
| 128 | `1.392e-3` | `4.839e-5` | `6.121e-5` | `7.658e-5` |

The second branch holds.  At N=64 all three corrected variants land in
`3.1e-5`--`4.0e-5`; at N=128 the combination is `1.58x` *worse* than PCG/DIC
alone.  Only at N=32, where the drop spans `6.4` cells, is the combination the
best of the four.  Both levers were therefore removing amplification of one and
the same underlying residual, which is now exposed at roughly `3.5e-5` at N=64 and
`7.7e-5` at N=128 -- still growing with refinement, and about `O(h^-1)` over that
pair.

The exposed floor is not algebraic.  PCG/DIC at tolerances `1e-9`, `1e-11` and
`1e-13` agrees to better than `0.4%` in every row (`1.901/1.914/1.914e-6`,
`3.510/3.509/3.509e-5`, `7.646/7.658/7.657e-5`).  Strict GAMG remains unusable
with the predictor off as well: floating-point exception at N=32 and N=128, and
`1.172e-3` at N=64.

### What is left

For an exact, constant curvature on a skewed mesh the residual is now insensitive
to the force form (product versus pressure potential, roundoff), the geometry
(analytic vertices, normals and tangent planes), the non-orthogonal corrector
count (converged at eight), `rAUf` variability, the momentum predictor, the
pressure solver, and the solver tolerance.  One operator is common to all four
cells of the table above and has never been removed: the `fvc::reconstruct` that
recovers the cell velocity from the face-normal flux difference in the pressure
correction.  It assumes the face values are samples of one smooth cell-centred
field and that `sum_f Sf (x) Sf` is well conditioned; at an interface cell in a
water/air system the first fails, because the velocity is only C0 there while
`grad(p)/rho` jumps by the density ratio, and on a skewed mesh the second fails
too, so the inverse mixes tangential into normal information.  That is consistent
with the observed signature: the same operator gives `1e-8` on uniform meshes and
`1.4e-3` on 10%-perturbed ones with exact curvature and analytic geometry.

The next gate separates source from gain and is a print statement.  Measure the
facewise residual `phig - p_rghEqn.flux()` *before* the common
`fvc::reconstruct`, and report it alongside the post-reconstruct maximum cell
velocity, on the frozen exact circle, uniform versus 10%-perturbed, N=32/64/128.
If the face residual sits at solver tolerance while the cell velocity is `1e-5`
or larger, the reconstruction is the gain and the fix is the velocity-recovery
operator -- an interface-aware or least-squares reconstruction, or keeping the
correction in flux space.  If the face residual is itself of that order, the
capillary and pressure operators do not cancel facewise on skewed meshes and the
reconstruction is innocent.

`momentumPredictor no` should not yet be adopted as a default.  It is `1.8x`
worse at N=32 uniform, both values being at tolerance, and with the predictor off
the viscous and convective terms act only through `HbyA`, which changes the
effective time integration.  That is benign for a frozen gate and a near-static
drop, but it must be reassessed before any oscillating-droplet frequency or
damping number is quoted.

Reproduction: `LEIA_MOMENTUM_PREDICTOR=no` with
`workflow/scripts/run_pressure_compatibility_gate.py`,
`workflow/scripts/run_pressure_solver_gate.py` and
`workflow/scripts/run_curvature_delivery_replay.py` (the last one without
`--fresh`, so the transported baseline is reused).  The predictor-off tables are
kept out of `docs/**/data/` so the committed tables continue to describe the
predictor-on production configuration.

## Face-centered curvature convergence gate (2026-08-06)

The replay gates measure how quiet a curvature delivery leaves the flow; this
gate measures how accurate the delivered quantity itself is.  The observable is
the face curvature the CSF force actually applies, `kappa_f` in
`G_sigma,f = sigma kappa_f snGrad(alpha) |Sf|`, on the active
`|snGrad(alpha)| > 0` faces of the static R = 1 mm circle initialised as an
exact signed distance, N = 32..512.  `leiaTestMeanCurvature` now assembles
`kappa_f` for every model exactly as its `surfaceTensionForce` path does and
writes tidy `leiaTestFaceCurvature.csv`; the study is
`config/faceCurvatureDroplet2D.yaml` (serial, seconds), the figure and order
table land in `docs/method-comparison/method-comparison-article/data/`.

| face curvature (`|Sf|`-weighted band L2) | order | L2 at N=512 [1/m] |
|---|---:|---:|
| quadratic cell centre, arithmetic (production CSF) | `h^1.13` | `11.35` |
| + stabilized foot point at the face | **`h^2.04`** | **`0.105`** |
| FVM div(grad psi/|grad psi|) | `h^1.16` | `11.36` |
| + stabilized foot point | `h^1.97` | `0.315` |
| Kang face interpolation | `h^1.14` | `8.42` |
| Newton-foot / closest-point extensions | `h^1.08` / `h^1.13` | `11.35` |
| foot-point height function | `h^1.75` | `0.56` |
| stabilized foot point, per-side variant | `h^1.40` | `2.02` |
| connected interface (defaults w=3, lambda=1) | `h^0.48` | `4.30` |
| interface mean (constant diagnostic) | `h^2.03` | `0.15` |
| FVM -div(grad alpha/|grad alpha|), 2x smoothed | diverges (`h^-0.4`) | `824` |

Three mechanism results.  First, the O(h) defect of the production face
curvature is the parallel-contour offset, not the estimator: the Newton-foot
and closest-point normal extensions leave the face error unchanged because a
quadratic fit's Hessian is constant per cell, so its curvature is pinned to the
stencil centre no matter where the foot lands.  Second, re-referencing the SAME
interpolated `kappa_f` to the interface through the stabilized foot point
(`footPointDistance` of the face centre, inverse-|psi| side weighting, then the
parallel-curve inverse `kappa_d/(1 - d_f kappa_d)`) restores clean second
order -- a 100x accuracy gain at N=512 for one extra scalar per active face.
Applying the same correction per side BEFORE the face combination is worse
(`h^1.4`): the correction must act on the face-combined value.  Third, the
CSF-support-weighted mean converges at `h^2`, so mean-curvature bias is second
order even when the facewise field is first order -- consistent with the replay
verdict that spatial variation, not bias, drives the parasitic velocity.

Pitfall found on the way, worth its own line: `leastSquaresPlaneCoeffs`
guarded its 4x4 normal-equations solve with `|det| < 1e-10 prod(diag)` on the
ABSOLUTE-coordinate matrix.  That ratio shrinks like `(h/|x|)^4` under
refinement, so at N = 512 on the 0.01 m box the guard silently returned the
flat-plane fallback for EVERY band cell: the geometric/Detrixhe-Aslam
indicators degraded to `sign(psi)` (the snGrad(alpha) band halved), and the
plane-based curvature models (foot-point height function, connectedInterface)
fell back to the symbolic cell value with no error or warning.  The guard now
runs on a centered, variance-normalised copy of the matrix -- a scale-free
shape test of the stencil cloud -- while the solve is untouched, so every
previously passing cell is byte-identical (verified: the refreshed
`curvatureDroplet2D` table changed in the N=512 row only).  The committed
`curvature_error.csv` N=512 row was biased LOW by the narrow corrupted band
(10.2 vs the honest 17.4 on the full band); both themes' tables are refreshed.

The connected-interface stall (`h^0.48` with the default half-width-3,
lambda=1 chain fit) is a static-accuracy counterpart to its replay quietness
and needs the manufactured-spectrum gate before any parameter is promoted.

Reproduction: `snakemake --workflow-profile profiles/local --configfile
config/faceCurvatureDroplet2D.yaml`; slide in
`docs/method-comparison/method-comparison-presentation/comparison.template.html`.
