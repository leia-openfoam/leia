# `fvc::grad(U)` was biased at every processor boundary — and the SL Taylor foot uses it

Handover note from the SDPLS / velocity-extension session (`github-51`) to the
semi-Lagrangian session. Written 2026-08-26. Everything below is measured or is
a file:line reference; nothing is estimated.

**Short version.** A defect in `velocityModel::setVelocity` biased `fvc::grad(U)`
by an O(1) — not O(h) — amount in every cell adjacent to a processor boundary,
for as long as it has existed. It is fixed (`30e6ba9`). The semi-Lagrangian
**Taylor foot integrator consumes `fvc::grad(U)`** in its second-order term, and
the affected solver feeds it a *prescribed* velocity, so **31 parallel kinematic
SL studies are contaminated, 15 of them with curated tables already in the
versioned docs.** An earlier audit of this bug's blast radius missed this path
and therefore under-reported it by an order of magnitude.

---

## 1. The defect

`src/leiaLevelSet/velocityModel/velocityModel.C`, `setVelocity()`, wrote

```cpp
UpatchField[faceI] = velocity(CfPatchField[faceI]);   // a FACE-CENTRE value
```

into **every** patch, including coupled (processor, cyclic) ones.

A coupled patch field does not hold the value at the face. It holds the value in
the **neighbouring cell centre** — that is what `processorFvPatchField::
patchNeighbourField()` returns, and what every coupled-patch operator
interpolates against. Gauss gradient forms `w·U_own + (1−w)·U_neighbour` and
expects the second term to be a *cell* value.

On a uniform mesh with `w = 1/2` the interpolation therefore returns
`0.5·U(x_P) + 0.5·U(x_f)` instead of `U(x_f)` — an error of
`0.5·(U(x_f) − U(x_P)) = a·h/4` for a linear field of strain `a`. Divided by the
cell volume in the Gauss sum that is an **O(a), not O(h), error in `grad(U)`**
in every processor-adjacent cell. It does not vanish under refinement.

**Fixed in `30e6ba9`**: coupled patches are skipped in the analytic loop and
filled by `U.correctBoundaryConditions()`, which returns the true neighbour cell
values. That is also strictly more accurate than the analytic face value would
have been, because a coupled face then interpolates exactly as an interior face.

### How it was found

`cases/1Dstretch` — a 1D uniaxial-stretch gate with a closed-form solution for
the whole ψ field (`v = a(x−x₀)e_x`, `a = 1 1/s`). Band mean `d(psi)/dx` at
`t = 1`, against the exact `1` (source `R`) and `e^{-1}` (sourceless):

| arm | serial | np=2 | np=4 | np=8 |
|---|---|---|---|---|
| noSource | 0.367872 | 0.367872 | 0.367872 | 0.367872 |
| R (before) | 1.000000 | 1.004681 | 1.005367 | **0.956834** |
| R (after) | 1.000000 | 1.000000 | 1.000000 | 1.000000 |

The error **grew with processor count** and the sourceless arm was clean at every
count — because advection uses `phi`, which `setVolumetricFlux` writes per-face
where the face value *is* the correct flux. Only consumers of `fvc::grad(U)` were
affected.

---

## 2. The path that was missed: the SL Taylor foot

`src/leiaLevelSet/semiLagrangian/pointValueScheme.C:353`

```cpp
const volTensorField gradU(fvc::grad(*trajectoryNew, "gradU"));
...
const vector accel = (uNew - uOld)/dt + (uNew & gradU[c]);      // :367
feet[c] = C[c] - uNew*dt + 0.5*accel*dt*dt;                     // :368
```

i.e. the `dt²/2` term of

```
x_d = x_c - u^{n+1} dt + 1/2 [ du/dt + (u^{n+1}.grad)u^{n+1} ] dt^2
```

is built from a velocity gradient.

Three facts make this the dominant exposure rather than a corner case:

1. **`pointValue` is the DEFAULT scheme** — `cases/default.parameter:134`,
   `SL_SCHEME pointValue`. Most studies never name it.
2. **`leiaSemiLagrangeLevelSetFoam` feeds a PRESCRIBED velocity** —
   its `createFields.H` calls `velocityModel->setVelocity(U)`, so the `U` whose
   gradient enters the foot trace is exactly the field the defect corrupted.
3. The two-phase SL solver is **not** affected: it takes `U` from the momentum
   solve, and `leiaSemiLagrangianLevelSetTwoPhaseFoam` contains no reference to
   `velocityModel` at all.

An earlier audit of this bug listed `normalProjectedScheme.C:161` and
`signedDistanceQuadraticWeightedLeastSquaresReconstruction.C:319` as the SL-side
`grad(U)` consumers, and therefore flagged only `npslConv2Dvortex` and
`sdCompare2D`. It did not check `pointValueScheme.C`. Because that is the
default, the true set is the whole parallel kinematic SL benchmark suite.

---

## 3. Affected studies

Criterion: `solver: leiaSemiLagrangeLevelSetFoam` **and** `np > 1`. Serial studies
are clean; every two-phase SL study is clean.

**31 studies. PUBLISHED = a curated `<study>_errors.csv` already exists under a
versioned `docs/*/data/tables/`.**

| study | np | SL_SCHEME | published |
|---|---|---|---|
| bench2Dvortex | 12 | default(pointValue) | |
| bench3Dshear | 12 | default(pointValue) | |
| bench3Ddeformation | 12 | default(pointValue) | |
| linearConv2Dvortex | 4 | default(pointValue) | **yes** |
| linearConv2DvortexClip | 4 | default(pointValue) | **yes** |
| linearConv3Dshear | 8 | default(pointValue) | **yes** |
| linearConv3DshearPoly | 8 | default(pointValue) | **yes** |
| linearConv3Ddeformation | 8 | default(pointValue) | **yes** |
| linearConv3DdeformationPoly | 8 | default(pointValue) | **yes** |
| uncachedConv2Dvortex | 4 | default(pointValue) | **yes** |
| uncachedConv3Dshear | 8 | default(pointValue) | **yes** |
| uncachedConv3DshearPoly | 8 | default(pointValue) | **yes** |
| uncachedConv3DdeformationPoly | 8 | default(pointValue) | **yes** |
| uncachedConv3Ddeformation | 8 | default(pointValue) | |
| uncachedCompare2D | 4 | default(pointValue) | |
| nslConv2Dvortex | 4 | pointValue, normalProjected | **yes** |
| nslConv3Ddeformation | 4 | pointValue, normalProjected | **yes** |
| nslConv3Dshear | 8 | pointValue, normalProjected | |
| npslConv2Dvortex | 4 | normalProjected | **yes** |
| sdCompare2D | 4 | default(pointValue) | **yes** |
| idwCompare2D | 4 | default(pointValue) | |
| idwCompare3Dshear | 8 | default(pointValue) | |
| idwCompare3Dshear128 | 8 | default(pointValue) | |
| idwCompare3Ddeformation | 8 | default(pointValue) | |
| idwCompare3Ddeformation128 | 8 | default(pointValue) | |
| 3DshearSL | 8 | default(pointValue) | |
| 3DdeformationSL | 8 | default(pointValue) | |
| bulkVortexSL | 4 | default(pointValue) | |
| bulkVortexSLHighRes | 4 | default(pointValue) | |
| capillaryEnvelopeKinematic | 4 | default(pointValue) | |
| vortexSLsmoke | 4 | default(pointValue) | |

The transport ground-truth table in `docs/plan-curvature-stabilization.md:20-30`
(2D reversed vortex shape 2.97 / volume 2.98; 3D shear 2.50 / 3.25; 3D
deformation 1.36 / 1.86; poly 1.46 / 2.51) is drawn from this set.

### What is NOT affected

- Anything **serial**.
- The **two-phase** SL solver and every study using it (`U` from the momentum
  solve, no `velocityModel`).
- **Sourceless / non-gradient** consumers: advection uses `phi`, which was always
  correct.
- The `nSL` normal-projected *value* path except through `grad(U)` — but
  `normalProjectedScheme.C:161` does consume it, so those arms are affected too.

### How large is the effect?

Unknown for these cases, and I have deliberately not guessed. What is measured is
the 1D gate above: at np=8 the error in a gradient-dependent quantity reached
**4.3e-02** where serial was exact. Whether the SL foot's `dt²/2` term — which is
one term of a second-order trajectory, not the whole thing — moves a fitted
convergence order by more than the fit noise is exactly what a re-run answers and
nothing else does.

---

## 4. Bearing on `IMPROVEMENTS.md` item 2

Item 2 ("Semi-Lagrangian foot points: use the converged `u^{n+1}`; add a genuine
RK2 integrator") is unaffected in its diagnosis and strengthened in its remedy.

**Problem A is independently corroborated, with a magnitude.** The item observes
that with `PSI_OUTER_CORRECTORS no` the foot trace is fed the AB2 extrapolation
`u^n + (dt/dt0)(u^n − u^{n−1})` — old time levels only — so a second-order
trajectory receives a first-order velocity.

The same defect *class* was measured this week in an unrelated code path: the
`closestPoint` velocity extension builds `Uext` at the start of the step and
holds it, and on the 1D gate the resulting interface error is exactly first order
in dt — `err/CFL` constant at **≈ −0.5** across CFL 0.5 / 0.25 / 0.125 / 0.0625.
So the price of "velocity lagged by one step level" in this code base is about
**0.5·CFL cells of interface displacement**; a quarter of a cell at the standard
CFL 0.5. Two independent instances of one defect class is much stronger evidence
than either alone, and the extension measurement supplies the size that the SL
analysis could not.

**Fix B (RK2) has a better argument than the item makes for it.** The item
motivates RK2 partly on order, but the existing Taylor foot is *already* second
order in the trajectory, so RK2 adds none. Its real advantage is structural:

> the Taylor foot needs a velocity **gradient**; RK2 needs only velocity
> **samples along the trajectory**.

Given `fvc::grad(U)` was silently wrong at every processor boundary for the
entire life of the code, an integrator that does not consume a velocity gradient
is materially more robust, not merely equal-order. That is the case I would make
for Fix B.

**Sequencing.** Implementing Fix A or Fix B *before* the re-runs means measuring
a new scheme against contaminated references. The baseline should be
re-established first.

---

## 5. What this session is and is not doing

**Not touching**, by standing agreement: `applications/solvers/
leiaSemiLagrangianLevelSetTwoPhaseFoam/`, `src/leiaLevelSet/semiLagrangian/`,
`docs/semi-lagrangian-level-set/`, `docs/method-comparison/`, and the curvature
gates (`leiaTestSignedDistanceEllipsoid`, `faceCurvature*3D*`). The
`pointValueScheme.C` reference above is a **read** for diagnosis only.

**Decisions that are yours**, not mine:

1. Whether the 31 studies are re-run, and in what order. My own contaminated set
   (29 SDPLS studies) is already tiered P0–P3 and blocked on Lichtenberg; these
   would merge into the same tiering, and the 15 published ones belong in P0.
2. Whether `PSI_OUTER_CORRECTORS` is defaulted to `yes` (item 2 Fix A). It changes
   every SL study's numbers, so it should be decided **together with** the re-run
   so the studies move once rather than twice.
3. Whether Fix B (RK2 foot integrator) is worth building, and on which argument.
4. Item 2's "related config drift" (three sources disagreeing on the curvature
   delivery: template `none` / METHOD.md `psiOverGradPsi` / STATUS.md
   `cellCentreInverse`) — entirely a curvature question and entirely yours.

**One correction to the record I would make regardless.** The commit message of
`30e6ba9` states that the SDPLS source is the only consumer of `fvc::grad(U)`.
That is wrong: `closestPoint.C:432`, `steadyUpwindLinear.C:75`,
`normalProjectedScheme.C:161`,
`signedDistanceQuadraticWeightedLeastSquaresReconstruction.C:319` and — the one
that matters most, because it is on the default path —
`pointValueScheme.C:353` all consume it.
