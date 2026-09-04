# STATUS — where the work stands

Living hand-off file. Written to be usable from a phone: every command below is
meant to be run **on Lichtenberg**, and nothing here needs a local OpenFOAM.

Last updated: 2026-09-04 (static local refinement around the interface: gates G-1/G0 passed, WB gate running; polyhedral 3D rung R/h = 13.8 completed).

Conventions this file assumes are already known: [CLAUDE.md](CLAUDE.md) (layout,
build, git discipline) and [CLUSTER.md](CLUSTER.md) (full verified cluster
workflow). This file is the *current state*, not the manual.

---

## 0. READ THIS FIRST — the translating droplet case was a CLOSED BOX (2026-09-02)

**Everything this file previously said about the translating droplet is VOID.** The
whole of the old section 0 — "rhoLENT + a matched momentum time scheme makes the
translating droplet run", the pairing table, the best-configuration table, the
droplet-leaves-the-domain retraction, the `div(rhoPhi,U)` scheme comparison — was read
off `maxMagUPrime`, `t_blow` and the travelled fraction on a case that had **no inlet
and no outlet**.

`cases/translatingDroplet2D/system/blockMeshDict` put left, right, top and bottom into
a single `walls` patch. Every field in `0.org` has always declared `inlet` and
`outlet`; OpenFOAM errors when a **mesh** patch is missing from a field but **silently
ignores a field entry that matches no mesh patch**, so
`inlet { type fixedValue; value uniform (0.05 0 0); }` was parsed and discarded on
every run this campaign ever did on that case.

So it was a closed slip-wall box. Slip imposes `U.n = 0`, and on the x-normal faces the
normal *is* x, so the uniform stream the case initialises with is incompatible with its
own boundaries. The pressure projection annihilated it on step 1:

| | closed box | repaired |
|---|---|---|
| first-step local continuity error | 1.00e-05 | **2.32e-20** |
| `max\|U-U0\|` at step 1 | 1.031e-01 | **2.159e-03** |
| `mean\|U-U0\|` | 5.011e-02 | **1.160e-05** |
| `maxMagU` | 0.05987 | **0.05207** = U0 |
| ambient mean Ux at t = 0.02 | -2.04e-03 | at U0 |
| droplet mean Ux | +6.14e-02 | — |
| whole-domain mean Ux | -7.3e-07 | — |

A domain mean of zero is what a closed box must have. `maxMagUPrime = max|U - (U0,0,0)|`
was therefore measuring the **annihilated free stream** at roughly 2*U0 from the first
step onward — not a spurious current, on any arm, in any study.

**Fixed** in commit 440107f: `blockMeshDict` and its template now split `inlet` (left),
`outlet` (right), `walls` (top+bottom, still free-slip). The repair gate is
`config/translatingFreeStreamGate2D` (200 steps, N=128), whose pre-registered criteria
are the middle column above; it passes. On the repaired mesh `max|U-U0|` settles around
5e-03 and `mean|U-U0|` drifts 1.2e-05 -> 7.6e-05 over 200 steps — a real parasitic
current roughly 50x below the artefact that was hiding it.

**Voided, not reused** (CLAUDE.md, "A wrong setup voids its data"):

- cluster study dirs renamed `studies/*_VOID_closedBox_20260902`;
- eight curated tables moved to
  `docs/method-comparison/method-comparison-article/data/tables/VOID_closedBox_20260902/`
  with a README — `alphaFTest_{donorPlane,averagedPlanes}`, `bestConfigTranslating2D`,
  `rhoBoundGate2D`, `rhoDdtGate2D`, `rhoLENTGate2D`, `translatingMap2D`,
  `wellBalancedTranslating2D`;
- the equal-density control that was mid-flight was **stopped by job id** rather than
  allowed to finish, although its isolation argument is formally boundary-independent.
  A setup wrong in one way is not assumed wrong in only that way.

**Re-running now.** `config/translatingRepaired2D` and
`config/translatingRepairedEqualRho2D` — twin 8-arm matrices, `MASS_FLUX` (rhoLENT,
geometricFaceDensity) crossed with `MOMENTUM_DIV_SCHEME` (upwind, limitedLinearV 1,
vanLeerV, linearUpwind gradU), at density ratio 838.824 and at ratio 1 with the density
sum held at 999.39 so both sit at the same fraction of the capillary limit. 16 jobs,
ids 54460353-60 and 54460363-70, in `.my_jobs`.

### The source is the curvature estimator; translation and density ratio are amplifiers

`config/kickOriginGate2D` (8 arms, 200 steps, N=128, repaired mesh) decomposed the
first-step disturbance. Two results, both pre-registered:

**1. The step-1 kick does not depend on the translation speed.**

| U0 | `L1(U-U0)/U_ref` @ step 1 | @ step 200 |
|---|---|---|
| 0 | 2.3338e-04 | 1.3097e-03 |
| 0.0125 | 2.3305e-04 | 1.1620e-03 |
| 0.025 | 2.3271e-04 | 1.1752e-03 |
| 0.05 | 2.3198e-04 | 1.5238e-03 |

0.4% variation across a 4x change in U0. The mass-momentum consistency term
`U0 * [ddt(rho) + div(rhoPhi)]` is identically zero at U0 = 0 and would scale linearly,
so the kick is neither a translation nor a consistency effect. **The stationary droplet
carries the same kick.** The 1e-5 long quoted for it is its SETTLED level; comparing that
against a first-step transient was never a valid comparison.

(U_ref = 0.05 m/s throughout, so the U0 = 0 row is on the same scale.) The one-step
disturbance is already FOUR PERCENT of the translation speed.

**2. Exact curvature removes it by six orders.** With `constantCurvatureSurfaceTension`
at kappa = 1/R = 1000 the force `sigma*kappa*snGrad(alpha)` is an exact discrete
gradient, the pressure equation absorbs it entirely, and the kick falls
**2.15e-03 -> 1.69e-09**, a factor of 1.27e+06. The curvature error is
`kErrL2Band = 70.6` against kappa = 1000, i.e. **7%** at R/h = 12.8, identical in all
eight arms because it is a t = 0 property of the initial circle.

**The decomposition.** In `max|U|(T) = u0(h) * exp(G(h))`:

- **u0, the source** = curvature-estimator error. U0-independent (L1 2.3338e-04 against
  2.3198e-04 relative to U_ref, i.e. 0.6%, across a fourfold change in U0) but NOT density-ratio independent: at ratio 1
  with matched dynamic viscosity the same first step gives L1 5.769e-06, a factor of 40.5
  smaller. RETRACTED 2026-09-03 -- "density-ratio independent" was asserted from a matrix
  that only ever varied U0 at fixed ratio 838.8, and amplifierGate2D falsifies it. The
  same curvature error divided by the light-phase density is a larger velocity, so the
  density ratio acts on BOTH factors.
- **G, the amplifier** = grows with U0 (over the full horizon the U0 = 0 and
  U0 = 0.05 arms separate by a factor of 67 in L1) and with the density ratio (the repaired 16-arm matrix: both ratios sit at
  the same order for 8000 steps, then only the ratio-838.8 arms go unstable).

Translation does not create the parasitic current. It stops it from relaxing.

Corollary for the matrix: read at step 5000 (t = 0.0375), where every arm is healthy,
`L1(U-U0)/U_ref` is 7.8e-03..9.8e-03 at ratio 838.8 against 5.7e-03..6.3e-03 at ratio 1 --
a factor of 1.4-1.6, NOT the 18x the first curation reported (which was L_inf, read inside
the blow-ups). That table was read at the
common step 8427, which was set by the first arm to die, so every ratio-838.8 row was
sampled inside its own blow-up. Corrected; the light-phase amplification argument that
predicted a factor of 420 is dead.

Open confounder: the equal-density control holds the KINEMATIC viscosities fixed, so the
dynamic viscosity ratio moved 54.8 -> 15.3. A matched-mu control is still owed.

### ANCHORED: the translating droplet reproduces Popinet's own benchmark (2026-09-04)

The translating droplet is Popinet's test case, introduced in JCP 228 (2009) 5838-5866
Sec. 6.2.2 for exactly this coupling -- he notes he was "not aware of any study" combining a
stationary droplet with surface tension and the advection of a circular interface.
`cases/popinetTranslating2D` + `config/popinet2D_*` reproduce it as a standalone study.

HIS SETUP, read from the PDF in the repo root rather than recalled: D = 0.4 in the unit
square at 64x64 (R/h = 12.8), periodic in x with symmetry top/bottom, We = 0.4, and
La = sigma D/(rho nu^2) over {120, 1200, 12000, infinity}. **Both ratios are 1** -- Sec 6.2.1
says "rho the constant density" and La carries a single rho and nu. His air/water ratios
(850, 55.72) appear ONLY in the Sec 6.3 capillary wave. His face properties use the
volume-fraction mixture with a face fraction that is "a simple average of the cell-centred
values", i.e. our `alg_lin`, which our own 36-arm ladder independently chose.

RESULT, in his convention (maximum over TIME of the spatial norm, velocity relative to U):

| N | R/h | L2 (his RMS) | Linf (his Max) | shape |
|---|---|---|---|---|
| 64 | 12.8 | 4.4e-03 | **4.86e-02** | 8.15e-04 |
| 128 | 25.6 | 2.5e-03 | 3.36e-02 | 2.08e-04 |
| 256 | 51.2 | 1.3e-03 | 2.47e-02 | 7.57e-05 |
| order p (R) | | **0.88 (0.999)** | **0.49 (0.999)** | 1.71 (0.996) |

At his own resolution we get **4.86% of U against his "of the order of 5%"**, and RMS
4.4e-03 against his figure's ~3e-03 peak. Both of his convergence statements reproduce:
"less than first-order for the maximum error and close to first-order for the RMS" -- we
measure 0.49 and 0.88, each with R = 0.999. Shape converges at 1.71, better than his
"roughly first order".

**This is the empirical basis for not reporting L_inf.** The half-order maximum against
near-first-order RMS reproduces across two unrelated discretisations, so it is a property of
the benchmark, not of our scheme.

Laplace sweep at N=64, Linf/U: 2.93e-02 (La=120), 3.96e-02 (1200), 4.86e-02 (12000),
5.88e-02 (inviscid). A factor of two across two decades of La plus the inviscid limit,
monotone -- viscosity DAMPS here, and the inviscid case is worst. He calls that "weakly
dependent on the Laplace number"; the factor is worth quoting rather than the adjective.
Note this is the constant-property case; with a viscosity jump the balance differs (the
viscosity gate finds removal of viscosity lowers L1 but raises L2).

Deviations that cannot be removed, stated in the config: no cyclic support in the SL
transport, so inlet/outlet replaces periodicity and the box is 2x1 to leave the droplet 57
cells clear of the outlet over a full t/T_U; and his adaptive quadtree VOF with
height-function curvature against this uniform-mesh level set, which is the point of the
comparison. Symmetry top/bottom is exact -- slip imposes U.n = 0.

Curated by `workflow/scripts/make_popinet_table.py` into
`docs/semi-lagrangian-level-set/sl-level-set-article/data/tables/popinetTranslating.csv`;
written up in the article as "Comparison with the reference translating-droplet benchmark".

### DECIDED: the instability does not survive a perfect capillary force

`config/amplifierGate{,EqualRho}2D`, 8 arms, 8000 steps, N=128, 2x2x2 over
{reconstructed, exact kappa=1/R} x {U0 = 0, 0.05} x {ratio 838.8, ratio 1 at MATCHED
dynamic viscosity}. All complete. Velocity errors are the volume-weighted L1 norm RELATIVE to the translation
speed, `L1(U-U0)/U_ref` with U_ref = 0.05 m/s held fixed so the U0 = 0 control sits on the
same scale as the rest of its column. **L_inf is not reported**: on a parasitic current the
maximum is a single-cell extremum set by the sub-cell interface position, it does not
converge under refinement, and a verdict built on it had to be retracted (2026-09-03).

| ratio | kappa | U0 | step 1 | step 1000 | step 5000 | step 8000 |
|---|---|---|---|---|---|---|
| 838.8 | exact | 0 | 1.08e-08 | 3.54e-11 | 1.32e-10 | 1.46e-11 |
| 838.8 | exact | 0.05 | 4.15e-09 | 5.09e-06 | 2.04e-05 | 3.00e-05 |
| 838.8 | reconstructed | 0 | 2.33e-04 | 1.77e-04 | 6.96e-05 | 4.45e-04 |
| 838.8 | reconstructed | 0.05 | 2.32e-04 | 4.32e-03 | 7.85e-03 | **2.98e-02** |
| 1 | exact | 0 | 4.13e-11 | 1.33e-13 | 1.54e-13 | 1.39e-12 |
| 1 | exact | 0.05 | 7.04e-11 | 4.74e-11 | 6.04e-11 | **7.77e-11** |
| 1 | reconstructed | 0 | 5.77e-06 | 3.05e-04 | 9.63e-05 | 1.86e-05 |
| 1 | reconstructed | 0.05 | 5.76e-06 | 2.04e-03 | 6.39e-03 | 5.72e-03 |

**The exact-kappa arms stay bounded at every combination.** The growth seen over 200 steps
in kickOriginGate2D was a transient; over the full horizon it plateaus at 3e-05 of the
translation speed at ratio 838.8 and sits at round-off, flat, at ratio 1. The
reconstructed-kappa arm at U0 = 0.05 reaches 3.0e-02 of U0, a factor of 1000 above it. The reconstructed-kappa arms at
U0 = 0.05 are the ones that destabilise (8427-9987 steps in the repaired matrix).

**Verdict: the curvature estimator is the whole story, and the variational capillary force
(article Section "Outlook: a variational capillary force") is the right target.** No
amplifier independent of the source needs attacking first.

With exact kappa the volume error (4.32e-03), shape error (6.9e-06) and travelled fraction
(0.9975) are IDENTICAL to three digits at both density ratios -- that residual is pure
semi-Lagrangian transport error, which cannot see rho. And the net propulsion is gone:
travel 0.9975 instead of 1.0257. The drift the animation shows is a curvature-error product.

The matched-mu control is now also closed: NU1/NU2 rescaled so mu1 = 9.982e-04 and
mu2 = 1.8207e-05 match the ratio-838.8 arms exactly, so only rho differs.

### What SURVIVES the retraction

These did not run on `translatingDroplet2D` and are unaffected:

- **`rhoLENT` is safe on the stationary droplet.** `config/rhoLENTStationary2D`, all 6
  arms complete, control reproducing the shared ladder to five digits; neutral-to-better
  (+1.0% / -22% / +0.1% at N=32/64/128), volume and shape matching to three digits.
  `rhoLENT` remains the shared `MASS_FLUX` default.
- **`geometricFaceDensity` carries a mass residual of 0.04-0.56 relative on the
  STATIONARY droplet** and does no harm there. That measurement holds the residual fixed
  and varies only `U0`, and is the reason to expect a `U0 x [ddt(rho) + div(rhoPhi)]`
  momentum source at all. It is now the *motivation* for the repaired matrix rather than
  a result about it.
- **The 2D stationary ladders on the shared configuration.** Unabsorbed capillary
  residual over the second half: `cellCentreInverse` wins at every rung, 4.60x / 3.81x /
  1.55x / 1.49x at N=32/64/128/256.
- **The kinematic transport gates** (prescribed velocity, no momentum solve, different
  case) — unaffected by pressure boundary conditions.

### The curvature chain is exonerated, four times over

`kinematicTranslation2D` (transport 2nd–3rd order, bounded-α error exactly 0);
`wellBalancedTranslating2D` (**exact constant curvature still diverged**);
`bestConfigTranslating2D` (the required extension does not help, hurts coarse rungs);
`rhoDdtGate2D` (density ddt is not the mechanism).

### Leading open defect

The droplet **inflates** monotonically, +5% to +8% over the horizon, and at matched
times the gain does **not** converge under refinement. Untouched by the mass-flux
model — α comes from ψ, so this is semi-Lagrangian transport under the computed
velocity field.

### New, runtime-selectable, default-off (nothing existing moves)

`massFlux { boundRho true; }` (clip ρ to `[rho2, rho1]`; the clip is **reported**),
`massFlux { massResidualDiagnostic true; }` (mass residual for every model),
`DROPLET_OFFSET_X` (droplet start placement). All three gated for bit-identity and on
4 MPI ranks before any cluster use.

---

## 1. What is being worked on

**Since 2026-09-04: static local mesh refinement around the interface.** The 3D
ladders are expensive because the whole 6R box carries the interface's cell size
(572k / 1.46M / 3.6M polyhedra at R/h = 13.8 / 18.9 / 25.6). The method reads and
writes the interface only in a narrow band, so the band is refined and the far field
stays coarse -- as a PRE-PROCESSING step, the solver untouched:
`workflow/scripts/leiaRefineHexMesh.py` (blockMesh + N passes of [0/ := 0.org;
leiaSetFields; topoSet seed 0 < alpha < 1 + face dilations; refineHexMesh], then
leiaSetFields on the final mesh) and `leiaRefinePolyMesh.py` (cfMesh re-meshed around
the psi = 0 iso-surface). New mesh kinds `hexRefined | polyRefined`, tokens
`REFINE_LEVELS / REFINE_BAND_CELLS / REFINE_SOURCE`, gates and configs in the
"RUNNING: static local refinement" block of section 4. The question the campaign
below left open -- the parasitic-current source and its amplifier -- is unchanged;
refinement is what makes R/h >= 20 in 3D affordable, on hex first and then polyhedra.

**Question:** where do the parasitic currents in the stationary-droplet test come
from, in the combination `leiaSemiLagrangianLevelSetTwoPhaseFoam` (quadratic
semi-Lagrangian level set) + balanced CSF surface tension?

**Where it stands.** Two hypotheses have been tested and one of them is dead:

1. *The spurious variation of the delivered face curvature drives it* —
   **falsified**. A delivery that removes the across-support variation of
   `kappa_f` (one inverted value per cut cell) made the instability **faster**,
   not slower. Static accuracy turns out to be **anti-correlated** with coupled
   stability across every delivery measured.
2. *The gain drives it* — the gain being
   `G = ||d kappa_f|| / ||d psi||` on the CSF force support, i.e. how far the
   delivered face curvature moves per unit perturbation of `psi`, reported as
   `G h^2`. This is a **strong but incomplete** predictor: it orders the extremes
   correctly and the two deliveries with equal gain have equal growth rate, but it
   does not order every pair.

The gain is measurable offline in seconds on a static circle, which is the main
practical gain from this thread: a candidate delivery can now be rejected before
any coupled run.

3. *A varying-curvature gate was missing* — **now built, and it retracted a
   decision.** On a 2:1 signed-distance ellipse both one-value-per-cut-cell
   deliveries **collapse to first order** (1.02 and 1.03 against 1.98 for the
   per-face inversion), while on the circle they looked like 2.00 and 2.05. Forcing
   `kappa_f` constant across a cut cell gives each face a value centred on the cell
   rather than on itself — an O(h) offset wherever curvature varies along the
   interface. Structural; better averaging cannot fix it.

**Current production delivery, and the only one that passes:**
`curvatureExtension stabilizedFootPointFace` — the per-face parallel-surface
inversion. Order 1.98 on the varying-curvature gate, `G h^2` = 0.647.

**Retired:** `cutCellFootPointFace` and `cellMeanFootPointFace`. The cell-mean arm
had the lowest gain and the longest coupled survival (t_blow 0.1049 s), but that was
measured on a **circle**, whose constant curvature hides the first-order defect
exactly as the circle static gate did. Kept in the code and in both gates as the
two extreme points that established the gain-versus-accuracy dissociation.

---

## 2. Where the writing lives

| what | path |
|---|---|
| The program document, all measured results, section by section | [docs/plan-curvature-stabilization.md](docs/plan-curvature-stabilization.md) |
| ... section 9 | WP8.1 — the across-support driver split |
| ... section 10 | WP8.2 — cut-cell delivery falsifies it; gain is the predictor |
| ... section 11 | WP8.3 — what the inverse assumes; the residual law; the cell-mean delivery |
| Slide deck (derivation slides + verdict slide) | [docs/semi-lagrangian-level-set/sl-level-set-presentation/quadratic-semi-lagrangian-level-set.template.html](docs/semi-lagrangian-level-set/sl-level-set-presentation/quadratic-semi-lagrangian-level-set.template.html) |
| Article (flux-space residual section) | [docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex](docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex) |

Decks are built with `bash docs/build-decks.sh`; only the `*.template.html` is
tracked, the built `*.html` is not.

## 3. Where the numbers live

All curated, all committed, all regenerated by committed scripts:

| what | path |
|---|---|
| Face-curvature convergence, 3 gates | `docs/method-comparison/method-comparison-article/data/tables/face_curvature_orders{,_3d,_ellipse}.{csv,tex}` |
| Gain `G h^2`, 3 gates | `docs/method-comparison/method-comparison-article/data/tables/curvature_gain{,_3d,_ellipse}.{csv,tex}` |
| Convergence figures | `docs/method-comparison/method-comparison-article/data/figures/face_curvature_convergence{,_3d,_ellipse}.png` |

The three gates are the circle (`faceCurvatureDroplet2D`), the sphere
(`faceCurvatureSphere3D`) and the **2:1 signed-distance ellipse**
(`faceCurvatureEllipse2D`) — the last is the only one with non-constant exact
curvature, and it is the one that decides accuracy.

Regenerate with:

```
python3 workflow/scripts/make_face_curvature_fig.py studies/faceCurvatureDroplet2D
python3 workflow/scripts/make_curvature_gain_table.py studies/faceCurvatureDroplet2D
```

(and the same two with `studies/faceCurvatureSphere3D` for the 3D artifacts).

## 4. State of the measurements

**The gate that decides accuracy** — 2:1 ellipse, curvature 250–2000 1/m,
order fitted on N >= 128:

| delivery | order | L2 at N=512 [1/m] | `G h^2` | order on the circle |
|---|---|---|---|---|
| arithmetic | 0.97 | 13.66 | 0.647 | 1.13 |
| **per-face inverse** | **1.98** | **0.2785** | 0.647 | 2.04 |
| cut-cell inverse | 1.02 | 5.848 | 0.818 | 2.00 |
| cell-mean inverse | 1.03 | 5.851 | 0.395 | 2.05 |
| interface mean (control) | −0.01 | 506 | — | 2.03 |

**The foliation gate** (2026-08-14, plan sec. 17,
`config/faceCurvatureEllipsoidPsi2D.yaml`): same ellipse, psi as TRUE signed
distance (fit error only) vs the QUADRATIC form (fit exact, inverse's foliation
bias only; beta varies 2x along the interface). CORRECTED 2026-08-14: the gate
row `stableFootPoint` is a GATE-ONLY foot-point-native variant, NOT the solver
delivery; the new `solverStabFootFace` row calls the shipped
stabilizedFootPointFace header and is bit-equal to the foot-corrected plain
quadratic — so PRODUCTION IS SECOND ORDER while the foliation is parallel
(0.278, order 1.98 on this ellipse; 2.04 circle, ~1.95 sphere) and first order
only where the foliation is non-parallel (8.69, order 0.94: the d*D bias); and
`quadraticNewtonFoot` INVERTS — 2nd order, 0.187 at N=512, when the fit is
exact, because evaluating AT the interface point needs no parallel-surface
correction at all. Curated:
`docs/method-comparison/.../tables/face_curvature_orders_foliation.csv`.

**The foot-point-EVALUATED delivery** (2026-08-15,
`curvatureExtension footPointEvaluatedFace`, plan sec. 17.4): each adjacent cell
finds the foot point of the FACE CENTRE on its own quadratic's zero set and
evaluates that fit's curvature there; the face takes the linear interpolation of
the two. No parallel-surface conversion, so no foliation assumption. It is the
exact complement of production — where one is second order the other is first:

| psi (geometry) | production | footPointEvaluatedFace |
|---|---|---|
| circle, exact SDF (N=512) | 0.105, order 2.0 | 10.5, order ~1.0 |
| 2:1 ellipse, SDF (N=512) | 0.278, order 1.98 | 12.7, order ~1.0 |
| 2:1 ellipse, quadratic form (N=512) | 8.69, order 0.94 | **4.8e-4 (exact)** |
| oscillating geom, a/b=1.21 (N=256) | 6.24, order ~1.07 | **3.3e-4 (exact)** |

`G h^2` is EQUAL to production within 3% (0.657 vs 0.642, circle at N=512), so a
coupled difference is an accuracy difference, not an amplification one.

**COUPLED VERDICT (2026-08-15): the prediction is FALSIFIED and the delivery does
NOT promote.** On the OSCILLATING case -- the discriminating one, where its static
curvature is EXACT (3.3e-4 vs 6.24, four orders better than production) -- the
foot-evaluated delivery blows up EARLIER at every resolution:

| N | t_blow production | t_blow footEval | ratio | volume error at common t |
|---|---|---|---|---|
| 64 | 0.02132 | 0.00772 | **0.36** | 2.41e-2 vs 3.36e-3 |
| 128 | 0.01454 | 0.00613 | **0.42** | 1.28e-2 vs 2.51e-3 |
| 256 | 0.01258 | 0.00398 | **0.32** | 6.65e-3 vs 4.36e-4 |

Volume error is 5-15x worse at every N; SHAPE error is comparable and at N=128
slightly better (2.59e-5 vs 4.36e-5), so the two disagree -- read together, never
singly. On the stationary CONTROL there is no consistent ordering at all
(t_blow ratio 0.95 / 1.12 / 0.50 at N=64/128/256).

This is the THIRD delivery for which static accuracy is anti-correlated with
coupled stability, and the second pair the gain `G h^2` fails to order (equal
within 3%, blow-up 3x apart). Curated:
`docs/method-comparison/.../tables/foot_evaluated_face_coupled.csv`.

NOTE ON SCORING: r(A2h) was the intended order parameter, but these runs blow up
within 400-9400 steps, so its fit window is too short to be stable (ratios swing
from 0.07 to 16.9 on the control). t_blow is the robust discriminator here --
monotone across all three resolutions on the test case -- and is what the verdict
rests on.

**The time/coupling axis** (2026-08-13/14, plan sec. 16): r(N,dt) = r0 + c·dt
measured at N = 64/128/256; GAMG, decomposition (seed only), adaptivity (now
disabled by construction — fixed capillary deltaT is the case standard), the
frozen-force lag (psiOuterCorrectors, at any Picard depth) and the x_d kernel
(leiaTestDeparturePoint: order 3, correct constants at omega·dt = 0.52) are all
exonerated by direct experiment. Mode-resolved: r(maxU) ≈ 2·r(A2h);
the corrugation rate r(A2h) is the order parameter and is nearly
dt-independent at N=256.

**THE COMBINATION WORKS (2026-08-18, 28 coupled arms). Curated:
`docs/method-comparison/.../tables/capillary_envelope_coupled.csv`.**

*cellCentreInverse + psiFilter biharmonicBand theta=0.2, 2D*: all of
N=64/128/256 reach t=0.1 and EVERY metric improves monotonically under
refinement --

    N        max|U|      volume      shape       min|grad psi|
    64       2.24e-3     4.96e-3     3.78e-5     0.9520
    128      4.15e-4     2.57e-5     5.24e-6     0.9955
    256      1.02e-4     7.33e-6     8.14e-8     0.9996
    order    +2.43,+2.02 +7.6,+1.8   +2.9,+6.0

This is the first configuration in the campaign in which a stationary droplet
with a LIVE surface-tension force survives at every resolution AND converges in
current, volume and shape simultaneously. At N=256 it beats production+filter
1.3x on current, 13x on VOLUME and 8.5x on SHAPE, with the profile essentially
undrifted (0.9996 vs 0.9931). Production+filter at the same theta is NOT
monotone (1.94e-4 -> 1.71e-3 -> 1.34e-4): it survives but does not converge.

*The theta=0.05 arm does not converge* for either delivery (production
5.6e-4/5.4e-3/3.4e-2; inverse 1.3e-2/4.9e-4/9.3e-3), so theta=0.2 is not merely
"more damping is better" -- theta=0.05 is below what the loop needs.

*3D, N=64 works and N=32 does not.* theta=0.2 at N=64 reaches t=0.1 (max|U|
5.89e-3, volume 1.77e-4, shape 4.41e-5, min|grad psi| 0.964) while BOTH filtered
N=32 arms die EARLIER than the unfiltered N=32 (0.0546/0.0548 vs 0.0613) -- and
the theta=0.2 one dies hardest (max|U| 22.5 vs 0.244). At N=32 in 3D the drop is
only R/h = 3.2, so the filter's band+1-ring support spans a large fraction of it
and its grid-scale injection is amplified rather than damped: the failure mode
recorded in the config header before the run. Working requirement so far:
**R/h >~ 6**.

*THE TRANSPORT-ORDER AXIS IS FLAT, which confirms the temporal-cap argument.*
Nine matched (N, theta) pairs of quadratic vs linearTaylor TRANSPORT with the
curvature fit pinned quadratic (GEOMETRY_FIT) agree to 3-4 significant digits in
every FILTERED arm -- e.g. N=256 theta=0.2: max|U| 1.335e-4 vs 1.341e-4, volume
9.65e-5 vs 9.73e-5, shape 6.94e-7 vs 6.99e-7. Only the unfiltered arms differ
at all, and mildly/inconsistently (N=128 linearTaylor dies later, N=256 earlier
-- within the 5-38% seed spread). This is the predicted consequence of the
coupled interface motion being capped at FIRST order in time by the Euler
momentum solve on a dt ~ h^1.5 capillary step (temporal error ~ h^1.5 dominating
the O(h^2) spatial error): in the coupled problem the transport order is not the
limiting factor. Combined with the kinematic arm, where linearTaylor's errors
GROW under refinement (volume 340% at N=256 on the reversed vortex — **number
suspect, see the 2026-08-26 gradU contamination notice below; np=4 kinematic
run**), the verdict is that dropping transport order buys NOTHING coupled and
costs everything kinematic -- the filter is the damping mechanism, not the
transport order.

**THE CELL-CENTRE INVERSE AND THE FILTER (2026-08-18).** Two constructions the
user proposed now hold the best coupled results on record; curated:
`docs/method-comparison/.../tables/cell_centre_inverse_coupled.csv`.

*curvatureExtension cellCentreInverse* -- the parallel-surface inverse applied
per CELL (kappa_c, K_c at the centre, d_c = signedOffset, same dimension-
seamless algebra), the transformed cell field delivered by plain interpolation.
Its non-gradient force content converges at order +2.01 over the full CSF
support where every other cell curvature is +0.09 (leiaTestRemainderTerm.csv).
Coupled, UNFILTERED: at N=64 the lowest parasitic floor ever measured with a
live force (max|U| 8.05e-6, volume 7.2e-6, shape 7.2e-7 at t=0.1); outlives
production at every N in 2D (0.100/0.078/0.036 vs 0.100/0.063/0.033); in 3D it
beats production 2.6x at N=64 AND inverts the finer-blows-first trend between
N=32 and N=64 (0.0613 -> 0.0786; two points, seed sensitivity 5-38%). In 2D the
fine meshes still blow, and at N=256 plain arithmetic outlives it (0.049 vs
0.036): the floor improved, the growth survived.

*psiFilter biharmonicBand theta=0.2, production delivery* (capillaryEnvelope,
quadratic half): ALL of N=64/128/256 reach t=0.1, and at N=256 the run ends
healthier than N=64 (max|U| 1.3e-4, volume 9.7e-5, shape 6.9e-7, min|grad psi|
0.993) -- the h-trend inverted with a live force. theta=0.05 survives everywhere
but the residual current grows ~10x per refinement level.

*The kinematic price arm settles the linear-transport idea*: on the reversed 2D
vortex linearTaylor's errors GROW under refinement (volume 340% at N=256, shape
diverging) -- its interpolation diffusion is only useful where the flow is
near-quiescent, so the theta filter is the damping mechanism of choice, not the
transport order. (**Numbers suspect** — see the 2026-08-26 gradU contamination
notice below. The comparison is partly common-mode, both arms sharing the same
contaminated foot, so the verdict may survive; the numbers are not defensible
until re-run.)

### Mesh alignment is EXONERATED for the R/h = 25 departure (2026-08-26)

The 2x2 probe pair designed for the question "is R/h = 25 a degenerate
pole-on-face alignment?" completed (`studies/alignmentProbe2D/`, quarter horizon
t = 0.025 s, filter off, np 1):

| arm | geometry | maxU 17 ms | maxU end | r [1/s], 17-25 ms | vol end | shape end |
|---|---|---|---|---|---|---|
| rh25 original | N=150 centred, poles ON faces | 3.99e-4 | 2.59e-2 | 598 | 5.55e-5 | 1.65e-6 |
| rh26_centered | N=156 centred, SAME class | 4.93e-4 | 3.21e-2 | 550 | 2.66e-5 | 4.32e-6 |
| rh25_shifted | N=150, centre +(h/3, h/7), poles OFF faces | 1.42e-2 | 1.56e-1 | 308 | 3.85e-4 | 4.82e-5 |

**Both probes depart, so per the pre-registered read-out the departure is
fineness, not alignment.** The same-class control (N = 156) reproduces the
original's trajectory (departure ~17 ms, r 550 vs 598), so nothing is special
about 25; and the shifted arm — the one that HAD to stay mild for alignment to
be load-bearing — is strictly WORSE at every sampled time (50x higher current
already at 11 ms, 7x worse volume, 29x worse shape at the end): breaking the
symmetry removes error cancellation instead of removing the pump. The
symmetric, pole-aligned setup was the FAVOURABLE configuration. Consequences:
no bug hides in the alignment degeneracy; the "white points" are not tied to
pole placement; the fix search stays on the amplifier itself (semi-implicit
capillary coupling / variational pairing), not on mesh positioning. Caveat:
both probes carry two mid-run restarts (WSL reboots; identical handling in both
arms, and the same-class control shows no restart inflation, so the comparative
reading stands); metrics CSVs cover 10-25 ms, joined records preserved as
`leiaSemiLagrangianLevelSetTwoPhaseFoam.joined.csv` per case.

### gradU contamination CLOSED for the 2D vortex arm: published orders stand (2026-08-27)

Forensic chain on `uncachedConv2Dvortex` N=256 CFL 0.5 (raw July run preserved in
`studies/uncachedConv2Dvortex.preGradUfix-20260826`, all new runs on pinned
binaries, `studies/gradUserialCheck/`):

1. **Bug effect on the endpoint: <= 2x, no order change.** July-10 code
   (`230b1dc`, isolated worktree build) run SERIAL — bug-free by construction —
   on the exact July inputs gives shape 1.105e-06 / volume 3.95e-05 at t=2,
   against the curated buggy-parallel 7.851e-07 / 8.01e-05: the seam bias moved
   shape 1.4x (flatteringly) and volume 2x. The published second-order
   transport table is a method property, not a bug artifact. (One arm, np=4;
   np=8+ 3D arms have more seam area — re-runs stay worthwhile, not urgent.)
2. **The transport has not changed since July: psi at t=2 is BIT-IDENTICAL**
   (max|diff| = 0.0 over 65k cells) between the July-10 build and today's, same
   inputs, serial.
3. **Everything alarming in the first re-run was the measurement layer:**
   (a) LOGGING REGRESSION — the solver now writes metric rows only at write
   times and never logs the final step, so the aggregation scored t=1.99970
   instead of t=2; the reversal cancellation collapses shape 20x in that last
   dt, and the missing gap ~ dt ~ h manufactures a fake order ~1. Also
   `L_INF_E_PSI` now reads 0 (metric dead) and the errors-CSV `method` column
   echoes the inert Eulerian `ADVECTION` dict entry for SL runs.
   (b) RETRACTED (2026-08-27, same day): the "indicator change" was the SAME
   gate. In the gated run alpha and the narrow band were last computed at
   t = 1.9997 with a band FROZEN at t ~ 1.0018 (maximum deformation), and the
   endTime field write dumped that stale alpha as 2/alpha — the wrong-side 0/1
   snaps in ~300 cut cells were stale-band artifacts, not algorithm changes.
   With the fixed writer (per-step alpha + band), HEAD's alpha(t=2) is
   referee-exact against analytic circle fractions (mean 2.3e-3, max 1.05e-2,
   zero cells off by >0.1 — indistinguishable from the July build), and every
   logged metric row matches the July-code serial run TO ALL DIGITS (t=1.9997:
   3.6890740384e-05 / 1.66685002663e-05 in both). The detrixheAslam indicator,
   922eb41 and a903db4 are all exonerated (1-step referee runs at both commits:
   mean 2.1e-3, max 4.5e-3). The coupled two-phase campaign was NEVER exposed —
   its solver computes alpha and the band every step in its own loop.
   nDefCorr 2 vs 3 and the grad(psi) macro: bit-identical, exonerated.
   One true kinematic-solver defect chain remains fixed in 84906ee: the gate.

### INVALIDATION: parallel kinematic SL baselines predate the gradU coupled-patch fix (2026-08-26)

`velocityModel::setVelocity` wrote FACE-centre values into coupled (processor)
patches — where every coupled-patch operator expects NEIGHBOUR-CELL values —
biasing `fvc::grad(U)` by O(1) in every processor-adjacent cell, for the life of
the code (fixed `30e6ba9`). The default SL foot integrator consumes that
gradient in its dt²/2 term (`pointValueScheme.C:353`), and the kinematic solver
`leiaSemiLagrangeLevelSetFoam` feeds it the prescribed velocity — so **all 31
parallel kinematic SL studies are contaminated, 10 of them with git-tracked
curated tables** (uncachedConv*, sdCompare2D under the SL theme; six linearConv*
under the linear-SL theme). The transport ground-truth table in
`docs/plan-curvature-stabilization.md` sec. 0 is drawn from this set and carries
its own notice. **NOT affected:** anything serial, and every two-phase SL study —
`leiaSemiLagrangianLevelSetTwoPhaseFoam` has no `velocityModel` reference and
takes U from the momentum solve, so the entire coupled parasitic-current
campaign above stands. Full analysis, affected-study list and the magnitude
measurement (1D gate: 4.3e-2 at np=8 where serial was exact):
`docs/gradU-coupled-patch-contamination.md`. Sequencing agreed 2026-08-26: first
a discriminating re-run of `uncachedConv2Dvortex` on the fixed binary diffed
against the tracked table; the RK2 foot integrator (gradient-free by
construction) and any `PSI_OUTER_CORRECTORS` default change wait until the
baseline is re-established, so studies move once.

*Regression found and fixed en route*: since 4c4f7f1 the shared createFields.H
constructed a standalone TRANSPORT-type reconstruction and called
meanCurvatureClosestPoint on it whenever the force is reconstructedCurvature --
an instant NotImplemented abort for every non-quadratic transport (all nine
linearTaylor envelope arms died at startup). Guarded on `kappa` not yet being
registered (2625a19); Eulerian solver byte-identical.

IN FLIGHT: capillaryEnvelope rerun on the fixed binary (54249071; its quadratic
arms double as the regression check) and the COMBINATION
cellCentreInverse x theta={0.05,0.2} in 2D (54249072) and 3D (54249073).

**THE SOURCE, MEASURED (2026-08-17, plan sec. 18).** New spectral harnesses
(`interface_mode_trajectory.py`, `mode_rate_vs_drift.py`) perturb the interface
by cos(m·theta) and fit growth + frequency per mode on the shipped solver, U/p
reset, validated by omega matching the capillary dispersion within 1.5%.
Result: the instability is the **m=2 capillary mode**, damped at N=64
(−11 1/s), unstable at N=128 already on a NEAR-CLEAN profile (+19 1/s at
minGradPsiBand 0.994) — that is the N=64→128 switch. Profile drift AMPLIFIES
(×3 by minGrad 0.71) but does not trigger; m=6 is damped everywhere. The dt
series pins the nature: r_2 = 18.8/12.7/8.0/4.0 at dt/1/2/4/8, Richardson
intercept **r0 = +0.03 1/s ≈ 0** (clean state; +2.1 upper bound on the drifted
one, where the coupling coefficient c rises 4.3e6 → 6.6e6). The semi-discrete
loop is NEUTRAL — the raw delta-psi map is exactly gauge-neutral (power
iteration: lambda = 1) — and the instability is the **explicit capillary
coupling**, rate c·dt. Supersedes sec. 16.1's finite-r0 reading (that fit
integrated over growing drift). Consequences: the semi-implicit capillary
force is un-shelved as the principal lever; SDF maintenance is the c-reducer
with a measured budget; delivery work is floor-only (the production face error
is one function of face-normal/interface-normal misalignment — the m=4 pole
bias — largest exactly on the highest-|snGrad(alpha)| faces). Curated:
`docs/method-comparison/.../tables/mode_rate_vs_drift.csv`,
`mode_rate_dt_series.csv`.

Constant-curvature gates — exact circle (N=32–512) and exact sphere (N=32–128).
These cannot rank accuracy; they are kept for the gain and for regression:

| delivery | 2D `G h^2` | 3D `G h^2` | 2D L2 (N=512) | 3D L2 (N=128) |
|---|---|---|---|---|
| arithmetic | 0.638 | 0.823 | 11.35 | 0.4191 |
| per-face inverse | 0.639 | 0.826 | 0.1049 | 0.005254 |
| cut-cell inverse | 0.829 | 1.005 | 0.07608 | 0.005176 |
| **cell-mean inverse** | **0.402** | **0.534** | 0.08422 | **0.004093** |

Coupled stationary droplet, np 8, capillary time step:

| arm | max\|U\| at t=0 | t_blow [s] |
|---|---|---|
| arithmetic N=128 | 2.31e-3 | 0.0803 |
| per-face N=128 | 2.83e-4 | 0.0668 |
| cut-cell N=128 | 1.48e-5 | 0.0202 |
| cell-mean N=128 | 1.13e-5 | **0.1049 — longest of the four** |
| per-face N=256 | 4.71e-5 | 0.0348 |
| cut-cell N=256 | 2.96e-6 | 0.0145 |
| cell-mean N=256 | 3.26e-6 | 0.0493 |

Growth rate of `max|U|` over matched physical windows, N=128 [1/s]:

| window [s] | arithmetic | per-face | cell-mean |
|---|---|---|---|
| 0.020–0.035 | −11 | 56 | 31 |
| 0.035–0.050 | 21 | 134 | 72 |
| 0.050–0.066 | 180 | 238 | 100 |
| 0.066–0.080 | 280 | 9529 | 59 |
| 0.080–0.100 | 386412 | dead | 69 |

The cell-mean arm is the only one that does not enter a runaway acceleration: its
rate stays at 59–100 1/s from t = 0.05 to 0.10. Its delivered curvature variation is
flat at 4–9 1/m through t = 0.05 and then rises three orders **in step with
`minGradPsiBand` leaving 1** (0.998 → 0.858 at t = 0.08 → 0.256 at the end), which
is the non-parallelism residual of plan section 11.2. One run, so a correlation —
but a testable one.

**Not established:** the cell-mean arm is *not* stable — it blows up, just later. And both static gates have constant exact curvature, so
the averaging cost `O(h^2 d^2 kappa)` and the non-parallelism term of section 11.2
are both identically zero there. A varying-curvature static gate is still missing.

### Domain size is now an axis, and 6R is certified equivalent (2026-08-18)

The 3D ladder that has to carry the dimensional confirmation is pinned at its
coarse end by `R/h >= 10`, so with cell-count doubling it spans only a **factor 2
in h** — error ratio 4 end to end, and a 10% endpoint perturbation moves a fitted
order by ~0.14. Widening it on the 0.01 m box is unaffordable (uniform
`h = 2.5e-5` is 6.4e7 cells, ~50 000 core-h).

The droplet fills 0.1% of that box in 3D and `g = (0 0 0)`, so
`p_in - p_out = sigma*2/R` does **not** depend on the box size. Truncating the far
field at fixed h therefore cuts cells as `(L/0.01)^3` while leaving h, dt, the
step count and every metric alone — on a **uniform** mesh, so `w = 1/2` at every
face (where the potential-form identity is exact) and the exact absorption of a
uniform kappa error by the pressure solve both hold by construction. Measured,
three box sizes x two matched cell sizes, all six arms reaching t = 0.1 with the
predicted step counts (9207 at `h = 1e-4`, 26042 at `h = 5e-5`):

| h | L | max&#124;U&#124; | vs 10R | volume | shape | min&#124;grad psi&#124; |
|---|---|---|---|---|---|---|
| 1e-4 | 10R | 1.56067e-05 | — | 4.2177e-04 | 3.2808e-07 | 0.99568 |
| 1e-4 | **6R** | **1.55901e-05** | **-0.11%** | 4.2191e-04 | 3.2859e-07 | 0.995679 |
| 1e-4 | 4R | 1.92075e-04 | **+1130%** | 2.8629e-04 | 3.7515e-07 | 0.996438 |
| 5e-5 | 10R | 1.80287e-05 | — | 2.1758e-05 | 5.2553e-08 | 0.999292 |
| 5e-5 | **6R** | **1.80234e-05** | **-0.03%** | 2.1756e-05 | 5.2553e-08 | 0.999292 |
| 5e-5 | 4R | 2.91047e-04 | **+1514%** | 8.4839e-07 | 1.0848e-06 | 0.998045 |

6R reproduces the reference box to **0.16% on every metric at both cell sizes**,
h-refinement ratio 0.865 vs 0.865. **4R fails by 12–16x and fails UPWARD** — the
near walls inject rather than clip, so this is the reflected-pressure-mode
failure and not the flattering one (wall damping would have masqueraded as an
improvement). At 4R the Laplace jump is still exact (72.7508 vs 72.7505 Pa) and
`min|grad psi|` is intact while `max|U|` and shape blow up and the volume error
*drops* 96%: the interface and the jump are fine, the outer velocity field is
not, which places the parasitic recirculation's extent **beyond 1R**. Its ratio
0.660 vs 0.865 also shows confinement alters the h-trend itself, not just its
level.

Curated: `docs/method-comparison/method-comparison-article/data/tables/domain_size_control.csv`,
scored by `workflow/scripts/domain_size_control.py` (gates: max&#124;U&#124; and
min&#124;grad psi&#124; 5%, volume and shape 15%, refinement-ratio spread 10%, plus the
monotone-in-L direction test — which is why the control needed three box sizes and
not one pair).

Consequence: the 3D ladder runs at **L = 6R** for **379 core-h against 1776** for
the identical four cell sizes on the reference box — 4.7x, uniform mesh,
equivalence measured rather than argued, and the same cost as the octree-refined
alternative, which would have bought hanging-node faces and a constant-kappa gate
for nothing.

The plumbing: `DOMAIN_LENGTH` is a case token (`DOMAIN_HALF_LENGTH` derived for
the droplet centre) and **the capillary step now follows h, not the cell count** —
materialize evaluates it at `nRef = CAPILLARY_REF_LENGTH/h`. The N-based form gave
a truncated box a dt `(0.01/L)^1.5` too large; measured before the fix,
`L = 4e-3, N = 40` got 4.29e-5 instead of 1.0861e-5. Verified bit-identical for
every pre-existing token shape, and `blockMesh` polyMesh output byte-identical at
`DOMAIN_LENGTH = 0.01` in 2D and 3D.

### The wide ladders (2026-08-19): 2D completes, 3D destabilises at R/h ~ 16

`cellCentreInverse` + biharmonicBand theta = 0.2, curated in
`docs/method-comparison/method-comparison-article/data/tables/wide_ladder_coupled.csv`
via `workflow/scripts/make_wide_ladder_table.py`.

**2D (L = 10R, np 8), now four points — h halves each level:**

| N | R/h | max&#124;U&#124; | volume | shape | min&#124;grad psi&#124; |
|---|---|---|---|---|---|
| 64 | 6.4 | 2.2400e-03 | 4.9600e-03 | 3.7800e-05 | 0.952 |
| 128 | 12.8 | 4.1500e-04 | 2.5700e-05 | 5.2400e-06 | 0.996 |
| 256 | 25.6 | 1.0200e-04 | 7.3300e-06 | 8.1400e-08 | 0.9996 |
| **512** | **51.2** | **5.0993e-05** | **1.9658e-07** | **1.4690e-08** | **0.999883** |

Orders: max&#124;U&#124; **+2.43, +2.02, +1.01**; volume +7.59, +1.81, +5.22; shape
+2.85, +6.01, +2.47. Laplace jump 72.7405 against the exact 72.7400 Pa (7e-6).

The fourth point **falsified prediction (b)**: max&#124;U&#124; was predicted at ~2.5e-5
(order 2 continued) and came in at 5.10e-5, **order 1.01**. Volume and shape keep
falling (nothing rose, so the joint-convergence failure mode did not occur), but
the current is now first order. With three points this looked like clean second
order — this is precisely what the fourth point was for.

**3D (L = 6R, np 32), three of four points; N_L = 120 rerunning (54304158):**

| N_L | R/h | max&#124;U&#124; | volume | shape | min&#124;grad psi&#124; |
|---|---|---|---|---|---|
| 60 | 10.0 | 1.0035e-03 | 1.6314e-04 | 4.9531e-06 | 0.9936 |
| 76 | 12.7 | 1.0574e-03 | 1.0514e-04 | 1.8848e-06 | 0.9965 |
| 95 | 15.8 | **7.0362e-02** | 6.1514e-04 | 1.2074e-05 | **0.8586** |

All three reached t = 0.1. **The delivery itself is fine and second order in 3D:**
the WP8.1 delivered-curvature variation at t = 0 falls 3.184 -> 1.942 -> 1.252,
order **+2.09 / +1.97**, and t = 0 max&#124;U&#124; falls 2.12e-4 -> 7.82e-5 -> 4.09e-5
(order +3.58 over the span). The cellCentreInverse property transfers to the
coupled solver in both dimensions (2D t = 0 driver orders +2.12 / +2.14 / +1.83).

**What breaks is amplification, not delivery.** N_L = 60 and 76 follow an almost
identical trajectory (dip to ~1.2e-4 at t = 0.05, rise to ~1.0e-3 at t = 0.1 —
resolution-independent). N_L = 95 departs at t ~ 0.04 and grows exponentially at
**r ~ 79 1/s** (x2.2 per 0.01 s), passing 5e-3 at t = 0.068; min&#124;grad psi&#124;
only leaves 0.95 at t = 0.087, i.e. the profile degradation **follows** the current,
as established.

**The WP0 band mode spectrum is FLAT in every arm, including the unstable one:**
A2h/A4h/A8h stay at 1.00x of their t = 0 values for the whole run. The theta = 0.2
filter holds the normal-direction corrugation pinned even while max&#124;U&#124; grows
70x. What grows instead is the delivered curvature variation: acrossSupport
**182x** and alongInterface **143x** at N_L = 95, against 1.5x/0.9x at N_L = 60 and
4.4x/3.2x at N_L = 76.

So the 3D instability is **not** fed by normal-direction psi roughness — the only
thing the WP0 chain measures. Candidates it leaves open, both dimension-specific:
the two tangential directions along which psi can develop structure the
normal-chain diagnostic is blind to (2D has one, 3D has two), and the Gaussian
curvature term K of the parallel-surface inverse, which is identically zero in 2D
and active in 3D. Neither is measured yet.

**Confinement is ruled out as the cause.** The wall sits 2R = 2e-3 m from the
interface at every level (L is fixed), i.e. 20, 25 and 32 cells at N_L = 60, 76,
95 — the wall gets *further* in cell units under refinement, so a confinement
effect would improve, not degrade. The 2D control measured 6R equivalent to 10R at
0.03–0.11% at h = 1e-4 and 5e-5, bracketing N_L = 95's h = 6.32e-5, and N_L = 60
and 76 sit in the same box and are stable. The definitive 3D cross-check (N = 158
on the 10R box, same h) is available at ~598 core-h / 19 h if wanted.

**Answer to the question the 3D case was built for: dimension matters, sharply.**
The same configuration converges monotonically in 2D through R/h = 51.2 and
destabilises in 3D between R/h = 12.7 and 15.8.

### K is EXONERATED, and it is load-bearing (2026-08-19)

`stationaryDroplet3DwideNoK` — the same 3D ladder with the Gaussian-curvature
term of the parallel-surface inverse switched off, i.e. 3D run through the algebra
2D uses. K is identically zero in 2D, so the switch was validated as
**bit-identical in 2D** and different in 3D before use.

| N_L | R/h | K | reached | max&#124;U&#124; | volume | shape | min&#124;grad psi&#124; | driver t=0 |
|---|---|---|---|---|---|---|---|---|
| 60 | 10.0 | yes | 0.1 | 1.0035e-03 | 1.6314e-04 | 4.9531e-06 | 0.9936 | 3.184 |
| 76 | 12.7 | yes | 0.1 | 1.0574e-03 | 1.0514e-04 | 1.8848e-06 | 0.9965 | 1.942 |
| 95 | 15.8 | yes | 0.1 | 7.0362e-02 | 6.1514e-04 | 1.2074e-05 | 0.8586 | 1.252 |
| 60 | 10.0 | **no** | 0.1 | **3.0343e-01** | **5.5017e-02** | 2.1370e-05 | 0.7798 | **962.6** |
| 76 | 12.7 | **no** | **DIED 0.0955** | 9.2909e-01 | 8.4954e-02 | 1.5291e-03 | 0.0000 | **964.1** |
| 95 | 15.8 | **no** | **DIED 0.0772** | 1.2781e+00 | 2.0926e-02 | 7.2732e-04 | 0.6334 | **959.2** |

Both K-off failures are genuine divergences, not wall-clock kills: arm 00002
logged `task 21: Floating point exception (core dumped)`, arm 00001 a task failure
with rank 19 absent from the terminated list.

**Removing K does not stabilise anything — it destroys the delivery.** The t = 0
delivered non-gradient content goes from 3.184 -> 1.942 -> 1.252 at order **+2.03**
with K, to 962.6 / 964.1 / 959.2 at order **+0.01** without it: 770x larger in
absolute terms and ZEROTH order. Every coupled arm degrades — even the surviving
`R/h = 10` one is 300x worse in current and 340x worse in volume — and the two
finer ones diverge where the K-aware form reached t = 0.1.

The zeroth order is analytically expected, which is worth stating because it means
the measurement is understood and not merely observed. On a sphere kappa = 2/R and
K = 1/R^2, and the inverse reduces exactly:

    with K:     (2/R - 2d/R^2)/(1 - 2d/R + d^2/R^2) = (2/R)/(1 - d/R) = 2/(R - d)
    with K = 0: (2/R)/(1 - 2d/R)                                      = 2/(R - 2d)

so dropping K makes a relative curvature error of `d/(R - 2d) ~ d/R = O(h)`, and
the non-gradient content — which differences kappa across a face — is then `O(h)/h
= O(1)`. Zeroth order, as measured.

**Verdict: the 3D instability at R/h ~ 16 happens DESPITE a second-order K-aware
delivery, and the K term is not the amplifier.** Of the two dimension-specific
candidates, K is eliminated and the **tangential structure** is the one left: 2D
has one tangential direction, 3D has two, and the WP0 chain measures only the
normal direction, which is why it reads flat at 1.00x while the current grows 70x.
Nothing currently instrumented can see it — that diagnostic does not exist yet and
is the next thing to build.

### INVALIDATION: every filtered result predates a psi-filter seam bug (2026-08-19)

An interFoam-vs-leia coupling audit found, and measurement confirmed, that the psi
filter was **decomposition-dependent**. Two defects in `psiFilterEqn.H`:

1. `Lpsi = psi - fvc::average(fvc::interpolate(psi))` inherits its patch types from
   the SECOND operand, and `fvc::average` carries `calculated` patch fields, so on a
   processor patch `Lpsi` was **uncoupled**: the second application of `L`
   interpolated a stored face value instead of the neighbour cell's `L(psi)`,
   injecting a seam-localised perturbation every step into the field whose second
   derivative the curvature model differentiates.
2. The band dilation looped internal faces only, so the filtered **cell set** itself
   depended on the decomposition.

Measured, 3D `N_L = 60` on L = 6R, 20 steps, serial against np=4, **filter OFF as
the control**:

| metric (filter ON) | before | after fix | filter OFF |
|---|---|---|---|
| max&#124;U&#124; | 2.43e-04 | **1.72e-06** | 1.55e-06 |
| A2hL2Band | 8.38e-07 | **3.72e-10** | 4.09e-09 |
| volumeRelError | 1.92e-05 | **6.55e-07** | 3.64e-07 |

53–205× worse than the unfiltered control before; matching it after, with the
control bit-unchanged. Fixed in `f83a1ab`; `config/seamConsistency3D{serial,par4}.yaml`
make it a one-command regression.

**Consequence: the 2D ladder (np=8, 4714–106689 steps) and the 3D ladder (np=32,
9207–18344 steps) both ran with this defect.** The convergence orders, the
`R/h ~ 16` onset and the K verdict all rest on filtered runs. Those numbers are
**not method properties** until re-run. Re-runs are in flight (section 5).

### The 3D "instability under refinement" may not be an h effect at all

| N_L | R/h | dt [s] | steps | r [1/s] | **r·dt per step** | e-folds |
|---|---|---|---|---|---|---|
| 60 | 10.0 | 1.0861e-05 | 9207 | 40.47 | **4.395e-04** | 4.05 |
| 76 | 12.7 | 7.6186e-06 | 13126 | 52.32 | **3.986e-04** | 5.23 |
| 95 | 15.8 | 5.4514e-06 | 18344 | 90.55 | **4.936e-04** | 9.05 |

The rate rises 2.2× but the amplification **per step is constant to 24%**, because
`dt = CAPILLARY_DT_COEFF/nRef^1.5` shrinks with h. So `r ~ 1/dt`, and the e-folds
accumulated by t = 0.1 grow only because a finer mesh takes more steps to get there.
That is an operator-splitting / force-lag signature rather than a spatial-truncation
one, consistent with the delivered non-gradient content converging at +2.03
(R² 0.9997) while the response diverges, and with the earlier capillary-dt sweep
measuring ~90% of the growth rate as dt-proportional at N = 128.

`stationaryDroplet3DwideFixedDt` tests it: all three meshes at one fixed
dt = 5.45144e-6 (the finest arm's own value, so that arm cross-validates both
studies), 18344 steps each. If the two coarse meshes — stable at their native dt —
destabilise purely by taking more, smaller steps, the mechanism is the once-per-step
frozen capillary force and "finer is less stable" is an artefact of tying dt to h.
Note `psiOuterCorrectors` defaults **off**, so the three outer correctors currently
converge momentum and pressure against a force that cannot change.

### Two further defects found by the audit and fixed

- **93.5% of cells were being inverted having never been filled.**
  `applyCellCentreInverseCurvature` looped all cells; for an unfilled one (kappa
  exactly 0) the inverse returns `-2Kd/(1+Kd²)` — identically zero in 2D, but in 3D
  a curvature manufactured from nothing, peaking at 1/R around d ~ R. At 3D
  `N_L = 60`: 201888 of 216000 cells. `nNoFit` read 0 throughout because it counts
  non-finite offsets, which `signedOffset` never returns. A/B (serial): the coupled
  metrics move ~1e-7, so it was benign for the FORCE — alpha is 0/1 that far out —
  but the written kappa field was 93.5% spurious in 3D, and it stops being benign
  the moment the force support widens.
- **The t=0 Young–Laplace solve was unreferenced and under-converged.** Singular
  system (zeroGradient p_rgh everywhere, `needReference()` true) solved without
  `setReference` — it converged only because a closed-domain flux divergence sums to
  zero — and a bare `solve()` inherited `relTol 0.01`, leaving the initial pressure
  ~1% off and handing the first predictor an impulse of order 0.01·sigma·kappa. Both
  fixed. The third asymmetry, unit Laplacian coefficient vs pEqn's rAUf-weighted one,
  is recorded but unfixable without 1/A, which does not exist before the first UEqn.

**The pressure–velocity coupling itself is clean.** Read independently in both
codebases and diffed: UEqn structure, viscous stress (both carry `grad(U)` and
`grad(U)ᵗ`), capillary force in the momentum RHS and `phig` but never the matrix,
`rAUf` cancelling exactly, `p_rghEqn.flux() = rAUf*magSf*snGrad(p_rgh)`, the same
velocity update line, the same `ddtCorr` multiplier outside the cancellation. Ours is
tighter than interFoam's in one respect: it evaluates `GSigma` once into a `tmp` at
pimple scope, where interFoam re-evaluates it in both places. Also corrected: the
"kappa and the differenced indicator come from different discrete objects" issue is
**not** unique to us — interFoam builds kappa from `interpolate(grad(alpha))` while
differencing the same field with `snGrad(alpha)`.

---

### The two-factor law, and what it does to the score (2026-08-20)

`filterOffAmplifier3D` (4/4) and `upwindConvection2D` (8/8) both completed. Full
tables and derivations: `docs/plan-shannon-parasitic-currents.md` sec. 0d.

**The filter is not the source.** With the psi filter removed entirely, the two
fine 3D arms still grow: `A` = +7.68e-04 (R/h = 15.8) and +1.16e-03 per step
(R/h = 20.0), i.e. 3.5 and 7.6 e-folds over a quarter horizon. `A` changes sign
between `R/h` = 12.7 and 15.8 — the stability boundary belongs to the amplifier.

**Every endpoint factorises exactly** as `max|U|(T) = u_0(h)*exp(G(h))`, `u_0` the
one-step kick and `G = A*T/dt` the e-fold count at fixed physical time. Measured:

| | u_0 (should be h^2 force x h^1.5 step = h^3.5) | G |
|---|---|---|
| 3D, filter off | h^+4.22, h^+2.91, h^+3.30 | **h^−3.27** |
| 2D, BDF2/upwind | h^+3.62, h^+3.51, h^+3.69 | **h^−0.95, h^−0.56** |

So refinement shrinks the kick ~11x per doubling and grows the exponent acting on
it ~1.5–2.1x. Currents keep converging only while `G < 3.5/p`, `p` the exponent of
`G` in `h`: budget 4.7 in 2D (spent at `G` = 3.69, `N` = 512 — hence max&#124;U&#124; still
falling, 1.80e-05, the ladder's smallest) and 1.07 in 3D (spent between `R/h` =
12.7 and 15.8 — past it, refining 1.26x costs 26x). **2D and 3D differ in `p` by
~4.3x, not in the sign of `A`. `p` is the number a fix must change** — and it is
scoreable at a quarter horizon, 4x cheaper than `t_blow`.

**The filter's damping flips sign.** Against the θ = 0.2 arms at matched `u_0`
(2.120e-04 vs 2.118e-04 and 4.085e-05 vs 4.061e-05, so same starting state):
`R/h` = 15.8 → the filter is **5.86x better** (G 3.52 → 1.76, half the rate
removed); `R/h` = 10.0 → **1.61x worse** (G −0.21 → +0.27, damping turned into
growth). The biharmonic band filter carries its own weak source, which dominates
wherever there is no corrugation to remove. **θ must therefore be
resolution-dependent** — not as tuning, but because the source term has to vanish
with the content it accompanies. It also explains the non-monotone θ sweep.

**Convective scheme: inert.** `upwind` vs `linearUpwind gradU` agree to 4–5
significant figures across the whole 2D ladder (`u_0` to 5 digits, max&#124;U&#124; ≤0.8% at
`N` = 64 and ≤0.06% for `N >= 128`, volume/shape ≤0.01%). Closes the ∇U hypothesis
for the *stationary* droplet only — it stands for translating and oscillating.

**BDF2 vs Euler, matched window** (1179/3333/9428 steps): Δg = +11.1% / +2.9% /
**−3.0%**, sign-flipping, ±3% for `N >= 128`; volume and shape within 1.2%. The
momentum time discretisation is not the amplifier. BDF2 stays the default on
formal grounds (a second-order foot-point trace should not be fed a first-order
velocity) at no measurable cost. *Retracted:* an unmatched comparison of these
same arms made BDF2 look 3.1x worse at `N` = 512 — that was the horizon, since
`gAvg` is an endpoint estimator and the rate decays. **Never compare `gAvg`
across unequal step counts.**

**2D volume error converges at fourth order**: 3.059e-03 / 1.442e-04 / 8.933e-06
/ 4.978e-07 → orders 4.41 / 4.01 / 4.17; shape 1.219e-06 / 1.520e-07 / 5.846e-08
/ 1.296e-08 → 3.00 / 1.38 / 2.17.

### Two process failures to not repeat (2026-08-20)

1. **`squeue` returning nothing is not an empty queue.** Lichtenberg's primary
   controller `mssd0001` was DOWN (both backups UP) and `squeue` returned zero
   rows, which I read as "all drivers finished". They were all running.
   `scontrol ping` distinguishes the two; the filesystem always does — step count
   from `grep -c '^Time = ' log.*` and the log's mtime.
2. **I deleted a live case directory.** A cleanup loop tested completeness as
   "does the arm have its metrics CSV?", but for interFoam that file is written by
   `post_solve` and so exists *only after* the run ends — every in-progress
   interFoam arm classified as dead. `interFoamDroplet2D_00003` was removed at
   90404/106689 steps (~16 h lost; it was minutes from its wall limit anyway).
   The leia studies were untouched because their solver appends to the CSV every
   step. **Liveness comes from log mtime, never from a declared output, and
   destructive cleanup must not run while any driver may be alive.**

### The full-horizon stability gate, and a RETRACTION of the t_blow baseline (2026-08-31)

`config/fullHorizonStability2D`: `{semi-implicit off, increment} x {cellCentred,
projectedFlux}` at N = 128, filters off, run to the FULL horizon t = 0.1 s
(13,334 steps) on four concurrent Lichtenberg jobs (driver 54434376). Every
arm COMPLETED; none took a floating-point exception.

**RETRACTION FIRST. The 0.078 s blow-up time for `cellCentreInverse` at N = 128
quoted in this file does not reproduce on the current code.** The reproduction
control (off + cellCentred) reached t = 0.1 alive, with max|U| = 4.65e-3 m/s.
It is unmistakably unstable -- the fitted rate over the last fifth of the run is
**+118 1/s**, an e-fold every 8.5 ms, accelerating -- so the instability itself
is intact and reproduces qualitatively; extrapolating that rate puts the blow-up
near t ~ 0.146 s. What has moved is the ARRIVAL TIME, by about 1.9x, which is
well outside the 5-38% seed sensitivity this campaign has measured for t_blow.
Until the cause is found, **no result may be compared against the older t_blow
table** (the 0.100/0.078/0.036 and 0.100/0.063/0.033 rows above). Comparisons
WITHIN the new study are unaffected: all four arms ran on one commit, one set of
binaries, one mesh and one dt law, simultaneously.

| arm | max\|U\| t=0 | final max\|U\| | r, last 20% | volume rel err | shape L2 | min\|grad psi\| |
|---|---|---|---|---|---|---|
| off + cellCentred (control) | 6.61e-05 | **4.65e-03** | **+118.3 1/s** | 9.86e-05 | 3.59e-06 | 0.9170 |
| off + projectedFlux | 6.61e-05 | 3.34e-06 | -52.0 1/s | 1.14e-06 | 1.96e-07 | 0.9980 |
| increment + cellCentred | 4.21e-05 | 7.53e-06 | +10.5 1/s | 8.18e-07 | 2.00e-07 | 0.9981 |
| increment + projectedFlux | 4.21e-05 | 4.55e-06 | -45.3 1/s | 1.19e-06 | 1.96e-07 | 0.9980 |

(shape L2 starts at 1.91e-07 and min\|grad psi\| at 0.9981, so the three
non-control arms end the run essentially where they started after 13,334 steps.)

**The result the campaign was after, and it is not the one the plan predicted.**
The semi-implicit capillary force was the whole point of the A->B->C plan; the
arm that wins carries NONE of it. `off + projectedFlux` -- no artificial
viscosity anywhere, no coefficient to tune -- is the best arm in the matrix on
every metric. Tracing the semi-Lagrangian foot point with the discretely
divergence-free face flux instead of the cell-centred extension is what turns
growth into decay. That is a fix with no tuned constant, which is what the
no-filtering rule actually demands.

Two secondary findings, neither load-bearing:
- The increment form's deferred correction does NOT reach round-off at 3 outer
  correctors. Its lagged term norm contracts ~5x per outer pass but still stands
  at 0.07% of the term's own magnitude on the last pass, so the "zero artificial
  dissipation at outer-loop convergence" property the form was designed around
  is not achieved at this loop setting. The loop axis was dropped on an
  equivalence measured with the term OFF, where the loop genuinely is inert;
  with the term ON the loop is what makes that property true or false.
- `spuriousDispPerStepOverH` in `writeDropletMetrics.H` is computed as
  `divUBandL2*hCell*deltaT/max(hCell,VSMALL)` -- `hCell` cancels, so the metric
  is exactly `divUBandL2*deltaT` and has no h-dependence despite its name. It
  carries no information independent of `divUBandL2` and must not be quoted as
  a second diagnostic. NOT yet fixed (fixing it changes that column mid-campaign).

### ANSWERED: the projectedFlux win is the RECONSTRUCT OPERATOR, not solenoidality (2026-08-31)

`config/traceAmplifierDt2D` (driver 54435583, 12 arms, all COMPLETED with exactly
doubling step counts 13,334 / 26,667 / 53,334 / 106,668). Order parameter on the
fixed window [0.02, 0.1] s, `dr = [G_x - G_projectedFlux]/(T - t_ref)`:

| trace | dt | G | r [1/s] | dr [1/s] |
|---|---|---|---|---|
| cellCentred | 7.50e-06 | **+3.120** | +38.99 | **85.85** |
| cellCentred | 3.75e-06 | **+0.880** | +11.00 | 52.75 |
| cellCentred | 1.87e-06 | -1.529 | -19.12 | 20.31 |
| cellCentred | 9.37e-07 | -2.640 | -32.99 | 4.92 |
| projectedFlux | 7.50e-06 -> 9.37e-07 | -3.748 -> -3.033 | -46.85 -> -37.91 | 0 (ref) |
| reconstructedU | 7.50e-06 | -1.664 | -20.80 | **26.05** |
| reconstructedU | 3.75e-06 | -1.331 | -16.64 | 25.10 |
| reconstructedU | 1.87e-06 | -1.697 | -21.21 | 18.22 |
| reconstructedU | 9.37e-07 | -2.333 | -29.16 | 8.75 |

**The hypothesis that solenoidality is the mechanism is FALSIFIED at the step
that matters.** `reconstructedU` -- smoothed, NOT divergence-free -- already
removes 70% of the cellCentred amplifier (85.85 -> 26.05) with no solenoidality
whatsoever. rho = 0.303 at the base step.

The second control (`reconstructedMomentum`) was then found to be a NULL STEP,
and the 4-rank gate established it for free: with `levelSet.velocityExtension.type
= none` -- the default, and what every stationary-droplet study uses -- the
extension is the IDENTITY, so `Uext == U` bitwise and the arm came back
BYTE-IDENTICAL to `reconstructedU`. The split is therefore fully determined and
proved by bit-identity rather than inferred from a rate:

| mechanism | share of the cellCentred amplifier | how established |
|---|---|---|
| **reconstruct operator** | **70%** | dr 85.85 -> 26.05 |
| velocity extension | **0%** | null step, bit-identity |
| solenoidality | **30%** | dr 26.05 -> 0 |

**This is not a violation of the no-filtering rule, and the earlier framing of it
as one was wrong.** `fvc::reconstruct(linearInterpolate(U) & Sf)` carries NO
tunable coefficient: it restricts the traced velocity to the discrete face-flux
space -- the same space in which the pressure projection enforces the divergence
constraint. The rule targets constants that need re-fitting per resolution (the
theta = 0.2 band filter: 5.86x better at R/h = 15.8, 1.61x WORSE at R/h = 10.0).
A parameter-free operator is a discretisation choice. The statement the data
supports is: **advect the interface with a velocity that lives in the space the
pressure projection actually controls.**

Two method notes worth keeping:
- rho is only meaningful where `dr_cc` is large. It reads 0.897 at dt/4 and 1.78
  at dt/8, but `dr_cc -> 0` as dt -> 0 by construction (all traces converge to the
  same physics), so those are ratios of small numbers -- and rho > 1 means
  reconstructedU is WORSE than cellCentred, which no version of the hypothesis
  predicts. Scored at the coarsest step the same data reads 0.303, the opposite
  verdict. The curation script now refuses to label the weak rungs.
- `dr_cellCentred` does NOT scale as dt: ratios per halving are 1.63 / 2.60 /
  4.13 against a predicted 2.00, i.e. the local exponent climbs 0.70 -> 1.38 ->
  2.05. So `r = r0 + c*dt` with constant c does not describe this. The
  pre-registered corollary DID hold: `G_cellCentred` changes sign between dt/2
  and dt/4, so the control stops growing once the step is small enough.
- `projectedFlux`'s own rate drifts -46.85 -> -37.91 across the sweep, so
  subtracting it as "the shared physical damping" is only approximately valid;
  the three traces are converging toward a common intercept but have NOT
  converged (spread 85.8 -> 8.75). Self-validation is partial.

### RUNNING: does it survive refinement, and does it damp the physical mode? (2026-08-31)

Everything above is N = 128, one case, stationary. Two refinement ladders now on
the cluster, each carrying its OWN matched `cellCentred` baseline because the
t_blow retraction means no historical number is a valid reference:

| study | driver | arms |
|---|---|---|
| `projFluxStationary2D` | **54436905** | 8 -- N = 32/64/128/256 x {cellCentred, projectedFlux} |
| `projFluxOscillating2D` | **54436906** | 6 -- N = 32/64/128 x {cellCentred, projectedFlux} |

Stationary asks whether **G stops growing under refinement** or projectedFlux
merely lowers the floor (the "floor improved, the growth survived" outcome
already on record for cellCentreInverse). Oscillating is the **moving-interface
gate** required before promotion, and its discriminator is **period and decay
rate against Lamb**, not max|U|: the reconstruct operator keeps only face-normal
components, so if it damps the physical capillary mode this is the case that
shows it.

Both are NEW configs, not edits: `config/stationaryDroplet2D` carries
`export_slides: true` and would have overwritten the curated paper table.

**The 4-rank gate caught a silent no-op here too.** The oscillating case's
`fvSolution.template` had NO `traceVelocity` entry, so every arm would have run
at the solver default whatever the config said -- a two-arm comparison in which
the arms differ in nothing, completing cleanly and producing a plausible,
meaningless table. It surfaced only because materialization reported
`axes={} -> 1 variations`. Tokens ported; verified afterwards on 4 ranks that the
axis expands to 2, the arms differ, and the cellCentred arm is BYTE-IDENTICAL to
the pre-change run.

### DECIDED: hanging nodes do not break the balanced force -- the constant-curvature gate (2026-09-04)

The objection `stationaryDroplet3Dwide` raised against octree refinement (w = 1/2 broken
at 2:1 transition faces, so the potential-form identity and the exact absorption of a
uniform kappa no longer hold) was put to the measurement it asked for:
`constantCurvatureSurfaceTension` with kappa = 2/R = 2000 on the refined N_fine = 60
meshes (interface band 51 640 cells / 1 688 transition polyhedra; ball control 48 784)
and on the uniform N = 60 mesh (216 000), 921 steps to t = 0.01, np 4 (+ np 1 for the
band). Read out with `make_refined_equivalence_table.py` (`max_t` rows and t = T):

| arm | ranks | max_t mean\|u\| | max_t L2\|u\| | mean\|u\|(T) | L2\|u\|(T) | Laplace jump |
|---|---|---|---|---|---|---|
| refined, interface band | 4 | 3.1e-10 | 3.8e-10 | 2.2e-11 | 1.3e-10 | 145.470 |
| refined, ball control | 4 | 3.5e-10 | 4.0e-10 | 2.3e-11 | 1.8e-10 | 145.470 |
| refined, interface band | 1 | 4.6e-10 | 5.3e-10 | 2.0e-11 | 1.4e-10 | 145.470 |
| uniform N = 60 | 4 | 7.2e-10 | 8.2e-10 | 3.7e-12 | 6.0e-12 | 145.470 |

Same Laplace jump to the digit (exact 145.48 Pa, -6.9e-5 relative) on refined and
uniform meshes; the spurious velocity is at the pressure solver's round-off floor in
every arm -- the refined PEAKS are below the uniform peak, the refined ENDPOINTS
(1.3-1.8e-10 m/s) are above the uniform's decayed 6e-12 -- all five orders below the
parasitic currents the real curvature produces on these meshes (~1e-5 m/s). The
pre-registered falsification (refined > 1e-9 while uniform at round-off) did NOT occur.
**I mis-set the pass line**: 1e-12 at every step is below the solver's floor and no arm
reaches it, uniform included; the decision rests on refined-vs-uniform, which was the
falsifiable half. Volume, shape and kErrL2Band agree to 1.2 % (the shape difference is the
initial discretisation of the sphere on two meshes). Consequence: the hex route proceeds
to the production-force ladder (G3/G4); the octree rejection in
`config/stationaryDroplet3Dwide.yaml` and in section 4 below is AMENDED by this block --
equal cost was its other argument, and at 4.2x fewer cells that one is gone too.

**G1 on this gate is blind by construction** (every dynamic column is round-off, so np 1
vs np 4 differ in reduction order only; the Laplace jump and the geometric metrics agree).
The decomposition test with the PRODUCTION force ran on the laptop:
`stationaryDroplet3DrefinedGate4` (np 4) vs `stationaryDroplet3DrefinedGate1` (np 1),
same refined mesh, 921 steps each, both COMPLETED. **Not the 1e-10 I pre-registered**:
at t = 0.01 the two agree to 1.7e-4 (mean|u'|, 2.838e-6 m/s), 5.8e-5 (L2|u'|, 1.172e-5),
3.6e-4 (volume error 2.167e-6), 7.9e-4 (shape 5.50e-7 m), 1.6e-8 (Laplace jump
145.628 Pa), 1.9e-4 (kErrL2Band) relative, and to four digits in the maxima over time
(1.657e-5 / 5.655e-5 m/s). That tolerance was unattainable by construction -- the
pressure solve converges to a RELATIVE tolerance and GAMG's agglomeration depends on the
decomposition -- and this repository has already measured a 30x tighter pressure solve to
move exactly these metrics by 2.4e-4 relative (CLAUDE.md, cancellation-dominated
problem). So 1e-4 is the solver's decomposition dependence, PROBABLY; the yardstick that
turns that into a verdict is the same np 1 / np 4 pair on the UNIFORM N = 60 mesh, running
now to t = 0.003 (`stationaryDroplet3DuniformGate4/Gate1`): PASS = the uniform pair's
relative differences are of the same order as the refined pair's at t = 0.003; a refined
pair 10x or more apart is a defect at the transition polyhedra. Second mis-set threshold
of the day, same cause: I wrote round-off-level lines for quantities the linear solver
limits at its tolerance.

### RUNNING: static local refinement around the interface (2026-09-04)

**Why.** R/h >= 10 on a uniform 3D mesh costs 216k (hex, N = 60) to 3.6M (poly,
R/h = 25.6) cells for a droplet that occupies 0.1 % of the box. The method touches the
interface only in a band: the phase indicator fills `0 < alpha < 1` cells (layer 1),
the capillary force acts on their faces (layer 2), each face curvature comes from a CPC
fit reading one cell-point-cell ring (layer 3), one semi-Lagrangian step and the
explicit viscous transpose term read one more (layer 4). **Four complete fine layers
each side is the stencil minimum; `REFINE_BAND_CELLS 6` is the margin.**

**Mechanism -- pre-processing only, solver untouched** (commit `ac77c4a`). Per pass:
`0/ := 0.org` -> `leiaSetFields` -> `topoSet` (seed = cells with `0 < alpha < 1`,
i.e. the `snGrad(alpha)` support, plus face dilations) -> `refineHexMesh refineCells
-overwrite` (hexRef8; `cellLevel`/`pointLevel` PERSIST, so 2:1 holds across passes --
OpenFOAM's plain `refineMesh` does not persist them and is never used). After the last
pass `0/ := 0.org` -> `leiaSetFields` on the FINAL mesh: **nothing mapped survives**
(a mapped alpha is a smeared alpha; verified `U`/`p_rgh` cmp-identical to `0.org`,
`leiaSetFields` idempotent on the final mesh). Poly: cfMesh has no in-place refiner and
v2512's only one (hexRef8) is hex-only, so `leiaRefinePolyMesh.py` re-meshes with the
psi = 0 iso-surface (STL) as a `surfaceMeshRefinement` surface -- not yet exercised.

**MEASURED, and it corrected the plan.** Face dilation is Manhattan growth: at N = 60,
3 / 4 / 5 dilations leave **4.6 / 5.8 / 7.0** complete fine layers at the worst point
(the sphere's diagonals), not the 6 / 8 / 10 an axis-aligned count predicts. The driver
therefore adds dilations until the psi it already has on the current mesh proves
`2|psi|/h - 1 >= REFINE_BAND_CELLS` for every unselected cell; psi verifies the width,
the criterion stays the alpha seed. Interface band at N_fine = 60: 5 dilations,
**51 640 cells against 216 000 uniform (4.2x)**, 6.97 layers, 1 688 transition
polyhedra, max non-orthogonality 25.2 deg. Ball control (whole interior fine): 3
dilations, 48 784 cells, 7.19 layers.

**Gates passed (laptop).** G-1 render: five re-rendered studies differ only in two
blockMeshDict comment lines, the new `curvature` entry (3D case) and the new tokens in
`case_params.json`. G-1 run: uniform N = 30, 65 steps, HEAD vs new templates --
**cmp-IDENTICAL** solver CSV. G0: both `refinedWB` arms PASS the band check, `N_CELLS`
pin exact (the capillary dt handle MUST be the fine spacing; the check fails otherwise).

**Running now.** The hanging-node question that `stationaryDroplet3Dwide` raised
(w = 1/2 broken at transition faces) is answered by the constant-curvature well-balanced
gate: `constantCurvatureSurfaceTension` with kappa = 2/R = 2000 on the refined mesh
(interface band, ball control), its serial twin (np 1 vs np 4) and the uniform N = 60
control. Pre-registered: refined arms at the uniform arm's round-off floor -> the
objection is answered; > 1e-9 m/s while uniform is at round-off -> it stands, stop hex,
go poly. Refined arms LANDED (laptop, np 4, 921 steps each): interface band max over
time mean|U'| 3.1e-10, L2 3.8e-10 m/s; ball control 3.5e-10 / 4.0e-10; Laplace jump
145.470 Pa against the exact 145.48. Round-off scale, above the pre-registered 1e-12 --
the uniform control's floor decides whether that is the pressure solver's floor or the
mesh. Laptop chain: `refinedWB` -> `refinedWBserial` -> `uniformWB` (profiles/local8).
**GC PASSED**: the same `refinedWB` ran as a `profiles/slurm` study (orchestrator
54477820, solves 54477827/28 in `.my_jobs`): the driver built both meshes in the serial
`case_pre` job with exactly the laptop's statistics, and the cluster solver CSVs are
cmp-IDENTICAL to the laptop's -- every one of 66 columns at every one of 921 steps
(`workflow/scripts/compare_metrics_csv.py`).

**P0 (polyhedral driver) smoke, laptop, `polySmoke3D` box at maxCellSize 1.5e-4.** The
loop runs end to end: iso-surface STL (4846 facets) from the coarse psi, cfMesh
`surfaceMeshRefinement` around it, re-initialised fields, band check PASS. Two things
MEASURED that the plan had assumed: (1) a `cellSize` request of half the base gave cells
TWO octree levels finer (2.98e-5 for 7.5e-5 asked) and 6x the cells -- cfMesh rounds a
size to its octree -- so the driver now requests whole levels
(`additionalRefinementLevels`, the meaning of REFINE_LEVELS): 190 243 -> 337 907 cells
(1.78x), interface cells at 5.95e-5 = one level below the 1.04e-4 base; (2) cfMesh
grades beyond the requested thickness: 11.8 complete fine layers for 6 asked. The size
profile vs distance (`sizeProfileOverHfine` in refinedBand.csv) is the honest picture
on polyhedra, where one octree level's dual cells scatter in size: median size in fine
cells, by distance from the interface in fine cells, `0-18: 1.00 | 18-24: 1.56 |
24-48: 2.0 | 48-72: 1.03` -- fine out to 18 for 7.6 requested, the base level beyond,
and fine AGAIN at the box edges and corners: cfMesh's feature refinement of the domain
boundary, present in every pMesh mesh including the uniform rungs (hex, for comparison:
`0-9: 1.00 | 9+: 2.00`, exactly the requested band). N_CELLS for a polyRefined study is
pinned from the built mesh (`N_CELLS_suggested`), as for the uniform rungs.

**First production-force pair LANDED (cluster, N_fine = 60, refined np 8 vs uniform np 32,
2302 steps each, t = 0.025).** Read with `make_refined_equivalence_table.py`:

| metric | t | refined | uniform | rel. diff | pre-registered |
|---|---|---|---|---|---|
| max_t mean\|u'\| | -- | 1.657e-5 | 1.660e-5 | -0.18 % | 10 % |
| max_t L2\|u'\| | -- | 5.655e-5 | 5.656e-5 | -0.013 % | 10 % |
| mean\|u'\| | 0.0125 / 0.025 | 9.10e-7 / 2.59e-7 | 9.10e-7 / 2.49e-7 | +0.011 % / +3.9 % | 10 % |
| L2\|u'\| | 0.0125 / 0.025 | 4.626e-6 / 1.155e-6 | 4.626e-6 / 1.154e-6 | -0.013 % / +0.13 % | 10 % |
| volume error | 0.0125 / 0.025 | 2.438e-6 / 6.19e-7 | 2.438e-6 / 6.21e-7 | +0.012 % / -0.24 % | 5 % |
| shape (zeroSetRadialL2, m) | 0.0125 / 0.025 | 5.45e-7 / 5.47e-7 | 5.41e-7 / 5.42e-7 | +0.81 % / +0.91 % | 5 % |
| Laplace jump rel. error | 0.025 | 1.010e-3 | 1.010e-3 | +1.4e-4 | 1 % |
| kErrL2Band | 0.025 | 2.200 | 2.203 | -0.13 % | 5 % |
| A2hL2Band | 0.025 | 8.073e-7 | 8.073e-7 | ~1e-6 | 5 % |
| core-hours | -- | 1.22 | 4.12 | 3.4x cheaper | recorded |

Every metric of the vector inside its tolerance; the two meshes are equivalent to the
level at which two DECOMPOSITIONS of the same mesh differ (1e-4, see G1 above -- refined
at np 8 against uniform at np 32 sits at exactly that level). Three more rungs and the
controls are running; the ORDERS wait for all four.

**Then (in order) -- SUBMITTED 2026-09-04 evening after the gate above, orchestrator ids
in `.my_jobs`.** G3/G4 on the cluster: `stationaryDroplet3Drefined` (fine N =
60/76/96/120, one level), controls `refinedL2` (two levels at 120) and `refinedBall`,
uniform twins `stationaryDroplet3Duniform` (60/76/96, RE-RUN on the current code -- the
wide ladder was filtered and predates alg_lin) and `uniform120` (1.73M cells, ~90
core-h, exactness over CPU hours). Pre-registered tolerances in the config headers:
volume, shape, kErrL2Band, A2hL2Band within 5 %, pLaplace within 1 %, L1/L2 of |U|
within 10 % with the same trend; orders within +-0.3 of the uniform ladder's. Then
P0-P2 for `polyDroplet3Drefined_r13p8`. Curation:
`make_refined_mesh_table.py`, `make_refined_equivalence_table.py`,
`make_refined_mesh_fig.py`.

**Polyhedral uniform ladder (the study this refinement work is for).** `polyDroplet3D_r13p8`
(572k cells, N_CELLS 84) COMPLETED its 3813 steps to t = 0.025 without blow-up:
mean|U'| 7.5e-6, L2|U'| 3.1e-5 m/s, volume error 2.6e-6, shape (zeroSetRadialL2)
4.9e-7 m = R/2060, Laplace jump 145.546 Pa against the exact 145.48 (+0.045 %),
kErrL2Band 1.68 of 2000 (0.08 %), min|grad psi| 0.997. `r18p9` at t = 0.0171 of 0.025
(mean 1.5e-5, L2 5.3e-5) and `r25p6` at t = 0.0054 (mean 1.3e-6, L2 5.7e-6) are still
running -- compare at EQUAL time when they land, not at their current endpoints.

## 5. Lichtenberg — what is running

Login: `ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de`
Repo on the parallel file system: `/work/scratch/tm83tomy/leia`
Account `special00004`. Every job **must** set `--mem-per-cpu`.

| job | what | limit | output |
|---|---|---|---|
| `54354379` `leia-curv` `long` | **interFoamDroplet2D** — re-running the `N` = 512 arm only (the other three are complete at the full 0.1 s horizon with `interFoam.csv` present). 106689 steps at the measured 1.52 steps/s = ~19.5 h. Replaces the arm lost to the cleanup bug above. | 28 h | `studies/interFoamDroplet2D/` |
| DONE | **filterOffAmplifier3D** 4/4, **upwindConvection2D** 8/8, **upwindConvection3D** 4/4, **filterThetaScaling3D** 6/6 — all analysed, section 4. | — | `studies/*/` |
| DONE | **stationaryDroplet3Dwide**, **cellCentreInverseFiltered512**, **domainSizeControl10R/6R/4R**, **psiOuterCorrectorsGain3D**, **ddtOrderGain3D** | — | `studies/*/` |
| `54477820` `leia-curv` `long` | **stationaryDroplet3DrefinedWB** GC (2026-09-04): the static-refinement mesh rule as a serial `case_pre` job + the constant-curvature well-balanced gate on the refined mesh, np 4, two arms x 921 steps. From `/work/scratch/tm83tomy/leia-curvature`. | 7 d (orchestrator) | `studies/stationaryDroplet3DrefinedWB/` |
| `54477891`-`54477895` `leia-curv` | **static-refinement ladder** (2026-09-04 evening, submitted after the hanging-node gate passed): orchestrators for `stationaryDroplet3Drefined` (fine N = 60/76/96/120, np 8), `stationaryDroplet3DrefinedL2`, `stationaryDroplet3DrefinedBall` (controls at 120), `stationaryDroplet3Duniform` (60/76/96, np 32) and `stationaryDroplet3Duniform120` (1.73M cells, np 64, ~90 core-h). Read with `make_refined_equivalence_table.py` when ALL arms of a pair are complete -- never before. | 7 d (orchestrators) | `studies/stationaryDroplet3D{refined*,uniform*}/` |
| running | **polyDroplet3D_r18p9** (1.46M cells, np 32) at t = 0.0171 / 0.025 and **polyDroplet3D_r25p6** (3.6M, np 64) at t = 0.0054 / 0.025; **polyDroplet3D_r13p8** DONE (3813 steps, section 4). | 23 h | `studies/polyDroplet3D_*/` |

Not from this work, do not cancel: the SDPLS session's jobs in
`/work/scratch/tm83tomy/leia` (`of-bo-cht-*`, `23fbd13d-*`, `leia-studies` on
`long`). **Curvature work runs only from `/work/scratch/tm83tomy/leia-curvature`
with `$HOME/OpenFOAM/curvature-v2512` binaries and job name `leia-curv`;** never
pull, stash, `Allwmake` or edit `profiles/*` in theirs.

> **First launch of these two (jobs 54140311/312) was thrown away.** Both timed
> out at 4 h having advanced 111 of 4714 steps — not physics: snakemake's
> jobstep wrapper was confining all 8 ranks to one CPU. Fixed in
> `profiles/slurm` (`--cpu-bind=none`, commit `aad3284`); measured 0.0077 ->
> 86 steps/s. Check `wall/CPU` (ClockTime/ExecutionTime) on any cluster study:
> `~np` means confined, `~1` means healthy.

Score the two footEval studies on **r(A2h)** — the corrugation growth rate, the
order parameter — not on t_blow, whose e-fold count varies 5.4–13.3 across the
matrix. Report **volume and shape error together**; the N=64 local smoke already
showed them disagreeing in sign (footEval 3.3x better on volume, 1.4x worse on
shape, on the stationary control).

Check what is yours at any time:

```
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "squeue -u tm83tomy"
```

Watch the parasitic-current trace grow (`TIME,maxMagU` is the first two columns):

```
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "tail -1 /work/scratch/tm83tomy/leia/studies/stationaryDropletFootEvalFace/stationaryDroplet2D_00000/leiaSemiLagrangianLevelSetTwoPhaseFoam.csv | cut -d, -f1,2"
```

A finished run either **blows up** — a floating-point exception with `maxMagU`
past 1e50 and the time step collapsed, which is a *result*, not a crash — or
reaches `endTime`. Both leave the full CSV.

**Other jobs on this account are not from this thread.** A separate study
orchestrator (`leia-studies`, and the snakemake-spawned hash-named jobs) runs the
transport/SDPLS studies. Do not cancel those.

## 6. Lichtenberg — how to run things

Pull the newest code first, then rebuild what you touched:

```
cd /work/scratch/tm83tomy/leia && git pull --rebase
source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake
```

`./Allwmake` builds the library, all solvers and all test apps. To rebuild one
target only, e.g. `wmake applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam`.

**Gotcha that has bitten twice:** sourcing the OpenFOAM `bashrc` returns 1 on its
last internal call, so any script with `set -e` dies silently at that line. Wrap
it as `set +eu; source ...; set -e` — the same thing `workflow/Snakefile`'s `sh()`
helper does.

**Static gates + gain (serial, minutes, no MPI).** These are the cheap screening
runs — a new delivery should be scored here before any coupled run:

```
cd /work/scratch/tm83tomy/leia
snakemake --workflow-profile profiles/local --nolock \
  --configfile config/faceCurvatureDroplet2D.yaml --config export_slides=False
snakemake --workflow-profile profiles/local --nolock \
  --configfile config/faceCurvatureSphere3D.yaml --config export_slides=False
```

Both apps run per case automatically: `leiaTestMeanCurvature` writes
`leiaTestFaceCurvature.csv` (accuracy) and `leiaTestCurvatureNoiseGain` writes
`leiaTestCurvatureNoiseGain.csv` (gain), side by side in each case directory.

**Coupled droplet runs (MPI).** Prepare on the login node, then submit the solve
as a standalone `sbatch` — this is the pattern that works:

```
snakemake --workflow-profile profiles/local --nolock \
  --configfile config/stationaryDropletCellMean.yaml --until decompose --cores 4
```

then submit; a ready-made script is on the cluster at
`/work/scratch/tm83tomy/leia/submit_cellmean.sh` (edit the config name, the job
name and the time limit at the top). Keep 2D runs at **8 ranks** with time limits
sized to the measured pace (~0.03 s of physical time per wall-clock hour at N=128
on 8 ranks) — do not pad to 12 h or 48 h.

**Do not** delete `processor*/` subdirectories with a filter like
`! -name 0`: that removes `processor*/constant` and with it the decomposed mesh.
Use `decomposePar -force` to rebuild if it happens.

**Bring results back** (raw output travels by rsync, never by git):

```
make pull-study STUDY=stationaryDropletCellMean
```

Available configs for this thread: `config/stationaryDropletStableFoot.yaml`
(per-face), `config/stationaryDropletCutCell.yaml` (cut-cell),
`config/stationaryDropletCellMean.yaml` (cell-mean),
`config/faceCurvatureDroplet2D.yaml` and `config/faceCurvatureSphere3D.yaml`
(static gates). `workflow/README.md` documents every study.

## 7. Next, in order

*(2026-08-14: items 1 and the psi-transform half of item 3 are DONE and measured
— symmetric face averaging lost the order (plan sec. 14), and the non-parallel
foliation is now gated directly by `config/faceCurvatureEllipsoidPsi2D.yaml`
(plan sec. 17). The time/coupling axis is closed: capillary-dt sweeps, GAMG/
decomposition/adaptivity exonerations, the psiOuterCorrectors falsification and
the x_d kernel verification are plan sec. 16. THE order parameter is now the
corrugation growth rate r(A2h) ~ 14–160 1/s (N = 64–256), nearly dt-independent
at N=256 — every successful lever so far bought prefactor only. Leads: the
closestPointNewton evaluation inverts to 2nd order with a 46x smaller constant
than the production delivery under the broken foliation, when the fit is exact
(sec. 17.3, noise gain unmeasured); production itself is CONFIRMED second order
wherever the foliation is parallel — sec. 17.1 corrected; and item 2 below
remains untouched.)*

1. **Done — superseded.** (Symmetric face averaging: measured order 1.10, plan
   sec. 14. The gain–order trade-off on the delivery axis is closed.)
2. **Fix the signed-distance assumption where it actually survives.** It is not in
   the curvature — it is in the phase indicator, which locates the interface with a
   *linear* least-squares fit and a first-order `psi/|grad psi|` offset
   (`applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H:154-161`
   and `src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C:130-140`), so it
   sits in `snGrad(alpha)`, the other factor of the CSF force. The Hessian-corrected
   root already exists as `offsetDistance` and costs zero extra derivatives. It
   moves `alpha`, so it must be re-gated on transport and volume conservation.
3. **Done — see section 4.** The varying-curvature gate now exists
   (`config/faceCurvatureEllipse2D.yaml`, `cases/ellipseDroplet2D`) and is what
   retracted the cell-mean adoption. Two follow-ons remain worth doing: a 3D
   companion (a torus — exact signed distance with non-constant mean curvature) to
   confirm the same ranking in 3D, and a ψ-transform arm
   (ψ → ψ·(1 + b cos mθ), zero set held fixed) which is the only way to exercise
   the non-parallel-foliation regime, since every signed-distance gate has D = 0.
4. **Test the mechanism on data already owned, no new runs.** ψ is written at
   t = 0, 0.02, …, 0.10 in the cell-mean N=128 case. Compute
   D = Δ_Γ(ln\|∇ψ\|) − \|∇_Γ ln\|∇ψ\|\|² there and check whether d·D accounts for the
   band curvature error going 70 → 187 → 1100 1/m as `minGradPsiBand` falls. Noise
   is irrelevant for a diagnostic, so the wide stencil that could never go in the
   solver is fine here.
5. Still open from the program document: T1.3 (frozen-band event characterization)
   and T2 (offline anisotropic-fit study).

**Acceptance criterion** for any future curvature delivery, checked in the same
table: `G h^2 <= 0.65` **and** fitted `E_L2` order `>= 1.9` **on the ellipse gate,
not on the circle**. A candidate that meets
only the second is, on the recorded evidence, a regression.

## 8. Repository sync

Hub is GitHub `leia-openfoam/leia`, branch `development`. Code moves by git; raw
simulation output (`studies/`, `runs/`) moves by rsync and is git-ignored, along
with `.snakemake/`, built deck `*.html` and article PDFs.

At the time of writing: laptop, GitHub and Lichtenberg are all on the same commit.

The laptop working tree also carries uncommitted regenerated study tables and
figures from a **different** thread (the linear/nSL semi-Lagrangian studies) and a
cluster stash `stash@{0}` with the same. Those are deliberately not committed
here — they belong to whoever is running those studies.
