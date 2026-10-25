# STATUS — where the work stands

Living hand-off file. Written to be usable from a phone: every command below is
meant to be run **on Lichtenberg**, and nothing here needs a local OpenFOAM.

Last updated: 2026-10-14 (section 11.27: the static test of the Eulerian source terms at zero velocity. The first cluster run is VOID: at u = 0 the psi matrix is diagonal, a small source update started below the solver tolerance, and PBiCGStab did nothing in 68 % of the psi solves (39 cases never changed psi); stage C was cancelled; fixed with minIter 1 and a solver-iteration check (S0); the laptop gate passes except the bit-inertness prediction G3, which was wrong about PBiCGStab (the change is round-off except for softWall on the whole domain); the re-run on Lichtenberg is next. Part 1 corrected: only softWall (and less cubicQ) on the whole domain is decomposition-dependent. section 11.26: the combined source F = a + lambda (1 - q): the centred q of the source step fails at zero velocity and its growth law is retracted; the one-sided q `upwindFit` passes gate E2.1 at order 3.6; the kinematic arms of the 2D method gate fail on the shape order for the combined source and both single terms; the band breaks the strain term in 1D; the 2D shape loss is open. section 11.25: the prescribed flux is projected to divergence-free on polyhedral meshes, after fixing two defects of correctFlux; no transport metric degrades by more than 1.1 %. Section 11.24: L_INF_E_PSI fixed in the shared error writer, which served three solvers; one writer for four solvers. Section 11.23: the kinematic SL solver computes the narrow band before alpha and stays, one kinematic solver per flow solver; the CFL 1.0 row of the SL 2D convergence table is retracted. Section 11.22: leiaRedistancedLevelSetFoam retired, leiaLevelSetFoam reproduces it bit for bit; the semi-Lagrangian solver is not yet reproducible by leiaLevelSetFoam. Section 11.21: Build Tests runs on OpenFOAM-v2606 in our own CI image and passes; Lichtenberg has OpenFOAM-v2606, and the new clone leia-dev runs the development line there through profiles/slurm. Section 11.20: the feature is merged into development and main, the knowledge-base site is live, two branches are pruned. Section 11.19. `2Dtranslation` was a REVERSED flow until today, so every earlier number of the case is void; re-run one-way: the cone bound is 14x worse than no bound at N = 64 and does not converge, the unbounded scheme converges at order 2.90, 2.16, 1.86. Fixed psi values on a patch make the SL update unstable. The production curvature cellCentreInverse is second order on the signed-distance ellipse and ellipsoid, and its gain equals the per-face inverse's (0.647 at N = 512). Earlier, 2026-09-28: FALSIFIED: the distance-cone bound is NOT a transport bound. The advection ladder settles it -- on pure advection with every arm on the IDENTICAL mesh it is WORSE than no bound on uniform translation, where its L = 1 is exactly valid (the 2.3-3.6x of that date was read from the reversed flow; one-way 14x), 9-19x worse in strained flow, and the stencil L mode LOSES THE PHASE on 3D shear hex (E_VOL_REL = 1.0000). The resolution ladder agrees: the eikonal error is 7.119e-03 at N = 64 and 7.066e-03 at N = 128, a FLOOR at order 0.01, and the centroid error reverses to +119 % worse than none. RETRACTED: the earlier 'every interface metric improves 25-68 %' holds only at N = 64 on Popinet's NEARLY UNIFORM velocity, the single regime where L = 1 is defensible. The unbounded POLYHEDRAL advection arm does diverge at step 198 and every bound prevents it, but the falsified monotone clip beats the cone bound on every advection gate. WHAT SURVIVES: the slValueBound family (none | stencilBounds | lipschitzCone as one gated, byte-inert study axis, the seat for the review's Rank 1) and the -mode growth instrument, self-validated against the power iteration. See METHOD.md 8.3.)

Conventions this file assumes are already known: [CLAUDE.md](CLAUDE.md) (layout,
build, git discipline) and [CLUSTER.md](CLUSTER.md) (full verified cluster
workflow). This file is the *current state*, not the manual. Since 2026-09-28 the concise, cross-linked record of every method line (the decisions, the retractions, why something failed, what is open) is the knowledge base `docs/knowledge-base/` (start at its `index.md`; section 11.18).

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
- MARKED 2026-09-28: the studies behind `cases/default.parameter` lines 824-876
  (`MASS_FLUX_ALPHAF_TIME_LEVEL`, `MASS_FLUX_PROJECT`, the +5 to +8 % inflation) also ran on the closed
  box: `alphaFTimeLevelTranslating2D` (5cd2c98, 08:31), `massFluxComparison2D` (86f8333, 09:38) and
  `translatingClearOutlet2D` (65fc13c, 00:07), all on 2026-09-02 before the fix 440107f (17:25). Their
  numbers are void; the block carries a VOID marker and stays as history.

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

RE-RUN IN SI UNITS 2026-09-08 (D = 1 mm, rho = 1000, nu = 1e-6, sigma = 0.012 N/m,
U = 0.0693 m/s -- see "the Popinet cases run in SI units" below). The whole ladder and the
Laplace sweep were re-run: N = 128 and N = 256 and the four sweep arms on the cluster
(54500576 / 54500577 / 54500578), N = 64 on the laptop. All 28 curated entries reproduce the
numbers above to at most 0.051 %, the step counts match (1563 / 4420 / 12501), and the fitted
orders are unchanged: L2 0.91 (R = 0.999), Linf 0.49 (R = 0.999), shape 1.71 (R = 0.996). The
shape column is now normalised by R, so the table reads 4.08e-3 / 1.04e-3 / 3.78e-4 for the
same measurement. The table is therefore regenerable from the committed SI configs.

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

RETRACTED 2026-09-27. Three of the four studies below (`wellBalancedTranslating2D`,
`bestConfigTranslating2D`, `rhoDdtGate2D`) ran on the closed box and are VOID (list above).
The conclusion also contradicts the amplifierGate2D verdict above: the curvature estimator is
the source. Do not cite this subsection.

`kinematicTranslation2D` (transport 2nd–3rd order, bounded-α error exactly 0);
`wellBalancedTranslating2D` (**exact constant curvature still diverged**);
`bestConfigTranslating2D` (the required extension does not help, hurts coarse rungs);
`rhoDdtGate2D` (density ddt is not the mechanism).

### Leading open defect

RETRACTED 2026-09-27: the +5 to +8 % inflation below comes from `translatingClearOutlet2D`,
which ran before the fix 440107f (closed box). It is not a measurement of the repaired case.

The droplet **inflates** monotonically, +5% to +8% over the horizon, and at matched
times the gain does **not** converge under refinement. Untouched by the mass-flux
model — α comes from ψ, so this is semi-Lagrangian transport under the computed
velocity field.

### New, runtime-selectable, default-off (nothing existing moves)

CORRECTED 2026-09-27: `boundRho` defaults to true since 28a1383 (2026-09-02), and with
`MASS_FLUX rhoLENT` (82ca995, the same day) it is active in every study that does not override it.

`massFlux { boundRho true; }` (clip ρ to `[rho2, rho1]`; the clip is **reported**),
`massFlux { massResidualDiagnostic true; }` (mass residual for every model),
`DROPLET_OFFSET_X` (droplet start placement). All three gated for bit-identity and on
4 MPI ranks before any cluster use.

`levelSet { semiLagrangian { valueBound ...; } }` — the bound on the reconstructed value,
runtime-selected: `none` | `stencilBounds` (the falsified monotone clip) | `lipschitzCone`
(the distance-cone bound). The default is the sentinel `fromClipSwitch`, which follows the
legacy `clipToStencilBounds`, so nothing existing moves. Its three cone entries
(`lipschitzMode`, `lipschitzConstant`, `onInadmissible`) are read by `lipschitzCone` only.
Tokens `SL_VALUE_BOUND`, `SL_CONE_L_MODE`, `SL_CONE_L`, `SL_CONE_INADMISSIBLE`. Gated for
bit-identity on 8 arms over 1563 steps at np = 4, plus an exact-field unit gate.

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

CORRECTED 2026-09-27: since c935883 (2026-09-01) the default is `cellCentreInverse`
(METHOD.md Sec. 4.1 and 8.1). The paragraph below is the state before that commit.

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

SUPERSEDED 2026-09-28 marker: the filtered coupled results below (psiFilter biharmonicBand theta = 0.2
and 0.05, and the transport-order comparison read in filtered arms) predate the psi-filter seam bug;
see the INVALIDATION of 2026-08-19 (section 4) and docs/knowledge-base/retractions/psi-filter-seam-bug.md.

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

SUPERSEDED 2026-09-28 marker: the filtered arms below (psiFilter biharmonicBand theta = 0.2 and 0.05)
predate the psi-filter seam bug; see the INVALIDATION of 2026-08-19 (section 4) and
docs/knowledge-base/retractions/psi-filter-seam-bug.md. The unfiltered cellCentreInverse arms are outside
that invalidation; every np > 1 coupled run before 2026-09-27 carries the parallel defects of section 11.14.

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

SUPERSEDED 2026-09-28 marker: both ladders ran with psiFilter biharmonicBand theta = 0.2 and predate the
psi-filter seam bug; see the INVALIDATION of 2026-08-19 (section 4), the two-factor law of 2026-08-20
(section 4) and docs/knowledge-base/retractions/psi-filter-seam-bug.md.

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
converge momentum and pressure against a force that cannot change. CORRECTED 2026-09-27:
`PSI_OUTER_CORRECTORS` defaults to yes since 5cbfaaa (2026-08-28).

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

**DECIDED (uniform yardstick landed, 276 steps to t = 0.003, both pairs read at the same
instant).** Relative np 1 / np 4 difference, uniform vs refined: mean|u'| 2.5e-6 vs
1.3e-5 (refined 5.4x), L2|u'| 2.6e-5 vs 2.3e-5 (0.86x), volume 9.3e-5 vs 8.1e-5 (0.87x),
shape 1.2e-2 vs 3.0e-4 (refined 40x SMALLER), Laplace jump 6e-9 vs 1.6e-8 (2.6x),
kErrL2Band 1.6e-3 vs 4.2e-3 (2.6x), A2hL2Band 2.6e-7 vs 1.7e-7 (0.66x); maxima over time
identical to four digits in both. Same order everywhere, nothing near the 10x line, one
metric better by 40x: the 1e-4 is the linear solver's decomposition dependence, present
identically on the hanging-node-free mesh, and the transition polyhedra add nothing
detectable. G1 PASSED on the refined mesh with the production force. (Aside worth
keeping: the UNIFORM pair's 1.2 % shape difference is the metric's own decomposition
sensitivity -- zeroSetRadialL2 is crossing-count weighted -- not the flow's.)

### DECIDED: static refinement around the interface is equivalent to the uniform mesh at matched h (2026-09-04, evening)

The production-force ladder is complete, four rungs, refined vs uniform at matched fine
spacing, capillary dt and step count on the 6R box (`stationaryDroplet3Drefined` vs
`stationaryDroplet3Duniform` + `uniform120`; `make_refined_equivalence_table.py`):

| N_fine (R/h) | cells refined / uniform | steps | peak L1 / L2 | L1 / L2 at T | volume | shape | Laplace jump | kErrL2Band | A2h | core-h refined / uniform |
|---|---|---|---|---|---|---|---|---|---|---|
| 60 (10) | 51 640 / 216 000 | 2302 | -0.18 % / -0.01 % | +3.9 % / +0.13 % | -0.24 % | +0.91 % | +1e-4 | -0.13 % | 0 | 1.2 / 4.1 (3.4x) |
| 76 (12.7) | 91 944 / 438 976 | 3281 | -0.20 % / -0.02 % | +3.8 % / +0.33 % | -0.37 % | -0.37 % | -1e-4 | +0.15 % | 0 | 2.7 / 10.6 (4.0x) |
| 96 (16) | 169 000 / 884 736 | 4659 | -0.17 % / -0.02 % | -2.4 % / -0.21 % | -0.11 % | +0.93 % | +1e-4 | -0.36 % | 0 | 7.0 / 25.1 (3.6x) |
| 120 (20) | 289 808 / 1 728 000 | 6511 | -0.20 % / -0.08 % | -0.79 % / -0.54 % | +1.7 % | +0.13 % | 0 | +0.11 % | 0 | 14.1 / 97.3 (6.9x) |
| 120, two levels | 125 448 / 1 728 000 | 6511 | -0.4 % / -0.09 % | -0.9 % / -0.6 % | +1.8 % | +0.65 % | +2e-4 | -0.6 % | 0 | 7.3 / 97.3 (**13.4x**) |
| 120, ball control (vs the band arm) | 296 752 / 289 808 | 6511 | +0.04 % / +0.02 % | +0.61 % / +0.43 % | -2.0 % | -0.03 % | +6e-5 | -0.16 % | 0 | 17.7 / 14.1 |

Pre-registered tolerances (10 % velocity L1/L2, 5 % volume/shape/kappa/A2h, 1 % Laplace jump)
met at every rung. Fitted orders at t = T over the four rungs, refined vs uniform
(`refined_convergence_orders.tex`): shape **2.08 / 2.08** (R 1.000), Laplace jump 2.30 / 2.30
(0.998), kErrL2Band 1.80 / 1.80 (0.999), A2hL2Band 1.99 / 1.99 (1.000), volume 2.45 / 2.48
(R 0.76, non-monotone), parasitic velocity at T **-1.03 / -1.12 (L1), -0.79 / -0.80 (L2)** --
never more than 0.09 apart (pre-registered +-0.3). The negative velocity order is the
uniform ladder's own: the end-of-run currents grow with refinement (the campaign's open
defect), and the refined mesh reproduces it as faithfully as the convergent metrics; the PEAK
currents converge at 2.2 (L1) / 2.15 (L2) on both meshes. The two-level control matches the
one-level N = 120 arm to 0.2 % (peaks) and <= 0.8 % (T) at half its cost. **The ball control
LANDED** (whole droplet interior fine, 296 752 cells, 6511 steps): against the interface-band
arm, peaks +0.04 % (L1) / +0.02 % (L2); at T L1 +0.61 %, L2 +0.43 %, shape -0.03 %, Laplace
jump +6e-5, kErrL2Band -0.16 %, A2hL2Band 5e-7; volume error -2.0 % at T and +5.7 % at T/2 --
the one nominal tolerance excursion of the campaign, on a value of 4.6e-8 (an absolute
difference of 3e-9), stated as such. The INNER transition, six fine layers inside the
interface, is as inert as the outer one. **The octree rejection of `stationaryDroplet3Dwide`
is retracted in full**: its hanging-node argument fell to the constant-curvature gate, its
equal-cost argument to 3.4-13.4x.

Where the CPU went: refined ladder 25 core-h + two-level 7 core-h against 137 core-h for the
uniform twins, on the same commit (`b72c7ef` stamps; drivers on disk `24223a2`, identical for
hex).

**P2 -- the POLYHEDRAL refined mesh is equivalent to its uniform twin too** (cluster, np 8 vs
np 16, N_CELLS 84 both, 3813 matched steps to t = 0.025; `polyDroplet3Drefined_r13p8` vs
`polyDroplet3D_r13p8`, interface cells 7.937e-5 in both, 150 362 vs 572 039 cells): peaks
-0.07 % (L1) / -0.014 % (L2); at T L1|u'| -0.09 % (7.485e-6 vs 7.491e-6), L2|u'| -0.014 %
(3.147e-5 vs 3.148e-5), volume +0.14 %, shape -0.50 %, Laplace jump -2e-4, kErrL2Band
+0.74 %, A2hL2Band -0.13 %; core-hours 6.4 vs 17.7 (**2.7x**; the poly ratio is smaller
than hex because cfMesh grades ~2x wider than asked, 11.8 fine layers for 6, and because
its 'uniform' twin already carries a half-size boundary layer). Same octree level at the
interface, same dt, same steps, same answer -- the graded cfMesh mesh is one more mesh to
the solver.

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

**P0 on the real box LANDED, and it CORRECTS the polyhedral ladder's quoted resolution.**
`polyDroplet3Drefined_r13p8` (maxCellSize 2.0e-4, one level): 88 929 -> 150 362 cells,
3.8x fewer than the twin's 572 039, interface cells at 7.937e-5 m, 5.9 complete fine
layers, `N_CELLS_suggested` 76 -- the band check failed the mesh rule on the config's 84,
as designed. The same tool on the uniform twin `polyDroplet3D_r13p8` finds ITS interface
cells at 7.937005e-5: IDENTICAL to seven digits (same octree level, same root cube). The
twin's N_CELLS 84 came from h_eff = (V/N)^(1/3) = 7.228e-5, an average pulled down by
cfMesh's boundary refinement -- its size profile reads `0-36: 1.00 | 36-72: 0.51-0.55`,
i.e. the polyhedral "uniform" meshes carry a half-size boundary layer. **Measured on all
three rungs, the interface cell size is EXACTLY 2^(-1/3) x maxCellSize** (the dual cell
holds half an octree box): 7.937e-5 / 5.556e-5 / 3.9685e-5 m, i.e. **R/h = 12.6 / 18.0 /
25.2 at the interface, not the 13.84 / 18.91 / 25.6 quoted from the average**; the h ratios
and hence the fitted orders are unchanged, the N_CELLS pins 84/114/156 (interface-honest
76/108/151) make the capillary dt 10 / 5.3 / 3.1 % conservative. Rung headers amended; the
rungs are kept as run. For a like-for-like P2 the refined arm runs at the twin's N_CELLS 84
(equal dt, equal 3813 steps) and DECLARES the deviation with the new inert token
`REFINE_ALLOW_PIN_MISMATCH` (the check still reports it). P1, the 4-rank smoke of the
graded polyhedral mesh (`polyDroplet3Drefined_r13p8_smoke4`, 458 steps), is running on the
laptop; the full study goes to the cluster on a PASS. **Its first attempt was a LAUNCH
FAILURE of my own making, not a result**: I had put `LD_PRELOAD=libjemalloc` (cfMesh's
glibc-2.39 workaround) in the study-global preamble, and the MPI solver segfaulted at
startup -- empty log, rc 139 -- exactly the trap CLAUDE.md records; the post-solve step
recorded it as "diverged". Voided; the second attempt lost its quotes in snakemake's
`--config` parsing (mesh rule failed before running anything); voided. Preload now a bare
path scoped to pMesh (`LEIA_PMESH_PRELOAD`). **P1 PASSED (third attempt, 458 steps, np 4):**
against the uniform twin at t = 0.003, mean|u'| 7.593e-6 vs 7.598e-6 (-0.07 %), L2|u'|
2.420e-5 vs 2.421e-5 (-0.02 %), volume +0.008 %, shape +0.30 %, Laplace jump +1.2e-8,
kErrL2Band +1.2 %, A2hL2Band +0.10 %, maxima over time within 0.09 % / 0.05 %. The full
study `polyDroplet3Drefined_r13p8` (np 8, 3813 steps, twin-matched N_CELLS 84) is
submitted to the cluster (orchestrator 54477977) -- **P2 PASSED, see the DECIDED block above.**

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
at np 8 against uniform at np 32 sits at exactly that level).

**N_fine = 76 pair LANDED (3281 matched steps): the same picture.** Maxima over time within
0.20 % (L1) and 0.019 % (L2); at t = 0.025 L1|u'| +3.8 %, L2|u'| +0.33 %, volume -0.37 %,
shape -0.37 %, Laplace jump -8.7e-5, kErrL2Band +0.15 %, A2hL2Band 2e-6; core-hours 2.66
vs 10.62 (4.0x cheaper; 91 944 vs 438 976 cells). The one metric that moves at T in both
rungs, L1|u'| by +3.8-3.9 %, is 2.2e-7 m/s against a peak of 1e-5 -- the decayed tail.

**N_fine = 96 pair LANDED (4659 matched steps): the same picture, and three rungs now fit.**
Maxima within 0.17 % (L1) and 0.024 % (L2); at t = 0.025 L1|u'| -2.4 %, L2|u'| -0.21 %,
volume -0.11 %, shape +0.93 %, Laplace jump +1.2e-4, kErrL2Band -0.36 %, A2h 4e-8;
core-hours 6.95 vs 25.09 (3.6x; 169 000 vs 884 736 cells). Provisional three-point orders
at t = 0.025, refined vs uniform: shape **2.05 vs 2.05** (R 0.9995 / 0.9999), Laplace jump
2.37 vs 2.37, kErrL2Band 1.75 vs 1.75, A2hL2Band 1.99 vs 1.99, volume 0.36 vs 0.36 (R 0.20,
non-monotone on both), and the parasitic velocity **-1.68 vs -1.81 (L1), -1.46 vs -1.47
(L2)** -- NEGATIVE on both meshes: at t = 0.025 the currents grow with refinement on the
uniform hex ladder too (the campaign's open defect), and the refined mesh reproduces it
order for order.

**uniform120 LANDED (1 728 000 cells, np 64, 6511 steps, 97.3 core-hours) and the TWO-LEVEL
control (30 -> 60 -> 120, 125 448 cells, 7.25 core-hours) reproduces it:** peaks within
0.38 % (L1) / 0.09 % (L2); at t = 0.025 L1|u'| -0.93 %, L2|u'| -0.62 %, volume +1.8 % (of
7.4e-8), shape +0.65 %, Laplace jump +2.4e-4, kErrL2Band -0.58 %, A2hL2Band 6e-6 --
**13.4x fewer core-hours** at R/h_fine = 20. The one-level N = 120 arm, the ball control,
r25p6 and the poly P2 study are still running; final orders when the one-level arm lands.

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
kErrL2Band 1.68 of 2000 (0.08 %), min|grad psi| 0.997.

**`r18p9` LANDED (6028 steps to t = 0.025, no blow-up) -- and the finer polyhedral rung is
WORSE in velocity and curvature while shape converges.** At the interface these rungs are
R/h = 12.6 and 18.0 (measured, above). At t = 0.025: mean|u'| 7.5e-6 -> 3.5e-5 m/s
(4.7x larger), L2|u'| 3.1e-5 -> 1.18e-4 (3.7x), kErrL2Band 1.68 -> 4.0 (2.4x), volume
error 2.6e-6 -> 3.4e-6; shape 4.85e-7 -> 2.46e-7 (order 1.9) and A2hL2Band 4.6e-7 ->
2.3e-7 (order 2.0) DO converge; Laplace jump 145.55 -> 145.53, min|grad psi| 0.997 ->
0.998. The time series says why: on r13p8 L2|u'| sits at 2.0-2.6e-5 for the whole run and
its L1 peaks (8.7e-6 at t = 0.02) then decays; on r18p9 both norms fall to a minimum near
t = 0.005 (L2 9.6e-6) and then GROW exponentially to the end -- e-folding rate of L2|u'|
over t = 0.0125-0.025: **29 1/s at R/h = 12.6, 124 1/s at R/h = 18.0** -- with kErrL2Band
climbing 0.67 -> 4.0 alongside. This is the campaign's open defect (section 4: the hex
ladders' growth rate rising 40 -> 52 -> 91 1/s under refinement, 'the 3D instability under
refinement'), now measured on polyhedra with the production chain and no filter.

**`r25p6` LANDED (2026-09-05, 9650 steps to t = 0.025, no blow-up) -- the prediction of a
faster growth was WRONG in its detail: the finest rung grows, but less than the middle one.**
The uniform polyhedral ladder at t = 0.025, h = the measured interface cell size:

| rung | R/h | steps | peak L1 / L2 | L1 / L2 (T) | volume | shape (m) | Laplace jump rel. | kErrL2Band | A2hL2Band | min\|grad psi\| |
|---|---|---|---|---|---|---|---|---|---|---|
| r13p8 | 12.6 | 3813 | 9.1e-6 / 3.1e-5 | 7.5e-6 / 3.1e-5 | 2.6e-6 | 4.85e-7 | 4.6e-4 | 1.68 | 4.64e-7 | 0.997 |
| r18p9 | 18.0 | 6028 | 3.5e-5 / 1.18e-4 | 3.5e-5 / 1.18e-4 | 3.4e-6 | 2.46e-7 | 3.6e-4 | 4.01 | 2.33e-7 | 0.998 |
| r25p6 | 25.2 | 9650 | 1.7e-5 / 5.9e-5 | 1.1e-5 / 4.4e-5 | 2.1e-6 | 1.72e-7 | 1.5e-4 | 2.42 | 1.19e-7 | 0.999 |

Equal-time L2|u'| (kErrL2Band in brackets), r13p8 | r18p9 | r25p6: t = 0.005: 2.2e-5 (1.35) |
9.6e-6 (0.79) | 4.3e-6 (0.39); t = 0.010: 2.1e-5 (1.34) | 1.7e-5 (1.05) | 1.1e-5 (0.43);
t = 0.015: 2.4e-5 (1.38) | 3.5e-5 (1.47) | 2.2e-5 (0.53); t = 0.020: 2.9e-5 (1.52) | 6.8e-5
(2.29) | 3.6e-5 (1.06); t = 0.025: 3.1e-5 (1.68) | 1.18e-4 (4.01) | 4.4e-5 (2.42). L2|u'|
e-folding rate over t = 0.0125-0.025: **29 / 124 / 82 1/s**; over r25p6's last 5 ms before
t = 0.0217 it was 209 1/s, and its L1 then PEAKED (1.7e-5) and fell to 1.1e-5 at T -- the
growth is not monotone in time either. Three-point orders at T against the interface h:
shape **+1.50** (R 0.987), A2hL2Band **+1.97** (R 1.000), Laplace jump ~1.6, volume +0.33
(R 0.44), min|grad psi| improving; the parasitic velocity **-0.61 / -0.53** (R -0.26, not a
power law: the middle rung is the worst) and kErrL2Band -0.54. So on polyhedra, as on hex:
everything geometric converges at 1.5-2, the Laplace jump converges, and the parasitic
current at a fixed physical time does not, because a resolution-driven growth sets in after
an initial decay, at a rate that is not even monotone in h (29 -> 124 -> 82). Every rung sits
far below blow-up (|u'| <= 1e-4 m/s at t = 0.025). Compare rungs at EQUAL time only. The
campaign's open defect stands, now with three polyhedral points and the growth-rate line.

### RUNNING: Popinet's translating droplet in 3D on polyhedral meshes (2026-09-05)

The 2D reproduction (`popinet2D_La12000_N64/128/256`, anchored to Popinet's 5 % at
R/h = 12.8) is taken to three dimensions on cfMesh polyhedra: `cases/popinetTranslating3D`
(2 x 1 x 1 box, D = 0.4, U = 1, equal densities and viscosities, La = 12000, We = 0.4, T = D/U
= 0.4; STL with named solids so the mesh has exactly `inlet`/`outlet`/`walls`), configs
`popinet3D_La12000_poly_r12p8/r19p2/r25p6` (R/h at the INTERFACE, N_CELLS pinned 64/96/128 via
the measured dual factor 2^(-1/3); ~0.7M / 2.4M / 5.6M cells; ~9 / 56 / 205 core-hours
expected). Pre-registered read-out: max over time of L1 and L2 of |u - U0|/U0 (Popinet's
convention), volume, shape, Laplace jump, kErrL2Band at T; compare rungs at equal time; a
rung that grows to blow-up is a result.

**The 4-rank smoke (`popinet3D_La12000_poly_smoke4`) BLEW UP at step 3, and the ladder is
NOT submitted.** The mesh is right (658 683 polyhedra, patches exactly inlet/outlet/walls,
interface cells 1.5624e-2 -> N_CELLS 64 exact to 0.01 %, checkMesh OK). Steps 1-2 are
healthy (continuity 1e-14, Courant 0.11); at step 3, in the FIRST outer iteration and
before any velocity solve, the phase volume jumps +5 % (0.01665 -> 0.01751) -- the
level-set transport moved the interface ~12x the physical U dt in one SL update -- then
the momentum residual is 1, the capillary force from the corrupted alpha/kappa drives
|u| to 3e4 m/s, Courant 3417, FPE at step 8. Not a launch artefact (preload scoped, log
complete). Two controls, both on four ranks, 78 steps:
- `popinet3D_La12000_hex_smoke4` -- the SAME 3D case on the hexahedral twin (524 288
  cells): COMPLETED, phase volume constant to 1e-8, L1|u-U0|/U0 7.3e-5, L2 5.4e-4,
  Laplace jump 10.044 (exact 10), kErrL2Band 0.71. The 3D translating setup is sound.
- `popinet3D_La12000_poly_smoke4_ccTrace` -- the polyhedral case with the cell-centred
  trace velocity: blew up identically (volume 0.0167 -> 0.0294 -> 0.21). The trace velocity
  is exonerated; the defect is polyhedra + two-phase flow at O(1) velocity (the kinematic
  polyhedral tests and the stationary polyhedral droplet, both healthy, never exercised
  that combination: one has no flow coupling, the other |u| ~ 1e-5).
A per-step field dump with three outer correctors shows only the aftermath (368k cells with
|dpsi| >> U dt, worst |dpsi| ~ 5 at cells 0.01-0.03 from the SIDE WALLS at the droplet's
streamwise station). Running now: the same dump with ONE outer corrector (localises the
first corruption) and the case with `SL_FIT householderQR` (the reconstruction's own
warning names near-degenerate stencils on the near-wall polyhedra; a normal-equations
quadratic fit can return huge coefficients there while passing the rank check).

**ROOT CAUSE FOUND (2026-09-05, morning): ill-conditioned quadratic stencils on cfMesh's
boundary-layer size-transition polyhedra.** Chain of measurements, each one run:
- QR fit: blows up IDENTICALLY (step-3 phase volume 0.017512193 vs 0.017512208) -> the
  fit's conditioning-in-arithmetic is not it. One outer corrector: same +5 % from a single
  SL update. ONE rank: identical (96 / 16 039 / 20 770 bad cells at t1/t2/t3) -> not a
  decomposition defect. The far-field re-extension belongs to another reconstruction class
  and never runs here.
- Against the exact translated signed distance: t0 exact everywhere (1e-8); after ONE step
  every interior cell is within U dt but 96 cells are wrong by up to 1.4 h -- ALL of them
  cfMesh boundary-layer cells (size 6.0e-3 = one octree level below the 1.56e-2 interior)
  in the OUTERMOST layer where they meet the interior, 0.024 from a wall, at arbitrary
  streamwise positions; 1 of 98 touches a processor patch. The error then multiplies ~15x
  per step through the stencils (2.2e-2 -> 0.36 -> 3.9), crosses zero in the boundary
  layer, the phase indicator paints a spurious interface along the WALLS (that is the +5 %
  phase volume -- the droplet is a bystander), and the capillary force on it blows the flow
  up.
- Offline replication from the mesh files: the offenders have **10 point-neighbours for the
  9-coefficient quadratic** (healthy interior cells 14-16), weighted normal matrices with
  condition numbers **3e7-6e12** (healthy ~1e2), quadratic coefficients **3e2-3e8**
  (healthy <= 1.5) while the fitted gradient is right to three digits. Such a fit is exact
  at d = 0 and wrong by coefficient x (U dt)^2 at the foot: 1e8 x 1e-20 = nothing on the
  stationary droplet (|u| ~ 1e-5), 7.5e3 x 6.5e-8 = 5e-4 per step at U = 1 offline, and
  worse in the code's own stencil. The kinematic polyhedral gates never saw it because
  their velocity vanishes on the walls, exactly where these cells live. The degeneracy test
  (Cholesky pivot <= 1e-15 ABSOLUTE, on entries of order h^4 ~ 6e-8) cannot see a
  condition number of 1e12; it caught 1 cell of 660 000.

**FIX (library, `uncachedQuadraticWeightedLeastSquaresReconstruction::build`, this
morning):** the quadratic fit's admissibility is decided ONCE per mesh from geometry
alone -- assemble the weighted normal matrix from stencil positions, Jacobi-scale it to
unit diagonal, take its smallest Cholesky pivot (~ 1/sqrt(condition)); below
`quadraticPivotTol` (default 1e-2 -- in the measured gap between the pathological ~1e-3 and
the healthy >= 2.5e-2 stencils) the cell uses the LINEAR
fit, whose coefficients are the first three of the same buffer (every consumer already
honours the per-cell coefficient count). A rank decision, not a filter, no shape or mesh
assumed; `quadraticPivotTol 0` reproduces the former behaviour exactly (the control arm),
`writeFitOrder true` writes the per-cell order as `slFitOrder`. Tokens
`SL_QUAD_PIVOT_TOL` / `SL_WRITE_FIT_ORDER` (default.parameter; entries in the two 3D
templates). This is a change to the LIBRARY, not to the refinement work: a defect of the
polyhedral reconstruction that any translating polyhedral case would hit.

**Gates for the fix (running):** hex bit-identity (`stationaryDroplet3DbitIdentity`, 65
steps, cmp against the pre-fix CSV -- the scaled pivots on hex stencils are ~0.1, so the
fallback must never trigger there), the polyhedral 4-step dump (phase volume constant,
boundary-layer errors ~1e-4, `slFitOrder` counting the fallback cells and where), the
polyhedral 78-step smoke, and the uniform polyhedral stationary smoke `polySmoke3D` against
its committed export (the boundary-layer garbage was 1e-12 there; the result must be the
same to round-off, else the polyhedral rungs are re-run). The Popinet 3D ladder is
submitted only after all four.

**Census of the fallback cells (`slFitOrder`, one-step serial runs on the fixed code):**
- hexRefined N_fine = 60 (51 640 cells, 1 688 hanging-node polyhedra): **0** fallback cells;
  the smallest scaled pivot over the whole mesh is 0.64. Hanging-node hex transitions are
  well-conditioned; the hex ladder is untouched by the fix.
- polyRefined r13p8 (150 362 cells): **3 555 fallback cells (2.4 %), every one in cfMesh's
  wall boundary layer** (within 1.3-1.8e-4 of the box, sizes 5.6-7.5e-5), the nearest
  14 fine cells from the interface, none within the band. The graded octree transition
  AROUND the droplet is well-conditioned; the polyhedral refinement results stand.
So the answer to 'is it only cells near walls' is: in every mesh this repository has, yes
-- cfMesh's boundary layer is the only size transition that starves the point-neighbour
stencil -- but the condition is topological, and the fix decides per cell from geometry
so it does not depend on that staying true.

**First gate pass with `quadraticPivotTol 1e-3` (recorded, then superseded).** The
polyhedral four-step dump: 29 255 of 658 683 cells (4.4 %) demoted to linear, phase
volume constant to 3e-9 (was +5 % at step 3), first-step boundary-layer error 2.2e-4 (was
2.2e-2, now at the interior level) -- but by step 3 two cells had grown to 2.4e-2 again,
and the 78-step smoke DIVERGED at step 12 (later than the original step 8, same route).
Offline replication of those cells: quadratic scaled pivots 1.04e-3 and 1.05e-3 -- a hair
ABOVE the tolerance, the same 10-neighbour boundary-layer stencils as the demoted ones;
their linear fits are perfectly conditioned (pivot 1.0). Every healthy cell measured has
a pivot >= 2.5e-2 (hex: 0.64). The pathological class straddles 1e-3 and there is a clean
gap up to 2.5e-2, so the default is set to **1e-2**, the count of quadratic cells within a
factor 3 above the tolerance is printed (the gap must stay visible), and a linear fit
that fails the same test cascades to the constant reconstruction. Rebuilt; the four gates
run again on that code (pipeline queued behind the 1e-3 runs).

**`polySmoke3D` is not a round-off yardstick -- its pre-fix reference was itself contaminated.**
Re-run on the 1e-3 code against the pre-fix run: L1|u'| +0.26 %, L2 +7.4 %, volume -1.1 %,
Laplace jump -2.7 %. That study is a coarse near-blow-up smoke (R/h 8.4 at the interface,
currents 0.1 m/s, volume error 4 %): its foot displacements were 2e-6, large enough for the
former garbage coefficients (1e8 x 4e-12 ~ 4e-4 per step) to have shaped the reference. The
committed `polySmoke3D_errors.csv` is therefore superseded by the fixed-code export. The
regression yardstick for the polyhedral rungs is the QUIET one: `polyDroplet3D_r13p8`
(|u| ~ 1e-5, displacement 6e-11, former contamination ~4e-13) re-run for 200 steps on the
final code (`polyDroplet3D_r13p8_regress200`) and compared to the cluster CSV row by row;
PASS = round-off (1e-8 relative). Queued behind the final gates.

**Gates on the 1e-2 + cascade code:** hex bit-identity cmp-IDENTICAL (0 demoted cells,
smallest pivot 0.64); polyhedral four-step dump CLEAN -- no cell beyond 0.1 h at any step,
phase volume constant to 3e-9, 34 514 cells demoted (5.2 %), none to constant; the only
growing errors are the eight box-corner cells, exactly +U dt per step (never updated; a
bounded nuisance, not an amplifier). `polySmoke3D` re-exported. **But the 78-step smoke
diverged again, at step 16** (8 -> 12 -> 16 as the tolerance went 0 -> 1e-3 -> 1e-2). The
dump says why: 3 840 kept-quadratic cells lie within a factor 3 above the tolerance
(pivots 1e-2..3e-2, condition ~1e3-1e4); their coefficients (~1e3) put ~6.5e-5 per step
into ~1 % of the cells (the count above 0.01 h grows 5 098 -> 11 873 over four steps), and a
quadratic cell of condition c amplifies a neighbour's error by ~c (d/h)^2 per step: at
c ~ 1e4 and d/h = 0.016 that is a gain of ~2.6, which turns the linear accumulation
exponential after ~8 steps -- step 16. So the admissible condition is set by the CFL,
c < (h/d)^2, not by the pathological class alone: at the scheme's design point d/h ~ 0.5
only stencils about as good as a hexahedral one (pivot ~0.5-0.6) are unconditionally
safe. The threshold is therefore taken from the PIVOT DISTRIBUTION of the meshes we run
(`slFitPivot`, written next to `slFitOrder`): interior polyhedra measured so far 0.7-0.9,
boundary-layer / transition cells 1e-3..0.12, hex 0.64 -- and the band cells' pivots are
checked on every mesh so that no second order is lost where the interface is.

**Pivot census (five meshes, one-step serial runs with `writeFitOrder true`;
`workflow/scripts/sl_fit_pivot_census.py`, table `sl_fit_pivot_census.csv` in the SL
article data). Classes: band |psi| < 6h, near 6-12h, far interior, far-small (cfMesh
boundary-layer / grading-transition cells).**

| mesh | cells | min band | min near | min far interior | far-small: n, min | demoted at 0.3 | cells in [0.3,0.4) |
|---|---|---|---|---|---|---|---|
| hex uniform N=60 | 216 000 | 0.789 | 0.789 | 0.637 | none | 0 | 0 |
| hex refined N=60, L=1 (1 688 hanging-node cells) | 51 640 | 0.773 | 0.769 | 0.637 | none | 0 | 0 |
| poly refined R/h 13.8, L=1 | 150 362 | 0.935 | 0.736 | 0.000 (grading transition) | 1 879, 3.3e-2 | 9 860 (6.6 %) | 35 |
| poly uniform R/h 13.8 | 572 039 | 0.949 | 0.949 | 0.615 | 129 888, 0 | 40 317 (7.1 %) | 40 |
| poly uniform Popinet-3D N=64 (feature-edge mesh) | 661 903 | 0.949 | 0.949 | 0.451 | 149 927, 0 | 46 849 (7.1 %) | 15 |

Reading: on every mesh the band cells sit at >= 0.74 (polyhedral >= 0.94), hex never
falls below 0.64 (hanging-node cells included), the polyhedral interior never below 0.50.
The degenerate stencils -- cfMesh's wall boundary layer on the uniform boxes, and on the
graded mesh ALSO the coarse-fine grading transition around the droplet (my size-based
"interior" class puts those there; they are 12h and more from the interface) -- form a tail
from 0 up to ~0.3, and the 0.3-0.4 bin is nearly empty on all three polyhedral meshes
(35 / 40 / 15 cells of 150k / 572k / 662k). Stability requires the admitted scaled condition
number c to satisfy c (U dt / h)^2 < 1: 1e-2 (c ~ 1e4) diverged at step 16 with
U dt / h = 0.016, so pivots up to 3e-2 amplify; **0.3 admits c <~ 100, is stable for
U dt / h < 0.1, sits in the gap, costs no second order within 12h of any interface, and
leaves every hex mesh bit-identical.** DECIDED: `quadraticPivotTol` default **0.3**
(code and `SL_QUAD_PIVOT_TOL`); the marginal-class report is `[tol, 1.5 tol)`.

**RETRACTED as a fix gate: `polyDroplet3D_r13p8_regress200` vs the cluster CSV.** The run
completed (200 steps, np 8, 31 196 cells demoted on the 1e-2 code) but the comparison is
confounded three ways: the cluster twin was MESHED by the cluster's own cfMesh build (a
different mesh to round-off), decomposed on 32 ranks (the linear-solver partition
dependence measured earlier), and run on the pre-fix code. Measured: 0.1-0.6 % on volume,
shape, curvature and A2h, 1e-2 on max|u'| (L_inf, not a verdict metric), 8e-9 on the
Laplace jump -- unattributable, so it decides nothing. Directory kept as
`polyDroplet3D_r13p8_regress200_tol1e-2_20260905` (a valid 200-step np-8 run of the
uniform polyhedral rung on the 1e-2 code, nothing more). The same-mesh, same-decomposition
gate replaces it: `polyDroplet3Drefined_r13p8_smoke4` (np 4, 458 steps, pre-fix CSV of
2026-09-04) re-run on ITS OWN mesh and decomposition with the 0.3 binary.

**Gates on the 0.3 code:** hex bit-identity `stationaryDroplet3DbitIdentity` (65 steps):
**cmp-IDENTICAL** to the pre-fix CSV, 0 cells demoted, smallest pivot 0.64.
Same-mesh, same-decomposition polyhedral regression (`polyDroplet3Drefined_r13p8_smoke4`,
np 4, 458 steps) -- FIRST ATTEMPT VOID, and the reason is a trap worth its own line: the
copy was taken from the finished study's `0/`, which is NOT the initial state. The solver
writes its projected initial pressure into time 0 of the processor directories at startup
and `reconstructPar -withZero` copies it back; the re-run therefore started from a converged
p_rgh (first solve: initial residual 6e-7 in 6 iterations, against 1 in 27 on 2026-09-04)
and differed from the study CSV from step 1 in the density/clipping diagnostics (verdict
metrics at round-off: volume 4e-15 absolute, shape and curvature 1e-10 relative, Laplace
jump identical). Attribution: `quadraticPivotTol 0` (the pre-fix algorithm) and 0.3 are
BIT-IDENTICAL on this mesh -- serial one-step fields (psi, U, p_rgh, kappa, alpha, phi)
identical to the byte, 20 steps on 4 ranks cmp-identical, two independent 4-rank launches
cmp-identical -- so the fix is inert here and the parallel path deterministic; only the
copied initial pressure differed. Corrected gate (0/ rebuilt from U, alpha, psi + p_rgh
from 0.org, same decomposition byte for byte): **cmp-IDENTICAL over all 458
steps** -- the 0.3 code reproduces the 2026-09-04 pre-fix study CSV byte for byte, first
pressure solve included (initial residual 1, 27 iterations). Two false starts on the way,
both in the INITIAL STATE and neither in the code: the finished `0/` holds the solver's
projected p_rgh AND its recomputed alpha.water (both rewritten into time 0 at startup and
copied back by `reconstructPar -withZero`; `psi` and `U` were not touched), so a re-run
copied from it differed from step 1 at 1e-15..1e-11 absolute. Regenerating `0/` from
`0.org` + `leiaSetFields -alphaName alpha.water`, as the workflow does, closed the gap
exactly. Pristine HEAD (worktree build, separate user dir) vs today's tolerance-0 path:
cmp-identical over 20 steps, so recompilation changed nothing either.
**The fix is bit-inert on every mesh it was measured on: hex (0 demoted cells) and the
refined polyhedral droplet (6.6 % demoted, none within 12h of the interface).** The
committed polyhedral `_errors.csv` exports therefore stand as produced.

Popinet-3D polyhedral 78-step smoke `popinet3D_La12000_poly_smoke4` (N = 64, np 4, 47 626
cells demoted = 7.2 %, none to constant): **COMPLETED all 78 steps** -- the first
polyhedral translating run to survive (0 -> step 8, 1e-3 -> 12, 1e-2 -> 16, 0.3 -> done).
Its vector at t = 0.02 (5 % of T) next to the hexahedral twin at equal steps: max over
time of mean|u'| 9.7e-5 (hex 7.3e-5), L2|u'| 5.3e-4 (hex 5.4e-4), max|u'| 9.5e-2 (hex
1.3e-2; L_inf, not a verdict metric -- to be located when the ladder lands), volume error
9.4e-6 (hex 2.5e-6), shape 9.2e-5 (hex 7.0e-5), centroid error 3.0e-5 (hex 5.4e-5),
Laplace jump 10.052 (hex 10.044, exact 10), band curvature error 0.81 (hex 0.71). L1/L2
parity with hex at matched N; the ladder decides convergence.
`polySmoke3D` re-run on 0.3 (17 273 demoted): max_t mean|u'| 1.200e-2 / L2 1.236e-1,
volume 4.44e-2, Laplace 150.3, kErr 2174 -- the four tolerance variants of this coarse
near-blow-up case (pre-fix 1.19e-2 / 1.246e-1 / 4.40e-2 / 153.7 / 2266; 1e-3; 1e-2) scatter
by a few percent; it is a smoke, not a yardstick, and its export is the 0.3 run's.

### RESULT, and a VOID: the polyhedral Popinet-3D MESH was defective -- four box edges were not feature edges (2026-09-05, night)

The user's single-phase question settled it. OpenFOAM's own `simpleFoam` (laminar, steady,
the Popinet BCs: inlet U = 1, slip walls, pressure outlet) on the pMesh used for every
polyhedral Popinet run does NOT keep the uniform stream: after 40 SIMPLE iterations from the
exact uniform state, |U - U0| = 8.3 % in the eight concave inlet-corner cells (our solver's
constant 0.0945), 3 % in the wall layer, 1 % in the interior, pressure -0.19..+0.22 where 0
is exact; a converged 302-iteration solve lands on the same perturbed field. A uniform stream
is an exact discrete solution whenever every wall face is perpendicular to it, so the mesh
was checked face by face: **1 984 of the 40 960 wall faces (4.8 %) have normals tilted into
the flow, up to 8 deg at the box corners.** Cause: `box2x1x1.stl` puts the four side walls in
ONE solid, so the four edges between them were never feature edges and cfMesh's Voronoi dual
wrapped faces around them; the stationary box `box6mm-edges.stl` (one solid per plane) has 0
tilted faces -- the stationary polyhedral ladder is unaffected. An explicit `boundaryLayers`
block does not cure it (8.6 %, p +-0.18). The cure is the feature-edge surface:
`surfaceFeatureEdges box2x1x1.stl box2x1x1.fms -angle 45` (12 edges) -> pMesh gives 0
tilted faces on walls, inlet and outlet, and `simpleFoam` then holds the uniform stream to
**|U - U0| <= 5.7e-16, p within 1e-14** after the same 40 iterations (its normalised residual
sits at 0.4 only because the exact solution makes the normalisation degenerate).

**VOID (wrong setup): every polyhedral Popinet-3D result on the STL mesh** -- the two ladder
rungs, the 78-step smokes and their exports, the field dumps, the corrector sweep's
translating arms, the census row "poly uniform Popinet-3D N = 64", and the two `inflowOnly`
horizon attempts. Directories renamed `_VOID_tiltedWallFaces_20260905` (laptop and cluster),
the committed `popinet3D_La12000_poly_smoke4_errors.csv` removed. What SURVIVES as method
findings, because each was measured on hex or on the stationary box as well: the quadratic-fit
admissibility (pivot 0.3, bit-inert everywhere), the boundary-face effect on the SL stencil
(8-50 % transport fractions on hex too; `inflowOnly` exact at the outlet in 2D), and "the
non-orthogonal correction converges in one pass" (stationary r13p8). What must be
RE-MEASURED on the corrected mesh: whether the outlet checkerboard amplifier exists at all --
it lived in the slab next to the tilted faces. The case now ships `box2x1x1.fms`
(`SURFACE_FILE` default in `popinetTranslating3D_poly.parameter`, all poly configs switched);
the 4-rank smoke on the corrected mesh is the gate, then the r12p8 horizon on the cluster.

**Gate on the corrected mesh (2026-09-07, 8 ranks, 78 steps): PASSED, and the polyhedral case
is now the BETTER one.** 0 tilted faces on walls/inlet/outlet (43 368 / 5 289 / 5 289 faces);
max_t max|u'| **8.6e-3** (tilted mesh 9.5e-2, hex 1.3e-2) and it sits at the droplet surface
(0.19 from the centre, R = 0.2) instead of the box corners; mean|u'| **4.1e-5** (hex 7.3e-5),
L2|u'| **2.8e-4** (hex 5.4e-4); interface metrics unchanged (volume 9.4e-6, shape 9.2e-5,
Laplace 10.052, band curvature 0.81); no cell above 0.03 U. The boundary-face transport
fractions are the same as before (inlet 0.082, outlet 0.240, walls 1.002): that defect is
the stencil's, not the mesh's, and stays selectable. r12p8 and r19p2 resubmitted on the
corrected mesh (ids in `.my_jobs`); the outlet amplifier question is answered by whether
r12p8 passes step 865.

**Corrected-mesh rungs (2026-09-07): BOTH DIVERGED, at steps 1154 (r12p8) and 1308 (r19p2) --
the same final steps as on the tilted mesh (1151 / 1301), but by a DIFFERENT route.** The
quiet phase is now better than hex (r12p8 to step 217: max|u'| 9.7e-3, mean 1.1e-4, L2 4.1e-4;
hex to step 815: 4.2e-2 / 5.6e-4 / 1.6e-3), then at step 267 (t = 0.069) max|u'| jumps
8.7e-3 -> 1.6 -> 9.5 in three steps while L2 stays 6e-4: ONE cell. It is the smallest cell of
the mesh -- 0.153 h on the box edge y = z = 0 at x = 1.749, 280x below a regular cell, the one
concave cell and the one face with volume ratio < 0.01 that `checkMesh` flags -- a sliver the
Voronoi dual makes where an edge point sits badly. Its velocity spike throws the SL foot
several cells, psi there reaches +41 / -38 (t = 0.2, alternating along the edge over x =
1.70-1.76, 65 cells with |psi| > 3, 1 939 cells with |u'| > 0.1), a fake interface forms, and
the contamination spreads until the droplet is destroyed (volume error 1.4e-4 at step 600,
1e-3 at 800, 0.15 at 1000). The coincident final step is the droplet's death by contamination,
not a geometric event. So: tilted faces (STL) gave the outlet-slab checkerboard; the
feature-edge dual gives one sliver on an edge; both are boundary-geometry artefacts of
cfMesh's pMesh, not of the level set. The extruded `boundaryLayers` on the `.fms` does NOT remove it (smallest cell 1.367e-8 against
1.378e-8; the concave cell stays). Sensitivity to the cell size (feature-edge mesh, `checkMesh
-allGeometry`): maxCellSize 0.019685 / 0.0190 / 0.0195 / 0.0200 give a smallest cell of 1.4e-8 /
6.7e-8 / 9.5e-8 / 3.8e-8 (h^3 = 3.8e-6, i.e. 0.15-0.29 h) -- small edge cells are systematic in
pMesh's dual at feature edges, but only 0.019685 also carries the inverted boundary point and
the volume-ratio-0.0095 face at x = 1.75 on the edge y = z = 0, exactly the seed. 0.0195 has
neither flag (smallest cell 0.29 h, non-orthogonality 58.9 deg). RUNNING (2026-09-07 night):
smoke at 0.0195 on the laptop as the 4-rank gate, then the r12p8 horizon on the cluster
(`popinet3D_La12000_poly_*_mcs0195.yaml`, pre-registered in the config). If it too is seeded on
an edge cell, the polyhedral Popinet benchmark needs a different boundary mesh (e.g. OpenFOAM
`polyDualMesh` of a tet mesh, or a hex boundary layer with a polyhedral core), not a cell size.

**0.0195 mesh (no flagged defects, smallest cell 0.29 h): gate PASSED (max|u'| 7.9e-3, L1 4.1e-5,
L2 2.7e-4 -- below hex), horizon DIVERGED at step 1026.** And this time the ORDER is measured:
the outlet-corner flow is perfect in the quiet phase (|u'| 3e-8 at t = 0.02; the field maximum
sits at the droplet), max|u'| grows only 0.5 %/step from 1e-2, the **fake zero set appears at
step 508 (t = 0.130) and the first velocity spike at step 530** -- the level set fails first,
the flow follows. At t = 0.2: 1 544 far cells with psi < 0, edge/corner psi error up to 117,
velocities of 12-14 U in the outlet-corner cells (x = 1.998, 0.33-0.35 h), while this mesh's
smallest cells (0.29 h, on the edge near x = 0.5-0.67) are quiet -- so not "the smallest cell"
but the outlet's small boundary cells. Three pMesh variants, three seeds, one class: **the
semi-Lagrangian far-field transport develops a growing alternating mode in cfMesh's small,
asymmetric boundary cells near the outlet (onsets 865 / 267 / 508), independent of wall tilt,
slivers, corrector count and cell size; the geometry sets where and when.** Pre-registered
verdict of the 0.0195 test: the polyhedral Popinet benchmark cannot be rescued by a cell size.
**Kinematics-before-dynamics control (SIGMA 0, 2026-09-08): the transport ALONE fails.** The
velocity stays an exact uniform stream for all 1563 steps (max|u'| 2.3e-14, L2 3.4e-15 -- the
flow solver on the feature-edge pMesh is perfect), yet the passive level set develops the fake
zero set at **step 508, the same step as with surface tension**, and by T = 0.4 the far field
is destroyed (outlet-layer |psi - psi_exact| 6.7e4, interior 9e3, 61 335 far cells across
zero, "volume error" 379 %). The capillary coupling played no part in the onset; it only
turned the level-set failure into a velocity blow-up afterwards. This is the missing
kinematic gate -- uniform translation THROUGH inlet and outlet on a polyhedral mesh -- and it
fails; the 3D polyhedral kinematic gates passed because u = 0 on their walls. Mechanism, now
without alternatives: a point-value semi-Lagrangian update with a least-squares fit amplifies
a checkerboard error at ~ (U dt / h) x stencil asymmetry per step; on cfMesh's small,
one-sided boundary cells at the outlet (0.33 h, local U dt / h = 0.05) that is ~0.5-1 %/step,
on symmetric hex stencils it cancels. The remedy every quasi-monotone SL scheme uses is the
clip of the reconstructed value to the stencil's bounds -- `clipToStencilBounds`, which this
repository already carries (token `SL_CLIP`) with the note "REQUIRED on general polyhedra"
from the KINEMATIC solver, while all three two-phase templates hardcoded it `false` (now the
token, inert default). **Measured in 2D N = 64 hex before adopting anything: the clip is NOT
inert at the interface** -- L1/L2 |u'| -6 % / -5 %, volume error +30 %, shape +5 %, centroid
+12 % -- so a global clip is a method change, and the clean form is a band-aware clip (far
field only), to be implemented only if the clip stabilises the far field at all. RUNNING
(2026-09-08): the SIGMA-0 passive case with the global clip (`..._sigma0_clip`, the scheme
test: PASS = no zero-set jump to T = 0.4) and the quiet stationary rung with the clip
(`polyDroplet3D_r13p8_clip`, 200 steps: what the clip does to a polyhedral band).

### RESULT: the Popinet-3D polyhedral ladder DIVERGES late, from an inlet/outlet level-set transport defect that is NOT the fit fix (2026-09-05, evening)

Two rungs landed, both DIVERGED, both on the right binary (`curvature-v2512`, `slFitPivot`
present), both with the expected demotion count (7 % of cells, wall boundary layer):

| rung | cells | ranks | demoted | blow-up | zero-set jump | max\|u'\| before onset |
|---|---|---|---|---|---|---|
| r12p8 (N 64) | 658 683 | 32 | 47 626 | step 1151, t = 0.29 | step 865, t = 0.222: 1.4e-4 -> 0.28 in ONE step | 0.09451 constant from step 20 |
| r19p2 (N 96) | 2.4M | 64 | 111 898 | step 1301, t = 0.18 | step 1002, t = 0.140: 5.7e-5 -> 0.34 in ONE step | 0.08833 constant from step 20 |

**Signature, identical on both:** the zero-set error jumps four orders of magnitude in a
single step while every velocity norm is still quiet (L2|u'| 1.5e-3 / 7e-4, growing only
0.1-0.16 % per step; max|u'| a CONSTANT boundary-cell value since step 20) and the Laplace
jump still reads 10.06 -- a fake piece of zero level set appears far from the droplet, its
radial error of 0.7-0.9 places it near the OUTLET (x ~ 1.8), the curvature there is garbage
(kErrL2 310 two rows later), the surface tension acts on it, and the velocity explodes 30-60
steps later. The finer mesh blows EARLIER IN TIME but at a similar STEP count: a per-step
amplifier, not a physical instability. The interface itself was at x = 0.72 / 0.64, more than
a diameter from every boundary.

**Where the level set was wrong before the onset (r12p8 field at t = 0.2, `psi_check`,
|psi - psi_exact| against the translated sphere):** bulk interior median 2e-5 (the transport
is fine); the first regular cells behind the outlet boundary layer (x = 1.97) at **1.73**
(psi 2.80 against 1.07, i.e. +0.12 h PER STEP, seven times the translation speed and in the
wrong direction) -- grown from 1.5e-3 at t = 0.02, an e-folding every ~100 steps; inlet
layer lagging by 0.16; corners frozen at 0.2 (= U t, the never-updated corner cells of the
4-step dump). No far cell was within h of zero at t = 0.2; twenty steps later one was.

**The defect underneath, measured on the two 78-step smokes (hex AND poly) at t = 0.02,
measured / exact change of psi in the first cell layer:**

| layer | polyhedral | hexahedral |
|---|---|---|
| side walls (u.n = 0) | 1.002 | 1.002 |
| inlet | **0.082** | **0.198** |
| outlet | **0.240** | **0.491** |

uniform across each layer to three digits. Tangential transport is exact; NORMAL to an
inflow/outflow boundary the level set moves at 8-50 % of the exact rate. Cause, in the
code: the CPC/CFC stencil (`extendedCentredCellToCellStencil`) treats every physical
boundary FACE as a data point (`cellToCellStencil::validBoundaryFaces` excludes only coupled
and empty patches; `slReconstruction::buildCentreTail` stores their centres, `collectData`
their patch values). With `psi zeroGradient` on inlet/outlet/walls each boundary face brings
the cell's OWN value at h/2 -- 0.19 h in cfMesh's small boundary cells -- with a 1/d^2 weight
four to thirty times a neighbour's, and the fit's boundary-normal gradient collapses toward
zero. The smaller the boundary cell, the closer the face, the worse: 8 % on polyhedra, 20 %
on hex. The polyhedral mesh then AMPLIFIES the boundary-layer error in the size-transition
cells behind it (the hex twin over the full horizon is running locally as the control:
does the boundary defect alone blow up, or does it need the transition?). This defect is
invisible on every closed-box case (walls, u.n = 0) and every kinematic gate (u = 0 on the
walls), exactly like the fit defect before it -- and the 2D Popinet run completes T = 0.4
with the same BC because a 2D hex mesh has no size transition.

**The fit fix stands:** same demotion counts as the laptop smoke, bit-inert on the quiet
cases, and the divergence route is a different one (far-field transport at the outlet, not
garbage coefficients: the transition cells there are linear at 0.3).

**Actions.** r25p6 (54482121) cancelled by id 1 h 15 min in -- it would only have measured
this defect -- and its 96-core solver child 54482215, which the driver's cancellation had
ORPHANED (53 minutes), cancelled by id once found in the driver's `.err` (lesson in
CLUSTER.md and the guides). r19p2 had already diverged. The ladder is VOID until the
boundary treatment is fixed; the study directories stay on the cluster as the record.

**Controls and the fix, measured (same evening).**
- **Hex twin over the FULL horizon (`popinet3D_La12000_hex_r12p8`, local np 8, 1563 steps):
  COMPLETED, no zero-set jump.** max_t mean|u'| 7.4e-4, L2 2.0e-3 (the 2D N = 64 reference:
  1.9e-3 / 4.4e-3), max 4.4e-2; final volume 1.4e-3, shape 3.0e-4, centroid 4.5e-4, Laplace
  jump 10.055 (exact 10), band curvature error 0.72. So the boundary defect alone is a
  bounded, linear error on hex; the AMPLIFIER that turns it into a blow-up is polyhedral --
  the size-transition cells behind cfMesh's outlet boundary layer.
- **`stencilBoundaryFaces exclude` (no boundary face as data) is NOT the fix: 2D N = 64
  DIVERGED at step 716.** The outlet became exact (error 0.000 at x >= 1.6 through t = 0.2)
  but the inlet -- whose departure foot lies OUTSIDE the domain -- extrapolated from the
  interior and grew an alternating-sign mode (+-0.1 at t = 0.15, +-1.55 at t = 0.20). With the
  faces INCLUDED the inlet layer is pinned by its zeroGradient value (transport fraction 0.009
  in 2D) and stable, and the outlet carries a stationary +-0.015 mode. The correct rule is per
  boundary face by the sign of its flux: an inflow or no-flow face IS data (the boundary
  condition speaks there), an outflow face is NOT (the foot lies inside, upwind).
- **`stencilBoundaryFaces inflowOnly` (implemented in `slReconstruction`, decided once from
  the flux at the first step so the stencil stays static for the pivot decision): 2D N = 64
  COMPLETED 1563 steps, outlet transport fraction 0.962 -> 1.000, inlet 0.009 (BC-pinned,
  unchanged), every droplet metric IDENTICAL to `include` to all printed digits** (mean
  1.8615e-3, L2 4.4332e-3, volume 1.213e-3, shape 8.152e-4, Laplace 5.0212): the far field
  never touches the droplet. Default stays `include` (bit-identity `A5` cmp-IDENTICAL on the
  final base class); token `SL_STENCIL_BOUNDARY_FACES` in the three Popinet/3D templates.
- **First `inflowOnly` horizon attempt on the polyhedral r12p8 case (laptop, np 8): DIVERGED
  at step 637, fake zero set at step 464 -- EARLIER than the unmodified stencil (865).** The
  outlet was fixed (transport fraction 0.240 -> 1.000 at t = 0.05 and 0.10) but the WALL and
  edge cells blew up instead (wall-layer error +1.05, edges +-0.67 at t = 0.10). Cause, measured
  on the written flux: slip walls carry a ROUND-OFF flux of either sign (1e-20, against 4e-4 on
  the outlet) and my outflow test `phi_f > 0` had no tolerance -- 9 719 of 40 960 wall faces
  (23.7 %) were classified outflow at random, so those cells lost part of their wall faces and
  got asymmetric one-sided stencils. Fixed: outflow = flux above 1e-6 of the largest boundary
  flux (a collective, called uniformly). Re-running the horizon; the 2D horizon re-checked
  first.
- **Non-orthogonal correction is sufficient -- measured (2026-09-05, evening).** The
  polyhedral meshes are strongly non-orthogonal (`checkMesh`: Popinet-3D N = 64 max 54.4 deg,
  average 4.2, skewness 1.7, 8 concave corner cells; stationary r13p8 max 61.8 deg, average
  5.7, skewness 2.0) and `nNonOrthogonalCorrectors` had been HARDCODED to 1 in every template
  (hex is orthogonal, so it never showed). Tokenised (`N_NON_ORTHOGONAL_CORRECTORS`, the
  materialiser's existing name, default 1) and swept 1 / 3 / 6 on the cluster (np 32):
  - stationary r13p8, 200 steps, all three arms valid (12 / 24 / 42 p_rgh solves per step,
    verified in the logs): mean|u'| differs by 5e-6 relative between 1 and 6 correctors, L2 by
    3e-7, the Laplace jump by 3e-9, band curvature by 8e-9, shape by 1e-9, volume error by
    5e-4 of a 2e-6 quantity. The correction converges in the first pass; the polyhedral
    parasitic current is NOT a Laplacian-correction artefact.
  - translating smoke, re-run with the stale arms removed (54482746), all three arms valid:
    1 / 3 / 6 correctors give max|u'| 9.4568e-2, mean 9.7169e-5, L2 5.3311e-4, volume
    9.434e-6, shape 9.224e-5, Laplace 10.0523 IDENTICAL to every printed digit, at 100 / 109 /
    119 s. The constant boundary-cell velocity defect (the 8 concave inlet-corner cells) is
    not a correction artefact either. Non-orthogonality is corrected sufficiently on both
    polyhedral cases; the count stays 1.
  - VOID: arms 00000/00001 of the translating sweep carried the DEFAULT tokens (U = 0.05,
    sigma = 0.07274, N 32, dt 7.2e-4) -- their `case_params.json` says commit bdd3f98,
    18:48:37: they were rendered by the FIRST, mis-specified driver (54482448, cancelled) before
    the cancel landed, and the corrected driver found the outputs present and did not re-render
    them (snakemake reuses existing outputs). A different case, dead at step 0 in the metric
    writer; not results. Lesson: after cancelling a driver, delete its study directory before
    resubmitting. Re-run with the stale arms removed (54482459 -> resubmitted).
- **Second `inflowOnly` horizon attempt (flux tolerance in place): DIVERGED at step 1301, fake
  zero set at step 839 (t = 0.215) -- the same step as the unmodified stencil (865 / 1151).**
  The outlet's first layer is now transported exactly (median fraction 1.000 at every write),
  yet the growing mode is STILL at the outlet: the size-transition cells at x = 1.97 behind
  cfMesh's boundary layer (error 0.08 at t = 0.10 -> 0.13 at 0.15 -> 0.66 at 0.20, an e-fold
  every ~200 steps) and the layer itself, where adjacent cells carry OPPOSITE errors (-0.226
  and +0.225 at (1.986, 0.648, 0.451/0.470)): an alternating-sign mode a minority of cells
  carries while the median fraction reads 1.00. Boundary-face data (this fix) and the
  non-orthogonal correction (the 1/3/6 sweep) are both excluded as its cause. What remains is
  the semi-Lagrangian update itself on systematically ASYMMETRIC stencils: the least-squares
  gradient of a checkerboard error cancels on a symmetric (hex) stencil and does not on the
  one-sided, size-graded stencils of the transition slab, and the update feeds it back at a
  rate of order (U dt/h) x asymmetry per step -- 0.016 x ~0.3 = the observed 0.5 %/step. The hex
  twin has no such mode (horizon completed). The slab cannot be removed by dictionary: cfMesh's
  Voronoi generator (`pMesh`) adds an explicit layer only when a `boundaryLayers` block is
  present, and ours has none -- the half-size boundary cells (0.38-0.45 h) and the size
  transition behind them are intrinsic to the dual construction. Discriminators that remain:
  a polyhedral mesh from a different generator (OpenFOAM `polyDualMesh` of a tet mesh, or
  cfMesh `tetMesh` itself, both without a graded boundary slab) to test "coherent stencil
  asymmetry"; and at the scheme level an update that cannot amplify a checkerboard -- a
  reconstruction that interpolates every stencil centre, or the quasi-monotone clip
  (`clipToStencilBounds`, present, off) restricted to the far field. Neither has run.
  **Answer to "can the 3D polyhedral Popinet droplet run to T = 0.4": not on this mesh with
  this scheme -- the far field, not the droplet, fails at step ~840 regardless of the outlet
  treatment.** `stencilBoundaryFaces` stays selectable, default `include`.
- Script: `workflow/scripts/sl_boundary_transport_check.py` (the numbers above).

#### The submission record (superseded by the result above)

`popinet3D_La12000_poly_r12p8` (N 64, ~0.66M cells, np 32) -> orchestrator **54482104**;
`popinet3D_La12000_poly_r19p2` (N 96, ~2.4M, np 64) -> **54482105**; `popinet3D_La12000_poly_r25p6`
(N 128, ~5.6M, np 96) -> **54482121**, submitted through `run-study-r25p6.sbatch` (a local copy of
the wrapper whose snakemake line carries the profile's `set-resources` list with the serial mesh
job at 32 GB / 240 min and pre/post-processing raised in proportion: the 3.6M-cell polyhedral mesh
step measured 11.7 GB and 71 min against the profile's 16 GB / 120 min; the user's pinned profile
is untouched). All three in `.my_jobs` (the ledger was cleaned of the submit plugin's stdout
chatter, see CLUSTER.md). Binary: `curvature-v2512` rebuilt 2026-09-05 at commit 0c5b9db,
`slFitPivot` symbol present, solvers relinked. Pre-registered read-out (`make_popinet_table.py`,
to be extended to 3D poly): max over time of L1 and L2 of |u - U0|/U0, volume, shape, Laplace
jump and band curvature error at T = 0.4; PASS = the L1/L2 disturbance decreases along the
ladder and sits at parity with the hex twin at matched N (the 78-step smoke: 9.7e-5 / 5.3e-4
vs 7.3e-5 / 5.4e-4). The L_inf column is reported, not judged.

### DECIDED: the quasi-monotone clip stops the polyhedral far-field failure (2026-09-08)

Two pre-registered gates ran on the cluster (jobs 54496854 and 54496855). Both COMPLETED.

1. **sigma = 0 passive transport, polyhedral mesh at maxCellSize 0.0195, N = 64, full horizon
   (1563 steps), clip ON against the clip OFF control of 2026-09-07.** The velocity field is
   an exact uniform stream in both runs (max|u'| 2.3e-14 and 2.5e-14), so only the transport
   is under test.

   | | clip OFF | clip ON |
   |---|---|---|
   | fake zero set | step 508 (t = 0.130) | NONE |
   | zero-set error at T | 8.66e-01 (max 1.009) | 3.98e-04 (max 3.99e-04) |
   | volume error at T | 3.79 (379 %) | 1.08e-06 |

   The two runs agree to every printed digit until step ~500 (zero-set error 9.1462e-05 at
   step 100 and 1.3261e-04 at step 400 in BOTH). The clip does nothing while the level set
   stays smooth. It acts only where the semi-Lagrangian update creates a new extremum, which
   is the checkerboard mode of the one-sided boundary stencils.

2. **Stationary polyhedral droplet R/h = 13.8, 200 steps, clip ON against the control.** Every
   metric is identical to seven digits: mean, L2 and max |u'|, volume error, shape error,
   centroid error, Laplace jump, band curvature error and A2h. The clip never fires, because
   the interface moves 1e-10 per step. This rung proves inertness ONLY for a quiet interface.

**What this does not settle.** The 2D hexahedral translating droplet measured the clip as NOT
inert at the interface (N = 64: L1 -6 %, L2 -5 %, volume error +30 %, shape +5 %, centroid
+12 %). A global clip is therefore a method change, not a switch to flip. The clean form is a
BAND-AWARE clip: apply it outside the narrow band only, where the level set carries no
interface. The interface metrics are then untouched by construction, and the form stays
mesh-agnostic, compact-stencil and MPI-decomposable. It is not implemented yet. `SL_CLIP`
stays runtime-selectable and default false.

### DECIDED: the Popinet cases run in SI units, and the conversion is exact (2026-09-08)

Popinet writes the benchmark in dimensionless form (rho = sigma = U = 1, D = 0.4 of the box
height). This repository ran it that way. Every dictionary therefore held a 0.4 m droplet
travelling at 1 m/s in a 2 m box, which is not a physical configuration. The cases now run the
DIMENSIONAL twin at the same dimensionless groups:

| quantity | SI value | fixed by |
|---|---|---|
| D | 1.0e-3 m (R = 5e-4) | chosen |
| rho, both phases | 1000 kg/m^3 | chosen (density ratio 1) |
| nu, both phases | 1.0e-6 m^2/s | chosen (viscosity ratio 1) |
| sigma | 0.012 N/m | La = sigma D/(rho nu^2) = 12000 |
| U | 0.0692820323 m/s | We = rho U^2 D/sigma = 0.4 |
| box | 5 x 2.5 x 2.5 mm | height = D/0.4, length 2 heights |
| T_U = D/U | 0.01443375673 s | one diameter of travel |
| dt at N = 64 | 9.23685e-06 s | 0.2323 of the Brackbill limit |

Five dimensional quantities carry two groups, so three are free. D, rho and nu are the choice
(water at 20 C); sigma and U follow. Re = sqrt(La We) = 69.28 and Oh = 9.13e-3 are unchanged,
so every number Popinet reports still applies.

**The similarity gate (the evidence the conversion is exact).** The SI 2D case at N = 64 was
run against the preserved dimensionless run
(`studies/popinet2D_La12000_N64_DIMENSIONLESS_20260908`). Both take 1563 steps. Normalised by
U and R, every metric agrees:

| metric | SI | dimensionless | difference |
|---|---|---|---|
| max_t L1\|u'\|/U | 1.861459e-03 | 1.861462e-03 | 0.0002 % |
| max_t L2\|u'\|/U | 4.433300e-03 | 4.433248e-03 | 0.0012 % |
| max_t Linf\|u'\|/U | 4.861838e-02 | 4.861816e-02 | 0.0004 % |
| max_t shape/R | 4.075464e-03 | 4.076008e-03 | 0.0133 % |
| final volume error | 1.212290e-03 | 1.213101e-03 | 0.0668 % |
| final centroid error/R | 4.190834e-03 | 4.192087e-03 | 0.0299 % |
| final pLaplace/(rho U^2) | 5.021250e+00 | 5.021236e+00 | 0.0003 % |

The worst difference is 0.067 %, which is round-off. **This also clears the absolute
linear-solver tolerances at the new scale** (p_rgh 1e-9, U 1e-8, psi 1e-12): the SI pressure
is 4.8 times larger and the SI velocity is 14.4 times smaller than in Popinet's units, and
nothing moved. The published 2D table therefore stands unchanged, in dimensionless form.

What changed: `cases/popinetTranslating{2D,3D,3D_poly}.parameter`, the Popinet block of
`cases/default.parameter`, both `blockMeshDict.template`s (POPINET_XLEN is now an ASPECT
RATIO and the box length is `#eval{ POPINET_XLEN * ylen }`; the 2D empty direction scales with
the height), all 20 `config/popinet*.yaml`, and the box surface
(`box5x2p5x2p5mm.stl` + `.fms`, made by `make_box_stl.py` and `surfaceFeatureEdges -angle
45`). `workflow/scripts/popinet_si.py` owns the set and checks a config against it:

    python3 workflow/scripts/popinet_si.py check config/popinet*.yaml

Two curation defects the conversion exposed, both fixed: `make_popinet_table.py` divided by a
hardcoded `U = 1` and `TU = 0.4` (it now reads TRANSLATION_SPEED, DROPLET_RADIUS and END_TIME
from `case_params.json`, and normalises the shape error by R); `aggregate.py` reported
`h = 1/N_CELLS` for every hexahedral study (it now reports `DOMAIN_LENGTH/N_CELLS`, identical
on the unit-length kinematic cases, and the true cell size on a dimensional one -- the
stationary-droplet tables carried 1/120 where the cell is 5e-5 m).

**The 3D gates, on 4 ranks locally (the gate before the cluster).** Both smokes COMPLETED at
78 steps and reproduce their preserved dimensionless twins:

| gate | worst normalised difference | note |
|---|---|---|
| 3D hex, 78 steps | 0.056 % (endpoint volume error) | velocity norms 0.005 % |
| 3D polyhedral (cfMesh, 674 493 cells), 78 steps | 0.33 % (band curvature) | velocity norms 0.014 % |

cfMesh built the same mesh at the new scale: 674 493 cells at maxCellSize 4.875e-05 m against
674 493 at 0.0195, band 2778 cells in both, interface cell 0.991 h_N in both, band check PASS
with a 0.96 % pin error. The mesher is not bit-reproducible under a coordinate scaling, so the
0.33 % on the polyhedral rung is mesh noise, not a unit effect.

**Absolute thresholds in the library, audited by hand.** The gate exercises the whole solver
path, so a scale-sensitive constant would have shown. A grep confirms why none did: every
literal in `src/leiaLevelSet` that a length, a velocity or a time is compared against is
either relative (`footPointTolRel` 1e-6, `1e-8*rTrust`, `1e-10*Uscale`, `1e-8/cbrt(V)`), or it
guards a dimensionless quantity (`alphaTol` 1e-8 on alpha, `trajectoryNormalEpsilon` 1e-12 on
|grad psi|, which is 1 for a signed distance), or it belongs to a model this case does not
select (velocity extension, SDPLS source, volume correction, band renormalisation). The
linear-solver tolerances ARE absolute, and the similarity gate clears them at this scale.

### MEASURED: the amplification bound of the fit, and why a smaller step cannot help (2026-09-08)

The question was whether the polyhedral failure is a Courant-number effect: the near-wall
cells are small, so their local Courant number is higher, and perhaps they need their own
step limit. Two measurements answer it.

**1. The local Courant number, from OpenFOAM's own `CourantNo` function object**, on the SI
smoke meshes (both N = 64, same dt):

| mesh | class | cells | median Co | max Co | min h/h_band |
|---|---|---|---|---|---|
| hex | whole mesh | 524 288 | 0.0164 | 0.0166 | 1.000 |
| poly | band | 27 027 | 0.0230 | 0.0231 | 1.000 |
| poly | far-interior | 495 118 | 0.0230 | 0.0593 | 0.702 |
| poly | far-small (cfMesh slab and edges) | 152 348 | 0.0431 | 0.1507 | 0.274 |

So the intuition is right in DIRECTION: the small cells carry 6.5 times the bulk Courant
number, and 23 % of the mesh is in that class. It is wrong in MAGNITUDE: the largest local
Courant number on the mesh is 0.15, and the time step is not Courant-controlled at all
(`adjustTimeStep no`; the step comes from the capillary law).

**2. The amplification bound of the reconstruction** -- the cheap static test, now in the
library as `fitProbeDisplacement` (inert at its default, writes `slFitAmplification`; see
workflow/README.md "The amplification bound of the fit"). The update is linear in the stencil
values, so the Lebesgue constant of the fit at the departure displacement,
`Lambda_c = |1 - sum_j g_j| + sum_j |g_j|`, bounds what one step can do to any mode.
Measured at `d = -U dt`:

| mesh | class | median | p99.9 | max | share > 1.10 |
|---|---|---|---|---|---|
| hex | band | 1.0164 | 1.0164 | 1.0164 | 0 |
| hex | far-interior | 1.0164 | 1.0430 | 1.0527 | 0 |
| poly | band | 1.0205 | 1.0205 | 1.0205 | 0 |
| poly | far-interior | 1.0205 | 1.0687 | 1.0963 | 0 |
| poly | far-small | 1.0329 | 1.2348 | 1.2608 | 0.39 % |

The bound is 5 times larger on the polyhedral mesh (1.26 against 1.05), and its whole excess
sits in the small one-sided cells: the five worst cells are 0.32-0.33 h, and the median
excess falls monotonically with cell size (0.039 below 0.4 h, 0.021 at h). Those are exactly
the cells where the sigma = 0 control grew its fake zero set. No cell has a singular normal
matrix, so the bound is defined everywhere.

**The decisive consequence.** `Lambda - 1` is PROPORTIONAL to the displacement, hence to the
local Courant number. At `|d|/U dt = 0.25 / 0.5 / 1 / 2 / 4` the maximum reads

    hex   1.0133  1.0265  1.0527  1.1041  1.2030
    poly  1.0684  1.1347  1.2608  1.4871  1.8413

so over a fixed physical time the accumulated bound is
`Lambda^(T/dt) = exp((Lambda - 1) T/dt) = exp(c U T)`, which does not contain dt. **A smaller
time step buys proportionally more steps at proportionally smaller growth per step, so it
cannot remove the failure.** `config/popinet3D_poly_sigma0_dtSweep.yaml` tests that
prediction directly: the sigma = 0 control at dt, dt/2 and dt/4, pre-registered to fail at
the SAME physical time (steps 508 / 1016 / 2032) if the analysis holds, and falsified if the
failure moves later as dt falls.

**What the bound does not say.** Lambda > 1 is necessary for growth, not sufficient: the
hexahedral mesh reaches 1.05 and runs the horizon, and the realised growth on the polyhedral
mesh (0.5 % per step) is a fiftieth of its bound. The bound proves the converse cleanly,
though: a reconstruction whose Lambda is 1 everywhere cannot create or grow a new extremum.
That is what a quasi-monotone clip enforces, and it is why the clip removed the failure
completely while leaving the first 500 steps identical.

### CONFIRMED: the polyhedral failure is not a time-step effect (2026-09-09)

The pre-registered sweep `config/popinet3D_poly_sigma0_dtSweep.yaml` ran on Lichtenberg as
job **54503136**. All three arms COMPLETED at np 32. Each arm is the sigma = 0 passive
transport control on the SI polyhedral mesh at N = 64. The arms differ in `FIXED_DELTA_T`
only: dt, dt/2 and dt/4, with dt = 9.23685e-06 s.

The read-out is the first step at which the far field grows a false zero set. Two columns
mark that step together: `zeroSetRadialL2/R` passes 1e-2, and `maxGradPsiBand` leaves 1.016.

| arm | dt [s] | steps to failure | failure time [s] | growth rate of zeroSetRadialL2/R [1/s] | the same rate per step [%] |
|---|---|---|---|---|---|
| dt | 9.23685e-06 | 528 | 4.877e-03 | 132.78 | 0.1227 |
| dt/2 | 4.618425e-06 | 1009 | 4.660e-03 | 128.59 | 0.0594 |
| dt/4 | 2.3092125e-06 | 1974 | 4.558e-03 | 126.05 | 0.0291 |
| reference run, Popinet's units, dt | 9.23685e-06 | 509 | 4.702e-03 | 132.94 | 0.1229 |

Read the table in this order:

1. **The number of steps to failure doubles and quadruples as the step falls.** The measured
   ratios are 1 : 1.91 : 3.74. The pre-registration predicted the steps 508 / 1016 / 2032 and
   the measurement gives 528 / 1009 / 1974.
2. **The failure TIME does not move.** The three times agree to 6.5 %, and the smallest step
   fails EARLIER, not later.
3. **The growth rate per unit physical time does not move.** It stays between 126 and 133 per
   second across a factor 4 in the step. The rate per STEP falls by that same factor.
4. **The pre-failure trajectories are one function of time.** `zeroSetRadialL2/R` at matched
   physical times agrees to 0.3 % in all three arms, from t = 1e-4 s to t = 4.5e-3 s.

This is what the amplification bound predicts. Lambda - 1 is proportional to the departure
displacement, hence to dt, so the bound over a fixed physical time is
`Lambda^(T/dt) = exp((Lambda - 1) T/dt) = exp(c U T)`, which does not contain dt.

**The falsifier did not occur.** The failure did not move later, and it occurred in every arm.
The near-wall cells therefore do NOT need their own step limit. A local time-step constraint
cannot repair this failure, and a global one only pays more steps for the same growth. The fix
must change the update operator, not the step.

**The same data confirms the SI conversion once more.** The reference run in Popinet's units
and the SI arm at the same step give the same growth rate (132.94 against 132.78 per second)
and the same failure time to 4 %. Their metric levels differ by exactly the length scale 400.

**A trap this measurement exposed.** `zeroSetRadialL2`, `zeroSetRadialLinf`, `centroidError`
and the `m2*` columns are ABSOLUTE LENGTHS in metres. In SI they are 400 times smaller than in
Popinet's units. A fixed threshold against them reads every SI run as clean, which is how the
first pass over this sweep misread all three arms as passing. Divide by `DROPLET_RADIUS`
before any threshold, any curated table and any comparison across unit systems. The committed
curation scripts use relative gates or divide by R already, so none of them carries the defect.

### DECIDED: the amplification defect is POLYHEDRAL CELLS, not cfMesh and not boundary layers (2026-09-09)

**The question.** Every polyhedral conclusion in this campaign rested on ONE mesher: the
workflow's `mesh: poly` is pMesh, and there was no cartesianMesh and no snappyHexMesh path. So
"is it only pMesh that is so bad?" was untestable by construction.

**The measurement, and it needs no coupled run.** The amplification bound
`Lambda_c = |1 - sum_j g_j| + sum_j |g_j|` depends on the stencil GEOMETRY alone, so one mesh
pass per mesher answers it. All four arms use the SAME box (5 x 2.5 x 2.5 mm), the SAME
nominal cell size and the SAME probe displacement `d = -U dt = -6.39948e-07 m`, which makes
them directly comparable to the recorded hex 1.0527 and pMesh 1.2608.

| mesher | cells | Lambda median | Lambda p99.9 | Lambda max | cells > 1.10 | cells demoted | smallest pivot |
|---|---|---|---|---|---|---|---|
| blockMesh (hex) | 524 288 | 1.0164 | 1.0430 | 1.0527 | 0 | 0 | 0.637 |
| cartesianMesh (cfMesh HEX) | 281 216 | 1.0131 | 1.0499 | 1.0549 | 0 | 0 | 0.629 |
| snappyHexMesh + 3 wall layers | 617 984 | 1.0164 | 1.0430 | 1.0566 | 0 | 0 | 0.410 |
| pMesh (cfMesh POLYHEDRAL) | 674 493 | 1.0205 | 1.0815 | **1.2608** | **588** | **47 622** | 1.4e-05 |

The excess above 1 is 0.053 / 0.055 / 0.057 / **0.261**. The three hexahedral-family meshers
cluster within 0.004 of each other; pMesh is 4.6 to 4.9 times worse.

**cartesianMesh is the decisive arm.** It is cfMesh's own HEXAHEDRAL mesher reading the SAME
`meshDict`, the SAME `box5x2p5x2p5mm.fms` and the SAME `maxCellSize 4.875e-05` as pMesh. It
lands at 1.0549, indistinguishable from plain blockMesh, with zero cells above 1.10 and zero
demoted. So the defect is not cfMesh, not the surface-based workflow and not the size grading:
it is the POLYHEDRAL cells.

**A RETRACTION of my own hypothesis.** I argued the danger was cfMesh's boundary-layer grading,
and predicted snappyHexMesh with layers would inherit it. The cell-size histogram refutes it:
cartesianMesh has 8 cells below 0.8 h and NONE below 0.6 h, snappyHexMesh with three wall
layers has 1536 cells in [0.4, 0.6) h and none below 0.4 h, while pMesh has 27 191 cells below
0.4 h and its five worst amplifiers are all 0.32-0.33 h. Layer addition on a hexahedral mesh
costs only +0.004 in the excess. The 0.33 h one-sided cells come from pMesh's polyhedral
generation.

**quadraticPivotTol and Lambda are essentially UNCORRELATED, and that closes the QR question.**
On pMesh the top 0.1 % of Lambda -- 674 cells -- are ALL still quadratic, and the worst
amplifier (cell 14736, Lambda 1.2608, size 0.333 h) has pivot **0.757**, six times the 0.3
tolerance and perfectly well conditioned. Conversely the 47 622 cells the pivot test DID demote
to a linear fit reach only Lambda 1.0511; they are harmless. So the pivot test measures
conditioning, Lambda measures amplification, and they select different cells.

That explains the 2026-09-05 result directly: `SL_FIT householderQR` blew up IDENTICALLY to the
normal equations (step-3 phase volume 0.017512193 against 0.017512208) because the amplifying
cells were never badly conditioned. There was nothing for better arithmetic to fix.

**Consequence for a truncated-SVD fit.** A classical truncated SVD keyed on the singular
spectrum would leave these cells untouched -- at pivot 0.76 there is no small singular value to
drop. Rank reduction still LOWERS Lambda (it shrinks the pseudo-inverse in the directions that
make `||g||_1` large, and the rank-1 limit is the weighted mean with every `g_j >= 0`,
`sum g = 1`, `Lambda = 1` exactly), and the headroom is 100 % of the top amplifiers. But the
TRIGGER must be Lambda, not the singular values: choose per cell the highest rank whose
`Lambda <= 1`. `Lambda <= 1` is the non-amplification condition, not a fitted coefficient.

**A trap this measurement produced, and it is the repository's own family.** The first
snappyHexMesh arm was VACUOUS and looked like a pass. `minVol 1e-13`, the OpenFOAM tutorial
default, is an ABSOLUTE volume in cubic metres sized for a metre-scale mesh; at h = 3.9e-05 m a
cell holds 5.9e-14 m^3, below it. snappyHexMesh flagged all 1 593 344 faces illegal, removed
every extrusion, reported "Added 0 out of 98304 cells (0%)" and EXITED 0 -- and the probe read
the untouched background blockMesh, reproducing its Lambda to eight digits (1.0527251). Only the
layer-addition line and the cell count reveal it. Verify the ARTEFACT, never the exit code, and
never a threshold that carries a length.

**Tools, all committed and regenerable.** `workflow/scripts/fit_amplification_probe.sh` builds
one mesh with any of the four meshers and runs the one-step probe;
`workflow/scripts/sl_fit_amplification_census.py` tabulates Lambda by cell class, by fit order
and by cell size, and writes
`docs/method-comparison/method-comparison-article/data/tables/sl_fit_amplification.csv`.
`cases/popinetTranslating3D/system/snappyHexMeshDict` is the layer-addition dictionary.

**What this means for the campaign.** A hexahedral-family unstructured workflow -- cfMesh
`cartesianMesh` or snappyHexMesh with boundary layers -- is clear of this defect. The clip
therefore leaves the critical path, and the open work for a spreading-droplet demonstration is
the two capabilities that do not exist yet: a wetting model (no `contactAngle` anywhere in
`src/`, `applications/`, `cases/`, `config/`) and, if dynamic refinement is wanted, cache
invalidation in `slReconstruction` on topology change (`haveCentres_`, `radius_`, `centreTail_`,
`keep_` and `phi_` are all sized or filled once, and `mesh.update()` is already called).

### G4 FALSIFIES THE CANDIDATE: the bound and the extremum exemption are in direct conflict (2026-09-09)

**The gate ran and its pre-registered falsifier occurred.** `popinet3D_poly_sigma0_clipGate`,
four arms, np 32, all COMPLETED at 650 steps on the polyhedral mesh, sigma = 0 so only the
transport is under test. Every arm ran the `curvature-v2512` binary (the `Exec` line is
checked). Job ids 54506624-27 and 54506640-43, in `.my_jobs`.

| arm | clip / keepExtrema | failure step | failure time [s] | zeroSet/R at T | volume err at T |
|---|---|---|---|---|---|
| 00000 | false / false | **527** | 4.8771e-03 | 3.4268e-01 | 1.0112e-04 |
| 00001 | false / true | **527** | 4.8771e-03 | 3.4268e-01 | 1.0112e-04 |
| 00002 | true / false | **NONE** | -- | 8.9732e-04 | 3.9322e-06 |
| 00003 | true / true | **506** | 4.6831e-03 | 1.2678e-01 | 2.2880e-05 |

Failure = the first step at which `zeroSetRadialL2/DROPLET_RADIUS` exceeds 1e-3, the fake zero
set in the far field.

**The controls hold, so the verdict is readable.** Arm 00000 reproduces the recorded failure
at step 528 to one step, which also proves the amplification diagnostic and the whole clip
code bit-inert on a POLYHEDRAL mesh. Arms 00000 and 00001 are BIT-IDENTICAL, so
`clipKeepExtrema` is unread when the clip is off. Arm 00002 re-establishes in SI that the
global clip removes the defect: no failure in 650 steps, the zero-set error suppressed by a
factor 380 and the volume error by 26.

**And the candidate fails at essentially the same step as no clip at all.** 506 against 527 is
4 %, inside the 5-38 % scatter this campaign documents for genuine instabilities. The extremum
exemption returns the growth the bound removed.

**Why, and it is not repairable by a better extremum test.** The exemption withheld 59.2 % of
the clip's firings on this mesh -- 2 680 917 cell-steps against 6 578 187 over 1950 corrector
calls. So **59 % of the cells the clip must bound on a polyhedral mesh ARE THEMSELVES stencil
extrema of psi.** That is the mechanism stated in the gate's own falsifier text: a growing
checkerboard has extrema at its own peaks, and any rule that exempts "a cell that is its
stencil's extremum" hands the defect exactly the cells it needs to grow in. The spurious
extremum and the genuine one are not separated by BEING an extremum.

The two requirements are therefore in direct conflict under this formulation:

| | hexahedral interface | polyhedral far field |
|---|---|---|
| clip, no exemption | +30.4 % volume error | defect REMOVED |
| clip + exemption | +1.03 % volume error | defect RETURNS at step 506 |

**RETRACTION of my own next step.** The repair I proposed and the author agreed to -- exempt a
cell whose own quadratic has a stationary point INSIDE the cell -- is now doubtful for a
measured reason, and it should not be built as stated. A checkerboard peak's own fit also has
a stationary point inside its cell, so that test would very likely withhold the same 59 % and
fail G4 the same way. The stationary-point test is a better DETECTOR of an extremum; the
measurement says the problem is not detecting extrema.

**What actually separates them is SCALE, and that is where to look next.** The level set's
genuine extrema are the medial axis of the signed distance -- the apex of the distance cone at
the droplet centre, the box corner farthest from the interface. They are wide: the apex sits
12.8 cells from the interface and the curvature of psi there is of order 1/R. A spurious
extremum is one cell wide, wavelength 2h. Two properties follow that a threshold-free rule
could use, and neither is tested yet:
- At a genuine extremum of a signed distance `|grad psi|` falls toward 0, because the field is
  smooth and stationary there. At a checkerboard peak `|grad psi|` is of order 1 or larger,
  because the oscillation steepens the gradient. `|grad psi| = 1` is the normalisation the
  signed distance already carries, so this is a statement about a dimensionless number, not a
  tuned length.
- A genuine extremum survives widening the stencil by one layer; a mesh-scale oscillation does
  not. The layer count is discrete geometry, not a fitted coefficient.

**The other honest option, and it must be stated.** Arm 00002 shows the global clip works on
polyhedra, and the hexahedral firing probe shows its whole cost on hex comes from SIX cells.
So a rule that names those six cells correctly would satisfy both requirements at once. The
above two properties are candidates for naming them. Until one is measured, the state is: the
clip removes the polyhedral defect and costs +30.4 % volume error on the hexahedral
translating droplet, and no formulation yet does both.

One caveat on arm 00002 that must not be lost: its zero-set error is still GROWING slowly --
4.46e-04 at step 1, 7.38e-04 at 500, 8.97e-04 at 649 -- and 8.97e-04 sits just below the 1e-3
threshold. The global clip suppresses the growth by a factor 380; it does not stop it. A longer
horizon may still cross. That is a measurement to make before the clip is called a fix, and it
is cheap: extend END_TIME on arm 00002 alone.

### RETRACTED AND REPLACED: the clip's damage is not the narrow band, it is the level set's own extrema (2026-09-09)

**What was believed.** The quasi-monotone clip stops the polyhedral far-field failure but is
not inert at the interface: on Popinet's 2D hexahedral translating droplet at N = 64 it costs
+30 % volume error, +12 % centroid and +5 % shape. The stated remedy was a BAND-AWARE clip --
apply it outside the sign-change narrow band only, where the level set carries no interface,
so the interface metrics are untouched by construction. That is section 2 of
[ROADMAP-poly3D.md](ROADMAP-poly3D.md), requirement 3, and it is now RETRACTED.

**What the measurement says.** The band-aware clip is implemented, runtime-selectable and
gated (`clipRegion outsideBand`, token `SL_CLIP_REGION`). It withholds the clip from 144 band
cells exactly as designed. It removes NONE of the damage. Every interface metric moved by
exactly what the global clip moved it, to every printed digit:

| | L1\|u'\| | L2\|u'\| | volume err | shape L2/R | centroid/R | kappa L2 band |
|---|---|---|---|---|---|---|
| baseline, clip off | 9.0723e-05 | 2.0419e-04 | 1.2123e-03 | 4.0755e-03 | 4.1908e-03 | 1.7498e+02 |
| clip, region all | +8.9 % | +3.3 % | +30.4 % | +4.9 % | +12.5 % | -3.1 % |
| clip, region outsideBand | +8.9 % | +3.3 % | +30.4 % | +4.9 % | +12.5 % | -3.1 % |

**The mechanism, located rather than argued.** The solver now writes `slClipFired`, so the
firing cells can be found. `config/popinet2D_clipFiringProbe.yaml` did that. The mesh has
exactly SIX cells where psi is the extremum of its own stencil -- the four box corners, which
are local maxima of the distance field, and the two cells at the droplet centre, which are the
apex of the distance cone and a local minimum. EVERY cell the clip fired in was one of those
six. Nothing else fired, at any step, and the `all` and `outsideBand` arms fired in the same
cells. The band was never where the damage came from.

A quasi-monotone bound CANNOT represent an extremum. Where psi_c is already the stencil
minimum, `lo` IS psi_c, so every reconstructed value below it is pulled back up to psi_c and
the extremum is flattened -- at every step, for the whole run. The fit undershoots at the apex
because a smooth quadratic cannot follow a non-differentiable minimum. So the clip fires on a
UNIFORM HEXAHEDRAL mesh, which the roadmap assumed it would not.

**The rule that follows, and it is classical.** Exempt a cell that is ALREADY its stencil's
extremum: the standard exemption of a monotonicity-preserving limiter at a smooth extremum.
It carries NO coefficient -- the test is `psi_c == lo` or `psi_c == hi`, exact in floating
point because lo and hi are taken over a stencil that contains psi_c. A flat stencil is not an
extremum, so the clip keeps enforcing exactness on a constant field. What stays bounded is the
SPURIOUS extremum, the one the fit creates in a cell that was not an extremum, which is the
polyhedral defect the clip exists for.

Measured, same case, same horizon (`config/popinet2D_clipRegionGate.yaml`, 8 arms, all
COMPLETED at 1563 steps): the exemption cuts the damage by a factor 30 and the clip's activity
from over 600 cell-steps to THREE.

| | L1\|u'\| | L2\|u'\| | volume err | shape L2/R | centroid/R | kappa L2 band | clip activity |
|---|---|---|---|---|---|---|---|
| clip, no exemption | +8.9 % | +3.3 % | +30.4 % | +4.9 % | +12.5 % | -3.1 % | 44 per step |
| clip + exemption | +0.3 % | +0.2 % | +1.03 % | +0.2 % | +0.4 % | +0.1 % | 3 cell-steps in total |

The band rule and the extremum rule are ORTHOGONAL and the matrix separates them: `all` and
`outsideBand` are identical in both rows, so on hex the band exclusion removes exactly zero
firings and earns nothing. Whether it earns its place must be argued on the polyhedral mesh or
the entry should be dropped rather than kept as a second knob.

**Inertness, gated.** All four clip-off arms of the 8-arm matrix are bit-identical to each
other, so both new tokens are inert. The G1 and G2 baselines were re-run at HEAD BEFORE any
clip code existed and came back bit-identical to the studies of 2026-09-08 -- which also
proves commit `7bf52eb`, the amplification diagnostic, bit-inert at its default
`fitProbeDisplacement (0 0 0)`. G0 shows the rendered case differing by the two token lines
only.

**The residual, characterised exactly, and NOT yet fixed.** The candidate is not bit-identical
on hex, so the roadmap's requirement 1 is not met. The residual is ONE cell: the droplet-centre
apex. The apex minimum is shared by a symmetric PAIR of cells, equal to the last digit at
t = 0. The exact equality test holds while that tie is exact and fails the moment a
rounding-level perturbation breaks it -- the survivor sits at
`(psi_c - lo)/(hi - lo) = 0.004`, a hair off the minimum, so it is clipped, and being clipped
pushes it further off. THREE cell-steps of clipping in the first steps then amplify to +1.03 %
volume error by step 1562: the candidate is identical to the baseline to ELEVEN digits through
step 100, its first byte difference is at step 20, and it grows from there. That amplification
of a far-field perturbation is the campaign's own central finding, arriving from a new
direction.

**A diagnosis that was tested and did NOT hold.** The full-horizon probe found the
pre-exemption survivors in a strip at `i = 126`, the single cell layer just inside the OUTLET,
every one a near-minimum at `(psi_c - lo)/(hi - lo) ~ 0.006`. `slReconstruction.H` records
that `stencilBoundaryFaces include` makes the outlet carry a stationary alternating error, and
`inflowOnly` is the gated remedy for it. It works on the FIRINGS and not on the DAMAGE:
`config/popinet2D_clipStencilGate.yaml` removed those eighteen firings -- the cumulative
counter goes to zero after the first steps -- and left the volume error at +1.027 %, against
+1.0 % with `include`. The `inflowOnly` baseline's own final metrics match the `include`
baseline to nine digits. So the outlet layer is not the residual, and that hypothesis is
retracted too.

**A process lesson, and it cost an hour.** A write-time SAMPLE of the per-step firing count is
not evidence. The `inflowOnly` arm reported 0 bounded cells at every one of its three write
times while its metric CSV had already left the baseline at step 20 -- the clip fires in a
burst in the first steps and then stops, and every write time sampled the quiet part. The
solver now accumulates the count every step and prints the total, the number of steps that
fired and the FIRST step that fired, and it writes `slClipFiredEver` so a firing between two
write times can still be located. Read the cumulative counter, never the per-step one.

**What is next, in one line: G4, and nothing else.**
`config/popinet3D_poly_sigma0_clipGate.yaml` asks whether the clip still REMOVES the
polyhedral far-field failure once the extrema are exempt. It can falsify the whole approach: a
growing checkerboard has extrema at its own peaks, so the exemption may hand the defect
exactly the cells it needs to grow in. Four arms, 650 steps, np 32, sigma = 0, at the first dt
of the sweep whose clip-off control failed at step 528. Do NOT build the apex detector before
that gate has run.

The apex detector, if G4 passes and requirement 1 still has to be met: exempt a cell whose OWN
quadratic has a stationary point INSIDE the cell. The gradient and the Hessian already exist
in the fit (the curvature path computes them), the test is whether the stationary point lies
within the cell, and the cell's own size is geometry rather than a tuned number. It is
coefficient-free, per cell, compact-stencil, and it is one implementation for 2D and 3D.

### The capillary time step on polyhedral meshes -- audited (2026-09-08)

The step follows `dt = CAPILLARY_DT_COEFF/nRef^1.5`, which reproduces 0.2323 of the Brackbill
limit `sqrt((rho1+rho2) h^3/(2 pi sigma))` at the spacing `h = DOMAIN_LENGTH/N_CELLS`. On a
polyhedral mesh `N_CELLS` is a PIN, so the question is which cell it describes. Measured with
`leia_refine.band_check` on the built meshes:

| case | interface cell | dt/Brackbill there | smallest cell anywhere | dt/Brackbill there |
|---|---|---|---|---|
| Popinet 3D poly, maxCellSize 0.0195, N = 64 | 0.991 h_N (all band cells equal) | 0.236 | 0.292 h_N (box edge) | 1.476 |
| stationary poly r13p8, N_CELLS 84 | 1.111 h_N | 0.198 | 0.247 h_N (boundary slab) | 1.889 |

Reading, in order:

1. **At the interface both cases are inside the limit**, at the intended safety factor. The
   capillary limit is a condition on cells that CARRY the capillary force, and the band cells
   of both meshes are uniform to a fraction of a percent.
2. **The far-field slivers are 1.5 to 1.9 times ABOVE the local limit.** They carry no surface
   tension, so no capillary wave is under-resolved there. It matters only for a case whose
   interface can reach a wall or a box edge; such a case must pin the step to the smallest
   BAND cell, not to the interface cell.
3. **The stationary polyhedral rung runs 16 % below its own law**: `N_CELLS = 84` describes
   h = 7.143e-05 m while cfMesh made 7.937e-05 m (the band check reports `pinRelError` 0.100,
   suggested pin 76). It is conservative, not wrong, but that rung is not dt-matched to its
   hexahedral twin. Re-pin to 76 before any dt-matched comparison.
4. The advective and viscous numbers are far from their limits everywhere: Co = 0.0165 and
   nu dt/h^2 = 0.0062 at the interface, Co = 0.056 in the smallest cell.

## 5. Lichtenberg — what is running

Login: `ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de`
Repo on the parallel file system: `/work/scratch/tm83tomy/leia`
Account `special00004`. Every job **must** set `--mem-per-cpu`.

| job | what | limit | output |
|---|---|---|---|
| `55048916` `leia-gate2D` `long` | **methodGate2D, FIXED** (2026-09-27 02:10, commit 1150e68, stamps g1150e68): the coupled-face density fix, the droplet-metric fix and `translatingSeamNp1` (section 11.14), all seven candidates, 56 arms. Artefacts checked before sbatch. The pre-fix run `55044205` and the translating re-run drivers `55047550`-`55047554` were stopped by id at 02:07, stragglers only; their studies are preserved as `methodGate2D_*_pre-20260927-*` and `methodGate2D_summary_pre-20260927-020856` (section 11.15). | 7 d (orchestrator) | `/work/scratch/tm83tomy/leia/studies/methodGate2D_*` |
| `54354379` `leia-curv` `long` | **interFoamDroplet2D** — re-running the `N` = 512 arm only (the other three are complete at the full 0.1 s horizon with `interFoam.csv` present). 106689 steps at the measured 1.52 steps/s = ~19.5 h. Replaces the arm lost to the cleanup bug above. | 28 h | `studies/interFoamDroplet2D/` |
| DONE | **filterOffAmplifier3D** 4/4, **upwindConvection2D** 8/8, **upwindConvection3D** 4/4, **filterThetaScaling3D** 6/6 — all analysed, section 4. | — | `studies/*/` |
| DONE | **stationaryDroplet3Dwide**, **cellCentreInverseFiltered512**, **domainSizeControl10R/6R/4R**, **psiOuterCorrectorsGain3D**, **ddtOrderGain3D** | — | `studies/*/` |
| `54477820` `leia-curv` `long` | **stationaryDroplet3DrefinedWB** GC (2026-09-04): the static-refinement mesh rule as a serial `case_pre` job + the constant-curvature well-balanced gate on the refined mesh, np 4, two arms x 921 steps. From `/work/scratch/tm83tomy/leia-curvature`. | 7 d (orchestrator) | `studies/stationaryDroplet3DrefinedWB/` |
| `54477891`-`54477895` `leia-curv` | **static-refinement ladder** (2026-09-04 evening, submitted after the hanging-node gate passed): orchestrators for `stationaryDroplet3Drefined` (fine N = 60/76/96/120, np 8), `stationaryDroplet3DrefinedL2`, `stationaryDroplet3DrefinedBall` (controls at 120), `stationaryDroplet3Duniform` (60/76/96, np 32) and `stationaryDroplet3Duniform120` (1.73M cells, np 64, ~90 core-h). Read with `make_refined_equivalence_table.py` when ALL arms of a pair are complete -- never before. | 7 d (orchestrators) | `studies/stationaryDroplet3D{refined*,uniform*}/` |
| DONE | **polyDroplet3D_r13p8 / r18p9 / r25p6** all complete to t = 0.025 (3813 / 6028 / 9650 steps; section 4) and the **static-refinement ladder + controls + polyDroplet3Drefined_r13p8** all complete (section 4, DECIDED). | — | `studies/polyDroplet3D_*`, `studies/stationaryDroplet3D{refined*,uniform*}` |
| (superseded) | **polyDroplet3D_r18p9** (1.46M cells, np 32) at t = 0.0171 / 0.025 and **polyDroplet3D_r25p6** (3.6M, np 64) at t = 0.0054 / 0.025; **polyDroplet3D_r13p8** DONE (3813 steps, section 4). | 23 h | `studies/polyDroplet3D_*/` |

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

## 9. SDPLS thread — sync state (2026-09-22)

The source-term line (SDPLS) is developed in a separate session from the curvature work
of sections 0-8. This section records the state of its three repositories and the actions
of 2026-09-22. Nothing in sections 0-8 was edited for it.

### 9.1 Where the three repositories stand

| location | at the start of 2026-09-22 | after this sitting |
|---|---|---|
| laptop `development` | 32a0d9b | the commits of 9.2, pushed |
| GitHub `origin/development` | 32a0d9b | same as the laptop |
| cluster `/work/scratch/tm83tomy/leia` (the SDPLS session's clone) | b3aa65e (2026-08-19), 221 commits behind | b16fc41 (fast-forwarded 2026-09-22 16:16), then the record commit of 9.5 |
| cluster `/work/scratch/tm83tomy/leia-curvature` (curvature session) | 81d04e7 | 32a0d9b (fast-forwarded 2026-09-22 11:20) |

The SDPLS code (`src/leiaLevelSet/sdplsSource/**`) has not changed since 4923800 (2026-08-26).
The SDPLS clone lacks five code commits of 2026-08-20 to 2026-08-26 (8d8467e, ef7341b, dc4d835,
e6ee5e1, 4923800), the configs `sdpls1Dstretch`, `sdplsRdivCoefficientError2Dvortex`,
`sdplsRdivProductRule2Dvortex` with `cases/1Dstretch`, and the article, deck and table commits
of 2026-08-20 to 2026-08-28 (a50d235, e346d50, 92a81ee). Its article is the 112 KB version;
upstream is 206 KB.

### 9.2 What this sitting changed on the laptop

- 246c96e: seventeen per-study result tables versioned under
  `docs/sdpls-level-set/sdpls-article/data/tables/` (13 copied byte for byte from the cluster
  clone, where they were untracked; 4 inputs of the generated convergence tables from the
  laptop's `studies/` tree), and one Reproducibility paragraph in `sdplsLevelSet.tex` that
  names them. PDF rebuilt (52 pages).
- 169172e: the negative-results deck caught up with the article: stacks 3 and 4 retracted in
  place (reversed-sign branch), five stacks added, 10 stacks, 26 slides, every number traced
  to an article line (245 numbers, 0 misses). Visual check in a browser pending.
- 7825401: the SDPLS planning briefs and the combined-source note (not implemented).
- the `workflow/README.md` SDPLS section; this section; the CLUSTER.md note of 9.3; the
  `run-studies.sbatch` hook of 9.4.

Laptop gate (before any commit): the six cases `1Dstretch_00000..00002` and `00009..00011` of
`studies/sdpls1Dstretch` (serial, run on this laptop at 093e3a4-dirty on 2026-08-26) were re-run
with the current binary (`libleiaLevelSet.so` of 2026-09-10, up to date with HEAD). Both
`leiaLevelSetFoam.csv` and `gradPsiError.csv` reproduce byte for byte: 12 comparisons with
`compare_metrics_csv.py --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`, 12 PASS.

### 9.3 Findings about the SDPLS clone, measured 2026-09-22

- **Four modified tracked tables** (`sdpls_convergence*.csv`, `convergence_orders*.tex`) are
  not data. The clone's `studies/sdplsConv3Dshear/sdplsConv3Dshear_errors.csv` carries an EMPTY
  `oscillation` column (an older `aggregate.py`), so `make_convergence_table.py` did not blank
  the t = T gradient columns of the 3D shear rows. Repair: re-aggregate that study in place
  (CLUSTER.md, "Aggregate in place"), regenerate, expect `git diff` empty.
- **Thirteen untracked tables** in `docs/sdpls-level-set/sdpls-article/data/tables/` are now
  tracked upstream with identical bytes (246c96e). Before the pull they are moved aside; after
  it they are compared with `cmp` and deleted.
- **Binary provenance.** `$HOME/OpenFOAM/tm83tomy-v2512/platforms/*/lib/libleiaLevelSet.so` was
  rebuilt on 2026-09-09 at 11:00:40 from the curvature clone (that clone pulled 0046961 at
  11:00:20); it contains the `compositeFlux` and `exponentialImplicit` symbols of 2026-08-26,
  while the SDPLS clone at b3aa65e has object files and build logs of 2026-08-19. Any SDPLS run
  launched from that clone after 2026-09-09 ran a library about 200 commits ahead of its own
  source. `sacct` for 2026-09-15 to 2026-09-22 shows no leia job; before trusting an SDPLS
  result of that period run
  `sacct -u tm83tomy -S 2026-09-09 -X --format=JobID,JobName%30,State,WorkDir%60 | grep '/leia$'`.
- **Modification-time sweep.** Between 10:59 and 11:12 on 2026-09-22 every file and directory
  under `/work/scratch/tm83tomy` received a fresh mtime, content unchanged where checked. File
  ages on that scratch mean nothing for that day: `foam_log_state.sh` `age=`, and Snakemake's
  `--rerun-triggers mtime`. Dry-run every study before launching it (9.5).

### 9.4 The pin, decided 2026-09-22

The SDPLS session pins its binaries to `$HOME/OpenFOAM/sdpls-v2512`. `run-studies.sbatch`
sources `$PWD/.leia_env` after `etc/bashrc` when that file exists, and does nothing otherwise.
The clone-local `/work/scratch/tm83tomy/leia/.leia_env` (unversioned; listed in
`.git/info/exclude`) exports `WM_PROJECT_USER_DIR=$HOME/OpenFOAM/sdpls-v2512` and the
`FOAM_USER_APPBIN`, `FOAM_USER_LIBBIN`, `PATH`, `LD_LIBRARY_PATH` overlay of SLURM.md 1.1. The
default dir `tm83tomy-v2512` is rebuilt from the pulled clone as a consistent fallback.

**MEASURED 2026-09-22 16:28 (job 54823248): the driver-side hook alone does NOT reach the
ranks.** A job submitted from the pinned driver shell, whose script runs the committed
profile's `env_preamble` (`module purge; module load ...; source etc/bashrc`), resolves
`leiaLevelSetFoam` and `FOAM_USER_LIBBIN` to `tm83tomy-v2512`: the preamble's `etc/bashrc`
resets `WM_PROJECT_USER_DIR` and puts the default dirs first again. So the SDPLS session's
solves run the DEFAULT dir today, which is consistent because 9.5 rebuilt it from b16fc41. The
pin becomes effective only when the job side sources the overlay too: one inert clause
`; [ -f "$PWD/.leia_env" ] && . "$PWD/.leia_env"` at the end of `profiles/slurm/config.yaml`
line 132. That line is the one the curvature clone edits locally, so committing the clause
makes that clone's next `git pull --ff-only` abort until its edit is stashed and re-applied
(or replaced by its own clone-local `.leia_env`). DECISION PENDING: commit the clause, with
that coordination, or stay on the default dir. After every launch, whichever way:
`grep -m1 '^Exec' studies/<study>/<case>/log.leiaLevelSetFoam` tells which dir ran.

### 9.5 Cluster procedure and verdicts (executed 2026-09-22 16:15-16:30 by this session, authorized)

Logs: `/work/scratch/tm83tomy/leia-sync-20260922/log/` (one file per step). Gate scratch:
`/work/scratch/tm83tomy/leia-gate-20260922/` (kept, small). Rollback anchor in the clone:
branch `sdpls-pre-pull-20260922` = b3aa65e. Job ids appended to `.my_jobs`: 54822859,
54823068, 54823248 (all `leia-sdpls-gate`, all COMPLETED). No SDPLS job ran between
2026-09-09 and 2026-09-22 (`sacct`), so no result was produced with the mismatched pair.

1. **Fast-forward** b3aa65e -> b16fc41 (228 commits). The 13 untracked tables were moved
   aside, compared with `cmp` after the pull (13 of 13 identical) and the aside copies
   removed. The clone has no modified tracked file left.
2. **The four regenerated tables.** Cause confirmed: `studies/sdplsConv3Dshear/
   sdplsConv3Dshear_errors.csv` had an empty `oscillation` column. Re-aggregated in place
   (`aggregate.build_database` over 36 cases; copies `*.pre-reaggregate-20260922` kept): the
   column is `on`, and the studies copy now carries four columns more than the tracked docs
   copy (`surfaceTensionForce`, `psiDdt`, `gradientErrorMax`, `gradientErrorMaxHalf`), every
   shared column identical. `make_convergence_table.py --method sdpls` then reproduced the
   committed tables exactly: `git diff` of the tables dir is EMPTY. Trap: the login shell's
   `python3` is miniconda 3.14 without numpy; the workflow's interpreter is the snakemake
   shebang `/shared/spapps_2026_01/gcc-11.5.0/python-3.11.14-*/bin/python3`.
3. **Builds** from the pulled clone: `sdpls-v2512` in 2 min 5 s (0 errors, 22 binaries,
   `compositeFlux` present), `tm83tomy-v2512` relinked in 8 s from the same objects;
   `curvature-v2512` untouched (2026-09-10 13:52).
4. **Gate job 54822859** (mpsd0001, 39 s, np = 4, both cases 129 steps, `Exec` names
   `sdpls-v2512` because the gate script sourced the overlay itself):
   - noSource control (`sdplsConv2Dvortex/2Dvortex_00000`, reference a2175d1-dirty of
     2026-08-12): `leiaLevelSetFoam.csv` identical at tolerance 0; `gradPsiError.csv`
     identical except the three strain-diagnostic columns `NARROW_MIN_R`, `NARROW_MEAN_R`,
     `NARROW_MAX_R`. `leiaSetFields.csv` identical.
   - R arm (`sdplsExpSource2Dvortex/2Dvortex_00000`, reference 40ecf59 of 2026-08-19):
     FAIL as pre-registered: 7 solver columns and 17 gradient columns differ, band gradient
     error L1 up to 1.2e-1 relative.
   - Attribution, not a regression of this sync: 30e6ba9 and 808d384 (2026-08-26,
     "setVelocity wrote FACE values into COUPLED patches: every parallel SDPLS run was
     wrong"; "The extension flux was the raw flux at every processor face") lie between the
     two references and the laptop reference 093e3a4 (2026-08-26), whose SERIAL 1D gate
     passed byte for byte (9.2). The control's signature -- transport identical, strain
     diagnostic changed -- is what a corrected `fvc::grad(U)` at processor patches produces.
     The 2026-08-19 parallel run is the defective one; the article already suspends those
     results (Section 7).
5. **Serial confirmation job 54823068** (np = 1, new binary, 9 s), serial against np = 4 of
   the SAME binary: noSource, every metric within 5.8e-6 (most below 1e-7). R arm: band
   gradient error L1 1.2e-3 and L2 2.4e-4, volume error 5.0e-4, shape 5.5e-6, whole-domain
   gradient extremes 0.94 (`MIN_MAG_GRAD_PSI`, serial 1.3e-3 against 2.2e-2 at t = 0.74) and
   0.19 (`MAX`). Against the OLD np = 4 reference the band L1 differs by 1.2e-1.
   **CORRECTED the same day (first written as "a residual decomposition dependence of the R
   path remains"; that reading was wrong).** Two measurements settle it. (a) Growth in time:
   at step 1 the serial-vs-np4 difference is 1.1e-10 in the R arm and 1.3e-10 in the
   sourceless arm, the same size and equal to the psi solver tolerance (PBiCGStab, DILU,
   1e-10); the sourceless arm stays at 1e-8 to 1e-9 for all 129 steps, the R arm grows to
   1e-7 at t = 0.17, 5.6e-3 (L2) at t = 0.42 and 3.6e-2 at t = 0.74. A seam operator would
   show at step 1; this is amplification of solver-level noise. (b) The single-variable test,
   job 54827753 (serial, same binary, same decomposition, only the psi tolerance changed from
   1e-10 to 1e-14): the R arm differs from its own tolerance-1e-10 run by 0.95 in
   `MIN_MAG_GRAD_PSI` (same step 107), 3.1e-2 in L2, 2.9e-2 in L1 and 5.7e-4 in volume, the
   same amounts as serial against np = 4; the sourceless arm by at most 4.1e-6. So the
   transport is seam-consistent, the R source path carries no remaining decomposition
   dependence, and the R arm at N = 32, T = 1 (non-reversed) amplifies solver-tolerance noise
   by seven orders over 130 steps, mostly where the far-field `|grad psi|` collapses toward
   zero (its minimum reaches 1e-3). Every `gradPsiErrorCSV` statistic is collective, so these
   are solution differences of that sensitive far field; the band metrics move 1e-3 and the
   shape 5e-6. Consequence for the SDPLS session: a seam gate on the R arm must tighten the
   psi tolerance (1e-14) or read the band metrics at a short horizon, otherwise solver noise
   masquerades as a seam defect; and the narrow-band parallel fixes of 2026-08-15 (4a73c71,
   6d144e9) are not involved, because this case's R path does not read the band
   (`GRAD_PSI fvc`, `MOLLIFIER none`) and the fixes are in both binaries compared.
6. **Pin probe job 54823248**: see 9.4. The ranks of a workflow-submitted job run the default
   dir. Consistent today; decision pending.
7. **Dry run** with the Makefile's flags (`--rerun-triggers mtime`): `sdplsConv2Dvortex` would
   re-run 6 of 42 cases, `sdplsExpSource2Dvortex` 2 of 14, because the 2026-09-22 mtime
   sweep gave some inputs a later second than their outputs, and `_solve` deletes the CSV
   first. Without `--rerun-triggers mtime` every case would re-run (params and code
   triggers). Before re-launching ANY existing SDPLS study on the cluster: `snakemake
   --workflow-profile profiles/slurm --configfile config/<study>.yaml --touch` (marks the
   outputs up to date without running), or preserve the study dir with a dated suffix.

### 9.6 Phone-runnable checks

```
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "squeue -u tm83tomy -n leia-studies"
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "tail -3 /work/scratch/tm83tomy/leia/.my_jobs"
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "cd /work/scratch/tm83tomy/leia && git status -sb | head -1"
ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de "grep -m1 '^Exec' /work/scratch/tm83tomy/leia/studies/<study>/<case>/log.leiaLevelSetFoam"
```

### 9.7 Not done, on purpose

- GitHub `main` is 440 commits behind `development` (one own commit, a `.gitignore` update).
- The briefs of 7825401 are not implemented (`interfaceDefects`, `sdplsCombined`).
- `benchVortexEulerT2/T8_errors.csv` tracked under the method-comparison theme are the
  2026-07-30 curation (15 rows, no commit column); the 2026-08-12 runs (30 rows, 6228f44-dirty)
  are not versioned. That theme's call.
- The deck's visual check in a browser.

## 10. Build policy migration: per-clone binaries, stamps, the library split (2026-09-22 ->)

Plan: `docs/plan-library-split-and-build-policy.md` (63ce3ed). Decisions in its section 0.
This section is the running record of its execution; one entry per work package.

### 10.1 WP1 -- `etc/leia-env.sh`: every clone runs its own binaries

**What changed (laptop, one commit).** `etc/leia-env.sh`, sourced AFTER OpenFOAM's
`etc/bashrc`, sets `WM_PROJECT_USER_DIR` to the clone root, so binaries install into
and run from `<clone>/platforms/<WM_OPTIONS>/{bin,lib}` (git-ignored), and it removes
every OpenFOAM USER dir from `PATH` and `LD_LIBRARY_PATH` while keeping OpenFOAM's own
installation and ThirdParty. Sourced by: `Allwmake`, `Allwclean`, the Snakefile shell
helper `sh()` (after the profile's or config's `env_preamble`, so every job, local or
SLURM), `run-studies.sbatch`, `run-sink-decay.sbatch`, `workflow/Snakefile.pressure-compatibility`,
and the seven probe and arm scripts under `workflow/scripts/` that launch a leia binary.
`run-studies.sbatch` now runs the clone it was submitted FROM (`cd $SLURM_SUBMIT_DIR`);
the job name is overridden on the command line (`sbatch -J leia-curv run-studies.sbatch`).
Not changed: `profiles/slurm/config.yaml` (the helper follows its preamble), the study
configs (their `[ -f $HOME/.leia_env ]` clauses become inert when that file is retired),
the legacy OpenFOAM-v2212 `cases/*/Allrun_*.sbatch`.

**MEASURED while building it.** The first path filter of the environment file removed
every entry matching `/OpenFOAM/<name>-v2512/platforms/`, which also matched OpenFOAM's
own `OpenFOAM-v2512` and `ThirdParty-v2512`. The link of every solver then failed with
undefined `Foam::UPstream::mpi_*` and `Foam::vtk::*` references: GNU ld resolves the
transitive dependencies of shared libraries through `LD_LIBRARY_PATH`. The filter now
keeps everything under `$WM_PROJECT_DIR` and `$WM_THIRD_PARTY_DIR`.

**Laptop gates, all PASS (2026-09-22 evening).**
1. Build: 22 binaries and 4 libraries in `<clone>/platforms/linux64GccDPInt32Opt/`;
   `which leiaLevelSetFoam` and `ldd` resolve to the clone; the account default
   `tmaric-v2512` is no longer on `PATH`; sourcing twice is idempotent (one entry).
2. Relocation: `sdpls1Dstretch` cases 00001 and 00010 re-run with the relocated binary
   reproduce the published `leiaLevelSetFoam.csv` and `gradPsiError.csv` at tolerance 0.
3. Through the workflow: `snakemake --workflow-profile profiles/local --configfile
   config/sdpls1Dstretch.yaml --config studies_dir=$PWD/tmp/gate-wp1 --until solve`
   (stopped before the report rule, so no data table is rewritten): the rendered job
   shell sources the environment file directly after the preamble; all 12 cases
   reproduce both CSVs at tolerance 0 (24 comparisons); `git status docs/` unchanged.

**Cluster (both clones), 2026-09-22 22:09-22:22, all PASS.** Logs under
`/work/scratch/tm83tomy/leia-sync-20260922/log/wp1_*`.
- Preconditions: no job of either clone running or pending.
- Curvature clone: its two private overlays (`profiles/slurm/config.yaml` line 132 and
  `run-studies.sbatch`: job name `leia-curv`, its own `cd`, the `curvature-v2512`
  exports) are in `stash@{0}` there, recoverable; then fast-forward to 019ed08; build
  2 min 2 s into `/work/scratch/tm83tomy/leia-curvature/platforms/`, 22 binaries,
  `which` and `ldd` resolve there. Probe job 54840842 (profile preamble, then the hook,
  as every workflow job does): `WM_PROJECT_USER_DIR` and `which
  leiaSemiLagrangianLevelSetTwoPhaseFoam` name the curvature clone. The curvature session
  submits with `sbatch -J leia-curv run-studies.sbatch` from its clone root from now on.
- SDPLS clone: the superseded clone-local `.leia_env` removed; fast-forward to 019ed08;
  build 8 s into `/work/scratch/tm83tomy/leia/platforms/`, 22 binaries.
  - Job 54840841 (np = 4, preamble then hook, `srun` as the profile does): `Exec`
    names `/work/scratch/tm83tomy/leia/platforms/linux64GccDPInt32Opt/bin/leiaLevelSetFoam
    -parallel`; both CSVs identical at tolerance 0 to the 16:22 gate output of the same
    code (`leia-gate-20260922/exp_00000`).
  - Job 54840840, the real launcher chain: `sbatch --export=ALL,STUDY=sdpls1Dstretch
    run-studies.sbatch` from the clone root, 5 min 1 s, 63 of 63 steps, 12 cases; all
    12 cases reproduce the laptop's published `leiaLevelSetFoam.csv` and
    `gradPsiError.csv` with max relative difference 0.0 (gcc 11.5 against gcc 13.3;
    this serial 1D case is bit-reproducible across the two machines). The report rule
    rewrote `docs/.../sdpls1Dstretch_errors.csv` with the new run's provenance columns
    only; restored with `git checkout`, the clone has 0 modified tracked files. The four
    regenerated theme tables stayed identical. `studies/sdpls1Dstretch` now exists on
    the cluster (12 cases, small).
- Old install dirs and `$HOME/.leia_env` untouched until WP6; `sdpls-v2512` and the
  default dir were last written 16:19 and 16:20, `curvature-v2512` 2026-09-10.

**Follow-up, MEASURED 2026-09-23 (WP3 gate).** The strip also removed the only `pMesh` a cluster
job could find: cfMesh lives in the dedicated install `$HOME/OpenFOAM/cfmesh/platforms/<opt>/`
(25 binaries and `libmeshLibrary.so`, 2026-09-01) and was on PATH only through the retired
`curvature-v2512` overlay; the first polyhedral mesh rule after WP1 failed with `pMesh: command
not found`. `etc/leia-env.sh` now appends `$LEIA_CFMESH_DIR` (default `$HOME/OpenFOAM/cfmesh`)
AFTER the clone's own directories. The laptop received the same layout from the cfMesh files of
`tmaric-v2512` (24 binaries plus the library, copied 2026-09-23), so `pMesh` resolves there too.

### 10.2 WP2 -- a version stamp in every library, printed in every banner and case

**What changed (laptop, one commit, 2026-09-23).** `etc/leia-stamp.sh <dir> <library>`
writes `<dir>/leiaStamp_<library>.C` (git-ignored) with one static registrar object
carrying `git describe --always`, plus `-dirty` when `src`, `applications`, `workflow`,
`cases` or `config` has uncommitted or untracked changes -- the same rule as the
`gitCommit` column of `materialize.py`. It rewrites the file only when the stamp changed,
so wmake recompiles exactly then; `Allwmake` calls it before `wmake src/leiaLevelSet`.
`src/leiaLevelSet/leiaVersionRegistry.{H,C}` collects the registered pairs; the five
solvers call `leia::reportVersions(Info)` and `leia::writeVersions(runTime)` right after
`createTime.H`, which prints `leia library <lib> : <stamp>` in the banner and writes the
same lines to `<case>/leia.version` (master rank, `runTime.globalPath()`).
`aggregate.py` reads that file into the column `libStamps` of `<study>_database.csv` and
`<study>_errors.csv`, next to `gitCommit`. A pulled-but-not-rebuilt clone then shows a
stamp that differs from `gitCommit` in every table and every log.

**Laptop gates, all PASS.** Build 1 min 6 s (the library recompiled for the new units);
the stamp file is generated and ignored; the stamp string is in the `.so`. Published case
`sdpls1Dstretch/1Dstretch_00001` re-run with the stamped binary: banner line present,
`leia.version` present, both CSVs identical at tolerance 0. `leiaTestSdplsSource`: 89
passed, 0 failed. Workflow run of `sdpls1Dstretch` up to the `aggregate` rule into a
scratch tree: `libStamps` present in the database and the errors table with the value
`libleiaLevelSet shared-method-config-2026-09-01-135-gf05d03c-dirty`, all 12 cases at
tolerance 0 against the published study, no data table touched. Note: on this laptop
the stamp reads `-dirty` as long as the curvature thread's two untracked scripts sit
under `workflow/`, exactly as `gitCommit` does; the stamp follows the same rule on
purpose. The mismatch test (clean binary, dirtied source, no rebuild) runs on the SDPLS
cluster clone, whose code paths are clean: see the cluster entry below.

**Cluster (both clones), 2026-09-23 11:09-11:13, all PASS.** Logs under
`/work/scratch/tm83tomy/leia-sync-20260922/log/wp2_*`. No job of either clone was running.
Both clones fast-forwarded to cc79df4; `Allwmake` generated the stamp and rebuilt in about
50 s each; both `.so` files carry `shared-method-config-2026-09-01-136-gcc79df4` -- no
`-dirty`, because both clones have clean code paths.
- SDPLS clone, np = 4 (job 54871894): `Exec` names the clone's `platforms/`; the banner
  prints `leia library libleiaLevelSet : shared-method-config-2026-09-01-136-gcc79df4`;
  `leia.version` lies in the case root only (not in the processor directories); both CSVs
  identical at tolerance 0 to the 2026-09-22 16:22 gate output of the same physics code.
- Curvature clone, serial (job 54871895): `which` resolves its own `platforms/`; the same
  banner line and case file; the 1D case `1Dstretch_00001` reproduces the SDPLS clone's run of
  2026-09-22 at tolerance 0 -- two clones, one commit, identical binaries in effect.
- Mismatch test in the SDPLS clone: an untracked probe file under `workflow/` makes
  `materialize._git_commit()` read `cc79df4-dirty` while the built binary's banner still reads
  `...-gcc79df4`; the stamp script would write `-dirty` on the next build. A pulled or
  edited clone that was not rebuilt is therefore visible in every log and, through
  `libStamps`, in every table. Probe removed; code paths clean again.


### 10.3 WP3 -- `libleiaLevelSet` split into `libleiaCore` and seven method libraries

**What changed (branch `wp3-library-split`, one commit, 2026-09-23 evening).** One `Make/`
per part under `src/leiaLevelSet/`: `libleiaCore` (`profile`, `narrowBand`, `phaseIndicator`,
`velocityModel`, the version registry), `libleiaSdplsSource`, `libleiaSemiLagrangian`,
`libleiaVelocityExtension` (links SemiLagrangian), `libleiaRedistancer`, `libleiaVolumeCorrection`,
`libleiaSurfaceTension` (with `fvOptions/semiImplicitCapillaryForce`, moved by `git mv` under
`surfaceTensionForce/`; no line of code changed), `libleiaAdvection` (links SdplsSource,
SemiLagrangian, VelocityExtension). `schemes/levelSetBlended` stays uncompiled, as it was in the
monolith. Every library links with `--no-undefined` and carries its own version stamp; `Allwmake`
builds the eight in dependency order, removes a stale `libleiaLevelSet.so` and the monolith's
leftover objects, and runs `etc/leia-check-deps.py` first, which prints the header includes
between the parts and refuses one that the link graph does not cover. Each solver links every
library it can select a model from (`README.md`, Build: the library table and the link matrix);
the 12 test applications and 2 utilities that linked the monolith link what their headers use.
Documents in the same commit: `CLAUDE.md`/`AGENTS.md` (build bullet), `README.md`,
`workflow/README.md` (Install), `CLUSTER.md` (stamp paragraph), `docs/IMPROVEMENTS.md` (the
ownership table names the library of every owned source; the index is now versioned), the
comments of `cases/sdplsSourceUnit/Allrun.sh`, `profiles/local8/config.yaml` and
`psiConservationCSV.H`. Not edited: `METHOD.md` line 35 still names `libleiaLevelSet` (thread C's
file; the owner decides).

**Laptop: build only.** The user's instruction of 2026-09-23 evening: no case runs on the laptop,
every rung runs on Lichtenberg. Build in a separate worktree (`/home/tmaric/OpenFOAM/repos/leia-wp3`,
its own `platforms/`, the main clone's binaries untouched): 2 min 10 s with 16 compile jobs; 8
libraries and 22 binaries; no link error at the first attempt, so the dependency table of the
plan's section 1 is complete; `etc/leia-check-deps.py` prints exactly its 9 edges and PASS; `ldd`
of the 22 binaries shows the link matrix and no `libleiaLevelSet`; every `.so` carries the stamp
`shared-method-config-2026-09-01-137-g5614bdc-dirty` (dirty: the uncommitted split itself).

**Gate (Lichtenberg, every rung, tolerance 0) -- pending.** Baselines run with the pre-split
binaries of the SDPLS clone (`cc79df4`, HEAD 5614bdc) into
`/work/scratch/tm83tomy/leia-wp3-20260923/baseline` (drivers 54887761, 54887767, 54887768, poly
arms 54887764-66); the split runs from its own clone `/work/scratch/tm83tomy/leia-wp3split` into
`.../split`; the polyhedral rung reuses the baseline meshes (`advect_bound_arm.sh`). The verdicts
follow in the next commit of this section.

**Gate (Lichtenberg, every rung, tolerance 0) -- PASS, 2026-09-23 19:22-22:23.** Scratch tree
`/work/scratch/tm83tomy/leia-wp3-20260923/{baseline,split}` (logs in `log/`, scripts `wp3_*.sbatch`,
`wp3_compare.py`); no case ran on the laptop (the user's instruction of that evening). Baseline:
the pre-split binaries of the SDPLS clone (`cc79df4`, HEAD 5614bdc), stamp
`shared-method-config-2026-09-01-136-gcc79df4`. Split: a third, temporary clone
`/work/scratch/tm83tomy/leia-wp3split` on branch `wp3-library-split` at a235d0a, built in 175 s,
stamp `...-138-ga235d0a` in all eight libraries. Both sets ran through the same driver into two
trees; every CSV in every case root was compared pairwise with `compare_metrics_csv.py --tol 0
--skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`. Rungs, all PASS:

| rung | cases (steps) | CSV pairs |
|---|---|---|
| unit `leiaTestSdplsSource` | 89 passed, 0 failed, in both trees | -- |
| exact 1D `sdpls1Dstretch`, serial | 12 (82/111/171/295) | 36 |
| seam `seamConsistency3Dserial`, `...par4` (np 4) | 2 + 2 (20) | 8 + 8 |
| coupled SL two-phase `stationaryDroplet3DbitIdentity`, serial | 1 (65) | 4 |
| `leiaLevelSetFoam` np 4: `sdplsExpSource2Dvortex` N 32 | 2 (both schemes) | 6 |
| coupled Eulerian two-phase `sdplsPsiBudgetDroplet2D` N 32, np 8 | 2 completed (1667) + 1 identical abort | 9 |
| `leiaRedistancedLevelSetFoam` np 4: `bulkVortexGRL` N 32, T 0.5, four redistancers | 4 (74) | 12 |
| hex 2D `advConv2Dtranslation`, np 8, 3 bounds x N 32-256 | 12 | 36 |
| hex 2D `advConv2Dvortex`, np 8, 3 bounds x N 32-256 | 12 (153/255/463/883) | 36 |
| hex 3D `advConv3DshearHex`, np 48, 3 bounds x N 32/64/128 | 9 (374/703/1362) | 27 |
| poly 3D `advConv3DshearPoly`, np 48, shared mesh per resolution: 49 911, 347 073, 2 389 233 cells | 3 bounds x 3 meshes (967/1989/3945) | 27 |
| loading: `ldd` of the 22 binaries in the split clone | the link matrix of `README.md`, no `libleiaLevelSet` | -- |

CORRECTED 2026-09-29: the `advConv2Dtranslation` row ran the REVERSED translation (11.19);
its 36 pairs are void, and the equivalence of the hex 2D rung rests on the vortex row.

209 CSV pairs compared, 0 differences. The polyhedral rung built each mesh ONCE (plain
`pMesh` jobs 54889779/82/85, digests 7e22acf1fbc9, 133ef0399684, 73ac6afbe5c8) and
`advect_bound_arm.sh` copied it into both trees' arms, so the cfMesh non-reproducibility never
enters the comparison. The split solvers print one stamp line per loaded library
(`leiaSemiLagrangeLevelSetFoam`: Core and SemiLagrangian; `leiaLevelSetFoam`: seven).

**Two findings about the current code, not about the split (both trees behave the same).**
1. `leiaLevelSetTwoPhaseFoam` writes its metrics to `leiaLevelSetFoam.csv`
   (`errorCalculation.H` line 37), but the workflow's solve rule waits for `<solver>.csv`, so every
   study of that solver ends in `MissingOutputException` after a completed run; `psiConservation.csv`,
   `capillaryFluxResidual.csv` and `leiaLevelSetFoam.csv` are written and were compared here.
2. The `Rdiv` arm of `sdplsPsiBudgetDroplet2D` aborts at step 1 in the coupled solver:
   `FOAM FATAL IO ERROR: Entry 'velocityModel' not found in dictionary "system/fvSolution"`
   (`sdplsRdiv.C` line 444, since ef7341b: the divergence form reads the prescribed velocity model,
   which a coupled case has not). The workflow records it as a diverged result. Both are for the
   SDPLS thread to decide; nothing in this work package changes them. **Both fixed 2026-09-24,
   on the user's instruction: see 10.4.**

**Three traps met on the way, all recorded.** (a) A command-line `--config` replaces the slurm
profile's `config:` list (CLUSTER.md, SLURM.md): the first baseline pass ran `mpirun` in one-task
steps and was voided (jobs cancelled by id, tree removed). (b) `pMesh` was on no job's PATH after
WP1 (10.1 follow-up above). (c) `foam_log_state.sh` printed "0 0" for a zero-step log (fixed,
4940440). The `advConv3DshearPoly` study lists each arm twice (`N_CELLS` 32 and 64 from
`cases/3Dshear_poly.parameter`); the gate meshed the N 32 case of each resolution.

**Merge and production clones (2026-09-23 23:20-23:37).** `development` fast-forwarded to 88e9480 on the laptop
(the split a235d0a plus four follow-up commits: the helper fix, the `--config` note, the
cfMesh environment line, the 10.1 follow-up) and pushed, then c6944e6 (`/agent-input/` in
`.gitignore`); both cluster clones pulled and rebuilt into their own `platforms/` (eight
libraries, no monolith, 22 binaries, `leia-check-deps` PASS, `pMesh` resolved; stamps
`...-142-g88e9480` in the SDPLS clone and `...-143-gc6944e6` in the curvature clone, the
commit each clone was built from), and each ran the probe set (unit, the 12 serial 1D
cases, the two np 4 cases) against the baseline tree: 42 CSV pairs each (12 serial 1D cases and the two np 4 cases, three CSVs per case)
identical at tolerance 0, unit 89 of 89 in both, no binary loads `libleiaLevelSet`
(drivers 54898114 and 54898205, trees `final-leia` and `final-leia-curvature`). The temporary split clone
was removed. Laptop: the main clone rebuilt after the merge; the worktree and the branch removed.

The scratch tree `/work/scratch/tm83tomy/leia-wp3-20260923` (baseline, split and probe results,
logs, driver scripts) stays as the gate's evidence until 2026-09-30; the temporary clone
`leia-wp3split` was removed on 2026-09-23. Next in the plan: WP4 (touch every finished study once,
then the launch guard) and WP6 (retire the old install folders after the first real runs).

### 10.4 WP4 -- study safety: touched once, guarded every time; two SDPLS fixes

**What changed (development, 2026-09-24 morning, three commits).**
1. `workflow/Snakefile`: the solve rule waits for `SOLVER_CSV`, the file the solver writes;
   for `leiaLevelSetTwoPhaseFoam` that is `leiaLevelSetFoam.csv`, for every other solver
   `<solver>.csv` as before. The dry-run DAG of every other study is unchanged by construction.
2. `src/leiaLevelSet/sdplsSource/sdplsRdiv.C`: the exact-reference diagnostic reads the
   `velocityModel` dictionary with `subOrEmptyDict`; a coupled case (no prescribed velocity
   model) now takes the documented skip instead of `FOAM FATAL IO ERROR` at the first assembly.
   Kinematic cases read the same dictionary as before.
3. `workflow/scripts/guard_finished_cases.py` (new, stdlib only) runs before every `snakemake`
   the Makefile starts (`GUARD` in `studies-one` and the seven study loops): it repeats the
   launch's dry run, and if a listed `solve` job targets a case whose log is COMPLETED or
   DIVERGED it prints the cases and exits 2. `LEIA_ALLOW_RERUN=1` turns that into a warning
   (preserve the directory first). Documents: `workflow/README.md` (Run and SDPLS sections),
   `SLURM.md` 4.2, `CLUSTER.md` ("Study safety").

**Laptop.** `Allwmake` rebuilds `libleiaSdplsSource` only (3 min 21 s in total); dry runs of
`sdplsPsiBudgetDroplet2D` (solve output `leiaLevelSetFoam.csv`, 45 jobs on an empty tree) and
`sdpls1Dstretch` (output unchanged) parse; the guard on `sdpls1Dstretch` against the laptop's
finished study: 0 solve jobs listed, exit 0. No case ran on the laptop.

**Cluster -- pending.** Both clones pull and rebuild; the guard is measured on
`sdplsConv2Dvortex` before the touch pass (expected: 6 of 42 cases refused) and after (exit 0);
`snakemake --touch` marks every study directory with a config up to date (measured 2026-09-24
on a scratch copy with one CSV removed: touch creates no file); the probe set (unit, 1D with
its Rdiv arms, np 4) is compared with the WP3 baseline tree at tolerance 0; and the N 32 arms of
`sdplsPsiBudgetDroplet2D` run through the fixed solve rule, the Rdiv arm included. Verdicts
follow in the next commit of this section.

**Cluster, 2026-09-24 08:39-08:58, all PASS** (logs `/work/scratch/tm83tomy/leia-wp3-20260923/log/wp4_*`,
`drv-54909721`, `drv-54909722`, `drv-54909812`). Both clones fast-forwarded to 8c36b90 and rebuilt
(only `libleiaSdplsSource` recompiled, 35 s; stamp `shared-method-config-2026-09-01-147-g8c36b90`).
- Guard, before the touch pass, SDPLS clone, `sdplsConv2Dvortex`: 6 of the 6 listed `solve` jobs
  targeted finished cases (`2Dvortex_00036` to `00041`, COMPLETED, 700 steps, CSV present); exit 2.
  That is the 2026-09-22 hazard, refused. After the touch pass: 0 `solve` jobs listed, exit 0.
- Touch pass (`snakemake ... --touch`, one dry run before and after each study): SDPLS clone 43
  studies touched, 3 already up to date, 0 dry-run failures; curvature clone 98 touched, 1 up to
  date, 0 failures. `git status docs/` is empty in the SDPLS clone; the curvature clone's 27 lines
  (5 modified and 22 untracked tables, every mtime 2026-09-22 10:59) are the curvature session's
  own state from before the pass (the 5 modified tables are the same 5 it has modified on the
  laptop). Code paths clean in both. Studies that still list jobs after the touch, all for
  outputs that never existed: 18 in the SDPLS clone, 6 in the curvature clone -- mostly
  `reconstruct`, `aggregate`, `report` of pipelines that never ran to the end, plus `solve` jobs
  over cases without a CSV (`benchVortexGRLfrozen` 2, `benchVortexVET2` 2,
  `oscillatingDropletFootEvalFace` 6, `stationaryDropletFootEvalFace` 5, `stationaryDropletRenorm`
  1; curvature clone `advConv3DshearHex` 1, `filterThetaScaling3D` 1). The guard permits those
  unless the case log says COMPLETED or DIVERGED.
- Probe sets in both clones (unit, the 12 serial 1D cases with their Rdiv arms, the two np 4
  cases): 42 CSV pairs each identical to the WP3 baseline at tolerance 0, unit 89 of 89; the 1D
  Rdiv arms print `sdplsRdiv exactRef: SKIPPED (... velocityModel 'uniaxialStrain' ...)` as before.
- Coupled Eulerian study `sdplsPsiBudgetDroplet2D`, N 32 arms through the fixed solve rule, np 8:
  all three cases COMPLETED with 1667 steps and no `MissingOutputException`; `noSource` and `R`
  identical to the WP3 baseline (8 CSV pairs); `Rdiv`, which aborted at step 1 before, runs to the
  end and prints `sdplsRdiv exactRef: SKIPPED (implicitSurface 'implicitSphere', velocityModel
  'none', parRun true)`; `aggregate.build_database` over the three cases gives 3 rows with 160
  populated columns each and the eight stamps in `libStamps`. What the Rdiv arm's numbers mean is
  for the SDPLS thread to read; this gate shows only that the arm runs and aggregates.
- The first attempt at that gate targeted the old CSV name through my scratch target selector
  (nothing ran); the second attempt, with the mapped name, is the one recorded above.

WP4 is closed. Open in the plan: WP6 (rename `$HOME/OpenFOAM/tm83tomy-v2512`, `curvature-v2512`,
`sdpls-v2512` and `$HOME/.leia_env` with a dated suffix after the first real runs from the new
folders; delete a week later; the laptop's `tmaric-v2512` keeps only its cfMesh copy, now also in
`$HOME/OpenFOAM/cfmesh`).

### 10.5 WP6 -- the old install folders retired (2026-09-24), deletion due 2026-10-01

**Condition met first.** Since WP1 every case that ran from either clone -- the WP2 probes, the
full WP3 gate (both trees), the WP4 gates and probes -- printed the clone's `platforms/` path
in its `Exec` line and the clone's stamps in its banner; the launcher chain itself
(`sbatch run-studies.sbatch` -> `make studies-one` -> guard -> snakemake) ran from each clone
today (jobs 54916369, 54916370: guard 0 solve jobs, "Nothing to be done", `ALL-STUDIES-DONE`).

**What was retired, Lichtenberg, 2026-09-24 12:27 (log `wp6_retire.log`).**
- `$HOME/OpenFOAM/tm83tomy-v2512` STAYS: its `platforms/linux64GccDPInt32Opt/lib` also holds
  `libweiFunctionObjects.so`, a library of another project on the shared account (its jobs
  `openfoam-bayesopt-cht/...` were running). Only the leia files moved, into
  `platforms/linux64GccDPInt32Opt/retired-leia-2026-09-24/{bin,lib}` (22 binaries, 4 libraries).
- `curvature-v2512` -> `curvature-v2512.retired-2026-09-24` (25 MB; its cfMesh copy is a duplicate
  of `$HOME/OpenFOAM/cfmesh`), `sdpls-v2512` -> `sdpls-v2512.retired-2026-09-24` (18 MB),
  `~/.leia_env` -> `~/.leia_env.retired-2026-09-24` (the `[ -f $HOME/.leia_env ]` clauses in the
  study configs are inert now). A plain OpenFOAM shell finds no `leiaLevelSetFoam` any more; with a
  clone's `etc/leia-env.sh` it finds the clone's, and `pMesh` in `$HOME/OpenFOAM/cfmesh`.
- Laptop, the same rule: `$HOME/OpenFOAM/tmaric-v2512` stays (36 binaries and 18 libraries of
  other projects, plus cfMesh and `run/`); its 22 leia binaries and 4 leia libraries moved into
  `platforms/linux64GccDPInt32Opt/retired-leia-2026-09-24/`. A plain shell finds no leia binary.
- Found on the way: `run-sink-decay.sbatch` sourced `$HOME/OpenFOAM/curvature-v2512/etc/bashrc`,
  a file that never existed, so it sourced nothing and stopped at `etc/leia-env.sh`; it now
  sources OpenFOAM's `etc/bashrc` (this commit).

**Verification after the retirement.** Probe drivers 54916367 (SDPLS clone) and
54916368 (curvature clone), 12:27-12:36: unit 89 of 89, `ldd` shows no monolith, the 12 serial 1D
cases and the two np 4 cases give 42 CSV pairs each identical to the WP3 baseline at tolerance 0;
every `Exec` line names the clone's `platforms/`, every banner its seven stamps, and no retired
folder is on any job's PATH. The launcher chain (jobs 54916369, 54916370) passed the guard and
found nothing to do on the two up-to-date studies. My chain-job output files were removed from
the clone roots (my ids only).

**Deletion.** On or after 2026-10-01, if both sessions have run a study from the new folders by
then, delete `curvature-v2512.retired-2026-09-24`, `sdpls-v2512.retired-2026-09-24`,
`~/.leia_env.retired-2026-09-24`, the two `retired-leia-2026-09-24` subfolders (cluster and
laptop); record the date here. Until then they are the rollback of this work package.

### 10.6 OpenFOAM-v2606 on the laptop; the davof workspace; the TwoPhaseFlow hook (2026-09-28)

**Why.** The DAVOF (dual area/volume-of-fluid) work will be developed in leia and benchmarked
against TwoPhaseFlow's plicRDF/isoAdvector through the `interFlow` reference arm. The newest
OpenFOAM on the laptop is v2606 (`$HOME/OpenFOAM/OpenFOAM-v2606`, Opt only), and TwoPhaseFlow's
upstream master was ported to it on 2026-08-19. Decision: v2606 is the standard from today,
everywhere in the repository (the cluster files too; Lichtenberg has only v2512 and needs a
v2606 source build before its studies run again).

**What changed (three commits before this one).**
- `etc/leia-env.sh`: the strip rule no longer hard-codes `-v2512`. The account default that
  `etc/bashrc` put first on the paths is captured before `WM_PROJECT_USER_DIR` is overridden and
  removed by its exact prefix, for whatever version is sourced (fallback: the bashrc's own
  formula); the `<name>-vNNNN/platforms` pattern stays for retired user dirs. New: a version
  stamp, `platforms/<WM_OPTIONS>/.openfoam-version`, written by `Allwmake`, removed by
  `Allwclean` (before sourcing), refused by `leia-env.sh` on mismatch, because v2512 and v2606
  share `WM_OPTIONS=linux64GccDPInt32Opt` and the loader cannot tell the two apart. cfMesh:
  `$HOME/OpenFOAM/cfmesh-$WM_PROJECT_VERSION` when it exists, else `$HOME/OpenFOAM/cfmesh`, with
  a warning if a stamped cfMesh belongs to another version. TwoPhaseFlow: `$LEIA_TPF_DIR`
  (default the `TwoPhaseFlow` clone next to this one) appended after leia's own directories, the
  same shape as the cfMesh hook, so the `interFlow` arm needs no study-config reference.
- One `sed` over 363 files: `OpenFOAM-v2512/etc/bashrc` -> `OpenFOAM-v2606/etc/bashrc` in 349
  study configs, `config/config.yaml`, the slurm and local8 profiles, both sbatch drivers, and
  the workflow scripts (`local-mesh-wrappers/pMesh` included).
- Docs: README, AGENTS/CLAUDE (byte-identical), requirements, Makefile, CLUSTER, SLURM, METHOD.

**The davof workspace, `$HOME/OpenFOAM/repos/davof/`.** Fresh clones, both building into their
own `platforms/`: `leia` (this branch) and `TwoPhaseFlow` on a local branch `davof` = the fork's
`feature/parasitic-currents` (the `Foam::twoPhaseFlow` inline-namespace fix for the
libgeometricVoF collision, 5aa7e9a) merged with upstream `ba10702`; no conflicts. TwoPhaseFlow got
an `etc/davof-env.sh` of the same shape as ours (account default stripped, clone-local install,
version stamp), sourced by its `Allwmake`/`Allwclean`; no `Make/files` or `Make/options` changed.
Not pushed. An untracked `davof/env.sh` sources OpenFOAM, then TwoPhaseFlow, then leia.
`$HOME/OpenFOAM/cfmesh-v2606` was built from OpenFOAM-v2606's own `plugins/cfmesh` with
`./Allwmake -prefix=<dir>/platforms/<WM_OPTIONS>` (25 binaries, `libmeshLibrary.so`, stamped);
the SourceForge-based `$HOME/OpenFOAM/cfmesh` stays the v2512 build.

**MEASURED, laptop, 2026-09-28.**
- leia on v2606: `./Allwmake` clean without any source change; `leia-check-deps: PASS`, eight
  stamps written, 22 binaries and 11 libraries in `platforms/linux64GccDPInt32Opt`, stamp `v2606`.
- Bit-identity across OpenFOAM versions: `cases/2Dtranslation` (367 steps) run with the v2512
  build of `repos/leia` (c431677) and with this v2606 build; `compare_metrics_csv.py --tol 0`
  on `leiaSemiLagrangeLevelSetFoam.csv`: 368 rows, 12 columns, max relative difference 0 in
  every column, PASS. Same result at 1e-12 and 1e-8, as it must be.
- The environment: after `etc/leia-env.sh` under v2606, zero `tmaric-v2606` entries on PATH and
  LD_LIBRARY_PATH (that directory holds a stale `libVoF.so`, `libalphaFieldFunctions.so` and
  `setAlphaField` of an unrelated build), eight OpenFOAM-v2606 entries kept, leia's directories
  first, then OpenFOAM and ThirdParty, then `cfmesh-v2606`, then the TwoPhaseFlow clone. Sourcing
  twice adds nothing. Under a v2512 shell both env files refuse with the stamp message (rc 1).
- TwoPhaseFlow on v2606: builds clean (7 executables, 9 libraries). Collision checks on
  `libVoF.so`: 41 defined symbols under `Foam::twoPhaseFlow::reconstructionSchemes`, 0 under
  `Foam::reconstructionSchemes`, 0 strong symbols shared with `libgeometricVoF.so`;
  `run/parasiticCurrents/stationaryDroplet2D` (function objects loaded) exits 0 after 20 steps,
  0 "Duplicate entry"; `interFlow -lib geometricVoF` on the same case also exits 0 with 0
  duplicates. `ldd interFlow` resolves `libVoF`/`libsurfaceForces` in the clone and
  `libalphaFieldFunctions` in OpenFOAM-v2606, nothing in `tmaric-v2606`.
- The interFlow arm through the real launch path: a temporary study (copy of
  `config/interFlowDroplet2D.yaml`; N 64, RDF, np 4, END_TIME 0.002) via
  `make studies-one`: 10 of 10 Snakemake steps, `interFlow -parallel` 94 steps to `End`, 0
  "Duplicate entry", `interFlow.csv` written by the post-solve script. The only `interFlow` on
  a v2606 path is the davof clone's (`command -v` in the driver shell). The report step's
  `make_convergence_table.py` prints "no sdpls studies found" and rc 1 for this arm, as before
  (the report rule tolerates it). Study, config and the copied errors table were deleted.
- Account defaults untouched: `ls -la --time-style=full-iso` of `tmaric-v2606` and
  `tmaric-v2512` `platforms/*/{bin,lib}` before and after are identical, and `find -newer` the
  before-listing finds 0 files (`davof/logs/acct-default-*-{before,after}.txt`).

**Open.** Lichtenberg: build OpenFOAM-v2606 there before the next cluster study (every preamble
now sources v2606). The cfMesh codebases differ (community plugin on v2606, SourceForge on
v2512), so polyhedral meshes may differ slightly between the two; measure before comparing
poly results across the versions. The `feature/eulerian-rholent` branch edits `Allwmake` and
`etc/leia-env.sh` too; its merge onto this `leia-env.sh` is manual.

### 10.7 DAVOF static gate 1: the Gauss-identity normal on a sphere, measured (2026-09-28)

**Why.** The DFG proposal's reviewers asked for a preliminary result: initialise the cell
volume fractions `alpha_c` and the face liquid-area fractions `alpha_f` of a sphere, recover
the interface area normal of every cell from the DAVOF identity `m_c = -sum_f alpha_f S_f^out`
(`|m_c|` the interface area in the cell, the direction the normal out of the liquid), and show
its convergence next to the best PLIC reconstructions (TwoPhaseFlow's plicRDF, gradAlpha,
isoAlpha). Expectation before the run: DAVOF second order, the PLIC schemes first.

**What exists now (commits 98bb59f, 9cfcd4b, this one).** `src/leiaLevelSet/davof`
(`libleiaDavof`, `davofState.H`): the DAVOF state with three sources of `(alpha_c, alpha_f)`:
`exactSphere` (closed-form circle-polygon face fractions, quadrature cell fractions, both
exact), `linearInterpolant` (the exact signed distance at the cell centre, the face centres
and the points, linear on the Detrixhe-Aslam tets and on the face fans with the same vertex
values: ONE C0 piecewise-linear surface, for which the identity is exact), and
`planePhaseIndicator` (leia's production plane-DA state). `applications/test/leiaTestDavofNormal`
scores `e_c = |m_c/|m_c| - n_exact(xS_c)|` at the reconstructed patch centroid over the
interface cells (`|m_c| > wispTol V_c^(2/3)`, wispTol 1e-3) in L1 (mean), L2 (rms), Linf,
area-weighted variants, the total area `sum_c |m_c|` against `4 pi R^2`, the tet/face
consistency and the closure, and runs OpenFOAM's geometricVoF plicRDF / gradAlpha / isoAlpha
(the algorithms TwoPhaseFlow uses) on the IDENTICAL alpha as the like-for-like cross-check.
Study `config/davof/sphereNormal3D.yaml` (R = 1 mm in a 4R box, centre off every symmetry
plane and grid line, hex N = 20/40/80/160, R/h = 5..40, serial, 12 cases in ~12 min),
`make studies-one STUDY=davof/sphereNormal3D`; the plane gate `cases/davof/planeNormal3D/Allrun`;
results in `docs/davof/davof-article/data/{tables,figures}` (the tpf rows merged with
`make_davof_normal_table.py --tpf-results`). The TwoPhaseFlow arm (davof clone, branch `davof`,
commits 5fdca12 and ac879d1, not pushed): `reconstructionError` gained six appended columns
(`LNormalDiff1/2/Inf`, `areaSum`, `areaExact`, `nInterface`; the first twelve byte-identical)
and `run/benchmark/reconstruction/sphereNormal3D` is the same ladder through its own
framework (genCases/runAll/getResults + pytest/oftest in `~/OpenFOAM/repos/davof/tpf-venv`).

**Pre-registered prediction 1 (config header, before run 1): FALSIFIED for the normal, held
for the area and the identity.** With the `linearInterpolant` state the DAVOF normal converges
at order ONE, not two: p(L2) = 1.01 (pairwise 1.05 / 1.01 / 0.99), p(L1) = 1.02, L2 at N = 160
2.43e-3; the area error is second order, p(A) = 2.00 (1.13e-4 at N = 160); the identity is exact,
MAX_CONSISTENCY = |m_c - sum of the tet patches|/h^2 <= 4.2e-14, closure |sum_c m_c|/A ~ 1e-15.
Mechanism: the face crossings of a linear interpolant of the exact distance are O(h^2)
accurate, so the face fractions carry O(h) errors as fractions, and the direction of a sum of
O(h^2) vectors weighted by O(h)-accurate fractions is O(h) accurate; the magnitude (area) is
not affected at leading order. The same holds for the plane-DA state (p(L2) = 1.13, p(A) =
2.00, L2 at N = 160 4.50e-3, Linf 0.058, larger by the C0 mismatch of neighbouring planes).

**Pre-registered prediction 2 (config header, before run 2): HELD.** With the exact face and
cell fractions (`exactSphere`) the DAVOF normal is second order: p(L2) = 1.98 (pairwise
2.16 / 2.03 / 1.73), p(L1) = 2.05 (2.14 / 2.05 / 1.96), area-weighted p(L1) = 2.01, p(L2) = 1.98,
p(A) = 2.00; L2 falls 6.96e-3 -> 1.56e-3 -> 3.83e-4 -> 1.15e-4 over N = 20/40/80/160, the area
error to 5.1e-5. Linf (p = 1.22, 5.5e-3 at N = 160) is set by the smallest corner clips just
above the wisp threshold and is reported, not gated. Interface cells 472 / 1863 / 7418 / 29739,
wisps 8 / 34 / 122 / 432. The DAVOF normal costs 0.02..3.0 s serial per rung.

**The cross-check on the identical alpha (OpenFOAM geometricVoF in the leia app).** plicRDF:
p(L2) = 1.06 on every state (exact 1.08 / 1.11 / 0.99 pairwise), L2 at N = 160 3.87e-3 -- on the
exact alpha 34x the DAVOF error, on the linear-interpolant alpha 1.6x. gradAlpha (the
centred-gradient normal): p(L2) = 0.07, L2 = 0.066 at every rung, 0.124 even on the exact
plane of the plane gate: the discrete gradient of a sharp indicator does not converge to the
normal (Pilliod & Puckett, JCP 199, 2004). isoAlpha: p(L2) = 0.52, 0.021 at N = 160. The
TwoPhaseFlow arm agrees to within a few percent on its own alpha and error measure: plicRDF
cutCellIso p(L2) = 1.08 (1.09 / 1.14 / 1.00), 3.97e-3 at N = 160, cutCellImpFunc 1.06, 4.01e-3;
gradAlpha 0.07, 0.066; isoAlpha 0.55 / 0.51, 0.0216 / 0.0214; pytest (N <= 80): 7 passed, the
order assertion is 0.9 for plicRDF (measured 1.07) and none for gradAlpha (0.11; the
pre-registered 0.9 for it was wrong and the test now says so). The "plicRDF is second order"
expectation from its paper is NOT reproduced in this metric on this ladder by either code.

**Gates and inertness.** Plane gate at tolerance 1e-12 for BOTH leia states (they are exact
for a plane): linearInterpolant normal 5.3e-14, consistency 2.1e-15, alpha vs the exact plane
cut 2.1e-15; planePhaseIndicator 9.5e-13, 3.4e-14, 2.7e-14 (the least-squares fit's
conditioning); plicRDF on that exact plane alpha 1.9e-4, gradAlpha/isoAlpha ~0.1. Parallel:
np 4 vs serial on the N = 40 sphere, all three sources, every norm within 5.5e-11 relative,
round-off quantities within 7.3e-13 absolute. The Snakefile edit (a `case:` may name a
sub-folder; per-case directories use the basename): the dry-run job listings of all 352
configs are identical before and after (sorted; Snakemake's job order is not deterministic);
the report rule's behaviour for every existing solver is unchanged (only the new app is
exempted from the generic plots). The plane-DA alpha of `leiaSetFields` is reproduced bit
for bit (MAX_ALPHA_DIFF_DA = 0): its narrow band is built in createFields.H on a uniform psi,
so it initialises EVERY cell from a plane, and cells 0.5..0.9 h away still clip a corner --
a fact worth knowing about the production initialiser.

**What it means for the proposal.** The identity is exact and its normal is as accurate as
the face fractions are as fractions: O(h^2) face fractions give a second-order normal that
beats plicRDF by more than an order of magnitude at equal input; O(h) face fractions (any
linearised initialisation, and any transport that only locates the interface to O(h^2)) give
a first-order normal, still 1.6x better than plicRDF on the same input. WP1's face update
therefore has to deliver `alpha_f` to O(h^2) as a fraction, i.e. the face crossing to third
order or the face liquid area by an exact-geometry evaluation, and that is the next gate: the
one-step face-fraction update the reviewers asked for. The total area is second order for
every state.

**Incidents.** Two `pkill -f` patterns matched their own `wsl.exe ... bash -c` command line
and killed it (exit 15) before the kill lines ran; the TwoPhaseFlow ladder was then started
four times over and its 4.1e6-cell cutCellImpFunc cases (~8 GB each) exited with rc 1 out of
memory. Stopped from a script whose invocation matches no pattern (`scratch/stop-tpf.sh`),
restarted once; `runAll.sh` now runs coarse to fine (an alphabetical glob starts with N160).
`pgrep -x` cannot match `reconstructInterface` (names > 15 characters). PyFoam needs a
FoamFile header in `simulationParameter` (added; OpenFOAM merges it on `#include`).

**Open.** Push `development` after review. The polyhedral (cfMesh) ladder of the same gate
(`mesh: poly`), the one-step face-fraction update test, and a wispTol sweep are one config
each away. `leiaSetFields`' every-cell plane initialisation deserves its own note in the
phase-indicator docs.

### 10.8 DAVOF static gate 1b: the explicit plane position, measured (2026-09-29)

**Why.** The 2020 area-of-fluid draft (`research/articles/2020-article-area-of-fluid`) holds the
explicit position formula the proposal alludes to. Re-derived from the Gauss theorem for x on the
liquid sub-volume, corrected in the draft (sign; one degree of freedom, the signed distance p_k of
the plane from the cell centroid, not a "centroid" vector; the face-centroid form is EXACT on
planar faces because (x - x_k).n_f is constant there; a remark on the d-fold amplification of a
volume defect) and added to the proposal as eq:avof-position; verified symbolically and in exact
rational arithmetic (`verification/aof_identities_sympy.py`: corner cuts of the cube and the
square, 157 rational planes through a cube in every configuration, zero residual, the 2020 sign
returns -p_k, a defect dV moves the plane by 3 dV/|A|). Implemented as
`davofState::planePosition()`:

    p_c = ( 3 alpha_c V_c - sum_f alpha_f (x_f - x_c).S_f^out ) / |m_c|,   plane {x : (x - x_c).n_c = p_c},

then the plane is cut through the cell on the Detrixhe-Aslam tets (exact for a plane): the
polygon centroid `xPlane.davof`, area `APlane.davof`, cut fraction `alphaPlane.davof`. Position
error per model: E_POS_L1/L2/LINF, the distance [m] of the polygon centroid (the geometricVoF
models' `centre_`) from the exact surface over the interface cells; realizability diagnostics
MAX_VOL_DIFF_PLANE = max |alphaPlane - alpha|, MAX_AREA_DIFF_PLANE = max |A n - m_c|/h^2; the
plane gate checks E_POS_LINF/h and both diagnostics at 1e-12 for both leia states.

**Pre-registered prediction 3 (config header): HELD.** p(E_POS_L2): exactSphere 2.00 (pairwise
2.00 / 2.00 / 2.00), linearInterpolant 2.03 (2.07 / 2.03 / 1.99), planePhaseIndicator 2.00;
E_POS_L2 at N = 160: 4.35e-8 m (1.7e-3 h), 1.69e-8 m, 2.02e-7 m; L2/h falls linearly (1.4e-2 ->
1.7e-3 for the exact state). Every plane-based position is second order, plicRDF, isoAlpha and
gradAlpha included (2.00 / 1.99 / 1.95 on the exact alpha; TwoPhaseFlow's own LCentre1: plicRDF
2.00). On the exact alpha the explicit position coincides with plicRDF's volume-matched one to
0.02 % (4.350e-8 vs 4.349e-8 m): the identity lands on the volume-matching plane when the state is
consistent to leading order. On the linear-interpolant alpha the explicit position is 4.7x more
accurate than the volume-matched plane of the same alpha (1.69e-8 vs 7.94e-8 m); on the
plane-indicator state, where the two planes of a face disagree, the volume-matched plane is
slightly better (1.74e-7 vs 2.02e-7 m). DAVOF's foot point x_c + p_c n_c: L2 4.6e-8 m at N = 160
for the exact state, the same order.

**One diagnostic prediction was WRONG.** MAX_VOL_DIFF_PLANE for the curved exact state was
pre-registered as O(h^2); measured 5.9e-2, 3.3e-2, 1.7e-2, 8.6e-3, i.e. O(h): the fraction cut by
a plane at the patch's mean offset differs from the exact fraction by (sagitta x patch area)/V_c =
O(h^2/R x h^2 / h^3) = O(h/R). The area diagnostic is O(h) as predicted, noisy in the max
(0.18, 0.099, 0.087, 0.056). Both are round-off on the plane gate (1.9e-15 / 1.5e-14 and 6.0e-14 /
9.6e-14 for the two states); E_POS_LINF there 1.0e-16 m and 8.9e-17 m.

**Where things are.** leia commit (this one): davofState.{H,C} (planePosition, four fields),
leiaTestDavofNormal (position norms, foot point, realizability, gate), make_davof_normal_table.py
(E_POS_* metrics, figure, tex column; TwoPhaseFlow's LCentre1/LCentreInf merged as E_POS_L1/LINF),
config header, README/docs, theme data (`davof_normal_e_pos_l2_sphereNormal3D.png` and the
tables). Article: `sections/area-of-fluid-method-interface-recon.tex` corrected in place (comments
mark every change), `sections/results.tex` written with the plane and sphere tables,
`verification/aof_identities_sympy.py` added. Proposal: eq:avof-position with the two properties
(fraction accuracy decides the order; 3-fold sensitivity to inconsistency) inserted after
eq:avof-normal. Not touched in the proposal, still flagged from the 2026-09-28 review: the face
update's denominator (|S_f| where the pre-image area belongs) and the curvature factor (H vs 2H).

**Next (agreed 2026-09-29).** Proposal: the order tables and TikZ schematics of both tests in the
preliminary-work section. leia: TwoPhaseFlow's PLIC VTK writer as a leia library, the plane and
sphere interfaces of gradAlpha, plicRDF and DAVOF at three resolutions rendered to PDF by a
Snakemake/Python job and included in the proposal automatically.

### 10.9 The PLIC surfaces as VTK, the proposal's tables and figures (2026-09-29)

**What.** `libleiaDavofInterface` (`src/leiaLevelSet/davofInterface`, ported from TwoPhaseFlow's
`src/postProcessing/interface`, Henning Scheufler, DLR): `plicInterfaceSurface(mesh, normalOut,
pointOnPlane)` cuts the plane of every interface cell through the cell with OpenFOAM's
`cutCellPLIC` into one polygon of a `meshedSurface` (`meshCells()` maps polygons to cells) and
writes ASCII legacy VTK polydata with per-polygon CELL_DATA from cell fields; the
`sampledSurface` `type davofInterface` (fields `normal m.davof`, `centre xPlane.davof`, optional
`wispTol`) is the function-object form for the `surfaces` FO of a solver. Decoupled from the
reconstructionSchemes registry on purpose: the caller names the fields, so the DAVOF plane and
a geometricVoF scheme (-normal_, centre_) are drawn by the same code. `leiaTestDavofNormal`
writes `postProcessing/davofInterface/<time>/plic.{davof,plicRDF,gradAlpha,isoAlpha}.vtk` with
`cellId`, `eNormal`, `ePos`. `workflow/scripts/render_davof_interfaces.py` (report rule, after
the table) parses the VTK itself and draws models x resolutions (the coarsest three) coloured by
the error on one log scale, orthographic projection with a painter's sort (mplot3d is broken on
this host: apt + pip matplotlib), head-on for a planar surface; `figures/davof_plic_<study>.{pdf,png}`.
`config/davof/planeNormal3DLadder.yaml` (the gate's plane at N = 8/16/32, seconds) gives the
plane renderings. `make_davof_normal_table.py` writes the compact
`tables/davof_normal_proposal_<study>.tex`; `export_davof_proposal.py --dest`
(`make davof-proposal DEST=...`) copies table and figures into the proposal, which `\input`s
and `\includegraphics` them (`figures/davof/`), so a study rerun plus the export updates the
proposal's preliminary-results section without editing it.

**Proposal (2026-DFG-AVOF, uncommitted there).** A paragraph "Preliminary results for the DAVOF
core (static gates, 2026)" before "Technical positioning": the two gates, the compact order table
(Table davof-prelim), TikZ schematics of both gates (the cut cell with alpha_f, the area normal
and p_k; the sphere in its box with the interface cells and a zoomed cell with the exact arc,
the recovered plane, e_k and d_k), the three findings, and the PLIC figure (plane and sphere,
DAVOF / plicRDF / gradAlpha x three resolutions). Builds clean with xelatex (Windows TeX Live
2025; WSL lacks Arial), 0 unresolved references.

**MEASURED, the plane ladder (N = 8/16/32, linear interpolant).** DAVOF normal L2 = 1e-14 at every
N (round-off; the fitted "order" is meaningless), plicRDF 1.3e-4, gradAlpha 8.8e-2, isoAlpha
7.7e-2 at N = 32: the comparators do not recover a plane from exact plane fractions. The sphere
study was rerun (run 4, identical metrics, now with the VTK surfaces).

**Incidents.** A sed replacement wrote a literal `\n` into the Allwmake loop (fixed by a
script); the report rule tolerates a failing table/render step (rc recorded in the log, the
study passes) -- watch for "FAILED (rc=1)" lines; the plane ladder's first run recorded no
DAVOF tokens because its fvSolution was not a template (case_params.json records only the
tokens the templates reference) -- fvSolution.template now, DOMAIN_LENGTH in the .parameter.

### 10.10 The order is the fraction order: quadratic faces, the arm renamed, the proposal's equations (2026-09-29)

**Why.** Tomislav's reading of the first runs: the initialisation should be the exact signed
distance at the cell centre, the face centres and the corner points in the Detrixhe-Aslam
formulas, and that should give a second-order normal without iteration. That IS the arm the study
had called `linearInterpolant` (a misleading name: nothing is interpolated from elsewhere; "linear"
is what the DA formulas assume inside each simplex), and it measured first order. Renamed
`detrixheAslam` (alias kept). The mechanism, now written down: the DA formulas integrate the chord
surface, whose crossings are O(h^2) off a curved interface; per face that is an O(h) error AS A
FRACTION although the integrated volume and area are second order; the Gauss sum's direction
inherits the per-face order. Hence two estimates that define WP1 (now equations in the proposal):
|dn| <= sum_f |d alpha_f| |S_f| / |A| = O(h^q) for |d alpha_f| = O(h^q), and |d alpha_f| = O(d_Sigma/h)
with d_Sigma the position error of the discrete interface used to clip the face: any per-cell
plane gives q = 1, O(h^3) gives q = 2.

**New arm `quadraticFaces`** (`davofQuadraticFaceGeometry.H`): the exact signed distance also at
the edge midpoints of every fan triangle, the quadratic Lagrange interpolant cut exactly (edge
roots in the numerically stable form; the sliver between chord and conic arc by 12-point
Gauss-Legendre, analytic integrand); cell fractions stay DA. Triangles with other than two
boundary crossings or an arc that is not a graph over its chord keep the DA value and are counted
(N_FACE_FALLBACK).

**MEASURED, pre-registered prediction 4: HELD.** quadraticFaces p(E_L2_N) = 1.95, p(E_L1_N) = 2.05,
E_L2_N within 4 % of exactSphere at every rung (1.19e-4 vs 1.15e-4 at N = 160), p(A) = 2.00,
p(E_POS_L2) = 2.01; fallbacks 4 / 2 / 0 / 4 faces on the four meshes (the "none at R/h >= 10"
clause missed by a handful). Plane gate exact for all three leia states (E_POS_LINF 9e-17 m). A
first version of the edge-root formula failed the plane gate at 1e-6: for a plane the quadratic
degenerates to a line and the textbook root formula cancels; the stable form fixed it. So: a
second-order normal from exact signed distances alone, non-iteratively, for any implicit surface,
with one polynomial degree more in the face cut. The transport must preserve exactly that
fraction accuracy, which is the WP1 question as now stated.

**Proposal (2026-DFG-AVOF, committed there).** After the face update: the two estimates as
displayed equations (eq:normal-propagation, eq:clipping-law); the Consistency Hypothesis restated
on the normal and position orders with q, the realizability defect a diagnostic relative to the
exact state's O(h) floor; Stage 1 names the two q = 2 candidates (curvature-carrying clipping
surface, sliver correction) and the self-consistency argument (a first-order curvature suffices,
second-order normals provide it), stability of the loop as Stage 4; success criteria split
(position >= 1.8 unconditionally, normal >= 1.8 with q = 2, the q = 1 outcome a stated result);
"compatible with, but not dependent on" replaced by the designed coupling with the iso-surface
fallback; "no stencil" replaced by "cell-local recovery, one ring for the curvature in the
update"; Table 1 caption names the cell set and the wisp threshold; Table 1 is nine rows (exact,
quadratic faces, Detrixhe-Aslam x DAVOF, plicRDF, gradAlpha), isoAlpha and the plane indicator
quoted in the text. Page budget kept at 17 by dropping the WP2 convergence plot (its numbers are
Table 2's), the repeated sentences in WP2, WP3 and the timeline, and shorter captions.

**Open (agreed).** The remaining question Tomislav wants to discuss next: the q = 2 update itself
(which candidate, and whether to show a one-step result before submission). Then the 25 September
items outside WP1.

### 10.11 DAVOF static gate 2: the curvature of the DAVOF state on an ellipsoid, measured (2026-09-29)

**What.** The curvature of the DAVOF state as run-time-selectable models inside
`libleiaDavof` (`src/leiaLevelSet/davof/curvature`, `fvSolution davof.curvature`,
nothing in the level-set libraries): `quadricFit` fits, per interface cell, a
quadric graph in the frame of `n_c = m_c/|m_c|` about the plane-polygon centroid
to the Hermite data of the point-neighbour ring (OpenFOAM's `CPCCellToCellStencil`,
gathered for the interface cells only, parallel-aware), two slope rows per
neighbour from its normal and one position row from its centroid weighted by
`positionWeight` (0 = the normals decide the second derivatives, the headline
model `normalsOnly`; 1 = the naive Hermite fit, `hermite`). Convention
`kappa = kappa_1 + kappa_2 = div(n)`; also `K`, `kappa1`, `kappa2`. Delivery
(`davofCurvatureDelivery.H`) with the parallel-surface closed form of the
level-set solver, because the interface centroid is neither the cell centre nor
a face centre: the contour-referenced cell field on the force band by the forward
map `kappa(d) = (kappa + 2dK)/(1 + d kappa + d^2 K)` at the cell centre's offset
(the stable quadratic root of the normal-projected SL note), the face field by
`interpolate(kappaCell)` and `parallelSurfaceInverse(kappa_f, d_f, K_f)` (copied
verbatim from `stabilizedFootPointFaceCurvature.H`), next to the models' own foot
value; packed models cross coupled patches by `syncTools::swapBoundaryFaceList`.
`leiaTestDavofNormal` scores every model against the exact geometry
(`signedDistanceEllipsoid`: total curvature at the closest point, the Gaussian
curvature from the axes; `implicitSphere` 2/R; `implicitPlane` 0), writes
`leiaTestDavofCurvature.csv` and the headline columns into the wide and the
models CSV; `-expectExact` gates `E_KAPPA_LINF h` at round-off. New study
`config/davof/ellipsoidNormal3D.yaml` (half-axes 1.0/0.8/0.6 mm, the sphere
gate's centre, N = 20..160, states quadraticFaces / detrixheAslam /
planePhaseIndicator, the geometricVoF cross-check), PREDICTION 5 in its header.
`make_davof_normal_table.py` carries the curvature metrics, the series
`davof:curv:<model>` and the compact proposal table with `p(kappa)` and the
relative error (the p(A) column dropped).

**Measured (the config header has the full record).**

| state | p(n) L2 / L1 | p(x) | p(kappa) L2, pairwise | rel. L2(kappa) at N = 160 | p(K) / cell / face |
|---|---|---|---|---|---|
| quadraticFaces | 1.78 / 1.96 | 2.00 | 1.55 (1.91, 1.54, 1.21) | 4.55e-3 | 1.59 / 1.36 / 1.63 |
| detrixheAslam | 1.00 / 1.05 | 2.04 | 1.31 (1.88, 1.34, 0.70) | 6.75e-3 | 1.23 / 1.25 / 1.38 |
| planePhaseIndicator | 1.12 / 1.14 | 2.01 | 1.28 (1.93, 1.36, 0.53) | 8.60e-3 | 1.28 / 1.18 / 1.35 |

plicRDF on the identical alpha: 1.10 (5.36e-3 at N = 160, 18x the quadratic-face
DAVOF error 3.01e-4), position 2.00; gradAlpha 0.10; isoAlpha 0.65. Area order
2.00 for every state. `hermite` within 3 % of `normalsOnly` everywhere.

**Against PREDICTION 5.** HELD: quadraticFaces normal second order in L1 and
position second order; curvature at least first order with the relative error
inside the predicted 2e-3..8e-3 at N = 160; K, cell and face deliveries at the
same orders; plane gate at round-off; np-4 = serial. MISSED by 0.02: the L2
order of the quadratic-face normal (1.78 against the 1.8 threshold) -- 47 of
18642 interface cells (0.25 %) carry e > 2e-3 and set the rms (1.72e-4 without
them); the L1 order is 1.96. ABOVE the predicted band: the curvature order
(1.55 against 1.0 +- 0.3, pre-asymptotic, the finest pair 1.21). FALSIFIED (two):
(a) the curvature of the q = 1 states converges (1.31 and 1.28, criterion
"> 0.5"), 1.5x less accurately than q = 2 at N = 160 and slowing at the finest
pair (0.70 against 1.21) -- the chord state's O(h) normal error is smooth along
the interface, so its surface gradient is O(h)/R, not O(1); the curvature
separates the fraction orders LESS sharply than the normal does (1.00 against
1.78), the opposite of the pre-registered claim; (b) the naive Hermite fit shows
no point-fit floor: in the joint least squares the two slope rows per neighbour
(leverage U) outweigh the one position row (leverage U^2/2), so the centroids
barely move the second derivatives; a floor would need positions alone.

**Operational.** With the geometricVoF cross-check a case needs 1.2 GB at
N = 80 and about 9 GB at N = 160 (0.62 GB at N = 80 without it): the ellipsoid
study runs ONE case at a time (`scratch/run-study-ellipsoid.sh`, `--jobs 1`);
the first run with three concurrent N = 160 cases was OOM-killed, and the ring
data were then still gathered for all 4.1e6 cells (6 GB) -- now the stencil
lists of non-interface cells are cleared before the map is built. WSL shuts the
distro down when no session is open, which killed a detached (`nohup &`) run
mid-ladder: long runs go through an attached session. The Snakefile edit is one
line (the curvature CSV in the `_preprocess` rm list; `git diff`), the sphere
study re-renders its tables unchanged apart from the new column layout.
Commits: the library and app, the study and docs, this record.

## 11. Halo-limited extension and gradient-control sources (2026-09-26 ->)

Plan: `docs/plan-halo-limited-gradient-control.md` (approved 2026-09-26). Branch
`feature/gradient-controlled-level-set`, laptop worktree `~/OpenFOAM/repos/leia-gcls`; the
never-edited reference worktree `~/OpenFOAM/repos/leia-gcls-base` holds the pre-change binaries
for every bit-identity gate: before each gated step it is checked out (detached) at the last commit
before that step and rebuilt (c431677 for B3; 0c7079c for Phase C). Theme: `docs/gradient-controlled-level-set/`.
Decisions of 2026-09-26: the SL line first; 3D h ratio >= 1.3; the dossiers stay out of git; ONE
2D and ONE 3D gate study test every method (no per-method study configs).

### 11.1 Phase A: rules, plan, docs theme (2026-09-26)

- CLAUDE.md and AGENTS.md: the section "Extension without modification" (byte-identical).
- `docs/plan-combined-source-terms.md` (WP2, WP3) and `docs/combined-source-terms/improvement-sdpls-combined.md`
  are SUBSUMED: the combined source is `gradientControl` with law `linearQ` and strain weight `full`.
- The theme skeleton builds (`make article-gcls`, the deck exports).

### 11.2 Phase B1: workflow repairs, and a retraction (2026-09-26)

**RETRACTED: every 2D advection order in METHOD.md section 8.3.7 was 3/2 of the true value.**
`advection_convergence_table.py` used h_eff = nCells^(-1/3) for every case, so a 2D case with
N^2 cells had h_eff = N^(-2/3). The corrected orders (h = 1/N) are in METHOD.md 8.3.7: the
translation `none` arm goes 3.80, 2.10, -0.54 (was 5.70, 3.15, -0.82), the vortex `none` arm
2.76, 3.23, 2.09 (was 4.15, 4.85, 3.14). The conclusions of 8.3.7 do not change. OPEN: regenerate
the curated `advConv2D*_convergence.csv` on Lichtenberg (the studies are only there).

Other repairs, one commit:
- `aggregate.py`: the token copy sat inside the `leia.version` branch since WP2 (cc79df4), so a
  case without a stamp lost N_CELLS, END_TIME and h, and had no completeness check. MEASURED on
  three existing studies, old against new aggregator: every changed column goes from blank to its
  value, except one diverged run (`alphaFTest_averagedPlanes`, DIVERGED at step 4634), whose
  mid-run row is no longer published as final. The method label of eight `sdpls1Dstretch` arms
  changes from `euler+div:...` to `euler+SDPLS:R/...`: the old label hid which arms had a source.
  No physical value changed.
- `aggregate.py`: completeness now follows OpenFOAM's stopping rule (t_last >= END_TIME - dt/2)
  instead of 0.99 END_TIME; `seamConsistency3Dpar4` (t_last = 98.7 % of END_TIME) is complete.
- The solve rule writes `<case>/.leia_launch` (np, SLURM job id, host, start, end, exit code,
  classifier line); `aggregate.py` reports `nRanks` (from the log header), `wallClockSolve`,
  `steps`, `logState`, `slurmJobId`. `foam_log_state.sh` prints an additive `nprocs=` field.
- DEFERRED, D-c (dims of the `*Droplet2D` cases recorded as 3): the only readers are the 2D field
  plots of `plots.py` and `marker_ref.py`; a fix would switch them on for every droplet study, and
  the gates take the dimension from the gate definition.

### 11.3 Phase B3: the droplet band gradient uses the unlimited metric scheme (2026-09-26)

`writeDropletMetrics.H` computed `minGradPsiBand`, `maxGradPsiBand`, `gradPsiL2ErrorBand` and
`gradPsiRatioBand` with `fvc::grad(psi)`, the advection scheme `cellLimited leastSquares 1`; it now
uses `gradPsiMetric` (unlimited `leastSquares`), as the templates require. Gate (laptop, np 4,
`stationaryDroplet2D` N = 64, 30 steps, reference binaries of c431677 against the new binary): all
four CSVs identical at tolerance 0, the four columns included. The limiter acts only at a local
extremum of psi, and a smooth band has none; the fix changes the columns only where the band
profile is corrupted.

### 11.4 Phase B4: the oscillating droplet's level-set surface is a token (2026-09-26)

`cases/oscillatingDroplet2D` initialised psi with `implicitEllipsoid`, which is ALGEBRAIC:
psi = sum (x_i - c_i)^2/a_i^2 - 1, so |grad psi| = 2/a_i is 1.8e3 to 2.2e3 at the interface. Its
gradient columns, and every band criterion in psi units, had no meaning in any oscillating-droplet
study. New token `DROPLET_SURFACE` (default `implicitEllipsoid`, bit-identical); the method gates
pin `signedDistanceEllipse`. Gate: `oscillatingLadder2Dshared` rendered before and after; only three
comment lines differ. OPEN, for the author: void the past `oscillatingDroplet2D` studies
(`_VOID_algebraicPsi_<date>`) or keep them with this caveat.

### 11.5 Phase B5 and B6: the method gates exist, and the 3D oscillating droplet (2026-09-26)

- `config/gates/methodGate{2D,3D}.yaml`, `config/candidates/baseline{,Eulerian}.yaml`,
  `workflow/scripts/render_gate_configs.py`, `workflow/Snakefile.gate`,
  `workflow/scripts/make_gate_summary.py`, `workflow/scripts/richardson.py`, Makefile targets
  `gate` and `studies-one-file`; CLAUDE.md and AGENTS.md section "Method gates".
- `richardson.py --self-test`: 64 checks PASS (orders p = 1, 2, 3 recovered to 1e-8 on the four
  integer gate ladders, extrapolated value to 1e-10, the three convergence types).
- 2D smoke, laptop, baseline (`SMOKE=1`, N 32/45/64, about 20 steps): all 14 cases COMPLETED at
  np 4; exact1D q error 1.6e-6 (PASS); seam np 1 and np 8 against np 4: column-scaled
  difference <= 1.9e-12 (PASS).
- 3D smoke, laptop, baseline (N 20/26, np 4): all 11 cases COMPLETED; seam np 8 against np 4:
  1.1e-12 (PASS). New case `cases/oscillatingDroplet3D` (mode-2 prolate spheroid, a = 1.1 R,
  b = c = R/sqrt(1.1), `signedDistanceEllipsoid`): initial volume error -3.23 % at N = 26,
  against -3.22 % for the sphere at the same N; initial mode-2 coefficient 7.4e-5 m (sphere
  2e-8 m), decreasing over the 20 steps.
- FOUND by the gate: every extension model with a linear solve died at step 1 in the eight
  droplet templates ("Entry 'UextFinal' not found"): their solver entry was `Uext`, not
  `"Uext.*"`. Fixed (inert for `velocityExtension none`, which never solves for Uext).
- FOUND by the gate's no-effect check: `VELOCITY_EXTENSION closestPoint` on the SL line gives
  CSVs identical to the baseline in every arm: the default `projectedFlux` trace ignores the
  extension. This is the second break named in CLAUDE.md "Extension without modification";
  step C4 of the plan removes it.
- FOUND: `writeDropletMetrics.H` sums the zero-set crossings over internal faces only, so the
  droplet shape columns miss the crossings on processor faces in a parallel run. Open. (CLOSED by b1798c3, 2026-09-27)

### 11.6 The pre-print's case figures (2026-09-26)

`docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_case_figures.py` writes one
TikZ figure per gate arm (exact 1D, 2D and 3D shear, the three 2D and three 3D droplets) with the
geometry, the initial condition and the boundary conditions taken from the case files. The 2D shear
streamlines are contours of the stream function; the 3D panels use an orthographic view (x right,
y into the page, z up) with computed hidden edges, and the camera keeps every front edge clear of
the droplet. Iterated on rendered pages: all half-width panels are at most 8.0 cm at scale 1 (two
per row in the 16.5 cm text width), labels on lines have white backgrounds. Section "Test cases"
of the article holds the three figures and the parameter table.

### 11.7 Phase C: the extension points are open, bit-identical (2026-09-26)

One commit, steps C1 to C5 of the plan. After it, a new source law, extension or trace flux needs
no solver edit.
- C1 `sdplsSource`: the selector tested the literal names `Rdiv` and `RdivStrictSp`; a virtual trait
  `usesDiscretization()` replaces them. `setTransportFlux()` gives a model the flux that advects
  psi; `eulerianAdvection`, `leiaRedistancedLevelSetFoam` and `leiaLevelSetTwoPhaseFoam` set it.
  The accessor `transportFlux()` stops with a fatal error when no flux was set.
- C2 new family `slSource` (base `none`, in `libleiaSemiLagrangian`, dictionary
  `levelSet.semiLagrangian.source`). `slAdvection::advect` calls it after the scheme step with the
  scheme's effective dt, and only when the model is active.
- C3 `leiaLevelSetTwoPhaseFoam` constructs a `velocityExtension` (`velocityExtensionFieldsEuler.H`)
  and transports psi with its flux in `div` and `Sp`. With `none` the flux is a copy of `phi`
  with the name `phi`, so the `div(phi,psi)` scheme still applies.
- C4 `leiaSemiLagrangianLevelSetTwoPhaseFoam`: `levelSet.semiLagrangian.traceFlux physical |
  extension` (default `physical`). `extension` traces reconstruct(phiExt). The solver warns when
  an extension is selected with `projectedFlux` and `physical`: that combination was the no-op
  that the gate found in 11.5.
- C5 `leiaSemiLagrangeLevelSetFoam`: the same `traceFlux` word and a `velocityExtension` root. For
  the reversed flow the root extends the t^n field (oscillation factor of t^n) for the old trace
  level and the t^(n+1) field for the new one.
- `Allwmake` refreshes the root lnInclude of `src/leiaLevelSet` before the libraries. MEASURED: a
  header in a new sub-directory (`semiLagrangian/source/slSource.H`) was not in it, and the
  solvers did not compile.

Gate, laptop: the reference worktree at 0c7079c (the state before Phase C, rebuilt) against the
new binaries. Each study runs up to `aggregate` in both trees; then every rendered dictionary,
every top-level CSV (tolerance 0, the two clock columns skipped) and every file of every written
time directory is compared byte for byte (`bitid2.sh`, `bitid_compare.py`).

| study (laptop sizes) | solver | covers | cases | result |
|---|---|---|---|---|
| `sdpls1Dstretch` (committed, serial, noSource/R/Rdiv, N 32-256) | `leiaLevelSetFoam` eulerian | C1 | 12 | identical |
| 1Dstretch np 4, N 32/64, `ADVECTION semiLagrangian` | `leiaLevelSetFoam` | C2 | 2 | identical |
| 1Dstretch np 4, N 32/64 | `leiaSemiLagrangeLevelSetFoam` | C2, C5 | 2 | identical |
| 2Dvortex np 4, R/beta | `leiaRedistancedLevelSetFoam` | C1 | 2 | identical |
| 2Dvortex np 4, noSource/R/Rdiv, N 32/64 | `leiaLevelSetFoam` eulerian | C1 | 6 | identical |
| stationaryDroplet2D np 4, noSource/R/Rdiv, N 32, 40 steps | `leiaLevelSetTwoPhaseFoam` | C1, C3 | 3 | identical |
| 2Dvortex np 4, reversed, cellCentred/projectedFlux, N 32/64 | `leiaSemiLagrangeLevelSetFoam` | C2, C5 | 4 | identical |
| 3Dshear np 4, N 24, T = 0.6 | `leiaSemiLagrangeLevelSetFoam` | C2, C5 | 1 | identical |
| stationaryDroplet2D np 4, N 64, 30 steps, projectedFlux/cellCentred x none/closestPoint | SL two-phase | C2, C4 | 4 | identical |
| translatingDroplet2D np 4, N 64, 100 steps | SL two-phase | C2, C4 | 1 | identical |
| stationaryDroplet3D np 4, N 30, 65 steps | SL two-phase | C2, C4 | 1 | identical |
| regression set, hex 2D: 2Dvortex and 2Dtranslation, N 32/64/128, np 4 | `leiaSemiLagrangeLevelSetFoam` | C2 | 6 | identical |
| regression set, hex 3D: 3Dshear N 16/24/32, T = 3, np 4 | same | C2 | 3 | identical |
| regression set, poly 3D: 3Dshear, 18 082 / 40 001 / 74 234 cells, T = 3, np 4, shared meshes | same | C2 | 6 | identical |

CORRECTED 2026-09-29: the 2Dtranslation arms of the hex 2D regression row ran the REVERSED
translation (11.19); they are void, and the row rests on its 2Dvortex arms.

53 cases, 0 differences. `leiaTestSdplsSource`: 89 passed in both trees, the same output. The
regression set ran the default arm at laptop sizes, not at the cluster sizes of the committed
`advConv*` studies: with `slSource none` the new call is skipped, so no arithmetic changes. The
polyhedral rung used `advect_bound_arm.sh` on one mesh per resolution (the same point digest in
both runs), so the mesher cannot enter the comparison.

### 11.8 Phase C6: the tokens of the new families, inert (2026-09-26)

- `cases/default.parameter`: `HL_RADIUS_CELLS 1`, `HL_M 2`, `HL_BETA 1`, `HL_DIRECTION levelSet`,
  `HL_SAMPLER stencilFit`, `SL_TRACE_FLUX physical`, `SL_SOURCE none`, `SL_SOURCE_BAND_CELLS 3`,
  `GC_LAW none`, `GC_STRAIN_WEIGHT none`, `GC_MU 1`, `GC_EPS 0.02`, `GC_SAT_C 1`, `SW_C_KAPPA 1.25`,
  `SW_DELTA_S 0.08`, `SW_P 5`, `SW_GAMMA 1.4722194895832` (artanh 0.9; the dossier's rounded 1.4722
  gives tanh = 0.89999), `SW_EPS_D 0`, `OMEGA_BETA 1`, `OMEGA_M 2`. The comment that listed
  `meshWaveExt` as a valid extension type now says `meshWave` (the registered name).
- The ten gate templates (1Dstretch, 2Dvortex, 2Dtranslation, 3Dshear, the three 2D and the three 3D
  droplets): the haloLimited entries at the end of `velocityExtension`, `traceFlux` after
  `traceVelocity`, the `source` block with its `law` at the end of `semiLagrangian`, and the `law`
  block in `sdplsSource`.
- `3Dshear` had no `velocityExtension` block, no `Uext` solver entries and no extension schemes, so
  no extension could be selected in the 3D shear arm. It has them now (copied from 2Dvortex).
- `1Dstretch` had no `traceVelocity` entry: its SL runs traced with the solver default
  `cellCentred`. The template carries the token; `cases/1Dstretch.parameter` keeps `cellCentred`,
  and both gates set the production trace `projectedFlux` as a line token of the SL line.

Gates:
1. Render diff: every committed config that uses one of the ten cases (250 configs), rendered in
   the reference tree (0c7079c) and in the new tree. Every change is an expected inert insertion or
   the `meshWave` comment: 250 of 250 PASS, 2 899 changed files, 126 487 added lines, 650 comment
   lines, no deleted or changed line. The classifier fails on two injected changes (a changed
   relTol, `traceFlux extension`). Thirteen configs first failed to render in one tree: a race of
   concurrent snakemake instances in the migration of `.snakemake`; they pass when rendered one at a
   time.
2. Bit identity, the Phase C set again on the new templates, new binaries (with D1 to D3 compiled
   in, not selected) against the reference: 53 of 53 cases, every CSV and every written field identical; every
   rendered dictionary differs from the reference only by the inert insertions. The polyhedral
   rung runs both binaries on the same C6 dictionaries (one mesh per resolution); the hex 3D
   rung covers the new 3Dshear dictionaries against the old ones.

### 11.9 Phases D1 to D3: the laws and the two source models (2026-09-26)

- D1 `libleiaGradientControl` (new library; `etc/leia-check-deps.py`, `Allwmake`, every solver's
  `EXE_LIBS`): the family `gradientControlLaw` (`none`, `linearQ`, `linearZ`, `cubicQ`, `cubicZ`,
  `twoThirdsZReg`, `saturatedLinearZ`, `softWall`) and its strategy `strainWeight` (`none`,
  `full`, `omega`). F = w(q) a + G(q, sigma). Every coefficient that a law reads is required;
  `softWall` refuses an even exponent. `boundedGradient` of the plan is not a class: the
  bounded-gradient equation with target 1 is `linearQ`, and a second name for the same law would
  break the one-name rule.
- D2 `sdplsGradientControl` (`sdplsSource` type `gradientControl`): F is the `nonLinearPart()`,
  so the discretization and mollifier strategies apply unchanged. It warns when a velocity
  extension and a strain weight both cancel the strain.
- D3 `slGradientControlSource` (`semiLagrangian.source` type `gradientControl`): psi <- psi
  exp(dt F) in the band |psi|/q <= bandCells h, F = 0 outside; |dt F| clamped at 30, counted and
  written (`slSourceF`, `slSourceClamp`). DEVIATION from the plan: q and n come from
  fvc::grad(psi, "gradPsiSource"), not from the geometry fit, because only the uncached quadratic
  fit provides fit derivatives; the case selects the scheme. The scheme keys come to the gate
  templates in the next commit, with their own render check.

Unit gates (laptop):
- `leiaTestGradientControlLaw` (new, no mesh): 109 passed. G(1) = 0 exactly for every law, the
  sign of G, dG/dq against a central difference, the equal q/z slopes, the soft-wall table of the
  dossier (3 to 10 percent, within 5e-4) and Theta* = 0.9, the strain weights, the SDPLS
  cancellation, RK4 against the logistic, z-logistic and implicit cubic solutions (RK4 order 4).
  Two injected law defects (a wrong linearZ slope, a flipped soft-wall sign) fail it.
- `leiaTestSdplsSource` (extended): 189 passed (89 before). gradientControl with law none and
  weight full equals R bit for bit; linearQ with mu 1 equals beta 1 (difference 0); every law and
  weight applies F psi at the exact affine values to 1e-10.
- `leiaTestSlSource` (new, `cases/slSourceUnit`, serial and np 4 with the band across the
  processor boundary x = 0): 49 passed on both. psi = psi^0 exp(dt F) in the band (240 cells) to
  1e-12, psi unchanged outside bit for bit, no sign change, the clamp at |dt F| = 50, the cell
  size, and current processor-patch values. A mutant without the halo update passes in serial
  and fails on 4 ranks.

### 11.10 Phase D4: the halo-limited extension, and the scheme keys (2026-09-26)

- `velocityExtension` type `haloLimited` (`velocityExtension/haloLimited/`) with four new strategy
  families, each runtime-selectable: `extensionTravel` (`capped`: S_R(d) = d [1 +
  (d/R)^(2m)]^(-1/(2m))), `extensionWeight` (`fractionReached`: w = (S_R/d)^beta),
  `extensionDirection` (`levelSet`: d = psi/sqrt(Q2), e = grad(psi)/sqrt(Q2), with the flat
  regularizer Q2 = q^2 + (1 - 4q^2)^3/4 for q < 0.5, zero above) and `extensionSampler`
  (`stencilFit`: the uncached quadratic value fit of the semi-Lagrangian line, one velocity
  component at a time). A linear sampler is not offered: at first order the extension is a
  weighted SDPLS source (dossier).
- Flux in the correction form phiExt = phi + w [u_h(Y) - u_h(x)].S_f, the mean of the two
  cells' models on a face; physical patches keep phi; on a coupled face both sides average their
  corrections in one orientation, so they carry exactly opposite values. R = radiusCells
  min(h_owner, h_neighbour), 0 < radiusCells <= 1, h from 1/deltaCoeffs over internal and coupled
  faces (no dependence on the decomposition). Written fields: hlWeight, hlDistance, hlReach,
  hlCorrection.
- The scheme keys of the new models in the ten gate templates: `gradPsiSource leastSquares`,
  `gradUSource Gauss linear` (D3), `gradPsiExtension leastSquares` (D4). Render diff against the
  reference with the C6 rules plus these three lines: 250 of 250 configs PASS.

Unit gate `leiaTestHaloLimited` (`cases/haloLimitedUnit`, serial and np 4, the planar interface
x = 0.01 next to the processor boundary x = 0, 240 coupled faces): 18 passed on both. S(0) = 0,
S'(0) = 1, S odd, |S| <= R up to |d| = 1e6 R, 1 - w = (beta/2m)(d/R)^(2m); for an affine velocity
the face correction equals -w S_R(d_f) (du/dx).S_f and Uext = U - w S_R du/dx to 1e-16 (the
quadratic models are exact); a uniform velocity leaves phi and U bit-identical; for the
divergence-free quadratic velocity of the dossier the correction equals w [u(Y) - u(x)].S_f with
the exact u to 2e-16; the correction field is exactly antisymmetric on the coupled faces. Two
floating-point defects found by the test and fixed: |S| exceeded R by one rounding step for
t > 1 (S is now sign(d) R (1 + t^(-2m))^(-1/(2m)) there), and (1 - w) U + w U was not U bit for
bit (Uext = U + w (u_h(Y) - U) now).

Final-state bit identity (all templates of C6 and D4, binaries with D1 to D4 compiled in and not
selected, against the reference 0c7079c): the Phase C set, 53 of 53 cases, every CSV and every
written field identical; the rendered dictionaries differ only by the inert insertions.

### 11.11 The 2D gate smoke with the six candidates, laptop (2026-09-26)

`SMOKE=1`, N 32/45/64, np 4, about 20 to 55 steps per case; baseline and S1, HL0, HL1q, HL1z,
HL2, FP0 on commit 952e593 (clean stamps): 49 studies, every case COMPLETED, no rule error. This is
the 4-rank check of CLAUDE.md before the cluster, not a measurement.
- The sources act (HL2 shear: max |dt F| up to 0.06 per step, no clamp) and every sample point
  stays inside R (max |Y - x|/R = 0.99999999).
- Seam check (np 1 and np 8 against np 4, column-scaled tolerance 1e-10): PASS for S1, HL0,
  HL1q, HL1z and HL2; FAIL for FP0 (E_VOL_ALPHA_REL differs by 1.4e-2 at np 1), as its candidate
  file predicted (closestPoint trusts halo Taylor data and falls back to a steady solve).
- The smoke verdicts are not results: N = 32 and 45 are below R/h = 10, and the horizons are 20
  to 55 steps. OBSERVATION for the gate to test: every extension candidate lowers the shear
  gradient band error at N = 64 (HL0 0.53, HL1q/HL1z 0.50, HL2 0.47, FP0 0.05 of the baseline) and
  raises the shear volume error about six times (HL0 6.8e-4 against 1.1e-4) with a lower order
  (0.54 against 1.27). S1 misses its target (0.83 > 0.8) and also raises the volume error.
- PRE-REGISTERED FOLLOW-UP (not run): with R = 1 h the extension velocity changes over one cell
  in the normal direction, so the trace error near |d| = R can be first order. Prediction: HL0
  with R = 2 h halves the finest-rung shear volume error of HL0 and keeps its gradient gain.
  It needs the cap radiusCells <= 1 raised (the local quadratic model is then extrapolated beyond
  its stencil). Falsified if the volume error does not fall by half.
- The 3D gate smoke (N 20/26, np 4; the seam arm at np 8), same candidates and binaries: 42
  studies, every case COMPLETED, no rule error; seam PASS for S1 and every HL candidate, FAIL for
  FP0 (as predicted). This is the 4-rank check for a later 3D launch, not a measurement.
- CORRECTED after this smoke: the gate summary scored `maxMagU` (the maximum of |U| over the cells,
  an L_inf norm) and gave it an order, against the rule that no L_inf norm is reported. It now
  scores `meanMagUPrime` (L1) and `l2MagUPrime` (L2) at T; the bounds check has no order. The
  smoke summaries were regenerated; the running cluster gate uses the old script, so its
  summaries are regenerated with the fixed script when every arm is done.
- The kinematic SL solver evaluated the extension three times per step (t^(n+1), t^n, t^(n+1));
  it now evaluates t^n first and t^(n+1) second. The HL0 shear smoke is bit-identical before and
  after (3 rungs, both CSVs at tolerance 0).

### 11.12 The 2D method gate runs on Lichtenberg (2026-09-26)

- New clone `/work/scratch/tm83tomy/leia-gcls` (branch `feature/gradient-controlled-level-set`,
  commit 4015404), its own binaries (`etc/leia-env.sh`); the other sessions' clones are not
  touched. Verified before the launch: the clone's solver on PATH, the stamps g4015404 without
  `-dirty`, the symbols of haloLimited, slGradientControlSource and softWall in the libraries,
  and a dry run of 70 jobs (49 arms, 7 exact1D checks, 7 summaries, 6 comparisons).
- Orchestrator `55042210` (`leia-gate2D`, partition `long`), `make gate GATE=methodGate2D
  CANDIDATES=S1+HL0+HL1q+HL1z+HL2+FP0 PROFILE=profiles/slurm` (baseline added). Ledger: the
  clone's `.my_jobs`; the child job ids are appended from the orchestrator's `.err`.
- Read-out when every arm is done: `studies/methodGate2D_summary/<candidate>/{summary,orders,
  vsBaseline}.csv` and `verdict.txt`; pull with `make pull-study STUDY=methodGate2D_summary`.
  The 3D gate runs only after a 2D PASS.
- FAILED FIRST LAUNCH: all 28 exact1D cases COMPLETED (4 rungs x 7 candidates), then every
  exact1d_check died on `import yaml`: the gate rules called the first python3 on PATH, which on
  Lichtenberg is miniconda without yaml (the laptop's has it, so the smoke passed). The gate
  stopped with nothing else to run (55042210 COMPLETED, rc 2). Fixed in f7067d9: the helpers run on
  snakemake's own interpreter, as the report rule of workflow/Snakefile already did. The clone
  pulled f7067d9 (only workflow, docs and STATUS changed since 4015404, so the binaries stay) and
  the gate was relaunched as `55043312`; the seven finished exact1D arms are not re-run.
- MOVED 2026-09-26 23:00, on the user's instruction to run in `/work/scratch/tm83tomy/leia` and not
  in `$HOME`. Checked first: every job of the account had WorkDir `/work/scratch/tm83tomy/leia-gcls`,
  and no file in `$HOME` had changed since the setup; `/work/scratch/tm83tomy/leia` was on `development`
  at c431677, idle since 2026-09-22, with no tracked change (67 untracked logs and folders). The
  `leia-gcls` gate (orchestrator `55043312`, 30 live children) was cancelled by id from its ledger
  (209 child ids). `/work/scratch/tm83tomy/leia` now has the branch `feature/gradient-controlled-level-set`
  checked out and its binaries are rebuilt from it; its untracked files and its `studies/` stay. To
  give the clone back to the development line: `git switch development && ./Allwmake`.
  Relaunched there as `55044205` (WorkDir `/work/scratch/tm83tomy/leia`); ledger: that clone's `.my_jobs`,
  lines tagged `[gcls session]`.
- 23:32: 21 of 49 arms done, no rule error; six exact1D checks PASS (baseline 2.5e-6, the same value
  as in `leia-gcls`). HL2's exact1D arm is complete but its check waits: the outer workflow runs at
  most 16 local jobs and 16 long `run_arm` jobs hold every slot, so HL2's other arms start late.
  Fixed for the next run in 778c06f (the checks and the summaries have priority); the running gate
  keeps the old order.
- Conda on the cluster: the base env (`~/miniconda3`, Python 3.14) has neither yaml nor numpy, and it
  is the first `python3` on PATH in the jobs. The gate helpers do not use it: they run on snakemake's
  own interpreter (module Python 3.11.14, yaml and numpy from `~/.local`).
- NOTE for the harvest: the report rule of every non-smoke arm copies its `<study>_errors.csv`
  into the theme's `data/tables`. Commit only the gate summary tables from there.

### 11.13 The best SL two-phase settings, the record corrected, the translating arm re-run (2026-09-27)

Request (user, 2026-09-27, before an 8-hour absence): find the best measured settings of the SL (the SL findings of 11.13 to 11.15 are summarised for the separate SL session in `docs/knowledge-base/sessions/sl-session-handover.md`, and the SL article now carries them; 2026-09-28)
method for the coupling to the two-phase Navier-Stokes equations (mass and momentum consistency
at high density ratio, the surface-tension force), apply them, and re-run the multiphase studies
on Lichtenberg.

What was done:

1. A read-only audit of METHOD.md, this file, `cases/default.parameter`, the case files and the
   study configs.
2. Both gates carry the best coupling explicitly (8a9b85a, the `twoPhaseCoupling` block). A render
   check proves that the block changes no case file of the stationary and oscillating arms: the
   generated case directories of `baseline` and `HL1z` are byte-identical before and after. The
   translating arm differs only in `endTime`.
3. The record is corrected: METHOD.md sections 4.1, 4.3, 6, 8.1 and 10, three comments in
   `cases/default.parameter` (the token set is unchanged, 144 of 144), and five places in this
   file, each marked CORRECTED or RETRACTED 2026-09-27.

Findings:

1. The closed-box VOID (section 0) removes most of the translating evidence:
   `matchedBDF2Translating2D` (the ddt pairing table), `massFluxComparison2D` (the rhoLENT
   translating numbers), `rhoDdtGate2D` and `volumeCorrectionTranslating2D` (the need for
   boundRho), `bestConfigTranslating2D` (the only none-against-cellCentreInverse comparison),
   `translatingClearOutlet2D`, and most likely `translatingLadder2D` (committed e8a6660 on the
   morning of the fix, no record of a run after 440107f, and its N = 128 tokens diverged after
   the fix in `translatingRepaired2D`).
2. No `cases/<case>.parameter` sets `CURVATURE_EXTENSION`. The global default is
   `cellCentreInverse` (c935883). The Popinet `none` is in each Popinet config, at density ratio 1.
3. `MASS_FLUX_BOUND_RHO` is active with rhoLENT, and its only measured basis is VOID.
4. METHOD.md sections 4.1 and 10 still described the configuration of 2026-07-31.

The best configuration as it stands (the `twoPhaseCoupling` block of both gates):

| token | value | evidence | status |
|---|---|---|---|
| `MASS_FLUX` | rhoLENT | `rhoLENTStationary2D`: residual +1.0 / -22 / +0.1 % at N = 32 / 64 / 128 | stationary only; no translating benefit shown after the fix |
| `MASS_FLUX_ALPHAF_SOURCE` | donorPlane | the author's instruction (2026-09-01) | no measurement |
| `MASS_FLUX_BOUND_RHO` | true | none valid | OPEN |
| `MOMENTUM_DDT_SCHEME`, `RHO_DDT_SCHEME` | backward, backward | the BDF2 rule; the matching argument | the pairing tables are VOID |
| `MOMENTUM_DIV_SCHEME` | upwind | inert on the 3D stationary droplet; the longest-lived arm of `translatingRepaired2D` (9987 steps) | every translating arm diverges |
| `VISCOSITY_FACE_MODEL` | alg_lin (2D) | the viscosity-contrast ladder at mu1/mu2 = 1000 | 3D templates have no token (OPEN) |
| `N_OUTER_CORRECTORS`, `PSI_OUTER_CORRECTORS` | 3, yes | 3 frozen equal 12 re-advected to 3-4 digits | measured |
| `CAPILLARY_FORCE_CENTRING` | endStep | the linear analysis; a 20-step smoke | no coupled measurement |
| `SURFACE_TENSION_FORCE`, `FACE_CURVATURE_SOURCE` | reconstructedCurvature, model | non-gradient force content at order +2.01 | constant curvature only |
| `CURVATURE_EXTENSION` | cellCentreInverse | the 2D stationary ladders: residual 4.60 / 3.81 / 1.55 / 1.49x lower at N = 32-256 | translating undecided; oscillating at N = 128 only |
| `CURVATURE_INVERSE_GAUSSIAN` | yes | h^1.95 against h^1.02 on the 3D sphere | measured |
| `CAPILLARY_DT_COEFF` | 0.010861 | 0.2323 of the Brackbill limit at every rung | measured |
| `PSI_FILTER`, `VOLUME_CORRECTION` | none, noVolumeCorrection | the repository rules | — |

The translating arm. The 0.1 s runs of `55044205` diverged in every candidate. The time of the
last CSV row, per rung (N = 100 / 142 / 200):

| candidate | N = 100 | N = 142 | N = 200 |
|---|---|---|---|
| baseline | 0.0904 | 0.0942 | 0.0772 |
| S1 | 0.0574 | 0.0433 | 0.0593 |
| HL0 | 0.0870 | 0.0734 | 0.0619 |
| HL1q | 0.0996 | 0.0724 | 0.0648 |
| HL1z | 0.0917 | 0.0718 | 0.0618 |
| HL2 | 0.0467 | 0.0365 | running (0.0185 at 01:00) |
| FP0 | 0.0511 | running (0.0432) | running (0.0264) |

1. These are valid results at the pre-registered horizon of 0.1 s, not a wrong setup. They show
   that the late instability is not a property of the baseline alone. S1 and HL2 diverge before
   0.05 s at one rung or more. The directories stay, renamed
   `methodGate2D_<candidate>_translating_endTime0p1_20260927`.
2. The gate's translating arm ends at 0.05 s with `cellCentreInverse`. The `none` of 8a9b85a is
   reverted, because its evidence is VOID. The arm then differs from the 0.1 s runs only in
   `endTime`, so each re-run must reproduce the first 0.05 s of its 0.1 s run byte for byte
   (`compare_metrics_csv.py`, the two timing columns skipped).
3. Stated as a change after seeing data: the horizon was pre-registered as 0.1 s and was changed
   to 0.05 s after the baseline diverged. The verdict at 0.05 s answers one question: does a
   candidate degrade the translating droplet while the production method is stable? The late
   instability stays open. Not tried after the fix: `cellCentreInverse` against `none` on a
   matched setup, a longer box, the RK2 foot integrator, midpoint force centring, and the
   semi-implicit capillary force.

### 11.14 Two parallel defects of the coupled SL solver; the Eulerian rhoLENT port (2026-09-27)

Found while preparing the Eulerian port, which reuses the SL mass-flux code. Every number below
is from the laptop, the gate's baseline translating case (N = 100, R/h = 10, END_TIME as stated),
binaries of the named commits.

**Defect 1, the solution: the face density on a coupled face (fixed in 28d13f0).**
`computeFaceAreaFractions` filled every boundary face from its local cell. On a processor face
each rank used its own cell's plane, whatever the flux direction, so the two sides disagreed on
alpha_f and rho_f*phi, and mass was not conserved across the seam where the interface crossed it.

1. Measured, N = 64, np 4, one step: 4 processor faces with alpha_f different on the two sides,
   rho_f up to 90 % apart; phi equal on both sides. After the fix: 0 faces differ.
2. np 4 against serial (column-scaled maximum difference, velocity metrics): 1e-5 to 5e-4 by
   t = 0.033 before the fix; 1e-8 to 1e-10 through step 5000 after it.
3. The per-rank mass residual (1e-10 to 1e-13) did not see the defect: each rank's auxiliary
   density balances its own flux.

**The late translating instability is not a decomposition artefact.** To t = 0.1:

| run | diverged at step | t [s] |
|---|---|---|
| serial | 7139 | 0.0775 |
| np 4, before the fix (the cluster baseline: the same step) | 8331 | 0.0904 |
| np 4, with the fix | 8000 | 0.0868 |

The runs agree until the growth starts near t = 0.06-0.065; from there any difference is
amplified.

**Defect 2, the diagnostics: internal-face loops in writeDropletMetrics.H (fixed in b1798c3).**
The zero-set sample, the curvature/gradient band and the face statistics used internal faces
only; on several ranks the processor faces were missing. np 4 against serial, to t = 0.05, with
defect 1 already fixed:

| column | before | after |
|---|---|---|
| zeroSetRadialL2 (gate: shape) | 3.42e-02 | 1.20e-07 |
| gradPsiL2ErrorBand (gate: gradient band) | 3.82e-02 | 1.24e-08 |
| kErrL2Band (gate: curvature) | 7.59e-02 | 2.11e-08 |
| m2Amplitude (gate: oscillating shape) | 6.31e-01 | 1.92e-07 |

Serial runs are unchanged by both fixes: bit-identical to the previous binary over 4604 steps,
every column (`compare_metrics_csv.py --tol 0`).

**Consequence.** Every SL two-phase result on more than one rank before 2026-09-27 carries both
defects: a decomposition error of about 1e-4 in the solution before any instability, and a
decomposition bias of several percent in the shape, gradient-band, curvature and mode-2 columns.
OPEN, author decision: which earlier parallel studies to re-run or void. The running 2D gate
(`55044205` and the five translating re-runs) uses the pre-fix binaries; it completes as the
pre-fix record, then the whole gate re-runs with the fixes (PRESERVE=1), and the two are
compared for verdict changes.

**The gate now checks the coupled solver's decomposition (aaae627).** `translatingSeamNp1` runs
the coarsest translating rung in serial and compares the error-vector columns of the droplet CSV
with the np 4 run at tol 1e-5 (3.5e-7 at most after the fixes). A reference that did not complete
is NOT_COMPARABLE, not a FAIL. The seam option now works on any arm.

**The Eulerian two-phase solver gets the SL mass flux (c094bd8).** It kept rho, rhoPhi and muf at
their t = 0 values for the whole run. The SL mass-flux code is now shared, moved verbatim into
`createMassFluxFields.H`, `updateFaceDensity.H` and `updateMassFlux.H` (the SL files reconstruct
byte for byte). The Eulerian solver reads the same `levelSet.massFlux` dictionary with the same
defaults, rebuilds the face density from the new interface, solves the rhoLENT auxiliary density
on every outer corrector, resets rho after the PIMPLE loop, and logs the density-bound numbers.
`PSI_OUTER_CORRECTORS yes` rebuilds its interface on every outer corrector. Laptop, np 4, 1000
steps (t = 0.0109 s):

| run | travelled fraction | L2 \|U-U0\| [m/s] | L1 \|U-U0\| [m/s] |
|---|---|---|---|
| Eulerian, new, ratio 840 | 1.00012 | 4.85e-4 | 2.54e-4 |
| Eulerian, previous binary, ratio 840 | 0.420 | 1.55e-1 | 3.05e-2 |
| Eulerian, new, ratio 1 | 1.00078 | 1.04e-4 | 3.59e-5 |
| Eulerian, previous binary, ratio 1 | 1.00339 | 2.18e-4 | 9.10e-5 |

np 4 against serial (new): U, alpha, psi and p_rgh equal to 1e-8 relative. Stationary droplet,
np 4: displacement 2.5e-13 m, L2 |U| 8.0e-5 m/s. Relative mass residual at ratio 840: 1.2e-10.
The previous Eulerian binary moves the heavy droplet at 42 % of the stream: every Eulerian
two-phase result with a moving interface and a density contrast before this commit is wrong
physics (the OPEN decision on the frozen-rho studies, plan item D-m). Two steps remain before
the Eulerian line can join the coupled gate arms:

1. The curvature pipeline. The Eulerian production curvature is the older one
   (`meanCurvatureClosestPoint` plus the `stabilizedFootPointFace` delivery,
   `reconstructedCurvatureFields.H`); it does not read `curvatureExtension`, so it cannot run
   `cellCentreInverse`, the SL best. It needs the same curvature-extension dispatch.
2. The droplet CSV (plan F2). `writeDropletMetrics.H` needs two small adapters: the CSV name
   (hard-coded to the SL solver) and the fit object of one research diagnostic
   (`slAdv->reconstruction()`); everything else it reads exists in the Eulerian solver now.

Also noted, not fixed: `rhoClipFraction` counts round-off clips in pure cells (37-48 % of the
cells at a clip L1 of 1e-12), so as an error metric it measures round-off; the research
diagnostics `A4h*`, `A8h*` and `driver*` still use internal faces only.

**A setup defect of the 3D translating case (fixed 2026-09-27).** `cases/translatingDroplet3D`
had no `dropletReferenceVelocity` entry, so the metrics writer used (0 0 0): the disturbance
metrics (`maxMagUPrime`, `meanMagUPrime`, `l2MagUPrime`) reported the translation itself and the
zero-set and centroid errors reported the displacement. Seen in the 3D gate smoke:
`meanMagUPrime = 0.0500 = U0`. The solution does not read the entry. The 2D translating case and
both Popinet cases have it. Now the 3D template renders `dropletReferenceVelocity (0.05 0 0)`.
VOID for these columns: `traceTranslating3Dhex` and `traceTranslating3Dpoly_r10p0/r12p7/r15p8`
(no conclusion in STATUS.md, METHOD.md or the docs cites them). The 3D gate has not run.

**The late translating instability: four one-change discriminators (laptop, 2026-09-27).** The
gate's baseline translating case, N = 100 (R/h = 10), np 4, to t = 0.1 s, binaries with both
parallel fixes. One entry changed per run. The read-out is the divergence time against the
decomposition scatter of the unchanged case (serial 0.0775 s, np 4 0.0868 s with the fixes,
0.0904 s without). A shift inside that scatter is not a result. ONE resolution only: these are
indicators for the choice of the next ladder, not results.

| change | result | t [s] | reading |
|---|---|---|---|
| none (the reference, np 4) | DIVERGED step 8000 | 0.0868 | — |
| `curvatureExtension none` | DIVERGED step 6404 | 0.0695 | 20 % earlier: cellCentreInverse delays the instability |
| `footIntegrator rk2` | DIVERGED step 7759 | 0.0842 | inside the scatter: no effect |
| `capillaryForceCentring midpoint` | DIVERGED step 5467 | 0.0593 | 32 % earlier: worse |
| density ratio 1 (499.695 both) | COMPLETED | 0.1 | the instability needs the density contrast |

Consistent with section 0 and `translatingRepaired2D` (`none` at N = 128 diverged at
t = 0.063-0.075 s, ratio 1 completed). The divergence onset (growth from t = 0.06 s) is when the
droplet centre is at 5.5 mm and its leading edge 3.5 mm (35 cells) from the outlet. Not yet
tested: a longer box (the distance to the outlet), and the semi-implicit capillary force.

### 11.15 The pre-fix verdicts, the outlet, and the fixed 2D gate (2026-09-27)

**The pre-fix 2D gate** (`55044205` plus the five translating re-runs, pre-fix binaries, summaries
in `methodGate2D_summary_pre-20260927-020856`): every candidate's verdict is FAIL. The kinematic
shear arm alone fails all of them. The coupled-solver defects of 11.14 do not touch the kinematic
solver, whose seam checks pass at 1e-12:

| candidate | shear shape error at N = 136 | ratio to baseline (6.19e-4) | shape order (baseline 3.07) |
|---|---|---|---|
| HL0 (halo-limited extension, no source) | 8.49e-3 | 14 | 1.39 |
| S1 (soft wall, no extension) | 0.254 | 410 | 0.06 |
| HL1q (extension, linear q law) | 0.107 | 170 | 0.02 |
| HL1z (extension, linear z law) | 1.65 | 2700 | — |

The extension and the sources degrade the transport of the interface itself. The coupled arms
add divergences where the baseline completes (for example HL1q and HL1z at stationary N = 142).
The fixed gate re-runs the kinematic arms with unchanged kinematic code; they must reproduce the
pre-fix ones byte for byte.

**The translating re-runs to 0.05 s** (pre-fix binaries): the baseline completes all three rungs,
so criterion 1 can hold. All 15 re-runs reproduce the first 0.05 s of their 0.1 s runs byte for
byte (`compare_metrics_csv.py --tol 0`, the pre-registered check of 11.13), S1 at N = 142 to the
identical divergence row (6741).

**The late translating instability: the outlet triggers the fast divergence; a slower interior growth remains** (laptop, N = 100, fixed binaries). CORRECTED 2026-09-27 03:40: this heading first said "the late translating instability is the outlet". I was wrong: the 20 mm box to t = 0.25 s diverges too, later (see the end of this paragraph). The same
case in a box twice as long in x (20 mm, the same h, the outlet 10 mm further away) COMPLETES
t = 0.1 s, and its spurious current decays:

| t [s] | 0.02 | 0.05 | 0.06 | 0.065 | 0.07 | 0.08 |
|---|---|---|---|---|---|---|
| L2 \|U-U0\|, 10 mm box (np 4) | 1.04e-3 | 1.95e-3 | 1.79e-3 | 4.90e-2 | 7.56e-2 | 1.35e-1 |
| L2 \|U-U0\|, 20 mm box (np 8) | 7.39e-4 | 1.38e-3 | 1.22e-3 | 1.09e-3 | 9.10e-4 | 6.73e-4 |

(The 20 mm box has twice the volume, so its volume-weighted L2 norm is about 1/sqrt(2) of the
10 mm value at equal local levels; the early ratio is 0.71.) In the 10 mm box the growth starts
when the leading edge is 3.5 mm (35 cells) from the outlet. This is the geometric event that the
section "Is the interface still inside the domain?" of CLAUDE.md warns about; it needs the density
contrast (ratio 1 completes). The gate's 0.05 s horizon ends before it. OPEN, author decision:
a longer translating box for the gates and ladders (a new length token with the current box as
its default).

Confirmed at N = 142 (R/h = 14.2, np 8, fixed binaries): the 10 mm box DIVERGES at step 14198
(t = 0.0911 s, volume error 0.48 at the last row); the 20 mm box COMPLETES t = 0.1 s (volume error
2.2e-3). L2 |U-U0|:

| t [s] | 0.02 | 0.05 | 0.06 | 0.07 | 0.08 | 0.09 | 0.1 |
|---|---|---|---|---|---|---|---|
| 10 mm box | 6.80e-4 | 1.26e-3 | 9.66e-4 | 7.04e-2 | 1.29e-1 | 6.46e-1 | — |
| 20 mm box | 4.81e-4 | 8.95e-4 | 6.85e-4 | 4.28e-4 | 9.64e-4 | 1.56e-3 | 1.80e-3 |

A slower growth remains in the 20 mm box at N = 142 after t = 0.07 s: 4.3e-4 to 1.8e-3 in 0.03 s
(about 46 1/s), with the droplet more than 11 mm from the outlet. That is an interior growth,
about 70 times slower than the outlet-driven one; at N = 100 the 20 mm box decays over the same
window. It is the question a long-box ladder (three rungs, 0.1 s or longer) must answer next.

The long-box ladder to t = 0.1 s (laptop, fixed binaries; np 8, 8, 16 at N = 100, 142, 200): all
three rungs COMPLETE. At t = 0.1 s:

| metric | N = 100 | N = 142 | N = 200 | pairwise orders |
|---|---|---|---|---|
| zeroSetRadialL2 [m] | 1.38e-4 | 7.08e-5 | 3.99e-5 | 1.90, 1.68 |
| centroidError [m] | 1.95e-4 | 1.00e-4 | 5.62e-5 | 1.89, 1.69 |
| phaseVolumeRelError | 1.30e-2 | 2.19e-3 | 2.20e-3 | 5.1, 0.0 |
| gradPsiL2ErrorBand | 7.70e-2 | 7.80e-2 | 4.53e-2 | -0.04, 1.58 |
| l2MagUPrime [m/s] | 5.14e-4 | 1.80e-3 | 7.36e-4 | -3.6, 2.6 |
| kErrL2Band [1/m] (exact 1000) | 15.4 | 21.6 | 20.2 | -1.0, 0.2 |

Shape and centroid converge near order 1.7-1.9. The spurious current stays bounded at every rung
but is not monotone in h (the N = 142 bump after t = 0.07 s); the curvature error does not converge
(1.5-2.2 %).

**The 20 mm box to t = 0.25 s, N = 100: DIVERGED at step 19590 (t = 0.2127 s).** L2 |U-U0| (and the
droplet centroid x):

| t [s] | 0.10 | 0.12 | 0.14 | 0.16 | 0.18 | 0.21 |
|---|---|---|---|---|---|---|
| L2 \|U-U0\| | 5.1e-4 | 7.7e-4 | 1.4e-3 | 3.4e-3 | 8.5e-2 | 5.9e-1 |
| centroid x [mm] | 7.7 | 8.7 | 9.8 | 10.9 | 13.4 | 14.7 |

The growth starts near t = 0.10 s at about 30 1/s, with the droplet more than 8 mm from the outlet,
and turns explosive after t = 0.16 s (the centroid then runs ahead of U0 t: the droplet breaks up).
So there are TWO phenomena: (1) a fast growth that the outlet triggers in the 10 mm box (a factor 27
in 0.005 s at t = 0.06 s, N = 100), and (2) a slower interior growth that the long box only
postpones. The gate's 0.05 s horizon ends before both. N = 142 in the 20 mm box: DIVERGED at step
34187 (t = 0.2194 s), the same pattern.

**A geometric signature (04:10).** The time at which L2 |U-U0| first exceeds 10 times its t = 0.03 s
level, and where the droplet is then:

| run | t [s] | centroid x [mm] | box centre [mm] | centroid - centre [mm] |
|---|---|---|---|---|
| 10 mm box, N = 100, serial | 0.0630 | 5.76 | 5.0 | +0.76 |
| 10 mm box, N = 100, np 4 | 0.0630 | 5.76 | 5.0 | +0.76 |
| 10 mm box, N = 142, np 8 | 0.0656 | 5.72 | 5.0 | +0.72 |
| 20 mm box, N = 100, np 4 | 0.1608 | 10.98 | 10.0 | +0.98 |
| 20 mm box, N = 142, np 8 | 0.1611 | 11.04 | 10.0 | +1.04 |

In every run the explosive growth starts when the droplet's TRAILING edge (centroid - R, R = 1 mm)
crosses the middle of the box, within 0.3 mm, in serial and in parallel, at both resolutions and
both box lengths. The jump time reproduces to 0.2 % between N = 100 and 142 in the 20 mm box: a
geometric event (CLAUDE.md, "Is the interface still inside the domain?"), not a growth that starts
at t = 0. No case file places anything at x = L/2 (the inlet fixes U and alpha, the outlet fixes
p_rgh, psi is zeroGradient everywhere). So the position in the box and the elapsed time are
confounded in these runs. The discriminator, with its prediction written before the result: the
20 mm box, N = 100, droplet starting at x = 5 mm instead of 2.5 mm. If the event is tied to the box
centre, the jump comes with the centroid near 11 mm, at t near 0.12 s; if it is tied to time or
travel, near t = 0.16 s with the centroid near 13 mm.

RESULT (04:40): the jump came at t = 0.1372 s with the centroid at 12.15 mm (7.15 mm travelled);
the run diverged at t = 0.1917 s. NEITHER prediction holds, so the box-centre coincidence of the
table above is NOT supported. Up to t = 0.12 s the two runs have IDENTICAL spurious-current
histories (7.4e-4, 8.0e-4, 1.2e-3, 6.7e-4, 5.1e-4, 7.7e-4 at t = 0.02 ... 0.12 s): the slow phase is a
function of time only, as translation invariance requires. The explosive phase then comes
earlier for the droplet that is further downstream (its leading edge 6.9 mm from the outlet at the
jump, against 8.0 mm for the later jump of the droplet started at 2.5 mm). The reading that fits
every run so far: a slow interior growth in time, and an amplification that grows as the droplet
approaches the outlet; the jump happens when the two together cross a threshold. In the 10 mm box
the outlet is close from the start, so the jump comes at t = 0.063 s. Not yet tested: a far longer
box (the interior growth alone), an outlet condition other than fixed p_rgh with zeroGradient U.

**The fixed 2D gate:** orchestrator `55048916` (section 5). While it runs, the cluster clone stays at 1150e68 (no pull: new cases would record a commit the binaries do not carry); its STATUS.md is refreshed from GitHub without a pull (`git show origin/<branch>:STATUS.md > STATUS.md`; STATUS.md is outside the paths that mark a case dirty). Before the next pull there: `git checkout -- STATUS.md`.

**Two corrections of the gate's scoring, made after the fixed gate started (so its orchestrator
still scores the old way; the summaries and verdicts are recomputed with the corrected script
after it finishes):**

1. `rhoClipFraction` is reported, not scored (4bff922): it counts round-off clips.
2. The oscillating arm no longer scores `meanMagUPrime` and `l2MagUPrime` (this commit). There
   |U - U_ref| is the physical oscillation, not a spurious current, and a candidate that damps
   the oscillation more scored as better. Plan section 5.4 pre-registered the oscillating arm's
   entry as the damping rate; the implementation had deviated from it. The two norms are now
   reported as `oscL2MagU` and `oscMeanMagU`. The login node `lcluster5` stopped
answering at about 02:00; `lcluster3` and `lcluster4` answered and were used.

**Fixed gate, first checks (02:52):** the kinematic arms (exact1D, shear, seamNp1, seamNp8) of the
fixed gate reproduce the pre-fix gate byte for byte: 177 CSV pairs identical at tolerance 0, as
predicted (the fixes do not touch the kinematic solver). The baseline's `translatingSeamNp1`
(serial) against its np 4 run, the error-vector columns of the droplet CSV: at most 3.47e-7
(`meanMagUPrime`), 1.2e-7 (`zeroSetRadialL2`), 2.1e-8 (`kErrL2Band`), the same numbers as on the
laptop to every printed digit; the check passes at its tolerance 1e-5. The baseline translating
arm completes all three rungs.

**Fixed against pre-fix baseline (the whole vector, the value at T; 03:40).** Most entries agree to
0.4 % or better; the kinematic shear arm to 0.0 %. The seam defects mattered where the processor
boundaries cross the droplet:

| arm, N | metric | pre-fix | fixed | change |
|---|---|---|---|---|
| stationary 100 | shapeError | 3.90e-4 | 4.27e-4 | +9.4 % |
| stationary 100 | meanMagUPrime | 1.48e-8 | 1.35e-8 | -8.7 % |
| oscillating 100 | l2MagUPrime (the physical oscillation) | 2.41e-3 | 1.22e-3 | -50 % |
| oscillating 200 | gradientBandError | 3.41 | 5.13 | +51 % |
| oscillating 200 | volumeError | 8.18e-3 | 3.04e-3 | -63 % |

The baseline's oscillating arm is unhealthy at N = 200: the gradient band error is 5.1 (0.10 and
0.14 at N = 100 and 142), so |grad psi| in the band is far from 1 by t = 0.1 s at the finest rung.
The oscillating arm's orders and verdicts must be read with that in mind.

**Mid-gate pull (04:37), and why it was safe.** The no-pull rule protects case provenance: a case
materialized after a pull records a commit the binaries do not carry. Checked first: every one of
the 56 arms had materialized its cases, and every case records 1150e68 (0 exceptions). Between
1150e68 and the pulled 35a8878 only `workflow/scripts/make_gate_summary.py` changed in `src`,
`applications`, `workflow`, `cases` and `config`. So the orchestrator's remaining summaries and
comparisons now use the corrected scoring (4bff922, 4ef98db), and the baseline summary was
recomputed with it (the oscillating rows carry `oscL2MagU`; `l2MagUPrime` is empty there).

**Finish estimates (04:36).** HL0, HL1q and HL1z within about 1 h; S1 by about 08:20; HL2 by about
11:00. FP0's stationary and oscillating N = 200 run at 0.19 and 0.14 steps/s and need 29 h and
41 h: they will time out at the 10 h job limit. FP0's verdict does not depend on them: its shear
arm (bit-identical to the pre-fix gate) fails at 34 times the baseline's shape error, and its own
seam check fails (closestPoint is decomposition-dependent, 0.45).

**The interior growth without the outlet: a 40 mm box (N = 100, laptop, np 8, to t = 0.3 s).** The
outlet stays more than 17 mm from the droplet. The run COMPLETES, but it degrades:

| t [s] | 0.10 | 0.14 | 0.18 | 0.22 | 0.26 | 0.30 |
|---|---|---|---|---|---|---|
| L2 \|U-U0\| (40 mm box; 1/sqrt(2) of the 20 mm scale) | 3.6e-4 | 1.0e-3 | 2.2e-3 | 4.4e-3 | 6.9e-3 | 8.3e-3 |
| phaseVolumeRelError | 1e-2 | 4e-3 | 1e-2 | 4e-2 | 8e-2 | 1e-1 |
| centroid x - (x0 + U0 t) [mm] | 0.2 | 0.3 | 0.6 | 1.3 | 2.5 | 3.9 |

The spurious current grows from t = 0.10 s at a slowing rate (about 26 1/s at first, 5 1/s at the
end), the droplet volume drifts by 10 % (CORRECTED 2026-09-29: 15 %, 1.5e-1 at t = 0.3 s in the second table below and 0.147 in the archived history; the table above rounds it) and the droplet runs 3.9 mm ahead of the stream. So the
interior degradation is real and needs no outlet; the outlet only turns it into an explosive
divergence. For the method this is the open defect of the translating droplet at the water/air
density ratio.

The same 40 mm box at N = 142 (np 16, to t = 0.3 s) COMPLETES with a smaller degradation:

| t [s] | 0.10 | 0.18 | 0.22 | 0.26 | 0.30 |
|---|---|---|---|---|---|
| L2 \|U-U0\|, N = 100 | 3.6e-4 | 2.2e-3 | 4.4e-3 | 6.9e-3 | 8.3e-3 |
| L2 \|U-U0\|, N = 142 | 1.3e-3 | 1.1e-3 | 1.6e-3 | 2.4e-3 | 3.6e-3 |
| volume error, N = 100 | 1.3e-2 | 1.2e-2 | 3.6e-2 | 8.1e-2 | 1.5e-1 |
| volume error, N = 142 | 2.2e-3 | 9.1e-3 | 1.9e-3 | 9.6e-3 | 2.8e-2 |
| centroid lead [mm], N = 100 | 0.19 | 0.63 | 1.32 | 2.48 | 3.86 |
| centroid lead [mm], N = 142 | 0.10 | 0.58 | 0.82 | 1.21 | 1.82 |

At t = 0.3 s the degradation is 2.3x (spurious current), 5.4x (volume) and 2.1x (lead) smaller at
N = 142 than at N = 100: it shrinks with h, as a discretization error does, not as a
resolution-independent instability. TWO rungs only: no order is stated. The third rung (N = 200,
160 000 cells, 78 000 steps) belongs on the cluster. With the outlet more than 17 mm away the
translating droplet is usable over 0.3 s at N = 142; near the outlet it is not.

**Fixed-gate verdicts so far (05:32).** HL1q: FAIL. HL1z: FAIL (43 FAIL lines). Both fail on every
arm: the shear arm (as before the fix), the stationary arm at N = 200 (HL1q shape error 3.8 R and
volume error 76 %; HL1z 3.7 R and 84 %: these runs COMPLETE, but with the droplet destroyed, so
the completion count alone would mislead), the translating and the oscillating arm. Waiting: HL0
and S1 (oscillating finest rungs), HL2, FP0.

**The baseline's oscillating arm goes unstable at the finest rung (fixed gate, 06:20).** The band
gradient error `gradPsiL2ErrorBand` (L2 of |grad psi| - 1 in the band) over time:

| t [s] | 0.01 | 0.03 | 0.05 | 0.07 | 0.08 | 0.09 | 0.10 |
|---|---|---|---|---|---|---|---|
| N = 142 | 0.036 | 0.066 | 0.076 | 0.099 | 0.111 | 0.123 | 0.140 |
| N = 200 | 0.054 | 0.087 | 0.159 | 0.386 | 0.581 | 1.24 | 5.13 |

At N = 142 the drift grows about linearly; at N = 200 it grows exponentially after t = 0.05 s, and
at t = 0.10 s the velocity norm jumps tenfold (8.1e-2). The run COMPLETES the 0.1 s horizon on the
edge of a divergence, so its values at T (period 8.06 ms against 10.0 and 9.82 ms at the coarser
rungs; the Celik classification of period and damping: divergent) cannot carry a verdict. This
drift of |grad psi| without redistancing is exactly the defect the gradient-control candidates
exist to remove; the candidates fail for other reasons (the shear transport). OPEN, author
decision: the oscillating arm's horizon (as for the translating arm), and a baseline check that
reads the gradient drift over time, not only at T.

**HL0 (halo-limited extension alone): FAIL, and it destroys the level set in the coupled arms
(fixed gate, 06:45).** HL0 is the dossier's core idea (protect |grad psi| by the extension, no
source). Values at T (baseline in brackets):

| arm | N = 100 | N = 142 | N = 200 |
|---|---|---|---|
| oscillating, gradient band error | 7.5 (0.10) | 2.2e6 (0.14) | 2.0e7 (5.1) |
| oscillating, volume error | 0.74 (7.6e-3) | 10.4 (7.3e-3) | 6.7 (3.0e-3) |
| translating, shape error / R | — | — | 0.196 (3.1e-3) |
| shear (kinematic), shape error | — | — | 8.5e-3 at N = 136 (6.2e-4), order 1.39 (3.07) |

The oscillating runs COMPLETE with garbage fields. The stationary arm is the exception: only
volume-error regressions at the 1e-7 level. The ladder of CLAUDE.md places the defect: the exact
unit gates of D4 pass (zero normal strain to 1e-12 on affine U, uniform U bit for bit), and the
first rung above them, KINEMATIC transport (shear, no force), already fails at order 1.39. So the
defect is in the kinematic path of the extension, before any coupling. First suspect, named as a
risk in the plan (section 9): the `stencilFit` sampler (a least-squares fit of U; exact only for
affine U). Next step, the cheapest discriminator: the shear arm with HL0 at two samplers or radii,
and HL0 on the exact-velocity 1D and 2D kinematic cases with the flux correction written out.

**S1 and HL2: FAIL (fixed gate, 12:14).** S1 (soft wall, no extension): 19 FAIL lines; the
translating N = 142 rung diverges before 0.05 s, and its `translatingSeamNp1` check FAILS: serial
and np 4 differ by 76 % in `l2MagUPrime` (the baseline: 3.5e-7). The soft-wall source makes the
coupled solution sensitive to round-off, so its result depends on the decomposition; this is
the mesh-noise-floor diagnostic of CLAUDE.md, applied to the decomposition. HL2 (extension plus
soft wall): 15 FAIL lines; all three translating rungs diverge; its seam check is NOT_COMPARABLE
(the np 4 reference diverged). Every candidate of the first campaign now FAILS the 2D gate: S1,
HL0, HL1q, HL1z, HL2, and FP0 (its shear arm and seam check). The 3D gate therefore does not run.
Remaining: FP0's oscillating arm and stationary N = 200, which time out at the 10 h limit.

**The fixed 2D gate is complete (orchestrator `55048916`, 12:57, rc 2: FP0's three TIMEOUT
solves, 10 h limit).** Every summary and verdict recomputed with the corrected script (13:00,
HEAD c718abb, no source change since the build 1150e68). Verdicts: every candidate FAILS.

| candidate | verdict | target ratio | regressions | orders | seam | completion |
|---|---|---|---|---|---|---|
| FP0 | FAIL | 0.939 | 3 | 8 | 5 | 5 |
| HL0 | FAIL | 0.828 | 16 | 13 | 0 | 0 |
| HL1q | FAIL | 0.204 | 22 | 18 | 0 | 0 |
| HL1z | FAIL | 6.490 | 23 | 19 | 0 | 0 |
| HL2 | FAIL | 0.754 | 7 | 5 | 0 | 3 |
| S1 | FAIL | 0.707 | 12 | 5 | 1 | 1 |

Coupled seam check (`translatingSeamNp1`, tol 1e-5): baseline 3.5e-7, HL0 4.7e-6, HL1q 4.9e-7, HL1z
5.3e-7 PASS; S1 0.76 FAIL; FP0 1.11 FAIL; HL2 NOT_COMPARABLE. The trade-off the gradient-control
sources show on the shear arm: HL1q lowers the band gradient error of the baseline (1.62, order
0.05: in a shear flow |grad psi| changes physically without redistancing) five-fold, and raises
the shape error 173-fold. Tables committed: `docs/gradient-controlled-level-set/gcls-level-set-
article/data/tables/methodGate2D_*` (the verdicts, one LaTeX table per arm, and the summary,
orders, vsBaseline and seam CSVs of every candidate; regenerate with `make_gate_summary.py
--docs` and `make_gate_tables.py`). The 3D gate does not run: no 2D pass.

### 11.16 The pre-print is written; its data are archived per software version (2026-09-27)

`docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex` is now a numerical
paper (21 pages): for each method the continuum model, the OpenFOAM discretisation and the gate
results, then a discussion that names, for each failure, the discrete mechanism it points to and
the experiment that can decide it (the source step moves the discrete zero set where F differs
between the two cells across the interface; the halo-limited extension loses order through the
face-to-cell reconstruction of the corrected flux; the sampler fits across the interface in the
coupled arms). Every number comes from ONE folder of secondary data,
`gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/` (116 files,
13 MB: the fixed gate, the pre-fix gate, the laptop runs; `README.md`, `MANIFEST.csv` with
SHA-256). `figures/make_archive.py` writes the archive from the raw output;
`figures/make_result_figures.py` writes the five result figures and eight tables from it; the gate
tables come from `make_gate_tables.py` on the archived summaries. A new pre-print version gets a
new folder; an older folder is never overwritten. The raw laptop runs are preserved in the
git-ignored `runs/gcls-laptop-20260927` (1.5 GB, 539 files), with a copy on Lichtenberg in
`/work/scratch/tm83tomy/leia/runs/`. `.gitignore` now exempts `docs/**/data/figures/*.pdf` and
`docs/**/data/archive/**/*.csv` (curated results).

### 11.17 Gate blind spots found on 2026-09-28 (OPEN repairs)

Two blind spots of the method gates, found while the record was corrected on 2026-09-28. Both
repairs are OPEN. The measurements are in the technical report
`docs/gradient-controlled-level-set/gcls-technical-report/`; the gates are described in
`docs/knowledge-base/concepts/method-gates.md`.

- (a) **The exact-1D closed form is vacuous for every candidate.** `workflow/scripts/make_gate_summary.py`
  lines 190-198 set `q_exact` to `None` for every candidate that is not plain (any of
  `VELOCITY_EXTENSION`, `SL_SOURCE` or `SDPLS_SOURCE` set), so `qError` is `None` and the criterion
  "target" is vacuous for the candidates. Repair: integrate the closed form of each candidate per
  band cell, dq/dt = q (F - alpha K(d/R)) and dd/dt = alpha d (1 - c^2), and compare `qBandMean`
  with it.
- (b) **The band-gradient metric is blind to the cell-scale mode.** `gradPsiMetric leastSquares` is
  centred, so the mode (-1)^(i+j) d, which destroyed the stationary droplet under HL1q and HL1z,
  does not change it. Repair: an L2 of the second difference of psi, or the min and max of a
  one-sided q, and a curvature column in the kinematic arms.

### 11.18 The knowledge base, the technical report and the paper split (2026-09-28/29)

The user's request of 2026-09-28: separate the semi-Lagrangian topic from the source-term topic,
build a cross-linked knowledge base in `docs/` from the slides and pre-prints that records the
reasoning and the decisions, and assess the first gradient-control campaign critically in a
technical report. Plan: `docs/plan-knowledge-base-2026-09-28.md` (approved
2026-09-28). Commits d1e3414 to the final commit of 2026-09-29 on `feature/gradient-controlled-level-set`.

1. **The knowledge base** `docs/knowledge-base/` (an Obsidian vault, built by Quartz v5):
   131 notes (7 hubs, 17 models, 55 concepts, 16 decisions, 13 retractions, 8 cases, 12 studies,
   3 sessions), the decision log and the retraction log. Start at `index.md`; the living handover
   is `sessions/current.md`; the handover to the SL session is `sessions/sl-session-handover.md`.
   Tools in `.quartz/`: `check_kb.py` (frontmatter, wikilinks, log coverage; the CI gate),
   `check_links.py` (every pinned GitHub link against the history), `build_graph.py` (the data of
   the 3D graph page `graph3d/graph.htm`), `missing_notes.py`, `build.sh` (`make kb`,
   `make kb-serve`, `make kb-graph`). The workflow `.github/workflows/knowledge-base.yml`
   replaces the dead `docs.yml`.
2. **Publishing is blocked by one setting (RESOLVED 2026-09-29: `main` carries the vault since the merge of 11.20, and its run published the site).** The `github-pages`
   environment allows deployments from `main` only (GitHub API, 2026-09-29). The workflow builds
   and checks the site on this branch; to publish from it, add the branch under Settings >
   Environments > github-pages > Deployment branches, or merge into `main`.
3. **The technical report** `docs/gradient-controlled-level-set/gcls-technical-report/`
   (`make report-gcls`, 11 pages, summary first, every claim marked MEASURED, DERIVED or
   HYPOTHESIS). CORRECTED, the pre-print's reading of the failures: the linear laws read no
   velocity, so the coupled loop cannot be the primary mechanism for HL1q/HL1z; the source step
   with a centred gradient of q is linearly unstable at a rate of order mu (DERIVED); the
   measured growth rates on the stationary droplet at N = 200 are about 290 1/s (curvature
   error) and 320 to 350 1/s (spurious current) over 10 to 20 ms against mu = 270 1/s. I was
   wrong to quote 245 and 276 1/s from the first reading of the histories; the archive gives the
   numbers above. The halo-limited extension at R = h relocates the normal strain into a shell
   at d = 1 to 2h (K = 1.27 at 1.5R); the exact-1D band means 0.46 to 0.59 measure it. The
   crossing shift is O(h) and cannot explain the O(1) shape errors. The experiment ladder E0.1
   to E3 is pre-registered in the report and in `concepts/gradient-control-next-experiments`.
4. **The paper split.** The four SL-baseline findings of the gcls pre-print moved into the SL
   article (the subsection "Consistency under domain decomposition", the subsection
   "Mass-momentum-consistent flux (rhoLENT)" with the Eulerian port, the paragraphs "The late
   instability: the outlet and the interior growth", the oscillating drift in the Limitations),
   with their tables and figures (`git mv`); `make_result_figures.py` writes them into the SL
   article's `data/`. The gcls abstract, contributions, discussion and conclusions follow the
   report. Both articles build without undefined references (52 and 19 pages).
5. **CLAUDE.md and AGENTS.md** (byte-identical): "The knowledge base is the point of
   reference" and "How we work: supervisor and expert developer".
6. **Record corrections** (the 13 audit flags of the plan): this file's date, SUPERSEDED
   markers at the filtered 2026-08-18/19 blocks, the VOID list of section 0 (the morning studies
   of 2026-09-02), 11.17; METHOD.md CORRECTED notes (lines 35, 353, 377, 756);
   `cases/default.parameter` (the clip, the pre-fix residual, the VOID block); the gate configs
   (the VOID `rhoDdtGate2D` basis, the pre-fix residual); `HL0.yaml` (the read-out conflated the
   interface with the band); `baselineEulerian.yaml`; the banner default of the SL solver
   (`geo_lin`); `workflow/README.md` (eight method libraries); `docs/IMPROVEMENTS.md`, the nPSL
   note, the library-split plan.
7. **Link integrity.** The note writers read the tree of d1e3414 and pinned many links to
   8867581, where the same line numbers show other text. `check_links.py` moved 65 pins on
   decisive evidence; a reviewer decided the 159 ambiguous ones by hand (125 OK, 14 moved to the other commit, 19 re-anchored, one pointed at the wrong file).
8. **The writers read the record against itself** and found inconsistencies that nobody had
   corrected; they are listed, with the corrections still to make at the source, in
   `docs/knowledge-base/sessions/sl-session-handover.md` (items 1 to 15). The two most important:
   the production curvature `cellCentreInverse` was never scored on the varying-curvature ellipse
   gate (the acceptance criterion of section 7 is not demonstrated for it; scored on 2026-09-29,
   11.19: the orders pass; the gain equals the per-face inverse's, and the per-rung form of the
   limit fails for every second-order delivery), and a row of the SL
   article's viscous-term table comes from the frozen-muf run of 4267d7b (re-run or retract: an
   author decision). Corrected at the source today: CLAUDE.md and AGENTS.md (a second clone would
   have loaded the library ~200 commits ahead; no job ran with it, 9.5), CLUSTER.md (the banner
   lists one stamp line per loaded library, up to nine), the 40 mm volume change of 11.15 (15 %).
   Open and not changed, and the most consequential: `2Dtranslation` is a REVERSED translation (no
   `oscillation` entry, the default is on, `tau` = `endTime`), so the circle moves at most 0.16 and
   returns at T, and the ladder of METHOD 8.3.7 read its errors at T (DERIVED from
   `velocityModel.C` lines 49-50; author decision: a wrong setup or a reversed-flow gate); the
   published polyhedral orders 3.28 and 1.46 ran with `SL_CLIP true` (the production default
   diverges at step 198 on the coarsest polyhedral 3D shear rung). Also open:
   `foam_log_state.sh` tests launch failures before death signatures on a
   zero-step log, the Snakefile the other way round; `make_archive.py` does not refuse an existing
   archive folder ("never overwritten" is a convention); the Phase C bit-identity scripts are not
   committed; every method-gate arm renders a hex mesh, so the polyhedral rung exists only in the
   regression set.
9. **Process: the spend limit stopped every agent three times**, each time before the end of a
   5-hour window of the individual limit (the windows that ended 2026-09-28 16:10 and 21:10 and
   2026-09-29 02:10). Since the second stop
   the work is restart-safe: the manifest `.quartz/manifest.md` fixes every slug and owner,
   `missing_notes.py` lists what is missing, each writer saves a source digest
   (`.quartz/digests/`, git-ignored) before drafting and every note at once, and the tree is
   committed after each batch (the incremental commits of 2026-09-28 and 2026-09-29).

### 11.19 The one-way 2Dtranslation, the production curvature on the ellipse and ellipsoid (2026-09-29)

The author's requests of 2026-09-29: `cases/2Dtranslation` must translate the droplet left to
right only; the production curvature `cellCentreInverse` must be tested on the signed-distance
ellipse and ellipsoid, after a check of what exists; two findings go to the SL session (its
hand-over note); then the merges into `development` and `main`. All runs below: laptop,
`leia-gcls` binaries (v2512), branch `feature/gradient-controlled-level-set`, commits f94e2a4a
to the last commit of this section.

1. **`2Dtranslation` was a REVERSED flow (retraction; every earlier number of the case is
   VOID).** Commit cd97e6da (2026-09-01) made the case from `2Dvortex` and dropped the line
   `oscillation @!OSCILLATION!@;`. `velocityModel.C:49` defaults `oscillation` to TRUE and `:50`
   defaults `tau` to endTime, so U(t) = cos(pi t/T) (1 0 0): the circle came back to its start at
   T, and the error metrics compared T with the initial fields (no `psiEnd`). The recorded data
   show it: U = (0.034 0 0) at t = 0.2446 and (-0.99991 0 0) at t = 0.4979; the centroid moves
   0.2500, 0.4082, 0.2502. Void: `kinematicTranslation2D`, `coneBoundMesh2Dtranslation` and its
   four arms (METHOD 8.3.4), `advConv2Dtranslation` (8.3.7), the 2Dtranslation rows of the
   regression and bit-identity gates (8.3.8, 10.3, 11.7, the C5 gate of the halo-limited plan).
   The trees are renamed `*_VOID_reversedTranslation_20260929` (laptop and Lichtenberg), the three
   curated tables are in `data/tables/VOID_reversedTranslation_20260929/` with a README.
   I was wrong in the METHOD statements "the unbounded translation saturates at N = 256" and
   "the cone bound is better at the coarsest rung": both were read from the reversed flow.

2. **The fix.** `OSCILLATION off` in `cases/2Dtranslation.parameter`. A new token
   `END_REFERENCE none|translate` (inert default `none` in `cases/default.parameter`): with
   `translate`, `workflow/scripts/write_end_reference.py` writes `psiEnd`/`alphaEnd` (the circle
   moved by U T, made by the same `leiaSetFields` from the pristine `0.org` inputs) before
   `decomposePar`; the Snakefile runs it from the template's content, `advect_bound_arm.sh` and
   `Allrun.sh` call it too. The three kinematic solvers print which error reference they use
   (the silent fallback to the initial fields hid the reversal for four weeks); inert,
   bit-identical CSVs on `2Dvortex` N = 32 np 4 for the SL and the Eulerian solver.
   `aggregate.py` keeps `shapeError` for a one-way row with an end reference. psi is
   zeroGradient on all four patches (item 3).

3. **Fixed psi values on a patch make the SL update unstable (MEASURED; a finding about the SL
   scheme, for the SL session).** The approved plan put the exact psi on all four patches; that
   was my recommendation, and I was wrong. The SL fit reads every physical patch value as a
   stencil datum (`SL_STENCIL_BOUNDARY_FACES include`, the default). Serial probe of four
   boundary variants (`workflow/scripts/translation_bc_probe.sh` and `_scan.py`),
   `E_GEOM_ALPHA_REL` at T:

   | psi on left (inflow) / right (outflow) / top, bottom | N = 128, CFL 0.5 | N = 256, CFL 0.5 | N = 128, CFL 1 |
   |---|---|---|---|
   | exact / exact / exact | 1.72e-03 | **1.98e-02** (outflow edge 4.5e+05) | **1.79** |
   | exact / zeroGradient / exact | 1.61e-03 | 4.44e-04 | **1.79** |
   | exact / zeroGradient / zeroGradient | 1.61e-03 | 4.44e-04 | **1.79** |
   | zeroGradient everywhere (committed) | 1.61e-03 | 4.44e-04 | 2.15e-03 |

   (a) The exact value on the OUTFLOW patch: the outflow-edge error grows by a factor of about
   1.1 per step (N = 64 to 256, CFL 0.5), reaches 4.5e5 at T at N = 256 and makes 96 false
   zero-set cells. It corrupted the `none` N = 256 rung of the first re-run (`E_GEOM_ALPHA_REL`
   order -3.53); those trees are `*_VOID_outflowDirichlet_20260929`. (b) The exact value on the
   INFLOW patch at CFL 1, where the departure point of the first cell column lies outside the
   domain: 2246 false zero-set cells in the inflow region; those trees are
   `*_VOID_inflowDirichlet_20260929`. With zeroGradient everywhere both are stable: outflow edge
   1.7e-03 at N = 256 (CFL 0.5), 8.8e-04 at N = 128 (CFL 1). No production case gives psi a
   fixed value, so no production run is affected; CLAUDE.md has the corollary. Consequence for
   this case: the inflow region x < t keeps psi0(0, y) >= 0.10 (no zero crossing, the interface
   metrics are unaffected), and the GLOBAL gradient metrics measure that frozen region.

4. **The repair gate** (`config/translationRepairGate2D{serial,par4}.yaml`, pre-registered; four
   runs, the first three void: a non-pristine end reference, then the exact outflow and the exact
   inflow value). RE-RUN 3 passes every line on both indicator arms: U = (1 0 0) at
   every write; 108 steps (the reversed case took 86); centroid at T (0.74972, 0.50000);
   `E_VOL_ALPHA_REL` at t = 0 1.65e-14; `E_GEOM_ALPHA_REL` 2.000000 at t = 0 (disjoint circles)
   and 7.210e-03 at T; every patch face equals its boundary cell; np 4 against serial 2.7e-12
   column-scaled; outflow edge 5.92e-03 at t = 0.2486 and 6.90e-03 at T; no false zero set.

5. **The translation studies, re-run** (all 34 runs COMPLETED, no false zero-set cell at any
   write time). `advConv2Dtranslation` (np 8) is the new preserved baseline of the hex 2D
   regression rung:

   | N | `none` E_GEOM | order | `stencilBounds` | `lipschitzCone` | order | cone / none |
   |---|---|---|---|---|---|---|
   | 32 | 5.364e-02 | — | 5.591e-02 | 1.323e-01 | — | 2.5 |
   | 64 | 7.210e-03 | 2.90 | 7.210e-03 | 1.019e-01 | 0.38 | 14.1 |
   | 128 | 1.609e-03 | 2.16 | 1.614e-03 | 9.664e-02 | 0.08 | 60.1 |
   | 256 | 4.442e-04 | 1.86 | 4.472e-04 | 7.535e-02 | 0.36 | 169.6 |

   `E_VOL_ALPHA_REL` of `none`: 2.156e-02, 4.557e-03, 5.798e-04, 2.359e-05 (orders 2.24, 2.97,
   4.62); of the cone: 7.7e-02 to 9.1e-02, no convergence. `E_BOUND_ALPHA` = 0 everywhere. The
   one-way errors are about 20 times the reversed ones (7.21e-03 against 3.41e-04 at N = 64),
   because nothing cancels at T. The four cone-bound arms on one N = 64 mesh (np 4, digest
   fd13ce29ac8e): `none` 7.21e-03, cone unity 1.02e-01 (14.1x), cone stencil 6.82e-02 (9.5x),
   `stencilBounds` 7.21e-03. The cone arm depends on the psi boundary value 0.1 away from the
   interface: 3.03e-02 with the exact outflow value, 1.02e-01 with zeroGradient (probe, N = 64);
   the mechanism is not measured.

   `kinematicTranslation2D` (18 arms, np 4), `E_GEOM_ALPHA_REL` at T, N = 32 / 64 / 128:

   | CFL | E_GEOM_ALPHA_REL | orders | E_VOL_ALPHA_REL | orders |
   |---|---|---|---|---|
   | 0.25 | 6.24e-02 / 8.62e-03 / 1.75e-03 | 2.86, 2.30 | 2.94e-02 / 7.81e-03 / 1.06e-03 | 1.91, 2.88 |
   | 0.5 | 5.36e-02 / 7.21e-03 / 1.61e-03 | 2.90, 2.16 | 2.16e-02 / 4.56e-03 / 5.80e-04 | 2.24, 2.97 |
   | 1 | 4.47e-02 / 2.42e-02 / 2.15e-03 | 0.88, 3.49 | 2.07e-03 / 1.79e-03 / 6.83e-06 | 0.21, 8.04 |

   The pre-registered read-out, on the geometric alpha error (never L_inf): second order at fixed
   CFL PASSES at CFL 0.25 and 0.5; at CFL 1 the orders are irregular (OPEN). The mechanism
   prediction (the error collapses at CFL 1) is FALSIFIED again: at CFL 1 the shape error is 3.4
   times worse at N = 64 and 1.34 times worse at N = 128 than at CFL 0.5; only the volume error
   falls, 85 times at N = 128, on one rung. The null control passes: projectedFlux and
   cellCentred agree to 5.0e-12 column-scaled in every column of every row.

6. **`cellCentreInverse` on the signed-distance ellipse and ellipsoid (MEASURED; the first
   scoring of the shipped delivery on varying curvature).** What existed: the ellipse gate of
   2026-08-12 to 08-15 (per-face inverse 1.98, cut-cell 1.02, cell mean 1.03, arithmetic 0.97)
   ran BEFORE `cellCentreInverse` existed (9bb4f9d5, 2026-08-18); the sphere ladder h^1.95
   against h^1.02 scored the per-face inverse; the ellipsoid gate ran on 2026-08-26 and was never
   curated; the static gates scored cCI only as a remainder term. What the author remembered is
   most likely those gates. Added (e089939f): the rows `solverCellCentreInverse` and
   `solverCellCentreInverseNoK` in `leiaTestMeanCurvature` and the models `cellCentreInverse` and
   `cellCentreInverseNoK` in `leiaTestCurvatureNoiseGain`, both calling the shipped function and
   the arithmetic face interpolation. Inert: every old row byte-identical on the circle (N = 128),
   the sphere (N = 32) and the ellipse (N = 128). The two EXISTING studies re-ran; no new study.
   Active-face L2 error [1/m]:

   | gate | N | cCI with K | cCI without K | per-face inverse |
   |---|---|---|---|---|
   | ellipse, signed distance | 64 / 128 / 256 / 512 | 11.2 / 3.28 / 0.831 / 0.208; orders 1.77, 1.98, 2.00 | identical | 14.2 / 4.32 / 1.12 / 0.279 |
   | ellipsoid, signed distance | 32 / 50 / 80 / 128 | 0.698 / 0.166 / 0.0571 / 0.0231; fit 50-128: 2.10 | 1.82 / 1.30 / 0.825 / 0.508; fit 1.00 | 0.643 / 0.193 / 0.0689 / 0.0286; fit 2.03 |
   | ellipsoid, implicit psi | 32 / 50 / 80 / 128 | 1.08 / 0.727 / 0.475 / 0.309; fit 0.91 | 1.84 / 1.48 / 1.10 / 0.685 | 1.19 / 0.728 / 0.468 / 0.308 |

   The pre-registered read-out (the two config headers): the orders PASS on both signed-distance
   gates; the K-off control is identical in 2D and first order in 3D (the Gaussian term carries
   the order); not falsified. The gain: G h^2 (linear regime) is 0.618 / 0.651 / 0.673 / 0.647
   at N = 64 to 512. The criterion of section 7 as PCS 12 applied it (the finest rung and the
   fitted order; the per-face inverse "passes" with 0.647 and 1.98) is MET with the same numbers.
   The stricter per-rung form that I pre-registered FAILS at N = 128 (0.6514) and 256 (0.673);
   the per-face inverse fails it at N = 256 (0.673; 0.6499 at N = 128 is below the limit by
   1e-4), the arithmetic delivery too (0.670 at N = 256), and only
   the first-order cell-mean (0.416) and symmetric face-mean (0.467) stay below 0.65. Correction
   of the same day: I first wrote "every delivery fails it at N = 256"; the knowledge-base
   writer found the two first-order exceptions. Ellipsoid, signed distance, N = 128: cCI
   0.838, per-face 0.836, arithmetic 0.834. The implicit psi is first order for every delivery
   (its foliation is not parallel); its gain per unit psi (0.060) is not comparable with the
   signed distance. Curated by hand (the finalize rule runs the same scripts):
   `face_curvature_ladder_ellipse.csv`, `face_curvature_orders_ellipse.*`, `curvature_gain_ellipse.*`,
   the `_ellipsoid3d_<surface>` tables and figures, `face_curvature_orders_ellipsoid3d_foliation.*`.
   `make_face_curvature_fig.py` and `make_curvature_gain_table.py` now write one set per
   `PSI_SURFACE` (they mixed the two psi of the ellipsoid gate in one fit), label the production
   rows, name the right geometry in the title and write the per-rung ladder. METHOD 4.1, 4.3, 8.1.

7. **Traps, for any session.**
   - `leiaSetFields` is not idempotent on a non-pristine `0/`: a second call on its own output
     wrote 0/1 indicator values at the new circle (`alphaEnd` off by 0.46 in a cell, 1.6 %
     volume). `write_end_reference.py` rebuilds `0/` from `0.org` for its call.
   - A verification study through the full workflow runs the finalize rule, which overwrote
     curated figures (`alpha_field_2Dvortex_hex.png`, three velocity-extension figures; restored
     with `git checkout`). Run verification studies with `--until solve`.
   - The per-value relative test of `compare_metrics_csv.py` fails on round-off columns
     (`E_BOUND_ALPHA` <= 1.6e-14, the t = 0 volume rows); read the column-scaled difference.

8. **Open.** (a) The SL treatment of Dirichlet psi data and of a departure point outside the
   domain (item 3; the SL session). (b) The CFL-1 orders of `kinematicTranslation2D` (0.88,
   3.49). (c) The cone bound's dependence on a boundary value 0.1 away (the bound is falsified
   anyway). (d) `3Dtranslation` keeps its own `init_End` with a hardcoded centre and symlinks
   that `materialize.py` breaks; no study runs it. (e) The circle and sphere face-curvature
   tables have no cCI rows yet: their studies were not re-run (the sphere inertness rung, N = 32:
   K-on 0.079, K-off 1.87).

### 11.20 The merges of 2026-09-29: development, main, the pruned branches

The author's request of 2026-09-29: merge everything into `development`, then into `main`, with
`feature/eulerian-rholent`, and prune the obsolete feature branches (plan steps E to G).

1. **`development` = 4f12ae71**, the merge of `feature/gradient-controlled-level-set` (537a6011,
   which contains `feature/eulerian-rholent`) onto the OpenFOAM-v2606 move (185e356a). The two
   predicted conflicts are resolved: STATUS.md keeps 10.6, then section 11; METHOD.md line 35
   keeps both notes (the eight libraries, v2606). The merge adds to the feature exactly
   development's own change: the same 375 files and the same line counts. The next commit,
   1d568a87, puts v2606 in the three places that the feature added with v2512 (the banner of the
   `oscillatingDroplet3D` fvOptions template, the launch recipe of the halo-limited plan 5.6, the
   default of `translation_bc_probe.sh`), makes the clone recipe of CLUSTER.md check out
   `development`, pins the two branch links of `concepts/method-gates` to a commit, and rewrites
   `sessions/current`.
2. **`main` = d89e3f56**: `git --no-replace-objects merge -s ours --no-ff origin/main` on 1d568a87,
   pushed as a fast-forward, no force. Its tree is development's (`git diff --quiet HEAD^1 HEAD`),
   and the old `main` 8f7b89fe (the history before the filter-repo rewrite of 2026-07-09) is its
   second parent. `development` does not take this commit. The local `refs/replace/*` map 8f7b89fe
   to its rewritten twin, so every command of this step ran with `--no-replace-objects`; with the
   replace refs active, git sees the old `main` as an ancestor of `development`. Later merges from
   `development` into `main` are normal merges with a recent base; this section reaches `main` by
   the first of them.
3. **GitHub Actions on `main`.** Knowledge base: build and deploy succeeded (run 36559552822, 11:27
   UTC); https://leia-openfoam.github.io/leia/ now serves the vault instead of the old Doxygen
   pages, so 11.18 item 2 is resolved. Build Tests (`.github/workflows/build.yml`, Docker image
   `tmaric/openfoam-v2206_ubuntu-focal`) FAILED in its compile step (run 36559552814): the code
   needs v2606, and the README badge is red. OPEN, author decision: a v2606 image or the
   retirement of the workflow.
4. **Branches.** Deleted on GitHub: `feature/eulerian-rholent` and `feature/velocity-extension`;
   both are contained in `main` (`git --no-replace-objects merge-base --is-ancestor`). KEPT:
   `feature/gradient-controlled-level-set`, because the Lichtenberg clones
   `/work/scratch/tm83tomy/leia` (40c39d7) and `leia-gcls` (f7067d9) track it and cannot move to
   `development` until OpenFOAM-v2606 is built on Lichtenberg (10.6); `leia-curvature` is on
   `development` (c431677). After that switch, delete it with
   `git push origin --delete feature/gradient-controlled-level-set`. The seven old branches of
   other authors stay (author decision, 2026-09-29). No tag was pushed. The laptop keeps its local
   branches; `leia-gcls-base` has `feature/eulerian-rholent` checked out, with no upstream.
5. **Not touched:** the laptop clone `~/OpenFOAM/repos/leia` of another session (`development` at
   c4316771, 101 commits behind, 51 uncommitted entries) and the Lichtenberg clones. The pull of
   the laptop clone, for its owner:
   1. Move aside the 8 untracked CSVs in
      `docs/semi-lagrangian-level-set/sl-level-set-presentation/data/tables/` that are
      byte-identical to development's (`npslConv2Dvortex`, `nslConv2Dvortex`,
      `nslConv3Ddeformation`, `nslKinematicGate`, `sl_fit_pivot_census`, `veSmoke`,
      `veSmokemeshWaveExt`, `veSmokesteadyUpwind`, each `*_errors.csv` or `.csv`).
   2. `git stash`.
   3. `git pull --ff-only`.
   4. `git stash pop`, then resolve `.gitignore`, `popinetTranslating.csv` and
      `stationaryDroplet3DbitIdentity_errors.csv`, which development also changed.
6. **To read this file on Lichtenberg** without a change to any working tree:

   ```bash
   cd /work/scratch/tm83tomy/leia-curvature && git fetch -q origin && git show origin/development:STATUS.md | less
   ```

### 11.21 CI on OpenFOAM-v2606; Lichtenberg on OpenFOAM-v2606 (2026-09-30)

The author's requests of 2026-09-30: a new Docker image, because CI must work; OpenFOAM-v2606 on
Lichtenberg, built in `$HOME/OpenFOAM` if needed; the input for the SL session delivered to it.

1. **The CI image.** `ghcr.io/leia-openfoam/openfoam-v2606_ubuntu-noble` (`:latest` and `:<sha>`),
   built and pushed by `.github/workflows/ci-image.yml` from `.github/docker/openfoam-v2606/Dockerfile`:
   Ubuntu 24.04, the openfoam.com package `openfoam2606-dev` (its dependencies bring the compiler and
   the OpenMPI development files), `openmpi-bin`, git, python3, bc, procps, make. The package is
   private to the organisation and linked to the repository; the build job pulls it with its token.
   Why our own image: the laptop has no Docker, and OpenCFD's newest development image is 2512
   (2606 exists only as the runtime image `opencfd/openfoam-run`, without compilers). The first image
   run pushed the image and failed its check on `foamVersion`, which the packages do not install
   (exit 127); the check now prints `WM_PROJECT_VERSION` and the second run passed.
2. **Build Tests** (`.github/workflows/build.yml`), on pushes to `main`, `development` and
   `releases/**`, on pull requests to `main` and `development`, and by hand. Its one step is
   `.github/scripts/ci-build-and-smoke.sh`, which runs the same way on a workstation
   (`.github/scripts/ci-build-and-smoke.sh $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc`):
   1. `./Allwmake`;
   2. `etc/leia-check-build.py` (new): every `EXE`/`LIB` target of every `src/**/Make/files` and
      `applications/**/Make/files` exists, because `wmake all` does not always fail on one broken
      application;
   3. `cases/2Dtranslation/Allrun.sh` (serial, N = 128, CFL 0.5), classified by `foam_log_state.sh`;
   4. its numbers: `E_GEOM_ALPHA_REL` = 2 at t = 0 (the exact end reference is in use), below 3e-3 at
      T; `E_VOL_ALPHA_REL` below 1e-12 at t = 0 and 1e-3 at T; `E_BOUND_ALPHA` = 0.
   MEASURED, laptop (OpenFOAM-v2606 source build, fbc5916b): 26 executables and 14 libraries, 0
   missing; 168 steps; `E_GEOM_ALPHA_REL` 2.000000 -> 1.6090e-03 and `E_VOL_ALPHA_REL` 2.27e-14 ->
   5.7984e-04, the `kinematicTranslation2D` values at N = 128, CFL 0.5 (11.19). GitHub: run
   36706765044 on `development` f010ed1a PASSED in 8 minutes (the build step 7.5 minutes). A failure
   writes an error annotation with the tail of its log, and a pass a notice with the numbers: the
   public check-runs API returns annotations without login, while the job logs need admin rights.
3. **Two traps, fixed in the script.** Sourcing OpenFOAM's `etc/bashrc` passes the script's own
   arguments to it, and it sources every readable file among them: a path to the bashrc in `$1`
   recursed until bash segfaulted (exit 139); the script runs `set --` first.
   `cases/2Dtranslation/Allrun.sh` has mode 644 in git (exit 126); the script runs `bash ./Allrun.sh`.
4. **Knowledge base.** Its deploy job runs on `main` only; on `development` the build job is the
   check, and the run ends green (it ended red at every development push before).
5. **Lichtenberg on OpenFOAM-v2606** (`$HOME/OpenFOAM/OpenFOAM-v2606`, `ThirdParty-v2606`,
   `cfmesh-v2606`; the recipe is in CLUSTER.md; every job id is in `$HOME/OpenFOAM/.my_jobs` and
   `/work/scratch/tm83tomy/leia-dev/.my_jobs`).
   1. Pass 1, job 55179573 (48 cores, 40:49, exit code 0): OpenFOAM plus ThirdParty, then cfMesh
      from the release's `plugins/cfmesh` (25 executables, `libmeshLibrary`, stamp `v2606`; `pMesh`,
      `cartesianMesh`, `tetMesh` present). But only 151 applications and 113 libraries: `icoFoam`,
      `interFoam`, `simpleFoam`, `potentialFoam`, `pimpleFoam`, `pisoFoam`, `scalarTransportFoam` and
      `setFields` were missing, and `foamInstallationTest` reported 1 critical error. Cause: the job
      sourced `etc/bashrc` before ThirdParty built FFTW, so FFTW's `lib/` was not on the loader path;
      the link of `noise` against `librandomProcesses.so` (which needs `libfftw3.so.3`) failed and
      make stopped before the standard solvers. The exit code 0 of the job was no evidence.
   2. The two leia jobs behind it (55179650 running for 3 minutes, 55179651 not started) were
      cancelled by id, and the partial leia build removed.
   3. Pass 2, job 55179978 (6:57): `etc/bashrc` sourced with FFTW present. 270 applications and 129
      libraries (v2512: 270 and 128); every standard solver present; `foamInstallationTest`
      "Critical systems ok"; `librandomProcesses.so` resolves `libfftw3.so.3` in ThirdParty-v2606.
   4. A new clone for the development line, `/work/scratch/tm83tomy/leia-dev` (development
      f010ed1; the other three clones are untouched), built by job 55179979 (32 cores, 1:16) with
      `.github/scripts/ci-build-and-smoke.sh`: PASS, 26 executables and 14 libraries, 0 missing,
      stamp `v2606`; 168 steps; `E_GEOM_ALPHA_REL` 2.000000 -> 1.6090e-03 and `E_VOL_ALPHA_REL`
      2.27e-14 -> 5.7984e-04, the laptop's values to every printed digit; `pMesh` and
      `cartesianMesh` resolve from `cfmesh-v2606`.
   5. The workflow end to end: `sbatch -J leia-dev-smoke --export=ALL,STUDY=translationRepairGate2Dpar4
      run-studies.sbatch` (driver 55179980, 8 child jobs 55180038 to 55180053 through
      `profiles/slurm`): rc 0; both arms COMPLETED, 108 steps, 4 ranks, the binary of the clone's own
      `platforms/`. Against the laptop's RE-RUN 3 (11.19 item 4: v2512, 4 ranks), every metric column
      is bit-identical except the four volume sums, which differ by at most 2.7e-12 column-scaled
      (1.1e-10 relative at one early row), the order of the parallel sum. End values on both:
      `E_GEOM_ALPHA_REL` 7.2095e-03, `E_VOL_ALPHA_REL` 4.5568e-03.
   So Lichtenberg runs the development line on OpenFOAM-v2606 through the committed profile. The
   study's two curated tables were removed from the clone.
6. **The SL session.** The full input (the tasks T1 to T3 and the results of 2026-09-29) is queued
   for "SLLS: Polyhedral stationary and translating" (a Remote Control session on another machine,
   offline on 2026-09-29 and 2026-09-30); delivery waits until that machine reconnects.
7. **Open.** The clones `leia` (the SDPLS session), `leia-gcls` and `leia-curvature` on Lichtenberg
   still run v2512 binaries; their owners move them (CLUSTER.md, "Where leia lives"); after that,
   `feature/gradient-controlled-level-set` can be deleted on GitHub. TwoPhaseFlow is not built on
   Lichtenberg for v2606 (the `interFlow` arm). A parallel smoke run and the seam check are not in
   the CI yet. The actions `checkout@v4`, `upload-artifact@v4` and the Docker actions target
   Node.js 20; GitHub runs them on Node.js 24 with a warning.

### 11.22 `leiaRedistancedLevelSetFoam` retired; the semi-Lagrangian pair is not equivalent yet (2026-09-30)

The author's question of 2026-09-30: `leiaLevelSetFoam` can run with no source term, no velocity
extension and a redistancer, so why a second solver? Then the request: retire it, and run the same
test on the semi-Lagrangian pair. Laptop, OpenFOAM-v2606 source build, a worktree of `development`
06e15ecf, 4 ranks everywhere.

1. **Why it was redundant (DERIVED, then MEASURED).** `leiaRedistancedLevelSetFoam` assembled
   `ddt(psi) + div(phi psi) - psi div(phi) = S_SDPLS(psi)` inline, with the deferred-correction
   loop, then called the gated redistancer. `leiaLevelSetFoam` hands the same step to its
   `eulerian` advection model, which builds the same equation with the extension flux; with
   `velocityExtension none` that flux is the prescribed one. Then it calls the same redistancer and
   a volume correction whose default does nothing. The committed case dictionaries already select
   these models, so only the executable changes. `leiaLevelSetFoam` dates from 2021 (a0d38c2e);
   the redistanced solver came in 2a53364b (2026-07-30), whose message says only that it "exercises
   the redistancer line on its own". I was wrong in my answer of the same day: I wrote that both
   solvers came in that one commit. That commit added the redistanced solver and changed
   `leiaLevelSetFoam`.
2. **The retirement gate (pre-registered: bit-identical in every metric column at every step, and
   the final fields byte-identical).** Three scratch studies through the workflow up to the solve
   rule, before (the old solver, the unchanged tree) and after (`leiaLevelSetFoam`, the changed tree
   with the old solver deleted and its binary removed): 2Dvortex N = 32, T = 0.5 with
   `noRedistancing`, `PDE`, `planeFootWave`, `anchoredEikonal`; 2Dvortex N = 32 with the SDPLS
   sources `R` and `beta`, each with `strictNegativeSpLinearImplicit` and `explicit`; 3Dshear hex
   N = 24, T = 0.75 with `PDE` and `planeFootWave`. PASS: in all 10 arms the solver CSV (14 columns)
   and `gradPsiError.csv` are identical in every column at every step (74 steps in 2D, 85 in 3D;
   `compare_metrics_csv.py --tol 0`, timing columns skipped), and the final psi and alpha files are
   byte-identical on all four ranks. Control: in the same gate, `PDE` and `planeFootWave` differ
   from each other in 7 columns, so the identity is not trivial. The aggregated errors table has the
   same values, with the solver column now `leiaLevelSetFoam` and the redistancer and trigger filled;
   `make_grl_fig.py` finds the new rows. An earlier test of the same day (the v2512 build of
   `leia-gcls`, `bulkVortexGRL` rendered at N = 32, `PDE` and `planeFootWave`) had given the same
   result.
3. **What changed.** The solver folder is deleted. The nine study configs that ran it select
   `leiaLevelSetFoam` with an explicit `theme: geometrically-redistanced-levelset` (the theme used to
   follow from the solver name). `aggregate.py` fills `redistancer` and `redistTrigger` for
   `leiaLevelSetFoam` rows too (the only reader is `make_grl_fig.py`); `make_grl_fig.py` accepts both
   solver names, because the curated tables keep `leiaRedistancedLevelSetFoam` in the rows that it
   produced; `paths.py` keeps its theme entry for those archived studies. Comments in the case
   templates, `cases/default.parameter` (token set unchanged, checked), three library files, the
   README, the Makefile, `workflow/README.md` and the two redistancing decks name `leiaLevelSetFoam`.
   The curated data tables, the plan documents and the older STATUS sections keep the old name: they
   record what ran. 25 executables and 14 libraries now; `etc/leia-check-build.py` finds none
   missing.
4. **Two traps, fixed.** (a) `git rm` leaves the untracked object folder `Make/<WM_OPTIONS>/` of the
   deleted solver in every clone that built it, and `wmake all` then stops on the missing
   `Make/files`: `Allwmake` rc 2 while every target exists. `Allwmake` now removes that folder (only
   when no `Make/files` is left) and the stale binary. (b) `make_grl_fig.py` crashed in its order
   table (division by zero when two rows of one group share h) on the committed tables alone, before
   this change, so the curation of every redistancing study failed there. A guard skips the order
   for equal h.
5. **The semi-Lagrangian pair: NOT equivalent (MEASURED).** `leiaSemiLagrangeLevelSetFoam` against
   `leiaLevelSetFoam` with `levelSet.advection.type semiLagrangian`, identical starting copies, the
   rest of the case dictionary unchanged: 2Dvortex N = 32, T = 0.5 (74 steps); 2Dtranslation N = 64,
   T = 0.5 (108 steps); 3Dshear hex N = 24, T = 0.75 (85 steps); each with `SL_TRACE_VELOCITY`
   `cellCentred` and `projectedFlux`, and the phase indicators `geometric` and `detrixheAslam`
   (12 arms, 24 runs, all COMPLETED).
   1. `cellCentred`: psi is byte-identical at T in all 6 arms, and every gradient metric is identical
      at every step. The transport is the same.
   2. `projectedFlux`, the production default: psi differs in every cell of 2Dvortex and 3Dshear (by
      up to 2.2e-05 and 6.9e-05 at T), and in one cell by 1e-15 in 2Dtranslation (a uniform
      velocity, where the two trace velocities agree to round-off). The `semiLagrangian` advection
      model has no trace-velocity option: it always traces with the cell velocity, and it has no
      flux-extension trace either.
   3. The alpha metrics differ even with `cellCentred`: `E_VOL_ALPHA` at 21 to 65 of the steps, the
      relative volume error by up to 5.6e-04 (2Dvortex), 7.5e-04 (2Dtranslation) and 7.7e-04
      (3Dshear); `E_GEOM_ALPHA` at 21 to 63 steps, the same order. That is not small against the
      error itself: at step 53 of 2Dvortex the SL solver has a relative volume error of 6.5e-04
      where `leiaLevelSetFoam` has 9.2e-05 (7.1x); at T of 3Dshear 2.14e-03 against 1.86e-03
      (+15 %); at T of the two 2D cases the two agree. The final alpha differs in 2 of 1024, 10 of
      4096 and 3 of 13824 cells, by up to 0.024, 0.082 and 0.058. Cause: both phase indicators read
      the `NarrowBand` field, and the SL solver computes alpha BEFORE `narrowBand->calc()` (with the
      band of psi^n), `leiaLevelSetFoam` after it (the band of psi^{n+1}). The order dates from the
      first SL versions (July 2026, ec160a46 and 4044f472) and carries no comment. The two-phase SL
      solver has had the consistent order since it was written (2a53364b; `slAlphaEqn.H` L167-L172:
      "previously it lagged one step, giving cells that had just entered the band a sign-based 0/1
      alpha for one step"); the kinematic SL solver kept the lag. OPEN, author decision: the
      per-step volume and shape errors of every SL kinematic study carry this stale-band effect,
      measured here only at N = 24 to 64, with an unknown order in h; a fix changes their metrics
      CSVs and needs the advection regression set at three resolutions.
   4. `L_INF_E_PSI` differs in every arm, also where psi is identical: the Eulerian
      `advectionErrors.H` of `leiaLevelSetFoam` still has the mis-parenthesised sign test
      `sign(psi == sign(psi0))`, which the SL copy fixed on 2026-08-27 (84906ee4). L_inf is never
      reported, so no reported number is affected. OPEN: port the fix.
   So `leiaSemiLagrangeLevelSetFoam` stays. To retire it: port the trace-velocity options into
   `semiLagrangianAdvection`, decide the phase-indicator order, port the sign-test fix, then repeat
   this test.

### 11.23 The kinematic SL solver computes the narrow band before alpha; it stays (2026-10-01)

The author's decisions of 2026-10-01: (a) fix the stale band of 11.22 item 5.3; (b) the kinematic
SL solver stays: "we need one kinematic solver for each flow solver if we use reversed prescribed
analytic velocity". A kinematic solver runs the interface step of its flow solver with a prescribed
velocity, so a reversed-flow test verifies the flow solver only if the two run the same interface
code. I withdraw the proposal of 11.22 to retire `leiaSemiLagrangeLevelSetFoam`: the
`semiLagrangian` model of `leiaLevelSetFoam` does not run the interface step of the SL two-phase
solver (no `projectedFlux` trace). Laptop, OpenFOAM-v2606, np 4, a worktree of `development`
1962df1f.

1. **The change (6f63418a).** `leiaSemiLagrangeLevelSetFoam` calls `narrowBand->calc()` after the
   advection in every step, then computes alpha on the report steps with the band of psi^{n+1}:
   the order of the SL two-phase solver (`slAlphaEqn.H` L167-L172). Before, it computed alpha with
   the band of psi^n and refreshed the band only on the report steps. Both phase indicators give a
   sign-based 0/1 alpha outside the band (`detrixheAslamPhaseIndicator.C` L151-L182). The Eulerian
   pair already has the consistent order (`alphaEqn.H` L52-L92 of `leiaLevelSetTwoPhaseFoam`,
   `leiaLevelSetFoam.C` L174-L181).
2. **Where the band reaches (code reading).** The default trajectory (`trajectoryVelocity input`)
   does not read the band; only `normalClosestPoint` does (`pointValueScheme.C` L205-L218), and only
   a two-phase study runs it (`config/transISTTrajectoryGate.yaml`). Every study that reports at
   write times only (`REPORT_WRITE_ONLY true`, the `benchVortex*` method comparison) runs
   `leiaLevelSetFoam`, which refreshes the band in every step. So the defect reached only the alpha
   metrics of `leiaSemiLagrangeLevelSetFoam`, by one step.
3. **The gate.** Each rung was pre-registered in a scratch config header before it ran
   (`slBandGate{2Dvortex,3DshearHex,3DshearPoly}.yaml` first, then
   `slBandGate{3DshearPolyClip,2DvortexPublished}.yaml`, each written after the earlier results and
   before its own run). One prepared and decomposed case per rung (the workflow up to `decompose`),
   two identical copies: the pre-fix solver (the source of 1962df1f, compiled under another name
   against the same libraries) and the fixed solver. The rungs:
   - 2Dvortex hex N = 32, 64, 128, 256, T = 2, CFL 0.5, the `none` arm of `config/advConv2Dvortex.yaml`;
   - 3Dshear hex N = 32, 48, 64, T = 3, CFL 0.3 (`advConv3DshearHex`);
   - 3Dshear poly MAX_CELL_SIZE 0.04, 0.03, 0.0225 (49,910 / 105,193 / 245,554 cells), CFL 0.3,
     `SL_CLIP false` (`advConv3DshearPoly`, the production default);
   - 3Dshear poly in the published configuration (`config/uncachedConv3DshearPoly.yaml`:
     `SL_CLIP true`, `SL_STENCIL face`, CFL 0.5), MAX_CELL_SIZE 0.03125, 0.0235, 0.0176
     (105,158 / 216,747 / 473,103 cells);
   - the published 2D table itself (`config/uncachedConv2Dvortex.yaml`): N = 32 to 256 in steps of
     sqrt(2), CFL 0.5 and 1.0.
   All 54 runs COMPLETED; pre and post have the same step count in every pair. The 3D rungs are
   under-resolved for orders (R = 0.15, R/h = 4.8 to 11.7; the method gates ask R/h >= 10 at the
   first rung): their orders bound the change on these meshes, they are no convergence result.
4. **P1, the transport is unchanged: PASS for psi.** psi is byte-identical at every written time on
   every rank in all 27 pairs, and every column of `gradPsiError.csv` is identical at every step.
   The read-out also listed `L_INF_E_PSI` as a psi column. That was wrong: the column zeroes the psi
   error in the cells where alpha is 0 or 1 (`advectionErrors.H` L21-L40), so it depends on alpha.
   It changes in 6 of the 13 pairs of the first four rungs. L_inf is never reported.
5. **P2, the fix acts: yes, except in the published polyhedral configuration.** The alpha columns
   differ at 58 of 154 to 845 of 884 steps (2Dvortex N = 32 to 256), 334 of 375 to 684 of 704
   (3Dshear hex) and 427 of 968 to 1192 of 1730 (3Dshear poly, production). In the published
   polyhedral configuration they are identical at every step (0 of 790, 1008, 1317); I have not
   found why. At T, the cells whose alpha differs lie within one cell of the exact interface on hex
   (2D N = 128: 5 cells at |psi0| = 0.61 to 0.69 h; 3D N = 32: 4 cells at 0.62 to 0.88 h), and 29 h
   away from it in the production polyhedral arm: a false zero set of the far-field defect (item 8).
6. **The size, per step, on the production ladders (CFL 0.5 in 2D, 0.3 in 3D).** The shape error
   moves little: the median relative change per step is at most 0.01 %, the maximum 1.7 % (2D) and
   0.19 % (3D hex). The volume error moves much more, and in
   2D the effect grows with N: the median relative change per step is 0.52 / 1.82 / 4.51 / 6.22 % at
   N = 32 / 64 / 128 / 256 (3D hex 0.27 / 0.19 / 0.10 %; poly production 0.00 %). Whether a read-out
   instant is hit depends on whether a cell enters the band in the last step before it.
   DERIVED: near T/2 the velocity factor cos(pi t/T) is near zero, so few cells enter the band
   there, and the T/2 read-out is protected; the T read-out is not.
7. **P3, the size at the read-out instants, and the decision.** Pre-registered threshold: an order
   change above 0.1, or a change of a finest-rung error above 5 %, means the curated table must be
   regenerated.
   - The production ladders: 2Dvortex E_GEOM at T -1.66 % at N = 128, unchanged elsewhere;
     E_VOL_ALPHA_REL at T/2 unchanged at every N, at T +4.80 % at N = 128 (2.233e-04 to 2.341e-04);
     consecutive-rung orders change by at most 0.068. 3Dshear hex: every change at most 0.35 %,
     orders by at most 0.011. Both poly configurations: unchanged to the printed digits.
   - **The published 2D table (`uncachedConv2Dvortex_errors.csv`): the threshold is MET.** The
     pre-fix runs reproduce the published numbers (fitted orders within 0.008, `volumeError` within
     0.8 %, `shapeError` within 1 % at N = 32). Fitted orders as in `make_convergence_table.py`:
     | CFL | | shapeOrder | volumeOrder | volumeHalfOrder |
     |---|---|---|---|---|
     | 0.5 | published | 2.843 | 3.303 | 1.766 |
     | 0.5 | pre / post | 2.839 / 2.841 | 3.299 / 3.294 | 1.763 / 1.763 |
     | 1.0 | published | 2.378 | 3.542 | 1.819 |
     | 1.0 | pre / post | 2.374 / **2.465** | 3.540 / **3.193** | 1.811 / 1.814 |
     At CFL 1.0 the shape error at T falls by 10.5, 14.8, 23.3, 11.0 and 15.5 % at N = 64, 90, 128,
     181 and 256 (the stale band inflated it), and the volume error at T moves by -10.5 % to +511 %
     (N = 128: 5.368e-05 to 3.282e-04). At CFL 0.5 the finest-rung errors are unchanged and the
     fitted orders move by at most 0.005. DERIVED: at CFL 1 the interface crosses up to one cell per
     step, so a cell that enters the band is more often cut by the interface.
   So the CFL 1.0 row of the SL article's 2D convergence table is RETRACTED as published (shape
   order 2.378, volume order 3.542); the measured values with the fixed solver are 2.465 and 3.193
   (laptop, v2606). The table is not regenerated yet: author decision (it belongs to the SL article
   of the SL session). 19 curated tables carry rows of this solver; 8 of them have arms at CFL 0.8
   or 1.0 (`uncachedConv2Dvortex`, `npslConv2Dvortex`, `nslConv2Dvortex`, `sdCompare2D`,
   `linearConv2Dvortex`, `linearConv2DvortexClip`, `linearConv3Dshear`, `kinematicTranslation2D`);
   none of them is re-measured.
8. **Two observations, not caused by the fix (pre and post agree).** (a) The production
   configuration fails on the polyhedral 3D shear rung: `E_VOL_ALPHA_REL` at T is 7.6, 11.6 and 18.7
   and grows with refinement. This agrees with the record (`value_bound_advection_3DshearPoly.csv`,
   2026-09-10: the `none` arm DIVERGED at step 186, 38.3) and with item T2 of the SL hand-over: the
   published polyhedral numbers need `SL_CLIP true`. (b) Some volume errors on these ladders are not
   monotone in h (2Dvortex at T/2: 5.21e-03 at N = 64, 9.49e-03 at N = 128; 3Dshear hex at T:
   2.77e-03 at N = 48, 3.31e-03 at N = 64). Recorded only.
9. **P4, the pair of 11.22 (`cellCentred`): now identical.** With the fix, the kinematic SL solver
   equals `leiaLevelSetFoam` `semiLagrangian` in every column except `L_INF_E_PSI` (the sign-test
   defect of the Eulerian writer, still open), and the final psi and alpha files are byte-identical:
   2Dvortex N = 32 and 3Dshear hex N = 24, two phase indicators each.
10. **CI.** `.github/scripts/ci-build-and-smoke.sh` on the fixed tree passes. The smoke numbers at T
    changed: `E_GEOM_ALPHA_REL` 1.6090e-03 to 1.5892e-03, `E_VOL_ALPHA_REL` 5.7984e-04 to
    5.6281e-04 (2Dtranslation N = 128, serial).
11. **Open.** (a) Regenerate the curated tables of item 7, the 2D ones on the laptop and the 3D ones
    on Lichtenberg, and correct the article and deck text that quotes them: author decision, and
    the SL article belongs to the SL session. (b) The `L_INF_E_PSI` sign test in `leiaLevelSetFoam`.
    (c) A shared library function for the interface step of each kinematic and flow solver pair,
    so that the two cannot drift again (a refactor with a bit-identity gate).

### 11.24 `L_INF_E_PSI`: the sign test fixed in the shared error writer; one writer for four solvers (2026-10-01)

The author's instruction of 2026-10-01: fix the `L_INF_E_PSI` defect of 11.22 item 5.4 first, then the
Helmholtz projection of the prescribed flux on polyhedral meshes (11.25). Laptop, OpenFOAM-v2606, np 4.

1. **The defect, wider than 11.22 said.** `applications/solvers/leiaLevelSetFoam/advectionErrors.H`
   tested `sign(psi == sign(psi0))`: the comparison is evaluated first and `sign(0) = 1`, so the
   guard was always true and the psi error was zeroed in every bulk cell, also where the sign had
   switched. The file serves THREE solvers: `leiaLevelSetFoam`, and both two-phase solvers
   (`leiaLevelSetTwoPhaseFoam/advectionErrors.H` is a git symlink to it, and
   `leiaSemiLagrangianLevelSetTwoPhaseFoam` includes that directory). 11.22 named only
   `leiaLevelSetFoam`; that was incomplete. The kinematic SL solver had its own copy, fixed on
   2026-08-27 (84906ee4); the two copies differed only in this line and in comments.
2. **The fix (60d3558e).** The test is `sign(psi) == sign(psi0)`, and the kinematic SL solver's copy
   is now a symlink to the same file: four solvers compile one writer, so the copies cannot drift
   again. The comments that cite line numbers of the file are updated (+8 lines).
3. **The gate**, pre-registered in the scratch config headers (`linfGate*.yaml`). Before = 44dff788,
   after = 60d3558e, both built; one prepared, decomposed case per arm, two identical copies, np 4.
   Arms: `leiaLevelSetFoam` 2Dvortex N = 32, T = 0.5 (`noRedistancing`, `PDE`) and 3Dshear N = 24,
   T = 0.75; `leiaSemiLagrangeLevelSetFoam` 2Dvortex N = 32 and 3Dshear N = 24;
   `leiaLevelSetTwoPhaseFoam` stationaryDroplet2D N = 32, 50 steps;
   `leiaSemiLagrangianLevelSetTwoPhaseFoam` translatingDroplet2D N = 32, 100 steps, and a 500-step
   arm added after the first result (its read-out written before its run). All COMPLETED.
   - P1 PASS: every column of every CSV is identical at every step except `L_INF_E_PSI`, and the
     final fields are byte-identical (the kinematic arms at T; the 500-step droplet at its write
     time 0.02 s, 72 files; the short two-phase arms write no time after 0).
   - P2: `L_INF_E_PSI` is never lower after the fix. It is higher at 14 and 17 of 75 steps (2D,
     the two redistancer arms), 28 of 86 (3D) and 145 of 501 (the 500-step droplet). The 100-step
     droplet did not change: it moved 2.9e-04 m, less than one cell (h = 3.1e-04 m), so no bulk
     cell switched sign. My prediction of a rise there was wrong for that horizon. The stationary
     droplet did not change either.
   - P3 PASS: the kinematic SL solver is identical in every column.
   - P4 PASS: the `cellCentred` pair of 11.22 (`leiaSemiLagrangeLevelSetFoam` against
     `leiaLevelSetFoam` `semiLagrangian`) is now identical in every column, `L_INF_E_PSI` included:
     2Dvortex N = 32 and 3Dshear N = 24, two phase indicators each.
4. **Reach.** No script reads `L_INF_E_PSI`, and L_inf is never reported: the raw CSVs of the three
   solvers carry the zeroed column, and no curated number changes. 11.23 item 11(b) is done.

### 11.25 The prescribed flux on polyhedral meshes: Helmholtz projection; two defects of `correctFlux` fixed (2026-10-01 to 02)

The author's instruction of 2026-10-01: on polyhedral meshes the face flux U(x_f) & S_f has a large
divergence, and the SL `projectedFlux` trace reconstructs its cell velocity from that flux, so the
flux must be projected to divergence-free with a Helmholtz projection, on polyhedral meshes only.
The same day the author stopped the large laptop runs (the RAM was full); the transport gate moved
to Lichtenberg (`/work/scratch/tm83tomy/leia-fluxproj`, branch `feature/flux-projection`).

1. **The divergence (G0, MEASURED).** A new log line at every start (`fluxDivergenceReport`): on the
   cfMesh 3D shear meshes the prescribed flux has max|div(phi)| = 2.3 to 4.0 1/s (the velocity
   gradient scale is about 2 pi), mean |div| 0.0016 to 0.019 1/s (falling about as h^1.3), and up to
   6 % of a cell's face flux uncancelled. The maximum does not fall with refinement: 3.97 / 3.98 /
   3.99 / 4.00 1/s at 105,158 to 4,832,366 cells. On uniform hex meshes it is round-off (1.8e-14 to
   2.7e-13 1/s): for this field the three face differences of a hex cell are 2, -1 and -1 times
   sin(2 pi x) sin(2 pi y) sin(2 pi z) sin(pi h), and they cancel exactly.
2. **Two defects of the existing `correctFlux`, fixed.**
   (a) The pure-Neumann problem was pinned with `setReference(0, 0)` on every rank: N ranks solved a
   problem with N pinned cells, so the projected flux kept a divergence at N - 1 cells and depended
   on the decomposition. Now the master rank alone pins its cell 0 (the serial choice in serial).
   (b) The projection cannot change a boundary flux (zeroGradient p), so a nonzero net boundary flux
   makes the Neumann problem incompatible, and the pinned cell absorbed all of it. The 3D shear top
   and bottom patches carry 0.136 m3/s each, but on cfMesh meshes their face-centre sums differ by
   -4.8e-05 to -2.8e-07 m3/s (on hex they cancel to 7e-49). After the projection the pinned cell kept
   |div| = 5 to 56 1/s, and the mean |div| equalled imbalance / volume to four digits. Now the
   outflow of the non-coupled patches is scaled first so that the net boundary flux vanishes, as
   `adjustPhi` does (factor 1.00002 to 1.00035 on the laptop meshes).
   With both fixes the projected flux has max|div| <= 2.8e-09 1/s on every mesh, laptop and cluster,
   and serial against np 4 agrees to 1.0e-15 m3/s (2.3e-12 of the largest face flux, the 0.03 mesh).
   The only study that ran the old code, `benchVortexSLimprovedPerturbed` (`-fluxCorrection` on a
   perturbed mesh, np 4), carried defect (a); it feeds no curated table.
3. **The change (9ccc9783, cfb1a58f).** `velocityModel.fluxProjection none | helmholtz` (default
   `none`): `setVolumetricFlux` reports the divergence and, with `helmholtz`, projects the flux once
   at the start; the cos(pi t/tau) scaling keeps it divergence-free. Token `FLUX_PROJECTION`: `none`
   in `cases/default.parameter`, `helmholtz` in the poly layer of 3Dshear, 3Ddeformation, 3Drotation
   and 3Dtranslation. Only `leiaSemiLagrangeLevelSetFoam` runs these cases on polyhedra (13
   configs); the two-phase cases and `leiaTestMeanCurvature` do not use the velocity model.
4. **G1, hex inert: PASS.** Fresh hex cases rendered with `fluxProjection none`, the old build
   (2b642398) against the new one, np 4: every CSV column and every field byte-identical (SL kinematic
   3Dshear N = 24 and 2Dvortex N = 32, `leiaLevelSetFoam` 3Dshear N = 24).
5. **G2, transport.** Pre-registered in a laptop file before the laptop runs and in the headers of
   `config/gates/fluxProjectionGate3DshearPoly{,Clip}.yaml` before the cluster runs; `none` against
   `helmholtz` on identical copies of each prepared mesh (`workflow/scripts/flux_projection_gate.sh`),
   trace `projectedFlux`, table by `workflow/scripts/flux_projection_gate_table.py`.
   - Laptop, before the stop (the finest rungs, R/h 9 to 12): published configuration, shape at T
     -0.3 / +0.2 %, volume at T +12.1 / -1.1 %, volume at T/2 and band gradient unchanged. The laptop
     prediction (a lower volume error at T at both rungs) was falsified there. Four smaller pairs were
     stopped: no result.
   - Lichtenberg, 14 runs, all COMPLETED at T = 3, np 16, wall/CPU = 1.0. Published configuration
     (`config/uncachedConv3DshearPoly.yaml`: `SL_CLIP true`, face stencil, CFL 0.5), 105,158 /
     347,073 / 1,272,989 / 4,832,366 cells (R/h 7 to 25):
     | quantity | none | helmholtz | change |
     |---|---|---|---|
     | shape error at T | 3.02e-02 / 1.08e-02 / 2.98e-03 / 8.11e-04 | 3.03e-02 / 1.08e-02 / 2.97e-03 / 8.20e-04 | +0.1 / +0.3 / -0.5 / +1.1 % |
     | its orders | 2.59 / 2.97 / 2.93 | 2.59 / 2.98 / 2.89 | |
     | volume error at T/2 | 1.03e-01 / 2.51e-02 / 3.18e-03 / 4.65e-03 | the same to four digits | 0.0 % |
     | volume error at T | 2.14e-03 / 1.36e-03 / 1.92e-04 / 3.28e-05 | 2.08e-03 / 1.27e-03 / 1.87e-04 / 7.21e-06 | -2.4 / -6.7 / -2.4 / -78.0 % |
     | band gradient error at T/2 | 2.29 / 2.45 / 2.61 / 2.64 | the same to three digits | 0.0 % |
     Q1 (a lower volume error at T at the two finest rungs and a higher volume order) passes as
     written. But the signed volume error at T changes sign across the rungs (-, +, +, -; the same
     sign in both arms at each rung): its consecutive orders are no convergence rate, and the -78 %
     at the finest rung is an absolute shift of 2.6e-05 toward zero, against a volume error of
     4.65e-03 at T/2. Q2 (shape and band gradient within 5 %) passes.
   - Production configuration (`SL_CLIP false`), 49,910 / 347,073 / 2,389,233 cells: fails with and
     without the projection (volume error at T 7.6 / 20.1 / 25.8, band gradient at T/2 up to 6e+24);
     the projection moves the failed metrics by +0.1 to +9.4 %. No rescue, as pre-registered.
6. **Reading.** The projection removes a real defect of the prescribed flux, one that does not shrink
   with refinement, for one Poisson solve at the start (217 to 912 PCG iterations). On the 3D shear it
   does not move the shape error, the T/2 volume error or the band gradient error up to 4.8M cells,
   and it shifts the T volume residual toward zero. It does not cure the far-field failure of
   `SL_CLIP false` on polyhedra, which is attributed to the fit amplification. Adopted on polyhedral
   meshes on the author's instruction, because no metric degrades by more than 1.1 %: METHOD 8.1 row
   `FLUX_PROJECTION`.
7. **Observations, the same in both arms.** The T/2 volume error of the published configuration is
   not monotone (3.18e-03 at 1.27M cells, 4.65e-03 at 4.83M cells), and the T volume error changes
   sign across the rungs.
8. **Open.** 3Ddeformation on polyhedra (3 configs) is not measured with the projection. Raw data:
   `/work/scratch/tm83tomy/leia-fluxproj/studies/fluxProjectionGate3DshearPoly{,Clip}` on Lichtenberg
   (job ids in that worktree's `.my_jobs`).

### 11.26 The combined source F = a + lambda (1 - q): the centred q fails at zero velocity, the one-sided q passes, the kinematic gate fails on the shape order (2026-10-13)

The author's direction of 2026-10-13: the SL work waits until the error amplification of the
quadratic reconstruction is understood (task T5, sent to the SL session), and this session pushes
the combined source of D. Bothe and K. Soga, F = a + lambda (1 - q) with a = n.D.n, and builds the
knowledge base into the shared memory of the agents (cb4026eb). Everything below ran on the laptop
(2D, serial or np <= 8, OpenFOAM-v2606, development cb4026eb plus this commit). The raw studies are
in `~/OpenFOAM/data/2026-10-13-combined-source/` (5.5 GB): `e1`, `e21`, `e21-np4`, `methodGate2D`,
`probes`. Every number in this section is MEASURED unless it is marked.

1. **Gate E1: the source step with the centred q at zero velocity.** `config/gates/sourceZeroFlow2D.yaml`
   quotes the pre-registration of 2026-09-28; driver `workflow/scripts/source_zero_flow_gate.sh`,
   diagnostic `workflow/scripts/source_checkerboard_diag.py`. N = 100, u = 0, linearQ, mu in
   {27, 270, 2700} 1/s at mu dt = 0.01, band 3 and 1000 cells, psi0 = d and 1.5 d, to mu t = 40.
   - E1.1 FALSIFIED as written: no growth at e^(mu t). The sawtooth amplitude (D2(D2 psi)/64 in
     |d| <= 3h) rises from 5.2e-8 to 2e-5 by mu t = 4, stays there to mu t of about 16, and then
     grows at 0.067 mu (1.87e-4 at mu t = 40). The interface cells flatten later: at mu t = 80 the
     median psi/d of the interface cells is 0.33 (OpenFOAM, the control arm of E2.1) and 0.18 (the
     numpy copy of item 2; the two agree to 2.6e-8 up to mu t = 40 and separate in the unstable
     regime).
   - A design flaw of E1, mine: dt = 1e-2/mu makes the three mu the same run in step space, so the
     rates 1:10:100 hold by construction (0.067 mu at all three). The numpy copy answers the
     dt question instead (item 2); the OpenFOAM run E1b was not needed.
   - E1.2 (psi0 = 1.5 d) PARTLY confirmed: the zero crossings move by O(0.5 h) at u = 0 (up to
     0.91 h at mu t = 40). The mechanism "one sublattice to zero within 5/mu" is FALSIFIED: q is
     restored first (band mean 1.499 to 0.993), and both sublattices flatten later (median psi/d
     0.38 at mu t = 40).
   - E1.3 (band 1000) FALSIFIED ("B disappears, A remains"): the source on the whole domain is
     worse (q band mean 4.5 and a displacement of 1.6 h at mu t = 40).
2. **The numpy copy of the source step** (`workflow/scripts/source_step_replica.py`, a research
   instrument). At u = 0 the SL step is the identity (the quadratic fit is anchored at the cell
   value and evaluated at the cell centre). The copy reproduces E1 to 2.6e-8 after 4000 steps
   (`source_step_replica.py validate`). Its tables (`source_step_replica.py table`):
   - The damage does not depend on mu dt: 0.068 mu at mu dt = 0.001 and 0.067 mu at 0.01
     (band 3, N = 100).
   - A wider band is worse: band 12 flattens the interface cells (median psi/d 0.0000) by
     mu t = 40.
   - A seeded mode 1e-6 d (-1)^(i+j) does not grow at e^(mu t).
   - THE MECHANISM (DERIVED, consistent with the three lines above). For a sawtooth with a smooth
     envelope, delta_c = (-1)^(i+j) A, the linearised centred step is A_t = mu d n.grad A: the
     envelope moves TOWARD the zero set with speed mu |d|, while smooth errors move away from it.
     This is the reversed group velocity of the sawtooth mode under a centred difference. The
     "eigenmode at rate mu", (-1)^(i+j) d, is the unbounded envelope arriving from larger |d|; a
     band cuts it off. RETRACTED: the growth law of the technical report
     (`docs/knowledge-base/retractions/centred-q-eigenmode-growth-rate.md`; a retraction paragraph
     in the report's summary and at sec:eigenmode). Consequence: the coupled rates near mu of HL1q
     (290 and 320 to 350 1/s against 270) are not explained by the source's own map; their cause
     is OPEN again (E1.4 decides the scaling with mu).
   - The q estimators, psi0 = d, band 3, mu t = 80, N = 50 / 100 / 200 (the maximum zero-crossing
     shift in h; the order of rms(psi - d)): centred 0.93 / 0.74 / 0.66, destroyed; the first-order
     one-sided q (Godunov, or least squares over the upwind face neighbours) 0.70 / 0.29 / 0.39,
     stable but first order; the same with the interface cells anchored 0 / 0 / 0, order 1.9 and
     2.0, but a wrong interface gradient stays (psi0 = 1.5 d keeps psi/d = 1.5 there); the cell's
     own quadratic fit 0.93 / 0.68 / 0.84, as centred; the single most upwind neighbour's fit
     0.22 / 0.058 / 0.20, no convergence; an upwind quadratic over the 5x5 ring 0.16 / 0.037 /
     0.0077, order 3.1 and 3.3, not compact; the upwind-weighted mean of the upwind face
     neighbours' fits 0.034 / 0.0043 / 0.0006, order 3.7 and 3.6. A smooth band taper with the
     centred q delays the damage and does not remove it.
3. **The change.** `slGradientControlSource`: `qScheme centred | upwindFit` (token `GC_Q_SCHEME`,
   default `centred`). `upwindFit`, per cell, two passes over face neighbours: (1) the SL fit is
   updated on the current psi, and its gradient g_u and Hessian H_u at every cell centre are
   exchanged across processor patches; (2) g_c = sum_u w_u (g_u + H_u (x_c - x_u)) / sum_u w_u over
   the face neighbours u upwind of c (closer to the zero set or across it),
   w_u = sign(psi_c) (psi_c - psi_u) > 0; q = |g_c|, n = g_c/q; F = 0 in a cell without an upwind
   neighbour. `quadraticWeightedLeastSquaresReconstruction` gets `fitDerivatives` (the uncached
   fit had it). `anchorInterface false | true` (token `GC_ANCHOR_INTERFACE`, default `false`) is a
   research probe: F = 0 in every cell with a face neighbour of the opposite sign. Both tokens are
   method tokens of `methodGate2D` and `methodGate3D`. Inertness: the default path re-ran E1 case
   00000_s1 byte-identically over 4001 steps (the CSV and psi), and the rendered `qScheme centred`
   path over 4000 steps (E2.1, P1); the null arm of the anchor probe (F = 0) is byte-identical to
   the baseline in 12 of 12 cases.
4. **Gate E2.1 PASS** (`config/gates/sourceZeroFlowUpwind2D.yaml`, pre-registered from the numpy
   copy before the run). `upwindFit` at u = 0, mu t = 80:

   | psi0 = d | N = 50 | N = 100 | N = 200 |
   |---|---|---|---|
   | zero-crossing shift max / h (copy) | 0.0341 (0.034) | 0.0043 (0.0043) | 0.0006 (0.0006) |
   | zero-crossing shift RMS (L2) / h | 0.01565 | 0.002104 | 0.000299 |
   | rms(psi - d)/h in abs(d) < 3h | 1.478e-2 | 2.271e-3 | 3.790e-4 |
   | sawtooth growth from mu t = 40 to 80 (limit 2x) | 1.69x | 1.32x | 1.16x |

   Correction 2026-10-14: the max row is an L_inf quantity. It stays in the table because P3 was
   pre-registered on it, but the L2 row is the one to quote: its orders are 2.89 and 2.82
   (`source_checkerboard_diag.py`, the same edges, on the preserved data in
   `~/OpenFOAM/data/2026-10-13-combined-source/e21/`). The estimator ranking in item 2 uses the
   max of the numpy copy; the copy's L2 values are not recorded.
   The orders of rms(psi - d) are 3.70 and 3.58. From psi0 = 1.5 d the median psi/d at the
   interface is 0.9978 (N = 100) and 0.9993 (N = 200). The centred control flattens (median psi/d
   0.04 at N = 50, 0.33 at N = 100). 4 ranks against serial at N = 100: psi agrees to 1e-13 at all
   161 writes, the metric columns to 1e-17 absolute (the reduction order). The C++ estimator
   reproduces the copy to three and four digits.
5. **The kinematic arms of `methodGate2D`: CS1, CS0s and CS0r FAIL** (`config/candidates/CS1.yaml`,
   `CS0s.yaml`, `CS0r.yaml`, pre-registered before the run). CS1 = linearQ plus strain weight full
   with `upwindFit`, M_mu = 1 (lambda = 1 1/s in the shear arm, band 3); CS0s = the strain term
   alone; CS0r = the restoring term alone (HL1q without the extension, with the one-sided q). The
   baseline re-ran on the same binaries. All 36 kinematic cases COMPLETED; the seams pass (np 1 and
   8 against 4 to 1e-12). Verdicts by `workflow/scripts/make_gate_summary.py`:

   | shear arm (2Dvortex, reversed, N = 68 / 96 / 136) | gradient band error at T/2, ratio to the baseline at 136 | shape error at T, order (baseline 3.25) | volume error at T, N = 136 (baseline 2.2e-4) |
   |---|---|---|---|
   | CS1 | 0.137 | 0.88 | 3.7e-2 |
   | CS0s | 0.575 | 1.24 | 1.2e-2 |
   | CS0r | 0.302 | 0.26 | 3.3e-2 |

   The exact1D arm (u = alpha x, alpha = 1, T = 1; the band mean q at T; the baseline gives
   e^(-1) = 0.368, exact): CS1 1.106 / 0.913 / 0.667 / 0.514 and CS0s 1.116 / 0.923 / 0.676 / 0.534
   at N = 32 / 64 / 128 / 256 FALSIFY the pre-registered q = 1. CS0r gives 0.546 / 0.514 / 0.523 /
   0.461, and its closed form 1/(1 + T) = 0.5 holds. CS0r hit the gradient target (0.30), against
   its prediction to miss it.
6. **Probes of the shape loss** (each config in `config/gates/` with its read-out written before
   the run).
   - The band, exact1D with `SL_SOURCE_BAND_CELLS 1000`: the strain term keeps q = 1.000119 /
     1.000040 / 1.000010 / 1.000003 (N = 32 to 256) on the whole domain, against 0.534 with the
     3-cell band. In 1D the band is the cause. (With a uniform a the whole-domain source is a uniform
     rescaling of the baseline field, so the zero set is the baseline's.)
   - The shear arm with the whole-domain source: the gradient band error is 0.17 to 0.33 of the
     baseline, but the shape error at T at N = 136 is still 1.76e-2 (CS0s), 3.99e-2 (CS1) and
     2.84e-2 (CS0r): 32x to 73x the baseline.
   - The anchor probe (`sourceAnchorProbe2D`): A1 FALSIFIED. With the interface cells anchored the
     strain arm's shape error at N = 136 is 1.29e-2, 23x the baseline (the limit was 10x): the
     direct zero-crossing shift of the cell-wise factor is not the cause.
   - The time-step probe (`sourceDtProbe2D`): D2, a spatial error. The strain arm's shape error is
     1.76e-2 / 1.80e-2 / 3.00e-2 at CFL 0.5 / 0.25 / 0.125 (N = 136); the null arm's is 5.5e-4 /
     2.6e-4 / 1.1e-4.
   - The deformation probe (`sourceDeformationProbe2D`): H1 not confirmed. The ratio strain/null
     at N = 136 is 18.6 at T = 1, 32.1 at T = 2 and 15.3 at T = 3 (H1 needed 3x between T = 1 and
     T = 2, and a larger ratio at T = 3). Even at T = 1 the strain arm's shape error converges at
     order 0.69 (8.9e-3 to 5.5e-3), the null arm's at 1.7.
   - The direction probe (`sourceDirectionProbe2D`): P1 FALSIFIED. With n from the centred
     least-squares gradient (`GC_Q_SCHEME centred`; with law none the centred q has no restoring
     term) the shape error at N = 136 is 34.3x the baseline (1.88e-2): the upwind direction is not
     the cause.
   - The fit probe (`sourceFitProbe2D`): C1 FALSIFIED by its clip part. With `SL_CLIP true` the strain
     arm has 20.2x the clipped null arm's shape error at N = 136 (1.41e-2 against 6.99e-4): the
     unbounded quadratic fit is not the cause. Its linear-fit arms are VOID (the 2Dvortex template
     has no geometryFit entry; the source asked the linear fit for derivatives and stopped at step 1).
   - The interface at T, N = 136, from the zero crossings of psi(T) against psi(0) on the cell edges
     (serial cell order): the baseline's crossings move by 0.009 / 0.010 h rms (x / y edges),
     smoothly along the interface (lag-1 correlation 0.85 / 0.97), and 2 + 2 edges change their
     crossing; the strain arm's move by 0.25 / 0.22 h rms with lag-1 correlation 0.06 / 0.36, and 74
     and 94 edges change their crossing. The shape loss is a CORRUGATION of the interface, not a
     displacement (mean shift -0.02 h). The applied F on the interface cells at T has an rms of 1.52
     and a cell-scale part D2(D2 F)/64 of rms 0.33. (A first version of this analysis, reported in the
     chat, sorted the cells on y rounded to 12 digits and mixed the rows; its numbers are void.)
   - The strain-direction probe (`sourceStrainNormalProbe2D`, a probe key `strainNormal psi | fixedX`,
     token `GC_STRAIN_NORMAL`, default `psi`): T1 CONFIRMED. With n fixed to e_x the source is
     F = D_xx, smooth, of both signs and of the size of n.D.n, and independent of psi: its shape
     error at T is 4.42e-3 / 1.49e-3 / 6.51e-4 (orders 3.16, 2.37), 1.2x the baseline at N = 136, and
     its volume error 2.06e-4 (baseline 2.16e-4). With n from psi it is 32.1x, as before. A smooth
     multiplicative F does not corrugate the interface; the dependence of F on psi through
     n = grad psi/|grad psi| does.
7. **Assessment.** (a) The one-sided `upwindFit` q is the repair of the source's q; it is ready
   (E2.1, 4 ranks, seams). (b) A 3-cell band breaks the strain term even in 1D; the whole domain
   repairs 1D. (c) In 2D the strain term alone, on the whole domain, loses the shape order (0.7 to
   1.2) through a spatial error per step: the interface corrugates, and F carries a cell-scale part.
   Falsified as the cause: the direct crossing shift (anchor), the time step, the deformation, the
   upwind direction n, the unbounded quadratic fit (clip). The cause (MEASURED by T1): the
   dependence of a = n.D.n on psi through n = grad psi/|grad psi|; a smooth psi-independent F of the
   same size keeps the shape at 1.2x the baseline. The mechanism (DERIVED, consistent with T1 and the
   six falsifications): the linearisation of psi n.D.n contains a transport of perturbations ALONG
   the level sets with speed 2 psi (n.D.t); both q schemes are centred in that tangential direction
   (upwindFit is upwind only along the normal), the E1 disease turned by 90 degrees, and the
   low-dissipation SL step does not damp it (the Eulerian SDPLS source R, which converges
   kinematically, is integrated inside a dissipative FV transport). Next: a strain term whose
   direction does not feed back the cell-scale content of psi: a normal-constant a from the zero
   set (E2.4: n and a at the foot point of the SL fit), then the same kinematic arms. No coupled arm
   runs before a kinematic arm passes.
8. **Records.** STATUS 11.26 (this section); the knowledge base: the retraction note, the negative
   log (the centred q, the taper alone, the first-order one-sided q, the cell's own fit, the single
   neighbour's fit, the band-limited strain term, the combined source on the kinematic gate, the
   anchor), `concepts/gradient-control-next-experiments` (the outcomes), `hubs/gradient-control`,
   `concepts/combined-source-note` (re-opened), `concepts/source-discretisation-defects`,
   `models/sl-source`, `sessions/current`. METHOD.md is not affected: `SL_SOURCE none` stays the
   default and nothing of this line is in production.

### 11.27 The static test of the Eulerian source terms: setup and the laptop smoke (2026-10-14)

The author's request of 2026-10-13/14, before the work returns to advection: test the level-set
source terms on their own. Only the Eulerian sources of `leiaLevelSetFoam` that the code has, those
of D. Bothe and K. Soga included; zero velocity, so the psi equation is `ddt(psi) == S(psi)` and a
run integrates the source alone in pseudo-time tau = t (every rate 1 1/s). A circle R = 0.25 with
the centre (0.5013, 0.4987) off the cell lattice, N = 64 to 512, eight initial fields with the exact
circle as zero set (the signed distance d, r^2 - R^2, its normalised form, the blends
d + eps d^2/(2R) with eps = (h/R)^2, h/R, 0.1, and d (1 + 0.1 cos(m theta)), m = 2, 4). The sources:
the seven gradient-control laws (linearQ = Soga's restoring term and the Bothe-Fricke-Soga
bounded-gradient equation; the Bothe-Maric-Soga dossier laws; Bothe's softWall with epsD 0.8, so
that its rate at u = 0 is 1 1/s), beta, and as null checks the strain sources R, Rdiv,
RdivStrictSp and the strain half of the combined source, which are identically zero at u = 0.
Source regions: the whole domain (MOLLIFIER none, the default) and a 3-layer band. tau = 80,
dt = 0.5 h. Plan and decisions: the plan of 2026-10-14; the read-out is in each header of
`config/staticSource/` (written before any run).

1. **The case and the fields** (ea7110fb). `cases/2DstaticCircle` (velocity (0 0 0),
   `adjustTimeStep no`, `deltaT #eval{ DT_COEFF/N_CELLS }`, writePrecision 17); psi0 by
   `setExprFields` from `system/initFields/<variant>`; the reference psiEnd = d and
   alphaEnd = alpha_DA(d) by `END_REFERENCE translate` at zero velocity. Three traps found and
   fixed before any study: (a) an `#include` inside the `expressions` list of setExprFieldsDict is
   expanded in the token stream and loses its `;`; (b) a list in the parent dictionary is not found
   by the expression driver, so the variant's list is copied into the entry with
   `initVariables $initVariables;`; (c) the base `config/config.yaml` sweeps PHASE_INDICATOR over
   two values, so every study pins it.
2. **The trial** (N = 64, tau = 2, linearQ in the band, serial): from r^2 - R^2 the band mean of
   |grad psi| is 0.733 at tau = 1 and 0.878 at tau = 2, against the logistic law of the zero set
   0.731 and 0.881; the Detrixhe-Aslam volume excess over that of d falls from 1.17e-3 to 6.3e-5.
3. **The laptop smoke** (`workflow/scripts/static_source_smoke.sh`; ea7110fb with clean library
   stamps; N = 64, tau = 2; stage A and B configs at np 4, Control, GcNone and GcBand also serial):
   296 cases, all COMPLETED. VOID 2026-10-14 (part 2): 60 of its 176 np-4 cases had psi solves
   without an iteration; the smoke did not check that.
   - check-init: 296 of 296 PASS (psi0 equals the analytic variant within 1e-14 in L1 and L2, the
     signs equal those of d, sdf equals psiEnd bit for bit, psiEnd and alphaEnd in every processor
     directory).
   - Null checks, 57 pairs, all bit-identical in both CSVs and in psi and alpha at every write:
     P0.3 (alpha of r^2 - R^2 equals alpha of its normalised form), P6.1 (beta equals linearQ),
     P7.1 (R, Rdiv, RdivStrictSp equal noSource), P8.1 (linearQ with strain weight full equals
     linearQ).
   - Seam, serial against np 4, the arms P9.1 judges (noSource and the band arms): psi agrees within
     2e-17 (L2), the gradient columns within 1.6e-12 (relative), the dimensionless volume columns
     within 2.5e-15 (absolute). P9.1 as written (1e-10 relative on every column) FAILS on the
     volume columns: they are differences of two nearly equal volumes, so their relative difference
     is the noise of near-zero numbers. The extremum columns (minima, maxima, L_INF_E_PSI) are not
     compared. The parallel path is correct.
   - The whole-domain arms (MOLLIFIER none; only reported): serial and np 4 differ by up to 5.2e-4
     in the volume columns and by 17 % in the gradient columns already at tau = 2, N = 64. The
     whole-domain sources amplify round-off far from the interface (the circle centre, the walls),
     so their numbers on the cluster belong to np 4 and are not decomposition-independent.
     CORRECTED 2026-10-14 (part 2, item 8): I was wrong to say this of every whole-domain source.
     Per law, with the fix and without it: linearQ, linearZ, saturatedLinearZ and twoThirdsZReg
     agree between serial and np 4 to round-off (volume columns 2.3e-15 absolute, gradient columns
     9e-13 relative); softWall differs by 7.5e-4 and 20 %, cubicQ by 4.6e-4 relative in the gradient
     columns, cubicZ by 8e-11. The decomposition dependence belongs to softWall, and less to cubicQ.
   - Nothing was written under docs/.
4. **The pace** (laptop, np 4): 7.8e-3 s per step at N = 128 (a full run, 160 N steps, 2.7 min),
   2.2e-2 to 2.9e-2 s per step at N = 256 (15 to 20 min); N = 512 about 2 h by extrapolation.
   solve_runtime: 30 min (null checks), 90 min (stage B), 480 min (Control and stage C).
5. **Next.** A Lichtenberg worktree at the commit of this section, the build, the drivers of stage A
   and B, stage C after stage A passes; then the metrics on the cluster and the results (part 2).

#### 11.27, part 2: the first cluster run is VOID; the psi solver did nothing below its tolerance (2026-10-14)

1. **What ran.** The worktree `/work/scratch/tm83tomy/leia-static` at 010b553a, built by the CI
   script (stamp `shared-method-config-2026-09-01-330-g010b553a`, `CI: PASS`). Stage A (Control,
   NullR, NullStrain; drivers 55435601/02/03) and stage B (GcNone, GcBand, Beta; 55435604/05/06):
   496 cases, all COMPLETED. The NullR case 00017 (RdivStrictSp, blendH1, N = 128) reached the
   30-minute limit on a slow node (0.168 s per step against a median of 7.7e-3) and was re-run by
   driver 55439898. The stage A null check (job 55440912) gave 114 bit-identical pairs. Stage C (the
   N = 512 rungs; drivers 55440959/60) started at 20:02.
2. **How the defect was found.** The first metrics job (55441345, metrics commit 4025143f) gave for
   linearQ in the band on d the zero-crossing shift against tau = 0 (shift2) at tau = 80: 0.44 h at
   N = 64, 0.30 h at N = 128, and exactly 0 at N = 256. At N = 256, psi at tau = 80 was
   bit-identical to psi0: the source had done nothing in 40960 steps.
3. **The mechanism.** At u = 0 the psi equation of `leiaLevelSetFoam` is
   `ddt(psi) == Sp(F, psi)`, so the matrix is diagonal. OpenFOAM normalises the solver residual by
   sum(|A psi - A psibar| + |b - A psibar|), here of the size of the whole field divided by dt. A
   source that acts only in a band, or that is of order h^2 on a signed distance, gives an initial
   residual below `tolerance 1e-10`, and PBiCGStab returns psi unchanged ("No Iterations 0"). This is
   the cancellation rule of CLAUDE.md seen from the other side: below the tolerance the solver does
   not under-converge, it does nothing.
4. **The size**, from the solver-performance lines of all logs (the count is now the command
   `static_source_metrics.py iterations`). Of the 416 cases with a nonzero source (gradient control
   and beta), 415 have psi solves without an iteration: 19.8 of 29.0 million (68 %; median per case
   0.70; band arms 0.59 / 0.67 / 0.74 and whole-domain arms 0.77 / 0.80 / 0.82 at N = 64 / 128 / 256).
   39 cases never changed psi: every cubicQ, cubicZ and softWall band case on d and on the h^2 blend
   at every rung (G'(1) = 0 makes F tiny near q = 1), and linearQ, linearZ, saturatedLinearZ and
   beta in the band on d and on the h^2 blend at N = 256. P4.1 (the flat laws keep the zero set on
   d) would have passed for the wrong reason, and the band arm of linearQ would have "converged".
5. **Consequences.** Under the rule "A wrong setup voids its data", every number of the static
   source studies at 010b553a is void, the null checks included. On Lichtenberg the directory is now
   `studies/staticSource_VOID_psiSolverDeadZone_2026-10-14` (with `VOID_README.txt` and
   `psi_zero_iteration_counts.txt`); stage C was cancelled at 20:28 with its 112 running jobs and
   both drivers, by id from the ledger. The laptop smoke of part 1 is void as well (item 3 of
   part 1). Retraction note: `docs/knowledge-base/retractions/static-source-solver-dead-zone.md`.
6. **The fix** (edb1d5bb). The psi solver of `cases/2DstaticCircle` gets `minIter`, token
   `PSI_MIN_ITER`, 1 in `cases/2DstaticCircle.parameter`. One PBiCGStab iteration with DILU solves a
   diagonal system to round-off. With F = 0 the residual is exactly zero and PBiCGStab stops at its
   singularity test before any update, so the null arms keep psi. The command `iterations` writes
   `staticSource_iterations.csv`, and `analyse` adds the row S0 (F != 0: no psi solve without an
   iteration; F = 0: no iteration at all). The smoke script runs the case that froze (linearQ, band,
   d, N = 256, tau = 1, np 4). Also: `analyse` reads its own tables without the stopped-run rule of
   `read_csv` (it dropped a last row with an empty field), and the null studies get 60 minutes.
7. **The gate**, pre-registered in the message of edb1d5bb, on the laptop with the binaries of the
   void smoke (no C++ change since ea7110fb), 301 cases:
   - G1 PASS. The dead-zone case: 1536 psi solves, none without an iteration; psi at tau = 1
     differs from psi0 (L2 of the change 3.5e-8).
   - G2 PASS. 301 of 301 COMPLETED; check-init 301/301; S0 301/301; null 57/57 bit-identical
     (P0.3, P6.1, P7.1, P8.1); `analyse` end to end; nothing written under docs/.
   - G3 FAIL as pre-registered. The prediction "every case whose psi solves all iterated in the void
     smoke is bit-identical to it" is false: 149 of those 190 cases differ (the 41 identical ones
     are the 40 cases with F = 0 and one band case, linearZ on the normalised r^2 - R^2, serial);
     the 110 cases with zero-iteration solves differ, as predicted. My model of the solver was wrong: PBiCGStab tests minIter also at its early exit after
     the first half-step (`PBiCGStab.C` line 223 of v2606), so with minIter 1 every solve does the
     full BiCGStab step, which adds a correction of round-off size. Measured size, L2 of
     psi_new - psi_void over L2 of psi at tau = 2: at most 4.3e-11 (median about 1e-15) in every arm
     except softWall on the whole domain, where it is 0.27 (median) to 0.59 (max). That arm turns
     round-off into an O(1) change within tau = 2; see item 8.
   - G4 as in part 1. For noSource and the band arms, serial against np 4: psi within 2.6e-17 (L2),
     the gradient columns pass the comparator, the dimensionless volume columns differ by at most
     2.3e-15 (absolute), which the comparator's relative test calls FAIL on numbers of size 1e-6,
     as in part 1. So the literal P9.1 still fails on the volume columns, at round-off.
   - The pace is unchanged: 2.34e-2 against 2.15e-2 s per step (linearQ, band, N = 256) and 2.93e-2
     against 2.93e-2 (softWall).
8. **The whole-domain seam, per law** (correction of part 1, item 3): serial against np 4 at
   N = 64, tau = 2, the same in the void and the new smoke: linearQ, linearZ, saturatedLinearZ and
   twoThirdsZReg agree to round-off (volume 2.3e-15 absolute, gradient 9e-13 relative); softWall
   differs by 7.5e-4 in the volume columns and 20 % in the gradient columns; cubicQ by 4.6e-4 in the
   gradient columns; cubicZ by 8e-11. The dead zone did not cause this. softWall on the whole domain
   amplifies round-off to O(1) (item 7, G3), so its cluster numbers belong to np 4 only.
9. **Next.** A fresh worktree `/work/scratch/tm83tomy/leia-static2` at the commit of this part, a
   fresh build, the drivers of stage A and B, stage C after stage A passes; the metrics from
   `leia-static-m` with `check-init`, `iterations`, `null` and `analyse`; S0 must PASS before any
   result is read.

### 10.12 The plane test on perturbed hexahedra: planar faces exact, warped faces need the wetted area vector (2026-10-02)

Asked by the panel review of the proposal (R3, R5: every static test so far ran on blockMesh hexahedra).
`cases/davof/planeNormal3D/Allrun_perturbed vertex|column <amplitude>`; pre-registered prediction:
round-off on any perturbed mesh, because the identity holds on any polyhedron. FALSIFIED in one half:

1. **Column-wise perturbation** (`workflow/scripts/perturb_columns.py`: the same (dx, dy) for all points
   above each other, boundary points fixed, 0.2 h and 0.3 h, skewness 0.41, every face planar):
   `-expectExact` PASSES for detrixheAslam and quadraticFaces at 1e-12 (normal, position, consistency,
   plane-cut alpha). planePhaseIndicator: normal 5e-12, consistency 5.6e-13, position 1e-12, but its
   `MAX_VOL_DIFF_PLANE` diagnostic reads 5.4 on the irregular cells and fails the check: a diagnostic
   defect of the plane-indicator path on non-cubic cells, not a DAVOF result (open).
2. **Vertex-wise perturbation** (`leiaPerturbMesh -alpha 0.1|0.2 -seed 0`, faces warp): the SCALAR face
   fraction loses the plane for every initialisation, E_LINF_N 0.080 (0.1 h) and 0.168 (0.2 h),
   `max |m_c - sum caps|/h^2` 0.049 and 0.102, position 0.038 h and 0.065 h, while
   `max |alpha - plane cut|` stays at 2e-15. Cause: on a warped face the wetted area vector
   sum_t a_t S_t is not parallel to S_f, so alpha_f S_f cannot represent it; (1) holds only for the
   vector. Fix measured the same day: `davofState::wettedSf_` (per-triangle wetted fractions times
   triangle area vectors, filled by quadraticFaces) and `mVec_ = -sum_f S_f^-`; the test prints
   `vector-state normal`: L2/Linf 6.6e-15 / 6.2e-14 (0.1 h), 1.1e-14 / 1.1e-13 (0.2 h), and
   9.5e-15 / 7.9e-14 on the column mesh (identical to the scalar there).

Consequence for the method (proposal, both versions, 2026-10-02): on non-planar faces the transported
face quantity is the wetted area vector and, for the explicit position (2), its first moment: four
scalars per face that reduce to alpha_f on a planar face. The first moment and the position on warped
faces are not yet measured (open); the Detrixhe-Aslam path does not fill `wettedSf_` yet (open).

Reworked the same day on the author's instruction (no per-face storage): the wetted face area is
computed from the barycentric fan, the triangles (x_f, p_i, p_{i+1}) picked up by circulating over the
face's own points with its centre x_f, each triangle's wetted fraction times its area vector summed
straight into the owner and neighbour cells during the initialisation (`davofState::mTri_`, both the
Detrixhe-Aslam path via `faceWettedAreaVector` and the quadratic path, which replaces the DA triangle
contribution of every face whose quadratic cut succeeds). `areaNormal()` uses that sum as m_c whenever
the state came from a surface; the scalar identity -sum_f alpha_f S_f is printed by the test as
`alpha_f S_f identity` (identical on planar faces). The per-face vector field of the first version is
gone.
Measured after the rework (vertex-perturbed 8^3 plane, `Allrun_perturbed vertex 0.1|0.2`): the davof
normal from the fan sum E_LINF_N 3.3e-14 / 6.2e-14 at 0.1 h (Detrixhe-Aslam / quadratic faces) and
1.8e-13 / 1.1e-13 at 0.2 h, consistency against the tet sum 1e-15, while the scalar identity reads
0.080 and 0.168; column 0.3 h and the Cartesian Allrun unchanged (3/3). `-expectExact` still fails on
the warped meshes through the POSITION: `planePosition()` uses sum_f alpha_f (x_f - x_c).S_f, which
is inexact on a warped face; with the exact normal the position error drops from 0.038 h to 0.0037 h
(the realizability diagnostics 0.004 / 0.0076 follow from it). Next: the first moment per fan
triangle (the wet polygon's centroid of each triangle) in the same loop, nothing stored.

Position (same day, on the author's instruction): the explicit position uses the same fan. Per triangle
the first moment (x_wet,t - x_c).S_wet,t with the centroid of the wet polygon
(`davof::triNegativePolygon`: a walk along the triangle boundary collecting the wet vertices and the
edge crossings; the quadratic path keeps its exact area fraction and takes the centroid from the
linear cut of the same triangle, exact for a plane and immaterial on a planar face, where
(x - x_c).n_f is constant), summed per cell on the fly (`momentTri_`); `planePosition()` uses it
whenever the state came from a surface. Measured: vertex-perturbed 0.1 h, both states PASS at 1e-12;
0.2 h: normal 1.8e-13 / 1.1e-13, position 1.5e-12 h, plane-cut alpha 2.6e-15, area 2.6e-13, curvature
2.3e-15; column 0.3 h: position 2.2e-12 h, the rest 1e-13 and below. The two position maxima sit in
corner-sliver cells (alpha 2e-5 and 0.99998, |m| 0.002 V^(2/3)), where p = (3 alpha V - M)/|m| divides
the difference of two O(V) numbers by a tiny area: round-off, so the position check of `-expectExact`
has a ten-fold allowance (1e-11 h). In such a cell on the warped mesh the old alpha_f moment put the
plane 0.004 h off. The ellipsoid at N = 20 with quadratic faces is unchanged to every printed digit
(normal 0.0111565466185 / 0.115032544554, identical plane offsets). The plane indicator keeps its
scalar path (its own fan sums are open), as is its realizability diagnostic on non-cubic cells.

### 10.13 DAVOF static gate 3: the zero-step regeneration, measured (2026-10-04)

What the proposal's Stage 1 pre-registers, run without a time step: from the
static state (alpha_k, {alpha_f}) every interface cell recovers its normal
(the Gauss identity), its explicit plane and, from the ring of recovered
normals, its curvature tensor (`quadricFit normalsOnly`); a clipping surface
is built from the three; every face of the cell is cut again with it; the
regenerated face fractions are compared with a reference AS FRACTIONS over
the (cell, face) pairs of the interface cells (each cell with its own view;
the owner-against-neighbour disagreement on shared faces is the C0
diagnostic); the normal and the position are recovered from the regenerated
data once more (`davofState::setRegenerated`). Code:
`src/leiaLevelSet/davof/davofRegeneration.H` (header-only: `cellSurface`,
`faceCut` with the quadratic cut on every fan triangle and a four-way split
where it falls back, `cellWettedArea`, `cellLiquidVolume[Refined]`,
`solveOffset`, `refFaceFraction`, `regenerate`), option
`leiaTestDavofNormal -regenerate <arm>` / `fvSolution davof.regenerate`
(token `DAVOF_REGEN_ARM`), output `leiaTestDavofRegen.csv`, the model row
`davofRegen` and the eight worst (cell, face) pairs in the log, tables and
figures `workflow/scripts/make_davof_regen_table.py` (`davof_regen_*`, run
by the report step for a study with that CSV; the compact table exported by
`make davof-proposal`), studies `config/davof/sphereRegeneration3D.yaml` and
`ellipsoidRegeneration3D.yaml` (PREDICTION 6 in their headers), gate
`cases/davof/planeNormal3D/Allrun_regenerate`.

Arms (the clipping surface of a cell; the paraboloid in the frame of the
plane polygon's centroid with the DAVOF normal, its second derivatives from
the fit's coefficients C3..C5 / h_c):

| arm | offset along the normal | q predicted |
|---|---|---|
| plane | the explicit plane (2) | 1 |
| paraboloidPlane | zero: vertex on the plane | 1 |
| paraboloidFaces | the cell's total wetted face area reproduced (one monotone scalar, bisection to 1e-13 h); THE METHOD'S anchoring (no volume matching, Tomislav 2026-10-03) | 2 |
| paraboloidVolume | alpha_k V_k with the Detrixhe-Aslam tet quadrature of the paraboloid | 1 on the chord cell fraction |
| paraboloidVolumeExact | alpha_k V_k with the linear tet cut on tets split twice (named "own cut volume" before the ladder; see the falsification) | 1 on the chord cell fraction, 2 on the exact state |

Reference: the exact face fractions on the `exactSphere` state; otherwise
the quadratic cut of the exact distance on every fan triangle subdivided
twice (its error 16 times below the quadraticFaces initialisation).

Plane gate: every arm regenerates the plane to round-off on both
surface-based states (alpha_f 1e-14, normal 1e-14, position 1e-16 m).

**Sphere ladder, N = 20 to 160 (R/h = 5 to 40), both states, five arms**
(orders over the three finest rungs, all four rungs in parentheses; the
face-fraction error as a fraction over the (cell, face) pairs of the
interface cells; `docs/davof/davof-article/data/tables/davof_regen_orders_sphereRegeneration3D.csv`):

| state | arm | p(alpha_f, L1) | p(alpha_f, L2) | p(mismatch, L1) | p(n, L1) | p(n, L2) | p(x, L2) | alpha_f L1 at N = 160 |
|---|---|---|---|---|---|---|---|---|
| exact sphere | plane | 1.00 (0.99) | 0.99 (0.99) | 1.00 (1.00) | 2.00 (2.05) | 1.87 (1.98) | 2.00 (2.00) | 3.49e-3 |
| exact sphere | paraboloidPlane | 0.98 (0.96) | 0.98 (0.96) | 1.00 (0.98) | 2.10 (2.11) | 2.35 (2.10) | 1.97 (1.92) | 1.19e-3 |
| exact sphere | paraboloidFaces | 1.98 (1.97) | 1.91 (1.91) | 1.97 (1.97) | 2.06 (2.06) | 2.34 (2.10) | 2.00 (2.00) | 1.12e-5 |
| exact sphere | paraboloidVolume | 1.00 (0.99) | 0.98 (0.98) | 1.13 (1.22) | 2.10 (2.10) | 2.31 (2.12) | 2.01 (2.04) | 1.05e-3 |
| exact sphere | paraboloidVolumeExact | 1.05 (1.15) | 1.20 (1.30) | 1.91 (1.92) | 2.07 (2.06) | 2.35 (2.10) | 2.01 (2.04) | 7.26e-5 |
| DA + quadratic faces | plane | 1.00 (1.00) | 0.98 (0.96) | 1.00 (1.00) | 1.97 (2.00) | 1.79 (1.87) | 2.00 (2.00) | 1.29e-3 |
| DA + quadratic faces | paraboloidPlane | 1.01 (1.02) | 0.98 (0.98) | 1.01 (0.99) | 1.98 (2.02) | 1.75 (1.95) | 2.02 (2.06) | 1.97e-3 |
| DA + quadratic faces | paraboloidFaces | 2.00 (2.07) | 1.18 (1.65) | 2.00 (2.03) | 1.99 (2.03) | 1.83 (1.90) | 2.00 (2.00) | 4.55e-6 |
| DA + quadratic faces | paraboloidVolume | 1.97 (2.03) | 1.49 (1.65) | 1.98 (2.03) | 1.99 (2.02) | 1.82 (1.88) | 2.00 (2.01) | 6.25e-6 |
| DA + quadratic faces | paraboloidVolumeExact | 1.00 (1.01) | 0.98 (0.99) | 1.04 (1.12) | 1.99 (2.03) | 1.85 (1.94) | 2.02 (2.05) | 9.81e-4 |

Pairwise L1 orders of the face-anchored arm: exact state 1.96, 1.96, 2.00;
chord state 2.23, 2.07, 1.92. At N = 160 its L1 error is 300 times below
the plane's on the chord state (4.6e-6 against 1.3e-3) and 310 times on the
exact state.

**Ellipsoid ladder (the proposal's state: DA cell fraction, quadratic
faces; half-axes 1.0/0.8/0.6 mm; N = 20 to 160), the reference the
twice-subdivided quadratic cut
(`davof_regen_orders_ellipsoidRegeneration3D.csv`):**

| arm | p(alpha_f, L1) | p(alpha_f, L2) | p(mismatch, L1) | p(n, L1) | p(n, L2) | p(x, L2) | alpha_f L1 at N = 160 |
|---|---|---|---|---|---|---|---|
| plane | 0.99 (0.98) | 0.97 (0.95) | 0.98 (0.99) | 1.93 (1.92) | 1.67 (1.69) | 2.00 (1.99) | 1.68e-3 |
| paraboloidPlane | 1.01 (1.04) | 0.98 (1.00) | 0.99 (1.01) | 1.96 (1.97) | 1.85 (1.80) | 2.11 (2.12) | 2.55e-3 |
| paraboloidFaces | 2.09 (2.11) | 2.16 (2.01) | 2.01 (2.01) | 1.93 (1.93) | 1.64 (1.70) | 1.99 (2.00) | 8.84e-6 |
| paraboloidVolume | 2.04 (2.05) | 1.83 (1.77) | 2.04 (2.07) | 1.92 (1.93) | 1.61 (1.68) | 1.99 (2.00) | 1.30e-5 |
| paraboloidVolumeExact | 1.02 (1.04) | 1.01 (1.01) | 1.16 (1.24) | 1.94 (1.93) | 1.66 (1.70) | 2.04 (2.10) | 1.27e-3 |

Pairwise L1 orders of the face-anchored arm 2.12, 2.18, 2.00; its L1 error
at N = 160 is 190 times below the plane's. The original state's normal on
this ladder: 1.97 (L1) / 1.77 (L2) over the three finest rungs, the numbers
of 10.11 (1.96 / 1.78) within the fit's rounding. Anchoring residual at or
below 2e-11 relative on every rung; triangle fallbacks 0 except
paraboloidVolume at N = 20 (6) and 40 (1); no cell without a bracket.

Reading. The face anchoring regenerates the face fractions at q = 2 on the
chord state (2.00 sphere, 2.09 ellipsoid in L1) and on the exact state
(1.98); the plane and the vertex-anchored paraboloid at q = 1 (0.98 to
1.01): the tilt of the explicit plane and the half-sagitta offset of the
vertex anchoring are the O(h) fraction errors of the derivation of
2026-10-03 (article folder, FINDINGS.md). The owner-neighbour mismatch
follows the same orders (2.00 / 2.01 face-anchored, 1.00 / 0.98 plane), so
the C0 defect of the cell-local regeneration is O(h^2) as a fraction for the
method's anchoring. The regenerated normal and position are second order
for every arm: in a zero-step test the regenerated fractions are the
clipping surface's own, and (1) returns that surface's mean normal (the
plane's exactly), so the fraction order, not the normal, carries q. The
regenerated normal's 1.99 / 1.83 (sphere) and 1.93 / 1.64 (ellipsoid) in
L1 / L2 satisfy the proposal's WP1 criterion (1.8 / 1.5) trivially for the
same reason; the result of this gate is the fraction order.

Predictions against the ladders (PREDICTION 6):

- plane q = 1: confirmed (1.00 / 1.00 / 0.99).
- paraboloidPlane q = 1: confirmed (0.98 / 1.01 / 1.01).
- paraboloidFaces q = 2: confirmed (1.98 / 2.00 / 2.09 in L1).
- paraboloidVolume q = 1 on the chord cell fraction: FALSIFIED (1.97 sphere,
  2.04 ellipsoid): the Detrixhe-Aslam tet quadrature of the paraboloid
  carries the chord defect of the cell fraction it is matched to, and the
  two defects cancel (the smoke run's finding holds over both ladders). On
  the exact state the same arm is q = 1 (1.00): there the quadrature's chord
  defect has no partner.
- paraboloidVolumeExact q = 2 on the exact state: FALSIFIED (1.05, pairwise
  1.38, 1.09, 1.01). Its "own cut volume" is the LINEAR tet cut on tets
  split twice: the chord defect is 1/16 of the Detrixhe-Aslam one, not zero,
  and the matched offset carries 1/16 of the chord shift; the finest L1
  error is 7.3e-5 against 1.05e-3 for the unrefined quadrature (ratio
  14.5). A truly exact paraboloid volume needs the quadratic cut. On the
  chord state the arm is q = 1 (1.00 / 1.02) as derived (15/16 of a chord
  shift off). The arm keeps its name with this note.

Consequence for the design (now over both ladders): a volume anchor is
consistent only in the quadrature that defined the cell fraction, and it
must be the chord quadrature for the chord state; the face anchor needs no
such pairing and is second order on both states and both surfaces.

**The chord state's L2 at N = 160 (sphere, face-anchored: 1.18 over the
three finest rungs, pairwise 2.25, 2.64, -0.29; Linf 1.4e-2 against 1.1e-3
at N = 80).** The worst-pair diagnostic names one internal face (6216215,
cells 2080726 / 2080886, alpha_k 0.013 / 0.983): both cells regenerate it
at 0.4650 / 0.4654, the state's own value is 0.4651, the twice-subdivided
reference says 0.4788. A grazing face at the sphere's pole: from alpha_k
the sphere enters the empty cell by about 0.01 h through this face, so the
in-plane circle has a radius of about 0.9 h and the quadratic cut of the
fan triangles misses 1.4e-2 of its area in the initialisation; the
face-anchored regeneration reproduces the state's value to 1e-4 and inherits
the defect. The next six pairs are corner slivers (alpha_k 0.99999) with
errors of 2e-3 to 3e-3 of the same kind. Such faces exist with an expected
count of a few h/R per rung (none at N = 80, one at N = 160), so their
weight in L1 is 1e-8 and L1 is clean (2.00); L2 and Linf at the finest rung
are theirs. The ellipsoid's L2 (2.16) has no such face. Open: the exact
fraction of a grazing face from the sphere's disc for the record; the
quadratic cut on a fan triangle whose in-plane curve radius is below the
triangle size (a split criterion by curvature, not by fallback).

Regression found and fixed: `computeFromExactSphere` let the Detrixhe-Aslam
fan sums (`mTri_`, 2026-10-02) stand in for the exact face fractions; the
exact state's normal was first order from 2026-10-02 to 2026-10-04 (2.0e-2 /
9.7e-3 at N = 20 / 40 instead of 7.0e-3 / 1.6e-3). The sphere ladder of
2026-09-29 (section 10.7) and the proposal's Table 1 predate the fan.
`hasMTri_` is false for the exact state now. On this ladder the exact
state's normal reads 2.00 (L1) / 1.88 (L2) over the three finest rungs
(pairwise L2 2.16, 2.03, 1.73), the position 2.00.

Running the ladders: `--jobs 1` (an N = 160 sphere case peaks at 5.1 GB;
a sphere N = 160 case runs 20 to 60 s, the paraboloidVolumeExact arm the
longest). The morning's six parallel jobs OOM-killed one solver, which
poisoned systemd's `init.scope` cgroup of the WSL VM so that every later
session died after 90 s until `wsl --shutdown` (AGENTS.md, "The DAVOF
zero-step regeneration"). Not yet in the proposal: the manuscript changes
only on Tomislav's decision (both versions, one sentence after finding
(iv) and the Stage 1 arms as results).

### 10.14 The iterated regeneration on Lichtenberg: q = 2 per step, q = 1 per period (2026-10-06)

The ladders of `config/davof/sphereRegenIterLong3D.yaml` (K = 5000, N = 20 and 40),
`sphereRegenIter3D.yaml` and `ellipsoidRegenIter3D.yaml` (K = 500, N = 20 to 160),
arms plane / paraboloidFaces / paraboloidVector, views own / cell, origin volume,
refit on, run in the clone `/work/scratch/tm83tomy/leia-davof` (jobs 55339952 to
55339954, resubmitted 2026-10-06 morning after the YAML boolean failure of the
overnight run: `DAVOF_REGEN_REFIT: [true]` became the token `True`; quote booleans).
Tables: `docs/davof/davof-article/data/tables/davof_regen_iter_<study>.csv` and
`davof_regen_iter_orders_<study>.csv`; per case `leiaTestDavofRegenIter.csv`.

Long ladder, completed 20:31 (L1 face-fraction error as a fraction; the per-period
accounting compares k proportional to 1/h: N = 20 at k = 40 against N = 40 at k = 80):

| arm / view | E(1) at N 20, 40 | E(2N) | E(4N) | E(5000) | per-period order |
|---|---|---|---|---|---|
| plane / cell | 3.92e-3, 1.98e-3 | unchanged | unchanged | unchanged | exact fixed point (P1 confirmed) |
| plane / own | 3.92e-3, 1.98e-3 | 2.48e-2, 1.52e-2 | 2.33e-2, 1.57e-2 | 2.32e-2, 1.60e-2 | saturates in 40 iterations at order 0.5; normal destroyed |
| paraboloidFaces / cell | 2.30e-4, 4.45e-5 | 5.52e-3, 2.98e-3 | 1.05e-2, 6.00e-3 | 9.37e-2, 8.97e-2 | 0.89 (2N), 0.81 (4N) |
| paraboloidFaces / own | 2.30e-4, 4.45e-5 | 6.25e-3, 1.63e-3 | 1.76e-2, 4.33e-3 | 0.37, 0.37 | diverges; 51 and 163 cells drop out |
| paraboloidVector / cell | 2.10e-4, 2.84e-5 | 8.21e-4, 4.80e-4 | 1.34e-3, 8.70e-4 | 8.0e-3, 7.5e-3 | 0.78 (2N), 0.62 (4N); constant 6 to 7 times lower |
| paraboloidVector / own | 2.10e-4, 2.84e-5 | 3.19e-2, 1.38e-2 | 6.21e-2, 4.14e-2 | 0.31, 0.36 | diverges; Newton fails in 60 and 277 cells |

Verdicts on the pre-registered lines of the config headers: P1 CONFIRMED (plane/cell
constant to every digit over 5000 iterations). P1' (plane/own increment first order):
the arm saturates, the increment is not the measure; the saturated level is order 0.5.
P2' FALSIFIED beyond the first few hundred iterations: at fixed K = 5000 the error of
paraboloidFaces is independent of h (9.4e-2 against 9.0e-2); the increment is O(h^2)
only at the start (smoke run, STATUS 10.13 addendum of 2026-10-05). P2'' half: the
vector arm's drift is a fifth to a seventh of the scalar arm's but not second order; the
Newton failures (11 to 15 % of the cells, scalar fallback) are the suspected residual.
The measured statement for the proposal: the zero-velocity regeneration with the
face-anchored paraboloid is q = 2 per step and q = 1 per period, the alpha_k-positioned
plane with cell views is the only fixed point, and a fixed point of the second-order
map needs a clipping surface that reproduces the cell's fractions (the proposal's
Stage 1 projection; a six-coefficient paraboloid fitted to the F fractions is the
identity of the map by construction for F <= 6 and is the next arm to test,
`paraboloidFull`, not yet implemented).

K = 500 ladders (sphere and ellipsoid, N = 20 to 160): at 23:30 at 94 % and 89 %, the
N = 160 cases running; their per-period orders over four rungs, the displaced-face cut
at k = 1 (P4) and the Newton-failure share against N are to be appended here when the
drivers end (`cluster_status2.sh` of the session record; the data return by git from
the clone). The proposal (commit 50451c2 of `research/projects`) carries the long
ladder's result in finding (v) and Stage 1.

**K = 500 ladders, completed 2026-10-06 23:39 (sphere) and 2026-10-07 00:24 (ellipsoid),
48 cases, no error; leia 0d2bffc carries the tables and figures.** L1 face-fraction error
as a fraction at k = 2N (one period at fixed Courant number, k proportional to 1/h) for
N = 20, 40, 80, 160 and the per-period order between successive rungs (4N lies inside
K = 500 only for N <= 80):

| arm / view | surface | E(2N) at N = 20, 40, 80, 160 | order at 2N | order at 4N |
|---|---|---|---|---|
| plane / cell | sphere | 3.92e-3, 1.98e-3, 9.92e-4, 4.96e-4 (= E(1), every digit) | exact fixed point | exact fixed point |
| plane / cell | ellipsoid | 5.48e-3, 2.82e-3, 1.40e-3, 7.05e-4 (= E(1)) | exact fixed point | exact fixed point |
| plane / own | sphere | 2.48e-2, 1.57e-2, 1.42e-2, 9.66e-3 | 0.66, 0.14, 0.56 | 0.57, 0.15 |
| plane / own | ellipsoid | 2.99e-2, 1.87e-2, 1.36e-2, 1.04e-2 | 0.68, 0.46, 0.38 | 0.63, 0.47 |
| paraboloidFaces / cell | sphere | 5.52e-3, 2.98e-3, 1.60e-3, 8.67e-4 | 0.89, 0.90, 0.89 | 0.81, 0.84 |
| paraboloidFaces / cell | ellipsoid | 1.29e-2, 6.30e-3, 3.40e-3, 1.85e-3 | 1.03, 0.89, 0.88 | 1.14, 0.81 |
| paraboloidFaces / own | sphere | 6.25e-3, 4.33e-3, 3.92e-3, 4.89e-3 | 0.53, 0.14, -0.32 | 0.56, -0.20 |
| paraboloidFaces / own | ellipsoid | 1.51e-2, 9.65e-3, 8.76e-3, 7.82e-3 | 0.64, 0.14, 0.16 | 0.40, 0.01 |
| paraboloidVector / cell | sphere | 8.21e-4, 4.80e-4, 2.50e-4, 1.26e-4 | 0.78, 0.94, 0.98 | 0.63, 0.86 |
| paraboloidVector / cell | ellipsoid | 2.34e-3, 1.12e-3, 5.58e-4, 2.95e-4 | 1.06, 1.00, 0.92 | 1.12, 0.94 |
| paraboloidVector / own | sphere | 3.19e-2, 4.14e-2, 9.36e-2, 1.47e-1 | diverges (-0.4 to -1.2) | diverges |
| paraboloidVector / own | ellipsoid | 2.75e-2, 4.82e-2, 1.01e-1, 1.69e-1 | diverges (-0.8 to -1.1) | diverges |

The report's own order table (`davof_regen_iter_orders_*.csv`, quantity EK = the error
at K = 500, three finest rungs) reads 1.71 (sphere) / 1.69 (ellipsoid) for
paraboloidFaces/cell and 1.76 / 1.83 for paraboloidVector/cell: at FIXED K the error
still falls with h because K = 500 is 1.6 periods at N = 160 and 12.5 at N = 20; the
per-period accounting above compares like with like and is the one to quote.

Verdicts. (1) The one-scalar face anchoring drifts at order 0.9 per period on both
surfaces over four rungs (the long ladder's 0.89 / 0.81 on the two coarsest meshes was
not a coarse-mesh artefact); (2) the area-vector anchoring drifts at 0.8 to 1.0 per
period with a constant 5.5 to 6.9 times lower (sphere 6.2 to 6.9, ellipsoid 5.5 to
6.3), tending to first order under refinement (0.98 and 0.92 on the finest pair), so
P2'' stays half-confirmed: a lower constant, not a higher order; (3) the own-view
(single-valued) arms do not converge at all: plane/own saturates, both paraboloid/own
arms diverge, the vector arm fastest; (4) the Newton failures of the vector arm
(N_NO_BRACKET at k = 1, scalar fallback) are a constant share of the interface cells,
60/472, 277/1863, 1155/7418, 4651/29739 = 12.7, 14.9, 15.6, 15.6 % on the sphere and
51/298, 208/1168, 851/4649, 3400/18642 = 17.1, 17.8, 18.3, 18.2 % on the ellipsoid;
they do not vanish under refinement, so they are a property of the anchoring (the
in-plane vertex position has no root in those cells), not the residual of a
second-order arm. (5) P4, the displaced-face cut at k = 1 (`leiaTestDavofRegen.csv`,
columns `E_ALPHAF_SHIFT_{M050,M025,P025,P050}_L1`, faces displaced by -0.5, -0.25,
+0.25, +0.5 h along the normal and cut against the k = 1 surface): the face-anchored
paraboloid gives 5.6e-4, 8.7e-5, 1.4e-5, 3.1e-6 at -0.5 h on the sphere (orders 2.7,
2.6, 2.2) and 1.7e-3, 2.6e-4, 4.7e-5, 1.0e-5 on the ellipsoid (2.7, 2.5, 2.2), the
vector arm the same within 10 %, the plane 1.05e-2, 4.8e-3, 2.2e-3, 1.06e-3 (1.1, 1.1,
1.0); the other three displacements agree within 30 %. So the ONE-STEP cut against
displaced faces (reviewer C's variant of the transport step) is q = 2 for the
paraboloid arms and q = 1 for the plane; it is the repetition, not the displacement,
that costs the order. P4 CONFIRMED.

Manuscript (2026-10-09, both versions, finding (v)): "first order (0.9 over four
refinements on the sphere and on the ellipsoid) for the face-anchored paraboloid and,
with a six times smaller constant, for the paraboloid anchored by the cell's wetted
area vector (0.8 to 1.0)"; the Detailed version names the ladders (sphere and
ellipsoid, N = 20 to 160, K = 500, and 5000 on the two coarsest sphere meshes).

### 10.15 The six-coefficient clipping surface (`paraboloidFull`): the identity of the zero-velocity map (2026-10-12)

The fifth proposal review (2026-10-12) read the Stage 1 projection as inert (a trigger
above the O(h^2) floor cannot act on an O(h^2) per-step drift) or undefined (target set,
degrees of freedom, one value per shared face). The record's answer, named in 10.14, is
the paraboloid whose six coefficients (offset, two slopes, the three tensor components,
in the frame of the recovered normal) are fitted to the cell's F face fractions, with the
ring tensor as prior: it reproduces the fractions by construction, so the zero-velocity
regeneration is the identity. Implemented 2026-10-12 as the arm `paraboloidFull` (leia
36a3667, a9e6798): `cellSurface` gains the slopes a1, a2 (value = w - c - a1 u - a2 v -
(b11 u^2 + 2 b12 u v + b22 v^2)/2, negative in the liquid; the other arms unchanged to
every digit, plane gate `Allrun_regenerate` all seven arms PASSED); unknowns
(c/h, a1, a2, b11 h, b12 h, b22 h) from the `paraboloidFaces` start, residual the
|S_f|-weighted difference between the face cut and the cell's own view of its fractions,
forward-difference Jacobian, normal equations with Marquardt damping relative to the
diagonal (the prior enters as damping, not as penalty rows: with four or five cut faces
the system is under-determined and a penalty row cannot reach 1e-10; the damped step
keeps the null-space components at the start, i.e. the ring tensor, zero slopes and the
bisection offset), step limited to h, halving up to six times, at most 50 iterations,
converged at a maximum face residual of 1e-10; on failure the best iterate is kept and
the cell counted in N_NO_BRACKET. New columns E0_ALPHAF_L1 / LINF in
`leiaTestDavofRegen.csv`: the input state's own fraction error against the reference
(k = 0).

Pre-registered (config headers, before the run): P7 zero-step identity (E at k = 1 at
round-off where the fit converges, chord and exact state); P8 non-converging share below
2 % and falling with refinement; P9 the displaced-face cut at k = 1 q = 2 with a constant
at or below the face-anchored arm's; P10 iterated to K = 500 the error stays at
round-off and the regenerated normal equals the state's.

Laptop smoke, sphere N = 20 (R/h = 5, 472 interface cells; feedback own, origin volume,
refit on; `studies/davof/paraboloidFullSmoke`, git-ignored):

| arm, state | E_ALPHAF_L1 (k = 1) | N_NO_BRACKET / 472 | MISMATCH_L1 | E_L2_N (k = 1) | SHIFT L1 at -0.5 / -0.25 / +0.25 / +0.5 h | E_ALPHAF_L1 (k = 50) |
|---|---|---|---|---|---|---|
| paraboloidFaces, chord | 2.304e-4 (= 10.13) | 0 | 4.6e-4 | 2.85e-2 | 5.6e-4 / 4.9e-4 / 4.3e-4 / 4.5e-4 | 8.54e-3 |
| paraboloidFull, chord | 9.616e-5 = E0 (9.617e-5) | 11 (2.3 %) | 1.8e-7 | 6.80e-3 | 7.4e-4 / 5.3e-4 / 1.04e-3 / 1.07e-3 | 9.626e-5 |
| paraboloidFaces, exact | 2.61e-4 | 0 | 6.0e-4 | 1.03e-2 | 4.5e-4 / 4.1e-4 / 4.3e-4 / 4.2e-4 | - |
| paraboloidFull, exact | 5.7e-7 | 13 (2.8 %) | 9.0e-7 | 6.95e-3 | 5.2e-4 / 3.9e-4 / 1.20e-3 / 1.33e-3 | 3.7e-6 |

Reading. P7 holds in its measurable content: on the exact state E(1) = 5.7e-7 comes
entirely from the 13 non-converged cells (Linf 4.2e-4 is one of them); on the chord
state E(1) equals E0 to four digits, because E_ALPHAF measures against the twice-refined
reference and not against the state, so "round-off" is unreachable there by
construction; the owner-neighbour mismatch is at 1e-7 (the other arms 5e-4). The
regenerated normal is the state's normal (E_L2_N 6.80e-3 = the davof row's 6.803e-3;
the face-anchored arm's 2.85e-2 is four times worse). P8 marginal at N = 20 (2.3 and
2.8 %; the failures are faces with tiny targets, 1e-7 to 1e-5, that the surface misses,
and six-cut-face cells with an irreducible least-squares residual of order 1e-5); the
ladder decides the trend. P9 FALSIFIED on the constant: the displaced-face errors are
1.3 to 2.4 times the face-anchored arm's (the fit is exact on the faces and worse away
from them; the order is for the ladder). P10 at K = 50: drift 1e-7 (chord) and 3e-6
(exact) in L1 from the failed cells feeding back through shared faces, the normal
constant to five digits; the identity of the map holds on every converged cell.

Ladders submitted 2026-10-12 (jobs 55408801 sphere, 55408802 ellipsoid; the K = 500
studies of 2026-10-06 moved to `*.k500-2026-10-06`; all four arms rerun): the
non-converging share against N, the displaced-face orders, E(K = 500) against E0 and
the mismatch over the rungs are to be appended here.

**Ladders done 2026-10-12 (ellipsoid 16:45, sphere 18:09; 64 cases, no error; leia
4e613f9 carries the tables and figures; all four arms rerun, the 6 October results of the
three old arms reproduced to every digit).** paraboloidFull, cell view, L1 face-fraction
error as a fraction; E0 = the input state's own error against the twice-refined reference;
E(k) after k applications; the fit's non-converged share N_NO_BRACKET / N_INTERFACE at
k = 1; the owner-neighbour mismatch at k = 1; the displaced-face cut at -0.5 h with the
paraboloidFaces value in parentheses:

| surface | N | E0 | E(2N) | E(500) | non-converged | mismatch | shift -0.5 h |
|---|---|---|---|---|---|---|---|
| sphere | 20 | 9.62e-5 | 9.99e-5 | 9.99e-5 | 3/472 (0.6 %) | 1.2e-5 | 7.4e-4 (5.6e-4) |
| sphere | 40 | 1.25e-5 | 1.28e-5 | 1.28e-5 | 1/1863 (0.1 %) | 7.1e-7 | 1.13e-4 (8.7e-5) |
| sphere | 80 | 2.03e-6 | 2.03e-6 | 2.03e-6 | 5/7418 (0.1 %) | 8.1e-9 | 1.67e-5 (1.42e-5) |
| sphere | 160 | 7.48e-7 | 7.47e-7 | 7.47e-7 | 15/29739 (0.1 %) | 6.5e-9 | 3.78e-6 (3.07e-6) |
| ellipsoid | 20 | 1.49e-4 | 1.51e-4 | 1.51e-4 | 0/298 (0.0 %) | 1.9e-6 | 1.96e-3 (1.69e-3) |
| ellipsoid | 40 | 3.30e-5 | 3.31e-5 | 3.31e-5 | 2/1168 (0.2 %) | 1.6e-7 | 3.39e-4 (2.62e-4) |
| ellipsoid | 80 | 5.82e-6 | 5.83e-6 | 5.83e-6 | 9/4649 (0.2 %) | 4.6e-8 | 6.76e-5 (4.65e-5) |
| ellipsoid | 160 | 1.38e-6 | 1.38e-6 | 1.38e-6 | 17/18642 (0.1 %) | 1.9e-8 | 1.50e-5 (1.03e-5) |

Orders of the displaced-face cut (-0.5 h, successive rungs): paraboloidFull 2.71 / 2.75 /
2.15 (sphere), 2.53 / 2.33 / 2.17 (ellipsoid); paraboloidFaces 2.69 / 2.62 / 2.21 and 2.69 /
2.49 / 2.17; paraboloidVector 2.95 / 2.99 / 2.41 and 2.80 / 2.59 / 2.12; plane 1.14 / 1.14 /
1.03 and 1.19 / 1.12 / 1.05. The own-view paraboloidFull has the same E(2N) and E(500) as
the cell view (9.63e-5 ... 7.48e-7 and 1.51e-4 ... 1.39e-6), i.e. the arm is the identity in
the single-valued view as well, where every other arm saturates or diverges (10.14).

Verdicts. P7 CONFIRMED in its measurable content: E(1) = E(2N) = E(500) = E0 to three
digits on every rung and both surfaces (the error against the reference cannot fall below
the input state's own error; the mismatch 1e-5 to 1e-8 is the fit's residual). P8
CONFIRMED: 0.0 to 0.6 % non-converged cells, 0.1 % on the three finer meshes (the laptop
smoke's 2.3 to 2.8 % at N = 20 came from its different fallback settings and is
superseded). P9 half: the displaced-face cut is q = 2 (2.2 to 2.8) but its constant is
1.2 to 1.5 times the paraboloidFaces arm's on every rung (the interpolating surface is exact
on the faces and worse away from them); the vector arm has the lowest displaced-face
constant of all. P10 CONFIRMED: over 500 applications the error does not move from E0,
and the regenerated normal equals the state's except in the non-converged cells (the N =
20 sphere L2 of 2.8e-2 is three such cells). Statement for the record and the manuscripts:
the paraboloid with six coefficients fitted to the cell's face fractions IS the fixed point
of the zero-velocity regeneration at every resolution tested; the transport question is
the remainder of the fit at the displaced pre-image, O(h^2) per step with a constant
slightly above the face-anchored arm's, and whether its sum over a period is O(h^2) or
O(h). Both proposal versions (finding (v)) and the preprint (Section 6.3) carry it since
2026-10-12 evening.

### 10.16 The static recovery on perturbed-hexahedral and polyhedral meshes (2026-10-12)

Part B of the plan of 2026-10-12: the static DAVOF recovery (normal, position, curvature)
of the plane, the sphere and the ellipsoid on three mesh classes, hexahedral (the ladders
of 10.7 to 10.11), vertex-perturbed hexahedral (`mesh: perturbed`: blockMesh, then
`leiaPerturbMesh` with its defaults, every point displaced by up to 0.1 h with the
deterministic seed 0, so the mesh rule needed no change; the faces warp) and polyhedral
(`mesh: poly`: cfMesh `pMesh` on the 4 mm box `box4mm-faces.stl`, one solid per box face
from `make_box_stl.py --solids faces`, so every box edge is a feature edge and the patches
are left/right/front/back/bottom/top, covered by one regex group in the `0.org` fields).
Configs `config/davof/{sphere,ellipsoid}Normal3D_{perturbed,poly}.yaml` and
`planeNormal3D_{hex,perturbed,poly}.yaml` (the plane on the N = 20/40/80/160 ladder of
the other two, case `planeNormal3DLadder`, whose `fvSolution.template` gained the
curvature block); `cases/davof/<case>/system/meshDict.template` and
`cases/davof/<case>_poly.parameter`; the gate `cases/davof/planeNormal3D/Allrun_poly`
(`system/meshDict`); alpha sources quadraticFaces and detrixheAslam on the non-hex meshes
(exactSphere is hex-only). Two derived tokens in `materialize.py`: `DOMAIN_CELL_SIZE` =
DOMAIN_LENGTH/N_CELLS (h) and `CFMESH_CELL_SIZE` = h/2^(-1/3), the cfMesh maxCellSize
whose dual cells at the interface measure h (the pin of the README's polyhedral droplet),
so the polyhedral studies run the SAME N_CELLS ladder as the hex ones. Measured on the
plane box at N = 20: with maxCellSize = h the mesh has 31471 cells and a mean cell size
(DELTA_X) of 0.63 h, with the pin 13017 cells, DELTA_X 0.85 h and 784 interface cells
against 772 on the hex mesh, i.e. the interface resolution of the hex ladder. The
Snakefile no longer passes `-fluxCorrection` to `leiaTestDavofNormal` on perturbed
meshes (an unknown option is a FatalError in argList). `make_davof_normal_table.py` keys
every series by (mesh, alphaSource, model) with the mesh class from `case_params.json`
(the per-study outputs are byte-identical to the old script's on every study present,
checked) and writes, for a study `<surface>Normal3D[_hex|_perturbed|_poly]`, the combined
table `tables/davof_normal_proposal_<surface>Normal3D_meshes.tex` (mesh column first, a
`\midrule` between the mesh blocks, comparator rows collapsed per block when the states
share alpha_k) and `davof_normal_orders_<surface>Normal3D_meshes.csv` from the sibling
study directories that exist (a missing one is noted); `make davof-proposal` exports them
with the other `davof_normal_proposal_*.tex`.

Pre-registered (the config headers, before any run). P11 (plane, every mesh class):
normal, position and curvature at round-off (1e-12) for both states with the fan sums.
P12 (perturbed 0.1 h and polyhedral): the quadratic-face normal converges at second order
with a larger constant than on the hexahedral mesh; the Detrixhe-Aslam normal at first
order; the position at second order for both; the one-ring curvature converges (order
between 1 and 2) with a smaller constant on polyhedra, where the point-neighbour ring is
larger. P13: plicRDF stays at order about 1.1 on every mesh class; gradAlpha does not
converge; isoAlpha 0.5 to 0.7.

The plane gate at N = 20 (`-expectExact`, tolerance 1e-12, position 1e-11 h), both
surface-based states (detrixheAslam / quadraticFaces):

| mesh | cells (interface) | max skewness | normal Linf | position Linf / h | consistency | curvature Linf h | area diagnostic | verdict |
|---|---|---|---|---|---|---|---|---|
| hexahedral | 8000 (772) | 0 | 5.6e-14 / 5.6e-14 | 1.3e-12 / 1.3e-12 | 2e-15 | 1e-13 | 6e-14 / 1e-13 | PASSED |
| perturbed 0.1 h | 8000 (772) | 0.13 | 2.2e-13 / 2.2e-13 | below 1e-11 | 1e-14 | 1e-13 | below 1e-12 | PASSED |
| polyhedral | 13017 (784) | 2.18 | 1.9e-13 / 1.9e-13 | 1.9e-12 / 6.8e-12 | 6.6e-15 / 8.2e-15 | 1.9e-14 / 2.5e-14 | 2.2e-12 / 1.1e-11 | area diagnostic above 1e-12 |

(The plane indicator is not exact on perturbed or polyhedral cells, normal Linf 0.098
and 0.17, as in 10.12.) P11 HOLDS in its stated content on all three mesh classes: the
normal, the position and the curvature of the plane are at round-off on polyhedra too,
and the Detrixhe-Aslam tets and the fan triangles needed no change for polyhedral cells
(no code change in `davofState`). The ONE quantity above the tolerance on the polyhedral
mesh is the area realizability diagnostic MAX_AREA_DIFF_PLANE, 2.2e-12 and 1.1e-11
against 1e-12: it divides by h^2 in the corner-sliver cells of a mesh with skewness 2.2,
the same floor for which the position already carries a ten-fold allowance;
`Allrun_poly` therefore exits 1, and whether the area diagnostic gets the position's
allowance is Tomislav's decision (the gate's tolerance was not changed).

Laptop smoke of the sphere ladders at N = 20 and 40 (`--jobs 1`; studies
`sphereNormal3D_{perturbed,poly}Smoke`, git-ignored; orders between the two rungs; the
hex row from the ladder of 10.10):

| mesh, state | L1(n) 20 to 40 (order) | L2(n) 20 to 40 (order) | Linf(n) at 40 | L2(x) order | plicRDF L2 (order) | gradAlpha / isoAlpha L2 order |
|---|---|---|---|---|---|---|
| perturbed, quadraticFaces | 2.90e-3 to 1.36e-3 (1.09) | 5.83e-3 to 2.84e-2 (-2.28) | 1.22 | 1.44 | 3.55e-2 to 1.76e-2 (1.01) | 0.11 / 0.39 |
| perturbed, detrixheAslam | 1.25e-2 to 6.27e-3 (1.00) | 1.93e-2 to 9.93e-3 (0.96) | 0.042 | 2.06 | same alpha_k | |
| polyhedral, quadraticFaces | 6.58e-4 to 1.26e-4 (2.38) | 2.22e-3 to 3.87e-4 (2.52) | 0.010 | 2.07 | 1.83e-2 to 8.97e-3 (1.03) | 0.12 / 0.16 |
| polyhedral, detrixheAslam | 2.78e-3 to 1.22e-3 (1.19) | 5.77e-3 to 2.53e-3 (1.19) | 0.019 | 2.01 | same alpha_k | |
| hexahedral, quadraticFaces | 3.03e-3 to 6.83e-4 (2.15) | 6.80e-3 to 1.56e-3 (2.13) | 0.025 | 2.01 | | |

Reading. The perturbed quadratic-face L2 at N = 40 is ONE cell: cell 28460 has
alpha_k = 0 (the Detrixhe-Aslam tets see no liquid) while the quadratic cut of its warped
fan triangles wets faces, so |m_c| passes the wisp threshold and the "normal" of a cell
without volume is scored (e = 1.22; N_FACE_FALLBACK 4 of 3.8e5 faces). Without that cell
L2 = 1.47e-3 and L1 = 7.06e-4, the hex constants (1.56e-3, 6.83e-4), and the pairwise
orders read 1.99 (L2) and 2.04 (L1): the perturbation of 0.1 h costs nothing in the
constant, and P12's "larger constant" is not seen; the ladder decides whether such
volume-less wetted cells recur (they are the realizability defect of the inconsistent
quadratic-face state in its sharpest form, and a candidate for the two-piece/defect rule
of the proposal). The polyhedral quadratic-face normal is second order with a constant
4.0 times BELOW the hex one at N = 40 (3.87e-4 against 1.56e-3; more faces per cell enter
the sum) and its position is second order: P12 FALSIFIED in the favourable direction for
the constant on the smoke. The Detrixhe-Aslam normal is first order on both mesh classes
(1.00, 1.19), its position second order (2.06, 2.01). P13 holds on the smoke: plicRDF
1.01 and 1.03 (its error 12 and 23 times the quadratic-face DAVOF error at N = 40),
gradAlpha 0.11 and 0.12, isoAlpha 0.39 and 0.16 (below the predicted 0.5 to 0.7 on two
rungs). The curvature part of P12 waits for the ellipsoid ladders (no curvature model on
the sphere).

PENDING: the ladders on Lichtenberg (`sbatch --parsable -J leia-davof-<cfg>
--export=ALL,STUDY=davof/<cfg> run-studies.sbatch` for sphereNormal3D_perturbed,
sphereNormal3D_poly, ellipsoidNormal3D_perturbed, ellipsoidNormal3D_poly,
planeNormal3D_hex, planeNormal3D_perturbed, planeNormal3D_poly; the laptop needs
`LEIA_PMESH_PRELOAD` for pMesh, the cluster does not); their orders over the three finest
rungs, the combined tables per surface and the verdicts on P12 and P13 are to be appended
here.

**Ladders on Lichtenberg (2026-10-13; leia 072de77 carries the tables; the four hex and
perturbed studies done, jobs 55421098, 55421100, 55421102, 55421103; the three polyhedral
studies resubmitted with `solve_mem_mb: 24000` after two OOM kills at N = 160, 4000 and
12000 MB, jobs 55422221-23, pending).** Orders over N = 20 to 160 (three finest in
parentheses), errors at N = 160:

| study | state / model | p1(n) | p2(n) | L2(n) | Linf(n) | p(x) | L2(d) [m] | p(kappa) | kappa err |
|---|---|---|---|---|---|---|---|---|---|
| sphere perturbed 0.1 h | DAVOF, quadratic faces | 2.12 (2.49) | 2.40 (4.11) | 9.49e-5 | 3.43e-3 | 2.07 | 7.77e-8 | - | - |
| | DAVOF, Detrixhe-Aslam | 0.99 (0.98) | 0.97 (0.98) | 2.57e-3 | 1.23e-2 | 2.04 | 1.77e-8 | - | - |
| | plicRDF | 0.90 (0.83) | 0.93 (0.87) | 5.25e-3 | 1.45e-1 | 1.50 | 2.37e-7 | - | - |
| | gradAlpha / isoAlpha | 0.06 / 0.17 | 0.05 / 0.20 | 6.71e-2 / 4.75e-2 | 1.45e-1 / 1.97e-1 | 1.51 / 1.57 | 2.40e-7 / 2.03e-7 | - | - |
| ellipsoid perturbed 0.1 h | DAVOF, quadratic faces | 2.02 (2.00) | 1.87 (1.91) | 2.81e-4 | 1.06e-2 | 1.99 | 1.04e-7 | 1.53 | 4.65e-3 |
| | DAVOF, Detrixhe-Aslam | 1.04 (1.03) | 1.01 (1.02) | 3.80e-3 | 3.01e-2 | 2.06 | 2.44e-8 | 1.27 | 7.23e-3 |
| | plicRDF | 1.06 (0.98) | 1.02 (1.01) | 6.43e-3 | 1.39e-1 | 1.64 | 2.39e-7 | - | - |
| | gradAlpha / isoAlpha | 0.10 / 0.29 | 0.11 / 0.31 | 6.86e-2 / 4.81e-2 | 1.45e-1 / 1.93e-1 | 1.63 / 1.78 | 2.46e-7 / 2.06e-7 | - | - |
| plane hex (ladder) | DAVOF, both states | round-off | round-off | 1.8e-14 / 9.8e-15 | 3.3e-13 | - | 4.6e-17 | - | - |
| | plicRDF / gradAlpha / isoAlpha | - | - | 2.51e-4 / 7.46e-2 / 4.18e-2 | 1.89e-2 / 0.40 / 0.54 | 0.12 / 1.31 / 1.45 | 8.3e-12 / 3.0e-8 / 3.6e-8 | - | - |
| plane perturbed 0.1 h | DAVOF, both states | round-off | round-off | 1.9e-14 / 1.2e-14 | 3.6e-13 | - | 4.6e-17 | - | - |
| | plicRDF / gradAlpha / isoAlpha | - | - | 5.26e-3 / 7.56e-2 / 6.08e-2 | 0.38 / 0.41 / 0.56 | 0.95 / 0.96 / 1.02 | 2.7e-7 / 2.6e-7 / 2.2e-7 | - | - |

Reading. P11 CONFIRMED on hex and perturbed (the plane at 1e-14 for every N; the
polyhedral ladder pending, the N = 20 gate passed). P12 half: the orders hold (quadratic
faces second order, Detrixhe-Aslam first, position second, the ellipsoid curvature 1.53
against 1.55 on hex), but the constant is NOT larger on the perturbed mesh (ellipsoid L2
2.81e-4 against 3.01e-4 on hex, and the L2 order rises from 1.78 to 1.87 because the
perturbation breaks the pole alignment that set the hex L2): the "larger constant" half
of P12 is falsified in the favourable direction. P13 CONFIRMED and sharpened: on the
perturbed meshes plicRDF's order falls to 0.90 / 0.93 (sphere) at a maximum normal error
of 0.145, and the comparators' position order falls to 1.5 to 1.8, while the DAVOF
position stays at 2.0 and its maximum normal error at 3.4e-3 (sphere) and 1.1e-2
(ellipsoid); on the perturbed plane plicRDF's normal error is 5.3e-3 with a maximum of
0.38 against 1e-14. The sphere's quadratic-face L2 order 2.40 (4.11 over the three finest)
is the N = 20 rung's volume-less cell of the smoke run (alpha_k = 0 with wetted warped fan
triangles) dropping out of the ladder as the meshes refine; L1 (2.12) is the order to
quote.

**Polyhedral ladders done 2026-10-13 11:20-11:27 (jobs 55422221-23, `solve_mem_mb` 24000;
leia af4c58f carries the tables and the combined `_meshes` tables).** Orders over
N = 20 to 160 (three finest in parentheses), errors at N = 160, cfMesh polyhedra with the
mean cell size of the hex ladder:

| study | state / model | p1(n) | p2(n) | L2(n) | Linf(n) | p(x) | L2(d) [m] | p(kappa) | kappa err |
|---|---|---|---|---|---|---|---|---|---|
| sphere poly | DAVOF, quadratic faces | 2.14 (2.04) | 2.19 (2.08) | 2.17e-5 | 1.63e-3 | 2.02 | 1.33e-8 | - | - |
| | DAVOF, Detrixhe-Aslam | 1.05 (0.99) | 1.04 (0.99) | 6.46e-4 | 5.07e-3 | 2.01 | 9.22e-9 | - | - |
| | plicRDF | 1.05 (1.04) | 1.06 (1.05) | 2.08e-3 | 1.07e-2 | 2.00 | 6.07e-8 | - | - |
| | gradAlpha / isoAlpha | 0.05 / 0.06 | 0.06 / 0.06 | 4.64e-2 / 5.37e-2 | 0.10 / 0.12 | 1.99 / 1.97 | 6.2e-8 / 6.4e-8 | - | - |
| ellipsoid poly | DAVOF, quadratic faces | 1.97 (1.97) | 1.79 (1.80) | 8.36e-5 | 4.23e-3 | 1.98 | 1.77e-8 | 1.83 | 2.33e-3 |
| | DAVOF, Detrixhe-Aslam | 1.26 (1.17) | 1.06 (1.02) | 8.89e-4 | 9.25e-3 | 2.01 | 1.22e-8 | 1.58 | 2.60e-3 |
| | plicRDF | 1.18 (1.12) | 1.11 (1.10) | 2.73e-3 | 2.02e-2 | 2.00 | 8.03e-8 | - | - |
| | gradAlpha / isoAlpha | 0.12 / 0.09 | 0.13 / 0.09 | 4.55e-2 / 5.36e-2 | 0.12 / 0.12 | 1.99 / 1.98 | 8.1e-8 / 8.3e-8 | - | - |
| plane poly | DAVOF, both states | round-off | round-off | 1.7e-14 / 8.1e-15 | 2.9e-13 / 2.1e-13 | - | 2.5e-17 | - | - |
| | plicRDF / gradAlpha / isoAlpha | - | - | 3.47e-3 / 6.99e-2 / 8.14e-2 | 0.47 / 0.55 / 1.95 | 1.47 / 1.58 / 1.40 | 5.9e-8 / 8.1e-8 / 1.4e-7 | - | - |

Reading. P11 CONFIRMED on all three mesh classes (the plane at 1e-14 on cfMesh polyhedra
for every N; plicRDF's maximum normal error there is 0.47, isoAlpha's 1.95). P12: the
orders hold on polyhedra (quadratic faces 2.14 / 1.97 in L1, Detrixhe-Aslam 1.05 / 1.26,
position 2.0, the ellipsoid curvature 1.83 against 1.55 on hex), and the constants are
SMALLER, not larger: the sphere normal L2 2.17e-5 against 1.19e-4 on hex (5.5x; the Detrixhe-Aslam state 6.46e-4 against 2.43e-3, 3.8x), the
ellipsoid 8.36e-5 against 3.01e-4 (3.6x), the curvature error 2.33e-3 against 4.55e-3
(2x); the polyhedral cells have more faces (the F fractions over-determine the normal)
and the point-neighbour ring is larger (the tensor fit sees more samples), which is the
mechanism the second half of P12 named for the curvature and did not expect for the
normal. P13 CONFIRMED: plicRDF 1.05 to 1.18 on polyhedra, gradAlpha and isoAlpha do not
converge at all there (0.05 to 0.13). Combined tables per surface
(`davof_normal_proposal_<surface>Normal3D_meshes.tex`): the sphere and ellipsoid ones lack
the hexahedral block because the hexahedral studies live on the laptop, not in the cluster
clone; the article stitches the three per-study bodies itself.

Verdicts of 10.16, all seven ladders in: P11 CONFIRMED (hex, perturbed, polyhedral). P12
CONFIRMED for the orders and FALSIFIED in the favourable direction for the constant: the
quadratic-face normal is second order on every mesh class with the hex constant on the
perturbed mesh and a 3.6 to 5.5 times smaller one on cfMesh polyhedra; the Detrixhe-Aslam
normal first order (1.0 to 1.3); the position second order everywhere; the one-ring
curvature 1.5 on hex and perturbed, 1.8 on polyhedra with half the error. P13 CONFIRMED
and sharpened: plicRDF 0.9 to 1.2 on every mesh class with maximum normal errors of 0.14
to 0.47 off the hex mesh, gradAlpha and isoAlpha below 0.3 and non-convergent on
polyhedra. Open: the volume-less wetted cells of the perturbed quadratic-face state (the
N = 20 rung of the sphere) as a case of the proposal's defect rule; the area diagnostic's
tolerance in `Allrun_poly` (Tomislav).

### 10.17 The transported plane and sphere tests through the TwoPhaseFlow fork (2026-10-14)

The first transport measurements of the DAVOF loop, run reproducibly through this workflow
with the fork's `advectorVoF` as the solver (configs `config/davof/planeTranslation3D_hex.yaml`
and `sphereTranslation3D_hex.yaml`, cases `cases/davof/{plane,sphere}Translation3D`, the row
token `TPF_ROW` expanded by `workflow/scripts/tpf_rows.py`, the converter
`tpf_advection_metrics.py`, the report `make_davof_advection_table.py`; fork at 8f098c1 with
the abort guard; cluster jobs 55428354/55 of 2026-10-13 23:08, interrupted by a NODE_FAIL of
node mpzd0001 at 00:31 that the executor does not resubmit, resumed as 55430286/87 with
`slurm-requeue` in the profile and the guard's parser fixed for snakemake 9's indented group
blocks; plane done 09:42, sphere 09:03; 224 and 63 cases, no solver failure; tables in
10f2f438). Seven rows on the identical initial alpha of `initDavofFields`: davof-cellLocal
(the own normal), davof-ringArea (the normal of the ring's union), davof-ringFit (the
weighted quadric through the ring's planes), davof-exact (the exact normal, control),
davof-plicRDFnormal (the loop driven by the plicRDF normal), plicRDF, isoAlpha. Plane:
detrixheAslam state, arms normal (U = n) and oblique (U = (1 0 0)), travel 0.2 L, N = 20 to
160, Co = 0.05 to 0.5, the gate at 1e-12 on the run maxima of L1(alphaf) and L1(n). Sphere:
quadraticFaces state, R = 0.15 in the 2 x 1 x 1 box, 3.3 diameters, N = 32 / 64 / 128,
Co = 0.1 / 0.2 / 0.5; errors at the end, orders least-squares over the three meshes.

Plane gate (32 cases per row: 2 arms x 4 N x 4 Co): davof-exact, davof-ringArea and
davof-ringFit PASS everywhere (max L1(alphaf) 1e-15 to 1e-14, max L1(n) 1e-15 to 1e-13);
davof-cellLocal passes 10 and fails 22: at Co 0.05 to 0.2 it is at round-off for N = 20 and
40 (and marginal, 1e-11 in L1(n), at N = 40 on the normal arm at Co 0.2), drifts to 1e-10
to 1e-7 at N = 80 and breaks at N = 160 (L1(alphaf) 3e-3 to 0.14), and at Co 0.5 it breaks
on the normal arm from N = 40 on (0.2 to 0.36); davof-plicRDFnormal fails everywhere (the
plicRDF normal of a plane is not exact at the corner cells: 1e-4 at N = 20 to 0.17 at
N = 160); plicRDF and isoAlpha carry no face fractions (gate "-"). P-A1 CONFIRMED for the
ring rows and the control; P-A2 FALSIFIED: the own normal does not pass at every N below
Co 0.5, the amplification through sliver cells grows with the number of steps and with
N, not with Co alone; P-A3 CONFIRMED.

Sphere, end of 3.3 diameters, L1(alphaf) as a fraction / L1(n) / L1(x) / Eshape at
N = 32, 64, 128 with the least-squares order (pairwise in parentheses):

| row, Co | L1(alphaf) | L1(n) | L1(x) | Eshape | bounded |
|---|---|---|---|---|---|
| ringArea, 0.1 | 2.33e-2, 1.22e-2, 7.72e-3: 0.80 (0.93 / 0.66) | 4.67e-2, 2.53e-2, 1.26e-2: 0.94 | 1.76e-3, 4.56e-4, 1.34e-4: 1.86 | 4.89e-4, 1.34e-4, 3.90e-5: 1.82 | 1e-13 |
| ringArea, 0.2 | 1.78e-2, 9.45e-3, 5.76e-3: 0.81 (0.92 / 0.71) | 4.06e-2, 2.23e-2, 1.09e-2: 0.95 | 1.35e-3, 3.56e-4, 1.02e-4: 1.86 | 3.85e-4, 1.08e-4, 3.05e-5: 1.83 | 1e-14 |
| ringArea, 0.5 | 1.19e-2, 6.06e-3, 2.96e-3: 1.00 (0.97 / 1.03) | 3.34e-2, 1.84e-2, 8.12e-3: 1.02 | 7.71e-4, 1.89e-4, 4.88e-5: 1.99 | 2.58e-4, 6.27e-5, 1.58e-5: 2.02 | 1e-14 |
| exact, 0.2 | 1.04e-2, 5.68e-3, 3.32e-3: 0.82 | 3.64e-3, 1.72e-3, 4.70e-4: 1.48 | 9.83e-4, 2.53e-4, 6.76e-5: 1.93 | 2.76e-4, 7.16e-5, 1.94e-5: 1.92 | 1e-14 |
| plicRDF, 0.2 | - | 3.89e-2, 1.91e-2, 9.20e-3: 1.04 | 1.53e-3, 3.80e-4, 1.04e-4: 1.94 | 4.32e-4, 1.14e-4, 3.08e-5: 1.91 | 1e-14 |
| isoAlpha, 0.2 | - | 5.66e-2, 3.46e-2, 2.03e-2: 0.74 | 1.60e-3, 4.25e-4, 1.20e-4: 1.87 | 4.67e-4, 1.28e-4, 3.53e-5: 1.86 | 1e-14 |
| plicRDFnormal, 0.2 | 2.07e-2, 1.17e-2, 1.68e-2: 0.15 | 4.03e-2, 3.32e-2, 4.10e-2: -0.01 | 1.58e-3, 4.68e-4, 3.09e-4: 1.18 | 4.31e-4, 1.15e-4, 6.26e-5: 1.39 | 1e-14 |
| ringFit, 0.2 | 7.93e-2, 1.04e-2, 6.02e-3 | 2.05e-1, 1.22e-2, 5.49e-3 | 1.32e-2, 3.83e-4, 1.05e-4 | 1.15e-3, 1.13e-4, 3.10e-5 | 1e-14 |
| cellLocal, every Co | aborted (interface cells above 10 x the initial count) at step 205 / 232 / 177 (Co 0.1), 111 / 117 / 97 (0.2), 46 / 56 / 45 (0.5) | | | | |

Reading and verdicts. P-B1 CONFIRMED: ringArea transports the face fractions at first
order (0.80 / 0.81 / 1.00 over the three Courant numbers, the clipping law's q = 1 for the
plane clipping, which the exact-normal control reproduces at 0.85 / 0.82 / 0.99: the order
is the update's, not the normal's), the normal at first order (0.94 to 1.02) at 1.1 to 1.35
times plicRDF's error (1.09e-2 against 9.20e-3 at N = 128 and Co 0.2; 8.12e-3 against
6.03e-3 at Co 0.5), the centre at second order (1.86 to 1.99) at plicRDF's error to three
digits (1.02e-4 against 1.04e-4), the shape error equal to plicRDF's (3.05e-5 against
3.08e-5), bounded to 1e-14 and conservative to round-off, at 1.16 times plicRDF's cost per
run (3797 against 3283 s at N = 128, Co 0.2, serial). P-B2 CONFIRMED: ringArea, exact,
plicRDF and isoAlpha bounded at every N and Co; cellLocal aborts at every N and Co (the
sliver amplification is a property of the own normal, not of the resolution); ringFit stays
bounded but is erratic (fit fallbacks 0 to 1702 per step; N = 64 at Co 0.5 ends at a normal
error of 0.54): unusable as it stands. P-B3 FALSIFIED in its stated direction: at fixed N
the ringArea errors RISE as Co falls (N = 64: L1(alphaf) 1.22e-2 / 9.45e-3 / 6.06e-3 for
Co 0.1 / 0.2 / 0.5; the normal 2.53e-2 / 2.23e-2 / 1.84e-2), because the per-step clipping
error accumulates over the 1/Co more steps; the exact control is flat in Co (6.57e-3 /
5.68e-3 / 7.27e-3) and plicRDF's own normal error falls with Co too (2.31e-2 / 1.91e-2 /
1.39e-2). NEW and open: the loop driven by the plicRDF normal (davof-plicRDFnormal) does not
converge on the sphere (normal error 4e-2 to 5e-2 at every N and Co, face fractions 1.7e-2
to 2.7e-2 at N = 128), although plicRDF alone converges at 1.04 on the same alpha; the
difference between the two rows is the DAVOF face update's flux rules on the faces beyond
the downwind faces (the second ring and the plane's flux where the upwind cell is a sliver
or is crossed) and the plane positioned by alpha along the helper's normal; which of the
two costs the convergence is not measured. On the plane the same row fails the gate because
the plicRDF normal is not exact at the corner cells. Conclusion for the design: the
one-ring area-weighted normal is the first stable DAVOF transport loop, first order in the
face fractions by the clipping law, equal to plicRDF in position and shape and within 35 %
in the normal, at a cost of 16 %; the paraboloid update is what q = 2 needs.

Appendix to 10.17 (2026-10-14 evening). (a) The perturbed twins (`planeTranslation3D_perturbed`,
`sphereTranslation3D_perturbed`, Co 0.2, 0.1 h vertex displacement; tables c9cc9ba): no row
passes the plane gate (the exact normal keeps the normal at round-off while the transported
scalar face fraction is off by 1 to 5 %), every DAVOF row diverges on the sphere, the
exact-normal control included (abort at step 144 / 217 / 315), the ring normal at 101 / 147 /
181; plicRDF and isoAlpha stay bounded with converging shape errors but non-convergent normals
(plicRDF 4.9e-2 / 5.5e-2 / 8.8e-2, ten times its static values). P-C1 FALSIFIED, P-C2 too
optimistic (divergence, not a floor), P-C3 half; P-D1 and P-D3 CONFIRMED, P-D2 FALSIFIED.
Reading: the scalar face fraction is not transportable on warped faces (the wetted area vector
and its first moment are the face quantities, 10.12), and the face update's plane-based flux
rules lose boundedness there; distorted-mesh transport is open for the framework as a whole.
(b) The plicRDF-normal control of the hex study: the "kept m" explanation (own normal where
the helper has no plane) is withdrawn, those cells receive no plane at all; the ring fallback
(`noNormalFallback ringArea`, fork 1019e09, row `davof-plicRDFnormal-ring`, f4e22198)
reproduces the control to all digits (the `*_hexControl` configs are not submitted). The
repeated-reconstruction hypothesis is falsified by `sphereTranslation3D_hexInterp`
(interpolateNormal false, e3d9ff2): the control's normal error 4.03e-2 / 3.45e-2 / 4.52e-2 at
N = 32 / 64 / 128 (switch on: 4.03e-2 / 3.32e-2 / 4.10e-2), standalone plicRDF within 4 % of
its values. P-G1 FALSIFIED, P-G2 CONFIRMED. The control's non-convergence stays open (the
candidate: the interaction of the plicRDF normal with the face update's flux rules beyond the
downwind faces); the ring-normal result does not depend on it.

Appendix to 10.17 (2026-10-15, cost). The study `sphereTranslation3D_hexTiming` (220c9b0:
rows ringArea, exact, plicRDF; N = 32 / 64 / 128; Co 0.2; the diagnostics at the end only via
the tokens TPF_AFE_WRITE_CONTROL onEnd, TPF_AFE_CONSISTENCY false, TPF_VFE_WRITE_CONTROL onEnd;
profiling on) ran on the fork at 33198a9 (the six bit-identical speedups of the day). P-T2
CONFIRMED: every end value equals the hex study's to all digits in all nine cases. P-T1
CONFIRMED at N = 32 (medians equal, 0.020 s per step) and N = 128 (1.02 against 1.19 s, the
run 799 against 870 s) and contaminated at N = 64 by node contention (the shared code paths
ran twice as slow in the ringArea run). Profile at N = 128: davof::reconstruct 94 s against
plicRDF::reconstruct 265 s; the face update about 45 s more than plicRDF's flux (two full-mesh
default passes 37 s, the second-ring sweep 7 s, the second-layer bounding 10 s). The 1.16x of
the hex study was the diagnostic's cost. The `TPF_PROFILING` block renders inert (false) in
every other study.

Appendix to 10.17 (2026-10-15 12:40, the second speedup round). Fork c1b3995 / 660d741 /
bfe0c00 (the default face fractions without the full-mesh passes, the area normal over the
published wet-face list in global face order, extendMarkedCells over the marked faces).
`planeTranslation3D_hexCheck` (16 cases) and `sphereTranslation3D_hexTiming2` (27 cases, three
repetitions through TPF_REPEAT; tables c104c1f): P-U1 CONFIRMED, every gate maximum and end
value equals the hex study's to all digits, the repetitions identical to each other. P-U2
CONFIRMED: medians of the per-step medians ringArea / plicRDF 0.02 / 0.02 (N = 32), 0.13 /
0.15 (N = 64), 1.31 / 1.72 s (N = 128); the same-node profile at N = 128 davof::reconstruct
106 s (20 % of advect) against plicRDF::reconstruct 398 s (49 %), formAreaNormal 4 s, the
default passes 12 s. Node noise: the same N = 128 case 639 to 1019 s across repetitions;
quote medians only.

Appendix to 10.17 (2026-10-16, the sagitta face clipping). The studies
`sphereTranslation3D_hexSagitta` (21 cases) and `planeTranslation3D_hexSagitta` (32 cases; fork
912fdb3, tables 8d538cc) FALSIFY P-S1, P-S2, P-S4 and confirm P-S5: the reference rows equal the
hex study's to all digits; with the exact normal the sagitta correction leaves L1(alphaf) at
1.028e-2 / 5.528e-3 / 3.233e-3 (N = 32 / 64 / 128, order 0.83) against the plane clipping's
1.039e-2 / 5.683e-3 / 3.322e-3 (0.82), applied on every cut face (860 / 3542 / 13570 per step, no
tensor fallback); the ring-normal loop is unchanged; the own normal with and without the sliver
fallback diverges at every N with thousands of tensor fallbacks per step; the plane passes at
1e-12 for the exact and the ring rows (8 of 8), the sliver variant in 6 of 8. Diagnosis in three
pre-registered steps: the formula with the exact tensor removes 92 to 97 % of the plane's face
error in pure geometry (the proposal folder's `sagitta_design_check.py`); the fitted tensor is
right on every interface cell (`sphereTranslation3D_hexSagittaDiag`, one step at N = 32 with the
tensor fields written, `davof_tensor_check.py`: eigenvalues 1.04 to 1.07 times 1/R, z_0 < 0 of
the expected size, P-D1 confirmed); and the per-face fields of fork 0a5207b (`alphafError`,
`sagittaCorrection`; `davof_sagitta_face_check.py`, `davof_sagitta_cell_check.py`) show the
correction tracking -e_plane in pattern (r = 0.86) with an offset of about -1.1e-2 per face,
the plane row's errors of mixed sign where the geometry gives one sign, and the state's volume
fractions below the exact sphere's on 384 of 412 cut cells (mean -8.5e-3 of V): the
quadraticFaces state carries the Detrixhe-Aslam chord cell fractions, the plane positioned by
them sits inward of the plane the correction's volume anchor z_0 = -(1/2) B:M assumes, and the
correction is off by that shift (the finding of 10.13 in transport: a volume anchor is
consistent only in its own quadrature). The anchor must come from the face fractions. Submitted
2026-10-16: `hexSagittaDiag2` / `hexSagittaDiag3` (one step on both states at Co 0.2 / 0.02,
P-D2, P-D3) and `hexSagittaExactState` (the ladder on the exact-sphere state, P-S7, P-S8), jobs
55459343-45; the article FINDINGS entry of 2026-10-16 has the full account.

Appendix to 10.17 (2026-10-16 01:00, the sagitta diagnostics on both states). `hexSagittaDiag2`
(Co 0.2) and `hexSagittaDiag3` (Co 0.02), one step at N = 32: P-D2 CONFIRMED, on the exact-sphere
state the plane's face error is one-signed (mean 1.28e-2 at Co 0.02, x and y/z faces alike) and
the correction removes 90 to 93 % of it (L1 6.34e-3 -> 4.40e-4 at Co 0.02, 5.79e-3 -> 5.86e-4 at
Co 0.2; r = 0.99, the storage identity to 1e-13); P-D3 CONFIRMED, on the quadraticFaces state
all 416 cut cells lie below the exact volume (mean -8.4e-3 of V), the plane's error is of mixed
sign and less than half the exact state's (the deficit's inward shift cancels most of the secant
overestimate), and the full correction makes it 2.2 times worse (2.78e-3 -> 6.22e-3). The x-face
asymmetry of the first diagnostic is a Courant-dependent sweep effect (gone at Co 0.02), which
the correction removes too. The exact-state ladder (job 55459345) decides P-S7 / P-S8.

Appendix to 10.17 (2026-10-16 11:00, the exact-state sagitta ladder). `sphereTranslation3D_
hexSagittaExactState` (12 cases, job 55459345, tables 0729fd5): P-S7 FALSIFIED, the sagitta row
with the exact normal ends at L1(alphaf) 8.42e-3 / 4.33e-3 / 2.58e-3 (order 0.85) against the
plane clipping's 1.15e-2 / 6.24e-3 / 3.64e-3 (0.83), a constant 27 to 30 % lower at the same first
order, where one step from the same state had removed 93 %; P-S8 nominally confirmed (the ring
loop's L1(n) 1 % lower, its face fractions 13 % lower, bounded, no tensor fallback) and void in
practice. Cause: the transported alpha is driven by the plane's volume fluxes and acquires a
one-signed O(h/R) deviation from the volume under the surface within the steps a cell holds the
interface, after which the volume-anchored offset overshoots as on the chord state. A
second-order face update under transport needs the offset anchored by the face fractions and
volume fluxes of the same clipping surface (Stage 1 of the proposal), not a correction of the
plane loop. The sagitta rows, configs and fields stay as the measured record; `faceClipping
plane` remains the default and the loop of 10.17 the baseline.

### 10.18 The dual advection with one curved surface per cell: the gate exact, the interpolating paraboloid unstable (2026-10-23)

The fork's `clippingSurface exactSphere | paraboloid` (isoAdvectionCurved.C) replaces the plane
sweep by one surface per cell for the volume flux (the liquid volume of the prism between the
pre-image and the face, adaptive Gauss over the face with exact segment fractions) and the face
update (the pre-image's quadratic face cut); the studies run through `cases/davof/*Translation3D`
with the row elements TPF_CLIPPING_SURFACE, TPF_FIT_ANCHOR and TPF_FIT_DOF derived from TPF_ROW
(`tpf_rows.py`), the converter columns of the "curved clipping" Info line and the report's run
maxima. Smokes 1 and 2 (`sphereTranslation3D_hexCurvedSmoke`, `hexCurvedSmoke2`; N = 32, Co 0.2):
the exact-surface gate gives exact face fractions (1e-14) from the second step, bounded and
conservative to round-off, carrying the initial state's volume unchanged (the chord deficit of
quadraticFaces, 8e-3 in L1; on the exact state a 2.4e-4 residual between the initialisation's
box-ball volumes and the ray quadrature, under measurement in `hexCurvedQuad`); the
six-coefficient surface interpolating the cell's face fractions is unstable from the first step on
both initial states and with every normal (first step 20 / 13 / 2.6 times better than the plane
loop's with the exact / own / ring normal, then growth 2 to 2.5 per step, stopped at step 10 to
21): the volume fraction never enters that surface and the faces of a cell come from different
surfaces, so the exact interpolation amplifies. Predictions P-S2b falsified as foreseen in the
reply to the proposal session; the quadrature cost 9 s per step at 1e-10 / 12 and 1.2 s at
1e-8 / 9 against 0.02 s for the plane loop. Running: smoke 3 (`fitAnchor volume`: the offset
shifted so that the volume under the surface is alpha_k V), the quadrature study, smoke 4
(`fitDof none`: the ring-tensor paraboloid positioned by the volume alone, and the one-scalar face
anchoring). The article FINDINGS entry of 2026-10-23 has the account; the ladders
(`hexExactSurface`, `hexParaboloid`, `planeTranslation3D_hexParaboloid`) wait for a stable
variant.

Appendix to 10.18 (2026-10-23 night). Smoke 3 (`hexCurvedSmoke3`, variant 2: the six-coefficient
face fit with the offset shifted afterwards to the volume alpha_k V, slopes and tensor as fitted):
diverges like the unanchored surface on both initial states with every normal. Quadrature study
(`hexCurvedQuad`): the ray quadrature of the cell volume under the exact sphere agrees with the
initialisation to 3e-6 in L1 at every tolerance 1e-8 to 1e-12 and depth 9 to 14 (cost 4.6 to 158 s
per step at N = 32); the 2.4e-4 of smoke 2 was the old tet-quadrature defect routine. Smoke 4
(`hexCurvedSmoke4`, variant 3, rows davof-paraboloidRing-*: the frame from the normal, the tensor
from the ring, zero slopes, the offset from the volume alone; fitDof none, fitAnchor volume):
STABLE with the exact and the ring normal over 160 steps at N = 32, Co 0.2 on the exact state
(bounded 1e-20, conservative, the defect at the anchor's tolerance): exact normal L1(alphaf) 5.7e-4
against the plane loop's 1.15e-2, own-faces normal 2.8e-4 against 3.7e-3, centroid 4.1e-4 against
7.7e-4; ring normal 5.9e-3 against 1.78e-2, normal 3.3e-2 (the ring's bias), centroid 5.3e-4
against 1.35e-3; the own-normal row diverges through the sliver cells (stopped at step 24); the
one-scalar face anchoring (paraboloidOne) drifts and is stopped at step 81. The anchor's Newton
with the flux's quadrature cost 28 s per step; fork 988fdf6 runs it at anchorTolRel 1e-6 with one
Newton and secant steps and the defect every tenth step. Queued: the gate ladder
`hexExactSurface` (55557458, running), smoke 5 (`hexCurvedSmoke5`: the interpolating surface at
Co 0.05 and at the tight quadrature, the proposal session's check), then the rebuild and the
ladders `sphereTranslation3D_hexParaboloid` (variant 3 rows, exact state, 36 cases; P-P1 to P-P5
restated in the config) and `planeTranslation3D_hexParaboloid` (48 cases).

Appendix to 10.18 (2026-10-24 early, smoke 5). The face-interpolating surface's instability is
independent of the quadrature tolerance (identical digits at 1e-8 / 9 and 1e-10 / 12, the same
abort step) and of the Courant number (at Co 0.05 stopped after 0.9 cell lengths against 2.8 at
Co 0.2: the amplification is per application of the map). P-S5a and P-S5b held.

Appendix to 10.18 (2026-10-24 morning). The plane test of variant 3 was not at round-off (1 to 76
anchor failures per case, errors to 4e-4, the own-normal rows diverging): the anchor's
Newton-plus-secant failed where V(c) is flat or the start far. Fork 2eb3496 brackets the root and
solves it by the Illinois rule; the plane test (55561277) and the sphere ladder (55561278) were
restarted on it after the rebuild (55561276); the partial runs are kept aside as
`*.secantAnchor-*`. The gate ladder (55557458) runs on, 14 of 27 cases done after three hours (the
N = 128 cases).

Appendix to 10.18 (2026-10-24 midday, the quadrature). The plane test with the bracketed anchor
stayed at 1e-5 to 4e-4: the adaptive six-point triangle rule's floor along the surface traces
(about 4e-6 of the cell volume at depth 9). Fork bcf7dd3 splits the triangle integrals at the
quadric traces (exact roots along sweep lines, Gauss-12 per piece, adaptive outer rule):
`hexCurvedQuad` rerun gives the defect 7e-9 in L1 at every tolerance (was 3e-6), 16 000 to 27 000
leaves per step (was 6.7 million), 1.3 to 1.6 s per step at N = 32. The plane test (55562658) and
the sphere ladder (55562659) restarted on it; the partial runs kept aside as `*.before-split-*`.

Appendix to 10.18 (2026-10-24 afternoon, the anchor's criterion and the two quadrature layers).
With the split quadrature the plane test stayed at 4e-6 to 1e-5 in L1(alphaf), with identical
digits across the anchor's changes: fork 9515ea9 refines the outer rule at least three levels
before accepting and runs the anchor and the flux at 1e-11 of the volume (templates 7e052553);
fork 5b8c613 closes the anchor on the offset (bracket width 1e-11 h) instead of the volume
residual, which left the offset unbounded in corner cells where V(c) is flat. The one-step
diagnostic `planeTranslation3D_hexCurvedDebug` (fork 90c0918, leia 35bc5190: the cell with the
largest anchor residual before the anchor, its pyramids against the exact linear cut of the
Detrixhe-Aslam tets; N = 40, Co 0.2, exact normal, both arms) located the error in the quadrature
of a cell whose centre lies a few 1e-3 h from the plane: the ray weight of the pyramid rule has a
layer of that width next to the trace (the rays cross the surface at s near zero on one side of
the trace and not at all on the other), the trace was an exact breakpoint, the layer was not, and
the 12-point rule on the piece missed it by up to 6.4e-5 of the cell per triangle (2.1e-4 in the
cell, 18 to 30 leaves). Fork cc8effa refines the inner integral adaptively per piece (1e-12, depth
20): the worst residual 2.1e-4 -> 3.9e-10 of the cell, the anchor's face shift 2e-4 -> 4e-10, the
one-step L1(alphaf) 4e-6 -> 4.1e-12 (normal arm) and 3.6e-12 (oblique), Linf 7e-10, 59 s per step
at N = 40 (the diagnostic's 927 000 leaves). The remaining 3.9e-10 sits in one triangle with 18
leaves: the kink of the outer sweep where the trace crosses the swept edge, which bisection
resolves at second order only; fork 272dc2c splits the outer sweep at those crossings (the same
exact roots along the edge). Build 55566595; the diagnostic (55566596), the plane test (55566597)
and the sphere ladder (55566598) rerun on it; the partial runs kept aside as `*.innerGauss-*`
and `*.innerOnly-*`. The gate ladder (55557458) runs on, 19 of 27 after seven hours.
The diagnostic on 272dc2c (55566596): the worst anchor residual before the anchor 2.0e-14 of
the cell (the quadrature against the exact linear cut: 0.8776474705299629 against
0.8776474705299703), the anchor's face shift max 8.8e-13 mean 2.2e-14, the surface defect max
1.0e-14, the one-step L1(alphaf) 2.4e-15 (normal arm, Linf 2.5e-14) and 6.5e-15 (oblique,
Linf 4.0e-12), 506 000 leaves (was 927 000), 44 and 38 s per step at N = 40 with the
diagnostics on. The split quadrature is at round-off on the plane; the plane test and the
sphere ladder of variant 3 run on this binary.
The cost: 44 s per step at N = 40 extrapolates to days for the N = 160 plane and the N = 128
sphere cases, and the three forced outer levels (9515ea9, against accidental acceptance at the
edge-crossing kinks) are not needed once those crossings are exact breakpoints. Fork 8e76295
forces one level: the diagnostic (55566844) repeats every digit of 272dc2c (worst cell 2.0e-14,
face shift 8.8e-13 / 2.2e-14, L1(alphaf) 2.4e-15 / 6.5e-15) with 128 400 leaves (was 506 000)
and 13.6 / 12.0 s per step at N = 40 with the diagnostics on (was 44 / 38). The first eleven
plane cases of 272dc2c (N = 20) were at round-off (run maxima of L1(alphaf) 1e-14, no anchor
failure) except the own normal on the normal arm at Co 0.5 (3e-8 after eight steps, the sliver
amplification of the plane loop). Build 55566843; the plane test (55566845) and the sphere
ladder (55566846) rerun on 8e76295, the partial runs kept aside as `*.minDepth3-*`. The
remaining cost is the anchor's cell-volume evaluations (24 pyramids per Illinois iteration);
the exact sphere-tetrahedron overlap (Strobl, Formella, Pöschel 2016) for the control and the
Evrard et al. (2023) closed forms for the paraboloid are the next levers if the N = 128 cases
do not fit the wall time.

Appendix to 10.18 (2026-10-24 evening, the run-time selectable volume rules). Tomislav's
decisions: Kromer and Bothe's method as the general volume rule, the exact sphere-tetrahedron
overlap for the sphere translation case, both run-time selectable; every implemented algorithm
cited in full with its DOI (and the journal version of every arXiv preprint). Fork 038f4d7:
`davofVolumeRule` (RTS, key `volumeRule` of the isoAdvector dictionary) with the rules
`splitQuadrature` (the reference, bit-identical to 8e76295), `kromerBothe` (J. Kromer, D.
Bothe, J. Comput. Phys. 475 (2023) 111840, doi:10.1016/j.jcp.2022.111840, arXiv:2111.01073;
the divergence theorem twice, the paraboloid patch as the line integral of H1 = c t1 - (b11
t1^3/6 + b12 t1^2 t2/2 + b22 t1 t2^2/2) along the conic traces, the wetted face areas by
Green's theorem; the conic in principal form (their Eqs. 26-34) only to order the crossings
and split long arcs, every arc integrated as the graph over its chord with the exact root
along the chord's normal and adaptive Gauss-Legendre 12; fallback per polyhedron to the split
quadrature, counted; `davofKromerBothe.H`) and `exactSphere` (S. Strobl, A. Formella, T.
Pöschel, J. Comput. Phys. 311 (2016) 158-172, doi:10.1016/j.jcp.2016.02.003, through their
MIT header-only reference implementation vendored as `src/VoF/davof/thirdParty/overlap`,
Eigen 3 from `EIGEN_INC`, set by `etc/davof-env.sh`; the cell as the centre-face-edge
tetrahedra, the prism as three tetrahedra per fan triangle). The anchor, the flux, the face
update and the defect all go through the selected rule; the Info line gains "; volume rule
<name> fallbacks N" (converter columns volumeRule, nRuleFallback; leia b9a1e935 + 877284a0:
the token TPF_VOLUME_RULE, default splitQuadrature, in both templates and parameter files).
Studies (predictions in the configs): `planeTranslation3D_hexCurvedDebugKB` (the one-step
plane diagnostic with kromerBothe: must reproduce the linear cut to round-off, no fallbacks,
cheaper than 13.6 / 12.0 s per step; 55573931), `sphereTranslation3D_hexRuleKromerBothe`
(variant 3 rows exact / ringArea at N = 32, Co 0.2, split against kromerBothe: P-KB1 agreement
to 1e-9 relative, P-KB2 fallbacks below 1 %, P-KB3 a third of the time; 55573932),
`sphereTranslation3D_hexRuleExactSphere` (the control, split against exactSphere: P-ES1 the
gate holds, P-ES2 the defect agrees, P-ES3 a third of the time; 55573933), all behind the
cluster build 55573930 (Eigen 3.4.0 downloaded to ~/OpenFOAM/thirdparty/eigen on the cluster
as on the laptop). The serial ladders (55557458, 55566845, 55566846) run on unchanged. The
MPI exchange of the surfaces follows the rules.
Results so far. `hexRuleExactSphere` (55573933): the exactSphere rule holds the gate with the
split quadrature's digits (L1(alphaf) 1.59e-10 at step 1, 2.7e-16 at step 2, 1.9e-14 at the
end, Eshape 3.5803e-4 in both, min alpha -9.7e-13 against -2.8e-13), P-ES1 and P-ES3 held
(0.090 s per step against 0.450, diagnostics included); P-ES2 falsified in the informative
direction: the surface defect is 6.6e-12 in L1 with the exact overlap against 2.7e-6 with the
split quadrature at the templates' depth 9, so the initialisation's cell fractions agree with
the exact sphere to 7e-12 and the residuals attributed to it (3e-6, 7e-9) were the ray
quadrature's floor at finite depth (preprint sentence corrected). `hexCurvedDebugKB`
(55573931): on the normal arm the kromerBothe rule reproduces the exact linear cut
(0.4575043255096169 against 0.4575043255096162), no fallbacks, L1(alphaf) 2.1e-15, 24 656
leaves and 0.25 s per step (was 13.6); on the oblique arm (U along x) the first step's fluxes
are wrong (alpha to 2.8 after bounding, 184 cells corrected): a prism configuration of the
axis-aligned displacement is mis-integrated; fork 955362b adds the per-prism comparison with
the split quadrature (key volumeRuleDebug) and the diagnostic reruns (55574157 after build
55574156). Fork 955362b also carries the parallel curved flux (the surfaces exchanged per
processor face, one side computes, syncProcPatches copies, counters reduced); its check is
`sphereTranslation3D_hexParallel` (np 4, both smokes' rows and rules, P-PAR1 round-off
agreement with the serial references).
The oblique-arm defect located and fixed (fork 1bb9724): the twelve mismatching prisms of the
diagnostic were all faces normal to y or z under the displacement (2e-5, 0, 0), whose d.Sf was
8.7e-28 instead of zero and passed the exact-zero test; the flat prism built for them is not a
closed polyhedron (top and bottom offset within one plane, the sides degenerate), so the face
sum of the divergence theorem does not cancel and gave up to 5e-14, five percent of the cell
volume, per face, while the split quadrature multiplied the same prisms by their zero
thickness. The sweep test is relative now (1e-14 of |d||Sf|) in the caller and in the rule,
and the rule falls back to the split quadrature, counted, where its volume leaves [0, V] by
more than 1e-9. Diagnostic on 1bb9724 (55574588): the normal arm unchanged (L1(alphaf)
2.1e-15), the oblique arm at round-off (L1(alphaf) 4.5e-15, Linf 2.7e-12, min alpha -8e-16,
7 fallbacks of about 6000 faces). The Kromer-Bothe sphere smoke reruns on 1bb9724 (55574594;
the run on 038f4d7 kept aside as `*.degenerate-*`). One more prism slipped through the 1e-14
sweep test (d.Sf 3e-27 of a face normal to the displacement: the zero components of the face
area vectors carry 1e-14 of round-off); the threshold is 1e-10 of |d||Sf| (fork 31028dd).
The four-processor check (`hexParallel`, 55574223, on 955362b): every case ran to the end,
conservative and bounded, the control's gate held (end L1(alphaf) 2.0e-14, 0.11 s per step
with the exactSphere rule), but the agreement with the serial runs was not at round-off: the
control's first step 2.33e-10 against 1.59e-10 and the variant 3 row 7 % off in L1 at step
64, with the norm over 1708 faces against 1668 and 3332 faces set against 3240. Mechanism:
under the translation along x the faces normal to y and z have phi = 0 exactly; both sides of
a processor patch saw themselves upwind by phiF >= 0, both computed the face fraction with
their own surface, the sync swapped the two values, and areaFractionError counted the
processor faces on both sides. Fork 164691d: a zero-flux processor face goes to the owning
side of the patch (the serial rule gives it to the owner cell), and the norm counts a
processor face once, on the owning side, adjacent when either cell is an interface cell (the
remote flag swapped across the patch). The check reruns (55576072 after build 55576071; the
first run kept aside as `*.zeroFluxTie-*`).
`hexRuleKromerBothe` (55574594, on 1bb9724): variant 3 at N = 32, Co 0.2 on the exact state,
160 steps, split quadrature (depth 9) against kromerBothe. Exact normal: L1(alphaf) at step 1
identical (2.581510e-4), at the end 5.173e-4 against 5.205e-4 (0.6 %), Eshape 3.5803e-4 to
twelve digits in both, the surface defect 5e-14 and 4e-14 (each rule consistent with its own
anchor), min alpha -5.6e-14 and -8.7e-19, no fallbacks, no linear cuts; ring normal: 5.6451e-3
against 5.6472e-3 (0.04 %). Cost: 37.6 s per step (split, 5334 s total) against 0.21 s
(kromerBothe, 36 s total), a factor 180; the split quadrature's cost on the sphere is the
adaptive resolution of the silhouette kinks of the pyramids, which the face-based rule does
not meet. P-KB1 (agreement to 1e-9) falsified, P-KB2 and P-KB3 held; whether the 0.6 % is the
split quadrature's depth-9 floor (2.7e-6 of V on the control) is decided by
`hexRuleDebug` (55577404, one step with the debug statistics at depth 9 and 14; prediction:
the difference falls with the depth). Decision: the split-rule sphere ladder (55566846; 8 of
36 cases after 4.3 h, the N = 128 cases unable to finish in 48 h at 600 s per step) cancelled
and kept aside as `*.splitRule-*`; the variant 3 ladder resubmitted with the kromerBothe rule
(55577406 after build 55577403; the N = 128, Co 0.1 case estimated at about 1.2 h).
`hexRuleDebug` first run (55577404): the split quadrature gives the same prism volumes at depth
14 as at depth 9 (converged by its own test) and the kromerBothe rule differs from it by up to
5.3e-6 of the prism volume (mean 1.3e-7) on a few x-normal prisms while the cell volumes
agree to 6.3e-9 (mean 2.5e-11); since the split quadrature and the exact overlap agreed to
6e-12 in the shape error over 160 steps of the control, the split quadrature is the reference
and the Kromer-Bothe PRISM is off at the 5e-6 level (the cells are right): the prediction of
the config falsified. Second run (55578132, fork 93149e3): the knobs kromerBotheMinSplit 1 | 8
(the arc quadrature's convergence) and kromerBotheFan true | false (planar faces whole) to
locate the defect in the side parallelograms or the arcs.
Parallel, second run (55576072, fork 164691d): the control agrees with the serial run at
round-off (the first step's 1.588e-10 identical, the end 1.928907e-14 against 1.928920e-14,
the largest relative difference 9e-4 at a step with L1 1.9e-15); the variant 3 row 0.3 % apart
at worst (step 98) and 6e-5 at the end, with two faces fewer processed: the second ring of a
processor face saw only the local cell's face neighbours, and a face between a remote
candidate and a local non-candidate was processed by nobody. Fork 91729ae exchanges per patch
face the candidate flag and the closest surface of the cell's face neighbours with its
distance, so the serial rule holds across the boundary; third run 55578374 after build
55578373 (the second kept aside as `*.localRing-*`).
Third run (55578374): the face counts match the serial run (3240 at the first step; 9 steps
with a different count), the control at round-off, but the variant 3 row still 0.29 % apart
at step 98; the per-step comparison shows the two runs apart by 2.8e-10 at step 1 (round-off:
the anchor's face shift differs in the eleventh digit), 1.1e-6 at step 2, 6.7e-5 at step 3 and
1e-5 to 2e-3 from then on without further growth, the interface-cell counts identical at
every step. A round-off difference amplified a million-fold in two steps and then bounded is
a discontinuous choice in the loop: the second-ring search takes the closest surface among
the face neighbours by `d < dMin`, and on the Cartesian mesh with the sphere centred on a
cell corner two surfaces are often at the same distance, so the first considered wins (the
local one in parallel, the lower cell index in serial); a different surface changes that
face's fraction by O(h^2). Fork c7b2058 breaks ties within 1e-9 by the surface origin
lexicographically (order-independent); the serial smoke (55582019) and the parallel check
(55582020) rerun on the same binary after build 55582018 (the earlier runs kept aside as
`*.ringTies-*`). The Kromer-Bothe ladder (55577406) started on 93149e3 and its later cases
load whatever library the cluster build last wrote (the version stamp per case records
which); it is a dry run and restarts on the final binary once the tie-break and the
brute-force comparison are in.
Fourth parallel run (55582020, the tie-break): no digit changed in either the parallel or the
serial run (both reproduce their earlier values exactly), so ties never occurred and the
difference is deterministic. Its signature, 2.8e-10 at the first step when every surface
still coincides with the sphere and 1e-4 to 1e-3 from the third, is a different surface on
some faces: the zero-flux faces (two thirds of the processor faces under the translation
along x) took the face-owner cell in serial, the lower cell index of the undecomposed mesh,
and the patch-owning side in parallel, an unrelated rank. Fork 2727892: in serial and in
parallel a zero-flux face goes to the cell with the lexicographically smaller centre (the
neighbour centres across a patch from patchNeighbourField); the serial smoke (55582697) and
the parallel check (55582698) rerun on this binary after build 55582696 (the earlier runs kept
aside as `*.tieBreak-*`). Fifth run: the convention changed the parallel result (end 5.17254e-4
against 5.17280e-4) and halved the worst difference (8.7e-4 at step 52) but left the pattern:
2.8e-10 at step 1 with the anchor's face shift in its eleventh digit, 1e-4 by step 3, 1e-5 to
1e-4 afterwards; the step-1 difference is the same as in every run, consistent with the ring
tensor's summation order across processors (a round-off difference that cannot be removed).
Whether the loop amplifies round-off a million-fold within two steps is decided in serial
by `sphereTranslation3D_hexNoise` (two serial runs differing only in curvedTolRel by 1 %;
prediction in the config): if they differ by 1e-5 or more from step 3, the loop has a
round-off noise floor and the parallel check is as good as serial reproducibility allows.

Appendix to 10.18 (2026-10-24 night, the harvest after a day's suspension of the session;
all SSH watches had died). (1) The split quadrature's defect, settled: in the self-test boxes
that contain the silhouette its integrand is right (the segment fraction against the exact
liquid length on a 60 x 60 grid agrees to the reference's own band width, 5e-7 of dx) but its
value depends on the sweep direction and the triangulation (fan 0.4904706, sweep reversed
0.4905202, the face's two diagonal triangles 0.4904690, exact 0.4905208 of the prism) and
it accepts with 38 to 60 leaves where the Kromer-Bothe rule, exact in all 72 boxes and
matching the exact per-face references to 1e-16, is the reference for the paraboloid; the
split quadrature stays the reference on the plane (no silhouette) and is cross-checked by
the exact overlap on the sphere (9e-12). Its early acceptance at the silhouette is open;
it is the fallback only. (2) The round-off test's first form was void: curvedTolRel 1e-11
against 9.9e-12 left the two serial runs bitwise identical (L1 3.710221e-4 both after 16
steps); the second form perturbs the sphere centre by 1e-13 (55608426). (3) The
Kromer-Bothe ladder dry run (32 of 36 cases, cancelled; kept aside as `*.dryRun-*`): the
exact and ring rows without fallbacks at 0.2 s per step (N = 32), 0.9 s (N = 64) and 5 to
8 s (N = 128; the N = 128, Co 0.1 cases in three hours); the own-normal row aborted at every
N and Co (steps 62 to 331) and the row with the ring normal below 0.1 h^2 at Co 0.5 (steps
57, 67, 79; at Co 0.1 and 0.2 it completed at N = 32 and 64), both with hundreds to 26 000
fallbacks per step to the split quadrature on their wild surfaces and single steps of up to
27 000 s; fork 7845dd2 caps the fallback (1e-8, depth 6, clamped to [0, V]). (4) The
split-rule control ladder (55557458): 26 of 27, the N = 128, Co 0.1 own-normal case stopped
by the 24 h limit at 121 s per step; the split-rule plane test (55566845): 44 of 48, the
exact and ring rows at round-off everywhere (run maxima of L1(alphaf) below 8e-14), the
own-normal row drifting (3e-8 to 8.5e-2 with 23 anchor failures), its four unfinished cases
own-normal cases at N = 80 and 160 stalled by the split quadrature on wild surfaces (kept
aside as `*.splitRule-*`). (5) The three ladders restarted on the final binary 7845dd2
after build 55608481: the plane test (55608482) and the sphere ladder (55608483) with
kromerBothe, the control (55608484) with exactSphere; configs' headers carry the split-rule
runs' results and the dry run's. (6) The round-off floor of the loop (`hexNoise`, 55608426,
two serial runs with the sphere centre 1e-13 apart, kromerBothe, 16 steps): 1.3e-10 of
L1(alphaf) at step 1, 2.0e-6 at step 2, 1.5e-6 at step 3, then 1e-9 to 2.4e-7, interface
counts identical. So the loop does amplify round-off a million-fold within two steps, to a
floor of about 1e-6, and the parallel runs' 1e-4 to 1e-3 were a hundred times above it: a
systematic difference remained, and the split quadrature's silhouette value, which depends
on the sweep direction, supplies it, because the fan of a processor face is traversed in
reversed order on the non-owner side; the two sides integrate the same prism to values
5e-6 of its volume apart. Sixth parallel run (55608669; a first submission, 55608551, went
out with the old config before the push and was cancelled) with the kromerBothe rule against
the ladder's serial case 00004 on the same binary; prediction P-PAR3 in the config: agreement
within the serial floor. RESULT (55608669): against the serial run of `hexNoise` (16 steps,
the same binary up to the fallback cap, which no fallback-free run touches) the four-processor
run differs by 2.8e-10 at step 1 and at most 3.6e-7 (step 4) afterwards, no step above 1e-6,
interface counts and mass identical, min alpha -1.8e-18: P-PAR3 holds, the parallel curved
flux agrees with serial within the loop's round-off floor (the 160-step comparison with the
ladder's case 00004 follows when that case is in). 160 steps in 16 s on four processors
(0.18 s per step, against 0.3 s serial: the N = 32 case is too small to scale). The sensitivity itself is a
property of the second-ring rule worth stating: runs that differ in round-off differ by up to
1e-3 of L1(alphaf), below the ladder's resolution of orders but above a bitwise comparison.
The `hexRuleDebug` second run (55578132): neither kromerBotheMinSplit 8 nor kromerBotheFan
false moves a digit of the 5.316e-6 prism maximum, so the arcs and the fan are
self-consistent; third run (55582004, fork 685a6e9) prints a brute-force midpoint-grid value
(n = 200, 400, 800, Richardson) of every printed mismatching prism as the third opinion:
inconclusive, the midpoint grid scatters at 1e-5 of the prism volume on an integrand with
kinks. Fourth run (55582053, fork 3503079): the split quadrature against the exact
sphere-tetrahedron overlap on the same prisms with the exact SPHERE as the surface agrees to
9.4e-12 of the prism volume over 556 prisms and to 3.0e-12 over the 36 equatorial ones (the
faces at the sphere's equatorial plane, where the segments are tangent to the surface and the
split rule's integrand has the silhouette kink it does not mark): the split quadrature is the
reference even there, and the Kromer-Bothe PRISM is off by up to 5.3e-6, one-signed, at the
equatorial prisms (the mismatching prisms are all at x = 0.5 under the displacement along
x). Fifth run (55582136, fork 2c8daf5): a self-test at the rule's construction against an
exact reference, thin axis-aligned boxes (dx = 0.2 h, h = 1/32) under the axis-aligned
paraboloid z = -7 (x^2 + y^2)/2, the equatorial configuration exactly (the traces on the
side faces are parabolas, on the top and bottom circles), the reference exact in x and
adaptive in y, with the per-face breakdown of every failing box. RESULT: the Kromer-Bothe
rule is exact (1e-15 of the prism) on all 36 boxes, the split quadrature at the run's
settings off by 2.3e-8 (box 15) and 8.7e-9 (box 32), i.e. an unmarked kink of the split
rule's sweep (a trace tangent to a sweep line) accepted early. The axis-aligned
configuration does not reproduce the 5e-6; the sixth run (55582183, fork 33836c8) adds a
tilted frame (e3 with 0.1 and 0.3 rad along the thin direction, 0.33 rad about it) and
b12 = 0.5, against the split quadrature at 1e-15 of the prism and depth 20. RESULT: all
cases with a tilt along the thin direction (a = 0.1, 0.3) agree to 5e-15, as do the a = 0
cases whose box sits at t1 = 0.02; the six cases with a = 0 AND the box straddling t1 = 0
(whichever b and zShift) have the Kromer-Bothe rule off by 3.2e-5 to 5.1e-5 of the prism
while the split quadrature agrees with its tightened version to 1e-10. So the defect needs
the frame axis exactly perpendicular to the thin direction and the box across the surface
origin's plane t1 = 0 (with b12 = 0.5 the parabolic traces on the faces normal to y then
carry their vertex inside the box). Seventh run (55582761, fork afc0008): exact per-face
references for the axis-aligned frame (the wetted area of every box face and the patch
integral int_S g, exact in one direction and adaptive in the other) next to the rule's
per-face terms, to name the defective term. RESULT (reversing the attribution): in the three
failing boxes the Kromer-Bothe wetted areas of all six faces, the face sum, the line integral
and the volume agree with the exact references to 1e-16 (e.g. V 2.99390157064e-06 both)
while the split quadrature, at the run's settings and tightened alike, is 5e-5 of the prism
below (2.99359526506e-06): the SPLIT quadrature is the one off, consistently, in prisms that
contain the silhouette of the surface along the sweep direction (the line where the segments
are tangent to the surface, at x = -b12 y/b11 inside the box only when the box straddles
t1 = 0), which it leaves unmarked; on the sphere control that line lay exactly on the
equatorial face plane x = 0.5 and the trace breakpoints caught it, hence the 9e-12 agreement
with the exact overlap there, while the paraboloid's frame at the PLIC centroid puts it
slightly inside the equatorial prisms. A value the split rule converges to with every
tightening is not a resolution problem; the integrand or the breakpoints are wrong near the
silhouette (under inspection in davofCurvedGeometry.H).

Appendix to 10.18 (2026-10-25, the control ladder with the exact overlap, 55608484, 70 min,
tables committed as leia 9e262ff). 24 of 27 cases complete: P-X1 HELD at every N and Co: the
face fractions exact to 1e-14 from the second step (L1 at the end 1.9e-14, 9.1e-14, 3.3e-13
for N = 32, 64, 128 at Co 0.2; the run maxima 1.6e-10, 1.9e-11, 2.4e-12 are the first-step
linear cuts of grazing faces, falling at third order, and absent at Co 0.5 where no face
grazes), bounded (min alpha -1e-12 to -1e-14), conservative, the surface defect against the
initialisation 1e-11 to 2.5e-10 in L1; the centroid second order (1.98) and the shape error
second order (1.98) identically in all three rows, as a flux that does not use the normal
must give. P-X2 FALSIFIED: the own-faces normal on the exactly transported fractions is
9.5e-4, 3.3e-4, 1.4e-4 at N = 32, 64, 128 (Co 0.2), least-squares order 1.40 (pairwise 1.53,
1.27), against the static recovery's second order on the same state (4.2e-5 at N = 160,
2.05, davof_normal_orders_sphereNormal3D.csv); the static value extrapolated to N = 32
agrees (1.1e-3) and the transported one is 2.1 times larger at N = 128. The two measures
differ in the cell set (the transport's wisp threshold) and possibly in the reference point
of the normal (the exact row's own L1(n) is 8e-5 to 1e-5, not zero); which explains the
factor is open and must be settled before the number enters the preprint. The ring normal
is first order as predicted (3.0e-2, 1.7e-2, 7.7e-3; 0.98). Cost 0.1 s per step at N = 32,
0.5 to 0.8 s at N = 64, 3 to 6 s at N = 128 (the split quadrature: 121 s). The three N = 128,
Co 0.1 cases aborted at step 29 on the overlap library's debug-only assertion ("negative
volume detected": a thin prism tetrahedron's overlap beyond [0, V] by more than sqrt(eps)
V); fork 3f32320 disables the assertions as a release build does and clamps the tetrahedron
overlap to [0, V]; the three cases rerun (55610001 after build 55610000; the aborted ones
kept aside as `*.overlapAssert-*`). RERUN DONE (2026-10-25 01:50, tables leia e65a334f, 27 of
27): the N = 128, Co 0.1 cases complete over 1280 steps at L1(alphaf) 7.8e-13 at the end
(run maximum 2.4e-12), min alpha -9.8e-15, defect 9.3e-10, 5.6 s per step; every number of
the other rows unchanged. The report script skips cases moved aside with a suffix inside a
study folder (leia c5653553), which the first re-score had counted as crashed duplicates.

Appendix to 10.18 (2026-10-25 01:10, the paraboloid sphere ladder, 55608483, PRELIMINARY: 23
of 36 cases; the ring and exact rows complete except exact N = 128 at Co 0.1; tables leia
61d6b0c8). The milestone's question is answered: the ring-tensor paraboloid positioned by
alpha_k with the EXACT normal transports the face fractions at order 2.84 per period at Co 0.2
(5.20e-4, 6.67e-5, 1.02e-5 over N = 32/64/128; Co 0.5: 2.87; Co 0.1 over two meshes: 2.97) and
2.76 after the first step (P-P1 HELD, P-P2 HELD above the hypothesis q >= 1.8): the per-step
error does not accumulate (end-to-first-step ratio 2.0, 1.9, 1.8 over N, stationary from the
second step at Co 0.5). Own-faces normal 2.05 (P-P5 HELD), centroid 1.99, shape 1.98, bounded
to -1.1e-13, no anchor failures, no rule fallbacks, anchor shift 1.75e-3 -> 4.0e-5 over N.
With the RING normal: 1.35 per period (5.65e-3, 2.12e-3, 8.71e-4; Co 0.1: 1.38, Co 0.5: 1.29),
1.05 after the first step (the ring normal's tilt sets the per-step order); normal 1.00 at
8.19e-3 (plicRDF 9.20e-3, plane loop 1.09e-2 at N = 128); centroid 2.09 at 2.85e-5 (plicRDF
1.04e-4); shape 2.30e-5 (3.08e-5); bounded to -1.5e-13; anchor failures at most 8 cells per
step (0.45 % at N = 64, Co 0.1: P-P4's 0.1 % clause FALSIFIED for the ring row, HELD for the
exact row); fallbacks to the split rule at most 3 per step. Courant dependence reversed
against the plane loop: at N = 128 the ring row's error rises with Co (8.52e-4, 8.71e-4,
1.22e-3), the exact row is flat. Own normal (P-P3): N = 32 at Co 0.1/0.2 completes without
tripping the guard but at L1(alphaf) 6e-2, normal error 1.15 and min alpha -0.1; Co 0.5
diverges at steps 62/64; the row with the ring normal below 0.1 h^2 at N = 32 completes at
5e-2 with normal errors 0.4 to 0.5; the remaining own-normal cases run with thousands of
capped fallbacks per step. Cost with the per-step diagnostics on: ring row 0.34 / 1.6 / 8.6 s
per step at N = 32 / 64 / 128, Co 0.2. Follow-ups submitted (leia 1e8b89d0, predictions in
the headers): the ring-QUADRIC-normal row of the ladder (sphereTranslation3D_hexParaboloidRingFit,
55611739; the proposal session's question whether a second-order computed normal keeps the
loop near q = 2) and the paraboloid timing study against plicRDF with the diagnostics off
(sphereTranslation3D_hexParaboloidTiming, 55611741). Preliminary scores sent to the proposal
session with the caveat that nothing enters the proposal before the final tables.

Appendix to 10.18 (2026-10-25 02:20, the paraboloid timing study, 55611741, 59 min, tables
leia on origin after the score): the ring-area row of the paraboloid loop against plicRDF on
the exact-sphere state, N = 32/64/128, Co 0.2, three repetitions, diagnostics at the end
only. Per-step medians (cpu/nSteps, includes the start-up): paraboloid 0.126 / 0.546 / 4.47 s,
plicRDF 0.021 / 0.155 / 1.47 s; factors 6.0 / 3.5 / 3.0. Step times of the logs at N = 128:
4.2 against 1.23 s (3.4). The prediction (2 to 6 times plicRDF at N = 128) HELD at 3.0; the
spread on the shared nodes is a factor 1.5 to 1.8 between repetitions (ring 2.43 to 4.47 s,
plicRDF 1.01 to 1.49 s at N = 128). The errors of the three repetitions reproduce the ladder's
digits exactly (L1(alphaf) 5.647061e-3 at N = 32 in all three and in the ladder). The plane
loop's earlier medians (hexTiming2: ring 1.69 s, plicRDF 1.96 s at N = 128) were on other
nodes; the comparison holds within a study only.

Appendix to 10.18 (2026-10-25 02:50, the ring-QUADRIC-normal row, 55611739, 9 of 9, tables on
origin after the score). All four pre-registered predictions HELD: P-RF1 bounded at every N
and Co (min alpha -2e-14); P-RF2 the transported normal 1.34e-3, 3.61e-4, 1.63e-4 at Co 0.2
(ring-area row 8.19e-3 at N = 128, plicRDF 9.2e-3: fifty times below), order 1.52 (pairwise
1.89 / 1.15, at the margin of the 1.5 clause and falling on the finer pair); P-RF3 q per
period 2.11 at Co 0.2 (8.75e-4, 1.37e-4, 4.73e-5; pairwise 2.68 / 1.53), 2.24 at Co 0.1, 1.93
at Co 0.5, first-step order 1.92 (end-to-first ratio 3.2, 2.8, 2.5 over N); P-RF4 no fit
fallback, no ring fallback, no anchor failure in any step of any case, rule fallbacks at most
2 per step. Centroid 2.00 at 2.57e-5 (the control's value), shape 1.98 at the control's
value, cost 4.5 s per step at N = 128 with the diagnostics on, as the ring-area row. Reading:
the closed loop's order is the order of the normal fed back; the ring quadric's normal is a
computed normal that keeps q at 2 over a period, against 1.35 with the ring-area normal and
2.84 with the exact one, and the pairwise orders say its own order falls toward one on the
finer pair, which the face fractions follow. The erratic behaviour of this normal in the
plane-clipped loop (fit failures in 1700 cells per step) is absent under the surface's own
flux. Preprint: paragraph written into Section 7 (untracked); the row is in CROWS of
make_table_advection.py. Sent to the proposal session.

Appendix to 10.18 (2026-10-25 03:40, the plane test with the kromerBothe rule, 55608482,
PARTIAL: 35 of 48 cases, N = 160 and the own-normal rows still running; gate table on origin,
e262e1a4). The gate is NOT held by the Kromer-Bothe run: the exact row passes at N = 20 and
40 (1e-14), reaches 1.4e-12 at N = 80 (Co 0.2) and 1.6e-6 at N = 160 (Co 0.5, run maximum;
3e-12 at the end); the ring row passes at N = 20 on both arms and at N = 40 on the normal arm
(2e-14), is at 1e-12 at N = 80 on the normal arm and at 6e-9 (N = 40), 7.5e-6 and 1.1e-5
(N = 80), 3.7e-6 (N = 160) on the oblique arm. The cause is in the counters: the rule falls
back to the capped split quadrature (1e-8, depth 6) in thousands of prisms per step on the
finer meshes (exact row: 0 at N = 20, 1 to 16 at N = 40, 1303 to 4010 at N = 80, 43779 at
N = 160; ring row oblique: 250 at N = 20, 2614 and 3307 at N = 40, 12365 and 13849 at N = 80,
69834 at N = 160), and the capped fallback is not exact on a plane while the uncapped split
rule is (the split-rule run of the same study, kept aside as *.splitRule-*, held the gate at
1e-14 up to N = 160 for both rows). On a plane the ring tensor is zero to round-off, every
trace is a line, and the rule should reduce to the exact polygon cut; which check rejects the
flat traces is not visible from the counters. Fork 01ee5d4 counts the fallbacks by reason
(conic, point, root, gradient, area, volume) and prints them per step; diagnostic study
planeTranslation3D_hexParaboloidDiagKB (55614649 after build 55614648; four cases, N = 80,
Co 0.5, exact and ring rows, both arms; prediction in the header: the chord-normal root).
Consequence for the records: the plane gate statement of the preprint and the proposal rests
on the split-rule run (valid: the same loop with the other volume rule, which is exact on
planes); the Kromer-Bothe run's plane table enters only after the fix, as the rule's own
test. The sphere ladder is unaffected (fallbacks 0 to 3 per step there) but the own-normal
rows' thousands of fallbacks per step may share the cause. Proposal session informed.
DIAGNOSED (55614649, 55615537; fork 01ee5d4, a9e20cc): every fallback is the AREA check, and
the printed faces (N = 40, Co 0.5: 326 and 2647 per step for the ring row, 5 for the exact
row's normal arm, none on its oblique arm) carry a LINE trace (33 to 40 of 40 printed; the
rest parabolas of round-off curvature): the ring tensor on the plane is 1e-10, the scaled
curvature 4e-19 against 1e-5 is classified as a line, and the chord-normal quadratic of the
arc integrand has a leading coefficient at round-off with one root on the chord and the
other at 1e9 face sizes; the root was chosen by the side of the chord the arc's midpoint lies
on, which for a flat arc is round-off, so the far root was taken in about half the faces and
the arc area came out at 1e11 to 1e13 times the face area, caught by the range check and
sent to the capped fallback. The first diagnostic's prediction (the root's discriminant)
named the right place for the wrong reason; the second's (an eccentric conic) was wrong: the
kind is the line. FIX (fork 748e401): the side decides only when the midpoint is off the
chord by more than 1e-9 of the face size, otherwise the root nearest to the chord is taken;
curved arcs (sagitta Ls^2/R, far above the threshold) are untouched, so the sphere ladder's
exact, ring and ring-quadric rows stand (their 0 to 3 fallbacks per step were short arc
pieces of the same kind, now resolved). RERUN: the Kromer-Bothe plane run and the second
diagnostic moved aside as *.kbArea-*, both resubmitted after the cluster build (the
diagnostic as the verification: zero fallbacks expected; the plane test with the prediction
of the split-rule run: the exact and ring rows at round-off at every N and Co). VERIFIED
(55616010, fork 748e401): zero fallbacks in all four cases over the whole run, the
quadrature leaves per step down from 1.6e6 to 2.4e7 to 2e3 to 6e3 (the capped fallbacks were
the cost); the gate at 6.9e-14 to 1.0e-13 on all four cases (the ring row's oblique arm was
3.1e-12 before). The plane test rerun is 55616011.

Appendix to 10.18 (2026-10-25 06:50, the cluster cleared and the sphere studies rerun on the
fix). The job list showed 96 own jobs: 71 solve jobs of the drivers of 2026-10-23 (the
split-rule and dry-run sphere ladders, studies long moved aside) still running after two
days (own-normal cases of the uncapped-fallback era, single steps costing hours), the old
Kromer-Bothe plane driver 55608482 with 11 children (its study moved aside at 05:13, so the
driver could not have collected them), and the sphere driver 55608483 with its five own-normal
children at t = 0.08 to 0.27 after 8.7 h (binary before the flat-arc fix). All 90 cancelled
by explicit job id (`cancel_stale.sh`; the watch scripts list drivers by name only, which is
how the 10-23 solves went unnoticed: a lesson, the watch now lists every own job). The sphere
ladder and the ring-quadric study moved aside as *.preFlatArc-* (their exact, ring and
ring-quadric rows are the scored ones above, 0 to 3 fallbacks per step) and both resubmitted
on the fixed binary (fork 748e401) so that the whole ladder is on one binary, with the
own-normal rows expected to run faster without the flat-arc fallbacks; the predictions of
the configs are unchanged, and the scored rows are predicted to reproduce their digits up to
the 0 to 3 resolved prisms per step.

Appendix to 10.18 (2026-10-25 07:30, the plane rerun on the flat-arc fix, 55616011, 43 of 48
at the score, own-normal N >= 80 still running; tables on origin). The exact row holds the
gate at every N, Co and arm (16 of 16; 7e-15 to 2.5e-13, the run maxima growing with N as
round-off does); the ring row holds it in 14 of 16 (8e-15 to 2.5e-13) and misses it on the
oblique arm at N = 160 (7.6e-10 at Co 0.2, 3.1e-10 at Co 0.5, against 3.7e-6 before the fix
and below 8e-14 in the split-rule run): there 3 to 4 faces per step in 4 to 5 steps fall
back on a DEGENERATE conic, the remaining classification gap: the ring tensor's round-off
curvature rotates the principal axes of a line trace onto its direction, one curvature and
the slope along that axis are both at round-off, and the conic B u^2 + 2 l u + K = 0 in the
other coordinate is two parallel lines, the second at -2 l/B, 1e10 face sizes away. Fork
9dbca10 classifies this as the near line when |B K| <= 1e-6 l^2 (the root in the stable
form); third diagnostic planeTranslation3D_hexParaboloidDiagKB3 (the two failing cases;
prediction: zero fallbacks, gate at 1e-13). The own normal fails the gate at every N but
N = 20, Co 0.2 (as in the split-rule run: 3e-8 at N = 20, Co 0.5, 8e-2 at N = 40, Co 0.5,
6e-3 at N = 80 on the oblique arm), with fallbacks of every reason on its wild surfaces.
Once the third diagnostic verifies, the final binary is 9dbca10 and the plane test, the
sphere ladder and the ring-quadric study are rerun on it once more, so that every table of
the preprint stands on one binary. VERIFIED (55618726, 2026-10-25 08:30): zero fallbacks in
both cases over all 160 and 64 steps (85000 quadrature leaves per step against 2.4e7 with
the fallbacks). FINAL RERUNS on fork 9dbca10 (cluster build 55618725): the intermediate runs
moved aside as *.on748e401-* (the plane run's own-normal N >= 80 cases cancelled unfinished),
plane test 55619126, ring-quadric study 55619127, sphere ladder 55619128; predictions as in
the configs, and the scored rows of the sphere ladder (exact 2.84, ring 1.35, ring quadric
2.11 per period) are predicted to reproduce their digits (their 0 to 3 fallbacks per step
were flat arc pieces, now resolved). Lesson from the cancellations: `ssh host bash -s "$X"`
splits a multi-id string into separate positional parameters (ssh joins the arguments with
spaces), so only $1 arrives; the id lists go into the command string (`cancel_ids.sh`).
BUT the third diagnostic's gate: with zero fallbacks the ring row's oblique arm at N = 160
still ends at 3.0579825799923463e-10 at Co 0.5, the SAME digits as with the four fallbacks of
the 748e401 binary, and at 2.6e-9 at Co 0.2 (7.6e-10 before), the transported ring normal at
8e-10 and 8e-9; the exact-normal row on the same binary is at 2.5e-13 and the split-rule run
of the same case below 8e-14. The residue is therefore not a fallback. Hypothesis (fourth
diagnostic planeTranslation3D_hexParaboloidDiagKB4, prediction in the header): the rule's
per-prism values differ from the split rule's at round-off (1e-14) and the closed loop with
the ring normal amplifies that difference at N = 160 (the sphere's measured round-off floor:
1e-9 to 2e-7 from a 1e-13 perturbation), which the exact-normal row cannot show because its
frame is prescribed; the alternative is a 1e-11 floor of the rule itself. The diagnostic runs
the failing case with volumeRuleDebug (|kb - split| per prism and cell, 8 steps) and its
split-rule twin over the full run. RESULT of the debug case (55619345, 8 steps, 131 s per
step with the comparison): every PRISM volume of the Kromer-Bothe rule is within 4.9e-14 to
6.3e-14 of the split quadrature's (mean 3.2e-15 over 28000 to 29000 prisms per step), but the
CELL volumes, which position the surface through the volume anchor, disagree in a few cells
per step: the maximum over the 25000 cells per step is 2.1e-10 of the cell in steps 3 and 5,
1.2e-11 in steps 1 and 7, and 3e-14 in the other four steps (mean 5e-15 to 2.4e-14). The
residue of the face fractions in this case (3e-10) is therefore the anchor's: a few cell
volumes off by 2e-10 shift the surface by 2e-10 h, and the ring normal carries it along.
The per-prism prediction of the fourth diagnostic held; the cell volumes name the defect's
place but not its cause (which face of the cell polyhedron; the prisms, built from the same
face routine on the same surface, are exact). Fork commit after 9dbca10 (debug branch only:
the first twelve cells per step above 1e-12 are printed with the per-face breakdown); the
diagnostic reruns. NOTE on the binaries: this rebuild happens while the final reruns
(55619126/27/28) run; the change is inside `if (debug_)`, so the cases that start on the new
library compute the same digits as the ones on 9dbca10, and the per-case tpfCommit stamp may
show either.

Appendix to 10.18 (2026-10-25 10:35, the FINAL ring-quadric study, 55619127, 9 of 9 on fork
9dbca10, one case stamped 8329439 as noted; tables on origin). Every case reproduces the
748e401-era digits exactly except N = 128 at Co 0.1 (5.06255e-5 against 5.06254e-5, a relative
2e-6: that case had 2 fallbacks per step before, now 0, so the loop's round-off floor moved
it); the orders are unchanged: q 2.24 / 2.11 / 1.93 per period, normal 1.67 / 1.52 / 1.31,
centroid 2.00, zero fallbacks, zero anchor failures, min alpha -2e-14. The ring-quadric row's
numbers in the preprint, the proposal and the messages to the proposal session stand.

Appendix to 10.18 (2026-10-25 10:50, the split-rule twin, 55619346): the ring row, oblique
arm, N = 160, Co 0.5 with the SPLIT quadrature on the same binary (9dbca10) ends at
L1(alphaf) 2.75e-13, normal 4.07e-13, over the 64 steps: the gate held. So the loop with the
ring normal does not amplify round-off to 3e-10 here; the Kromer-Bothe run's residue is the
rule's own: the few CELL volumes per step off by up to 2e-10 of the cell (fourth diagnostic),
which position the surface through the anchor. The fourth diagnostic's prediction is
therefore FALSIFIED in its conclusion (the per-prism floor and the twin came out as
predicted, but the cell volumes did not): a defect of the rule in the cell polyhedron, not
in the prisms built from the same face routine. The per-cell print (fork 8329439, study
rerunning as 55621064) names the face.

Appendix to 10.18 (2026-10-25 11:00, the FINAL sphere ladder and plane test on fork 9dbca10,
55619128 and 55619126, PARTIAL: the exact and ring rows complete in both, the own-normal rows
still running; tables on origin 5fdc7495, 685ce0aa). SPHERE: the exact row (9 of 9) and the
ring row (9 of 9) reproduce the scored digits (ring N = 64, Co 0.2: 2.118872e-3 against
2.118873e-3, the loop's floor where the 1 to 3 fallbacks per step are now resolved; exact
N = 128, Co 0.1 now 1.098e-5 as before) with ZERO rule fallbacks in every case; the orders
stand (exact 2.84, ring 1.35 per period). The own-normal rows now run at normal speed
(fallbacks in tens to hundreds of prisms per step, not thousands: the flat-arc defect was
most of their cost): the own normal completes at N = 32, Co 0.1 and 0.2 at L1(alphaf) 6e-2
with normal error 1.15 and diverges at step 266 (N = 64, Co 0.2), 208 (N = 128, Co 0.2), 63
and 65 (N = 32 and 64, Co 0.5); the row with the ring normal below 0.1 h^2 completes at
N = 32, Co 0.1 and 0.2 at 4.4e-2 and 5.2e-2 (normal 0.39 and 0.48); the rest pending. PLANE:
the exact row holds the gate in 16 of 16 (7e-15 to 2.5e-13), the ring row in 14 of 16, the
two misses being the oblique arm at N = 160 (2.6e-9 at Co 0.2, 3.1e-10 at Co 0.5: the
Kromer-Bothe cell-volume residue above, the split-rule twin of the Co 0.5 case at 2.75e-13);
the own-normal cases pending.

Appendix to 10.18 (2026-10-25 13:00, the PROFILE of the paraboloid loop on Lichtenberg, by a
profiling agent at Tomislav's request; OpenFOAM profiling trees of the timing study, perf with
DWARF call graphs, and chrono counters in a separate instrumented copy of libVoF under
/work/scratch/tm83tomy/davof/profile/, production library and studies untouched, nothing
committed; jobs 55623213 to 55623490). Ring row, N = 64 (1760 interface cells, 12850 faces set
per step, 4290 with flux), serial, one node: step 0.53 s, of which advect 0.48 (plicRDF 0.079):
the VOLUME ANCHOR 0.273 s = 57 % of advect (6.2 cell-volume evaluations per interface cell per
step: 1 start + 1.46 bracket + 3.73 Illinois; each evaluation 30 us = 24 fan triangles); the
anchor's two FACE-FRACTION DIAGNOSTICS (the "fit residual" before and the "face shift" after
the anchor) 0.067 s = 14 %; the PRISM volumes 0.042 s = 9 % (11 us per prism, 12 fan faces);
the FACE FRACTIONS 0.032 s = 7 % (2.8 us per face; two thirds of them on zero-flux faces);
second-ring search 0.2 %, ring normal 0.5 %, tensor fit 0.8 %, bounding 1.5 %, reconstruct
3 %. N = 128: step 4.34 s, advect 3.84 (plicRDF 0.89), anchor 2.08 (54 %), diagnostics 0.52,
prisms 0.35, fractions 0.26. Diagnostics ON (the ladder's setting) add 0.30 s/step at N = 64,
0.24 of it areaFractionError recomputing the exact sphere state every step. perf self time:
kbIntegrateFace's onArc lambda 17 %, kbIntegrateFace body 16 %, cellSurface::value 7.4 %,
kbFace::setFrame 3.9 %, the allocator family 13 % (Lists and std::function per face), the
adaptive recursion 3.3 %, conic2D::at 2.8 %, sincos 2.3 %. MEASURED remedy: whole planar faces
instead of the fan (`kromerBotheFan false`, the knob exists): step 0.67 -> 0.37 s at N = 64
(-45 %), 3.1 -> 1.9 s at N = 128 (-40 %), anchor 2.7 times cheaper, prisms 1.7 times, numerics
changed at round-off only (L1 0.00195336 -> 0.00195339), 0 fallbacks; the face fractions
ignore the knob (always fan). ESTIMATED remedies: the diagnostics behind a switch (-14 %); a
safeguarded Newton on the anchor with the exact derivative dV/dc = the projected wetted patch
area, one more accumulator in the arcs the rule already walks (3 instead of 6.2 evaluations,
anchor -50 %); the cell polyhedron built once per cell for all anchor evaluations and
stack arrays instead of Lists/std::function (-10 to -15 % of the rule's time); the conic's
rotation reused across the anchor iterations and skipped on uncut faces. Combined estimate:
advect 0.48 -> 0.15 to 0.18 s at N = 64, about twice plicRDF instead of six times. Nothing
of this is implemented; the ladders stand on the fan (kromerBotheFan true).

Appendix to 10.18 (2026-10-25 13:30, the per-cell print of the fourth diagnostic, 55621064,
fork 8329439). Twenty cells printed over the eight steps (diff/V 4.5e-12, 1.2e-11, 2.1e-10;
the same handful of alpha values recur, 0.0369, 0.1600, 0.5076, 0.8400, 0.9631: the same
configuration repeated across cells, as a plane on a Cartesian grid repeats it). Every
printed cell has a VERTEX within 1e-10 of the surface (phi = -1.2e-10 at the vertex against
1e-5 at the others, i.e. 1e-5 of the edge, far above the vertex-on-surface tolerance of
1e-14), so that one or two fan triangles carry a corner sliver of liquid of 2e-20 m^2
against a face area of 1.6e-10 (relative 2e-10). The discrepancy of the cell volume equals
one such sliver's contribution to the face sum (face 3 of cell 3096764: faceTerm 3.0e-25,
the cell's |kb - split| 1.9e-25). The sliver's wetted area by hand from the printed vertex
values (legs 5.4e-10 and 1.5e-10 at 45 degrees: 2.83e-20) is exactly the Kromer-Bothe value
(2.833e-20), so the rule's face value is right there; whether the split quadrature drops the
sliver (a breakpoint within its tolerance of the edge end) or another face of the same cell
is off in the rule is not decided by the log. Open, and bounded: 1e-11 to 2e-10 of the cell
volume in cells with a vertex grazing the surface, the loop's residue 3e-10 at N = 160 on
the ring row's oblique arm, nothing at the sphere ladder's level. The next step, when wanted,
is a unit test in kbSelfTest with a plane grazing a box vertex at 1e-10 against the exact
sliver volume, for both rules. Not pursued further now.

Appendix to 10.18 (2026-10-25 14:30, five agents launched on Tomislav's instruction "do all of
the other points, we have resources"; each reports to the main session, which records here).
Convention: agents that change code work on their own fork branch (`davof-perf`,
`davof-sliver`, `davof-warped`, pushed to tmaric only, never merged by them) in their own
clones under /work/scratch/tm83tomy/davof/agents/<name>/{TwoPhaseFlow,leia} (leia-env.sh
resolves the sibling fork, so their studies run their own build and the production clones
stay untouched); their leia configs carry unique suffixes (Perf, Sliver, Warped, N256) and
land in leia development; job ids go to .my_jobs. (1) P-X2: why the own-faces normal on the
exactly transported fractions of the control reads order 1.40 against the static 2.05, and
why the exact row's own normal error is not zero (cheapest check first: the diagnostic's
step-0 value against the static table; then cell set, reference point, state). (2) The fourth
mesh: the three paraboloid rows at N = 256 (33.6 million cells, parallel, diagnostics at the
end), Co 0.2, for the pairwise order 128 -> 256 of the ring-quadric row (production build, no
code change). (3) Cost: whole planar faces by a per-face planarity test (also for the face
fractions), the anchor's diagnostics behind a switch, a safeguarded Newton anchor on the
exact derivative (the projected wetted patch area), the cell polyhedron built once per cell
and stack arrays; verification by the plane gate, the sphere digits at N = 32/64 and the
timing study on its own build. (4) The vertex-grazing sliver: a kbSelfTest series with a plane
at 1e-6 to 1e-14 from a box vertex against the exact volume for both rules, the responsible
branch fixed, verified on the two failing plane cases. (5) Warped faces: the wetted face area
vector (and first moment) transported instead of the scalar fraction, behind a run-time
switch, verified by the planar reduction and the perturbed-hexahedra plane test with the
paraboloid rows. Merging any branch into davof is a decision for Tomislav after the reports.

Appendix to 10.18 (2026-10-25 15:30, P-X2 RESOLVED by the P-X2 agent; study
sphereTranslation3D_hexExactSurfaceStep0, leia 5fe40fb7, job 55624775; per-cell analysis
under /work/scratch/tm83tomy/davof/px2/). The own-faces normal on the exactly transported
fractions is as good as on the static exact state; the 1.40 is a property of the DIAGNOSTIC
and of the sphere's position on the mesh, not of the transport. (a) Identity: on a sphere
m_c = (A_c/R)(x*_c - c) with x*_c the exact patch centroid, so m_c/|m_c| IS the exact normal
at x*_c; every nonzero L1(n) on an exact state equals the mean tangential offset of the
chosen REFERENCE POINT from x*_c divided by R. The transport diagnostic (reconstructionError)
takes the exact normal at the centroid of the alpha-positioned PLIC polygon, the static test
(leiaTestDavofNormal) at the centroid of the Detrixhe-Aslam zero set; measured offsets 6.6e-3,
4.2e-3, 2.1e-3 h (PLIC) against 6.9e-3, 5.7e-3, 3.1e-3 h (static) at N = 32/64/128, and every
per-cell error equals offset/R to four digits. (b) The end value is Co-independent to nine
digits (the end state is a cell-by-cell translate of the initial one) and the per-step twin
reproduces it at every integer cell shift. (c) The static definition applied to the
TRANSPORTED fractions gives 1.387e-3, 6.707e-4, 1.641e-4 at R/h = 4.8, 9.6, 19.2; the static
study (R/h = 5, 10, 20, 40) gives 2.985e-3, 6.767e-4, 1.638e-4, 4.204e-5: agreement within
1 % at R/h 9.6 and 19.2, pairwise order 2.03 on the fine pair. The "2.1 times larger at
N = 128" of the earlier record compared at equal N instead of equal R/h (the static study
has R = 1 mm in a 4R box); at matched R/h the transport diagnostic is 0.32 to 0.83 of the
static value. (d) The ladder's sphere centre (0.5 0.5 0.5) sits on a mesh vertex at every N;
shifting it by 0.2 to 0.8 h changes the diagnostic's value by factors 0.82 to 1.9 and the
pairwise orders from 1.27 to 1.97: three-rung orders at R/h <= 19 from one symmetric position
carry +-0.5. (e) No O(1) outliers, no wisp effect; the diagnostic's cell set equals the static
threshold's. (f) The exact row's floor (8.2e-5, 1.9e-5, 1.1e-5): exactNormal evaluates the
exact normal at the centroid of the plane positioned along m_c, the scheme then repositions
the plane along that normal, and the error compares the exact normal at the two centroids
(0.4 theta h/R per cell, dominated by grazing cells); the exact normal is exact. Caveat for
per-cell analyses: the scheme's written fields at time t are the reconstruction at the START
of the last step (t - dt), alpha and alphaf the end state. Verdict on P-X2 as worded:
FALSIFIED by its diagnostic, HELD by the identity and under the static definition on the fine
pair; the number leaves the open list. Proposed by the agent and launched: the control ladder
with the sphere centre off the symmetry planes (TPF_SPHERE_CY 0.5137, CZ 0.4871), scored under
both definitions (prediction: static-definition L1(n) on the transported state within 5 % of
the static study at matched R/h, order 2.0 +- 0.2; the PLIC diagnostic 0.7 to 1.0 of it).
Code proposals for the merge round (not applied): exactNormal evaluated at the centroid of
the plane positioned along the exact normal itself (removes the exact row's floor);
reconstructionError to report the error against the exact patch centroid c + R m_c/A_c (zero
by identity: a round-off gate on the fractions) next to the PLIC-reference number and to
write the per-cell error; the table script to state the reference point and to compare
ladders by R/h.

Appendix to 10.18 (2026-10-25 17:00, the off-centre control, 55625094, leia b3d91618; P-X2
CONFIRMED, a new P-X1 finding). Sphere centre (0.5 0.5137 0.4871), N = 32/64/128, Co 0.2, rows
exact and cellLocal. The own-faces normal on the transported fractions is second order at a
generic position under both definitions: static definition 2.70e-3, 6.17e-4, 1.81e-4 (order
1.95, pairwise 2.13/1.77; 0.83, 0.84, 1.02 of the static study scaled to matched R/h), PLIC
diagnostic 1.77e-3, 4.25e-4, 1.14e-4 (1.98, pairwise 2.06/1.90; 0.63 to 0.69 of the static
definition). The agent corrects its earlier "within 1 %": scaled by (R/h)^-2 the control's
single-position values are 0.91 to 0.92 of the static study; the accurate statement is
agreement to 10 to 16 % in magnitude at one position (position scatter of the measure, the
state being exact to 1e-13) and agreement in order (2.03, 1.95 against 2.05). Exact row's floor
0.11 to 0.15 of the cellLocal value, carried by grazing cells (worst 1 % hold 60 %). NEW:
P-X1 fails at this position: L1(alphaf) exceeds 1e-12 on 128 of 320 steps at N = 64 (max
1.44e-8) and 256 of 640 at N = 128 (max 3.3e-9) with ZERO linear cuts, clips, skips and
fallbacks, the end state clean at 1e-13 and mass/boundedness at round-off; the symmetric
position showed none (2e-14 at N = 32 here too). Periodic in the fractional cell shift;
faces not located; presumably the vertex-grazing configuration (STATUS above, the sliver
agent); forwarded to the sliver agent as a test case. N = 256 skipped (16 h, 100 GB).

Appendix to 10.18 (2026-10-25 17:00, the FOURTH MESH, 55624859, leia 7c213fce/90011e91/81374eb;
np 96, 1280 steps, 62 to 66 min per solve, fork 8329439 = 9dbca10 numerically). N = 256 end
values: ringFit L1(alphaf) 2.10e-5, normal 7.22e-5; exact 2.06e-6, 6.97e-7; ringArea 3.90e-4,
3.98e-3; centroid 6.27e-6 (7.16e-6 ringArea), shape 5.74e-6 in all rows; bounded (-1.4e-13),
conservative, zero fallbacks and anchor failures. Pairwise orders 32->64 / 64->128 / 128->256:
ringFit face fractions 2.68 / 1.53 / 1.17 and normal 1.89 / 1.15 / 1.18; exact 2.96 / 2.71 /
2.30 (normal not a converging quantity: 1.34 then 4.01, the reference-point effect); ringArea
1.41 / 1.28 / 1.16 and 0.89 / 1.11 / 1.04; centroid 2.04 and shape 2.00 everywhere. VERDICT:
the computed-normal loop does NOT keep second order: on the finest pair the ring-quadric row
is first order in both the face fractions and the normal, with a constant 18 times below the
ring-area row's (ratio 18.4 at N = 128, 18.6 at N = 256); a fit e = c1 h + c2 h^2 to the two
finest meshes puts 88 % of the N = 256 error in the first-order term. P-N1's "second order
kept" FALSIFIED; P-N3 to P-N6 held; the exact row's 2.30 lies below its predicted 2.5 to 3.0.
CAVEAT: this parallel run skipped 7 to 87 faces per step (4e-4 of the curved faces; a face is
skipped when neither its upwind cell nor its first or second ring carries a surface, which
in parallel happens at processor patches: the serial runs skipped 0, the four-processor
validation 0 to 2), so the digits carry an unquantified parallel effect; the control (the
ringFit row at N = 128 on np 96 against the serial 4.73e-5 / 1.63e-4) is launched before the
N = 256 numbers enter the preprint. Infrastructure: profiles/slurm-fat (64 to 96 GB for the
serial mesh, initialisation, decomposition and reconstruction; the pre-processing peaked at
25 GB), run-studies.sbatch reads PROFILE from the environment, the Makefile lifts the local
task cap; 3000 MB per rank requested, 2.4 GB used.

Appendix to 10.18 (2026-10-25 17:00, WARPED FACES, the agent's branch davof-warped on top of
davof 8a1c5e4: a03dd87, ee31d97, 73a3f57, 60ae03c; leia e6e18aec, 2b2c2c5b; configs
planeTranslation3D_hexParaboloidWarped(.yaml, Ref.yaml), planeTranslation3D_perturbedWarped
(.yaml, Apex.yaml), sphereTranslation3D_perturbedWarped.yaml; tokens TPF_FACE_STATE,
TPF_TENSOR_CENTROID, TPF_KB_PRISM_APEX). DESIGN: the wetted area vector a_f = sum over the fan
triangles of A_wet n_t stored per face (surfaceVectorField wettedAreaDavof) next to alphaf;
no separate first moment is needed because every fan triangle contains the face centre, so
the first moment is (x_f - x_c).a_f exactly (three scalars per face, not six); planar faces
reduce to alpha_f S_f; the update takes the pre-image's wetted area per fan triangle under
the upwind surface; faces no surface sets take alpha_up S_f; the vector is exchanged with
alphaf across processor patches with its sign flipped; the reconstruction sums the vectors.
Two further pieces proved necessary: `tensorCentroid anchored` (cutCellPLIC closes the wet
part of a warped face with a chord, so the alpha-positioned centroids are off the plane and
the ring Monge fit returns a spurious tensor; re-anchoring each centroid by the volume rule
before the fit removes it) and `kromerBothePrismApex faceCentre` (the prism's top and bottom
fanned about the face centre, the cell's own fan, not the vertex mean). Two flux defects of
davof on warped faces found and fixed unconditionally (the prism oriented by the sign of
d.S_f with backward-swept triangles signed in the split fallback; the net liquid volume of a
warped face admitted in [-V_back, V_fwd]): before them the oblique arm failed at 6e-3 to 3e-2
in EVERY combination including davof as it stood, with 220 prism volumes per step rejected.
PLANAR REDUCTION: scalar byte-identical to the reference build in all 8 hexahedral cases;
vector at round-off (7e-15 to 6.5e-14). PERTURBED PLANE (vertices displaced up to 0.1 h,
N = 20/40/80, both arms, Co 0.2): exact row at round-off with anchored centroids for both face
states (1.5e-5 with cut centroids: so the 1 to 5 % failure of the earlier record belonged to
the plane sweep, the dual advection's fan fractions were already right); ring row at round-off
with vector + anchored (scalar: 2e-3 in the normal, flat in N) in 5 of 6 arm/N combinations,
the exception the oblique arm at N = 80 (1.1e-8 / 2.4e-8: single-face jumps at steps 17, 36,
62, 71, the first on a face with a vertex 2.2e-4 h from the plane: vertex grazing again);
apex control: vertexMean fails at 2e-5 to 3.5e-5 in every combination. PERTURBED SPHERE
(N = 32/64, three rows, 12 cases, no abort): NOT bounded with either face state (min alpha
-5e-4 to -3.3e-3 at N = 64, anchor failures in 7 to 1163 cells per run; hexahedra: -8e-13);
centroid 2 to 2.3 times the hexahedral value, shape error below it; the vector state matters
for the ring-quadric normal (9.1e-4 against 2.6e-3 at N = 64; two-rung slopes 1.25 against
0.45). OPEN: the vertex-grazing events (debug copy at agents/warped/debug/c47), the perturbed
sphere's boundedness and anchor failures, the parallel exchange untested beyond one
processor, the exactSphere rule's prism path on backward-swept triangles. Merging is for
Tomislav.

Appendix to 10.18 (2026-10-25 18:30, COST, the agent's branch davof-perf: ef142fe whole faces
by a per-face planarity test (`kbFace::warp` > 1e-10 L_s fans; `kromerBotheFan auto|true|false`,
also in the face fractions and prism tops), 8a1c5e4 the anchor's diagnostics behind
`anchorDiagnostics` (default off), decb99d the safeguarded Newton anchor on dV/dc = oint t1 dt2
(`kbSums::patchSum`; rtsafe of Numerical Recipes 9.4 with a Halley correction, Scavo and Thoo
1995, doi:10.2307/2975033; Illinois fallback; `anchorMethod illinois|newton`), a0460dc the
cell faces cached per (cell, time index), stack lists and plain recursion, f6090c4 the
uncut-face skip by the quadratic variation bound, 1e824cf the default `anchorMethod illinois`
and `anchorMethod compare`; leia cd450914 (converter: nAnchorEvals) and a4b1b309 (token
TPF_ANCHOR_METHOD, default illinois); README updated with the citations). VERIFIED: steps 2,
4 and 5 bit-identical; `kromerBotheFan true` + illinois reproduces the ladder's ring N = 32
digits and production's two plane misses bit for bit. DEFAULT BUILD (whole faces, illinois):
sphere ring 5.647047e-3 against 5.647061e-3 (2.5e-6), 2.118887e-3 against 2.118872e-3 (7e-6),
exact rows equal, everything else at round-off, anchor failures per step 2 against 3; plane
gate 30 of 32 below 2.5e-13; the two vertex-grazing misses move with the round-off of the
path (Co 0.2: 1.67e-9 from 2.63e-9; Co 0.5: 2.59e-9 from 3.06e-10, i.e. WORSE, and 2.8e-13 with
newton + fans: the "not worse" clause fails for that case until the sliver fix). TIMING
(medians of cpu per step, N = 32/64/128, Co 0.2): default build ring 0.046 / 0.333 / 1.35 s
against plicRDF 0.030 / 0.164 / 1.35 in the same study, ratios 1.5 / 2.0 / 1.0 (reference 6.0 /
3.5 / 3.0): raw speed-up 2.8 / 1.6 / 3.3, normalised 3.9 / 1.7 / 3.0; with newton 0.041 / 0.214 /
1.71 (ratios 1.4 / 1.3 / 1.3; the N = 128 figure confounded by node scatter: plicRDF 1.05 to
1.76 s). Anchor evaluations per cell and step: illinois 6.1 to 6.6, newton 3.1 to 4.0. NEWTON
DOES NOT VERIFY on the ring row (0.5 % lower at N = 32, 2.5 % at N = 64, with 0 anchor failures
against 70 and 148): the cause, found with `anchorMethod compare`, is a METHOD DEFINITION GAP:
cells whose alpha is exactly 0 or 1 but which are interface cells by |m_c| > wispTol have a
flat V(c) over an interval (every offset there is a root); Illinois stops at the first bracket
end of its doubling (0.75 h or 1.75 h from the start), Newton at the interval's edge (0.72 to
0.80 h), and because that surface also cuts the face prisms and second-ring faces outside
the cell the two loops diverge: THE LADDER'S RING DIGITS DEPEND ON AN ARBITRARY BRACKET END
in those cells. Decision for Tomislav: define the offset of such a cell (the interval edge,
the surface touching the cell, is the natural choice), make newton the default, rerun the
ladders. Not done: step 5(a) (conic rotation reuse, would change digits), parallel and
polyhedral tests of the planarity test, a new profile. Note: the shared laptop clone's libVoF
was the agent's step-2 build between 13:04 and 14:27 (no laptop runs happened); it is rebuilt
from davof 8329439 and the clone is on davof with a clean tree.

Appendix to 10.18 (2026-10-25 19:00, the PARALLEL CONTROL of the skipped faces, 55628348, leia
41fd73fb and e317c31): the ring-quadric and exact rows at N = 128, Co 0.2, on 96 ranks against
the serial ladder. End values differ by 4.5e-5 (ringFit L1(alphaf)), 4.1e-5 (its normal), 6.0e-4
(exact L1(alphaf)), 3.6e-4 (its normal), 1e-5 (centroid, shape) relative; bounded and
conservative as serial; 0.65 s per step on 96 ranks against 4.9 serial (7.5 min per solve).
Skipped faces per step at N = 128 on 96 ranks: 0 / 13 / 40 (min / median / max), nonzero in 638 of
640 steps; the serial twin 0 in every step; at N = 256 7 / 42 / 87 (2.1e-4 of the curved faces);
the skipped count is IDENTICAL between the ringFit and exact rows at every step (a property of
the candidate set and the decomposition, not of the surfaces). P-C1 (harmless at 1e-3) HELD
for both rows (the exact row's Linf at 1.3e-3 the only value above 1e-3); P-C2 (a decomposition
effect) HELD; P-C3 (cost) HELD. Consequence: a perturbation of 6e-4 moves the pairwise orders
128 -> 256 by 1e-3: the N = 256 verdict (the computed-normal loop first order) STANDS and its
digits enter the preprint. The mechanism (the agent's reading of the code, not demonstrated
per face): a candidate face of a non-interface cell whose nearest interface cell lies beyond
the one remote second-ring surface the processor exchange passes (corners of three parts,
two patches away) finds no surface and keeps the flux set before the curved pass. OPEN for the
merge round: extend the processor exchange (two layers of remote surfaces, or exchange the
second-ring choice itself) so that parallel runs skip none; the serial N = 256 twin (14 h on
one core) was not run.

Appendix to 10.18 (2026-10-25 19:30, the SLIVER defects SETTLED, the agent's branch davof-sliver
on davof 8329439: 5991916 self-test series and anchor event print, 09ab419 split-rule fix,
58cc582 Kromer-Bothe fix, 45436d4/2741c52 self-test additions; leia c7e1a845, 3280e7c9,
b2410d4a with the configs *Sliver*). TWO defects, neither the one I guessed. (1) The printed
"cell mismatch" of 2e-10 to 4.5e-12 was the SPLIT rule's error, not a sliver lost by
Kromer-Bothe: the exact plane volume of cell 3096764 from the printed vertices (mpmath, two
formulas to 2e-28) puts Kromer-Bothe 3.8e-15 and the split rule 1.23e-11 off; cause: a grazing
vertex as the sweep apex of a pyramid triangle puts the sweep-line root within 7.5e-6 of the
outer parameter, a layer the breakpoints do not mark and depth 9 does not resolve
(triangleIntegralSplit via pyramidLiquidVolume/prismLiquidVolume); fix: the sweep apex
rotated to the vertex farthest from the traces when the given apex is within 1e-3 of the
largest distance (split cells 1.2e-11 -> 8.2e-15 in the self-test). (2) The plane gate's two
misses were single-step anchor events (L1 at 2e-13, a jump at step 12 (Co 0.5) and 102 (Co 0.2)
with Linf 4e-6 and 4e-5): in one cell the Kromer-Bothe cell volume toggled by 1.2e-5 V between
offsets 1e-16 apart, so the anchor landed on a wrong offset; cause: the plane's ring tensor is
round-off (b ~ 1e-8 1/m), B/l in the face's scaled coordinates between 1e-13 and 1e-9 sits just
ABOVE the LINE threshold of 1e-13, the trace is classified as an ellipse or hyperbola with its
centre 1e12 Ls away, its parameter carries noise of 1e-4 Ls, and a corner sliver arc with a
chord of 3e-4 Ls is dropped (by sb > sa or by the inside test at the parameter midpoint); the
fan triangle's wet area ends 1.5e-4 off, inside [0, A], so no fallback; the same mechanism with
b = 0 on faces nearly parallel to the plane. Fix: a trace with max|B| <= 1e-3 |l| is ORDERED as a
line (arcs still integrated with the exact quadric), an arc between two exact crossings is
kept when their parameters coincide (its chord length decides), and the inside test of a
one-piece arc is made at the arc's point over the chord midpoint. Debug knobs kbGrazingFixOn,
splitApexRotationOn (default on). SELF-TEST (2 boxes x 4 normals x 3 grazed points x 2 sides,
d = 1e-6 to 1e-14): Kromer-Bothe cells exact (8e-16) at b = 0 before and after; under weak
curvature (b L = -1e-10, 1e-12) before 4.5e-7 to 4.5e-11 of V scaling with d (the dropped
sliver's signature), after <= 8e-16; event-cell replicas scanning 201 offsets within 1e-15:
1.6e-5 and 1.1e-5 V before, 1.4e-14 and 6.2e-15 after. VERIFIED on the agent's build: the two
failing plane cases PASS (2.30e-13 and 2.83e-13, from 2.63e-9 and 3.06e-10; 0 fallbacks, 0 anchor
events); no regression on 16 cases (exact and ring rows, N = 40/80, both arms, Co 0.2/0.5; all
<= 1.7e-13; exact rows digit-identical, ring rows at 1e-14); the sphere ring row at N = 32
bit-identical over 160 steps. The OFF-CENTRE control's excursions (1.44e-8, 128 steps) are
bit-identical before and after: that path (tetOverlap + the quadratic cut faceCutDisplaced)
calls none of the changed code; a grazing problem of the quadratic cut is the likely analogue,
not shown. OPEN: fallback (unresolved) cases when d ~ |b L| with three vertices grazed at once
and for the axis-near normal (1, 1e-7, 3e-7) (counted fallbacks, 0 in the plane runs); a
Kromer-Bothe face error of 4.6e-11 on faces tilted 3e-7 against the plane (conditioning of the
trace position); one axis-near case where Kromer-Bothe and the quadratic cut disagree by
1.2e-3 on a face with no valid reference; the fix untested in parallel; the apex rotation may
move split-rule values within the quadrature tolerance on curved surfaces. Merging is for
Tomislav. Correction to the records above: "which rule drops the sliver is open" is answered
(the split rule, in its sweep-apex layer), and the plane gate's residue was a separate
Kromer-Bothe arc loss, now fixed on the branch.

Appendix to 10.18 (2026-10-25 20:30, Tomislav's DECISIONS and the merge). Decisions: (1) the
offset of a cell whose alpha is exactly 0 or 1 while it counts as an interface cell is DEFINED
as the edge of the flat interval, the surface tangent to the cell; (2) the anchor becomes the
Consecutive Cubic Spline positioning (Marić, J. Comput. Phys.: X 11 (2021) 100093,
doi:10.1016/j.jcpx.2021.100093) applied to the paraboloid: between vertex crossings the iterate
is the root of the cubic Hermite interpolant through the bracket ends with the exact V and the
exact dV/dc = A_p, re-bracketed, bisection-safeguarded; `anchorMethod ccs`, the default;
(3) merge order sliver, perf, warped, then one cluster build and the reruns of the three
ladders and the plane test on the merged binary. DONE: davof-sliver merged into davof as
f133781 (no conflicts, 4 files; built on the laptop, pushed to tmaric). The cost agent is
rebasing davof-perf onto f133781 and implementing the CCS anchor with the acceptance criteria
revised for the new definition (exact row at round-off; the ring row's change traced to the
flat-interval cells with `anchorMethod compare`, expected 0.5 to 2.5 %; plane gate 16 of 16;
evaluations per cell 2 to 3; timing). The production library on the cluster stays at
9dbca10/8329439 until the final own-normal cases of the last reruns are done and all three
branches are merged.

Appendix to 10.18 (2026-10-25 21:30, davof-warped REBASED onto f133781: tmaric/davof-warped =
94dddd3, five commits; the pre-rebase head kept as the agent's local davof-warped-pre-rebase;
its base had carried two perf commits that f133781 lacks, dropped and returning with the perf
merge; conflicts in kbFanFaces/kbOrientedFan/kbFaceWettedArea resolved to f133781's
kromerBotheFan bool, to be switched to kbWholeFace after perf). Verified on the agent's clones
against a build of f133781 (jobs 55631503, 55631514, 55631702/03/07): scalar byte-identical in
all 8 hexahedral cases, vector at round-off (6.3e-14); the perturbed plane (48 cases rerun): the
twelfth case (ring row, oblique arm, N = 80) now at 7.0e-14 / 2.2e-16 (from 1.1e-8 / 2.4e-8: its
jumps were the dropped sliver arc), so vector + anchored centroids + face-centre apex is at
round-off in both rows, both arms, every N (P-W2 holds in full); the other combinations
unchanged. leia c132e854 (results in the config headers).

Appendix to 10.18 (2026-10-25 23:30, davof-perf MERGED as davof 1f914f9 after the agent's second
report: head 3f5270b = e5b6205 merge of f133781, 46db35e the edge definition and the CCS anchor,
3f5270b compare classes; leia 30d68880 TPF_ANCHOR_METHOD default ccs, converter column
nAnchorEdge, *Perf3 configs). The CCS anchor (Marić 2021, doi:10.1016/j.jcpx.2021.100093,
arXiv:2012.02227): Newton step from the start until bracketed, then the root of the cubic
Hermite interpolant through the bracket ends with exact V and dV/dc, re-bracketed, bisection
when the iterate leaves the bracket or the bracket is not halved in two steps; the edge rule
(kbOffsetRange: the exact range of phi + c over the cell from vertices and edge- and
face-plane stationary points; a target within 1e-14 V of 0 or V takes cMin or cMax moved
1e-12 h onto the flat side, no volume evaluation, counted as `edge N`). VERIFIED on the
agent's build: exact row digits equal to the ladder (L1(n) to 6e-7); ring row -0.53 % (N = 32)
and -2.47 % (N = 64) against the ladder, EVERY differing cell traced by `anchorMethod compare`
to the edge definition (1294 of 1302, 1322 of 1519) or to differences below the anchor's
resolution (|dc| A <= 1e-13 V; nearly full cells with corner slivers) or to Illinois failures,
none other; anchor failures 0 (ladder 70 and 148), fallbacks 0, mass 8e-17, bound 9e-14. Plane
gate 32 of 32 at round-off (max 2.75e-13; f133781 gives 2.83e-13 on that case). Evaluations
per cell and step 3.4 (ccs) against 6.3 (illinois) and 3.5 (newton): start, one Newton step,
the Hermite root and a confirming evaluation for the offset criterion. TIMING (medians of
three, Co 0.2): ring 0.051 / 0.248 / 1.457 s against plicRDF 0.030 / 0.204 / 1.345 in the same
study, ratios 1.70 / 1.21 / 1.08 (reference 6.0 / 3.5 / 3.0); raw 2.5 / 2.2 / 3.1 times faster
than the ladder's binary. Self-test bit-identical to f133781 with fans over all 534 lines, at
round-off with whole faces. Notes: the sliver fix's anchor event print now needs
anchorDiagnostics true; builds before the merge reject `kromerBotheFan auto`; leia's
sphereTranslation3D.parameter had TPF_KB_FAN ( true ), set to auto now so that the reruns use
whole faces. The warped agent rebases onto 1f914f9 next; then the warped merge, the cluster
build and the four reruns.

Appendix to 10.18 (2026-10-26 00:30, davof-warped MERGED: tmaric/davof-warped bbd12e4 rebased
onto 1f914f9 (scalar byte-identical to a 1f914f9 build in all 8 hexahedral cases, vector at
round-off, the perturbed plane at 5e-15 to 6.3e-14 in every exact and ring run with vector +
anchored + faceCentre; leia 36a5d971), merged as davof 7c06c57 (16 files; the VoF library
built on the laptop). ALL THREE BRANCHES ARE IN: davof 7c06c57 is the final binary of the
preprint. Next: the cluster build of 7c06c57 and the reruns of the plane test, the sphere
ladder, the ring-quadric study, the timing study and the N = 256 study (with the fat profile)
on it; the superseded own-normal cases of the pre-merge reruns are cancelled by id and their
studies kept aside as *.preMerge-*. SUBMITTED (2026-10-26 00:45): build 55634206; plane test
55634207, ring-quadric study 55634208, sphere ladder 55634209, timing study 55634210, N = 256
55634211 (afterok the build; the N = 256 with PROFILE=profiles/slurm-fat); the five pre-merge
study directories moved aside as *.preMerge-10251913; the pre-merge drivers 55619126/55619128
and their nine remaining own-normal solves cancelled by id. The settings of the reruns: the
merged binary's defaults (kromerBotheFan auto, anchorMethod ccs, faceState scalar,
anchorDiagnostics off), so the ring row's digits move by up to 2.5 % against the pre-merge
tables (the edge definition) and the plane table should read 16 of 16 for exact and ring.

Appendix to 10.18 (2026-10-26 01:20, the FINAL PLANE TEST on davof 7c06c57, 55634207, 48 of 48
in 24 min; tables on origin, every case stamped 7c06c57). Exact row 16 of 16 at round-off (max
L1(alphaf) 2.4e-13), ring row 16 of 16 (max 2.75e-13; the N = 160, Co 0.5 normal-arm entry of the
preprint's table 9.8e-14), own-normal row 1 of 16 (at round-off only at N = 20, Co 0.2; 1.7e-8
at N = 20, Co 0.5; 7e-2 to 1.5e-1 at N >= 40 on the normal arm at Co 0.5 and at N = 160 on both
arms; 3e-11 to 8e-7 in between; up to 10211 rule fallbacks per step on its wild surfaces, 12.9 s
per step at most, which is why these cases finish now where they crawled for hours on the
fan binary with the flat-arc fallbacks). The plane gate of the dual advection with the
Kromer-Bothe rule is therefore held by the exact and the ring normal at every N, Co and arm on
the final binary, and the own normal fails it beyond N = 20 as in the plane-clipped loop.

Appendix to 10.18 (2026-10-26 01:50, the FINAL TIMING on davof 7c06c57, 55634210, 18 of 18, every
case stamped 7c06c57; tables on origin). Medians of three per step (Co 0.2, diagnostics at the
end): ring row 0.053 / 0.358 / 2.117 s against plicRDF 0.031 / 0.184 / 1.829 s in the same study
at N = 32 / 64 / 128, ratios 1.7 / 1.9 / 1.16 (the N = 128 pair on a slower node than the agent's
run: 1.46 against 1.35 there; the ratio is what carries). The ring row's end digits on the final
binary: 5.616950e-3 (N = 32), 2.066471e-3 (N = 64), 8.709252e-4 (N = 128), identical in the three
repetitions, i.e. -0.53 %, -2.47 % and -0.003 % against the fan/Illinois ladder (the edge
definition's cells are coarse-mesh cases); zero fallbacks, zero anchor failures.

Appendix to 10.18 (2026-10-26 02:20, the FINAL N = 256 on davof 7c06c57, 55634211, 3 of 3 in 55
min with the fat profile, np 96; tables on origin). End values: ringFit L1(alphaf) 2.1015e-5,
normal 7.2214e-5; exact 2.0648e-6, 6.9690e-7; ringArea 3.8992e-4, 3.9805e-3; centroid 6.27e-6
(7.16e-6 ringArea), shape 5.738e-6; min alpha -1.9e-12; zero fallbacks and anchor failures. These
reproduce the pre-merge N = 256 run to three or four digits (2.10e-5 / 7.22e-5, 2.06e-6 /
6.97e-7, 3.90e-4 / 3.98e-3), so the fourth-mesh orders (ringFit 1.17 / 1.18, exact 2.30, ringArea
1.16 / 1.04 on 128 -> 256) stand on the final binary; the preprint's printed digits need no
change. The step now costs about half (55 min against 62 to 66 for the solve including
pre-processing on the merged binary).

Appendix to 10.18 (2026-10-26 03:00, the FINAL SPHERE LADDER on davof 7c06c57, 55634209, 36 of 36
in 1 h 41 min, every case stamped 7c06c57; tables on origin). Exact row: 5.2050e-4, 6.6668e-5,
1.0183e-5 at Co 0.2 (orders 2.84 / 2.87 / 2.87 at Co 0.2 / 0.1 / 0.5), digits equal to the
pre-merge ladder; own-faces normal 2.05; zero fallbacks and anchor failures. Ring row:
5.6169e-3, 2.0665e-3, 8.7093e-4 at Co 0.2 (order 1.34; Co 0.1: 5.7508e-3, 2.0047e-3, 8.5191e-4,
1.38; Co 0.5: 7.2302e-3, 2.5182e-3, 1.2179e-3, 1.28), normal 1.00 at 8.19e-3, bounded to -8e-13,
zero fallbacks and anchor failures (the edge definition replaces the 70 to 148 failures of the
Illinois ladder). OWN-NORMAL ROWS (P-P3 settled): the own normal DIVERGES at every N and Co on
the final binary, aborted at steps 66 / 76 / 82 (Co 0.1), 46 / 63 / 68 (Co 0.2), 45 / 51 / 54
(Co 0.5) for N = 32 / 64 / 128 with alpha reaching the abort bound, where on the pre-merge
binaries two N = 32 cases had limped to the end at 6e-2 without tripping the guard; the row
with the ring normal below 0.1 h^2 diverges everywhere except N = 32, Co 0.2 (completed at
3.85e-2, normal 0.35), aborted at steps 41 to 317. So P-P3's first clause (own normal diverges
at every N) HELD, its third (the sliver rule bounded at N = 64 and 128) FALSIFIED: the ring
normal below 0.1 h^2 is not enough, the whole loop needs the filtered normal. Hundreds to
thousands of rule fallbacks per step in both own-normal rows (wild surfaces), none in the
others.
