# plan-curvature-stabilization.md

Curvature-error / parasitic-current stabilization for the semi-Lagrangian level
set (leia, `leiaSemiLagrangianLevelSetTwoPhaseFoam`). Working document for a
coding agent — **v0.2** (v0.1 + the 2026-08-07 review amendments; deltas marked
**[v0.2]**).

Scope: the curvature estimation → delivery → CSF force path and the transported-ψ
profile feedback. Out of scope here (separate plans): band-restricted advection,
reverse-mode LENT reset design beyond the WP6 gate, DED application runs.

---

## 0. Ground truth — measured facts the agent must not re-litigate

From the quadratic-SL deck (`quadratic-semi-lagrangian-level-set.template.html`),
the comparison deck (`level-set-method-comparison.template.html`), and the
curated CSVs under `docs/*/data/tables/`. All new work is scored against these
numbers; none of them may regress.

**Transport (do not touch):**

| case | mesh | shape order | volume order |
|---|---|---|---|
| 2D reversed vortex, T=2, CFL ½ | hex 32→256 | 2.97 | 2.98 |
| 2D reversed vortex, CFL 1 | hex | 2.59 | 2.94 |
| 3D shear, T=3 | hex 32→128 | 2.50 | 3.25 |
| 3D deformation (LeVeque) | hex 32→160 | 1.36 | 1.86 |
| 3D deformation | poly (cfMesh) 1/32→1/128 | 1.46 | 2.51 |

Serial ↔ np=4 agreement ~1e-12 (halo exchange test). Cached and uncached
QWLSR agree cell-for-cell to ~1e-13.

**Curvature statics (exact SDF, |Sf|-weighted L2 on active faces):**
- 2D circle: every classical delivery ≈ O(h^1.1). Cell-centred κ: O(h^1.07),
  rel. L2 35.4% → 1.7%.
- 2D stabilized foot-point face re-referencing: **O(h^2.04)**,
  11.35 → 0.105 m⁻¹ at N=512 (108×).
- Per-side-then-combine inverse ordering: O(h^1.40) — order lost.
- **[v0.2] 3D sphere (faceCurvatureSphere3D, N=32..128):** raw delivery
  O(h^0.98); Gaussian-curvature-aware inverse **O(h^1.95)**
  (0.419 → 5.25e-3 m⁻¹ at N=128, 80×); per-side variant O(h^1.98);
  the K-less 2D scalar inverse as a CONTROL measures **O(h^1.02) at raw
  magnitude** (0.421 vs 0.419) — the Gaussian term is load-bearing.

**Coupled (static droplet, water/air, R=1 mm):**
- t_blow [s] vs N (arithmetic delivery): 0.47 (32), 0.44 (64), 0.11 (128),
  0.03 (256). **[v0.2]** With the h² foot-point delivery: onset (max|U|>0.1)
  0.033 (256), **0.010 (512)** — refit the exponent against the full onset set
  before using it as the WP7 yardstick (closer to N^(−1.7) at the fine end
  than the earlier N^(−3/2)).
- Instability is h-driven, transport-coupled (freezing ψ removes runaway).
- Ruled out by instrumentation: force imbalance, κ formula, capillary CFL,
  |∇ψ|-degeneracy (|∇ψ| ≈ 0.85–1.1 at spikes), density ratio.
- Mechanism: parasitic strain pumps a 2h normal-profile mode in ψ; on diagonal
  stencils the h/√2 sampling is incommensurate with the mode; LSQ misfit leaks
  into the tangential Hessian blocks (H_tt, H_nt) which the κ formula does NOT
  cancel; measured κ ∈ [−883, 2744] vs exact 953 at ε = 0.5. Grid-aligned poles
  cancel (error lands in H_nn, formula-invariant): ±10% at the same ε.
- Foot-point re-referencing: seed 30× lower, onset (exponent) unchanged.
- **[v0.2] Biharmonic band ψ-filter (θ = 0.05): a DELAY device, not a closure.**
  It carries N = 64–256 through the t = 0.1 gate horizon with the band |∇ψ|
  pinned and the κ-error flat, but the extended-horizon run (N=256, post
  seam-fix, np 8) blows up at **t = 0.167** (onset >0.1 m/s at 0.113,
  e-folding ~27 ms) — a ~5× fuse extension over the unfiltered 0.035. This is
  the paired-combination confirmation of the filter-alone lesson (feedback
  migrates to longer undamped wavelengths). WP3 characterizes it; it is no
  longer a closure candidate.
- Offline refit of measured 3×3 stencils reproduces every solver κ spike to
  0.00% (100% sign match): spikes are the fit's true response, not a bug.
- **[v0.2] Processor-seam consistency:** face values computed from rank-local
  fits MUST be synchronized across processor seams (swap-and-average,
  syncTools). The unsynchronized one-sided correction made np-8 runs 7–10×
  noisier than np-4 at identical (N, t); post-fix np4 ↔ np8 agree to 1.1% at
  t = 5e-3 (N=128 filtered). Every np-8 verdict taken pre-fix was discarded
  and re-measured.

**Preconditions currently enforced by the solver:**
- `stabilizedFootPointFace` requires `offsetCorrection none` (else double-correct;
  solver refuses) and a quadratic reconstruction (only the uncached QWLSR
  implements `footPointDistance`/`fitDerivatives`).
- Selected by `levelSet.curvatureExtension { type stabilizedFootPointFace; }` +
  `surfaceTensionForce { faceCurvatureSource registered; }` (registered face field
  `kappaStableFootFace`).

---

## 1. Problem decomposition

- **A — static seed:** SOLVED (foot-point re-referencing, O(h^2.04) 2D /
  O(h^1.95) 3D **[v0.2]**). No further work.
- **B — feedback exponent:** OPEN. The 2h mode → tangential-Hessian aliasing →
  κ speckle → force → strain loop — **[v0.2]** with the caveat that the
  filtered blow-up implies the loop can also run on longer-wavelength content;
  WP0 measures which. All leverage is UPSTREAM of the fitted H. WPs 0–3 and
  the promoted WP6 attack this.
- **C — pinch-off singularity:** |∇ψ| → 0 at genuine topology change. Localized,
  transient. Needs guards, not fixes. WP4's guard generalization covers it.

Rule for the agent: any proposed change must state which of A/B/C it targets.
Changes targeting A are rejected by default (measured: seed reduction does not
move the onset).

---

## 2. Work packages

### WP0 — band mode-SPECTRUM diagnostic (do first; everything downstream scores on it)

`minGradPsiBand` measurably does NOT lead the instability (|∇ψ| ≈ 0.9 at the
spikes). Add the quantity that does. **[v0.2] Measure a small spectrum
(2h, 4h, 8h), not the 2h mode alone**: the filter damps 2h near-maximally yet
the extended-horizon run still blew up — a 2h-only indicator could sit quiet
while the run dies, destroying its value as the WP6 trigger and the WP7 score.

Tasks:
1. Per band cell c: n̂ from the existing LLS plane (halo-only, fit-free;
   fallback: max |n̂·(x_N − x_c)| alignment); build a ±n̂ face-neighbour chain
   3 cells deep per side. Amplitudes from divided differences along the chain:
   `A2h_c = |ψ_up − 2 ψ_c + ψ_dn| / 4` (pure 2h mode of amplitude a gives
   A2h_c = a); A4h/A8h from the wider stations, one calibration constant per
   mode fixed against the WP1 harness injection. Optionally report a detrended
   variant (subtract the smooth h²κ|∇ψ|/4 background); decide from calibration.
2. Band norms `A{2h,4h,8h}_max`, `A{2h,4h,8h}_2` (volume-weighted), written to
   the curated CSVs next to `minGradPsiBand(Half)` every write time. Halo-only,
   fit-free.
3. Retrodiction: (a) N=128 arithmetic onset case with the diagnostic on — does
   the spectrum rise ahead of t_blow while minGradPsiBand stays flat?
   (b) **[v0.2]** N=256 filtered extended-horizon case — which mode carries the
   filtered runaway (2h pinned + 4h/8h growing ⇒ migration; 2h growing ⇒ θ
   insufficiency)? The answer retargets WP3/WP6.
4. Acceptance: cheap (<1% step cost), parallel-consistent (serial ↔ np=4 to
   ~1e-12 AND np4 ↔ np8 **[v0.2]**), leads t_blow by a usable margin on (a),
   carrier identified on (b).

### WP1 — Offline estimator-robustness study (anisotropic geometry fit)

**Background and one correction.** A plain rotation of the quadratic basis into
the interface frame is a NO-OP: the full quadratic space is rotation-invariant,
so the fitted polynomial (hence g, H, κ) is identical in any basis. (Proposed in
discussion; checked; vacuous. Kept below only as a harness sanity arm.)
The non-vacuous, mechanism-informed version: the leak is LSQ **misfit** of the
under-resolved normal mode spilling into tangential components. On the 2D
diagonal 3×3 stencil the ζ-stations (ζ = distance along n̂) are exactly
{−√2h, −h/√2, 0, h/√2, √2h} — five values. A polynomial **quartic in ζ**
interpolates any pure function of ζ at five stations exactly → zero misfit →
nothing left to spill into τ², τζ. So the candidate fix is an
**anisotropic-degree fit in the interface frame**: degree 4 along ζ, degree 2
along τ, cross terms {τζ, τζ²} — a genuine model change, not a basis change.

Tasks:
1. Extend the existing offline harness (`workflow/scripts/make_profile_mode_fig.py`
   machinery — the one that reproduced solver spikes to 0.00%) with a synthetic
   injector: ψ = d(x) + ε·h·sin(π d/h + φ), circle R = 1 mm in the 0.01 m box,
   matching the deck's ε convention; ε ∈ {0.05, 0.1, 0.2, 0.5}, φ swept.
2. Stencil placements: grid-aligned pole and 45° diagonal (2D), N-equivalent
   h for N ∈ {64, 128, 256}.
3. Arms:
   a. baseline quadratic WLS (current production) — reproduces deck numbers;
   b. rotated-basis quadratic — **sanity check: must equal (a) to machine
      precision**; if not, the harness is buggy;
   c. anisotropic ζ-quartic × τ-quadratic in the interface frame
      (normal source: exact; then LLS-plane normal; then lagged g/|g|);
   d. anisotropic IDW weights (tighter tangentially) on the plain quadratic —
      cheap cousin, one arm only.
4. Measure per arm: κ error at the cell centre; block decomposition of the
   fitted-H error (H_nn vs H_tt/H_nt); scaling of the diagonal error with ε
   (expect ~ε/h for (a); ~ε²/h for (c) with perturbed normal — the normal's own
   O(ε) error re-leaks a fraction).
5. Conditioning: record the LSQ condition number per arm; on the diagonal the
   ζ-stations are near-Nyquist — Runge/ill-conditioning risk for degree 4.
   Define the fallback ladder now: ζ-quartic → quadratic → linear → constant,
   with a **centered, scale-free** degeneracy test (see §4 conventions — the
   absolute-coordinate determinant guard already burned us once at N=512).
6. Acceptance to proceed to WP2: arm (c) reduces the diagonal κ error at
   ε = 0.5 by ≥ 1 order of magnitude vs (a), with condition numbers bounded
   and the pole case not degraded.
7. 3D note (open question, do not block on it): along a body diagonal the
   3×3×3 stencil has 7 ζ-stations; degree-6 is not realistic. Options to study
   later: degree-4 + regularization, or rely on WP3's filter to carry more of
   the load in 3D. Record findings; decide after 2D results.

### WP2 — Solver integration of the winning geometry fit

Only after WP1 acceptance.

1. New RTS model (e.g. `interfaceAlignedAnisotropicFit`) in the existing
   pluggable hierarchy — used ONLY by the geometry path (κ estimation, foot
   Newton, d_f). **Two-object discipline:** the transport reconstruction
   (`quadraticWeightedLeastSquares` / uncached) is not touched, not shared,
   not limited by this model.
2. Band-only: fit in the α ∈ (0,1) + k face-ring band. Caching: the rotation
   makes the basis per-cell/per-step → cached pseudo-inverse invalid → use the
   uncached (per-cell Cholesky) path. Cost bound: band is O(N^(d−1)) cells;
   report measured step-cost delta (expect small; uncached full-domain was
   1.4× CPU, this is band-only).
3. Consumers unchanged in interface: model returns (value eval, g, H) at
   requested points; foot Newton and seed formula
   d0 = R/|g| + H_nn R²/(2|g|³) work as-is (zero set no longer a conic —
   Newton unaffected).
4. Verification: leiaTestMeanCurvature circle ladder (statics must stay
   ≥ O(h^2.04) with re-referencing on) **[v0.2] and the 3D sphere ladder
   ≈ O(h^1.95)**; transport suite untouched by construction but re-run as
   non-regression anyway; serial ↔ np=4 **and np4 ↔ np8 [v0.2]**.

### WP3 — Biharmonic band ψ-filter: characterize the DELAY device **[v0.2 reframed]**

The filter (θ = 0.05) is measured as a ~5× fuse extension at N=256, not a
closure (blow-up t = 0.167 on the extended horizon). Characterize it as the
reference delay arm:

1. Convergence signature: run the full transport suite (2D vortex, 3D shear,
   3D deformation; hex + poly) with the filter ON. Acceptance: shape/volume
   orders unchanged within fit tolerance (it acts on the band profile only,
   but measure, don't assume). Also record its effect on
   `gradientErrorBandHalf` — expected to improve.
2. Scale study: θ ∈ {0.025, 0.05, 0.1} on the N=256 extended-horizon case,
   **A-spectrum instrumented [v0.2]** — the deliverable is delay-vs-θ scaling
   and the carrier mode, input to WP7's interpretation.
3. 3D smoke test on the shear case.
4. Interaction: filter and WP2 fit compose (filter reduces the mode, fit is
   robust to the remainder); the arm comparison happens in WP7, not here.

### WP4 — 3D parallel-surface inverse (Gaussian term) + guards **[v0.2: DONE, one residual]**

Landed (commit b7e947e) exactly as specified: κ^Γ = (κ_d − 2 d K_d)/(1 − d κ_d
+ d² K_d), K_d = (gᵀ adj(H) g)/|g|⁴, combine-first ordering, generalized
|denom| > ½ guard, K ≡ +0 in 2D (bit-identical, full-ladder regression), 3D
sphere gate O(h^1.95) with the K-less control failing at O(h^1.02). The
|g|-floor + linear fallback remains labeled a pinch-off/topology guard only.

Residual task: **analytic unit test on the torus** (κ1 ≠ κ2 with sign change,
saddle side K < 0), plus cylinder and sphere, both signs of d, exactness given
exact inputs, guard-branch activation, K = 0 reduction identity — as a
`leiaTest*` utility exiting nonzero on mismatch.

### WP5 — Normal-constancy: STATUS + guardrails (no new force-path work)

Status: on the faces the force actually uses (|snGrad α| > 0), the stabilized
foot-point re-referencing ALREADY delivers the required property — both active
faces in a normal column are mapped to the same interface value, i.e. κ_f is
single-valued along the normal within the ∇_f α support, to the delivery's
order O(h²). This is the architecture requirement, met.

Precision that the code and docs must keep (a coding agent will implement what
the words say):
- The mechanism is **measure the DISTANCE d_f at the foot, then algebraically
  invert the parallel-contour relation on the face-combined κ_f**.
- It is NOT "fetch κ at the projected point." Fetching κ at the foot is the
  normal extension — measured dead: O(h^1.08)/O(h^1.13) statics (constant
  Hessian pins the fit's curvature to the stencil centre), and 40–85× κ
  amplification in the coupled loop. Distance from the foot: yes. Curvature
  from the foot: no.

Residuals and where they are handled:
- Remaining κ_f variation along a column is O(h²) estimator scatter (part of
  the surviving cos4θ seed) → attacked upstream by WP1–WP3, not by more
  delivery machinery.
- In 3D the property held only with the Gaussian term → WP4 (done), not this WP.

DO NOT: parent/column delivery (every band cell inherits a parent interface
cell's κ). Statically it equilibrates; in the coupled loop it is measured
unstable (parent switching → force discontinuous in time; ray-distance parents
and temporal κ relaxation both fail — scope: measured in the parent-delivery
context **[v0.2 scope note]**).

Optional (low priority): a cell-centred re-referenced field κ_c^Γ using the
per-cell foot distance d_c from the cell's own fit — as a DIAGNOSTIC (map of
normal variation of κ across the support, before/after WP1–3) and as the
delivery layer if a v2 integral estimator is ever built. Not wired into the
force.

### WP6 — Gated profile reset (event-driven redistancing) — **[v0.2: PROMOTED, frozen-band variant runs in parallel with WP1]**

Rationale: the reset attacks the transported-ψ profile feedback directly (the
exponent), which no delivery-accuracy lever moved — and post-horizon-result it
is one of only two levers with a plausible exponent claim. **[v0.2] The
frozen-band variant does NOT need the displacement-sign gate**: with the
transported sign-change band + guard layers kept as Dirichlet anchors
(REDIST_FREEZE + REDIST_ANCHOR_LAYERS — the gated Eikonal already achieved
E_geom = E_vol = 0 exactly per event with this pattern) the interface is
untouched by construction; only the bulk/fringe profile that feeds the stencils
is refilled. The sign-histogram gate below stays reserved for any future
IN-BAND rewrite variant.

**Gate experiment for in-band variants only (unchanged from v0.1):** one reset
applied to an exact SDF: ψ ← d(feet from each cell's own fit). Measure zero-set
displacement vs h AND the sign histogram, plane-based vs quadratic-based.
- Sign-symmetric O(h³) (quadratic): per-step cadence licensed. One-signed at
  ANY order: per-step dead; event-driven only if per-event × count converges
  below transport order. One-signed O(h²) (expected for plane): confirms the
  retirement.

**[v0.3] BINDING METHODOLOGICAL RULE (user, 2026-08-07):** psi advection is
never modified, suspended, or specialized for the static test case — the full
semi-Lagrangian transport runs every step exactly as for a moving interface.
Any reset/stabilization candidate must be machinery that works identically for
a deforming, translating interface. PROMOTION of any candidate therefore
requires, beyond the static-droplet score: (1) the kinematic transport suite
(reversed vortex, 3D shear/deformation) with the candidate ACTIVE and the §0
transport orders unchanged; (2) a coupled gate in which the interface moves
through cells under capillary coupling (translating droplet,
transISTKinematic-type) with no degradation relative to the static verdict.
This applies to the frozen-band reset, the psi-filter, and every WP7 arm.

Frozen-band arm (runs after the WP0 diagnostic + a cheap static event
characterization):
1. Event trigger: WP0's A-spectrum leaving a healthy window (NOT a fixed
   cadence, NOT minGradPsiBand). Band cells with |ψ| ≲ h keep transported
   values (frozen-anchor pattern); the fill rewrites bulk/fringe only.
2. Reset values are geometric distances from the foot Newton (seed + 1–2
   re-projections). NEVER route the reset through the parallel-curve inverse
   (in 3D the scalar-inverse error would be baked into the master field).
3. Score in WP7 against the filter on t_blow(N) scaling — the claim worth
   having is a changed exponent, not a delayed onset.

### WP7 — Evaluation matrix and non-regression gates

1. **Offline (WP1):** arm table as specified; decision point for WP2.
2. **Static droplet N=128, arms:** baseline / WP2-fit / frozen-reset /
   fit+reset / WP3-filter (reference). Scored on t_blow and on the A-spectrum
   (WP0), not on minGradPsiBand. **[v0.2] Horizon 0.3, not 0.1** — the 0.1
   gate is measured too short to distinguish delay from closure.
3. **Exponent claim:** t_blow(N), N = 64–512, for the winning arm(s). Compare
   against the refit onset law (§0). A changed exponent is the publishable
   result; a shifted prefactor is a tuning result.
4. **Non-regression, every WP:** transport-order table of §0 unchanged;
   statics ≥ O(h^2.04) (2D) and ≈ O(h^1.95) (3D); serial ↔ np=4 to ~1e-12;
   **np4 ↔ np8 coupled agreement [v0.2]**; `offsetCorrection none`
   precondition still enforced.

---

## 3. Measured dead ends — DO NOT implement (with reasons)

1. Normal extension of κ (fetch κ at Newton/closest-point foot): no static
   gain (O(h^1.08)/O(h^1.13)); coupled 40–85× amplification.
2. Harmonic (Laplace) κ smoothing: correlates signed errors across faces into
   a coherent force — worse.
3. Kang/GFM interface-weighted κ_f in the coupled loop: statics 58× better,
   blows earlier (0.07 vs 0.44).
4. CST integral model as-is: first equilibrium at N=64, N=128 blows at 0.047.
   General lesson: integral reformulations smooth corrupted geometry; they do
   not create missing signal — the inputs are still fitted from ψ.
5. Per-side-then-combine parallel-curve inverse: O(h^1.40) (2D); confirmed
   ordering-sensitive. Combine first.
6. |g|-floor / linear fallback as an instability fix: the denominator never
   became small at the spikes. Keep only as pinch-off guard, labeled as such.
7. Plane-based band rewrite redistancing (any cadence): O(h²κ) one-signed
   compounding — retired with measurements. Quadratic IN-BAND rewrite only via
   the WP6 sign-histogram gate; the frozen-band variant is exempt (no in-band
   rewrite) **[v0.2]**.
8. Basis rotation of the full quadratic fit: mathematically a no-op
   (rotation-invariant function space). Harness sanity arm only.
9. Slope limiting (Venkatakrishnan/BJ) of the GEOMETRY fit: limiters kill H
   exactly where κ is needed. Geometry fit stays unlimited with a
   conditioning fallback ladder; limiting, if ever, is a transport-object
   question only.
10. Velocity extension **for pure advection** [v0.2 scoped]: worst accuracy at
    12–27× cost. (Its potential role as coupled-loop PROFILE MAINTENANCE is
    unmeasured — not a dead end, but out of scope unless WP1 + WP6 both fail.)
11. Any `offsetCorrection` other than `none` combined with
    `stabilizedFootPointFace`: double-correction (solver already refuses —
    keep the refusal).
12. **[v0.2]** The K-less scalar inverse in 3D: measured O(h^1.02) at raw
    magnitude on the sphere gate — never reintroduce it as a "simplification".

---

## 4. Repo conventions the agent must follow

- New models are runtime-selectable (RTS) classes in the existing pluggable
  hierarchies; never modify existing model classes to add behavior.
- All runs config-driven (`config/*.yaml`); every figure/table regenerated from
  curated CSVs by script; write-time-only reporting cadence preserved so wall
  clocks stay comparable.
- Standalone physics tests as `leiaTest*` utilities (pattern:
  `leiaTestMeanCurvature`, `leiaSetFields`) before coupled runs.
- Degeneracy/conditioning guards must be **centered and scale-free** (the
  absolute-coordinate determinant guard silently zeroed the N=512 band once).
  Add a regression test for any new guard.
- **[v0.2] Seam consistency:** any face value computed from rank-local fits
  must be synchronized across processor seams (swap-and-average, syncTools) —
  the unsynchronized variant cost 7–10× decomposition-dependent noise. The
  regression is np4 ↔ np8 agreement of a coupled run, NOT only serial ↔ np4.
- Every new field/diagnostic: serial ↔ np=4 agreement check.
- Two-object discipline everywhere: transport reconstruction vs geometry fit
  are separate objects with separate requirements (transport: single-valued,
  possibly limited; geometry: unlimited, conditioning-monitored).

---

## 5. Open questions (tracked, not blocking)

1. WP1.7 — 3D anisotropic degree along ζ: what is affordable and conditioned
   on 3×3×3 and on polyhedral stencils? Fallback: filter carries more load in 3D.
2. Normal source for the interface frame: LLS plane vs lagged g/|g| —
   WP1 arm (c) decides; tie-break on conditioning and on the ε²-scaling.
3. WP0 calibration constants on unstructured meshes (divided-difference form).
4. Whether the A-spectrum should also gate the geometry-fit fallback ladder
   (switch fits per-cell) — revisit after WP2 data.
5. κ_c^Γ diagnostic field (WP5 optional): build only if WP7 arm analysis needs
   the normal-variation map.
6. **[v0.2]** Which wavelength carries the FILTERED runaway (WP0 retrodiction
   (b)) — the single most retargeting-relevant unknown in the program.

---

## 6. WP0 retrodiction (a) — measured 2026-08-07, N=128 arithmetic baseline

Run: stationaryDroplet2D, arithmetic delivery, np 4, WP0 spectrum on;
FPE blow-up at t = 0.0803 (10 823 steps; the historical serial measurement was
0.105–0.11 — chaotic variance across decompositions).

Onsets (first 10x-floor crossing; floors from the first 1000 steps):

| observable | onset t |
|---|---:|
| minGradPsiBand < 0.95 (slow drift) | **0.028** |
| max|U| | 0.0655 |
| A2h band L2 | 0.0763 |
| kErrL2Band | 0.0784 |
| A4h band L2 | 0.0791 |
| A8h band L2 | never (floor 5.6e-6) |
| blow-up (FPE) | 0.0803 |

VERDICT vs the WP0 hypothesis:
1. **A_band is NOT a leading indicator of the growth phase** — the velocity
   onset precedes A2h by ~11 ms. The exponential growth (t ≈ 0.04–0.065,
   max|U| 3e-4 → 1e-2) happens with the 2h/4h/8h band amplitudes AT THEIR
   FLOORS and kErrL2Band nearly flat (79 → 104 1/m): the loop's early carrier
   is low-mode/coherent (consistent with the m=1-like drift seen in the
   filtered runs), and the grid-scale mode explosion is the ENDGAME
   (t > 0.07), not the driver. max|U| itself is the most sensitive integrator
   of the feedback.
2. A2h still leads the BLOW-UP by ~4 ms (~600 steps at N=128) — usable as a
   last-resort event trigger, but a reset fired there may be too late.
3. The ground-truth statement "minGradPsiBand does not lead" needs refinement:
   it does not predict the κ spikes locally (|∇ψ| ≈ 0.9 at spike cells), but
   its band-min DRIFT leads everything (0.96 → 0.91 over t = 0.02–0.06) — as a
   slow profile-health trigger for the FROZEN-BAND reset it is the earliest
   available signal, firing decades of steps before the endgame.

CONSEQUENCES for the program:
- WP6 trigger design inverts: primary trigger = the slow profile drift
  (minGradPsiBand or gradPsiL2ErrorBand window), NOT the A-spectrum;
  the A-spectrum stays as the endgame/carrier diagnostic and the WP3/WP7
  carrier-identification instrument (retrodiction (b) pending).
- WP1's target is unchanged (the diagonal-stencil aliasing is measured and
  real) but its expected effect shifts: it should weaken the ENDGAME
  amplification; whether the low-mode growth phase also feeds through the
  fitted H at sub-floor mode amplitudes is exactly what the WP7 fit arm will
  reveal.
- Figure: docs/method-comparison/.../figures/wp0_retrodiction_N128.png
  (workflow/scripts/make_wp0_retrodiction_fig.py); trace in
  runs/wp0Retro128 (gitignored).
