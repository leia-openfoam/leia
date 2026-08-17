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

## 7. WP0 retrodiction (b) — measured 2026-08-07, N=256 filtered, t -> 0.3

Run: stationaryDroplet2D, stabilized foot-point delivery + psiFilter
biharmonicBand theta = 0.05, np 8 (seam-synchronized), WP0 spectrum on.
Blow-up at t = 0.167 — reproducing the pre-instrumentation rerun to three
digits (deterministic).

Onsets (first 10x-floor crossing; floors = median of the first 2000 steps):

| observable | onset t |
|---|---:|
| max|U| | **0.0677** |
| minGradPsiBand < 0.95 | 0.0980 |
| kErrL2Band | 0.1159 |
| A2h band L2 | 0.1270 |
| A4h band L2 | 0.1417 |
| A8h band L2 | 0.1536 |
| blow-up (FPE) | 0.1672 |

Also: the filtered N=128 fuse datum (seam-fixed rerun, wall-limited at
t = 0.191): onset max|U| > 0.1 at t = 0.163. Filtered onsets therefore
still shrink under refinement (0.163 at N=128 vs 0.113 at N=256, ~1.4x per
doubling — milder than the unfiltered ~3.3x, but the trend is unchanged).

VERDICT:
1. Same ordering as the unfiltered N=128 case: max|U| leads EVERYTHING.
   During the whole growth phase (t = 0.05 -> 0.13, max|U| 6e-4 -> 1.7e-1)
   the band-profile second differences sit at their floors (A2h
   8.7e-8 -> 8e-7 m) and the band |grad psi| is healthy (0.993 at t = 0.05,
   0.945 at t = 0.10). The grid-scale profile corrugation and the |grad psi|
   collapse are consequences of the strong late flow, not its cause.
2. Carrier identification: NO wavelength of the band psi-profile carries the
   filtered growth phase — the mode amplitudes explode only in the endgame,
   in the order 2h, then 4h, then 8h. The earlier "migration to longer
   wavelengths" hypothesis is not supported for the paired combination.
3. The hard consequence: the growth phase is driven by the RESIDUAL SMOOTH
   spatial variation of the delivered kappa_f — present even at h^2 accuracy
   on a healthy signed-distance field. With the loop gain rising under
   refinement, no fixed-order accuracy improvement can close it for all N;
   this elevates the deck's Defect-1 argument (a face-constant kappa_f is
   unreachable for any interface shape, so pressure can never absorb the
   force error exactly) from statics to the measured driver of the coupled
   growth phase.

PROGRAM RETARGETING:
- WP1 remains justified as the ENDGAME breaker (the sawtooth amplification is
  real and kills the runs), but it cannot be expected to change the
  growth-phase exponent.
- The frozen-band reset (WP6) attacks bulk profile drift — measured here to
  LAG the velocity onset in the filtered pairing; it may extend the fuse but
  the same reasoning as for the filter applies. Keep, but with damped
  expectations; the moving-interface gates (v0.3 rule) still apply.
- The lever class that now matters for the EXPONENT: force-side
  reformulations that let the pressure absorb the residual kappa_f variation
  exactly (discrete-gradient / potential-form capillary flux), or per-N loop
  damping that scales with the gain. This is a NEW work package to be
  designed (WP8 candidate) — not covered by WPs 0-7.
- Figure: docs/method-comparison/.../figures/wp0_retrodiction_N256filtered.png.

## 8. WP8.0 — measured 2026-08-07: the projection is converged; the residual is structural

Instrumentation (`applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H`, shared
by both two-phase solvers): at the final pressure corrector, the residual
capillary face flux R_f = phig - p_rghEqn.flux() is measured on the ACTIVE
faces (|snGrad(alpha)| > 0) against phig itself, together with the pressure
solver's iterations and initial/final residual and the velocity residual
fvc::reconstruct produces from it. Per-corrector rows in
`capillaryFluxResidual.csv`.

Four arms, N=128 stationary droplet to t = 0.02, np 4:
delivery {arithmetic, stabilized foot point} x pressure algebra {production
GAMG (1e-9, relTol 0.01, 4-17 iterations/solve), strict DICPCG (1e-11,
relTol 0, ~265 iterations/solve)}.

| arm | R_f/phig at t=0 | median over run | max\|U_res\| at t=0 | max\|U\| peak |
|---|---:|---:|---:|---:|
| arithmetic, GAMG | 1.49e-3 | 4.13e-5 | 2.49e-3 | 4.44e-3 |
| arithmetic, PCG 1e-11 | 1.49e-3 | 4.13e-5 | 2.49e-3 | 4.44e-3 |
| foot point, GAMG | 1.80e-4 | 3.89e-5 | 3.04e-4 | 1.10e-3 |
| foot point, PCG 1e-11 | 1.80e-4 | 3.89e-5 | 3.04e-4 | 1.10e-3 |

FINDINGS:
1. **The pressure projection is not the constraint.** Tightening the pressure
   solve by ~30x in final residual (2.9e-10 -> 8.6e-12) and ~20x in iteration
   count changes the non-absorbable flux fraction, the velocity residual and
   max|U| by at most 2.4e-4 RELATIVE over the whole run. The residual R_f is
   therefore the genuine non-gradient component of the capillary flux, not an
   unconverged projection. (This does not contradict the 2026-07-28
   perturbed-mesh gate where strict PCG gained 18-29x: that was a NON-
   ORTHOGONAL mesh, where the Laplacian is harder and the solver residual can
   exceed the structural residual. On uniform orthogonal hex it cannot.)
2. **The pressure solve absorbs 99.85-99.98% of the capillary flux**, and
   fvc::reconstruct is proportional, not amplifying (velocity-residual
   fraction 1.4e-3 vs flux-residual fraction 1.5e-3 at t=0). The parasitic
   velocity is large only because the absorbed quantity is large: the
   capillary predictor velocity rAU*reconstruct(phig/rAUf) is ~1.7 m/s at
   N=128, so 1.5e-3 of it is 2.5e-3 m/s.
3. **The delivery's static advantage is transient in the coupled run.** The
   foot-point delivery starts 8.3x better in the non-absorbable fraction
   (1.80e-4 vs 1.49e-3) but by t = 0.02 the gap is 1.78x, and the median over
   the run is only 6% (3.89e-5 vs 4.13e-5). Once the flow perturbs psi, both
   deliveries' kappa_f error is dominated by the response to the perturbed
   field, not by the offset the re-referencing removes. This is independent
   confirmation of the WP0 retrodiction verdict and further deprioritizes
   accuracy work as an exponent lever.

CONSEQUENCES:
- Mechanism (i) of the WP8.0 hypothesis set is confirmed: the driver is the
  non-gradient part of G_sigma, irreducible by solver work or by delivery
  accuracy. The remaining lever classes are (a) force-side structural -- make
  more of the flux lie in the range of snGrad (element-DOF delivery with one
  curvature per interface element, or a ghost-fluid pressure-jump assembly);
  (b) gain-side -- semi-implicit capillarity.
- Pressure algebra can be left at production GAMG for orthogonal-mesh droplet
  studies (no accuracy penalty, ~20x cheaper per solve); the strict-PCG
  requirement stands only for non-orthogonal meshes.

## 9. WP8.1 — measured 2026-08-08: the across-support variation regrows and dominates

Instrumentation (`capillaryDriverSplit.H`, per step in the droplet metrics):
kappa_f is recovered from the delivered flux itself
(Q_f = G_sigma,f/(snGrad(alpha)_f |Sf|)), the interface normal is taken from
the SAME quadratic fit the curvature comes from (g/|g| via fitDerivatives --
the linear-plane normal would carry an orientation-correlated O(h) error into
a threshold-based classification; it is retained only as a fallback and as a
measured disagreement), and for every pair of active faces of a cell the
difference quotient is binned by the alignment of the pair separation with
that normal (> 0.9 ACROSS the support, < 0.3 ALONG the interface).

N=128 stationary droplet, np 4, both arms to blow-up:

| t [s] | arithmetic across / along [1/m] | foot point across / along [1/m] |
|---|---|---|
| 0.000 | 967 / 179 (5.4) | **17.0 / 12.8 (1.3)** |
| 0.005 | 968 / 276 (3.5) | 154 / 31.8 (4.9) |
| 0.010 | 990 / 247 (4.0) | 280 / 46.8 (6.0) |
| 0.020 | 1076 / 252 (4.3) | 491 / 92.0 (5.3) |
| 0.040 | 1201 / 276 (4.4) | 733 / 155 (4.7) |
| 0.060 | 1346 / 291 (4.6) | 2969 / 730 (4.1) |
| blow-up | 0.0803 | 0.0668 |

FINDINGS:
1. **Outcome 1 of the WP8.1 hypothesis set.** The across-support variation of
   the delivered curvature -- the component that is unphysical by
   construction, since the exact curvature is constant along interface normals
   within the support -- dominates the along-interface component by a factor
   4-6 for BOTH deliveries throughout the growth phase. The two components
   only equalise in the final decade before FPE, where everything diverges.
2. **The foot-point delivery's structural advantage is destroyed by the
   coupling, fast.** It starts 57x below the arithmetic delivery in the
   across-support component (17.0 vs 967) and reaches the same magnitude
   within ~0.05 s: x9 in the first 5 ms, x43 by t = 0.04. The re-referencing
   removes the parallel-contour offset of a CLEAN field; it does not prevent
   across-support variation from regenerating once the transported psi
   deviates from a signed distance, because every ingredient (foot distances,
   fit gradients, the |1-d kappa| guard branch) is recomputed from that psi.
3. The normal used for the classification is not a source of doubt: the
   quadratic-fit normal and the linear-plane normal agree to 0.04-0.07 degrees
   through the entire quiet phase, diverging to ~4 degrees only in the last
   milliseconds. The correction was methodologically necessary but numerically
   minor here.

CONSEQUENCE (this settles the WP8 ordering):
A delivery that attaches ONE curvature degree of freedom to each interface
element and extends it along that element's normal sets the across-support
component to zero IDENTICALLY, at every time step, regardless of how psi
degrades -- it is a structural property of the delivery, not an accuracy
property of the estimator. Since that component is measured to be the
dominant part of the driver in the coupled state, and since the accuracy-based
route demonstrably loses its advantage within milliseconds, the element-DOF
delivery is now the justified next construction. The machinery exists
(connectedInterfaceCurvature fills kappaInterfaceFace, consumed via
faceCurvatureSource registered); what it needs is the accurate estimator
(foot-point re-referenced quadratic) in place of its tangential chain fit,
and the moving-interface promotion gates of the v0.3 rule.

Figure: docs/method-comparison/.../figures/wp81_driver_split_N128.png
(workflow/scripts/make_driver_split_fig.py).

## 10. WP8.2 — measured 2026-08-08: the cut-cell delivery falsifies the WP8.1 consequence

The construction §9 called for was built (`curvatureExtension cutCellFootPointFace`,
commit 44c9a84): one curvature per CUT cell (0 < alpha < 1) — that cell's own
quadratic-fit curvature at its centre, inverted to the interface with that
cell's own stabilized foot-point distance through the parallel-surface inverse
(Gaussian term included, so 2D and 3D are one code path) — assigned to every
active face of the cell, the mean where both sides are cut. No `geometricD`
branching anywhere: cut cells, `fitDerivatives`, `cof(H)` and the foot search
are dimension-independent, and K is exactly +0 on a pseudo-2D fit.

STATIC GATES (Lichtenberg 54005348, serial). It is the most accurate delivery
in the table:

| gate | N ladder | order (L2) | L2 at the finest N | best competitor |
|---|---|---|---|---|
| 2D circle | 32-512 | h^2.00 | 0.0761 1/m (N=512) | per-face inverse 0.105, h^2.04 |
| 3D sphere | 32-128 | h^1.95 | 0.00518 1/m (N=128) | per-face inverse 0.00518, h^1.96 |

It beats the per-face inversion at EVERY 2D resolution (1.3-1.4x lower L2) and
ties it in 3D, with zero foot-search failures and zero active faces left
without a cut-cell side at every resolution in both dimensions.

FORCE STRUCTURE (coupled, N=128, np 8). It does what §9 asked for, and keeps
doing it:

| quantity | per-face inverse | cut-cell inverse |
|---|---|---|
| across-support driver L2 at t=0 [1/m] | 16.96 | 2.35 (7.2x lower) |
| across/along ratio at t=0 | 1.33 | 0.41 |
| across/along ratio at 0.5 t_blow | 5.34 | 0.97 |
| max abs U at the first step [m/s] | 2.83e-4 | 1.48e-5 (19x lower) |

The across-support component never regains dominance — the structural claim of
§9 holds for the whole run, not only at t=0.

COUPLED VERDICT (Lichtenberg 54005349/54005350). The droplet blows up SOONER:

| delivery | E_L2 at N=128 [1/m] | max abs U(t=0) | t_blow N=128 | t_blow N=256 | late growth rate [1/s] |
|---|---|---|---|---|---|
| arithmetic | 45.8 | 2.31e-3 | 0.0803 | — | 186 |
| per-face inverse | 1.67 | 2.83e-4 | 0.0668 | 0.0348 | 189 |
| cut-cell inverse | 1.18 | 1.48e-5 | 0.0202 | 0.0145 | 454 |

Both cut-cell runs ended in a floating-point exception with max abs U past
1e50 and the time step collapsed — a genuine blow-up, not a code fault. Growth
is smooth exponential in every arm (median step-to-step change of max abs U
1.4%, no step above 50%), so this is not a switching artefact of the
membership of the cut-cell set.

CONSEQUENCE 1 — the WP8.1 consequence is FALSIFIED. Removing the across-support
variation of kappa_f, which §9 measured to be the dominant part of the driver,
made the instability faster, not slower. The variation of the delivered face
curvature is therefore NOT what sets the growth rate: at t = 0.0103 the
cut-cell run carries roughly TEN times less total delivered-curvature variation
than the per-face run at the same physical time and has twice the velocity.
No further work should be justified by "it reduces the spurious component of
the capillary driver" — that lever is now measured and it is the wrong one.

CONSEQUENCE 2 — accuracy on an exact psi is ANTI-correlated with coupled
stability. Over the three deliveries above, static L2 falls 45.8 -> 1.67 ->
1.18 while t_blow falls 0.0803 -> 0.0668 -> 0.0202. Every static gate in this
program has been scoring deliveries on a number that does not order them the
way the coupled run does.

THE NUMBER THAT DOES TRACK IT: the gain, how far kappa_f moves per unit
perturbation of psi on the force support (new utility
`applications/test/leiaTestCurvatureNoiseGain`, run by the gate study on every
case; curated table `data/tables/curvature_gain.csv`). Reported as G h^2, the
amplification relative to the h^-2 that any curvature operator must pay:

| delivery | E_L2 (N=128) | G h^2 | late growth rate [1/s] |
|---|---|---|---|
| arithmetic | 45.8 | 0.644 | 186 |
| per-face inverse | 1.67 | 0.643 | 189 |
| cut-cell inverse | 1.18 | 0.821 | 454 |

The gain is constant to ~3% over the N = 64-256 ladder and linear in the
perturbation amplitude over a 100x range (eps = 1e-3 to 1e-1), so it is a
property of the delivery, not of the mesh or of the probe. The two deliveries
with equal gain (0.2% apart) have equal growth rate (1.6% apart) despite a 27x
difference in accuracy; the delivery with 28% more gain grows 2.4x faster. The
gain therefore ORDERS the deliveries correctly but does not predict the
magnitude — a 28% gain increase buying a 2.4x rate increase is the signature of
a loop operating near its stability threshold, and that nonlinearity is
recorded as an observation, not a model.

The mechanism the numbers support: the per-face inversion computes an
independent inversion per face, so a cell's four (2D) active faces average four
partly independent foot-distance errors; the cut-cell delivery puts the whole
cell's force on ONE inversion of ONE foot distance d_c, whose error is
amplified by the inverse's 1/(1 - kappa_d d + K d^2) and then applied to every
face with full weight. On an exact signed distance field d_c is exact and the
concentration is pure gain in accuracy; under transport it is pure gain in
noise.

WHAT THIS OPENS: the cheap offline predictor the program has been missing. Any
candidate delivery can now be scored in seconds, serially, on a static circle,
by two numbers — its convergence order AND its gain — and only those that do
not raise the gain are worth a coupled run. The next constructions should be
selected to LOWER G h^2 below 0.643 (more averaging, not less), which is the
opposite of the direction WP8.1 pointed.

## 11. WP8.3 — measured 2026-08-10: what the inverse assumes, and the cell-mean delivery

### 11.1 The parallel-surface inverse does NOT assume a signed distance

Verified in code and by derivation. With g = grad psi, H = Hess psi, beta = |g|,
n = g/beta, P = I - n n:

  kappa_d = (tr(H)|g|^2 - g.(H.g))/|g|^3 = tr(S),  S = P H P/beta
  K_d     = (g.cof(H).g)/|g|^4           = det(S)
  d       = root-found Euclidean distance to the fitted zero surface

None of the three uses |grad psi| = 1. footPointDistance
(src/leiaLevelSet/.../uncachedQuadraticWeightedLeastSquaresReconstruction.C:865-983)
Newton-projects onto the cell's own fitted zero surface and returns the
query-to-foot vector on the UNIT model normal, so |d| is a length in metres, not
|psi| and not |psi|/|grad psi|. Under psi -> f(psi), f(0) = 0: K_d is exactly
invariant, kappa_d is equivariant (flips sign with f'), d flips too, and the
inverse is exactly equivariant -- so the delivered kappa_Gamma is invariant,
PROVIDED f' neither vanishes nor changes sign anywhere in the band. Checked
numerically on a torus (R = 1.0, r = 0.35) with f(t) = 2.5t + 7t^2 + 40t^3, so
|grad psi| in [2.5, 4]: the inverse recovers kappa = 1/r + cos(th)/(R + r cos(th))
to 7 digits at both K > 0 and K < 0.

WEAKEST SUFFICIENT HYPOTHESIS: the level sets in the band are the PARALLEL
(offset) surfaces of Gamma -- equivalently psi = f(phi) with phi the exact signed
distance and f' > 0; equivalently a := (n.grad)n = 0; equivalently beta is
constant ON each level set (not equal to 1). That is strictly weaker than a
signed distance and no rescaling of psi can violate it.

Caveat: the invariance is a CONTINUUM statement. kappa_d and K_d are read off a
weighted least-squares quadratic of raw psi differences, and LS projection does
not commute with a nonlinear f. Measured (sphere, cell centre half a cell out,
N=128): delivered kappa_Gamma = 22.795 / 22.770 / 22.580 for phi, phi + 2phi^2,
phi + 20phi^2 against exact 22.727. The ORDER is reparametrization-robust (~h^2);
the CONSTANT is not. Fit conditioning, not the algebra, is why a near-signed-
distance psi is preferred.

### 11.2 The exact residual, and why no free version of it exists

Along an integral curve of n, for ANY smooth psi with g != 0:

  dkappa/ds = -(kappa^2 - 2K) + D,   D := div(a) = Lap_Gamma(ln beta)
                                            - |grad_Gamma ln beta|^2
                                          = -beta Lap_Gamma(1/beta)

The inverse integrates the Riccati part EXACTLY and drops D, so

  kappa_Gamma(inverse) - kappa_Gamma(true) = int_0^d D ds = d D|_Gamma + O(d^2).

Verified: psi = (|x|-R)(1 + 0.3 z/|x|) on R=1 gives err/d = -0.3928 for
d = 0.08 ... 0.005 (clean first order) against D = -0.3536 - 0.0392 = -0.3928.

CORRECTION TO AN EARLIER STATEMENT IN THIS THREAD: the residual was first written
as +Lap_Gamma(d), which has the WRONG SIGN and omits the -|grad_Gamma ln beta|^2
term. A correction built from that form would add the error instead of removing
it.

Three consequences:

(a) The residual is FIRST order in d. Non-parallelism does not degrade the
    correction gracefully -- it returns kappa_f to the uncorrected h^1.13/h^1.03
    regime.

(b) THE INVERSE CAN BE WORSE THAN NOT CORRECTING. Uncorrected error is
    -(kappa^2 - 2K - D)d; corrected error is +D d. So it improves kappa iff
    |D| < |kappa^2 - 2K - D|. On the 2:1 ellipse psi = x^2/A^2 + y^2/B^2 - 1 at
    the major-axis vertex (kappa = 2, D = 3, kappa^2 - 2K = 4): at d = 1e-3 the
    uncorrected error is 9.995e-4 and the corrected error 3.005e-3 -- the inverse
    is 3.01x WORSE than doing nothing there. At the co-vertex it still helps.
    Every static gate in this program is a circle or a sphere, where D = 0
    identically, so this failure mode has never been exercised.

(c) D contains n.grad(Lap psi) and n_i n_j n_k psi_ijk -- genuine THIRD
    derivatives. Two fields with identical (psi, g, H) at x_f but different third
    derivatives give different kappa_Gamma at O(d).

THE FREE-LOOKING VERSION IS DEGENERATE. For a quadratic model H is constant, so
D_fit is closed-form in (g, H) at no extra derivative order. But for an exact
signed distance H n = 0, whence

  D_fit = tr(H^2) - n.H^2.n = k1^2 + k2^2 = kappa^2 - 2K,   while  D_true = 0,

so subtracting d D_fit cancels the Riccati term exactly and returns the
UNCORRECTED value. This was measured independently on the exact circle before the
derivation: the projected Laplace-Beltrami of 1/|grad psi| computed from the fit
returns -kappa^2 where the truth is 0, converging to a ratio of -1.000 over
N = 32-512 at both grid-aligned and diagonal placement. The same objection kills
every same-fit two-point construction, including re-evaluating the fit's kappa at
a second station along its own normal.

A CUBIC FIT IS NOT MERELY EXPENSIVE, IT IS RANK-DEFICIENT ON THE CURRENT STENCIL:
with offsets in {-h, 0, h}, delta^3 = h^2 delta identically, so the 2D design
matrix is 8x9 with rank 7. It needs a stencil that does not exist here.

### 11.3 Where a signed-distance assumption DOES survive in the delivered force

Not in the curvature. The CSF face force is sigma kappa_f snGrad(alpha)_f, and the
alpha side locates the interface with a LINEAR least-squares fit and the
first-order offset d = pc[3]/nmag, i.e. psi/|grad psi|:

  applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H:154-161
  src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C:130-140

That estimator carries e_d = -(1/2) b_n d^2 + O(d^3) with b_n = (n.H.n)/beta --
exactly the offset error the curvature channel removed by root-finding. The
Hessian-corrected root d = 2 psi_c/(|g| + sqrt(|g|^2 - 2 psi_c h_nn)) already
exists as offsetDistance, and h_nn is already inside kappa_d, so the fix is ZERO
extra derivative order. It moves alpha, hence the force support and the cut-cell
set, so it must be re-gated on transport and volume conservation. Not yet
measured; recorded as the next structural candidate.

Minor, free: the Kang blend weights in the per-face path use |psi| as a length
(stabilizedFootPointFaceCurvature.H:130-131). Replacing them with the already
computed |d_O|, |d_N| makes the no-signed-distance claim exactly true. Predicted
effect on any measured number: none (the two foot distances agree to fit
truncation, moving kappa_Gamma by O(h^2-h^3)).

### 11.4 The gain mechanism, corrected: averaging versus concentrating

The earlier reading of the ladder -- that gain tracks derivative count -- is NOT
supported by the measurements:

  arithmetic -> perFaceInverse ADDS the cof(H) contraction and the foot-point
    Newton solve, and G h^2 is unchanged (0.6443 -> 0.6428, 0.2%) while static
    L2 changes 27x;
  cutCellInverse has IDENTICAL derivative content to perFaceInverse and 28%
    more gain.

The measured discriminator is how many independent estimates the delivered value
averages, not how many derivatives it takes. Concentrating one cell's force on a
single inversion of a single foot distance raises the gain; averaging n per-face
inversions lowers it.

### 11.5 The cell-mean delivery: the first arm to move the coupled behaviour

curvatureExtension cellMeanFootPointFace -- the per-face inversions, averaged over
each cut cell's active faces and assigned back to them.

| delivery | 2D G h^2 | 3D G h^2 | 2D L2 (N=512) | 3D L2 (N=128) |
|---|---|---|---|---|
| arithmetic | 0.638 | 0.823 | 11.35 | 0.4191 |
| per-face inverse | 0.639 | 0.826 | 0.1049 | 0.005254 |
| cut-cell inverse | 0.829 | 1.005 | 0.07608 | 0.005176 |
| cell-mean inverse | 0.402 | 0.534 | 0.08422 | 0.004093 |

Coupled stationary droplet (np 8, capillary time step):

| arm | max abs U(t=0) | t_blow [s] |
|---|---|---|
| arithmetic N=128 | 2.31e-3 | 0.0803 |
| per-face N=128 | 2.83e-4 | 0.0668 |
| cut-cell N=128 | 1.48e-5 | 0.0202 |
| cell-mean N=128 | 1.13e-5 | 0.1049 |
| per-face N=256 | 4.71e-5 | 0.0348 |
| cut-cell N=256 | 2.96e-6 | 0.0145 |
| cell-mean N=256 | 3.26e-6 | 0.0493 |

Growth rate of max abs U over MATCHED physical windows at N=128 [1/s]
(nan = that arm had already blown up):

| window [s] | arithmetic | per-face | cut-cell | cell-mean |
|---|---|---|---|---|
| 0.002-0.010 | -157 | -25 | 104 | -13 |
| 0.010-0.020 | -36 | -9 | 560 | 7 |
| 0.020-0.035 | -11 | 56 | 295937 | 31 |
| 0.035-0.050 | 21 | 134 | nan | 72 |
| 0.050-0.066 | 180 | 238 | nan | 100 |
| 0.066-0.080 | 280 | 9529 | nan | 59 |
| 0.080-0.100 | 386412 | nan | nan | 69 |

The cell-mean delivery is the LONGEST-SURVIVING of the four at N=128 (1.57x the
per-face inversion, 1.31x arithmetic, 5.2x cut-cell) and survives 1.42x longer than
the per-face inversion at N=256. It grows more slowly than every other arm in every
matched window from t = 0.02 on, and -- the qualitative difference -- it does not
enter the runaway acceleration the others do: its rate stays at 59-100 1/s from
t = 0.05 to 0.10 while per-face goes 238 -> 9529 and arithmetic 280 -> 386412.

MECHANISM, OBSERVED IN THIS RUN. The delivered curvature variation stays tiny for
most of the run and then grows together with the loss of the parallel foliation:

| t [s] | max abs U | across L2 | along L2 | kErrL2Band | minGradPsiBand |
|---|---|---|---|---|---|
| 0.005 | 1.75e-4 | 4.07 | 7.41 | 70.2 | 0.9975 |
| 0.020 | 4.45e-4 | 6.01 | 9.10 | 70.4 | 0.9983 |
| 0.050 | 2.73e-3 | 5.83 | 4.24 | 70.1 | 0.9881 |
| 0.080 | 2.13e-2 | 203 | 297 | 187 | 0.858 |
| 0.100 | 0.199 | 1721 | 1720 | 1100 | 0.724 |
| 0.1045 | 0.535 | 8095 | 10850 | 3308 | 0.256 |

The driver is flat at 4-9 1/m through t = 0.05 and then rises by three orders,
and it does so in step with minGradPsiBand leaving 1. That is exactly the
non-parallelism of sec. 11.2: the inverse's residual is the integral of
D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2 along the normal, which is zero
only while the level sets stay parallel. It is one run and therefore a correlation,
not a demonstrated causal chain -- but it is a testable one, and the measurement
that would settle it is the Pi = |Lap_Gamma ln beta|/(kappa^2 - 2K) sweep of
sec. 11.2(b) run on the coupled trace.

WHAT THIS DOES NOT ESTABLISH. It is not stability: the arm blows up, just later.
And the gain does not order everything -- arithmetic and per-face have the same
gain to 0.2%, yet arithmetic grows more slowly in every early window and outlives
it. Gain is a strong predictor, not a complete one.

A rate extrapolation from the two-point exponent d ln(rate)/d ln(G h^2) = 3.586
predicted a late rate of ~34 1/s for the cell-mean arm; the measured rate is
59-100 1/s. Direction right, magnitude over-predicted ~2-3x -- the exponent is a
local two-point fit and must not be extrapolated. The t_blow extrapolation made
from the truncated run (~0.12 s) came out 13% above the measured 0.1049 s.

OPEN, AND REQUIRED BEFORE ANY ACCURACY CLAIM IS GENERAL: the circle and the sphere
both have constant exact curvature, so the O(h^2 d^2 kappa) cost of the averaging
is identically zero on both static gates, and D = 0 there so 11.2(b) is never
exercised. A varying-curvature static gate (the perturbed circle
r = R[1 + eps cos(m theta)] harness, whose exact curvature is not the constant the
gate app assumes) is needed for both.

ACCEPTANCE CRITERION for any future "signed-distance-free" correction, checked in
the same table: G h^2 <= 0.65 AND fitted E_L2 order >= 1.9 at the Pi of the
coupled run. A candidate that meets only the second is, on this evidence, a
regression.

## 12. WP8.4 — measured 2026-08-12: the varying-curvature gate kills the cut-cell family

New gate: `config/faceCurvatureEllipse2D.yaml`, case `cases/ellipseDroplet2D`, a
2:1 signed-distance ellipse (semi-axes 2 mm x 1 mm in the 0.01 m box, N=64-512).
`signedDistanceEllipse` is a TRUE signed distance -- level sets stay parallel and
the foliation residual D of sec. 11.2 is still zero -- so this gate isolates ONE
thing: the exact curvature now VARIES along the interface,
kappa(theta) = ab/(a^2 sin^2 theta + b^2 cos^2 theta)^{3/2}, from 250 1/m at the
minor-axis vertex to 2000 1/m at the major-axis vertex, a ratio of 8. The gate app
takes every error against the LOCAL exact value queried from the surface at the
face centre. implicitSphere is left on the radius-based constant, since its
curvature() returns 1/R in 3D too; the circle and sphere gates are verified
BIT-IDENTICAL after the change.

RESULT — active-face L2 of kappa_f [1/m], and the fitted order:

| delivery | N=64 | N=128 | N=256 | N=512 | order (N>=128) | circle order |
|---|---|---|---|---|---|---|
| arithmetic | 102.7 | 52.38 | 28.57 | 13.66 | 0.97 | 1.13 |
| per-face inverse | 14.17 | 4.323 | 1.116 | 0.2785 | **1.98** | 2.04 |
| cut-cell inverse | 41.64 | 24.13 | 11.79 | 5.848 | **1.02** | 2.00 |
| cell-mean inverse | 43.44 | 24.27 | 11.81 | 5.851 | **1.03** | 2.05 |
| interface mean (control) | 472.2 | 497.4 | 509.1 | 505.6 | -0.01 | 2.03 |

BOTH one-value-per-cut-cell deliveries COLLAPSE TO FIRST ORDER, and are 21x less
accurate than the per-face inversion at N=512. The per-face inversion keeps
h^1.98. The circle gate showed the exact opposite ranking (2.00 and 2.05 against
2.04), and the interface-mean control -- exact on a circle by symmetry -- is
revealed as zeroth order, which is the same defect in its extreme form.

WHY, AND WHY IT IS STRUCTURAL. Assigning ONE value to every active face of a cut
cell forces kappa_f to be constant across a cell. Where the exact curvature
varies along the interface, the variation within one cell is O(h dkappa/ds), so a
face at the edge of the cell receives a value centred on the CELL, not on itself:
an O(h) offset, first order by construction. It cannot be repaired by averaging
better, because the defect is the lumping, not the estimator. Note that this is
NOT the same as symmetric smoothing: a symmetric average about each FACE stays
second order (its error is O(h^2 d^2 kappa)); it is the one-sided, cell-centred
assignment that costs an order.

The gain is unchanged by the geometry -- cell-mean 0.370-0.416, per-face
0.612-0.673, cut-cell 0.779-0.849 across the ellipse ladder, the same values as
on the circle -- confirming that G h^2 is a property of the delivery and not of
the shape.

CONSEQUENCE. The cell-mean delivery, adopted in sec. 11.5 on the strength of the
longest coupled survival (t_blow 0.1049 s against 0.0668 s), is FIRST ORDER on any
interface whose curvature varies and therefore cannot be used in a general case.
Its coupled advantage was measured on a CIRCLE, whose constant curvature hides the
defect exactly as the static circle gate did. The stationary-droplet test is blind
to this by construction. Both cut-cell-family deliveries are hereby retired from
consideration as production deliveries; they remain useful as the two extreme
points that established the gain-versus-accuracy dissociation.

THE ACCEPTANCE CRITERION OF sec. 11.5 IS AMENDED: G h^2 <= 0.65 AND fitted E_L2
order >= 1.9 ON THE VARYING-CURVATURE GATE, not on the circle. Under the amended
criterion the only delivery that passes today is the per-face inverse
(order 1.98, G h^2 0.647) -- which is the production default.

WHAT THIS OPENS. The direction that remains is to lower the gain WITHOUT lumping:
average each face's inversion symmetrically about THAT FACE (its interface-adjacent
face neighbours), so every face keeps a value centred on itself. Predicted to keep
h^2 while reducing the gain, and it is now cheap to falsify -- both gates plus the
gain are one command each. The cut-cell result says nothing against this variant;
it says the failure was the cell-centred assignment.

METHODOLOGICAL NOTE. This is the standing rule of the programme paying off: a
static test case must not be simpler than the interfaces the method is meant to
run on, or the method that wins is the one that exploits the simplification.
Sec. 10 said accuracy does not predict stability; sec. 12 says stability on a
circle does not predict accuracy on anything else. Both gates are needed, and the
coupled droplet needs a non-circular companion before any delivery is promoted on
coupled evidence alone.

## 13. WP7 refit — measured 2026-08-12: no delivery has changed the scaling

The programme's scoring criterion (sec. WP7) is whether an arm changes the
t_blow(N) EXPONENT, not its prefactor. Refit on the three deliveries that now have
two resolutions each (t_blow ~ N^p, capillary time step throughout):

| delivery | t_blow N=128 [s] | t_blow N=256 [s] | p |
|---|---|---|---|
| per-face inverse | 0.0668 | 0.0348 | -0.94 |
| cut-cell inverse | 0.0202 | 0.0145 | -0.48 |
| cell-mean inverse | 0.1049 | 0.0493 | -1.09 |

Every arm blows up SOONER in physical time as the mesh refines, and the cell-mean
arm's exponent is WORSE than the per-face one it beat on prefactor. Extrapolating
on the fitted exponents, its advantage decays:

| N | cell-mean [s] | per-face [s] | ratio |
|---|---|---|---|
| 128 | 0.1049 | 0.0668 | 1.57 |
| 256 | 0.0493 | 0.0348 | 1.42 |
| 512 | 0.0232 | 0.0181 | 1.28 |
| 1024 | 0.0109 | 0.0094 | 1.15 |
| 2048 | 0.0051 | 0.0049 | 1.04 |

So the cell-mean result was a PREFACTOR win that vanishes under refinement -- which
is what the WP7 criterion exists to catch, and a second independent reason (beside
its first-order accuracy, sec. 12) not to promote it. Solving for the resolution at
which each arm would survive to the 0.3 s horizon gives N < 50 for all three, i.e.
coarser than every mesh run: no delivery is on a path to the horizon, they all get
worse with refinement.

CONSEQUENCE FOR THE PROGRAMME. The curvature-DELIVERY lever is, on this evidence,
exhausted apart from one untested candidate: symmetric averaging about each face
(sec. 12), which is the only construction that could lower the gain without the
cell-centred lumping that costs an order. Everything else in the delivery space has
been measured. After that candidate is scored on the ellipse gate and the gain, the
programme should move to the two levers never yet touched:
  - the psi side: keep |grad psi| constant along level sets, so the parallel-
    foliation hypothesis the inverse needs is TRUE rather than corrected for;
  - the alpha side: the phase indicator's first-order psi/|grad psi| interface
    offset (sec. 11.3), which sits in snGrad(alpha), the other factor of the CSF
    face force, and has never been touched by any arm.
The foliation-residual diagnostic of sec. 11.2 run on the ALREADY SAVED psi fields
decides between those two at zero simulation cost, and should be done before either
is built.

SCORING NOTE, BINDING FROM HERE: an arm is not to be promoted on t_blow at a single
resolution. Two resolutions give a two-point exponent, which is the minimum, and
three (N = 128, 256, 512) are needed for the exponent to be a claim rather than an
estimate. That is the real cost of a decision and it is affordable for one or two
arms, not for a sweep -- which is exactly why the offline gates (order + gain,
minutes, serial) exist to reject candidates first.

## 14. WP8.5 — measured 2026-08-12: the delivery lever is closed, and psi degradation is the driver

### 14.1 Symmetric face averaging: the last delivery candidate, and it fails too

Sec. 12 left one construction untested: keep the per-face inversion and smooth it
over each face's OWN owner+neighbour active-face set, which is symmetric about that
face (owner and neighbour are mirror images through it). The prediction was that a
symmetric mean annihilates linear fields, so the smoothing error would be
O(h^2 d^2 kappa) and second order would survive while the gain still fell.

Measured on the 2:1 ellipse gate (order fitted on N >= 128):

| delivery | L2 at N=512 [1/m] | order | G h^2 |
|---|---|---|---|
| per-face inverse | 0.2785 | 1.98 | 0.647 |
| symmetric face-mean, theta = 0.5 | 2.014 | 1.10 | 0.445 |
| symmetric face-mean, theta = 1.0 | 4.009 | 1.07 | -- |
| cut-cell inverse | 5.848 | 1.02 | 0.818 |
| cell-mean inverse | 5.851 | 1.03 | 0.395 |

THE PREDICTION IS WRONG, and the reason is worth recording. The stencil is
symmetric GEOMETRICALLY but it is restricted to ACTIVE faces, and the active mask
depends on where the interface sits inside the cell. That asymmetry is O(1) in
which faces are included, so the mean is offset at O(h) and the order is lost
anyway. The variant is 3x better in absolute error than the cut-cell family and
does cut the gain by 31%, but it does not pass the amended criterion.

On the circle every delivery is second order (1.97-2.00), the third independent
confirmation that a constant-curvature gate cannot rank deliveries.

CONCLUSION. Within face-value averaging, gain and order trade against each other:
every construction that lowered G h^2 below the per-face value lost an order on a
varying-curvature interface. The curvature-DELIVERY lever is closed. The production
per-face inversion (order 1.98, G h^2 0.647) stands. One structurally different
variant remains unmeasured -- smoothing a CELL field, whose stencil does not depend
on the interface position -- and it is recorded as low expected value beside the
psi-side result below.

### 14.2 The foliation residual, measured on the saved coupled fields

`applications/test/leiaTestFoliationResidual` evaluates
D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2 (beta = |grad psi|) on written psi
fields through the volumetric-plus-projection identity, band-restricted to
|psi| < 2h. Run on the cell-mean N=128 coupled case (the arm that survived
longest, t_blow = 0.1049 s), against the band curvature error the solver itself
recorded. The decision rule was written into the app header BEFORE the run.

| t [s] | beta range in band | \|a_T\| L2 [1/m] | D L2 [1/m^2] | bias at h/2 [1/m] | bias at psi/beta [1/m] | measured kErrL2Band [1/m] |
|---|---|---|---|---|---|---|
| 0.00 | [1.000, 1.000] | 2.28 | 8.99e3 | 0.35 | 0.77 | ~70 |
| 0.02 | [0.9983, 0.9999] | 2.20 | 1.22e4 | 0.48 | 0.97 | 70.4 |
| 0.04 | [0.9966, 1.0036] | 4.27 | 2.63e4 | 1.03 | 2.38 | ~70 |
| 0.06 | [0.974, 1.025] | 30.6 | 2.24e5 | 8.73 | 22.7 | ~75 |
| 0.08 | [0.858, 1.089] | 164 | 1.39e6 | 54.5 | 131 | 187 |
| 0.10 | [0.724, 1.373] | 583 | 5.44e6 | 212 | 411 | 1100 |

READING, BY THE PRE-REGISTERED RULE. The rule was: same order as the measured
error and growing with it => the foliation residual explains the growth; one or
more orders below => spurious. The outcome is the first, with an important
qualification:

- EARLY the residual is negligible: 0.35-1.0 1/m against a 70 1/m error, under 1.5%.
  The early curvature error is the ordinary O(h^2) delivery error, exactly the
  static-gate value, and has nothing to do with the foliation.
- The residual then grows 600x (0.35 -> 212) while the total error grows only 16x
  (70 -> 1100). It is by far the FASTEST-GROWING term.
- By blow-up it is order-parity with the total: 212-411 against 1100, i.e. 19-37%.
  At t = 0.08 it is 29-70% of the 187 that the solver measured.

So the foliation residual accounts for roughly a quarter of the error at blow-up
and none of it at the start -- it is the term that takes over, not the term that
was always there. Decomposing the GROWTH rather than the total: of the +1030 1/m
increase from t = 0.02 to 0.10, the residual supplies +212, about 21%. The other
~80% is the fit error itself degrading, which is driven by the SAME cause: beta
spreading from [0.998, 1.000] to [0.724, 1.373] corrupts the least-squares
quadratic, not only the inverse that follows it.

CONSEQUENCE, AND THE NEXT LEVER. psi losing its signed-distance character is the
driver of the curvature-error growth, by two mechanisms that a single fix
addresses: the inverse's parallel-foliation hypothesis stops holding (measured,
21% of the growth) and the fit that feeds it degrades (the remainder). The lever is
therefore the OPERAND, not the operator: hold |grad psi| constant along the level
sets in the band. That is the only remaining candidate with a mechanism, a measured
magnitude, and a plausible route to the t_blow EXPONENT rather than its prefactor
(sec. 13) -- because it removes a term that grows 600x over a run rather than
rescaling one that does not.

CAVEAT, STATED PLAINLY. This is one run, one resolution, one delivery. The
attribution is a magnitude accounting, not a controlled experiment. The controlled
version is to hold the foliation parallel and re-measure: if the curvature error
stops growing and t_blow moves, the chain is demonstrated; if the error still
grows, the fit degradation is the whole story and the band renormalization only
buys the 21%.

## 15. WP9 — measured 2026-08-12: the psi-side lever closes too, and what that means

### 15.1 Band renormalization, coupled: 3.4x WORSE

Arm: `psiRenormalization { type footPointGradient; }` with the divisor smoothed
over each band cell's face neighbours and the operation gated on the band L2 of
the dimensionless non-parallelism |P H n|/beta * h (threshold 5e-3). Delivery held
at the production per-face inversion, so the psi side is isolated. N=128, np 8.

| arm | t_blow [s] |
|---|---|
| per-face reference (no renormalization) | 0.0668 |
| + band renormalization | **0.0195** |

The gate behaved exactly as designed -- while the foliation is healthy the run is
BYTE-IDENTICAL to the control (verified at N=64 on max abs U, min|grad psi| and the
band curvature error at every sampled time). It then fired 799 times, each firing
rescaling psi by up to 2% in a single step, and the droplet died at t = 0.0195.
The N=256 twin was cancelled at t = 0.0072 rather than burning 8 h on a
configuration already 3.4x worse on the prefactor; no exponent is reported for it.

### 15.2 Relaxation does not rescue it, and the reason is quantitative

Predicted remedy: spread the same total correction over many steps
(psi <- psi/beta_Gamma^omega), so the SMOOTH drift being corrected accumulates
coherently (omega*n) while the grid-scale fit error rides along incoherently
(omega*sqrt(n)) -- a sqrt(n) gain in correction-to-noise. Measured at N=64, gate
forced open, over 2 ms:

| arm | max psi kick | max abs U at 2 ms | min\|grad psi\| at 2 ms | kErrL2Band |
|---|---|---|---|---|
| control (no renormalization) | -- | 1.63e-3 | 0.9877 | 141 |
| omega = 1.00 | 5.96e-3 | 5.36e-3 | 0.9716 | 153 |
| omega = 0.20 | 1.19e-3 | 4.08e-3 | 0.9751 | 149 |
| omega = 0.05 | 2.97e-4 | 3.33e-3 | 0.9800 | 145 |

The kick scales EXACTLY linearly with omega, as designed. The harm does not: a 20x
smaller kick removes only 2.2x of the excess velocity. And the decisive number is
min|grad psi| -- the operation whose entire purpose is to hold |grad psi| near 1
makes it drift FASTER than doing nothing, at EVERY relaxation level tested.

WHY, quantitatively. beta_Gamma is estimated from the same quadratic fit, and its
error at N=64 is 0.3-0.6% (measured directly: on the EXACT signed-distance initial
field, where beta is identically 1, the fit returns beta_Gamma in [0.9936, 0.9972]).
The drift being corrected over this window is ~0.1%. The correction is therefore
dominated by the error of the estimator that computes it, at any omega. It can only
pay off once the drift greatly exceeds 0.5% -- which is exactly the late phase where
the loop amplifies any injection fastest.

### 15.3 The finding, across four levers

Four structurally different attempts have now been measured out:

| lever | what it did | outcome |
|---|---|---|
| cut-cell / cell-mean delivery | remove the across-support kappa_f variation | earlier blow-up (sec. 10); first order on varying curvature (sec. 12) |
| symmetric face averaging | lower the gain without lumping | first order anyway, active-mask asymmetry (sec. 14.1) |
| band renormalization | restore the parallel foliation | 3.4x earlier blow-up; drift made worse at every relaxation (this section) |
| accuracy of the inverse itself | remove the O(h) offset | 27x accuracy for 0.2% gain change; no stability effect (sec. 10) |

The pattern is one statement: THIS SYSTEM IS NOT LIMITED BY SYSTEMATIC ERROR IN ANY
OPERATOR. IT IS LIMITED BY ITS SENSITIVITY TO PERTURBATIONS OF psi. Every corrective
built from the same quadratic fit necessarily injects that fit's error into the
field, and the injection is amplified faster than the systematic error it removes.
Accuracy improvements are free (they do not change the gain) and stability
improvements bought by averaging cost an order. No arm has changed the t_blow
exponent (sec. 13).

CONSEQUENCE FOR WHAT TO TRY NEXT. If neither operator nor operand can be corrected
without injecting noise, the remaining target is the AMPLIFICATION PATH itself:
force -> velocity -> interface displacement -> psi -> fit -> force. The one element
of that loop never touched is the velocity the level set is advected BY. Advecting
psi with a normal-extended, interface-derived velocity instead of the raw cell
velocity removes the route by which grid-scale velocity noise reaches psi at all,
while leaving the momentum velocity untouched. The programme's own non-goal list
made this conditional on WP1 and WP6 both failing; both have now failed, so it is
licensed.

STATUS OF THAT ARM: the machinery is fully wired in this solver
(`velocityExtension::New`, `slVelExt->correct()`, and `slAdv->advect(psi, UextStar,
slVelExt->Uext())`), the default is `none`, and the token plus dictionary block are
now in the droplet cases. But all three models tried -- closestPoint, steadyUpwind,
meshWaveExt -- abort immediately with "failed lookup of alpha": they assume the
phase field is named `alpha`, while the two-phase cases name it `alpha.water`. That
is precisely the "unmeasured-for-maintenance" status the programme recorded for
velocity extension in the two-phase solver. It is a small fix -- make the phase-field
name configurable in the extension models -- but it is a CODE change, not a config
change, and it must be done before the arm can be measured.

## 16. The loop model and the time axis -- measured 2026-08-13/14: the growth rate decomposes, and every coupling lever is exonerated

### 16.1 The assembled loop model and the capillary-dt sweep

Assembling the programme's measured numbers into one causal chain -- delta_psi
(grid scale, amplitude eps) -> delta_kappa_f = 0.647*eps/h^2 (the measured gain)
-> CSF flux -> one explicit step -> interface displacement -> delta_psi -- gives
a per-step loop gain gamma ~ 0.647*sigma*dt^2/(rho*h^3), h-INDEPENDENT under the
capillary-scaled step. The discriminating test is a capillary-dt sweep at fixed
N (config/stationaryDropletDtSweep.yaml, N=128, coeff = dt_sigma/{4,8,16}; then
config/stationaryDropletDtScaling.yaml, N=64/256). Growth rates r are fitted
directly from ln(maxMagU) vs t in the window [3e-3, 3e-1] -- t_blow is NOT a
valid proxy (incubation and e-fold count vary; the e-fold count K = r*t_blow
runs 5.4..13.3 across arms).

    r(N, dt) = r0(N) + c(N)*dt          (per-N fits, <= 6% residual)

    N=64 :  r = 30.0/30.4/32.5  -> r0 ~ 31,  c ~ 0 (NO dt term)
    N=128:  r = 202.9/104.4/70.0 -> r0 = 20.8, c = 2.40e7 (dt term = 90%)
    N=256:  r = 249.5/181.8/125.6 -> r0 = 91.7, c = 6.07e7 (dt term = 65%)

The first halving at N=128 halves r within 3% -- the dt-linearity is real. But
the crude model's c ~ N^3 is falsified (measured N^1.3 from 128->256, absent at
64: the mechanism SWITCHES ON between N=64 and N=128, and omega_2h*dt is
N-independent at fixed coeff, so no naive per-mode explicit threshold explains
that). The floor r0 does not vanish and rises toward fine grids. Volume and
shape error at matched times IMPROVE monotonically with smaller dt (dt/16 at
N=128: volume 73x, shape 33x better) -- smaller dt is pure gain except cost.

Bonus falsifications: smaller dt = 2-4x MORE remaps and fits per unit time yet
LESS growth, so per-step noise injection is not the driver; at N=64 the maxU
amplitude at matched time scales LINEARLY with dt while r is flat -- there the
dt-dependence is all SEED (the quasi-static balanced-CSF residual ~ dt), not
rate.

### 16.2 Exonerations (each by direct experiment)

- ADAPTIVITY: fixed-dt re-run of the full 9-cell matrix (deltaT = coeff/N^1.5,
  adjustTimeStep no, writeControl runTime -- now the case standard per
  cases/stationaryDroplet2D/system/controlDict.template): 7/9 cells within
  +-3% on r, no systematic shift. The two +13/+21% outliers calibrate the
  seed-sensitivity of single-cell r fits.
- GAMG: p_rgh swapped to PCG+DIC at np4 reproduces the trajectory to 4-5
  significant digits in t_blow after ~9000 steps THROUGH a chaotic amplifier
  (probes p128_pcg/p64_pcg). The linear solver's internal path is irrelevant at
  our tolerances. (The docs' "parallel inconsistent" flag concerns GAMG as a
  PRECONDITIONER inside PCG -- a code path we do not use.)
- DECOMPOSITION: serial runs blow ~2x later at BOTH N (0.118 vs 0.065 at 128;
  0.367 vs 0.181 at 64) with r nearly unchanged at 128 (-5%) -- decomposition
  acts through the SEED, not the rate. It is NOT the linear solver (PCG
  reproduces the parallel trajectory exactly). Caveat: at N=64 serial r is +38%
  -- floor-rate fits carry that much decomposition/seed sensitivity.
- FORCE TIME-LAG (the implicit-Euler consistency argument): implemented as
  PIMPLE.psiOuterCorrectors (dictionary-driven, default no = byte-identical;
  createTransportFields.H + slAlphaEqn.H): every outer corrector restores
  psi^n and re-advects with the CURRENT velocity iterate, re-running band ->
  indicator -> kappa -> face delivery, so sigma*kappa_f*snGrad(alpha) converges
  to its t^{n+1} value within the step. Gates: switch-off bit-identical to the
  fixed-dt matrix (all rows, through blow-up); switch-on verified 3.00
  pipelines/step. RESULT: r = 202.8 vs 198.9 (+2%, noise). Coupling-depth
  bracket (config/stationaryDropletCouplingDepth.yaml): nOuter = 1/3/6, frozen
  or iterated -> r in [185, 203]. The ENTIRE outer-coupling axis is flat; the
  frozen-force lag is not the c*dt term. Six Picard passes bound any
  convergence shortfall out.
- THE x_d KERNEL: applications/test/leiaTestDeparturePoint drives the shipped
  slAdvection::advect with u = cos(omega t)*Omega ez^(x-c) (exact backward
  characteristic) and reads the feet off advected psi = x, y. Per-step foot
  error: order 3.00 at omega = 0; rel error = 0.07*(omega*dt)^2 (exact-supply)
  and 0.35*(omega*dt)^2 (AB2 substitution) -- clean quadratics with the
  theoretical constants. At the capillary operating point (omega_max*dt = 0.52
  at EVERY N) the AB2 path mislocates the stiffest mode's foot by 9.5% of its
  per-step displacement -- real, quantified, and NOT the growth carrier (the
  iterated arm used the true velocity iterate on passes >= 2; r unchanged).

### 16.3 The mode-resolved reframing

Fitting the SAME windows on the WP0 band-mode amplitudes and the curvature
error (columns already in the per-step CSV):

    r(maxU) ~ 1.5-2x r(A2h) at every cell; r(kErr) tracks r(A2h), not r(maxU).
    r(A2h) is far less dt-sensitive than r(maxU): at N=256 nearly
    dt-INDEPENDENT (160 -> 140 across 4x dt) while r(maxU) halves twice;
    r(A2h) = 22.6/15.1/13.9 at N=64, 99.7/68.8/51.4 at N=128.

The velocity observable grows at ~2x the rate of the psi corrugation
(superlinear response), and the "c*dt term" lives largely in maxU, much less
in the corrugation itself. THE ORDER PARAMETER OF THE INSTABILITY IS THE
CORRUGATION GROWTH RATE r(A2h) -- ~14-160 1/s across N=64..256 -- and it is
the quantity every so-far-successful lever (dt, filter, seed) only bought
prefactor against. The spatial loop psi -> fit -> kappa_f -> quasi-static
velocity response -> displacement, per unit physical time, is the one
mechanism left standing.

## 17. The foliation gate -- measured 2026-08-14: each ellipse arm isolates one term of the kappa_f error budget

config/faceCurvatureEllipsoidPsi2D.yaml runs the 2:1 ellipse with psi
initialized as the TRUE signed distance (parallel foliation, D == 0 -- the fit
error alone) and as the QUADRATIC FORM (x/a)^2+(y/b)^2-1 (beta varies by
a/b = 2 along the interface; psi exactly quadratic so the fit is EXACT and the
error is the inverse's foliation bias alone). Reference in both arms: the
zero-set ellipse's curvature at the closest point. Curated:
docs/method-comparison/method-comparison-article/data/tables/
face_curvature_orders_foliation.csv. Orders fitted on N >= 128 (L2, active
faces):

    model (foot-corrected)          SDF psi (D=0)      quadratic-form psi (D!=0)
    solverStabFootFace (PRODUCTION,
      the shipped header, bit-equal
      to quadraticCellCentre+foot)  0.278, order 1.98   8.69, order 0.94
    trHessian + foot                0.276, order 1.98   7.1e5, diverging
    stableFootPoint (GATE-ONLY
      foot-point-native variant)    4.95,  order 0.93   7.12, order 0.95
    quadraticNewtonFoot             13.7,  order 0.97   0.187, order 1.95

Findings:

1. CORRECTED 2026-08-14 (the first write-up of this section mislabeled the
   rows): the gate's `stableFootPoint` row is a GATE-ONLY experimental variant
   (per-side inversion at the cell centres, Kang-combined) that never calls
   the solver; the SHIPPED curvatureExtension stabilizedFootPointFace is
   interpolate(kappa) + ONE inversion at the face centre's own foot, and the
   new `solverStabFootFace` row -- which calls the actual solver header --
   is bit-equal to `quadraticCellCentre + foot` in both arms. THE PRODUCTION
   DELIVERY IS SECOND ORDER wherever the foliation is parallel (2.04 circle,
   ~1.95 sphere, 1.98 varying-curvature SDF ellipse) and FIRST order only
   where the foliation itself is non-parallel (0.94: the d*D model bias).
   The gate-only native variant is first order in both arms (tangential O(h)
   mislocation: its two per-side values belong to two different interface
   points); every coupled stability number of this programme was produced by
   the second-order delivery.
2. The error budget separates exactly as sec. 11 predicts: with D = 0 the
   foot-corrected plain quadratic is SECOND order (fit error only); with the
   fit exact and D != 0 it drops to FIRST order with the D*d bias alone --
   the inverse's model error measured in isolation for the first time.
3. quadraticNewtonFoot INVERTS between the arms: first order under the SDF
   (fit Hessian error dominates at its foot evaluation) but SECOND order,
   with a 46x smaller constant than the production delivery, under the
   quadratic form (exact fit + evaluation AT the interface point needs no
   parallel-surface correction at all). Evaluating at the foot beats
   correcting to the foot precisely when the foliation is the problem --
   the curvatureExtension closestPointNewton path already exists in the
   solver and becomes the natural candidate the moment the fit can be made
   exact enough (WP1/WP2 territory), with its noise gain still unmeasured.


### 17.4 footPointEvaluatedFace -- measured 2026-08-14: the foot-evaluated delivery on four geometries

Implemented as curvatureExtension footPointEvaluatedFace (dictionary-selected;
slReconstruction::curvatureAtFootPoint + computeFootPointEvaluatedFaceCurvature):
each adjacent cell finds the foot point OF THE FACE CENTRE on its own
quadratic's zero set and evaluates the fit curvature AT that point -- the value
belongs to the interface point itself, no parallel-surface conversion, no
level-set-spacing assumption; LINEAR face interpolation of the two per-cell
values (both estimate the SAME interface point); one-sided fallback; seam
sync. Gated through the solver header on every geometry (the shipped code
path, not an in-app variant). E_L2 over active faces at the finest N:

    geometry (psi)                       production      footPointEvaluatedFace
    circle, exact SDF (N=512)            0.105 (2.0)     10.5  (order ~1.0)
    2:1 ellipse, SDF (N=512)             0.278 (1.98)    12.7  (order ~1.0)
    2:1 ellipse, quadratic form (N=512)  8.69  (0.94)    4.8e-4 (EXACT: fit
                                                          reproduces psi, foot =
                                                          true closest point;
                                                          residual = Newton tol)
    oscillating geometry, a/b = 1.21,
    quadratic form (N=256)               6.24  (~1.07)   3.3e-4 (EXACT)

Noise gain: footEvalFace = perFaceInverse within 3% on both the circle
(0.657 vs 0.642 at N=512) and the oscillating geometry (ratio 1.03; on the
quadratic-form psi the *_DIMLESS column is scaled by |grad psi| ~ 2000 --
compare models within a case, not across psi scalings). By the equal-gain =>
equal-growth finding, the construction carries NO measured stability penalty
and NO stability gain.

VERDICT: the exact mirror of the production delivery, now measured on the
shipped path across four geometries. Evaluation at the foot is EXACT whenever
the fit reproduces psi (both quadratic-form initializations, including the
oscillating droplet's own t=0 field) and fit-error-limited (first order,
~100x the production constant) on true signed distances, where production is
second order. Neither dominates; the complementarity is now fully
quantified, and a per-face selection or blend between the two -- both already
computed from the same fit -- is the natural next construction, with the
arc-regression D-correction (sec. 17.3) as the alternative that upgrades the
conversion instead of switching the evaluation point.

## 18. The source, measured -- 2026-08-17: a neutral semi-discrete loop, destabilized by the explicit capillary coupling of its m=2 mode

Sections 15-16 closed every operator lever and left "the spatial loop per unit
physical time" as the suspect by elimination. This section replaces the
elimination argument with direct spectral measurements on the shipped solver
(new harnesses: `workflow/scripts/loop_spectrum.py`,
`interface_mode_trajectory.py`, `mode_rate_vs_drift.py`; curated:
`docs/method-comparison/.../tables/mode_rate_vs_drift.csv`,
`mode_rate_dt_series.csv`, figure `mode_rate_vs_drift.png`).

### 18.1 The delivered error is ONE function, and its shape explains the statics

leiaTestMeanCurvature now dumps every delivery per active face
(`leiaTestFaceCurvatureField.csv`). Production at N=128 on the exact-SDF
circle: the signed error collapses onto a single-valued function of the angle
between the CSF face normal and the interface normal (bin sd 0.10-0.43 against
an overall 1.53) -- +3.33 at alignment falling to -0.99 at 45-60 deg. The mesh
offers two face orientations while the interface normal rotates, so the
misalignment sweeps 0-90 four times per revolution: that IS the cos(4 theta)
pole bias (m=4 carries 99.7% of the smooth fluctuation, 20 cells/wave at
N=128 -- resolved, and mesh-locked, so refinement only shrinks its amplitude at
h^2). The narrow peak at alignment sits on the faces with the LARGEST
|snGrad(alpha)| (force-weighted mean error +0.97 vs plain +0.68), and the
peak/trough asymmetry is the whole skew of the distribution. Resolvability
lesson applied retroactively: interface-mode fits at N=64 are valid for
m <= 4 only (40 cells around the circumference).

### 18.2 The level-set gauge is exactly neutral -- measured, not asserted

Power iteration on the raw one-step delta-psi map converges to lambda = 1 with
|Mv - v| -> 0: perturbations that leave the zero set in place (a subspace
containing delta-psi ~ psi) change alpha, force, velocity NOT AT ALL, and are
neither damped nor amplified. This is the stability inequality's "zero damping
for delta-psi" as an operator property -- and it is the reservoir the
|grad psi| drift accumulates in. Corollary for measurement: the instability is
invisible to power iteration in raw delta-psi; the physical perturbation space
is INTERFACE DISPLACEMENT.

### 18.3 r_m(drift): the mode spectrum along the run's own death trajectory

Harness: perturb psi by eps*h*cos(m theta), restart the shipped solver from
snapshot base states with U and p_rgh RESET (isolating the psi-profile's
effect), fit a_m(t) = A0 exp(r_m t) cos(omega t + phi); validity = fit
residual AND omega against the capillary dispersion relation (on-plateau
m=2 fits: residual 0.002-0.005, omega within 1.5%). Base states: an N=128
serial production run to its own FPE at t = 0.1158 (the deck's serial
t_blow ~ 0.118 reproduced), snapshots labelled by their own band min|grad psi|.

    minGrad   t_snap    r_2      r_4      r_6    [1/s]
    0.9942    0.002    +18.8    +18.1    -45.3
    0.9773    0.009    +12.2     +7.0    -47.8
    0.9589    0.021    +10.5    +13.9    -51.2
    0.9568    0.073     +9.4     +2.9    -51.6   (same level, 50 ms older)
    0.9503    0.090    +19.8    +10.1    -58.2
    0.9371    0.096     +7.6    +15.3    -78.8
    0.8771    0.101    +29.2    -27.2    -84.9
    0.7076    0.107    +53.4   (+72.4)*  -29.4   (* fit residual 0.166)

Findings. (1) NO zero crossing in drift: m=2 is unstable ALREADY on the
near-clean profile at N=128 (+19 1/s) where the same mode at N=64 is damped
(-11 1/s) -- the N=64->128 switch of sec. 16.1, located: it is the m=2
capillary mode crossing neutral. (2) Drift AMPLIFIES (x3 by minGrad 0.71) but
does not trigger; m=6 stays damped everywhere, so the unstable content is
low-m SHAPE, not grid-scale corrugation. (3) Level-not-age holds for m=2
(+10.5 vs +9.4 at the same drift 50 ms apart). (4) The coupled trajectory's
quasi-stable plateau (minGrad 0.95-0.96 for 50 ms, serial) matches the weak
phase's e-folding time 1/r_2; the acceleration knee sits at the same LEVEL
(~0.94-0.95) in serial and np8 runs that leave it 2.3x apart in time.

### 18.4 The dt series: r_0 -> 0. The instability IS the explicit coupling

r_2 at fixed horizon, deltaT halved three times (clean state, minGrad 0.994):

    dt        7.5e-6   3.75e-6   1.875e-6   9.375e-7
    r_2 [1/s]  18.8      12.7      8.01       4.02      (omega identical, 1.001-1.006 of theory)

The finest halving halves the rate exactly; Richardson r_0 = 2*4.02 - 8.01 =
+0.03 1/s -- ZERO within noise against an operating rate of 18.8, with
c = r/dt consistent to 0.5% between the two finest points (4.29e6 vs 4.27e6).
On the drifted state (minGrad 0.877): 29.2 / 22.65 / 12.39 at dt/1/2/4,
two-point Richardson r_0 = +2.1 (an upper bound by the clean-state precedent),
c = 6.6e6: DRIFT RAISES THE COUPLING GAIN ~55%, NOT THE INTERCEPT. The
softened omega (0.92-0.93 of theory at every dt) is a physical property of the
drifted profile -- the delivered restoring force weakens as the foliation
degrades.

THE STATEMENT: the semi-discrete psi-loop is NEUTRAL (18.2 for the gauge
directions, r_0 ~ 0 for the interface modes); the parasitic-current
instability at N=128 is the EXPLICIT TIME COUPLING of the capillary force to
the m=2 interface mode, rate c*dt, with c raised by profile drift. This
supersedes sec. 16.1's "r0 does not vanish and rises toward fine grids" --
that intercept was fitted from max|U| over windows in which the drift (hence
c) was itself growing; at the mode level, on a held profile, the intercept is
zero. It also explains, at last, why L-stable Euler was load-bearing
(negative-results deck) and why N=64 is stable (the mode's physical/numerical
damping exceeds c*dt there).

CONSEQUENCES FOR THE PROGRAMME. (a) The semi-implicit / deferred-correction
capillary force -- shelved as task 24, argued against while r_0 looked finite
-- is now the PRINCIPAL licensed lever: the instability lives exactly in the
channel an implicit treatment closes. Dictionary-driven if revived, per the
standing rule. (b) SDF maintenance (WP6 frozen-band) is the c-reducer, not the
cure: budget = the measured c(drift) curve. (c) Delivery work is floor-only
(18.1): it sets the standing-current amplitude, not the exponent.
CAVEATS: one geometry (stationary droplet), N=128 primary; drifted-state r_0
from a two-point Richardson; the sec. 17.4 note "NO stability penalty" for
footEvalFace predates the coupled falsification (see STATUS sec. 4) and is
superseded.
