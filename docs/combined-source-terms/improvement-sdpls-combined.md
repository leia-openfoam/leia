# Combined source terms: `sdplsCombined` model and note corrections

Implementation brief for Claude Code. Repo root: `~/OpenFOAM/repos/leia` (this repo).
Read `docs/IMPROVEMENTS.md` first for global rules, notation, and file ownership.

## Task statement

The combined-source note (`docs/combined-source-terms/levelset_combined_source_note.tex`)
proposes the source `F = â + λ(1 − |∇ψ|)` — normal-strain cancellation (the SDPLS
`R` term) plus a logistic relaxation of the interfacial gradient toward 1 — with
the closed-form interfacial law `D_t^Σ g = λg(1−g)` and an invariant tube
`|g−1| ≤ η/λ` under bounded strain-estimate error `|δ| ≤ η < λ`. It is analysis
only; nothing is implemented. Two improvements:

1. **Code** — implement `sdplsCombined` as one new runtime-selectable source model
   in the existing `sdplsSource` hierarchy, plus the note's stage-0 unit gate.
   The calculus note supplies the geometric payoff to verify: interfacial
   |∇ψ| → 1 removes the tangential footpoint drift (`docs/IMPROVEMENTS.md`
   notation `T`), which is exactly what the `interfaceDefects` diagnostic from
   the SDPLS brief measures.
2. **Docs** — two corrections to the combined-source note where it borrows the
   footpoint projection without error control, plus symbol-collision hygiene.

## Files owned by this brief

Create:
- `src/leiaLevelSet/sdplsSource/sdplsCombined.{H,C}`
- one unit case (e.g. `cases/sdplsCombinedUnit/` with `Allrun.sh`), and/or a
  study config following the existing `sdpls*` YAML patterns.

Edit:
- `src/leiaLevelSet/Make/files` (register the new source file)
- `docs/combined-source-terms/levelset_combined_source_note.tex`

Read-only: `src/leiaLevelSet/sdplsSource/{sdplsSource,sdplsR,sdplsBeta}.{H,C}`,
`sdplsSource/discretization/*`, `sdplsSource/mollifier/*`, `sdplsSource/gradPsi/*`,
`applications/test/leiaTestSdplsSource/` — patterns to follow, not to modify.

## Steps — code

1. **Read** `sdplsR.{H,C}` and `sdplsBeta.{H,C}` end to end. `sdplsCombined` is
   their sum with one parameter; it must reuse, not duplicate:
   - the strain part via the same computation `sdplsR` uses
     (`(∇U & n̂) & n̂` with the dedicated `gradPsiSdpls`/`gradUSdpls` scheme
     keywords and the `gradPsi` model hierarchy),
   - the relaxation part following `sdplsBeta`'s dimensional handling
     (`sdplsBeta` returns `1[1/s]·(β − |∇ψ|)`; here the deficit is `(1 − |∇ψ|)`
     and the rate is a user parameter).
2. **Implement** `sdplsCombined : public sdplsSource` with
   `nonLinearPart(psi, U) = a(psi, U) + lambda_·(1 − mag(gradPsi(psi)))`, where
   `lambda_` is a required dictionary entry (`lambda`, dimensions 1/s — no silent
   default; the note's guidance is to choose it from the desired correction
   timescale τ_r = 1/λ, and a wrong implicit default would hide that decision).
   Support the existing `mollifier` and `discretization` selections unchanged —
   the Sc/Sp splitting (`explicit`, `simpleLinearImplicit`,
   `strictNegativeSpLinearImplicit`, `exponential`) already acts on whatever
   `nonLinearPart` returns. Register with `addToRunTimeSelectionTable`; add the
   source file to `src/leiaLevelSet/Make/files`. Header comment: the PDE, the
   interfacial logistic law, the invariant-tube property, one citation each to
   the combined-source note and (for the geometric payoff) the calculus note.
   KISS: this class should be of the same size as `sdplsBeta` (a few dozen
   lines).
3. **Stage-0 unit gate** (from the note's own verification program): v = 0,
   ψ₀ = g₀·x with g₀ ∈ {0.5, 2.0}, source-only evolution
   `∂_t ψ = λψ(1 − |∇ψ|)`. Exact solution: ψ stays linear with slope
   `g(t) = 1/(1 + (g₀⁻¹ − 1)e^{−λt})`. Run on a small 1D/2D mesh with
   `leiaLevelSetFoam` (velocity model off / zero), assert the band slope against
   the closed form at several times (script does the assertion, exit non-zero on
   failure). Also assert the sign-preservation property of the chosen
   discretization (no zero crossings created, Δt·max|r| within the resolution
   criterion the note states).
4. **Stage-1/2 study arm** (measurement, not pass/fail): clone one 2D-vortex
   SDPLS study config to add a `combined` arm (R vs beta vs combined, same
   discretization — the repo's one-discretization-per-study rule is checked by
   `make check-discretization`). Report band `||∇ψ|−1|`, shape, and volume orders
   alongside the existing arms. If the SDPLS brief's `interfaceDefects` function
   object is already merged, add it to the arm's `functions` (config-only usage —
   no source-file overlap) and report `T·h`: the combined source should hold the
   band `T` at or below the R arm's level and additionally bound pre-existing
   `|g−1|` deviation by the invariant tube. Expected per the theory: the ~O(h^1.2)
   band-gradient ceiling from transport truncation persists — the combined term
   adds robustness to strain-estimate error δ (equilibrium `g* = 1 + δ/λ`), not
   spatial order. State this expectation in the PR so the result is read
   correctly.

## Steps — docs (`levelset_combined_source_note.tex`)

1. **The Π₁ error statement.** The note's velocity-extension equation
   (`eq:velocity-extension-nonlinear`) uses
   `x − φ∇φ/|∇φ|²` as an exact stand-in for the metric projection p_t. Add a
   remark immediately after it: for a general regular φ this operator equals
   `p_t(x) + d²[½(nᵀHn/|∇φ|)n − ∇_Γ log|∇φ|] + O(d³)` — second-order accurate,
   with a normal defect part and a tangential drift part; the tangential part
   vanishes when |∇φ| is constant along the interface, which is precisely the
   regime the combined source drives the field toward (interfacial g → 1). Cite
   the calculus note. This turns an unstated idealization into a quantified
   modelling error and links the two documents.
2. **Symbol hygiene.** The note's `β` (constant of the simpler source F = β − g)
   and `q` (∇φ) collide with the calculus note's `β = nᵀHn` and `q = ∇_Γ|∇φ|`,
   and its `a` (scale factor φ/d) with leia's `a` (normal strain). Add a short
   notation-warning box (or rename the simpler-source constant to `β_s`) so the
   documents can be bundled without ambiguity. Keep edits minimal — this is a
   hygiene pass, not a rewrite.
3. **Implementation cross-reference.** Point the note's "discrete design"
   section at the new `sdplsCombined` class and the stage-0 gate, so analysis and
   code stay linked.

## Verification gates

1. Build (`./Allwmake`); `leiaTestSdplsSource` still passes (the discretization
   identity test is model-agnostic — run it with `sdplsCombined` selected as well
   if its dict permits, since all discretizations must agree at ψ = ψⁿ).
2. Stage-0 unit gate passes for both g₀ values and at least two λ values
   (λ·T_end ≈ 3 and ≈ 10), serial and np4.
3. Bit-identity of untouched models: an existing R-arm case reruns bit-identical
   (the new class must not perturb `sdplsR`/`sdplsBeta` code paths — it doesn't
   touch them, so this is a build-level regression check).
4. Study arm completes on the ladder; results reported, no promotion claims.

## Do not

- Modify `sdplsR`, `sdplsBeta`, the `Rdiv` variants (withdrawn with
  measurements — do not resurrect), any `discretization/`, `mollifier/`, or
  `gradPsi/` class.
- Touch `functionObjects/*` (SDPLS brief), `velocityExtension/*`,
  `phaseIndicator/*`, `semiLagrangian/*`, `levelSetImplicitSurfaces/*`, or
  solver headers (see the `docs/IMPROVEMENTS.md` matrix).
- Give `lambda` a default value.
- Use the combined source on the coupled stationary droplet in this PR — the
  coupled campaign has its own scoring rules (`docs/IMPROVEMENTS.md` rule 7) and
  belongs to a separate, gated effort.
