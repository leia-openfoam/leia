# Semi-Lagrangian level set: second-order interface offset in the phase-indicator chain

Implementation brief for Claude Code. Repo root: `~/OpenFOAM/repos/leia` (this repo).
Read `docs/IMPROVEMENTS.md` first for global rules, notation, and file ownership.

## Task statement

Replace the first-order interface offset `d = ψ/|∇ψ|` used by the phase-indicator
chain of `leiaSemiLagrangianLevelSetTwoPhaseFoam` with the Hessian-corrected root
`d₂`, behind a dictionary switch whose default is bit-identical to today. This is
the improvement leia's own `STATUS.md` §7 lists as *"Fix the signed-distance
assumption where it actually survives"* — the named next structural candidate,
untouched as of 2026-08-25.

## Files owned by this brief

Edit:
- `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H` (anchor: lines 154–161)
- `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createSLFields.H` (only if a registered field is needed; see Step 3)
- `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H` (only if the per-step refresh must be hooked there)
- `src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C` and `.H` (anchor: lines 130–140 in the `.C`)

Create:
- one study/config entry (YAML under `config/` or a case under `cases/`) for the
  verification gate below, following the existing `faceCurvatureEllipsoidPsi2D.yaml`
  pattern.

Read-only (do NOT edit): `src/leiaLevelSet/semiLagrangian/slReconstruction.H`
(the `signedOffset`/`fitDerivatives` API, lines ~293–332),
`uncachedQuadraticWeightedLeastSquaresReconstruction.{C,H}` (`offsetDistance`
inline in the `.H` ~lines 178–196; **byte-frozen**, call only),
`src/leiaLevelSet/phaseIndicator/levelSetPlaneReconstruction.C`.

## Background (why this is correct and cheap)

The CSF face force is `σ κ_f snGrad(α) |S_f|`. The κ factor was already made
signed-distance-free (parallel-surface inversion). The α factor still locates the
interface with a linear LSQ plane plus the **first-order** offset `ψ/|∇ψ|`. Its
leading error is

```
e_d = −½ (h_nn / G) d² + O(d³),      h_nn = n_ψᵀ H n_ψ,  G = |∇ψ|
```

(the calculus note's normal defect; leia's `plan-curvature-stabilization.md`
§11.3 derives the same coefficient independently). The fix is the corrected
distance

```
d₂ = 2ψ / (G + sqrt(G² − 2 ψ h_nn))        [= ψ/G + ψ² h_nn/(2G³) + O(d⁴), exact when ψ is quadratic along the ray]
```

which **already exists** as `offsetDistance` / `signedOffset(c, fallback)` on the
quadratic reconstruction, including its guard `D = G² − 2ψh_nn ≥ β·G²` (β = 0.25)
with a first-order fallback. Zero extra derivative order: `h_nn` is already
computed for curvature. Properties to preserve: exact zero-preservation
(ψ_c = 0 ⇒ d₂ = 0) and sign(d₂) = sign(ψ_c).

Note the scope: `d₂` corrects the **normal** distance to O(d³). It does not (and
cannot) correct the tangential footpoint drift `−d²∇_Γ log G` — irrelevant here,
because the offset is consumed as a scalar distance along the local normal.

## Steps

1. **Read** all owned and read-only files; confirm the anchors. Identify exactly
   which phase-indicator object the SL two-phase solver instantiates at runtime
   (check a shipped `fvSolution` from `cases/stationaryDroplet2D` and
   `createSLFields.H`) — the correction must land on the live code path(s):
   `faceAreaFraction.H` (face liquid-area fractions for the mass flux) and, if
   live, `geometricPhaseIndicator.C` (cell α by clipping).

2. **faceAreaFraction.H.** The solver has the reconstruction object (`recon`) in
   scope. Where the code currently forms the plane offset from `ψ/|∇ψ|`
   (lines ~154–161), add a switch — dictionary entry `interfaceOffset` in the
   solver's `levelSet` subdict, values `linear` (default, bit-identical) and
   `quadraticRoot`. Under `quadraticRoot`, obtain the offset from
   `recon.signedOffset(cell, fallback)`; when `fallback == true` or the fit order
   is < 2, use the existing linear expression unchanged. Keep the change local: a
   small helper lambda or inline function, one call site per offset use, no
   copy-pasted guard logic (Clean Code: one source of truth for the branch).

3. **geometricPhaseIndicator.C.** This class sits in `src/` and only has the
   linear plane fit — no Hessian. Follow the same dictionary pattern
   (`interfaceOffset linear|quadraticRoot`, default `linear`). For
   `quadraticRoot` it needs `h_nn` per band cell; the leia-conformant route is a
   registered `volScalarField` (e.g. `hnnPsi`) written by the solver each step
   from the fit (`fitDerivatives(c, g, H)` → `n·(H·n)`), looked up via
   `mesh.cfindObject<volScalarField>("hnnPsi")`; fall back to `linear` per cell
   when the field is absent or the guard trips. Register/refresh the field in
   `createSLFields.H` / `slAlphaEqn.H` (owned by this brief). If reading the code
   shows the SL solver never instantiates `geometricPhaseIndicator` on its live
   path, implement Step 2 only and record that finding in the PR description —
   do not add dead machinery (KISS).

4. **Shift, don't rebuild.** In both sites the correction is a shift of the
   existing plane constant along the plane normal by `(d₂ − d₁)`; the plane
   normal, the clipping, and the Detrixhe–Aslam fill formulas stay untouched.

## Verification gates (run in order; stop and report on failure)

1. **Bit-identity (default off).** Build (`./Allwmake`), run one existing smoke
   case (e.g. the `benchVortex` smoke or `stationaryDroplet2D` first steps) with
   `interfaceOffset linear`: α and ψ fields bit-identical to the pre-change build.
2. **Unit gate — Benchmark C of the note.** Circle R = 0.25, ψ = F(d) with
   F(r) = r + (c/2)r², c·h = O(1) at the coarsest level. Exact α is known from the
   circle. Ladder N = 32…256: with `linear`, α/volume error saturates at the
   −½(h_nn/G)d² = −½c·d² offset bias; with `quadraticRoot` the bias term must
   vanish to fit accuracy and the volume error order improve toward 2. Reuse the
   `pyFoamStudy` ladder machinery; add the config file owned by this brief.
3. **Varying-curvature arm.** Rerun the existing ellipse-as-quadratic-form
   indicator arm (`faceCurvatureEllipsoidPsi2D`-style init where β = |∇ψ| varies
   2× along Γ): the gap between the SDF arm and the quadratic-form arm in α/volume
   error must close with `quadraticRoot`.
4. **Transport re-gate.** 2D vortex ladder (CFL 0.5): shape and volume orders must
   not degrade (offset moves α, hence the force support and cut-cell set — this is
   leia's own stated gate condition for this change).
5. **Moving-interface rule.** Include one translating case; the switch must not be
   gated on anything static-only.
6. **Parallel.** serial ↔ np4 agreement at the repo's standing ~1e-12 level.

Only after 1–6 pass may a coupled `stationaryDroplet2D` comparison be run, scored
per `docs/IMPROVEMENTS.md` rule 7 (exponent, not prefactor). Do not claim
stability improvements from this change; claim removal of a quantified
first-order bias in `snGrad(α)`.

## Do not

- Edit `offsetDistance`/`footPointDistance` or anything in
  `uncachedQuadraticWeightedLeastSquaresReconstruction.*`.
- Touch `src/leiaLevelSet/velocityExtension/*`, `src/functionObjects/*`,
  `src/levelSetImplicitSurfaces/*`, `src/leiaLevelSet/sdplsSource/*`
  (owned by the other briefs — see `docs/IMPROVEMENTS.md`).
- Change the Detrixhe–Aslam tetrahedron fill or the plane LSQ itself.
- Flip the default to `quadraticRoot` in this PR.
