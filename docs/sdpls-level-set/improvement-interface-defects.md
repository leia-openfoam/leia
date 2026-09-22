# SDPLS: interface-defect diagnostics (T, N) and the missing geometric justification

Implementation brief for Claude Code. Repo root: `~/OpenFOAM/repos/leia` (this repo).
Read `docs/IMPROVEMENTS.md` first for global rules, notation, and file ownership.

## Task statement

Two improvements:

1. **Code** — a new function object `interfaceDefects` that measures, per time
   step over the interface band, the two first-order defect fields of a general
   level-set function:
   - tangential drift magnitude `T = |(I − n̂n̂ᵀ)∇log|∇ψ||`
   - normal defect `N = n̂·∇log|∇ψ|`

   This closes the SDPLS article's explicitly open item (its §source-error
   analysis, lines ~544–551: the band metric `||∇ψ|−1|` senses the *normal*
   component of the gradient perturbation while the source error δa senses only
   the *tangential* one — "the estimate as written does not separate them") and
   supplies the tangential-structure diagnostic the curvature program says does
   not exist yet.
2. **Docs** — a short "geometric consequences" subsection in the SDPLS article
   giving the method its missing justification theorem.

## Files owned by this brief

Create: `src/functionObjects/interfaceDefects/interfaceDefects.{H,C}`
(model the class layout, IO, and CSV conventions on the existing
`gradPsiError` / `gradPsiErrorCSV` / `psiConservationCSV` function objects in the
same library).

Edit: `src/functionObjects/Make/files` (add the new source).

Edit (docs): `docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex`
(and `refs.bib` if a citation is added).

Read-only: everything else. In particular do **not** touch
`src/leiaLevelSet/sdplsSource/*` (owned by the combined-source brief) or
`src/leiaLevelSet/narrowBand/*`.

## Background

Why these two scalars matter for SDPLS specifically. The SDPLS source
`f = n̂ᵀ∇v n̂` keeps `|∇ψ| = 1` **on the interface only** (the article's own
controlling identity `Dg/Dt = g(f − a)` holds at ψ = 0). The calculus note proves
what that buys geometrically:

> If |∇ψ| is constant on Γ then the tangential drift T vanishes on Γ, the
> normalized gradient satisfies `n_ψ(x) = n(Π(x)) + O(d²)` throughout the band
> (not just on Γ), and the cheap footpoint Π₁'s leading error becomes purely
> normal.

So `T` is exactly the quantity the R source should suppress and `N` the one it
deliberately leaves free. Today neither is measured — the article reports only
`||∇ψ|−1|` norms, which mix the two. Both fields need **one** extra gradient of
one scalar: with `g = |∇ψ|`, `w = ∇(log g)`,

```
N = n̂ · w,          T = |w − (n̂·w) n̂|,        n̂ = ∇ψ/(g + SMALL)
```

— no Hessian, no fit. Noise is acceptable: this is a diagnostic, not a solver
operator (the repo's own rule for diagnostics).

## Steps — function object

1. **Read** `src/functionObjects/gradPsiErrorCSV.*` and one more neighbor for the
   established pattern: construction from dict, `execute()/write()` split,
   CSV file handling on the master rank, field registration.
2. **Implement `interfaceDefects`:**
   - Inputs (dict, with defaults): `psiName` (default `psi`), `gradScheme`
     keyword to use for both gradients (default the mesh's `grad(psi)`),
     `bandWidthCells k` (default 3): band = cells with `|ψ|/max(g,SMALL) ≤ k·h_c`
     (iso-agnostic — do not use raw |ψ|, which fails on non-SDF fields; this is
     the same lesson as the curvature band gate).
   - Per `execute()`: compute `g`, `w = ∇log(max(g, SMALL))`, then per band cell
     `N` and `T` as above. Optionally (dict switch `writeFields`, default off)
     register dimensionless per-cell fields `TdefectH = T·h_c`, `NdefectH = N·h_c`
     for visualization.
   - Per `write()`: one CSV row per time — time, band cell count, L2 and L∞ of
     `T·h` and `N·h` over the band, plus L2 of `|g−1|` over the same band so the
     three metrics are directly comparable on identical support.
   - Parallel: the norms are global reductions (`gSum`/`gMax`); band membership is
     cell-local; no seam handling needed. Master-only CSV writes, as in the
     existing FOs.
   - Clean Code: no state beyond the output stream; small pure helper for the
     band predicate; header documents the two formulas and their meaning in two
     sentences each.
3. **Register** the source file in `src/functionObjects/Make/files`; build.

## Verification gates — function object

1. **Analytic unit checks** (a small utility case or an `Allrun` script; both
   fields are trivial to initialize with `setExprFields` or a short `codeStream`):
   - Benchmark A field `ψ = (1 + αx)·y` (flat interface, tangential scale
     variation): expect `N = 0`, `T = α/(1+αx)` at y = 0, converging at the
     gradient scheme's order.
   - Benchmark B field `ψ = y + (c/2)y²`: expect `T = 0`, `N = c` at y = 0.
   Assert both to tolerance in the script (exit non-zero on failure).
2. **SDPLS study replay.** Add the function object to the existing 2D-vortex
   SDPLS study configs (`controlDict.functions` include — config files, not
   source) for three arms: no source, R, beta. Expected signature (this is the
   article's missing measurement, so record whatever is found): the R arm should
   hold band `T·h` markedly below the no-source arm while `N·h` behaves like the
   existing `||∇ψ|−1|` metric; beta should sit in between. No pass/fail here —
   it is a measurement, reported in the PR.
3. **Overhead**: one gradient + reductions per step; confirm wall-clock change on
   the N=128 vortex arm is below ~2%.

## Steps — article

In `sdplsLevelSet.tex`, add a subsection at the end of the method analysis
(near the `Dg/Dt = g(f−a)` identity):

1. State the hierarchy theorem (calculus note, translated to the article's
   notation ψ, g = |∇ψ|): interfacial `g ≡ 1` ⇒ tangential drift
   `∇_Γ log g = 0` on Γ ⇒ `n_ψ(x) = n(Π(x)) + O(d²)` in the band and the cheap
   projection `x − ψ∇ψ/|∇ψ|²` has a purely **normal** O(d²) error. One short
   proof sketch or a citation to the note; do not reproduce the full derivation.
2. Connect it to the article's open split: identify the tangential component of
   the gradient perturbation with `T` and the normal with `N`; note that the
   source can control the interfacial evolution of `g` (hence T on Γ) but not N,
   and that this is the precise sense in which "SDPLS is exactly local at the
   interface".
3. Reference the new `interfaceDefects` measurements from gate 2 (one small table
   or one sentence per arm — keep it lean; the article is already long).
4. Symbol hygiene: the article's `a` (normal strain) vs the note's `a` (scale
   factor) — keep the article's symbols, spell the imported quantities as
   `∇_Γ log|∇ψ|` and `nᵀHn/|∇ψ|` rather than importing `T_φ/N_φ` macros.

## Do not

- Touch `src/leiaLevelSet/sdplsSource/*` (combined-source brief owns it),
  `velocityExtension/*`, `phaseIndicator/*`, `semiLagrangian/*`,
  `levelSetImplicitSurfaces/*`, or any solver header
  (see the `docs/IMPROVEMENTS.md` matrix).
- Put the diagnostic inside the solver loop or a library class — it is a
  function object only.
- Use raw `|ψ| < k·h` band tests (breaks on non-SDF fields).
- Add smoothing/filtering to the diagnostic fields.
