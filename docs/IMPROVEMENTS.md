# Level-Set Calculus → leia: implementation briefs

Prepared 2026-08-25. Each improvement brief lives in the `docs/` method folder it
improves and is a **self-contained implementation brief for Claude Code**, mapping
results from the technical note *Calculus of Level-Set Geometry in a Tubular
Neighborhood* onto one leia level-set method. Source note:
`level_set_calculus_note.tex` in
`C:\Users\TomislavMaric\Documents\research\dieter\input\2026-LevelSet-Calculus`
(from WSL: `/mnt/c/Users/TomislavMaric/Documents/research/dieter/input/2026-LevelSet-Calculus/`).
Full background analysis: the artifact
[Tubular Calculus in leia](https://claude.ai/code/artifact/12a1d8af-a0f1-4513-bb55-91b45583b07c).

Run Claude Code **inside WSL in the repo root** (`~/OpenFOAM/repos/leia`) so
builds and tests work. Read the repo's `CLAUDE.md`, `METHOD.md`, and `STATUS.md`
before any edit.

## Briefs and file ownership

The briefs are **pairwise disjoint in the source files they modify**, so they can be
executed independently and in parallel (e.g. one git branch or worktree each).
A brief may *read* any file; it may only *edit or create* files it owns. Since the
library split of 2026-09-23 (`docs/plan-library-split-and-build-policy.md`, WP3) the
owned source lives in one library each -- named in the table -- and a new model is
added to that library's `Make/files`, never to another one's.

| Brief | Improvement | Owns (edits/creates) |
|---|---|---|
| [semi-lagrangian-level-set/improvement-interface-offset.md](semi-lagrangian-level-set/improvement-interface-offset.md) | Second-order interface offset (d₂) in the phase-indicator chain | `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceAreaFraction.H`, `.../createSLFields.H`, `.../slAlphaEqn.H`; `src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.{C,H}` (libleiaCore); new study config |
| [velocity-extension/improvement-metric-footpoint.md](velocity-extension/improvement-metric-footpoint.md) | Metric-consistent footpoint descent; phase-name unblock; orthogonality residual; article theory | `src/leiaLevelSet/velocityExtension/*` (libleiaVelocityExtension); `docs/velocity-extension/*` |
| [sdpls-level-set/improvement-interface-defects.md](sdpls-level-set/improvement-interface-defects.md) | Interface-defect diagnostics (T_φ, N_φ); article justification section | `src/functionObjects/*` (new function object + `Make/files`); `docs/sdpls-level-set/*` |
| [normal-projected-semi-lagrangian/improvement-drift-gate.md](normal-projected-semi-lagrangian/improvement-drift-gate.md) | Footpoint-drift benchmark gate; documentation corrections | `src/levelSetImplicitSurfaces/*` (new surfaces + `Make/files`); new `applications/test/leiaTestFootpointDrift/`; `docs/normal-projected-semi-lagrangian/*` |
| [combined-source-terms/improvement-sdpls-combined.md](combined-source-terms/improvement-sdpls-combined.md) | `sdplsCombined` source model (SUBSUMED 2026-09-26 by `plan-halo-limited-gradient-control.md`: law `linearQ` + strain weight `full` of the `gradientControl` source); note corrections (Π₁ error statement, symbol collisions) | new `src/leiaLevelSet/sdplsSource/sdplsCombined.{H,C}` and its line in `src/leiaLevelSet/sdplsSource/Make/files` (libleiaSdplsSource); `docs/combined-source-terms/*`; new unit case |

Not covered on purpose: **GRL** (leia's record retires plane-based band rewrites —
dead end #7 in `docs/plan-curvature-stabilization.md`; nothing in the note overturns
that measurement) and the **method-comparison** article (a benchmark, not a method).

## Global engineering rules (apply to every brief)

1. **Read before edit.** Every brief lists its anchor files with line numbers as of
   2026-08-25; verify them against the working tree first — the repo moves.
2. **Bit-identity by default.** New behavior goes behind a runtime-selectable
   dictionary entry whose default reproduces current behavior exactly. leia runs
   bit-identity regression gates; switching a default is a separate, gated decision.
3. **Frozen code.** `footPointDistance` in
   `src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C`
   is byte-frozen by an in-code contract (comment near lines 994–999). Call it;
   never edit it. The same holds for `offsetCorrection none` paths documented as
   bitwise-unchanged.
4. **leia design idioms.** Runtime-selectable class hierarchies
   (`addToRunTimeSelectionTable`), dictionary-driven configuration with documented
   defaults, cell-local geometry (no field access inside per-cell geometry
   routines), exact zero-set preservation for anything multiplying or rewriting ψ,
   and processor-seam consistency (face values computed from rank-local data must be
   `syncTools`-synchronized; serial ↔ np4 agreement ~1e-12 is the standing gate).
5. **KISS / Clean Code.** Smallest change that does the job; no new base classes
   when an existing hierarchy fits; intention-revealing names; comments only for
   constraints the code cannot express; no dead switches; one concern per class.
6. **Binding methodological rule** (repo, v0.3): no candidate may be specialized
   for a static test. Whatever you change must run identically on deforming,
   translating interfaces, and promotion requires the kinematic transport suite
   plus a moving coupled gate.
7. **Scoring coupled claims.** Static accuracy is measured in this repo to be
   anti-correlated with coupled stability (three falsified deliveries). Any change
   touching the coupled solver is scored on the noise gain `G h² ≤ 0.65`, fitted
   error order ≥ 1.9 on the **varying-curvature ellipse gate** (not the circle),
   and the blow-up exponent over ≥ 3 resolutions — never on static prefactors.
8. **Sign conventions.** The calculus note uses Prüss–Simonett `κ = −div n`
   (sphere: −(m−1)/R); leia uses `κ = +div(∇ψ/|∇ψ|)`. All formulas quoted in the
   briefs are κ-sign-free or already stated in leia's convention. The note's
   symbols also collide with leia's (`β`, `q`, `a` mean different things); the
   briefs spell every quantity out.

## Common notation used in the briefs

For a level-set field ψ with interface Γ = {ψ = 0}, signed distance d, metric
(closest-point) projection Π, and G = |∇ψ|, n_ψ = ∇ψ/G, H = Hessian(ψ),
h_nn = n_ψᵀ H n_ψ:

- Cheap footpoint:      `Π₁(x) = x − ψ ∇ψ / |∇ψ|²`, error `Π₁ − Π = d² [ ½(h_nn/G) n − ∇_Γ log G ] + O(d³)`
- Corrected distance:   `d₂ = ψ/G + ψ² h_nn / (2G³) = d + O(d³)`; closed root form `d₂ = 2ψ / (G + sqrt(G² − 2 ψ h_nn))`
- Tangential drift:     `T = ∇_Γ log|∇ψ|` (footpoint drifts by −d²·T; **curvature-independent**, invisible on any signed-distance field)
- Normal defect:        `N = h_nn / G` (distance error of ψ/G is −½·N·d²)
- Orthogonality residual for a candidate foot y of query x: `R_tan = |(I − n̂n̂ᵀ)(x − y)|`, n̂ at y. Zero for the metric foot; `|ψ(y)| small` does **not** imply `R_tan` small.
