---
title: "The linear semi-Lagrangian pre-print"
description: "The map of the linear semi-Lagrangian pre-print (679 lines, last edited 2026-08-06) and its deck; the nestedLSQ orders, the linearTaylor instability, and the naming and provenance questions that stand against it."
aliases: []
kind: study
status: settled
part: advection
tags: [study, part/advection]
date: 2026-09-28
code: [docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex, docs/linear-semi-lagrangian-level-set/lsl-level-set-presentation/linear-semi-lagrangian-level-set.template.html, workflow/Snakefile.sl-linear, src/leiaLevelSet/semiLagrangian]
sources: [LSL article, LSL deck, RM lines 18-28, gradU post-mortem section 3, kb-raw B5]
---
# The linear semi-Lagrangian pre-print

> **Verdict (2026-09-28).** The pre-print `linearSemiLagrangianLevelSet.tex` documents the linear reconstruction line as a written, measured companion of the quadratic line ([title, line 43](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L43); last commit d2b50f3, 2026-08-06). Its measured result stands: the least-squares linear fit converges at first to second order at CFL 1/2, the raw Taylor extrapolation diverges in pure advection, and the stability is conditional and flow-dependent ([`sec:conv2d`, line 421](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L421); [`sec:conv3d`, line 469](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L469)). Three things stand against it: all six of its studies are on the gradU contamination list and the article predates the fix ([[retractions/gradu-coupled-patch-contamination]]); the name `nestedLSQ` is not a code name, and the roadmap describes the selected model `quadraticTaylor` as a quadratic built from twice-differentiated psi ([RM lines 18-28](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L18-L28)); and the "two-phase workhorse" role of `linearTaylor` is superseded by the production decision for the uncached quadratic fit ([[decisions/sl-reconstruction-uncached-qwls]]).

## What it is

**The article.** `docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex`, 679 lines, `elsarticle`, one author. Pre-print PDF: https://leia-openfoam.github.io/leia/preprints/linearSemiLagrangianLevelSet.pdf. Sections: Introduction ([`sec:intro`, 94](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L94)); Method with the shared foot ([`sec:foot`, 176](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L176)), the linear least-squares reconstruction ([`sec:recon`, 195](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L195)), the `linearTaylor` instability ([`sec:lintaylor`, 241](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L241)), the safeguards ([`sec:safeguards`, 300](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L300)) and the indicator ([`sec:indicator`, 314](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L314)); Software ([`sec:software`, 328](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L328)); Verification with metrics ([`sec:metrics`, 356](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L356)), 2D ([`sec:conv2d`, 421](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L421)), 3D ([`sec:conv3d`, 469](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L469)), the fields ([`sec:fields`, 531](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L531)), the 3D interfaces ([`sec:iso3d`, 555](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L555)) and the comparison of the three reconstructions ([`sec:comparison`, 589](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L589)); Limitations ([623](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L623)); Conclusions ([647](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L647)).

**The deck.** https://leia-openfoam.github.io/leia/decks/linear-semi-lagrangian-level-set.html, 34 sections in 7 horizontal groups: `#/1` roadmap; `#/2` method (7 vertical slides, including "linearTaylor is unstable in pure advection" and "linearTaylor's redemption: the two-phase workhorse"); `#/3` 2D verification (4); `#/4` 3D verification (7); `#/5` three reconstructions, two design points, one cautionary tale (2); `#/6` conclusions and reproduction (2).

**The data.** `lsl-level-set-article/data/tables/` (10 entries: `linearConv2Dvortex_errors.csv`, `linearConv2DvortexClip_errors.csv`, `linearConv3Dshear_errors.csv`, `linearConv3DshearPoly_errors.csv`, `linearConv3Ddeformation_errors.csv`, `linearConv3DdeformationPoly_errors.csv`, `lsl_convergence.csv`, `lsl_convergence_orders.csv`, `convergence_orders.tex`, `convergence_orders_extended.tex`) and `data/figures/` (11 PNG files). The presentation folder mirrors the same 10 tables and 11 figures. No `data/archive/` folder.

**How to build.** `make article-lsl`; `make sl-linear` (`workflow/Snakefile.sl-linear` runs the five `config/linearConv*.yaml` studies and rebuilds deck and article, [Makefile lines 217-219 and 263-265](https://github.com/leia-openfoam/leia/blob/8867581/Makefile#L217-L219)); `make decks`.

## Why it matters

It is the only document that measures the two linear reconstructions against each other and records why the single-cell Taylor extrapolation must not be used for reinitialisation-free advection. It also fixes the shared machinery (foot, indicator, band) that the quadratic line uses ([[concepts/linear-semi-lagrangian]]).

## Where in the code

`src/leiaLevelSet/semiLagrangian/` (the `slReconstruction` family: `linearTaylor`, `linearWeightedLeastSquares`, `quadraticTaylor`, and the production `uncachedQuadraticWeightedLeastSquares`); `config/linearConv*.yaml`; `workflow/Snakefile.sl-linear`.

## Evidence

The orders below are least-squares slopes of the generated table [`data/tables/convergence_orders.tex`](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/data/tables/convergence_orders.tex) and [`lsl_convergence_orders.csv`](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/data/tables/lsl_convergence_orders.csv).

| claim | number | where |
|---|---|---|
| 2D reversed vortex, CFL 1/2, seven rungs | shape 1.070, volume 1.452, band gradient 2.082 (global gradient 0.516) | [`sec:conv2d`, 421](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L421), MEASURED |
| 2D, CFL 1 | stable envelope of 4 rungs (to N=90), shape 0.966; destabilised from N=128 | same section, MEASURED |
| `linearTaylor` on the same vortex | band gradient defect 1.4e9 at N=128 and 2.0e21 at N=256; shape order -0.9; `nestedLSQ` 1.5e-3 and 4.8e-4 at the same rungs | [`sec:lintaylor`, 241, Table `tab:lintaylor`](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L241), MEASURED |
| 3D deformation, CFL 1/2 | shape 2.282 hexahedral (5 rungs to N=160), 1.521 polyhedral; band gradient 1.547 hexahedral | [`sec:conv3d`, 469](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L469), MEASURED |
| 3D shear, hexahedral | stable envelope of 2 rungs (N=32, 50), shape 1.180; first destabilised rung N=80; gradient defect O(1e5) and shape error 3.7e-2 at N=128 | same section and the CSV column `hLimit=0.0125`, MEASURED; the abstract quotes the limit as "past N=80" and the section as "past N=50", which are the same envelope |
| 3D shear, polyhedral | shape 1.775, band gradient 2.098, converges over the full ladder | same section, MEASURED |
| Stencil clip | "required on polyhedral meshes" | [`sec:safeguards`, 300](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L300), superseded: the clip is falsified as a fix, [[decisions/sl-clip-and-value-bound-off]] |

## Retracted or superseded inside it

1. **Provenance.** The criterion "solver leiaSemiLagrangeLevelSetFoam and np > 1" puts all six `linearConv*` studies (np 4 and 8, all published) on the affected list of the coupled-patch post-mortem ([section 3, lines 102-146](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L102-L146)). The article was last edited on 2026-08-06, twenty days before the fix. Its fitted orders are provisional until the studies are re-run; a re-run is not recorded in the sources of this note ([[retractions/gradu-coupled-patch-contamination]]).
2. **The name of the method.** The article calls its method `nestedLSQ` and describes it as a linear least-squares fit ([`sec:recon`, 195](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L195)). Its reproduction entry selects `reconstruction quadraticTaylor` ([`sec:software`, 328](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L328)). The roadmap states that `quadraticTaylor` was formerly named `nestedLSQ`, builds a quadratic expansion from twice-differentiated psi, and needs the stencil clip to stay bounded; it adds that no genuinely linear reconstruction measured there is both stable and convergent ([RM lines 18-28](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L18-L28)). The code map lists no `nestedLSQ` member ([[models/sl-reconstruction]]). OPEN: which model the published orders belong to.
3. **The two-phase role of `linearTaylor`.** The article states that `linearTaylor` serves the two-phase solver and outlives the quadratic pipeline ([`sec:comparison`, 589](https://github.com/leia-openfoam/leia/blob/8867581/docs/linear-semi-lagrangian-level-set/lsl-level-set-article/linearSemiLagrangianLevelSet.tex#L589)). That was the "consistent piecewise-linear pipeline" of the SL negative-results deck (`#/6/1`, a partial positive). The production two-phase reconstruction is `uncachedQuadraticWeightedLeastSquares` ([[decisions/sl-reconstruction-uncached-qwls]]).
4. **The clip.** Item 7 of the evidence table: SL_CLIP is false and the value bounds are `none` ([[concepts/value-bounds-and-clips]]).

## What it does not cover

Two-phase coupling (measured elsewhere, [[studies/sl-quadratic-pre-print]]), polyhedral 2D, the value-bound family, and any re-run after 2026-08-26.

## Open questions

- Which reconstruction class produced the published `nestedLSQ` numbers (item 2 above).
- Whether the six studies were re-run on the fixed binaries (item 1 above).

## Related

[[hubs/advection]], [[hubs/method-lines]]. [[models/sl-reconstruction]], [[models/sl-value-bound]], [[concepts/linear-semi-lagrangian]], [[concepts/departure-foot-ab2-centring]], [[concepts/idec-defect-correction-failure]], [[cases/kinematic-advection-cases]]. Siblings: [[studies/sl-quadratic-pre-print]], [[studies/npsl-design]], [[studies/method-comparison]].

## Log

### 2026-09-28
Created from the article at 8867581, the deck template, the data folder and the post-mortem list.
