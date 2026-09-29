---
title: "The 2D advection orders of METHOD 8.3.7 were 3/2 too high (2026-09-26)"
description: "CORRECTED 2026-09-26 - every 2D order of the converged advection ladders (translation none 5.70, 3.15, -0.82; vortex none 4.15, 4.85, 3.14) came from h_eff = nCells^(-1/3) applied to a 2D case, so the spacing was N^(-2/3) instead of 1/N; the corrected orders are 3.80, 2.10, -0.54 and 2.76, 3.23, 2.09; no conclusion changes; the curated vortex CSV is not yet regenerated; the translation numbers are void since 2026-09-29 (a reversed flow), and one way the translation orders are 2.90, 2.16, 1.86"
aliases: [3/2 orders retraction, advection_convergence_table h_eff, advConv2D orders]
kind: retraction
status: retracted
part: advection
tags: [retraction, part/advection]
date: 2026-09-28
code: [workflow/scripts/advection_convergence_table.py, config/advConv2Dtranslation.yaml, config/advConv2Dvortex.yaml]
sources: [STATUS 11.2, METHOD 8.3.7, commit f8aaa8e, advConv2D*_convergence.csv, STATUS 11.19]
---
# The 2D advection orders of METHOD 8.3.7 were 3/2 too high (2026-09-26)

> CORRECTED 2026-09-26 (commit [f8aaa8e](https://github.com/leia-openfoam/leia/commit/f8aaa8e)). The claim, published on 2026-09-10 (commit 32a0d9b), was the set of 2D orders of the converged advection ladders: for uniform translation with `none` 5.70, 3.15 and -0.82 between the rungs N = 32 / 64 / 128 / 256, and for the reversed vortex with `none` 4.15, 4.85 and 3.14 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L627-L634), [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3242)). The measurement is a reading of the script: `workflow/scripts/advection_convergence_table.py` computed h_eff = nCells^(-1/3) for every case. A 2D case has N^2 cells, so its spacing was N^(-2/3) instead of 1/N, and every 2D order was 3/2 of the true value ([the script's own note](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/advection_convergence_table.py#L14-L20)). Against h = 1/N the translation `none` arm goes 3.80, 2.10, -0.54 and the vortex `none` arm 2.76, 3.23, 2.09; the `lipschitzCone` arms go 1.70, 1.87, 0.65 and 1.63, 1.21, 0.79 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L641-L662)). Scope: the orders only. The errors, the ratios between the arms (0.85, 3.6, 4.3, 1.9 on translation; 8.6, 18.9, 76.7, 189.7 on the vortex) and every conclusion of 8.3.7 are unchanged, because a ratio of two errors at the same N does not depend on the spacing ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L632-L635)). The 3D rungs use the cube root correctly. No data is void. CORRECTED 2026-09-29: the translation ladder ran a reversed flow, so its errors, ratios and orders are void for that reason, while the 3/2 correction of the script stands ([[retractions/reversed-2dtranslation]]). One way, the translation `none` arm goes 2.90, 2.16, 1.86 and the `lipschitzCone` arm 0.38, 0.08, 0.36 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697)).

## The claim, and where it lived

- [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L627-L680), `METHOD.md`: rewritten in place; its first paragraph records "RETRACTED AND CORRECTED 2026-09-26" and both sets of orders, and the tables now carry the corrected ones.
- [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3242), `STATUS.md`: the retraction, with the OPEN item to regenerate the curated CSVs.
- The curated tables `advConv2Dvortex_convergence.csv` and `advConv2Dtranslation_convergence.csv` ([the vortex table](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/advConv2Dvortex_convergence.csv), [the translation table](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/advConv2Dtranslation_convergence.csv)), last committed in 32a0d9b: their `hEff` column reads 0.0992 at N = 32 (the wrong spacing) and their `order` column holds the old values (4.1467 and 4.8518 for the vortex `none` arm, 5.6993 and 3.1512 for the translation `none` arm). On 2026-09-29 the reversed translation table moved to `VOID_reversedTranslation_20260929/`, and the one-way [translation table](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/docs/method-comparison/method-comparison-article/data/tables/advConv2Dtranslation_convergence.csv) at the old path carries h = 1/N (0.03125 at N = 32).
- The [distance-cone retraction](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L657-L667) quotes the vortex ladder: its ratios (8.6x to 189.7x) stand, its orders were the wrong ones until 2026-09-26, see [[retractions/distance-cone-bound-as-transport-bound]].

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The script's spacing is the mesh spacing. | h_eff = (V_domain / nCells)^(1/d) with d = 3 for every case; for a 2D case with N^2 cells that is N^(-2/3); h_eff read 9.92e-02 at N = 32 where h = 3.125e-02. | MEASURED, [script](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/advection_convergence_table.py#L14-L20), [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L629-L634) |
| The factor on the order. | An order p fitted against N^(-2/3) is 3/2 of the order against 1/N, because log(e1/e2) / log(h1/h2) scales with the exponent of N in h. | DERIVED, [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3239) |
| The corrected orders. | Translation `none` 3.80, 2.10, -0.54 (was 5.70, 3.15, -0.82); vortex `none` 2.76, 3.23, 2.09 (was 4.15, 4.85, 3.14); cone 1.70, 1.87, 0.65 and 1.63, 1.21, 0.79. The translation values are VOID since 2026-09-29: the flow was reversed. | MEASURED, [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L641-L662) |
| The conclusions change. | They do not: the unbounded translation still saturates at N = 256 (the error rises from 7.947e-05 to 1.159e-04) and the cone bound is still worse on the vortex by a factor that grows under refinement (8.6x to 189.7x). RETRACTED 2026-09-29 for the translation: the saturation came from a reversed flow; one way, the error falls to 4.442e-04 at N = 256 at order 1.86. | MEASURED, [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L634-L635), [L650-L669](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L650-L669), [METHOD 8.3.7 one way](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697) |
| The script now reads the dimension. | `dims` from `case_params.json`, default 3; h_eff = nCells^(-1/dims). | MEASURED, [script](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/advection_convergence_table.py#L106-L113) |

## What survives

1. Every conclusion of [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L627-L680) and of the [distance-cone retraction](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L490-L496): the ratios between arms on the vortex. RETRACTED 2026-09-29: the translation ratios and the saturation of the unbounded translation at N = 256 (order -0.54) came from a reversed flow ([[retractions/reversed-2dtranslation]]).
2. The 3D rungs and the [advection regression of 8.3.8](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L682-L692) (bit-identity at tolerance 0 does not read an order).
3. The vortex orders 2.76 to 3.23, consistent with the 2D vortex orders of the SL article and of the gradU re-run (2.84 at CFL 0.5, [[retractions/gradu-coupled-patch-contamination]]). The translation order 2.10 is void since 2026-09-29; one way, the orders are 2.90, 2.16, 1.86 ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/METHOD.md#L692-L697)).
4. The rule that an order is reported next to its error and its spacing ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L510-L537)), see [[concepts/richardson-ladders-and-orders]].

## Propagation (checklist, same commit)

Done:

- [x] `METHOD.md`: [8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L627-L636) rewritten with the retraction first and the corrected tables.
- [x] `STATUS.md`: [11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3235-L3242).
- [x] The script: the [RETRACTED note](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/advection_convergence_table.py#L14-L20) and the dimension-aware spacing ([L106-L113](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/advection_convergence_table.py#L106-L113)).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The curated `advConv2D*_convergence.csv` still hold the old `hEff` and `order` columns (last commit 32a0d9b, 2026-09-10); [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3241-L3242) lists the regeneration on Lichtenberg as OPEN (the studies live only there). [METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L635-L636) says the tables "are regenerated"; the repository does not hold them yet. CORRECTED 2026-09-29: the translation table is now the one-way table with h = 1/N (commit a8b444d8); the [vortex table](https://github.com/leia-openfoam/leia/blob/aaa0a7dd/docs/method-comparison/method-comparison-article/data/tables/advConv2Dvortex_convergence.csv) still holds the old columns (0.0992 at N = 32).
- [ ] The 3D hexahedral and polyhedral converged rungs were "RUNNING, not reported" ([METHOD 8.3.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L677-L678)); no landing is recorded.
- [ ] Whether any other 2D table was built with the same script before the fix: not recorded.

## Related

Hubs: [[hubs/advection]], [[hubs/verification]]. Siblings: [[concepts/richardson-ladders-and-orders]], [[concepts/advection-regression-set]], [[concepts/bit-identity-and-inertness-gates]], [[models/sl-value-bound]], [[concepts/value-bounds-and-clips]], [[retractions/distance-cone-bound-as-transport-bound]], [[retractions/gradu-coupled-patch-contamination]], [[retractions/reversed-2dtranslation]].

## Log

### 2026-09-28
Written from STATUS 11.2, METHOD 8.3.7 and the script. Corrected 2026-09-26. Entered in [[retraction-log#2026-09]].

### 2026-09-29
The translation statements are CORRECTED: the translation ladder ran a reversed flow and is void (STATUS 11.19, [[retractions/reversed-2dtranslation]]); the 3/2 correction of the script stands. The one-way orders (2.90, 2.16, 1.86; cone 0.38, 0.08, 0.36) are in the verdict; the saturation is RETRACTED; the curated translation table is regenerated one way, the vortex table is not.
