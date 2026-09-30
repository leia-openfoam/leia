---
title: "\"The translating droplet ran\": the case was a closed box (2026-09-02)"
description: "VOIDED 2026-09-02 - every result of cases/translatingDroplet2D before commit 440107f: the mesh had no inlet and no outlet, the pressure projection annihilated the free stream on step 1, and the disturbance metric measured the dead stream at about 2 U0"
aliases: [closed box void, translatingDroplet2D closed box, VOID_closedBox_20260902]
kind: retraction
status: voided
part: mass-flux
tags: [retraction, part/mass-flux]
date: 2026-09-28
code: [cases/translatingDroplet2D/system/blockMeshDict, config/translatingFreeStreamGate2D.yaml, config/translatingRepaired2D.yaml, config/translatingRepairedEqualRho2D.yaml]
sources: [STATUS 0, STATUS 11.13, CLAUDE wrong-setup rule, METHOD 6 and 8.1 and 10, DP MASS_FLUX and MASS_FLUX_BOUND_RHO and RHO_DDT_SCHEME comments, SL article sec:translating]
---
# "The translating droplet ran": the case was a closed box (2026-09-02)

> VOIDED 2026-09-02. The claim was "`rhoLENT` plus a momentum time scheme matched to the density equation makes the translating droplet run the full horizon at N = 128", with a pairing table, a best-configuration table, a droplet-leaves-the-domain retraction and a `div(rhoPhi,U)` scheme comparison ([the old STATUS section 0, commit 3dece5a](https://github.com/leia-openfoam/leia/blob/3dece5a/STATUS.md#L14-L75), `STATUS.md`). The measurement: `cases/translatingDroplet2D` put all four sides of its `blockMeshDict` into one `walls` patch. The mesh had no `inlet` and no `outlet`, and OpenFOAM silently discarded the two field entries that matched no mesh patch ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L32)). The case was a closed slip box, and the pressure projection annihilated the uniform stream on step 1: the first-step local continuity error was 1.00e-05 against 2.32e-20 on the repaired mesh, `max|U-U0|` at step 1 was 1.031e-01 against 2.159e-03, `mean|U-U0|` was 5.011e-02 against 1.160e-05, and the whole-domain mean of Ux was -7.3e-07, as a closed box must have ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L34-L46)). So `maxMagUPrime = max|U - (U0,0,0)|` measured the annihilated free stream at about 2 U0 from the first step, not a spurious current. Fixed in commit [440107f](https://github.com/leia-openfoam/leia/commit/440107f) (inlet left, outlet right, walls top and bottom); the pre-registered repair gate `config/translatingFreeStreamGate2D` passes ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L48-L53), [config header](https://github.com/leia-openfoam/leia/blob/8867581/config/translatingFreeStreamGate2D.yaml#L1-L21)). Scope of the void: every number the case produced before the fix, on every arm, in every study, including the metrics that look unaffected by the boundary condition ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L661-L673)).

## The claim, and where it lived

- [The old STATUS section 0](https://github.com/leia-openfoam/leia/blob/3dece5a/STATUS.md#L14-L75), `STATUS.md` at commit 3dece5a (the parent of the fix): "the translating droplet now runs", the pairing table, the "best configuration measured" table, "why" (the mass residual times U0).
- [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L65), `STATUS.md`: rewritten in place as "READ THIS FIRST". Two later subsections still carried closed-box results and are marked RETRACTED 2026-09-27: [the curvature chain is exonerated](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L248-L259) and [leading open defect](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L260-L269).
- [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3566-L3572): the list of the studies the void removes.
- `cases/default.parameter`: the [`MASS_FLUX` rationale](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L668-L672) ("WHY IT WORKS", the 4.8e-02 residual), the [`MASS_FLUX_BOUND_RHO` basis](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L796-L802) (rho to -72.28) and the [`RHO_DDT_SCHEME` pairing tables](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L938-L947).
- `METHOD.md`: [section 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L289-L299), the 8.1 rows [`RHO_DDT_SCHEME`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L372), [`MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L387), [`CURVATURE_EXTENSION`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391) and [`MASS_FLUX_BOUND_RHO`](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L392), and [section 10](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L807-L815).
- Eight curated tables of the method-comparison theme, now in [`VOID_closedBox_20260902/`](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/VOID_closedBox_20260902/README.md) with a README: `alphaFTest_{donorPlane,averagedPlanes}`, `bestConfigTranslating2D`, `rhoBoundGate2D`, `rhoDdtGate2D`, `rhoLENTGate2D`, `translatingMap2D`, `wellBalancedTranslating2D`.
- The SL article, `docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`, [sec:translating, the free-stream footnote](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1931-L1942) `sec:translating`: records the defect and states that every result of the section was re-run.

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| The uniform stream is a discrete steady state of the case. | First-step local continuity error 1.00e-05 (closed box) against 2.32e-20 (repaired); `max\|U-U0\|` at step 1 1.031e-01 against 2.159e-03; `maxMagU` 0.05987 against 0.05207 = U0. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L34-L42) |
| The ambient carries the droplet at U0. | Ambient mean Ux at t = 0.02 s: -2.04e-03; droplet mean Ux +6.14e-02; whole-domain mean -7.3e-07, the zero net flux of a closed box. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L40-L46) |
| `maxMagUPrime` measures a spurious current. | It measured the annihilated free stream at about 2 U0 from step 1 on every arm of every study on the case. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L44-L46), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L684-L688) |
| The solver would have refused a field with a patch the mesh does not have. | OpenFOAM errors when a MESH patch is missing from a field; it ignores a FIELD entry that matches no mesh patch. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L23-L28) |
| The equal-density control was still valid, because `ddt(rho) = 0` and `div(rhoPhi) = rho div(phi) = 0` at ratio 1. | The control was stopped by job id. A setup wrong in one way is not assumed wrong in only that way. | Author decision, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L63-L65), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L690-L694) |
| The repaired case has a real parasitic current. | `max\|U-U0\|` settles near 5e-03; `mean\|U-U0\|` drifts 1.2e-05 to 7.6e-05 over 200 steps, about 50x below the artefact that hid it. | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L50-L53) |

## What survives

Studies that did not run on `translatingDroplet2D` ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L229-L246)):

1. `rhoLENT` on the stationary droplet: `rhoLENTStationary2D`, six arms, +1.0 / -22 / +0.1 % at N = 32 / 64 / 128, volume and shape to three digits. `rhoLENT` stays the `MASS_FLUX` default, see [[decisions/mass-flux-rholent]].
2. `geometricFaceDensity` carries a mass residual of 0.04 to 0.56 relative on the stationary droplet and does no harm there. That measurement is now the motivation for the repaired matrix, not a result about it.
3. The 2D stationary ladders on the shared configuration: `cellCentreInverse` lowers the unabsorbed capillary residual 4.60 / 3.81 / 1.55 / 1.49x at N = 32 to 256.
4. The kinematic transport gates (prescribed velocity, no momentum solve).

What replaced the voided results, on the repaired mesh:

1. `kickOriginGate2D`: the step-1 kick is independent of U0 to 0.4 %, and exact curvature removes it by a factor 1.27e+06 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L74-L102)).
2. `amplifierGate2D` and `amplifierGateEqualRho2D`: the exact-kappa arms stay bounded at every combination of U0 and density ratio ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L190-L227)).
3. `translatingRepaired2D` and `translatingRepairedEqualRho2D`, 16 arms ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L67-L72)): the density ratio is the mechanism, see [[retractions/mass-momentum-consistency-dominant-term]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: [section 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L14-L65) rewritten with the void first; the two [RETRACTED 2026-09-27 subsections](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L248-L269); [11.13 finding 1](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3566-L3572); the [MARKED 2026-09-28 bullet](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L66-L70) (commit d1e3414) for `alphaFTimeLevelTranslating2D`, `massFluxComparison2D` and `translatingClearOutlet2D`, all run on 2026-09-02 before the fix.
- [x] `METHOD.md`: [6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L289-L299), the [8.1 rows](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L372) (`RHO_DDT_SCHEME`, `MASS_FLUX`, `CURVATURE_EXTENSION`, `MASS_FLUX_BOUND_RHO`) and [10](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L807-L815), each marked CORRECTED or RETRACTED 2026-09-27.
- [x] `cases/default.parameter`: CORRECTED 2026-09-27 at [`MASS_FLUX`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L680-L684), [`MASS_FLUX_BOUND_RHO`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L796-L802) and [`RHO_DDT_SCHEME`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L942-L947); CORRECTED 2026-09-28 at the [4.8e-02 residual](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L675-L676) and the [VOID marker](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L821-L823) of the time-level and projection block (commit d1e3414).
- [x] `CLAUDE.md`: the rule ["A wrong setup voids its data"](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L661-L705) and the one-line mesh-patch check, see [[concepts/wrong-setup-voids]].
- [x] The study directories renamed `_VOID_closedBox_20260902` on the cluster; the eight tables moved to the [VOID folder](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/data/tables/VOID_closedBox_20260902/README.md); the equal-density control stopped by job id.
- [x] The SL article: the [free-stream footnote](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1931-L1942).
- [x] The line in [[retraction-log]].

Still missing:

- [ ] `translatingLadder2D` is "most likely" void (committed on the morning of the fix, no record of a run after 440107f), but it is neither renamed nor marked ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3570-L3572)).
- [ ] `MASS_FLUX_BOUND_RHO true` has no valid measured basis after the fix ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L392)), see [[decisions/mass-flux-bound-rho]] and [[concepts/bound-rho]].
- [ ] `CURVATURE_EXTENSION` on the translating droplet is undecided: the only `none` against `cellCentreInverse` comparison was `bestConfigTranslating2D` ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391)).
- [ ] The [`MASS_RESIDUAL_DIAGNOSTIC` comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L805-L812) still describes the U0 times residual mechanism without a marker.

## Related

Hubs: [[hubs/mass-flux]], [[hubs/verification]]. Siblings: [[cases/translating-droplet]], [[concepts/wrong-setup-voids]], [[concepts/rholent-mass-flux]], [[concepts/bound-rho]], [[concepts/alphaf-source-donor-plane]], [[concepts/ddt-scheme-pairing-bdf2]], [[concepts/density-ratio-amplifier]], [[decisions/mass-flux-rholent]], [[decisions/mass-flux-bound-rho]], [[retractions/mass-momentum-consistency-dominant-term]], [[retractions/late-translating-instability-is-the-outlet]], [[retractions/polyhedral-popinet-3d-mesh-defect]].

## Log

### 2026-09-28
Written from STATUS 0 and 11.13, CLAUDE.md, cases/default.parameter and the commit 440107f. Voided 2026-09-02. Entered in [[retraction-log#2026-09]].
