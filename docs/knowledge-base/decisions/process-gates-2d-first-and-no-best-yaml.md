---
title: "Process: the 2D gate before 3D, no config/best.yaml, METHOD.md in the same commit"
description: "Three process rules settled 2026-09-09 and 2026-09-26: any method change runs the 2D method gate before any other coupled study and the 3D gate only after a 2D pass; the best configuration lives in METHOD.md and the .parameter layers and never in a config/best.yaml base file, because Snakemake merges config files shallowly; a gate that changes a setting updates METHOD.md in the same commit"
aliases: [2D gate before 3D, no best.yaml]
kind: decision
status: settled
part: verification
tags: [decision, part/verification]
date: 2026-09-28
date_settled: 2026-09-26
decided_by: [config/gates/methodGate2D.yaml, config/gates/methodGate3D.yaml, "author decision 2026-09-09", "author decision 2026-09-26"]
code: [config/gates/methodGate2D.yaml, config/gates/methodGate3D.yaml, workflow/scripts/render_gate_configs.py, workflow/Snakefile.gate, workflow/scripts/make_gate_summary.py, cases/default.parameter]
sources: ["CLAUDE best-configuration section (L150-L187)", "CLAUDE method gates (L647-L687)", "METHOD L1-L29", "STATUS 11.5 (L3300-L3306)", "STATUS 11.13 (L3563-L3575)"]
---
# Process: the 2D gate before 3D, no config/best.yaml, METHOD.md in the same commit

> Three rules of process, each written from an incident. (1) Any change to the level-set method (a new model, a changed default, a composition root that is not inert) runs the 2D method gate before any other coupled study; the 3D gate runs only after a 2D pass; the laptop runs the unit tests and the 4-rank smoke ([CLAUDE L647-L687](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L647-L687), d91db67, 2026-09-26). (2) The best configuration lives in METHOD.md section 8.1 and in its executable form, the `.parameter` layers plus the study's `axes_override`; a `config/best.yaml` base file is forbidden, because Snakemake merges multiple `--configfile` arguments shallowly: measured on 2026-09-09, a base setting `SL_FIT` plus a study setting only `N_CELLS` rendered a case that had silently lost `SL_FIT` ([CLAUDE L184-L187](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L184-L187), [METHOD L25-L29](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L25-L29), 40b1093). (3) A gate that changes a setting updates METHOD.md in the same commit, and every row names the config and the number that decided it ([CLAUDE L150-L157](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L150-L157)). The guide itself holds no case-specific method settings since 2026-09-27, when the translating arm of the gate inherited a global default that the guide had described as case-specific ([CLAUDE L172-L178](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L172-L178)).

## The question

How does a number become a decision, and how does a decision stay executable? Three failure modes drove the rules: a method that gave good advection and bad hydrodynamics, which only the coupled arms show ([CLAUDE L657-L659](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L657-L659)); a base config that drops axes without warning ([METHOD L25-L29](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L25-L29)); and settings that nobody could trace to a measurement, folklore ([METHOD L3-L7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L3-L7)). The 2D gate is one study for every method: a candidate enters as a set of case tokens, `config/candidates/<name>.yaml` with its pre-registered read-out, and runs the exact 1D arm first, then the 2D shear, stationary, translating and oscillating arms at np 4 against the `baseline` candidate on the same commit and binaries ([methodGate2D.yaml L1-L30](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L1-L30), [[concepts/method-gates]]).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| two `--configfile` arguments, base with `SL_FIT`, study with `N_CELLS` only | rendered case | `SL_FIT` silently lost: the study's `axes_override` replaced the base's entirely | MEASURED 2026-09-09, [METHOD L25-L29](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L25-L29) |
| `VELOCITY_EXTENSION closestPoint` as a candidate on the semi-Lagrangian line | effect | a no-op, because the `projectedFlux` trace ignores the extension; the verdict fails a candidate whose every CSV equals the baseline's | MEASURED 2026-09-26, [CLAUDE L681-L685](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L681-L685) |
| the gate's translating arm, 2026-09-27 | `CURVATURE_EXTENSION` | inherited the global default `cellCentreInverse` while the guide called the setting case-specific; the gate now sets it per arm | MEASURED, [CLAUDE L172-L178](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L172-L178), [methodGate2D.yaml L57-L63](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L57-L63) |
| the first 2D gate campaign, 2026-09-26 to 27 | verdicts | every candidate FAIL on the kinematic shear arm alone before the fixes; the kinematic arms byte-identical (177 pairs) after them | MEASURED, [STATUS L3762-L3772](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3762-L3772), [[studies/method-gate-2d-campaign-2026-09]] |
| `richardson.py --self-test` | checks | 64 PASS: orders 1, 2, 3 recovered to 1e-8, the extrapolated value to 1e-10 | MEASURED, [STATUS L3300-L3306](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3300-L3306) |

Pre-registered verdict of the gate ([methodGate2D.yaml L19-L26](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L19-L26)): every baseline case completes; no regression beyond the tolerance in any vector metric at the finest rung and no order below the baseline's minus the tolerance; the candidate's own target; a candidate whose every CSV equals the baseline's fails; the exact 1D arm passes first. The ladder is a Richardson ladder: cells x2 per rung in 2D, h ratio at least 1.3 in 3D, three rungs at least, R/h >= 10 at the first rung, the same parity of N ([CLAUDE L670-L676](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L670-L676), [[concepts/richardson-ladders-and-orders]]).

## What it does not cover

1. The 3D gate has not run ([STATUS L3739](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3739)).
2. Two blind spots of the 2D gate found on 2026-09-28: the exact 1D closed form is `None` for every candidate, and the centred band metric is blind to the cell-scale mode ([[concepts/method-gates]], [[hubs/verification]]).
3. The record correction of 2026-09-27 (METHOD 4.1, 4.3, 6, 8.1 and 10; three token comments; five STATUS places) shows the same-commit rule was not kept between 2026-07-31 and 2026-09-27 ([STATUS L3571-L3575](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3571-L3575)); this knowledge base is the layer that now mirrors section 8.1.
4. The gate's coupling block copies the measured best values explicitly on the user's instruction of 2026-09-27; it does not replace the `.parameter` layers ([methodGate2D.yaml L57-L63](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L57-L63)).

## Related

[[hubs/verification]] - [[concepts/method-gates]] - [[concepts/richardson-ladders-and-orders]] - [[concepts/bit-identity-and-inertness-gates]] - [[concepts/error-vector-and-read-out-instants]] - [[concepts/seam-checks-and-decomposition-invariance]] - [[decisions/curvature-extension-cell-centre-inverse]] - [[decisions/mesh-family-hexahedral]] - [[cases/exact-1d-stretch]] - [[studies/method-gate-2d-campaign-2026-09]] - [[conventions]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED: no best.yaml and METHOD.md in the same commit on 2026-09-09; the gates on 2026-09-26. Entered in [[decision-log#2026-09]].
