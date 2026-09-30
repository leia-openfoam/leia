---
title: "Bit-identity and inertness gates"
description: "A change that claims to change nothing runs against the pre-change state and must reproduce every metric CSV at tolerance 0 with the two clock columns skipped; a residualControl block with tolerance 0 was not inert (2026-08-27), and a finished case's 0/ is not its initial state (2026-09-05)"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-29
date_settled: 2026-09-10
decided_by: [author decision 2026-08-27, cases/stationaryDroplet2D/system/residualControl.off, config/stationaryDroplet3DbitIdentity.yaml, docs/plan-library-split-and-build-policy.md]
code: [workflow/scripts/compare_metrics_csv.py, workflow/scripts/make_gate_summary.py, cases/default.parameter, cases/stationaryDroplet2D/system/residualControl.off, config/stationaryDroplet3DbitIdentity.yaml]
sources: [CLAUDE nothing-is-inert section, CLAUDE provenance section, CLAUDE regression set section, STATUS 4 poly regression 2026-09-05, STATUS 10.1, STATUS 10.3, STATUS 11.7, STATUS 11.8, STATUS 11.14, STATUS 11.15, plan-library-split WP3 gate]
---
# Bit-identity and inertness gates

> Verdict (2026-09-29). Nothing is inert until a measurement says so. Each new token, refactor, composition root or changed default runs against the pre-change state on a committed case, before any physics ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L764-L773), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L89-L91), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L178-L180)). The comparison is `compare_metrics_csv.py --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`. A `cmp` of the whole file reports a difference for two runs that are identical in every physical quantity ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L637-L644)). The rule has a measured basis. A `residualControl` block with tolerance 0 changed which corrector PIMPLE treats as final. It moved a 30-step N = 128 droplet run in the 8th significant digit (2026-08-27, [`residualControl.off`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/residualControl.off#L3-L8)). A re-run started from a finished case's `0/` differed from step 1 (2026-09-05). The reason: the solver writes its projected pressure and its recomputed phase indicator back into time 0 ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L822-L832)). The large gates of this repository passed with 0 differences. The library split reproduced 209 CSV pairs ([STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3065-L3094)). The Phase C extension points reproduced 53 cases ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3361-L3387)). The kinematic arms of the fixed 2D gate reproduced 177 CSV pairs ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3897-L3899)).

## What it is

The procedure ([plan WP3 gate](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-library-split-and-build-policy.md#L261-L281)):

1. Preserve the pre-change study, or build the pre-change binaries in a separate worktree or clone. Phase C used a reference worktree at 0c7079c ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3361-L3364)).
2. Run the same committed cases with both sets of binaries, on the same mesh and the same decomposition.
3. Compare every case CSV pairwise with the comparator at tolerance 0. Never compare only the aggregated tables.
4. PASS means every pair is identical. Any difference stops the push: "a refactor that changes one number is not a refactor".

The comparator `workflow/scripts/compare_metrics_csv.py` ([L1-L12](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/compare_metrics_csv.py#L1-L12)):

- For every numeric column of both files it computes the maximum over the rows of `|a - b| / max(|a|, |b|, floor)`, with `floor = 1e-14` by default.
- The default tolerance is 1e-10 ([L29-L32](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/compare_metrics_csv.py#L29-L32)). With `--tol 0`, any nonzero difference in any numeric column fails (DERIVED from [L55](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/compare_metrics_csv.py#L55)).
- Unequal row counts stop the comparison with exit code 1, because runs with unequal step counts are not comparable ([L36-L38](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/compare_metrics_csv.py#L36-L38)).
- Exit codes: 0 PASS, 2 FAIL, 1 structural mismatch.
- The two skipped columns are wall-clock times from [`advectionErrorsCsv.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/advectionErrorsCsv.H#L52-L53); they always differ.

A tolerance-0 CSV comparison sees the digits that the solver writes, not the fields in memory (DERIVED). The Phase C gate therefore also compared every rendered dictionary and every file of every written time directory byte for byte ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3361-L3364)).

The inertness rules:

- A new template token ships with its inert default in the same commit. `SDPLS_FLUX_BLEND` came without a default in 4923800, and the materialization of every study that touches `cases/2Dvortex` stopped ([`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L577-L583)).
- A new model is runtime-selectable and inert by default ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L696-L698)).
- A new composition root in a solver ships with the inert default, in its own commit, with a bit-identity gate ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L89-L91)).
- A change of the best configuration has a bit-identity gate, because every study that does not override that axis inherits it ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L178-L180)).
- A dictionary entry whose presence changes the algorithm is opt-in. `PIMPLE_RESIDUAL_CONTROL_FILE residualControl.off` keeps the PIMPLE dictionary byte-equal to its historical form ([`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L597-L608)).

The inverse check. The method gates fail a candidate whose every CSV equals the baseline's at tolerance 0: its tokens were not consumed ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L350-L367), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L435-L437)). On 2026-09-26 `VELOCITY_EXTENSION closestPoint` on the semi-Lagrangian line was such a no-op, because the default `projectedFlux` trace ignores the extension ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L682-L686)).

The initial state. A bit-identity re-run regenerates `0/` exactly as the workflow does: `0.org` plus the pre-processing (`leiaSetFields`). It never copies the finished `0/`. At startup the two-phase solver writes its projected initial pressure and its recomputed phase indicator into time 0 of the processor directories. `reconstructPar -withZero` then copies both into the serial `0/`; `psi` and `U` stay unchanged ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L822-L832)).

## Why it matters

- "This default changes nothing" is a claim, not a fact. The mere presence of the `residualControl` block changed which solver settings the last pass used (`U`/`p_rgh` against `*Final`) ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L768-L771)).
- A difference from the wrong initial state looks like a code change. The re-run from a finished `0/` was nearly blamed on a library change, which a tolerance-0 control then proved bit-inert ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L825-L829)).
- One solver cannot stand in for another. The `slValueBound` refactor was bit-identical on the coupled solver over 1563 steps and 8 arms. Nobody ran the advection solver, although it reaches the same rewritten function ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L600-L604), [[concepts/advection-regression-set]]).
- Bit identity separates the code from the machine. The serial 1D case gave a maximum relative difference of 0.0 between gcc 11.5 on the cluster and gcc 13.3 on the laptop ([STATUS 10.1](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2960-L2964)).

## Where in the code

- `workflow/scripts/compare_metrics_csv.py`: the comparator ([L1-L67](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/compare_metrics_csv.py#L1-L67)), created in [24223a2](https://github.com/leia-openfoam/leia/commit/24223a2) (2026-09-04).
- `workflow/scripts/make_gate_summary.py`: the comparator and the clock columns ([L46-L47](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L46-L47)); the no-effect check ([L350-L367](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L350-L367)).
- `cases/*/system/residualControl.off` and the token `PIMPLE_RESIDUAL_CONTROL_FILE` ([`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L597-L608)).
- `config/stationaryDroplet3DbitIdentity.yaml`: the committed bit-identity case of the semi-Lagrangian two-phase solver, N = 30 on the 6R box, 27 000 cells, 65 steps ([L1-L7](https://github.com/leia-openfoam/leia/blob/d1e3414/config/stationaryDroplet3DbitIdentity.yaml#L1-L7)).
- The Phase C field and dictionary comparison used `bitid2.sh` and `bitid_compare.py` ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3364)). Neither script is in the tree of d1e3414.

## Evidence

| claim | number | where |
|---|---|---|
| a `residualControl` block with tolerance 0 is not inert | a 30-step N = 128 droplet run moved in the 8th significant digit (2026-08-27) | [`residualControl.off`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/residualControl.off#L3-L8), MEASURED |
| a valueless token stops the materialization | every study that touches `cases/2Dvortex` (commit 4923800) | [`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L577-L583), MEASURED |
| a finished `0/` is not the initial state | first pressure solve: initial residual 6e-7 in 6 iterations against 1 in 27; with the pressure restored, alpha still moved the velocities at 1e-15; from a regenerated `0/`, byte for byte over 458 steps (2026-09-05) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L822-L832), [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1713-L1737), MEASURED |
| the committed SL two-phase case | cmp-identical over 65 steps against the pre-fix CSV (`quadraticPivotTol` 0.3, 2026-09-05) | [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1713-L1714), MEASURED |
| the `slValueBound` family | bit-identical on 8 arms over 1563 steps at np 4, coupled solver only | [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L285-L291), MEASURED |
| per-clone binaries (WP1) | 12 cases, 24 CSV comparisons at tolerance 0 | [STATUS 10.1](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2935-L2941), MEASURED |
| the library split (WP3) | 209 CSV pairs, 0 differences; unit test 89 of 89 | [STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3065-L3094), MEASURED |
| the WP4 and WP6 probe sets | 42 CSV pairs per clone identical to the WP3 baseline | [STATUS 10.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3180-L3181), [STATUS 10.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3223-L3226), MEASURED |
| the Phase C extension points | 53 cases, 0 differences in CSVs, fields and dictionaries | [STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3361-L3387), MEASURED |
| the Phase C6 tokens | render diff 250 of 250 configs, 2 899 changed files, 126 487 added lines, no deleted or changed line; bit identity 53 of 53 | [STATUS 11.8](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3407-L3419), MEASURED |
| the render diff detects a change | it fails on two injected changes (a changed `relTol`, `traceFlux extension`) | [STATUS 11.8](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3410-L3412), MEASURED |
| the two parallel fixes leave serial runs unchanged | bit-identical over 4604 steps, every column (2026-09-27) | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3683-L3684), MEASURED |
| the translating re-runs to 0.05 s | 15 of 15 reproduce the first 0.05 s of their 0.1 s runs; S1 at N = 142 to the same divergence row, 6741 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3782-L3785), MEASURED |
| the kinematic arms of the fixed 2D gate | 177 CSV pairs identical at tolerance 0, as predicted | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3897-L3899), MEASURED |
| a no-op candidate | `VELOCITY_EXTENSION closestPoint` on the SL line: every CSV equal to the baseline's (2026-09-26) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L682-L686), MEASURED |

## Decisions

- 2026-08-27: the presence of a `residualControl` block is opt-in, by a file-level switch ([`residualControl.off`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/stationaryDroplet2D/system/residualControl.off#L1-L8)).
- 2026-09-10: metric CSVs are compared with the committed comparator, never with `cmp` of the whole file ([30cb989](https://github.com/leia-openfoam/leia/commit/30cb989), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L637-L644)).
- 2026-09-05: a bit-identity re-run regenerates `0/` from `0.org` and the pre-processing ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L829-L832)).
- 2026-09-22: a refactor passes only when every pair is identical ([plan](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/plan-library-split-and-build-policy.md#L280-L281)).

## Open questions

1. The field and dictionary comparison of Phase C used two scripts that are not committed. A later gate cannot repeat that part of the comparison from the repository ([STATUS 11.7](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3364)).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/advection-regression-set]] (the rung set that an inert advection change must reproduce), [[concepts/cluster-provenance-and-binaries]] (the build-policy gates, the version stamps), [[concepts/method-gates]] (the no-effect FAIL), [[concepts/log-classifier-and-waiters]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/wrong-setup-voids]].
- Models: [[models/sl-value-bound]] (the family whose gate ran on one solver only), [[models/semi-implicit-capillary-force]] (`off` verified bit-identical, [`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L585-L593)).

## Log

### 2026-09-29
Created.
