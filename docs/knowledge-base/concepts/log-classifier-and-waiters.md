---
title: "Reading solver logs and waiting on them"
description: "One classifier reads every OpenFOAM solver log (foam_log_state.sh, 2026-08-26); a bare grep for 'Floating point exception' matches the trapFpe banner of every healthy run and failed three times; every waiter has a timeout and an exit condition that can occur (2026-08-31)"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-29
date_settled: 2026-08-31
decided_by: [author decision 2026-08-26, author decision 2026-08-31]
code: [workflow/scripts/foam_log_state.sh, workflow/Snakefile, workflow/scripts/aggregate.py, workflow/scripts/make_gate_summary.py, workflow/scripts/guard_finished_cases.py]
sources: [CLAUDE log classifier section, CLAUDE waiter section, STATUS 4 process failures 2026-08-20, STATUS 9.3, STATUS 10.4, STATUS 11.2, SLURM.md log section]
---
# Reading solver logs and waiting on them

> Verdict (2026-09-29). Nobody writes a new grep against a solver log. Every script and every interactive check calls `workflow/scripts/foam_log_state.sh`, which returns one of six states, each with its own exit code ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L320-L330), [`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L18-L34)). The rule exists because every OpenFOAM solver prints `trapFpe: Floating point exception trapping enabled` in its startup banner ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L309-L311)). A grep for a bare `Floating point exception` therefore marks every healthy run as diverged. That false positive occurred three times: in the first Snakefile classifier, in a dt-bisection probe, and in a completion waiter. The waiter declared three live studies finished ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L313-L318)). Completeness comes from the terminating `End` line; liveness comes from the log mtime and the step count, never from `squeue` alone ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L332-L339)). A related mistake, completeness read from a post-processing output, deleted a live interFoam arm at step 90404 of 106689, about 16 h of work ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1059-L1066)). Every hand-written waiter has a `timeout`, an exit condition that can occur, and a pattern that its own command line does not contain ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L341-L372)).

## What it is

The classifier reads one log and prints one line: the state, then `steps=`, `age=` and `nprocs=`. Consumers read the first word and the exit code; the fields after the state are additive ([`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L18-L38)).

| state | exit code | condition in the log | meaning |
|---|---|---|---|
| COMPLETED | 0 | a line `End` (pattern `^End$`) | the run finished |
| RUNNING | 1 | none of the other states, and the log grew in the last `--stall` seconds | the run is alive |
| STALLED | 2 | none of the other states, and the log did not grow for `--stall` seconds (default 900) | examine the process before any action |
| DIVERGED | 3 | `sigFpe::sigHandler`, `sigSegv::sigHandler`, `Floating point exception (core dumped)`, `Segmentation fault (core dumped)` or `FOAM FATAL` | the solver died; the blow-up is a result |
| LAUNCH_FAILURE | 4 | an MPI or SLURM launch error, and zero time steps | the run never happened; it is never recorded as a result |
| MISSING | 5 | no log file | nothing to read |

The script tests the states in this order: COMPLETED, LAUNCH_FAILURE, DIVERGED, STALLED, RUNNING ([`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L60-L82)). It accepts LAUNCH_FAILURE only for a log with zero steps, because an aborted launch also writes MPI noise. The fields:

- `steps=` is the number of `Time = ` lines.
- `age=` is the time in seconds since the last write to the log, from its mtime.
- `nprocs=` is the rank count of the log header `nProcs : N`, and 1 for a serial log.

`--wait` repeats the classification every `--poll` seconds (default 30) while the state is RUNNING. It exits with the code of the first other state ([`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L50), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L84-L90)). A log that stops growing becomes STALLED after `--stall` seconds, so `--wait` ends when the solver stops writing (DERIVED from the loop).

The patterns are the same strings as in the solve rule of `workflow/Snakefile`. A change goes into both files in the same commit ([`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L40-L46), [`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L277-L288)).

The waiter rules ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L341-L372)):

1. Use `foam_log_state.sh --wait` in place of a hand-written loop.
2. Give every hand-written poll loop a `timeout`. Without one, `until <cond>; do sleep N; done` polls until the session ends if the condition never occurs.
3. Never match a pattern that your own command line contains. `! pgrep -f "configfile config/<study>"` matched the shell that ran it, so the negation stayed false. `pkill -f <pattern>` kills its own shell (exit 144) for the same reason.
4. Make sure that the pattern can occur. Two waiters waited for text that ssh did not deliver ([[concepts/cluster-provenance-and-binaries]]).
5. Identify an orphan from a ledger of what you started, by id. Never sweep by pattern or by user: other sessions share the laptop and the cluster account.

One more trap: `pgrep -x` compares only the first 15 characters of the process name ([pgrep(1)](https://man7.org/linux/man-pages/man1/pgrep.1.html), NOTES). `pgrep -x leiaSemiLagrangianLevelSetTwoPhaseFoam` therefore matches no process. A liveness check with the full name then reports a live solver as absent (DERIVED; the repository record has no incident for it).

## Why it matters

- A check that makes a live run look dead or finished lets a cleanup delete it. On 2026-08-20 a completeness test on a post-processing file classified every running interFoam arm as dead ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1059-L1066)).
- A false LAUNCH_FAILURE discards physics. Until 2026-08-15 the solve rule matched a bare `srun: error`, which srun also prints when the solver dies. Snakemake then deleted the partial time series of eleven of twelve footEval runs ([`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L265-L276)).
- The opposite error records an infrastructure fault as physics. SLURM spread an np 8 allocation over two nodes, and the rule logged six "divergences" ([`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L243-L251)).
- Until 2026-08-19 the rule labelled every failure "diverged". This included a 3D arm killed by its solve limit at 85 % of the run and a heap corruption of interFlow at teardown ([`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L253-L263)).
- An empty `squeue` is not an empty queue. On 2026-08-20 the primary controller `mssd0001` was down, `squeue` returned zero rows, and every driver was running; `scontrol ping` separates the two cases ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1054-L1058)).
- `age=` comes from the mtime. A sweep gave every file under `/work/scratch/tm83tomy` a new mtime between 10:59 and 11:12 on 2026-09-22. Every age on that scratch meant nothing for that day ([STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2781-L2784)).

## Where in the code

- `workflow/scripts/foam_log_state.sh`: the classifier ([L1-L91](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L1-L91)). Created in [d6bccd0](https://github.com/leia-openfoam/leia/commit/d6bccd0) (2026-08-26); the zero-step count fixed in [4940440](https://github.com/leia-openfoam/leia/commit/4940440) (2026-09-23); `nprocs=` added in [f8aaa8e](https://github.com/leia-openfoam/leia/commit/f8aaa8e) (2026-09-26, [STATUS 11.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3272-L3274)).
- `workflow/Snakefile`, the solve rule: a diverged run is a result ([L209-L222](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L209-L222)); the launch record `<case>/.leia_launch` stores the classifier line of the finished log ([L223-L241](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L223-L241)).
- `workflow/scripts/aggregate.py` reads the launch record into the columns `logState`, `steps` and `nRanks` ([L97-L125](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/aggregate.py#L97-L125)).
- `workflow/scripts/make_gate_summary.py`: a gate rung whose baseline case is not COMPLETED is INVALID (criterion 1). A candidate that does not complete where the baseline completed fails ([L388-L402](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L388-L402), [[concepts/method-gates]]).
- `workflow/scripts/guard_finished_cases.py`: refuses a launch whose `solve` job targets a case with a COMPLETED or DIVERGED log, exit 2 ([L1-L31](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/guard_finished_cases.py#L1-L31)).
- Other callers: `advection_convergence_table.py`, `value_bound_advection_census.py`, `solve_decomposed_arm.sbatch`, `advection_poly_ladder.sbatch`, `advect_bound_arm.sh`.
- The site-independent copy of the rules: [SLURM.md](https://github.com/leia-openfoam/leia/blob/d1e3414/SLURM.md#L312-L347).

## Evidence

| claim | number | where |
|---|---|---|
| the trapFpe false positive recurred | three times: the Snakefile classifier, a dt-bisection probe, a completion waiter that declared three live studies finished | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L313-L318), MEASURED |
| completeness from a post-processing output deleted a live arm | `interFoamDroplet2D_00003` removed at step 90404 of 106689, about 16 h lost (2026-08-20) | [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1059-L1066), MEASURED |
| an empty `squeue` during a controller outage | zero rows while every driver was running (2026-08-20) | [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1054-L1058), MEASURED |
| a loose launch-failure pattern discarded physics | 11 of 12 footEval runs lost (pattern in use until 2026-08-15) | [`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L265-L276), MEASURED |
| a launch failure recorded as physics | 6 false "divergences" from one np 8 allocation over two nodes | [`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L243-L251), MEASURED |
| unbounded waiters | three still polling hours after their work finished (2026-08-31) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L347-L352), MEASURED |
| the zero-step count | `grep -c` printed 0 twice for a launch-failure log, and the LAUNCH_FAILURE test stopped with an error (2026-09-23) | [`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L63-L65), [STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3111), MEASURED |
| the launch guard reads the classifier | `sdplsConv2Dvortex`: 6 of 6 listed `solve` jobs targeted COMPLETED cases (700 steps), exit 2; after the touch pass 0 jobs, exit 0 (2026-09-24) | [STATUS 10.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3165-L3167), MEASURED |
| `--wait` ends when the log stops growing | STALLED after `--stall` seconds, default 900 | [`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L78-L90), DERIVED |
| `pgrep -x` with a long solver name matches nothing | the name comparison uses 15 characters | [pgrep(1)](https://man7.org/linux/man-pages/man1/pgrep.1.html), DERIVED |

## Why it failed, or why we think so

Our reading of the incidents above: each check read a secondary sign in place of the state of the run. A banner line is not a death signature. A file that exists only after the run is not a completeness test. An empty `squeue` is not an empty queue. The classifier reads the three primary signs: the `End` line, the death signatures, and the growth of the log.

## Decisions

- 2026-08-26: one classifier for solver logs, and the rule that mandates it ([d6bccd0](https://github.com/leia-openfoam/leia/commit/d6bccd0), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L307-L339)).
- 2026-08-31: the waiter rules ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L341-L372)).
- A DIVERGED run is a result; a LAUNCH_FAILURE is not ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L339), [`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L209-L222)).

## Open questions

1. The helper and the Snakefile use the same patterns in a different order. The Snakefile tests the death signatures first; the helper tests the launch errors first when the log has zero steps. For a zero-step log with both kinds of line, the helper returns LAUNCH_FAILURE and the Snakefile records a divergence (DERIVED from [`foam_log_state.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/foam_log_state.sh#L71-L77) and [`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L277-L288)). No such log is recorded.

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/cluster-provenance-and-binaries]] (ssh output, the job ledger, the touch sweep), [[concepts/seam-checks-and-decomposition-invariance]] (a deadlock leaves a live job with a truncated log), [[concepts/method-gates]] (criterion 1 reads the classifier), [[concepts/bit-identity-and-inertness-gates]], [[concepts/wrong-setup-voids]].

## Log

### 2026-09-29
Created.
