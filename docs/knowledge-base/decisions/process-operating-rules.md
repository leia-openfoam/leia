---
title: "Process: the operating rules that came from incidents"
description: "Fourteen operating rules that lived only in one machine's session memory, now shared through the vault: L2 and L1 never L_inf, the full advection in static cases, the STATUS hand-off, clone ownership, large runs on Lichtenberg, the shell and snakemake traps, the rebuild rule, incremental agents, the replace refs, the exact criterion form, CI without gh, the writing and media preferences; each with the incident that produced it"
aliases: [operating rules, agent rules, session memory rules]
kind: decision
status: settled
part: all
tags: [decision, part/all]
date: 2026-10-13
date_settled: 2026-10-13
decided_by: ["author decisions 2026-08-07 to 2026-10-01 (the date of each rule is in its row)"]
code: [CLAUDE.md, CLUSTER.md, workflow/scripts/foam_log_state.sh, docs/knowledge-base/.quartz/kb.py]
sources: ["the laptop's Claude session memory, 2026-08 to 2026-10", "STATUS 11.20"]
---
# Process: the operating rules that came from incidents

> Settled 2026-10-13. These rules were kept in the session memory of one machine, where no other agent, no cluster session and no other tool can read them. The vault is the permanent memory that every agent shares through git, so they live here now. Each row names the incident that produced the rule; the rules that `CLAUDE.md` already states in full (BDF2, the log classifier, the account-wide `scancel`, the 4-rank gate, the unstructured-FVM constraint, no filtering in production, the research loop) are not repeated here.

## The rules

| rule | why: the incident | since | where else |
|---|---|---|---|
| Report L2 and L1 of every metric, never L_inf, in tables, prose, commits, plots and curated CSVs. | L_inf is a single-cell extremum set by the sub-cell position of the interface; it does not converge even in simple cases on unstructured meshes, and a verdict built on it had to be retracted. | 2026-09-03 | [[concepts/error-vector-and-read-out-instants]] |
| A static case (the stationary droplet) runs the full advection every step, exactly as a moving one; a stabiliser must work identically for a deforming interface. | A shortcut that exploits a static interface makes the static test meaningless and produces models that cannot be used. | 2026-08-07 | [[cases/stationary-droplet]] |
| `STATUS.md` is the hand-off for sessions run from a phone, which use Lichtenberg only: keep it current, with what is running, what was falsified, and where the data are; keep the clones synced. | A phone session cannot inspect a working tree; what is not written in the repository is invisible there. | 2026-08 | [[sessions/current]] |
| Every cluster clone has one owner. Never pull, stash or build in another session's clone; identify your jobs by your own ledger. | Sessions share the `tm83tomy` account and one laptop; a cleanup aimed at one session's work hit another's. | 2026-08-15 | [[concepts/cluster-provenance-and-binaries]], `CLUSTER.md` |
| Large runs go to Lichtenberg; the laptop (24 cores, 15 GB) is for builds, smoke runs and small 4-rank gates, one small pair at a time. | Two polyhedral pairs (473k and 245k cells) filled the laptop's memory while the author worked on it. | 2026-10-01 | [[decisions/flux-projection-polyhedral]] |
| `pgrep -x` needs the 15-character truncated name (`leiaSemiLagrang`); `pgrep -f` matches its own shell. A script with `set -e` or `set -u` dies silently when it sources OpenFOAM's `etc/bashrc`. | A rebuild and five census runs "completed" with empty output. | 2026-09-05 | [[concepts/log-classifier-and-waiters]] |
| Start concurrent snakemake runs in one worktree a few seconds apart, and retry a start that fails on `.snakemake/incomplete/migration_underway`. | 13 of 250 parallel render jobs and one of three bit-identity queues failed on this race; each passed alone. | 2026-09-26 | — |
| After a header change, above all of an inline or virtual member, build from the root (`./Allwmake`), never only the library. | GCC speculative devirtualisation kept the old inlined body in a stale test binary; the test failed against correct code. | 2026-09-26 | — |
| Delegated agents save per unit: a source digest before the draft, each unit to disk when done, a manifest and a "what is missing" script, a commit per batch. | The API spend limit stopped all eight writer agents twice in one day, and those that had drafted nothing lost everything. | 2026-09-28 | `.quartz/manifest.md`, `.quartz/missing_notes.py` |
| Every command that relates `main` to `development` (merge, merge-base, push) runs with `git --no-replace-objects`. | 446 local replace refs hide the history split of 2026-07-09; with them git reports "up to date" where GitHub sees a 2022 merge base and 112 conflicts. `main` was superseded on 2026-09-29. | 2026-09-29 | STATUS 11.20 |
| A pre-registered read-out that invokes an existing criterion quotes how that criterion was applied: the rung, the fit window, the amplitude. | A stricter paraphrase of the ellipse gain criterion nearly gave a false FAIL. | 2026-09-29 | [[decisions/curvature-extension-cell-centre-inverse]] |
| The laptop has no `gh`: read CI states through the public check-runs and annotations API with curl (60 requests per hour), in a bounded poll. | The job logs need admin rights; the annotations carry the CI script's PASS and FAIL lines. | 2026-09-30 | [[concepts/ci-build-tests]] |
| Write fully, not tersely, in ASD-STE100; plain language and a formula for every new quantity; no analogies. | Dense tables and short verdicts cost the author more time to unpack. | 2026-09 | `CLAUDE.md`, communication |
| Animations for the author are GIF, with explicit colour-bar axes and fonts readable in a GIF; mp4 stays for the decks. | mp4 downloads fail for the author. | 2026-09 | — |

## What it does not cover

1. The rules of `CLAUDE.md` and `AGENTS.md`, which every agent loads at the start; this note holds what they did not.
2. One machine's session memory may still hold rules that are not here. A rule that must outlive a session is written in the vault, not only in a memory.

## Related

[[decisions/process-gates-2d-first-and-no-best-yaml]] - [[decisions/kinematic-solver-per-flow-solver]] - [[concepts/log-classifier-and-waiters]] - [[concepts/cluster-provenance-and-binaries]] - [[sessions/current]] - [[conventions]] - [[decision-log]]

## Log

### 2026-10-13
SETTLED: moved from the laptop's session memory into the vault, so that every agent can read it. Entered in [[decision-log#2026-10]].
