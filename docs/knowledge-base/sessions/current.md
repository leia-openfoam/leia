---
title: "Current handover"
description: "The living handover, rewritten each sitting: branch, commit, what is being worked on, what is open in order, the traps, where the numbers live"
kind: session
status: open
part: all
tags: [session]
date: 2026-09-28
---
# Current handover

> Branch `feature/gradient-controlled-level-set`, laptop clone `~/OpenFOAM/repos/leia-gcls`, last
> data commit 8867581 (the gcls pre-print), gate data stamp
> `shared-method-config-2026-09-01-192-g1150e68`. Last dated handover:
> [[sessions/2026-09-27-gcls-first-campaign]]. Read [[decision-log#2026-09]] and
> [[retraction-log#2026-09]] from 2026-09-26.

## What is being worked on

1. 2026-09-28: this knowledge base (Quartz site, three.js graph), the separation of the
   semi-Lagrangian topic from the source-term topic (four sections of the gcls pre-print moved to
   the SL article; [[sessions/sl-session-handover]]), the technical report on the first
   gradient-control campaign
   ([gclsTechnicalReport.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex)),
   the two rule sections of CLAUDE.md (the knowledge base as the point of reference; supervisor
   and expert developer), and the record corrections of STATUS.md, METHOD.md and the configs.
2. Gradient control: the next experiments are pre-registered in
   [[concepts/gradient-control-next-experiments]]; none has run.

## What is open, in order

1. Author decisions on the gate repairs ([[concepts/method-gates]]), the next gradient-control
   experiments ([[concepts/gradient-control-open-decisions]]), `boundRho`
   ([[decisions/mass-flux-bound-rho]]), which parallel studies to re-run
   ([[concepts/coupled-face-density-defect]]), the void of the five frozen-density Eulerian studies
   ([[concepts/eulerian-solver-mass-flux-port]]), a box-length token and the oscillating horizon
   ([[cases/translating-droplet]], [[cases/oscillating-droplet]]).
2. The third rung of the 40 mm translating box ($N = 200$) on the cluster.
3. Regenerate the `advConv2D*` convergence CSVs on Lichtenberg ([[retractions/advection-orders-3-2-factor]]).

## Build progress of the vault (incremental, restart-safe)

The vault is written in increments by note writers (agents) that a budget limit can stop at
any time. Nothing is lost between increments: every note is written to disk as soon as it is
done, its source digest goes to `.quartz/digests/<slug>.md` before the draft, and the state is
committed after each batch. To finish or resume:

1. `python3 docs/knowledge-base/.quartz/missing_notes.py docs/knowledge-base` lists the notes
   of `.quartz/manifest.md` that do not exist yet, per owner.
2. Relaunch one writer per owner, or per four notes, with `.quartz/writing-brief.md`,
   `.quartz/manifest.md` and `.quartz/raw-material-2026-09-28.md`; a writer reads the digest
   of a slug if it exists, else the sources; it never rewrites an existing note.
3. `python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base` must print 0
   problems before the final commit; until then the commits are marked incremental.

State on 2026-09-28 22:45: hubs 7/7, sessions 3/3, models 17/17, studies 12/12, concepts
12 of 55, retractions 4 of 13, decisions 0 of 16, cases 0 of 8 (the writers for the missing
groups are running).

## Traps that cost time

- Never `scancel -u`; cancel by id from `.my_jobs` ([[concepts/cluster-provenance-and-binaries]]).
- Never grep a solver log; use `workflow/scripts/foam_log_state.sh` ([[concepts/log-classifier-and-waiters]]).
- A finished `0/` is not the initial state; regenerate from `0.org` ([[concepts/bit-identity-and-inertness-gates]]).
- Check `constant/polyMesh/boundary` before reading any droplet metric ([[concepts/wrong-setup-voids]]).
- Every SL two-phase result on more than one rank before 2026-09-27 carries the two seam defects ([[concepts/coupled-face-density-defect]]).
- The cluster clone `/work/scratch/tm83tomy/leia` is on the feature branch at 8867581; no job of this session runs there.

## Where the numbers live

- The gate and the laptop runs: `docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/` (README, MANIFEST) ([[concepts/data-archive-per-version]]).
- The lab notebook: [STATUS 11](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3218-L4041); the best configuration: [METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L341-L500).
- Raw runs (git-ignored): `runs/gcls-laptop-20260927` on the laptop and on Lichtenberg under `/work/scratch/tm83tomy/leia/runs/`.

## Log

### 2026-09-28
Created with the vault.
