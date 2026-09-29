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
   ([gclsTechnicalReport.tex](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex)),
   the two rule sections of CLAUDE.md (the knowledge base as the point of reference; supervisor
   and expert developer), and the record corrections of STATUS.md, METHOD.md and the configs.
2. Gradient control: the next experiments are pre-registered in
   [[concepts/gradient-control-next-experiments]]; none has run.

## What is open, in order

0. Publishing: the `github-pages` environment allows deployments from `main` only; add this
   branch under Settings > Environments > github-pages > Deployment branches, or merge into `main`.
1. Author decisions on the gate repairs ([[concepts/method-gates]]), the next gradient-control
   experiments ([[concepts/gradient-control-open-decisions]]), `boundRho`
   ([[decisions/mass-flux-bound-rho]]), which parallel studies to re-run
   ([[concepts/coupled-face-density-defect]]), the void of the five frozen-density Eulerian studies
   ([[concepts/eulerian-solver-mass-flux-port]]), a box-length token and the oscillating horizon
   ([[cases/translating-droplet]], [[cases/oscillating-droplet]]), and the record inconsistencies the note
   writers found, listed with the corrections to make at the source in [[sessions/sl-session-handover]]
   (among them: the production curvature was never scored on the ellipse gate; a row of the SL
   article's viscous table comes from the frozen-muf run).
2. The third rung of the 40 mm translating box ($N = 200$) on the cluster.
3. Regenerate the `advConv2D*` convergence CSVs on Lichtenberg ([[retractions/advection-orders-3-2-factor]]).

## Build progress of the vault (incremental, restart-safe)

Complete on 2026-09-29: 131 notes (7 hubs, 17 models, 55 concepts, 16 decisions, 13 retractions,
8 cases, 12 studies, 3 sessions), both logs, `check_kb.py` 0 problems. The work stayed
restart-safe through three stops by the spend limit: the manifest `.quartz/manifest.md` fixes every
slug and owner, `.quartz/missing_notes.py` lists what is missing, each writer saves a source digest
in `.quartz/digests/` (git-ignored) before it drafts and each note at once, and the tree is committed
after each batch. To add a batch of notes later:

1. Add the rows to `.quartz/manifest.md` (slug, title, scope, sources, owner).
2. Launch one writer per three to five notes with `.quartz/writing-brief.md`; pin every link to the
   commit whose tree the writer reads.
3. Run `python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base` and
   `python3 docs/knowledge-base/.quartz/check_links.py docs/knowledge-base --words` before the commit.

Link integrity: the first writers read the tree of d1e3414 and pinned many links to 8867581;
`check_links.py` moved 65 pins on decisive evidence and a reviewer decided 159 ambiguous links by
hand (`.quartz/link-audit-2026-09-29.txt`).

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

### 2026-09-29
The vault is complete (131 notes); the link review is done; publishing waits for the Pages setting.
