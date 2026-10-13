---
title: "leia knowledge base"
description: "The front door - the moving parts of the leia level-set methods, what was decided, what was retracted, and why"
kind: index
status: settled
part: all
tags: [index]
date: 2026-09-28
---
# leia knowledge base

The concise, cross-linked record of the leia level-set methods for two-phase flow in OpenFOAM:
what each sub-algorithm is, what was decided and on which measurement, what was retracted, why
something failed or why we think so, and what is open. Start with the hub of the part you work
on, then follow its links. The details live in the pre-prints, the decks, `STATUS.md` and the
code; every note links to them.

## The moving parts

| part | hub | the question it answers |
|---|---|---|
| Interface advection and the phase indicator | [[hubs/advection]] | How is the level set transported, and how is the phase indicator built from it? |
| Viscosity | [[hubs/viscosity]] | Which face viscosity is consistent with the interface representation? |
| Surface tension | [[hubs/surface-tension]] | How is the capillary force built, and why does the stationary droplet not stay at rest? |
| Mass flux and density-ratio consistency | [[hubs/mass-flux]] | How is the mass flux made consistent with the momentum flux at a high density ratio? |
| Gradient control | [[hubs/gradient-control]] | Can source terms or a velocity extension keep the level set a signed distance without a reinitialisation? |
| Verification | [[hubs/verification]] | How is a claim measured, gated, retracted and preserved here? |

[[hubs/method-lines]] lists the seven method lines with their status.

## The record

- [[decision-log]]: one line per settled decision, in time order.
- [[retraction-log]]: one line per retracted, voided or corrected claim, in time order.
- [[negative-log]]: one line per route that was tried and failed, with the condition under which
  it may be re-opened.
- [[decisions/process-operating-rules]]: the operating rules that came from incidents.
- [[sessions/current]]: the living handover; the dated notes in `sessions/` are frozen.
- [[conventions]]: how to write here.
- The graph, interactive in 3D: [graph3d/graph.htm](graph3d/graph.htm) (the same links as the graph view in the sidebar).

## Query it from a terminal

The vault is the permanent memory of every agent and every session; it travels with git. From the
repository root (standard library only):

```bash
python3 docs/knowledge-base/.quartz/kb.py search gradient control checkerboard
python3 docs/knowledge-base/.quartz/kb.py negative clip
python3 docs/knowledge-base/.quartz/kb.py show decisions/psi-filter-none
python3 docs/knowledge-base/.quartz/kb.py since 2026-10-01
```

## Where the details live

- Pre-prints and decks: `docs/<theme>/` in the repository; the site serves them under
  `preprints/` and `decks/`.
- `STATUS.md` (the lab notebook), `METHOD.md` (the best configuration and its evidence),
  `CLAUDE.md` (the rules), the plan documents under `docs/`.
- Code: `src/leiaLevelSet/` (one core library and eight method libraries), `applications/solvers/`.

## Log

### 2026-09-28
Created.

### 2026-10-13
Linked the negative log and the operating rules; the query tool.
