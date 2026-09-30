---
title: "Gradient control: the open decisions for the supervisor"
description: "What the supervisor decides before the next campaign: the abandon list, the order of the experiments, the gate repairs, the reference arms, the Tier-I tests, the sequencing of S2, ALG and V1"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [technical report sec 6]
---
# Gradient control: the open decisions for the supervisor

> Open (2026-09-28). The technical report proposes; the supervisor decides; this note records.

1. **Abandon as formulated** (proposed): the halo-limited extension at $R = h$ as a gradient-control
   device; HL2; the soft wall as parameterised ($p = 5$, $\delta_s = 0.08$); FP0 as a production
   candidate (keep as a reference arm); the strain-weight (Taylor / weighted-SDPLS) route until the
   curvature-error seed is smaller ([[concepts/why-the-candidates-failed]]).
2. **The order of the next experiments** ([[concepts/gradient-control-next-experiments]]): E0.1,
   then E1.1 to E1.3, then E2.1 and E2.2 with the source-CFL guard, then the stationary droplet at
   three $M_\mu$; E1.4 to E1.7 in parallel as dictionary-only runs.
3. **The gate repairs** ([[concepts/method-gates]]): the candidates' closed forms in the 1D arm; a
   checkerboard-sensitive band diagnostic (L2 of the second difference of $\psi$, or min/max of a
   one-sided $q$); a curvature column in the kinematic arms; the zero-flow static arm; the
   extension's time level on outer passes ≥ 2 ([[concepts/halo-limited-extension]]).
4. **The Tier-I cases of the dossiers** (uniform translation, rigid rotation, planar strain with
   $q_0 = 0.92 / 1.00 / 1.08$, perturbed scaling) as gate arms before any coupled run.
5. **Sequencing of the routes**: the velocity-free linear law first; $F$ at the algebraic foot
   (E2.4) if the O(h) shift remains; S2 only after, with its $a_h/\sigma_h$ term off first; ALG only
   for droplets; V1 last.
6. **Shared decisions with the SL line**: a box-length token for the translating arm, the
   oscillating horizon, and the re-run list after the parallel fixes ([[sessions/sl-session-handover]]).

## Related

[[hubs/gradient-control]], [[concepts/gradient-control-overview]], [[sessions/current]].

## Log

### 2026-09-28
Created.
