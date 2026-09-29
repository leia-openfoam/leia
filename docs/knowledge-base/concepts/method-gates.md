---
title: "The method gates (2D and 3D)"
description: "One 2D and one 3D study test every level-set method as a set of case tokens against the production baseline on five arms with a four-criterion verdict; the first campaign failed every candidate; two blind spots were found on 2026-09-28 and their repairs are open"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-28
date_settled: 2026-09-26
decided_by: [config/gates/methodGate2D.yaml, config/gates/methodGate3D.yaml]
code: [config/gates/methodGate2D.yaml, config/candidates, workflow/scripts/render_gate_configs.py, workflow/Snakefile.gate, workflow/scripts/make_gate_summary.py, workflow/scripts/richardson.py, workflow/scripts/make_gate_tables.py]
sources: [CLAUDE method gates section, G2, workflow/README method gates, PHL 5, STATUS 11.5, STATUS 11.11-11.15, STATUS 11.17, gcls article sec:gate]
---
# The method gates (2D and 3D)

> Verdict (2026-09-28). Any change to the level-set method runs the 2D method gate before any other coupled study, and the 3D gate only after a 2D pass ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L603-L643), [[decisions/process-gates-2d-first-and-no-best-yaml]]). A method enters only as a candidate, a set of case tokens in `config/candidates/<name>.yaml` with its pre-registered read-out; the baseline takes its settings from the `.parameter` layering, never from a copy ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L1-L28)). The 2D gate runs five arms on 4 ranks: the exact 1D stretch first, then the reversed 2D shear, the stationary, the translating and the oscillating droplet, each on a three-rung Richardson ladder, plus seam checks on the coarsest shear and translating rungs ([README](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L192-L235)). The first campaign (2026-09-26 to 27) failed all six candidates; the kinematic shear arm alone failed all of them ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L4002-L4023), [[studies/method-gate-2d-campaign-2026-09]]). Two blind spots were found on 2026-09-28: the exact-1D closed form is `None` for every candidate, so criterion 3 is vacuous there, and the centred band-gradient metric cannot see the cell-scale mode that destroyed the stationary droplet under HL1q and HL1z ([STATUS 11.17](https://github.com/leia-openfoam/leia/blob/feature/gradient-controlled-level-set/STATUS.md#1117-gate-blind-spots-found-on-2026-09-28-open-repairs), working tree). Both repairs are open.

## What it is

Files ([README](https://github.com/leia-openfoam/leia/blob/8867581/workflow/README.md#L196-L204)): the gate definition holds the arms, the solver per line, the line tokens a candidate cannot set, the method tokens a candidate may set, and the verdict thresholds; `render_gate_configs.py` turns gate plus candidate into one study config per arm with explicit refusals; `Snakefile.gate` runs the arms through the ordinary workflow, exact1D first ([`Snakefile.gate`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/Snakefile.gate#L109-L122)); `make_gate_summary.py` writes `summary.csv`, `orders.csv`, `seam.csv`, `vsBaseline.csv` and `verdict.txt`; `richardson.py` computes the orders ([[concepts/richardson-ladders-and-orders]]).

The arms of `methodGate2D` ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L86-L166)):

| arm | case | ladder | horizon | what it tests |
|---|---|---|---|---|
| exact1D | `1Dstretch`, `u = alpha x`, closed form `q = exp(-alpha t)` | N 32, 64, 128, 256 | 1 s | the only true error; runs first, the others start after it passes ([[cases/exact-1d-stretch]]) |
| shear | `2Dvortex`, reversed, R 0.15 | N 68, 96, 136 (R/h 10.2 to 20.4) | 2 s | transport with no force; gradient at T/2, shape at T; seam at np 1 and np 8 against np 4 |
| stationary | `stationaryDroplet2D`, R 1 mm, L 10 mm | N 100, 142, 200 | 0.1 s | the parasitic current, the pressure jump, the curvature |
| translating | `translatingDroplet2D`, U0 0.05 m/s, start 2.5 mm | N 100, 142, 200 | 0.05 s | the Galilean test; the coupled seam check `translatingSeamNp1` at tol 1e-5 |
| oscillating | `oscillatingDroplet2D`, signed-distance ellipse a = 1.1 R | N 100, 142, 200 | 0.1 s | period and damping of mode 2 against Lamb; the band gradient drift |

Every droplet arm merges the `twoPhaseCoupling` block, the explicit copy of the measured best coupling ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L61-L84), [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3578-L3595)); `CURVATURE_EXTENSION` is set per arm. The droplet time step is `dt = 0.010861 N^-1.5`, 0.2323 of the Brackbill limit, at every rung ([[concepts/capillary-time-step]]). The 3D gate (`methodGate3D`, np 32) runs the same 1D arm, the 3D shear at N 68, 90, 118, a seam at np 16, and the three droplets in the 6R box at N 60, 78, 102 ([PHL 5.3](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-halo-limited-gradient-control.md#L521-L537)).

The verdict, pre-registered ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L20-L28), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L382-L447)):

1. Completion: a candidate case that DIVERGED where the baseline COMPLETED is a FAIL, and a result. A baseline rung that did not complete makes the run INVALID for that rung instead of a silent pass (added after the translating baseline diverged at every 0.1 s rung).
2. No regression: every vector metric at the finest rung at most 1.10 times the baseline, every least-squares order at least the baseline order minus 0.3.
3. The candidate's own target, from its file (for the first campaign: the shear band gradient error at the finest rung at most 0.8 times the baseline).
4. No effect: a candidate whose every CSV equals the baseline's is a FAIL; its tokens were not consumed. MEASURED 2026-09-26: `VELOCITY_EXTENSION closestPoint` on the SL line was such a no-op, because the `projectedFlux` trace ignores the extension ([STATUS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3281-L3306)).
5. The seam checks pass: the shear rung at 1e-10, the translating rung at 1e-5 on nine columns; a reference that did not complete is NOT_COMPARABLE ([[concepts/seam-checks-and-decomposition-invariance]]).

The vector is L2 and L1 only, read at the gate's instants ([[concepts/error-vector-and-read-out-instants]]). `rhoClipFraction` is reported, not scored (4bff922); the oscillating arm's velocity is the physical oscillation and is reported as `oscL2MagU`, not scored (4ef98db) ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3867-L3878)).

## Why it matters

A method has given good advection and bad hydrodynamics here before, and only the coupled arms show that ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L603-L616)). The gate is the executable form of the research loop: cheapest discriminator first, one variable plus the baseline on the same commit and binaries, the whole vector, the order next to the error.

## Evidence

| claim | number | where |
|---|---|---|
| the baseline smoke of the 2D gate | 14 cases COMPLETED at np 4; exact1D q error 1.6e-6; seam np 1 and np 8 against np 4 at most 1.9e-12 | [STATUS 11.5](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3281-L3306), MEASURED |
| the first campaign, pre-fix, on the shear arm alone | shape error at N = 136 against the baseline's 6.19e-4: HL0 14x (order 1.39 against 3.07), S1 410x (0.06), HL1q 170x (0.02), HL1z 2700x | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3747-L3763), MEASURED |
| the fixed gate's verdicts (target ratio, regressions, orders, seam, completion) | FP0 FAIL 0.939 / 3 / 8 / 5 / 5; HL0 FAIL 0.828 / 16 / 13 / 0 / 0; HL1q FAIL 0.204 / 22 / 18 / 0 / 0; HL1z FAIL 6.490 / 23 / 19 / 0 / 0; HL2 FAIL 0.754 / 7 / 5 / 0 / 3; S1 FAIL 0.707 / 12 / 5 / 1 / 1 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L4002-L4023), MEASURED |
| the coupled seam check | baseline 3.5e-7, HL0 4.7e-6, HL1q 4.9e-7, HL1z 5.3e-7 PASS; S1 0.76 and FP0 1.11 FAIL; HL2 NOT_COMPARABLE | same, MEASURED |
| the kinematic arms are unchanged by the parallel fixes | 177 CSV pairs identical at tolerance 0 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3879-L3886), MEASURED |
| the translating horizon | the 0.1 s runs diverged in every candidate (baseline at 0.0904 / 0.0942 / 0.0772 s); the arm ends at 0.05 s, a change stated after seeing data | [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3596-L3623), MEASURED |
| the baseline's oscillating arm at the finest rung | band gradient error 5.13 at N = 200 against 0.10 and 0.14; period 8.06 ms against 10.0 and 9.82 ms; Celik: divergent | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3955-L3971), MEASURED |
| the 1D closed form is vacuous for candidates | `q_exact = None` unless `VELOCITY_EXTENSION none`, `SL_SOURCE none` and `SDPLS_SOURCE noSource`; the check then passes with "oracle pending" | [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L190-L198), [same](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L478-L490), MEASURED (reading of the script) |
| the centred metric is blind to the cell-scale mode | at 4.3 ms the HL1q band gradient metric is 3.2e-4 against the baseline's 4.4e-4 while the curvature error is 10x; at 10 ms 3.9x against 175x | [[concepts/why-the-candidates-failed]], the technical report `docs/gradient-controlled-level-set/gcls-technical-report/` (uncommitted at 8867581), MEASURED; the invisibility of a pure checkerboard to a centred difference is DERIVED |

## The two blind spots and their repairs (open)

1. The exact-1D closed form. The summary sets `q_exact` to `None` for any candidate that has an extension or a source, so `qError` is `None`, criterion 3 is vacuous in the 1D arm, and the arm can only report completion for candidates. Repair: integrate the closed form of each candidate per band cell, `dq/dt = q (F - alpha K(d/R))` and `dd/dt = alpha d (1 - c^2)`, and compare `qBandMean` with it (experiment E0.1 of the technical report; predicted HL0 band mean about 0.49 against the measured 0.46 to 0.59).
2. The band-gradient metric. `gradPsiMetric leastSquares` is a centred difference, so the mode `(-1)^(i+j) d` does not change it. Repair: an L2 norm of the second difference of psi, or the minimum and maximum of a one-sided `q`, and a curvature column in the kinematic arms.

Both from [STATUS 11.17](https://github.com/leia-openfoam/leia/blob/feature/gradient-controlled-level-set/STATUS.md#1117-gate-blind-spots-found-on-2026-09-28-open-repairs) (working tree, not at 8867581) and [[retractions/gcls-coupled-loop-reading]]. A third repair is named in the same place: a zero-flow static arm ([[concepts/gradient-control-next-experiments]]).

## Decisions

- 2D before 3D, no `config/best.yaml`, METHOD.md in the same commit: [[decisions/process-gates-2d-first-and-no-best-yaml]].
- The explicit `twoPhaseCoupling` block (8a9b85a, user instruction 2026-09-27): the block changes no case file of the stationary and oscillating arms; the translating arm differs only in `endTime` ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3545-L3577)).
- The horizon of the translating arm, 0.05 s: pre-registered as 0.1 s and changed after the baseline diverged; recorded as such ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3617-L3623)).

## Open questions

1. The two repairs above, and the zero-flow static arm.
2. The horizon of the oscillating arm and a baseline check that reads the gradient drift over time ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3955-L3971)).
3. A box-length token for the translating arm ([[cases/translating-droplet]]).
4. The Eulerian line's coupled arms wait for its curvature dispatch and droplet CSV ([[concepts/eulerian-solver-mass-flux-port]]).
5. The 3D gate has not run: no 2D pass, and the 3D translating case had no reference velocity until 2026-09-27 ([[concepts/wrong-setup-voids]]).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/error-vector-and-read-out-instants]], [[concepts/richardson-ladders-and-orders]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/bit-identity-and-inertness-gates]], [[concepts/data-archive-per-version]], [[concepts/wrong-setup-voids]].
- Cases: [[cases/exact-1d-stretch]], [[cases/benchmark-cases]], [[cases/translating-droplet]], [[cases/stationary-droplet]], [[cases/oscillating-droplet]], [[cases/kinematic-advection-cases]].
- Campaign and readings: [[studies/method-gate-2d-campaign-2026-09]], [[concepts/why-the-candidates-failed]], [[studies/gcls-pre-print]], [[retractions/gcls-coupled-loop-reading]].
- Decisions: [[decisions/process-gates-2d-first-and-no-best-yaml]].

## Log

### 2026-09-28
Created.
