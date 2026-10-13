---
title: "The centred-q mode (-1)^(i+j) d grows at the rate mu"
description: "RETRACTED 2026-10-13 - at zero velocity the source step with the centred q shows no growth at e^(mu t): a bounded rise, a plateau, then a slow growth at 0.067 mu that does not depend on mu dt; the centred q is still the defect, through the reversed group velocity of the sawtooth mode; the coupled rates near mu are not explained by it"
aliases: [centred-q eigenmode, mechanism A growth law, eigenmode at rate mu]
kind: retraction
status: retracted
part: gradient-control
tags: [retraction, part/gradient-control]
date: 2026-10-13
code: [src/leiaLevelSet/semiLagrangian/source/slGradientControlSource.C, workflow/scripts/source_step_replica.py, workflow/scripts/source_checkerboard_diag.py]
sources: ["STATUS 11.26", "technical report sec 3 (sec:eigenmode) and appendix B", config/gates/sourceZeroFlow2D.yaml]
---
# The centred-q mode (-1)^(i+j) d grows at the rate mu

> RETRACTED 2026-10-13. The claim was: "for the distance-modulated checkerboard $(-1)^{i+j}E\,d$
> the source step $\psi_c \leftarrow \psi_c \exp(-\mu\Delta t(q_c - 1))$ with the centred $q$ is an
> exact eigenmode with growth rate $\mu$, independent of $\Delta t$, $h$ and the flow"
> ([[concepts/source-discretisation-defects]], section A; the technical report, `sec:eigenmode`).
> The measurement (gate E1, STATUS 11.26, $u = 0$, $N = 100$, band 3 cells): the sawtooth amplitude
> rises from 5.2e-8 to 2e-5 by $\mu t = 4$, stays there to $\mu t \approx 16$, then grows at
> 0.067 $\mu$; the three $\mu$ give the same numbers in step space. A seeded mode of amplitude 1e-6
> does not grow at $e^{\mu t}$ either (the numpy copy of the step, which reproduces E1 to 2.6e-8).
> Scope: the growth law, the predictions built on it (E1.1, the rates 1:10:100, the rate $\mu$ of
> E1.4), and the report's use of it to explain the coupled rates near $\mu$. The defect (the
> centred $q$) stands.

## The claim, and where it lived

- [[concepts/source-discretisation-defects]], section A ("an exact eigenmode with growth rate $\mu$").
- The technical report, `docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex`,
  `sec:eigenmode` and appendix B, and its summary table (HL1q row).
- [[concepts/why-the-candidates-failed]], the first row of the readings table and evidence item 1:
  the measured coupled rates (about 290 and 320 to 350 1/s against $\mu = 270$ 1/s) were read as this
  eigenmode.
- [[concepts/gradient-control-next-experiments]], rows E1.1 and E2.1.

## Why it was wrong

The linearisation is right, and the mode it names is not bounded: $(-1)^{i+j}E\,d$ grows with
$\lvert d\rvert$, so in a band of finite width it is not an eigenmode. Write the perturbation as a
sawtooth with an envelope, $\delta_c = (-1)^{i+j}A(\mathbf{x}_c)$, with $A$ smooth. The centred
difference of the sawtooth part is zero, and of the product it is $-(-1)^{i+j}\nabla_h A$. The
linearised step then reads
$$A_t = \mu\, d\, \mathbf{n}\cdot\nabla A,$$
a transport of the envelope TOWARD the zero set with speed $\mu\lvert d\rvert$, with $A$ constant
along its paths. Smooth errors move the other way: the characteristics of
$\psi_t = -\mu\psi(\lvert\nabla\psi\rvert - 1)$ leave the zero set with speed $\mu\lvert d\rvert$.
This is the reversed group velocity of the sawtooth mode under a centred difference (Trefethen,
"Group velocity in finite difference schemes", SIAM Review 24, 1982). The growth $e^{\mu t}$ of
$E\,d$ is the arrival of envelope from larger $\lvert d\rvert$; a band cuts that off. What follows,
MEASURED:

1. The rise saturates at the amplitude that the band edge supplies (2e-5 by $\mu t = 4$ at
   $N = 100$, OpenFOAM).
2. A wider band gathers more (numpy copy, $\mu t = 80$, $N = 100$: the median $\psi/d$ of the
   interface cells is 0.18 at band 3 and 0.0000 at band 12).
3. The damage does not depend on $\mu\Delta t$ between 0.001 and 0.01 (numpy copy: 0.068 and
   0.067 $\mu$): it is a property of the spatial operator, not of the time step.

## What survives

- The centred $q$ is the defect of the source step and must not be used: at $u = 0$ it flattens the
  interface cells (median $\psi/d$ 0.18 at $\mu t = 80$, $N = 100$) and moves the zero crossings by
  up to 0.74 h. [[negative-log]] has the line.
- A one-sided $q$ removes it. The compact two-pass `qScheme upwindFit` is the repair under test
  (gate E2.1, STATUS 11.26).
- At $u \approx 0$ the SL step is the identity (MEASURED: the copy without an SL step reproduces E1
  to 2.6e-8 after 4000 steps).
- OPEN, re-opened by this retraction: the cause of the coupled rates near $\mu$ (HL1q at
  $N = 200$). The source's own map does not grow at $\mu$, so the flow, the velocity extension or
  the SL step with a small velocity must carry it. The pre-print's coupled-loop reading
  ([[retractions/gcls-coupled-loop-reading]]) is not reinstated by this; E1.4 (the scaling with
  $\mu$) still decides.

## Propagation (checklist, same commit)

- [x] STATUS.md 11.26
- [x] METHOD.md: not affected (`SL_SOURCE none` is the default; the source is not in production)
- [x] the technical report: a retraction paragraph at `sec:eigenmode`
- [x] [[concepts/source-discretisation-defects]], [[concepts/why-the-candidates-failed]],
  [[concepts/gradient-control-next-experiments]], [[hubs/gradient-control]],
  [[retractions/gcls-coupled-loop-reading]]
- [x] the line in [[retraction-log]]

## Related

[[concepts/source-discretisation-defects]] - [[concepts/why-the-candidates-failed]] -
[[concepts/gradient-control-next-experiments]] - [[hubs/gradient-control]] - [[models/sl-source]]

## Log

### 2026-10-13
Retracted on gate E1 and the numpy copy of the source step. Entered in [[retraction-log#2026-10]].
