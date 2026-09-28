---
title: "Gradient control: source terms and velocity extension"
description: "Hub: keeping |grad psi| near one without a reinitialisation, by source terms (SDPLS, the SL source step) or by a velocity extension; the first campaign (2026-09-26/27) failed on every candidate, and the technical report of 2026-09-28 says why and what to do next"
kind: hub
status: open
part: gradient-control
tags: [hub, part/gradient-control]
date: 2026-09-28
---
# Gradient control: source terms and velocity extension

[[index]] <- back

## The question this part answers

Can a source term in the level-set equation, or an extension of the velocity off the interface, keep $q = |\nabla\psi|$ near one in the narrow band without a reinitialisation, without moving the zero set, and without losing the third-order transport of the baseline?

## Current verdict (2026-09-28)

Not yet. The first campaign tested six candidates against the production method in the fixed 2D method gate ([[studies/method-gate-2d-campaign-2026-09]]): the soft-wall source alone (S1), the halo-limited extension alone (HL0), the extension with a linear law in $q$ (HL1q) or in $z = q^2$ (HL1z), the extension with the soft wall (HL2), and the closest-point extension (FP0). Every candidate FAILS ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)): the kinematic shear arm alone fails all of them (HL0 raises the shape error 14 times at order 1.39 against 3.07; S1 410 times at order 0.06; HL1q 173 times at order 0.02), and HL1q and HL1z destroy the stationary droplet at $N = 200$ (shape error 3.8 R, volume error 76 %). The pre-print of 2026-09-27 ([[studies/gcls-pre-print]]) named the crossing shift of the source step, the reconstruct of the corrected flux and a coupled loop through the velocity; the technical report of 2026-09-28 ([gclsTechnicalReport.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex)) corrects that reading ([[retractions/gcls-coupled-loop-reading]]): the linear laws read no velocity, and the measured growth rates (about 290 and 320 to 350 1/s over 10 to 20 ms) are of the order of $\mu = 1/T_\mathrm{ref} = 270$ 1/s and far below the mesh-scale capillary rate, so the primary mechanism is the source's own discretisation, a centred gradient of $q$ that makes the cell-scale mode $(-1)^{i+j} d$ grow at rate $\mu$ ([[concepts/source-discretisation-defects]], DERIVED, discriminator: the zero-flow static test). The halo-limited extension at $R = h$ relocates the normal strain into a shell at $d = 1$ to $2h$ with a 27 % overshoot; it does not remove it ([[concepts/extension-strain-relocation]], DERIVED and MEASURED on the 1D arm). The gate has two blind spots: the exact-1D closed form is vacuous for candidates, and the centred band metric cannot see the mode ([[concepts/method-gates]]). Earlier lines of the same family: the Eulerian SDPLS source R converges kinematically and fails under coupling through the strain reading ([[concepts/sdpls-source-eulerian]]); the legacy velocity extensions are dominated on pure advection ([[models/velocity-extension]]); the normal-projected SL is closed ([[concepts/normal-projected-sl]]).

## Map

| kind | note | status | one line |
|---|---|---|---|
| concept | [[concepts/gradient-control-overview]] | open | the direction, the dossiers, the plan, the pre-print, the report |
| model | [[models/gradient-control-law]] | settled | the eight laws and the strain weights |
| model | [[models/sl-source]] | candidate | the SL source step (exponential update, band, clamp) |
| model | [[models/sdpls-source]] | settled | noSource, R, beta, Rdiv, RdivStrictSp, gradientControl |
| model | [[models/velocity-extension]] | settled | the eight extensions and their verdicts |
| concept | [[concepts/sdpls-source-eulerian]] | settled | R converges kinematically; coupled it amplifies the mode-4 current 260 times |
| concept | [[concepts/sl-source-step]] | candidate | the update, the zero-set shift, the gate results |
| concept | [[concepts/combined-source-note]] | retracted | subsumed by gradientControl linearQ with strain weight full |
| concept | [[concepts/halo-limited-extension]] | retracted | the design, the unit gates, HL0's failure |
| concept | [[concepts/closest-point-extension]] | retracted | FP0: seam-dependent, same order loss |
| concept | [[concepts/material-form-transport]] | open | the Sp(div u) correction and what Rdiv taught |
| concept | [[concepts/why-the-candidates-failed]] | open | the mechanisms, their marks, the discriminators |
| concept | [[concepts/source-discretisation-defects]] | open | the centred-q eigenmode, the band edge, the soft-wall slope |
| concept | [[concepts/extension-strain-relocation]] | open | the K(t) profile and the shell |
| concept | [[concepts/gradient-control-next-experiments]] | open | E0.1, E1.1 to E1.7, E2.1 to E2.4, E3 |
| concept | [[concepts/gradient-control-open-decisions]] | open | what the supervisor decides |
| retraction | [[retractions/gcls-coupled-loop-reading]] | retracted | the pre-print's loop reading for HL1q/HL1z |
| case | [[cases/exact-1d-stretch]] | settled | the 1D arm and its vacuous closed form |
| study | [[studies/gcls-pre-print]], [[studies/sdpls-pre-print]], [[studies/velocity-extension-pre-print]], [[studies/method-gate-2d-campaign-2026-09]] | open | the documents and the campaign |
| session | [[sessions/2026-09-27-gcls-first-campaign]] | settled | what the session built |

## Open, in order

1. E0.1: the candidates' closed forms in the exact-1D arm, so the gate's target criterion is real ([[concepts/gradient-control-next-experiments]]).
2. E1.1 to E1.3: the zero-flow static test of the source at $\mu \in \{27, 270, 2700\}$ 1/s, with $\psi_0 = d$ and $1.5 d$, band 3h and no band: decides the eigenmode and the band-edge mechanism.
3. E2.1 and E2.2: a one-sided (Rouy-Tourin) $q$ and a smooth band taper in the source step, with a source-CFL guard; then the stationary droplet at $M_\mu \in \{0.1, 1, 10\}$.
4. The velocity-free linear law as the only law family with a clean linearisation; $F$ at the algebraic foot (E2.4) if the O(h) shift remains.
5. Abandoned as formulated: HL at $R = h$, HL2, the soft wall as parameterised, FP0 as production, the strain-weight (Taylor/weighted-SDPLS) route.

## Log

### 2026-09-28
Created with the technical report.
