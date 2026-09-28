---
title: "Gradient control: the next experiments, cheapest first, with pre-registered predictions"
description: "E0.1 (closed forms in the 1D arm), E1.1 to E1.7 (dictionary changes), E2.1 to E2.4 (small code changes), E3 (the extended sampler and S2): each with its prediction and the outcome that falsifies it"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [technical report sec 6]
---
# Gradient control: the next experiments, cheapest first, with pre-registered predictions

> Open (2026-09-28): none has run. The rule of the technical report stands: no candidate re-enters
> the coupled gate before the kinematic ladder (E1.1 to E1.3, the dossiers' Tier-I cases) passes
> at the baseline's order. Each row is written before its run; the prediction is the read-out.

| id (cost) | experiment | prediction | what it falsifies |
|---|---|---|---|
| E0.1 (hours) | the candidates' closed forms in the 1D arm: $\dot q = q(F - \alpha K(d/R))$, $\dot d = \alpha d(1 - c^2)$ per band cell | HL0 band mean ≈ 0.49 (measured 0.46 to 0.59); S1 → 0.925 (measured 0.998 → 0.697); HL1q with $\mu = \alpha$: band mean → about 0.5 (measured 0.80 → 0.53) | nothing yet; it makes the 1D gate real ([[cases/exact-1d-stretch]]) |
| E1.1 (dictionary) | zero-flow static test: kinematic solver, $u = 0$, circle SDF, $N = 100$, linearQ, bandCells 3, $\mu \in \{27, 270, 2700\}$ 1/s, $\mu\Delta t \le 10^{-2}$ | $\max\lvert\psi_c - d_c\rvert$ grows as $e^{\mu t}$ with the pattern $(-1)^{i+j}d$; the fit curvature error grows at $\mu$; the centred metric stays O(seed) until the amplitude is O(h); rates 1:10:100 | mechanism A; if nothing grows at $u = 0$ the source map is stable and the flow-loop reading returns |
| E1.2 | E1.1 with $\psi_0 = 1.5d$ (dossier B1, Tier-I test 4) | $\lvert\psi\rvert \to 0$ in the interface cells on one sublattice within about $5/\mu$; zero-set displacement O(0.5h); with $c = 1.08$: 0.1 to 0.2h | mechanism B, if the zero set stays put to O($h^2$) |
| E1.3 | E1.1 and E1.2 with bandCells 1000 | B disappears, A remains | separates A from B |
| E1.4 | HL1q stationary at `GC_M_MU` 0.1 and 10 | survives 0.1 s / dies by about 3.5 ms; rate ∝ $\mu$ | A as the primary mechanism |
| E1.5 | HL0 translating with `psiOuterCorrectors no` | the departure from the baseline moves later or vanishes if the time-level mismatch (D) matters | D |
| E1.6 | HL0 translating and oscillating with the extension sampling `reconstruct(phi)` instead of `U` (`velocityExtension.C` lines 49 to 50) | translating improves, oscillating does not; then density and viscosity ratio 1 for the oscillating arm | the unprojected-injection reading against the jump reading |
| E1.7 | S1 with $\delta_s = 0.3$, $p = 1$ (slope about $6\sigma$) | seam ≤ 1e-5, shear shape error ≥ 10x lower, band-error gain shrinks | the bang-bang reading, if the seam still fails |
| E2.1 (days) | a one-sided (Rouy–Tourin) $q$ in `slGradientControlSource::apply` with the coupled-patch neighbours, plus the source-CFL guard $\mu\,3h\,\Delta t/h < 1/2$ | E1.1's mode damped at $2\mu d/h$; HL1q stationary at the baseline level for 0.1 s | A |
| E2.2 | a smooth band taper $F\,\tau(\lvert\psi\rvert/(3hq))$ with a $C^2$ cutoff, or no band | HL1q and S1 shear shape orders rise from about 0 to ≥ 1 | B |
| E2.3 | the kinematic `cellCentred` trace with $\mathbf{u}_\mathrm{ext}$ (`leiaSemiLagrangeLevelSetFoam.C` line 218) | the HL0 shear constant drops, the order stays below 2 (the shell remains); order 3 would mean the trace alone was the cause | the relocation reading |
| E2.4 | $F$ evaluated at the algebraic foot $x_c - d\,\mathbf{n}$ and copied along the normal | the O(h) crossing shift vanishes at first order; A and B remain unless E2.1 and E2.2 are applied | the crossing-shift contribution |
| E3 (weeks) | the extended sampler with `findCell` and a 3-layer halo for $R = 3h$; S2 (scalar memory) | only if HL survives E2.3 | — |

Then the first coupled gate of the next campaign: the stationary droplet with the repaired
source at $M_\mu \in \{0.1, 1, 10\}$, prediction "no growth at rate $\mu$"; a slower residual growth
is the force-side seed ([[concepts/parasitic-current-mechanism]]) and no source fixes it.

## Related

[[hubs/gradient-control]], [[concepts/why-the-candidates-failed]],
[[concepts/source-discretisation-defects]], [[concepts/extension-strain-relocation]],
[[concepts/gradient-control-open-decisions]], [[concepts/method-gates]].

## Log

### 2026-09-28
Created with the technical report.
