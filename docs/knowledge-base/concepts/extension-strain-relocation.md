---
title: "The halo-limited extension relocates the normal strain into a shell"
description: "K(t) = 1 - (1 - t^4)(1 + t^4)^(-3/2) is the normal strain the capped/fractionReached extension transmits at d = tR: zero on the interface, 1.27 at 1.5R; with R = h the shell sits inside the band and the stencils; measured on the 1D arm"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [technical report sec 4 and appendix B]
---
# The halo-limited extension relocates the normal strain into a shell

> Open (DERIVED 2026-09-28; MEASURED on the exact-1D arm). For travel `capped` with $m = 2$ and
> weight `fractionReached` with $\beta = 1$: $S = c\,d$, $w = c$, $c = (1 + t^4)^{-1/4}$, $t = d/R$.
> The extended velocity $\mathbf{u}_H = \mathbf{u} - wS\,\partial_n\mathbf{u}$ transmits the normal strain
> $K(t) = 1 - d(wS)/dd = 1 - (1 - t^4)(1 + t^4)^{-3/2}$ = 0, 0.14, 1.00, 1.27, 1.21, 1.11 at
> $t$ = 0, 0.5, 1, 1.5, 2, 3, and $\int_0^D (1 - K)\,dd = wS(D) \to R^2/D$: the strain is displaced
> by about $1.5h$ into a shell with a 27 % overshoot, not removed. The dossier's own requirement,
> that the transition be broad and outside the geometry-evaluation neighbourhood, is violated by
> $R = h$ (band $3h$, stencil $1.4h$).

## Evidence

| claim | number | where |
|---|---|---|
| exact-1D `qBandMean` at $N = 32, 64, 128, 256$: HL0 against the K-profile prediction (band average of $e^{-K}$ over $d = 0.5h, 1.5h, 2.5h$) | 0.589, 0.499, 0.510, 0.459 against 0.49; baseline $e^{-1} = 0.368$ at every rung; FP0 0.99, 1.09, 1.08, 0.96 | archive `gate/summaries/*/summary.csv`, MEASURED |
| shear arm at $t = 0.02$ s (five steps), band error at $N = 136$ | HL0 0.0200, baseline 0.0354, FP0 0.0005 | archive `gate/histories_shear.csv`, MEASURED |
| the q profile's variation over one cell and the reconstruction error at the feet | $\psi''' \sim 0.6/h^2$; error $\sim 0.1h$ against $\sim 4\times10^{-4}h$ for the baseline | DERIVED |
| HL0/baseline shape-error ratio over the three shear rungs | 4.3, 8.4, 14 (orders 1.39 against 3.07) | archive `summaries/HL0/orders.csv`, MEASURED |
| HL0's excess volume error appears before any q-kink | 2.3e-4 against 5.5e-5 at $t = 0.02$ s | archive histories, MEASURED |
| FP0 has no shell and loses the same order | 1.38 | archive verdict, MEASURED |
| HL0 on the translating droplet: tracks the baseline to 0.02 s (curvature error 77 against 67 1/m), then curvature error and velocity grow together | about 235 and 210 1/s over 30 to 45 ms | archive `histories_translating.csv`, MEASURED |
| the oscillating droplet's air-side boundary layer at $N = 200$ | $\delta \approx \sqrt{\nu T} \approx 8h$, $\partial_n u_t \approx 250$ 1/s (DERIVED); HL0's band gradient error grows at about 480 1/s over 10 to 20 ms (0.29 to 29.5) | archive `histories_oscillating.csv`, MEASURED |
| the sample point can lie outside the stencil hull | $Y_f$ up to $1.5h$ from a cell centre; `evaluateRaw` has no trust region | code, MEASURED |
| the R ladder is blocked | `radiusCells > 1` is a `FatalIOError` | [haloLimited.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension/haloLimited/haloLimited.C#L108-L114) |

## Why it matters

The pre-print's read-out for HL0 ("zero normal strain on the interface, so the band error falls")
conflated the interface with the band. The reconstruct argument of the pre-print is one
O($h\lvert\nabla\mathbf{u}\rvert$) contribution, neither necessary nor sufficient. The translating
droplet has no physical gradient to fit across (rigid translation, $D = 0$) and HL0 fails 64x there:
the sampler reads the cell-scale spurious current $u'$ ($\sim 6\times10^{-4}$ m/s) and injects
$w[\mathbf{u}_h(Y) - \mathbf{u}_h(x)] \sim u'$ into the trace flux UNPROJECTED, while the baseline traces
`reconstruct(phi)`. The ratio-1 test proposed by the pre-print cannot separate the two readings
(it shrinks $u'$ and the kink at once); sampling `reconstruct(phi)` can (E1.6).

## Decisions

Abandon HL at $R = h$ as a gradient-control device. A mesh-resolved extension ($R \ge 3h$, the
dossier's K taper with $\ell_0 = 3h$, $\ell_1 = 6h$, $K_p \approx 1.8$) needs a sampler that locates
the cell containing $Y$ with a 3-layer halo: the ALG route, after the source repair, for droplets
only (it does not help the vortex). Discriminators: E1.5 (psiOuterCorrectors no), E1.6
(reconstruct(phi) sampling), E2.3 (trace with $\mathbf{u}_H(x_c)$), E3 (the extended sampler).

## Related

[[hubs/gradient-control]], [[concepts/halo-limited-extension]], [[concepts/closest-point-extension]],
[[concepts/material-form-transport]], [[concepts/why-the-candidates-failed]],
[[concepts/gradient-control-next-experiments]].

## Log

### 2026-09-28
Created with the technical report.
