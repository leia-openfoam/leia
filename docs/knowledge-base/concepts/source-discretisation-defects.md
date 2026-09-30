---
title: "The source discretisation: the centred-q eigenmode, the band edge, the soft-wall slope"
description: "Three defects of the SL source step's discretisation, derived on 2026-09-28: a linear instability at rate mu with a centred gradient of q, an O(1) jump of |psi| at the hard band edge, and a nearly discontinuous per-cell map for the soft wall; the clamp is never reached"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/semiLagrangian/source/slGradientControlSource.C]
sources: [technical report sec 3 and appendix B]
---
# The source discretisation: the centred-q eigenmode, the band edge, the soft-wall slope

> Open (DERIVED 2026-09-28; the discriminators have not run). The source step reads $q$ from
> `fvc::grad(psi, "gradPsiSource")`, which every case template sets to `leastSquares`, a CENTRED
> gradient ([slGradientControlSource.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/semiLagrangian/source/slGradientControlSource.C#L146);
> [fvSchemes.template](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSchemes.template#L112-L113)), switches
> the source off with a hard test at $\lvert\psi_c\rvert = 3h_c q_c$ (line 170), and clamps
> $\lvert\Delta t F\rvert$ at 30 (lines 185 to 194). Three consequences follow.

## A. The centred-q eigenmode (rate mu)

Linearise $\psi_c \leftarrow \psi_c\exp(-\mu\Delta t\,(q_c - 1))$ about $\psi = d$:
$\delta_c \leftarrow \delta_c - \mu\Delta t\,d_c\,(\mathbf{n}\cdot\nabla_h\delta)_c$. For the
distance-modulated checkerboard $\delta_c = (-1)^c E\,d_c$ the central difference gives
$(\mathbf{n}\cdot\nabla_h\delta)_c = (-1)^{c+1}E$, so $\delta_c \leftarrow \delta_c\,(1 + \mu\Delta t)$:
an exact eigenmode with growth rate $\mu$, independent of $\Delta t$, $h$ and the flow (in 2D
$(-1)^{i+j}E\,d$). A one-sided (Godunov, Rouy–Tourin) $q$ damps the same mode at $2\mu d/h$. A pure
checkerboard without the $d$ modulation is invisible and neutral, which is why the centred band
metric of the gate cannot see the growing mode. Why the stationary droplet exposes it: at
$u \approx 0$ the SL step is the identity (`evaluateRaw` at the cell centre returns $\psi_c$), so
nothing damps the mode; at CFL 0.5 the quadratic fit halves a checkerboard per step, which hides
it in the shear arm. Why SDPLS never saw it: `R` reads $a$, not $q$, and `beta` is linearQ at
$\mu = 1$ 1/s. Evidence: [[concepts/why-the-candidates-failed]]. Discriminators: E1.1 (zero-flow
static test at three $\mu$), E1.4 ($M_\mu$ = 0.1 and 10), E2.1 (the upwind $q$).

## B. The band edge under an accumulated O(1) factor

A law that reaches its target multiplies $\lvert\psi\rvert$ inside the band by
$e^{\int F\,dt} = q_\mathrm{far}/q_\mathrm{target} \approx 2$ to $3$ in the shear arm and by 1
outside, so $\lvert\psi\rvert$ becomes non-monotone at $3h$ (about $7.5h$ inside against $3.5h$
outside), independent of $h$, and spurious crossings follow (MEASURED: HL1z band error 10.5 and
volume error 1.6 at $T$; HL1q band error 0.44 at $t = 1.9$ s and 6.58 at $T$). With a centred
gradient the fixed point "$q_c = 1$ in the band with outside data $c\,d$" splits into two
sublattices; for $c = 1.5$ the odd sublattice's fixed point at the interface cells has the wrong
sign ($\psi(+0.5h) \to -1.25h$), which the sign-preserving update turns into $\lvert\psi\rvert \to 0$
there: the zero set moves by O(0.5h) with $u = 0$. The z-law adds an unbounded $-\mu q^2/2$ for
the large $q$ this produces (its runaway at $t \approx 0.7$ to $1.0$ s in the shear arm).
Discriminators: E1.2 ($\psi_0 = 1.5d$), E1.3 (no band), E2.2 (a smooth taper).

## C. The soft wall's slope, and the clamp

The slope $\lvert dG/dq\rvert$ of the soft wall at $C_\kappa = 1.25$, $\delta_s = 0.08$, $p = 5$,
$\gamma = \mathrm{artanh}\,0.9$ reaches about $40\sigma$ at $\lvert q - 1\rvert \approx 0.07$
([gradientControlLaws.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/gradientControl/gradientControlLaws.C#L210-L220)), so
$\Delta t\,dG/dq \approx 0.5$ per step in the shear arm ($\sigma \approx 3$ 1/s, $\Delta t = 4.1\times10^{-3}$ s):
a nearly discontinuous per-cell map. Its signature is the 76 % seam failure of S1, which HL1q and
HL1z do not share (4.9e-7 and 5.3e-7). The clamp at $\lvert\Delta t F\rvert = 30$ is never reached
(HL1q at $q = 190$ has $\Delta t F = -0.2$); for linearQ $q \to 0$ gives $F \to \mu$, bounded, so
"the clamp acts where $q$ is small" does not hold. Discriminator: E1.7 ($\delta_s = 0.3$, $p = 1$).

## Decisions

Repair before any law is judged: a one-sided $q$ on unstructured meshes with the coupled-patch
neighbours, a smooth band taper or no band, a source-CFL guard $\mu\,3h\,\Delta t/h < 1/2$, verified
by E1.1 to E1.3 and the non-vacuous 1D closed forms (E0.1).

## Related

[[hubs/gradient-control]], [[models/sl-source]], [[concepts/sl-source-step]],
[[concepts/why-the-candidates-failed]], [[concepts/gradient-control-next-experiments]].

## Log

### 2026-09-28
Created with the technical report.
