---
title: "gradientControlLaw: the eight laws and the strain weights"
description: "The runtime-selectable law family F = w(q^2) a + G(q^2, sigma) that both source models read; the linear laws read no velocity, the soft wall reads the strain magnitude"
kind: model
status: settled
part: gradient-control
tags: [model, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/gradientControl/gradientControlLaw.H, src/leiaLevelSet/gradientControl/gradientControlLaws.C, src/leiaLevelSet/gradientControl/strainWeight.C]
sources: [STATUS 11.9, gcls article sec:laws]
---
# gradientControlLaw: the eight laws and the strain weights

> Settled 2026-09-26 (implemented, 109 unit checks, inert by default). The law family is
> $F = w(q^2)\,a + G(q^2, \sigma)$ with $q = |\nabla\psi|$, $a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n}$ the normal strain rate and $\sigma = |\mathbf{D}|$ the strain magnitude
> ([gradientControlLaw.H](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/gradientControl/gradientControlLaw.H#L170-L185)). Every
> law has $F(q = 1, \sigma = 0) = 0$. With strain weight `none` and a linear law, $F$ reads NO velocity:
> this fact decides the reading of the coupled failures ([[retractions/gcls-coupled-loop-reading]]).

## What it is

The family lives in `libleiaGradientControl` and is read by both source models, the semi-Lagrangian
source step ([[models/sl-source]]) and the Eulerian SDPLS source ([[models/sdpls-source]]). The
dictionary words are the `TypeName` strings of
[gradientControlLaws.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/gradientControl/gradientControlLaws.C) and
[strainWeight.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/gradientControl/strainWeight.C); the tokens are `GC_LAW`,
`GC_M_MU` ($\mu = M_\mu / T_\mathrm{ref}$ per arm), `GC_STRAIN_WEIGHT` and the soft-wall parameters
`SW_*` ([STATUS 11.8](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3371-L3402)).

| member | dictionary word | reads | status | one line |
|---|---|---|---|---|
| no rate | `none` | nothing | settled | $G = 0$; with weight `full` the law is the SDPLS strain reading $F = a$ |
| linear in q | `linearQ` | $q$ | candidate | $G = -\mu (q - 1)$; fixed point with $w = 1$: $q^* = 1 - a/\mu$ |
| linear in z | `linearZ` | $q$ | candidate | linear in $z = q^2$ with the same slope at $q = 1$; the term $-\mu q^2/2$ is unbounded for large $q$ |
| cubic in q | `cubicQ` | $q$ | candidate | untested in a gate |
| cubic in z | `cubicZ` | $q$ | candidate | untested in a gate |
| regularised two-thirds z | `twoThirdsZReg` | $q$ | candidate | untested in a gate |
| saturated linear z | `saturatedLinearZ` | $q$ | candidate | untested in a gate |
| soft wall | `softWall` | $q$, $\sigma$ | retracted as parameterised | $G = -C_\kappa\sigma\,\mathrm{sign}(e)\tanh(\gamma (\lvert e\rvert/\delta_s)^p)$ with $e = q - 1$; at $C_\kappa = 1.25$, $\delta_s = 0.08$, $p = 5$, $\gamma = \mathrm{artanh}\,0.9$ its slope reaches about $40\sigma$ at $\lvert e\rvert \approx 0.069$: a nearly discontinuous per-cell map ([[concepts/source-discretisation-defects]]) |

| strain weight | dictionary word | $w(q^2)$ | one line |
|---|---|---|---|
| none | `none` | 0 | the law reads no strain; the gate's candidates |
| full | `full` | 1 | the SDPLS reading; `sdplsGradientControl` with `none` + `full` equals `R` bit for bit (189 checks) |
| omega | `omega` | a function of $q^2$ | the dossier's weighted variant |

## Why it matters

The law decides what the source can see. A law without $\sigma$ and without $w$ closes the loop
$\psi \to q \to F \to \psi$ with no flow in it; that is why the HL1q/HL1z destruction of the
stationary droplet at a rate equal to $\mu$ is a property of the source discretisation and not of
the capillary coupling ([[concepts/why-the-candidates-failed]]).

## Evidence

| claim | number | where |
|---|---|---|
| unit checks of the laws and weights | 109 pass | [STATUS 11.9](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3403-L3436) MEASURED |
| the soft-wall slope | $1.25\sigma \cdot 2.57/\delta_s \approx 40\sigma$ at $\lvert q-1\rvert \approx 0.069$ | technical report, DERIVED |
| HL1q and HL1z agree to three digits until $t = 0.02$ s on the stationary droplet | same slope at $q = 1$ | archive `gate/histories_stationary.csv`, MEASURED |

## Decisions

No law is production. The candidates of the first campaign used `linearQ`, `linearZ` and
`softWall` with weight `none` ([[studies/method-gate-2d-campaign-2026-09]]).

## Open questions

The cubic, regularised and saturated members have never entered a gate; the technical report says
no law is worth a gate before the source discretisation is repaired.

## Related

[[hubs/gradient-control]], [[models/sl-source]], [[models/sdpls-source]], [[concepts/sl-source-step]],
[[concepts/sdpls-source-eulerian]], [[concepts/source-discretisation-defects]].

## Log

### 2026-09-28
Created.
