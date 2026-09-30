---
title: "The SL source step: the exponential update and the discrete zero set"
description: "psi <- psi exp(dt F) after the characteristic step: invariant zero set in the continuum, a moved zero crossing on the mesh (an O(h) effect), and the gate results of S1, HL1q, HL1z, HL2"
kind: concept
status: candidate
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [gcls article sec:source-step and sec:disc-source, technical report sec 3]
---
# The SL source step: the exponential update and the discrete zero set

> Candidate (2026-09-26); every gate arm with it FAILS. The update $\psi_c \leftarrow \psi_c\exp(\Delta t F_c)$
> keeps every cell's sign. Two neighbours across the interface at $\psi_A = -a$ and $\psi_B = b$
> hold the linear crossing at the fraction $a/(a+b)$ from $A$; after one step with $F_A \ne F_B$
> it moves by $\delta = h\,ab/(a+b)^2\,\Delta t\,(F_A - F_B) + O(\Delta t^2)$
> ([gcls article, discussion](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradient-controlled-level-set/gcls-level-set-article/gclsLevelSet.tex#L853-L880), DERIVED).
> This shift is first order in $h$: summed over the shear run it is at most about $0.2h$ for the
> linear laws, while the measured shape errors are O(1) in $h$ (orders 0.02 to 0.06). The technical
> report therefore demotes it to a secondary contribution; the primary mechanism is the
> discretisation of the step ([[concepts/source-discretisation-defects]]).

## What it is

See [[models/sl-source]] for the implementation: the band $\lvert\psi_c\rvert \le 3 h_c q_c$ (a hard
switch), the centred $q$, the clamp $\lvert\Delta t F\rvert \le 30$, once per outer pass on the
restored $\psi^n$.

## Evidence

| claim | number | where |
|---|---|---|
| S1: identical to the baseline where $\sigma = 0$ (stationary droplet, four digits); shear shape error | 0.254 against 6.19e-4 (410x), order 0.06 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023) MEASURED |
| S1: oscillating volume error | 21.3 against 3.0e-3 | archive verdict, MEASURED |
| S1: coupled seam check | 0.76 FAIL (baseline 3.5e-7) | same, MEASURED |
| HL1q: shear band gradient ratio / shape ratio | 0.204 / 173, shape order 0.02 | same, MEASURED |
| HL1z: shear band gradient error / shape error at $N = 136$ | 10.5 / 1.65; volume error 1.6 | same, MEASURED |
| HL2: translating arm | DIVERGED at all three rungs; oscillating volume error 14.6 | same, MEASURED |
| the summed crossing shift for the linear laws in the shear arm | $\le (h/4)\,\mu T\,\max\lvert\Delta q\rvert \approx 0.2h$ | technical report, DERIVED |

## Why it failed, or why we think so

Three mechanisms, in the order the report holds them responsible: (A) the centred-$q$ eigenmode
at rate $\mu$ (DERIVED; E1.1 decides), (B) the band edge under the accumulated factor
$e^{\int F\,dt} = q_\mathrm{far}/q_\mathrm{target} \approx 2$ to $3$ (DERIVED; E1.2, E1.3 decide),
(C) the crossing shift, O(h) (DERIVED; E2.4 removes it). For the soft wall the 40σ slope makes the
map bang-bang; its 76 % seam failure is the signature ([[concepts/source-discretisation-defects]]).

## Decisions

Repair the discretisation before any law is judged ([[concepts/gradient-control-next-experiments]]).

## Related

[[hubs/gradient-control]], [[models/sl-source]], [[models/gradient-control-law]],
[[concepts/why-the-candidates-failed]], [[concepts/halo-limited-extension]].

## Log

### 2026-09-28
Created.
