---
title: "The halo-limited directional extension (HL0)"
description: "The design (capped travel, fractionReached weight, level-set direction, stencilFit sampler, flux-correction form, R = h), the unit gates that pass, and the gate that fails: order 1.39 kinematic, destroyed coupled arms"
kind: concept
status: retracted
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
code: [src/leiaLevelSet/velocityExtension/haloLimited/haloLimited.C, src/leiaLevelSet/velocityExtension/haloLimited/haloLimitedStrategies.C, applications/test/leiaTestHaloLimited]
sources: [STATUS 11.10, STATUS 11.15, gcls article sec:hl and sec:disc-hl, technical report sec 4]
---
# The halo-limited directional extension (HL0)

> Retracted at $R = h$ (2026-09-28). The extension samples the velocity at a point $Y$ moved
> from the face towards the interface by a capped travel $S_R(d)$ along the level-set direction
> and blends it with a weight $w = (S/d)^\beta$ (`fractionReached`): $\mathbf{u}_H = (1-w)\mathbf{u} + w\,\mathbf{u}_h(Y)$,
> applied as a flux correction on internal faces with one quadratic model per cell (`stencilFit`)
> ([haloLimited.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension/haloLimited/haloLimited.C#L269-L357)).
> The unit gates pass to $2\times10^{-16}$ (affine and divergence-free quadratic velocities), and
> the first rung above them, the kinematic shear flow, fails at order 1.39 with 14 times the shape
> error ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)). The report's reading: the extension relocates
> the normal strain into a shell at $d = 1$ to $2h$ ([[concepts/extension-strain-relocation]]).

## What it is

Strategies (dictionary words of [haloLimitedStrategies.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension/haloLimited/haloLimitedStrategies.C)):
travel `capped` ($S = d\,(1 + (d/R)^m)^{-1/m}$, $m = 2$), weight `fractionReached` ($\beta = 1$),
direction `levelSet` (the fitted normal with the Q2 regulariser), sampler `stencilFit` (the CPC
quadratic WLS of each velocity component, evaluated at $Y$ with the cell's own model; no trust
region). `radiusCells` must be at most 1 (a `FatalIOError` otherwise), because the sampler
evaluates the face's own cells' models
([haloLimited.C](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/velocityExtension/haloLimited/haloLimited.C#L108-L114)). The
SL solvers trace with `fvc::reconstruct` of the corrected flux (`SL_TRACE_FLUX extension`); the
kinematic `cellCentred` branch bypasses the extension.

## Evidence

| claim | number | where |
|---|---|---|
| unit gates (D4) | zero normal strain to 1e-12 on affine U; uniform U bit for bit | [STATUS 11.10](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3437-L3471) MEASURED |
| shear arm, shape error at $N = 136$, order | 8.49e-3 against 6.19e-4 (14x); 1.39 against 3.07 | [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023) MEASURED |
| shear arm, band gradient ratio (target 0.8) | 0.828 | same, MEASURED |
| exact-1D `qBandMean` at $N = 32$ to $256$ (baseline $e^{-1} = 0.368$; FP0 0.96 to 1.09) | 0.589, 0.499, 0.510, 0.459 | archive `gate/summaries/HL0/summary.csv`, MEASURED |
| oscillating droplet at $N = 200$: band gradient error / volume error | 2.0e7 / 6.7 (baseline 5.1 / 3.0e-3) | STATUS 11.15 MEASURED |
| translating droplet at $N = 200$: shape error | 0.196 against 3.06e-3 (64x) | same, MEASURED |
| stationary droplet: only the volume error regresses | 2.0e-7 against 2.2e-8 | same, MEASURED |
| coupled seam check | 4.7e-6 PASS | same, MEASURED |
| the extension geometry on outer passes 2 and 3 comes from $\psi^{n+1,(k-1)}$ while $\psi^n$ is transported | `slVelExt->correct()` before the $\psi^n$ restore | [slAlphaEqn.H](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L65-L75) MEASURED (code) |

## Why it failed, or why we think so

[[concepts/extension-strain-relocation]]: the transmitted normal strain $K(t)$, $t = d/R$, is 0 on
the interface and 1.27 at $1.5R$ (DERIVED), so with $R = h$ the shell sits inside the band and the
stencils; the 1D arm measures the band mean predicted by the profile; the O(1) variation of $q$
over one cell makes the transport first order whatever the trace velocity. The pre-print's
`reconstruct` reading is one O($h\lvert\nabla\mathbf{u}\rvert$) contribution, neither necessary
(FP0) nor sufficient (HL0's early volume error). In the translating arm the sampler injects the
spurious current unprojected; in the oscillating arm it fits across the boundary-layer kink
(HYPOTHESIS; E1.6 decides).

## Decisions

Abandoned as a gradient-control device at $R = h$; a mesh-resolved variant (the ALG route,
$R \ge 3h$, the dossier's taper) only after the source repair and only for droplets.

## Related

[[hubs/gradient-control]], [[models/velocity-extension]], [[concepts/closest-point-extension]],
[[concepts/extension-strain-relocation]], [[concepts/why-the-candidates-failed]],
[[concepts/trace-velocity-projected-flux]].

## Log

### 2026-09-28
Created.
