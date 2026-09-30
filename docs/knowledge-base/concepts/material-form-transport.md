---
title: "Material-form transport and the Sp(div u) correction"
description: "What the material form of the source (R + div u, the Sp(div Phi) correction) taught: Rdiv diverges kinematically and is withdrawn; a normal-constant extension shears the band tangentially"
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-09-28
sources: [SDPLS article sec Rdiv, technical report sec 4]
---
# Material-form transport and the Sp(div u) correction

> Open (2026-09-28). Two facts stand. (1) The material form of the SDPLS source, `Rdiv` with the
> $\nabla\cdot\mathbf{u}$ term (the `Sp(div phi)` correction), diverges kinematically (orders −3.8,
> −1.3, −1.5) and is withdrawn ([SDPLS article](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L1263), MEASURED).
> (2) Any velocity that is constant along the normals is not a rigid motion of the band: level
> sets at distance $d$ from a circle of radius $R$ rotate at $\Omega(d) = \omega R/(R+d)$ under a
> closest-point extension and at $\omega(R + d(1 - c^2))/(R + d)$ under the halo-limited one, so a
> circle is invariant and a slotted disc is not (DERIVED; the Tier-I rigid-rotation test decides it).

## Why it matters

A source or an extension that is exact for the interface still transports the band with its own
kinematics. The band is what the quadratic fits read, so its distortion enters the curvature and
the phase indicator before it enters the zero set.

## Open questions

The Tier-I rigid-rotation and planar-strain tests with $q_0 \ne 1$ (never run); the non-solenoidal
part of the extension flux, which the `pointValue` scheme ignores on the zero set but the
`reconstruct` averages into the trace velocity.

## Related

[[hubs/gradient-control]], [[models/sdpls-source]], [[concepts/extension-strain-relocation]],
[[concepts/halo-limited-extension]].

## Log

### 2026-09-28
Created.
