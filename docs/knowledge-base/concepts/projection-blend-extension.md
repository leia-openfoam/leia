---
title: "The blended projection velocity extension on the Eulerian line (D. Bothe, use B)"
description: "Open 2026-10-31 - Uext = delta U + (1 - delta) U(x - dist n): on the 2D vortex and the 3D shear it halves the band gradient error at T/2 but makes the shape error at T 1.1 to 5 times and the volume error at T 8 to 11 times larger (B1 FAIL); the 3D deformation runs again after a setup error"
aliases: [projectionBlend, VELOCITY_EXTENSION projectionBlend, use B, blended extension]
kind: concept
status: open
part: gradient-control
tags: [concept, part/gradient-control]
date: 2026-10-29
code: [src/leiaLevelSet/velocityExtension/projectionBlend.H, config/projectionBlend/projectionBlend2Dvortex.yaml, config/projectionBlend/projectionBlend3Dshear.yaml, config/projectionBlend/projectionBlend3Ddeformation.yaml, workflow/scripts/projection_blend_metrics.py]
sources: [STATUS 11.28 item 9, gcls-static-source-report section 6, D. Bothe's notes of 2026-10]
---
# The blended projection velocity extension on the Eulerian line (D. Bothe, use B)

> 2026-10-31 (STATUS 11.28 item 9). MEASURED, `leiaLevelSetFoam`, reversed flows, 2D vortex
> N = 64/128/256 and 3D shear N = 32/64/128: B1 FAILS. The extension (second distance) halves the
> band gradient error at T/2 (ratio 0.34 to 0.46), but the shape error at T is 4.5x/5.1x (2D) and
> 1.10x/1.56x (3D) of the baseline, with the shape order 2.53 against 3.23 (2D) and 1.28 against
> 1.87 (3D); the volume error at T is 8x (2D) and 11x (3D) at the finest rung. The band gradient
> error converges in no arm. With the small source (linearQ, mu 0.1) the gradient error falls
> further (B3 PASS as written) and the shape and the volume become worse again. The 3D deformation
> study is void (setup error) and runs again.

## What it is

$P(x) = x - \mathrm{dist}\,n$, $U_{ext} = \delta U(x) + (1-\delta)U(P(x))$, with dist from the first
or the second-order rescaled distance ([[concepts/algebraic-rescaling-eulerian]]), $U(P)$ by a
first-order Taylor step from the cell centre, $\delta$ a raised cosine from 0 at $|dist| \le 1.5h$
to 1 at $|dist| \ge 4h$. It is the construction of haloLimited (HL0, failed 2026-09-27) without
the one-cell cap on the travel.

## Why it matters

D. Bothe's proposal: a velocity that is constant along the normals near the interface keeps
$|\nabla\psi| = 1$ there, so a small source can hold the profile.

## Where in the code

`VELOCITY_EXTENSION projectionBlend` (the velocity-extension family), tokens `PB_DISTANCE`,
`PB_W1`, `PB_W2`; the read-out B1 to B4 in the headers of `config/projectionBlend/*.yaml`;
`workflow/scripts/projection_blend_metrics.py`.

## Evidence

| claim | number | where |
|---|---|---|
| band gradient error at T/2, ext / base | 0.34, 0.46 (2D); 0.40, 0.43 (3D) | `projectionBlend_verdicts.csv` B1 |
| shape at T, ext / base | 4.45, 5.10 (2D); 1.10, 1.56 (3D) | same |
| shape order, ext against base | 2.53 / 3.23 (2D); 1.28 / 1.87 (3D) | `projectionBlend_orders.csv` |
| volume at T, finest rung | 1.47e-3 / 1.85e-4 (2D); 5.1e-2 / 4.7e-3 (3D) | same |
| band mean q at T/2 | 2.4 (2D), 3.2 (3D) without; 1.45 to 1.9 with | same |
| null pairs, states | 3/3 bit-identical, all COMPLETED | B4 |

## Why it failed, or why we think so

HYPOTHESIS: the baseline's high shape order at T comes from the time reversal of a linear
transport with a prescribed velocity; an extension that depends on psi breaks that cancellation.
In 2D the extension reduces the volume error at T/2 16 times, which points the same way. A flow
that does not reverse separates the two. Earlier hypothesis (laptop smoke): the extended field is
not divergence-free and the advective form mixes a limited upwind and a centred term.

## Decisions

## Open questions

- The 3D deformation re-run (driver 55814176).
- A non-reversing flow (`2Dtranslation`, a rotation) to test the hypothesis above.

## Related

[[concepts/algebraic-rescaling-eulerian]], [[concepts/static-source-gate-eulerian]],
[[hubs/gradient-control]]

## Log

### 2026-10-29
Created with the model and the pre-registered read-out B1 to B4.

### 2026-10-31
2D vortex and 3D shear measured (B1 FAIL, B2 recorded, B3 PASS as written, B4 PASS). The 3D
deformation study had no `law` block in its case template: every gradientControl arm stopped at
step 0. It is renamed `_VOID_noLawBlock_2026-10-31`, the template is fixed (9d483075), and the
study runs again.
