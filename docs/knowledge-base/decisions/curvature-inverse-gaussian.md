---
title: "CURVATURE_INVERSE_GAUSSIAN: the K-aware parallel-surface inverse"
description: "CURVATURE_INVERSE_GAUSSIAN yes, decided by config/faceCurvatureSphere3D.yaml on 2026-08-07: the inverse with the Gaussian curvature converges at h^1.95 on the 3D sphere against h^1.02 without K; with K off the coupled 3D arms at R/h = 12.7 and 15.8 diverge"
aliases: [CURVATURE_INVERSE_GAUSSIAN yes]
kind: decision
status: settled
part: surface-tension
tags: [decision, part/surface-tension]
date: 2026-09-28
date_settled: 2026-08-07
decided_by: [config/faceCurvatureSphere3D.yaml, config/stationaryDroplet3DwideNoK.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/stabilizedFootPointFaceCurvature.H, cases/default.parameter]
sources: ["METHOD 4.3 (L222-L241, CORRECTED)", "METHOD 8.1 row CURVATURE_INVERSE_GAUSSIAN (L398)", "RM L1362-L1396", "STATUS 4 (L846-L890)", "DP L296-L303"]
---
# CURVATURE_INVERSE_GAUSSIAN: the K-aware parallel-surface inverse

> `CURVATURE_INVERSE_GAUSSIAN yes` in the global default; switchable since e6edf3d (2026-08-19) ([DP L296-L303](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L296-L303), [METHOD 8.1 L398](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L398)). The parallel-surface inverse that converts the curvature of the level set at offset `d` into the interface curvature is `kappa^Gamma = (kappa - 2 K d) / (1 - d kappa + K d^2)`, with the Gaussian curvature `K = (g . cof(H) . g) / abs(g)^4` from the same fit; `K` is identically zero in 2D, so every 2D result is byte-identical ([METHOD 4.3 L235-L241](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L235-L241), [RM L1362-L1372](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1362-L1372)). Decided by the static 3D sphere gate of 2026-08-07: the K-aware inverse converges at h^1.95 (L2 5.25e-3 at N = 128) against h^1.02 (0.421, equal to the raw value) for the 2D scalar inverse without K ([RM L1373-L1396](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1373-L1396)). The coupled control of 2026-08-19 confirms it: with K off the delivered non-gradient content is 770x larger at zeroth order and two of three 3D arms diverge ([STATUS L846-L890](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L846-L890)).

## The question

The 2D parallel-curve inverse `kappa = kappa_d / (1 - d kappa_d)` is first-order wrong in 3D, where the exact parallel-surface relation `kappa_d = (kappa + 2 d K) / (1 + d kappa + d^2 K)` needs the product of the principal curvatures `K = kappa_1 kappa_2` ([METHOD 4.3 L222-L234](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L222-L234)). On a sphere the inverse with K reduces exactly to `2 / (R - d)`; with K = 0 it gives `2 / (R - 2d)`, a relative curvature error `d / (R - 2d) = O(h)`, and the non-gradient content, which differences `kappa` across a face, is then `O(h) / h = O(1)`: zeroth order ([STATUS L870-L878](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L870-L878), DERIVED). The question was whether the K term is load-bearing, and later whether its fit noise feeds the 3D instability at R/h = 15.8 ([stationaryDroplet3DwideNoK.yaml L23-L34](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3DwideNoK.yaml#L23-L34)).

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| static sphere, exact SDF, serial, N = 32 / 50 / 80 / 128: raw cell-centre curvature, arithmetic face value | active-face L2 order; L2 at N = 128 | h^0.98; 0.419 1/m | MEASURED, [RM L1373-L1396](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1373-L1396) |
| the same, stabilized foot point with the K-aware inverse | order; L2 | h^1.95; 5.25e-3 (80x) | MEASURED, same |
| the same, FVM `div(grad psi / abs(grad psi))` raw / corrected | order | h^0.97 / h^1.91 | MEASURED, same |
| the same, 2D scalar inverse control (no K) | order; L2 | h^1.02; 0.421 = raw: the correction buys nothing without K | MEASURED, same |
| 2D cases with the switch on and off | every metric | byte-identical (K = +0 exactly on pseudo-2D fits) | MEASURED, [RM L1371-L1373](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1371-L1373), [config L36-L38](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3DwideNoK.yaml#L36-L38) |
| coupled 3D ladder, L = 6R, K on, N_L = 60 / 76 / 95 (R/h 10.0 / 12.7 / 15.8) | reached; max abs U; delivered non-gradient content at t = 0 | 0.1 s in all three; 1.0035e-03 / 1.0574e-03 / 7.0362e-02 m/s; 3.184 / 1.942 / 1.252, order +2.03 | MEASURED, [STATUS L852-L858](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L852-L858) |
| the same with K off | the same | 0.1 s / DIED 0.0955 s / DIED 0.0772 s; 3.0343e-01 / 9.2909e-01 / 1.2781e+00 m/s; 962.6 / 964.1 / 959.2, order +0.01 | MEASURED, [STATUS L852-L868](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L852-L868) |

Pre-registered read-out ([stationaryDroplet3DwideNoK.yaml L23-L34](https://github.com/leia-openfoam/leia/blob/8867581/config/stationaryDroplet3DwideNoK.yaml#L23-L34)): a. R/h = 15.8 stable with K off means K's noise feeds the growth; b. still unstable means K is exonerated and the tangential structure remains; c. the t = 0 delivery order drops to about 1 either way. Outcome b and c occurred: K is not the amplifier, and removing it destroys the delivery ([STATUS L880-L890](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L880-L890)).

## What it does not cover

1. The sphere gate scored the K-aware inverse inside the face delivery `stabilizedFootPointFace`; the production token now acts inside `cellCentreInverse` (sub-key `gaussianCurvature`) ([[models/curvature-extension]]). The coupled control of 2026-08-19 ran that cell delivery.
2. Second order is measured on constant curvature (the sphere) only; a varying-K surface (a torus) is not tested ([METHOD 4.3 L240-L241](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L240-L241), [[concepts/cell-centre-inverse-curvature]]).
3. The inverse assumes a parallel foliation, `psi = f(signed distance)`; on a drifted level set that premise fails, and the 3D instability at R/h about 16 happens despite the second-order delivery ([STATUS L880-L890](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L880-L890)).
4. The 3D offset correction of the reconstruction (`SL_OFFSET_CORRECTION`) is a different object and stays `none` in the droplet templates ([METHOD 4.1 L161-L167](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L161-L167)).

## Related

[[hubs/surface-tension]] - [[models/curvature-extension]] - [[concepts/cell-centre-inverse-curvature]] - [[concepts/curvature-from-the-fit]] - [[concepts/face-curvature-deliveries]] - [[concepts/parasitic-current-mechanism]] - [[decisions/curvature-extension-cell-centre-inverse]] - [[decisions/surface-tension-reconstructed-curvature]] - [[cases/curvature-static-gates]] - [[cases/stationary-droplet]] - [[studies/poly3d-roadmap]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED on the sphere gate of 2026-08-07, confirmed by the K-off control of 2026-08-19. Entered in [[decision-log#2026-08]].
