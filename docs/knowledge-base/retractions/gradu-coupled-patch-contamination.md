---
title: "31 parallel kinematic studies contaminated by gradU on coupled patches (2026-08-26)"
description: "VOIDED 2026-08-26 - every parallel run of the kinematic semi-Lagrangian solver: setVelocity wrote face values into processor patches, fvc::grad(U) was biased O(1) at every seam, and the default Taylor foot consumes that gradient; 31 studies, 15 with curated tables; the 2D vortex re-run moved the endpoint at most 2x and kept second order"
aliases: [gradU contamination, setVelocity coupled patches, 30e6ba9]
kind: retraction
status: voided
part: advection
tags: [retraction, part/advection]
date: 2026-09-28
code: [src/leiaLevelSet/velocityModel/velocityModel.C, src/leiaLevelSet/semiLagrangian/pointValueScheme.C, cases/1Dstretch, config/uncachedConv2Dvortex.yaml]
sources: [docs/gradU-coupled-patch-contamination.md, STATUS 4 (2026-08-26 and 2026-08-27), PCS 0 contamination notice, METHOD 8, commit 30e6ba9]
---
# 31 parallel kinematic studies contaminated by gradU on coupled patches (2026-08-26)

> VOIDED 2026-08-26. The claims were the parallel kinematic semi-Lagrangian baselines: the transport ground truth of the plan (2D reversed vortex shape order 2.97 and volume 2.98 at CFL 1/2, 2.59 / 2.94 at CFL 1; 3D shear 2.50 / 3.25; 3D deformation 1.36 / 1.86; polyhedral deformation 1.46 / 2.51, [plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L21-L32)) and the same numbers in [METHOD 8](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347). The measurement: `velocityModel::setVelocity` wrote a face-centre value into every patch, coupled patches included. A coupled patch field holds the neighbour cell value, so `fvc::grad(U)` carried an error of order a (the strain), not order h, in every processor-adjacent cell, for the life of the code ([handover note, section 1](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L18-L43)). On `cases/1Dstretch`, whose solution is closed-form, the band mean of dpsi/dx with the source `R` read 1.004681 / 1.005367 / 0.956834 at np 2 / 4 / 8 against the exact 1.000000, an error of 4.3e-02 where serial was exact; after the fix every count reads 1.000000 ([section 1](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L47-L60)). The default foot integrator consumes `fvc::grad(U)` in its dt^2/2 term (`pointValueScheme.C:353`), and the kinematic solver feeds it the prescribed velocity, so all 31 parallel kinematic SL studies are contaminated, 15 of them with curated tables ([section 2 and 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L64-L146), [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L619-L640)). Fixed in commit [30e6ba9](https://github.com/leia-openfoam/leia/commit/30e6ba9). Scope of the void: `solver: leiaSemiLagrangeLevelSetFoam` and np > 1. Serial runs and every two-phase SL study are outside it: the two-phase solver takes U from the momentum solve and has no `velocityModel` ([section 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L148-L156)).

## The claim, and where it lived

- [plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L21-L61), `docs/plan-curvature-stabilization.md`: the transport table "do not touch", now under a contamination notice and a re-established 2D row.
- [METHOD 8](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347), `METHOD.md`: "Transport (prescribed velocity, settled)", 2.97 / 2.59, 2.95 / 3.28, 1.36 / 1.46. No marker.
- The SL article, `docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`, [the 2D convergence section](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1278-L1285): shape order 2.97 at CFL 1/2 and 2.59 at CFL 1.
- The curated tables `uncachedConv*`, `linearConv*`, `nslConv*`, `npslConv2Dvortex`, `sdCompare2D` under the SL, linear-SL and nPSL themes ([the list](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L110-L142)).
- [STATUS 4, 2026-08-18](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L513-L518) and [L543-L550](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L543-L550): the `linearTaylor` kinematic price arm (volume 340 % at N = 256), "numbers suspect".
- The SDPLS article, `docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex`, [sec:coupledpatch](https://github.com/leia-openfoam/leia/blob/8867581/docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex#L2678-L2679) `sec:coupledpatch`: the same defect on the Eulerian SDPLS side (29 studies, tiered for re-run).

## Why it was wrong, or why we think so

| claim | number | where |
|---|---|---|
| A prescribed analytic velocity is exact on every patch. | On a coupled patch the slot holds the neighbour CELL value; the face value gives 0.5 U(x_P) + 0.5 U(x_f) instead of U(x_f), an error a h/4 per face and O(a) in the Gauss gradient. | DERIVED, [handover note 1](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L28-L38) |
| The kinematic solver is np-invariant. | `1Dstretch`, source `R`, band mean dpsi/dx at t = 1: 1.000000 serial, 1.004681 np 2, 1.005367 np 4, 0.956834 np 8; the sourceless arm 0.367872 at every count. After the fix: 1.000000 at every count (max error 4.0e-10). | MEASURED, [handover note 1](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L47-L60), [commit 30e6ba9](https://github.com/leia-openfoam/leia/commit/30e6ba9) |
| Only the SDPLS source consumes `fvc::grad(U)` (the commit message). | `pointValueScheme.C:353` (the default foot), `closestPoint.C:432`, `steadyUpwindLinear.C:75`, `normalProjectedScheme.C:161` and the signed-distance QWLS reconstruction consume it too. | MEASURED, [handover note 5](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L228-L234) |
| The 1e-12 serial-to-np-4 agreement proved seam cleanliness. | It tested the reconstruction exchange, not the `grad(U)` path. | DERIVED, [plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L45-L48) |
| The size of the effect on the published orders. | 2D vortex N = 256, CFL 0.5, July code run serial on the July inputs: shape 1.105e-06 and volume 3.95e-05 against the curated buggy-parallel 7.851e-07 and 8.01e-05, so the seam bias moved shape 1.4x and volume 2x, with no order change. The re-established 14-arm orders: 2.84 / 3.30 at CFL 0.5, 2.38 / 3.54 at CFL 1. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L585-L592), [plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61) |
| The first re-run's 117x collapse was the algorithm. | It was the measurement layer: the metrics writer scored t = 1.9997 instead of t = 2, and a band frozen at t ~ 1.0018 was dumped as the endTime alpha. Retracted the same day; psi at t = 2 is bit-identical between the July build and HEAD, serial. | MEASURED, [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L593-L617) |

## What survives

1. Every serial kinematic result, every two-phase SL study, and every consumer of `phi`: the flux was written per face and was always exact ([handover note 3](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L148-L156)).
2. Second-order shape transport of the 2D reversed vortex as a method property, re-established on the fixed binaries at 2.84 / 3.30 (CFL 0.5), with the N = 256 endpoint equal to the bug-free serial run to all digits ([plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L52-L61)); the table `uncachedConv2Dvortex_errors.csv` was re-baselined in commit [8693326](https://github.com/leia-openfoam/leia/commit/8693326).
3. The argument for an RK2 foot integrator: it needs velocity samples, not a velocity gradient ([handover note 4](https://github.com/leia-openfoam/leia/blob/8867581/docs/gradU-coupled-patch-contamination.md#L189-L199)), see [[concepts/departure-foot-ab2-centring]].
4. The defect class entered the [4-rank rule](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L216-L223), see [[concepts/seam-checks-and-decomposition-invariance]].

## Propagation (checklist, same commit)

Done:

- [x] `STATUS.md`: the [INVALIDATION](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L619-L640) and the [closure for the 2D vortex](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L580-L617); the two "numbers suspect" flags ([L513-L518](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L513-L518), [L543-L550](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L543-L550)).
- [x] The plan document: the [contamination notice and the re-established row](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L34-L61).
- [x] The handover note `docs/gradU-coupled-patch-contamination.md` with the 31-study list.
- [x] `CLAUDE.md`: the [class list](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L216-L223).
- [x] The curated 2D vortex table re-baselined (commit 8693326, 2026-08-27); the raw July run preserved as `studies/uncachedConv2Dvortex.preGradUfix-20260826`.
- [x] The line in [[retraction-log]].

Still missing:

- [ ] The 3D shear, 3D deformation and polyhedral rows are "pending re-run" ([plan-curvature-stabilization 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L61)); the record has no entry of their landing, and their curated tables (`uncachedConv3D*`, `linearConv3D*`) carry no marker.
- [ ] [METHOD 8](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L343-L347) and the [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1278-L1285) quote the pre-fix 2.97 at CFL 1/2; the re-established value is 2.84. Which table the article's `tab:orders` reads is not checked here.
- [ ] The outcome of the SDPLS re-run tiering (29 studies, P0 to P3) is not recorded in this vault.

## Related

Hubs: [[hubs/advection]], [[hubs/verification]]. Siblings: [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/departure-foot-ab2-centring]], [[models/sl-scheme]], [[concepts/psi-outer-correctors]], [[cases/kinematic-advection-cases]], [[concepts/advection-regression-set]], [[models/sdpls-source]], [[studies/sdpls-pre-print]], [[retractions/psi-filter-seam-bug]].

## Log

### 2026-09-28
Written from docs/gradU-coupled-patch-contamination.md, STATUS 4 (2026-08-26 and 2026-08-27) and the plan's contamination notice. Voided 2026-08-26. Entered in [[retraction-log#2026-08]].
