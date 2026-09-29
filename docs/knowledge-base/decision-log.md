---
title: "Decision log"
description: "One line per settled decision, in time order - what was decided, the number, the gate"
kind: index
status: settled
part: all
tags: [index]
date: 2026-09-28
---
# Decision log

[[index]] <- back

One line per settled decision: `date · [[decisions/note]] · SETTLED claim: number · gate or config`.
Lines are never edited; a correction is a new line (`CORRECTED`, `REOPENED`) that links the old
one. Newest month at the bottom.

## 2026-07

- 2026-07-31 · [[decisions/surface-tension-reconstructed-curvature]] · SETTLED `SURFACE_TENSION_FORCE reconstructedCurvature`, `FACE_CURVATURE_SOURCE model`: constant curvature absorbed to 3.770e-9 to 1.039e-8 m/s on uniform meshes; Laplace jump 145.470 Pa in every well-balanced arm · `workflow/scripts/run_pressure_compatibility_gate.py`, `config/stationaryDroplet3DrefinedWB.yaml`

## 2026-08

- 2026-08-07 · [[decisions/curvature-inverse-gaussian]] · SETTLED `CURVATURE_INVERSE_GAUSSIAN yes`: K-aware inverse h^1.95 (5.25e-3 at N = 128) against h^1.02 without K on the 3D sphere; K off, two of three coupled 3D arms diverge · `config/faceCurvatureSphere3D.yaml`, `config/stationaryDroplet3DwideNoK.yaml`
- 2026-08-20 · [[decisions/momentum-schemes-bdf2-upwind]] · SETTLED `MOMENTUM_DDT_SCHEME backward`, `RHO_DDT_SCHEME backward`, `MOMENTUM_DIV_SCHEME upwind`: BDF2 vs Euler gain +11.1 / +2.9 / -3.0 % (noise); upwind inert to four significant figures on the stationary droplet · `config/upwindConvection2D.yaml`, `config/upwindConvection3D.yaml`
- 2026-08-25 · [[decisions/psi-filter-none]] · SETTLED `PSI_FILTER none`: the theta = 0.2 band filter is 5.86x better at R/h = 15.8 and 1.61x worse at R/h = 10.0; score every candidate with all filters off · `config/filterOffAmplifier3D.yaml`
- 2026-08-27 · [[decisions/sl-reconstruction-uncached-qwls]] · SETTLED `SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares`: shape order 2.84 (2D vortex, CFL 0.5), 2.95 / 3.28 (3D shear, hex / poly) · `config/advConv2Dvortex.yaml`
- 2026-08-31 · [[decisions/sl-trace-velocity-projected-flux]] · SETTLED `SL_TRACE_VELOCITY projectedFlux`: rate -52.0 1/s against +118.3 1/s for cellCentred (N = 128, full horizon); reconstruct operator 70 %, solenoidality 30 % · `config/fullHorizonStability2D.yaml`

## 2026-09

- 2026-09-01 · [[decisions/mass-flux-alphaf-donor-plane]] · SETTLED `MASS_FLUX_ALPHAF_SOURCE donorPlane`: author instruction, "only Gauss upwind worked in rhoLENT"; no measurement · author decision 2026-09-01
- 2026-09-01 · [[decisions/phase-indicator-detrixhe-aslam]] · SETTLED `PHASE_INDICATOR detrixheAslam`: static-circle volume error order about 2.0, equal to the geometric clip to eight digits; the default change is not inert · `config/phaseIndicatorConvergence.yaml`
- 2026-09-01 · [[decisions/curvature-extension-cell-centre-inverse]] · SETTLED `CURVATURE_EXTENSION cellCentreInverse` (case-dependent): stationary 2D residual 4.60 / 3.81 / 1.55 / 1.49x lower than none at N = 32 to 256; translating undecided, Popinet runs none · `config/stationaryLadder2Dshared.yaml`
- 2026-09-02 · [[decisions/mass-flux-rholent]] · SETTLED `MASS_FLUX rhoLENT`: capillary residual +1.0 / -22 / +0.1 % at N = 32 / 64 / 128, volume and shape to three digits (stationary only; translating rationale void) · `config/rhoLENTStationary2D.yaml`
- 2026-09-03 · [[decisions/viscosity-face-model-alg-lin]] · SETTLED `VISCOSITY_FACE_MODEL alg_lin` (2D; 3D open): the only model with positive orders in all three metrics at both ratios; L1 1.400e-3 at N = 512, order 1.10 at ratio 1000 · `config/mufGrid2D_jump1000.yaml`
- 2026-09-05 · [[decisions/sl-fit-normal-equations]] · SETTLED `SL_FIT normalEquations`: householderQR blows up identically (step-3 phase volume 0.017512193 vs 0.017512208); worst amplifier pivot 0.757 · `config/popinet3D_La12000_poly_dump4_qr.yaml`
- 2026-09-09 · [[decisions/mesh-family-hexahedral]] · SETTLED mesh family hexahedral: Lambda_max 1.0527 / 1.0549 / 1.0566 on blockMesh / cartesianMesh / snappy+layers against 1.2608 on pMesh; rho(B) - 1 4.41e-03 against 1.03e-02 · `workflow/scripts/fit_amplification_probe.sh`, `config/popinet3D_poly_sigma0_dtSweep.yaml`
- 2026-09-10 · [[decisions/sl-clip-and-value-bound-off]] · SETTLED `SL_CLIP false`, `SL_VALUE_BOUND none`: clip with exemption fails at step 506 vs 527; global clip +30.4 % volume; cone bound 3.6x to 189.7x worse on pure advection · `config/popinet3D_poly_sigma0_clipGate.yaml`, `config/advConv2Dvortex.yaml`
- 2026-09-26 · [[decisions/process-gates-2d-first-and-no-best-yaml]] · SETTLED process: the 2D method gate before any coupled study, 3D on a pass; no config/best.yaml (shallow merge lost SL_FIT, 2026-09-09); METHOD.md updated in the same commit · `config/gates/methodGate2D.yaml`
- 2026-09-27 · [[decisions/mass-flux-bound-rho]] · REOPENED `MASS_FLUX_BOUND_RHO true`: the basis rhoDdtGate2D (rho to -72.28 without the bound) ran on the closed box and is void; the bound stays active and unmeasured on the repaired case · `config/rhoDdtGate2D.yaml` (VOID)
