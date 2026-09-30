---
title: "VISCOSITY_FACE_MODEL: alg_lin, decided in 2D; 3D open"
description: "VISCOSITY_FACE_MODEL alg_lin in the global default, decided 2026-09-03 on the 36-arm mufGrid2D ladder: the only model with a positive order in all three metrics at both viscosity ratios, L1 1.400e-3 at N = 512 with order 1.10 at ratio 1000; the 3D droplet templates carry no token and run the solver default geo_lin"
aliases: [VISCOSITY_FACE_MODEL alg_lin]
kind: decision
status: settled
part: viscosity
tags: [decision, part/viscosity]
date: 2026-09-28
date_settled: 2026-09-03
decided_by: [config/mufGrid2D_jump1000.yaml, config/mufGrid2D_jump1000_N128.yaml, config/mufGrid2D_jump1000_N512.yaml, config/mufGrid2D_prod.yaml, config/mufGrid2D_prod_N128.yaml, config/mufGrid2D_prod_N512.yaml]
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/faceViscosity.H, applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H, cases/default.parameter]
sources: ["METHOD 8.1 row VISCOSITY_FACE_MODEL (L397)", "DP L705-L767", "STATUS 0 (L139-L143)", "STATUS 11.13 (L3602)", "SL article sec:viscous (L760-L900)"]
---
# VISCOSITY_FACE_MODEL: alg_lin, decided in 2D; 3D open

> `VISCOSITY_FACE_MODEL alg_lin` in the global default, decided on 2026-09-03 on a 36-arm ladder: six models, N = 128 / 256 / 512, viscosity ratios 54.83 and 1000, all with the explicit viscous term ([DP L729-L752](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L729-L752), [METHOD 8.1 L397](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L397)). `alg_lin`, the linear interpolation of the cell indicator mixed linearly into `mu_f`, is the only model with a positive order in all three metrics at both ratios, has the lowest error at N = 512 at both, and the travelled fraction closest to 1; at ratio 1000 its L1 error is 1.400e-3 with order 1.10 and correlation 0.998 ([DP L742-L752](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L752)). A brief switch to `geo_lin` is retracted: the algebraic arms had run with a face viscosity frozen at t = 0, a bug fixed in 39e59b3 ([DP L754-L759](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L754-L759)). Popinet's face properties are the same simple average ([STATUS L139-L143](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L139-L143)). OPEN in 3D: the four 3D droplet templates carry no `viscosityFaceModel` entry and run the solver's code default `geo_lin` ([METHOD 8.1 L397](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L397), [STATUS L3587](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3587)).

## The question

The viscous term carries one face coefficient `mu_f`. Two independent choices build it: where `alpha_f` comes from (`alg_`: the cell indicator interpolated linearly; `geo_`: the geometric face area fraction that also builds `rho_f`), and how `mu_f` is mixed (`_lin`: `alpha_f mu1 + (1 - alpha_f) mu2`; `_harm`: `1/mu_f = alpha_f/mu1 + (1 - alpha_f)/mu2`; `_blend`: `w mu_harm + (1 - w) mu_lin` with `w = abs(n . S_f) / abs(S_f)`, no free parameter) ([DP L705-L728](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L705-L728)). An earlier token with three values mixed the two axes and left `alg_harm` untested ([mufGrid2D_jump1000.yaml L1-L17](https://github.com/leia-openfoam/leia/blob/8867581/config/mufGrid2D_jump1000.yaml#L1-L17)). The ladder crossed them.

## The measurement that decided it

| arm | metric | value | where |
|---|---|---|---|
| ratio 1000, N = 128 / 256 / 512, translating droplet 2D, t = 0.02 s: `alg_lin` | L1 at N = 512; p(L1), R; p(L2); p(shape); travelled fraction | 1.400e-03; 1.10, 0.998; 1.15; +0.71; 1.0024 | MEASURED, [DP L742-L747](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L747) |
| ratio 1000: `alg_blend` / `geo_lin` / `geo_harm` | L1 at N = 512; p(L1); p(shape) | 3.086e-03, 0.83, +0.60 / 1.834e-03, 0.46, -0.37 / 4.670e-03, 0.15, -0.56 | MEASURED, [DP L742-L747](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L742-L747) |
| ratio 54.83: `alg_lin` / `geo_lin` / `alg_harm` / `geo_harm` / `alg_blend` / `geo_blend` | L1 at N = 512; p(L1) | 3.494e-03, 0.77 / 4.953e-03, 0.15 / 6.479e-03, 0.16 / 6.027e-03, 0.08 / 6.247e-03, 0.26 / 5.669e-03, 0.11 | MEASURED, [DP L733-L740](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L733-L740) |
| shape order across the grid | sign | negative for every model except `alg_lin`, and `alg_blend` at the large ratio | MEASURED, [DP L749-L752](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L749-L752) |
| the retracted `geo_lin` switch | cause | `muf` constructed once and updated only in an `if (geometric)` branch; fixed in 39e59b3 | MEASURED, [DP L754-L759](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L754-L759) |

Pre-registered read-out ([mufGrid2D_jump1000.yaml L29-L34](https://github.com/leia-openfoam/leia/blob/8867581/config/mufGrid2D_jump1000.yaml#L29-L34)): (i) does the `alpha_f` source or the mixing model dominate; (ii) does the blend beat both of its limits. Neither harmonic nor blend beat linear anywhere ([DP L761-L762](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L761-L762)).

The ladder ran on the translating droplet case after the closed-box fix 440107f (the configs date from 2026-09-03), at np 4 to 16 before the coupled-face fixes of 2026-09-27 ([mufGrid2D_jump1000.yaml L35-L52](https://github.com/leia-openfoam/leia/blob/8867581/config/mufGrid2D_jump1000.yaml#L35-L52), [[concepts/coupled-face-density-defect]]).

## What it does not cover

1. 3D, open. The templates of the three 3D droplet cases `stationaryDroplet3D`, `translatingDroplet3D` and `oscillatingDroplet3D` contain no `viscosityFaceModel` entry (`ellipsoidDroplet3D` is a static curvature case and needs none); the code default is `geo_lin` ([createFields.H L92](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L92)), and the banner prints `alg_lin (default)` ([printMethodBanner.H L139](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/printMethodBanner.H#L139)), so a 3D run reports a wrong default. The 3D gate needs the token before its verdict means anything ([[concepts/viscosity-open-items]]).
2. The parallel provenance of the ladder: it ran on 4 to 16 ranks before 28d13f0 and b1798c3; the decision holds until a re-run says otherwise ([[hubs/viscosity]]).
3. The Eulerian two-phase solver kept `muf` frozen until c094bd8 (2026-09-27); its earlier two-phase results carry that defect ([[concepts/eulerian-solver-mass-flux-port]]).
4. Why the harmonic and the blended models lose is not explained by the series-parallel argument ([SL article L845-L852](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L845-L852)).

## Related

[[hubs/viscosity]] - [[models/viscosity-face-model]] - [[concepts/viscosity-open-items]] - [[concepts/coupled-face-density-defect]] - [[concepts/eulerian-solver-mass-flux-port]] - [[concepts/rholent-mass-flux]] - [[decisions/mass-flux-rholent]] - [[cases/translating-droplet]] - [[cases/popinet-translating-droplet]] - [[studies/sl-quadratic-pre-print]] - [[decision-log]]

## Log

### 2026-09-28
SETTLED in 2D on 2026-09-03; 3D open. Entered in [[decision-log#2026-09]].

### 2026-09-29
CORRECTED: three 3D droplet templates lack the token, not four (ellipsoidDroplet3D is static).
