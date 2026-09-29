---
title: "The pressure projection and the linear solvers"
description: "The rAUf-weighted Poisson projection with fvc::reconstruct absorbs the capillary flux; on orthogonal hexahedra the linear solver is exonerated (a 30 times tighter solve moves the result by 2.4e-4), on 10 percent perturbed meshes strict PCG gains 18 to 29 times and a common floor of 3.5e-5 to 7.7e-5 m/s remains that only the reconstruct operator can explain (2026-09-28)."
aliases: [pressure projection, linear solver tolerance, GAMG against PCG, cancellation-dominated convergence, rAUf]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-08-07
decided_by: [docs/plan-curvature-stabilization.md section 8, workflow/Snakefile.pressure-compatibility, cases/stationaryDroplet2D/system/fvSolution.template]
code: [applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H, applications/solvers/leiaLevelSetTwoPhaseFoam/YoungLaplaceEqn.H, cases/stationaryDroplet2D/system/fvSolution.template]
sources: [METHOD 5, CLAUDE solver convergence, PCS 8, PCS 16.2, SL article sec:fluxresidual, SL article sec:droplet, RM pressure gates 2026-07-28 and 2026-07-30, STATUS 4 hanging-node gate, STATUS 4 Young-Laplace fix, SL deck 6/7, PSH 0g]
---
# The pressure projection and the linear solvers

> Verdict (2026-09-28). The projection solves `laplacian(rAUf, p_rgh) == div(phiHbyA)` with the capillary flux inside `phiHbyA` and recovers the velocity from the face difference through `fvc::reconstruct` ([METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L253-L257), [`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L58-L108)). The balance is a cancellation problem, and OpenFOAM normalises the residual by the large quantity, so a relative tolerance is safe only when it is far below the fraction of the source that carries the signal ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L731-L743)). On uniform orthogonal hexahedra that margin is six orders and the exoneration is direct: a strict `DICPCG` at `1e-11`, `relTol 0` (about 265 iterations per solve against 4 to 17 for the production GAMG) moves the non-absorbable fraction, the velocity residual and `max|U|` by at most `2.4e-4` relative over a run ([PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L556-L577)), and a swap of GAMG for PCG reproduces the blow-up time to 4 to 5 digits through a chaotic amplifier ([PCS 16.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1342-L1346)). That does not transfer to non-orthogonal meshes: on 10 percent perturbed meshes with exact curvature the strict solve gains 18.2 and 28.8 times at N = 64 and 128 ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1132-L1137)), and after the solver and the momentum predictor have each removed their share a common floor of about `3.5e-5` and `7.7e-5` m/s remains that is not algebraic and points at `fvc::reconstruct` ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1226-L1269), [METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L794-L799)).

## What it is

The production settings of the droplet cases: `p_rgh` GAMG with DIC smoother, tolerance `1e-9`, `relTol 0.01`; `p_rghFinal` with `relTol 0`; `U` a symmetric Gauss-Seidel smoother at `1e-8`; two correctors, one non-orthogonal corrector, `pRefCell 0` ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L67-L92), [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L186-L190)). OpenFOAM's normalisation factor is `sum(|A psi - A psibar| + |b - A psibar|)`, which at convergence is `2 sum|b|`: the spread of the two sides about the uniform field, so that a solver gets no credit for doing nothing ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L737-L742), [SL deck 6/7](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/6/7)). In a balanced-force problem the net force is a small difference of two large fluxes, 0.02 to 0.15 percent of a 1.7 m/s capillary predictor at N = 128, so the reference is the large quantity ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L733-L736), [[concepts/balanced-force-csf-flux]]). Only the final corrector of the final outer iteration uses `relTol 0`; with `psi` re-advected inside the outer loop, an under-converged intermediate velocity moves `psi` and `kappa` before any tight solve ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L752-L754), [[concepts/psi-outer-correctors]]).

## Why it matters

Two decisions depend on the margin being measured, not assumed. First, the linear solver is not a lever on orthogonal meshes: tightening `p_rgh` from `1e-9` to `1e-12` leaves the settled currents unchanged to four significant figures (`5.709e-4` and `1.250e-5` m/s at N = 32 and 64) at about 13 times the iterations, which is why the fast tolerance is the default ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L71-L77), [SL article `sec:droplet`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1679-L1683)). Second, the decomposition dependence of a parallel run, about `1e-4` relative between one and four ranks, is the solver's: the pressure solve converges to a relative tolerance and the GAMG agglomeration follows the partition, identically on a uniform mesh and on a refined one with hanging nodes ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1236-L1258), [[concepts/seam-checks-and-decomposition-invariance]]). A diverging smoother destroys diagnosability: a momentum smoother reaching `1e+98` at 1000 iterations makes a physical blow-up indistinguishable from a linear-algebra one in the log ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L755-L759), [[concepts/log-classifier-and-waiters]]).

## Where in the code

- The projection and its diagnostics: [`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L58-L108); the flux-space residual per final corrector ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L170-L215)).
- The t = 0 Young-Laplace solve was unreferenced on a singular system and inherited `relTol 0.01`, leaving the initial pressure about 1 percent off; both fixed. Its unit Laplacian coefficient against the `rAUf`-weighted one of `pEqn` (a factor of about 839 across the interface) cannot be fixed before the first momentum equation exists ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L960-L966), [PSH 0](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L66-L69)).
- The mass-flux projection solves a pure Poisson problem with PCG at `1e-12`, `relTol 0`, because there the correction is the small quantity ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L20-L31), [[concepts/mass-flux-projection]]).
- The pressure gates run from `workflow/Snakefile.pressure-compatibility` ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1107-L1117)).

## Evidence

| claim | number | where |
|---|---|---|
| the solver is exonerated on orthogonal hexahedra | GAMG (1e-9, relTol 0.01, 4 to 17 iterations) against DICPCG (1e-11, relTol 0, about 265): R_f/phig, U_res and max abs U within 2.4e-4 relative, N = 128 to t = 0.02 | [PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L556-L577), MEASURED |
| the settled currents are curvature-limited | 1e-9 against 1e-12: 5.709e-4 and 1.250e-5 m/s at N = 32 and 64 either way | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1679-L1683), MEASURED |
| the frozen-interface current is solver-independent | unchanged to eight significant figures between relTol 1e-2 and absolute 1e-12 | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1736-L1739), MEASURED |
| GAMG against PCG through the amplifier | t_blow reproduced to 4 to 5 significant digits after about 9000 steps | [PCS 16.2](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L1342-L1346), MEASURED |
| the non-orthogonal corrector count on perturbed meshes | 0 / 1 / 8 corrections at N = 128: 1.183e-2 / 1.146e-3 / 1.392e-3 m/s; eight suffice; convergence does not remove the current | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1070-L1086), MEASURED |
| strict PCG on perturbed meshes | 9.050e-6 / 4.029e-5 / 4.839e-5 against GAMG 8.584e-6 / 7.328e-4 / 1.392e-3 at N = 32 / 64 / 128 (18.2x and 28.8x); 1e-13 unchanged; strict GAMG takes an FPE | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1132-L1139), MEASURED |
| the momentum predictor is a second amplifier of the mesh-induced residual | predictor off: perturbed 8.584e-6 to 2.725e-6, 7.328e-4 to 3.137e-5, 1.392e-3 to 6.121e-5 (23x); variable-curvature replay unchanged to four figures | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1181-L1224), MEASURED |
| the two levers land on a common floor | N = 64: 3.1e-5 to 4.0e-5 for all three corrected variants; N = 128: the combination is 1.58x worse than PCG alone; the floor agrees to 0.4 percent at 1e-9, 1e-11, 1e-13 | [RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1226-L1252), MEASURED |
| the decomposition dependence is the solver's | np 1 against np 4 at 1e-4 relative on the uniform and the refined mesh; maxima over time equal to four digits | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1248-L1258), MEASURED |
| the first-step floor with exact curvature | 1.09e-8 m/s here; SAAMPLE's residual-driven pressure iteration reaches 1e-13 from step 1 at 7 to 19 corrections | [PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L755-L771), MEASURED in the reference |

## Why it failed, or why we think so

On skewed meshes the exact, constant curvature leaves a residual that is insensitive to the force form, the geometry, the corrector count, `rAUf`, the predictor, the solver and its tolerance. The one operator common to every arm and never removed is the `fvc::reconstruct` of the velocity correction: it assumes that the face values sample one smooth cell field, which fails at an interface cell where `grad(p)/rho` jumps by the density ratio, and on a skewed mesh it mixes tangential into normal information ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1254-L1269)). The same operator is the output stage of the amplifier on orthogonal meshes: SAAMPLE shows it diverges at the interface for fields with a gradient jump ([PSH 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L791-L805), [[concepts/trace-velocity-projected-flux]]).

## Decisions

- Production GAMG at `1e-9`, `relTol 0.01` for orthogonal droplet studies; strict PCG/DIC only on non-orthogonal meshes ([PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L600-L602)).
- At least eight non-orthogonal corrections on `mesh: perturbed`, enforced by the materialisation and the oracle runners ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1088-L1094)); the token default stays 1, and on the polyhedral rungs 1 / 3 / 6 correctors change the velocity norms by at most `5e-6` relative ([`cases/default.parameter`](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L1100-L1105)).
- `momentumPredictor no` is not a default: it is 1.8 times worse at N = 32 on the uniform mesh and changes the effective time integration ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1282-L1287)).

## Open questions

1. The skew-mesh floor and the facewise gate that separates source from gain, `phig - p_rghEqn.flux()` before the reconstruct on the frozen exact circle ([RM](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1271-L1280), [METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L794-L799)).
2. The margin on polyhedral and strongly non-orthogonal meshes must be re-established before quoting solver-converged results there ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L746-L750), [[decisions/mesh-family-hexahedral]]).
3. Residual-driven outer iteration (`OUTER_U_TOL`, `OUTER_P_TOL`) as the route to SAAMPLE's first-step floor ([[models/semi-implicit-capillary-force]]).

## Related

[[hubs/surface-tension]], [[concepts/balanced-force-csf-flux]], [[concepts/well-balanced-exact-curvature-gate]], [[concepts/parasitic-current-mechanism]], [[concepts/psi-outer-correctors]], [[concepts/trace-velocity-projected-flux]], [[concepts/static-local-refinement]], [[concepts/mass-flux-projection]], [[concepts/seam-checks-and-decomposition-invariance]], [[concepts/log-classifier-and-waiters]], [[models/semi-implicit-capillary-force]], [[decisions/mesh-family-hexahedral]], [[studies/poly3d-roadmap]].

## Log

### 2026-09-28
Created from METHOD 5, CLAUDE's solver-convergence section, PCS 8 and 16.2, the roadmap's pressure gates of 2026-07-28 and 2026-07-30 and STATUS 4.
