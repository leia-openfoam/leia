---
title: "The balanced-force capillary flux"
description: "The capillary force is one integrated scalar face flux in the pressure-flux space; a constant curvature is absorbed to round-off, the projection absorbs 99.85 to 99.98 percent of the flux, and the remainder is the alpha-weighted face gradient of the curvature (2026-09-28)."
aliases: [balanced-force CSF, capillary face flux, G_sigma]
kind: concept
status: settled
part: surface-tension
tags: [concept, part/surface-tension]
date: 2026-09-28
date_settled: 2026-09-04
decided_by: [config/kickOriginGate2D.yaml, config/amplifierGate2D.yaml, docs/plan-curvature-stabilization.md section 8]
code: [src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.H, applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H, applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H]
sources: [METHOD 5, PCS 8, SL article sec:surften, SL article sec:fluxresidual, RM frozen-circle gate, STATUS 0, STATUS hanging-node gate, SL deck 4/9]
---
# The balanced-force capillary flux

> Verdict (2026-09-28). The capillary force is kept as one integrated scalar face flux, `G_sigma,f = sigma kappa_f snGrad(alpha) |S_f|`, and it enters the pressure equation on the faces where the pressure gradient acts ([METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L243-L272)). A constant curvature is a discrete gradient that the pressure absorbs to round-off: `max|U|` is about `3e-11` m/s on a frozen circle ([METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L259-L262)), and the exact curvature cuts the step-1 kick of the translating droplet from `2.15e-3` to `1.69e-9` ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102)). With the reconstructed curvature the projection still absorbs 99.85 to 99.98 percent of the flux, and a 30x tighter pressure solve moves the remainder by at most `2.4e-4` relative ([PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L568-L583)). The remainder is structural: on a uniform mesh it is exactly the alpha-weighted face-normal gradient of the curvature field ([SL article `sec:fluxresidual`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2200-L2228)). The delivery, the density interpolation and the projection are therefore exonerated; the curvature error is the source ([[concepts/parasitic-current-mechanism]]).

## What it is

The continuum surface force is `f_sigma = sigma kappa grad(alpha)` ([SL article `sec:surften`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L919-L931)). The solver never builds it as a cell vector. Every member of [[models/surface-tension-force]] returns the owner-oriented face flux `G_sigma,f` ([`surfaceTensionForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.H#L29-L106)). The pressure equation then reads

1. `phig = (G_sigma - g_h snGrad(rho) |S_f|) rAUf`, added to `phiHbyA` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L58-L67));
2. `laplacian(rAUf, p_rgh) == div(phiHbyA)` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L81-L84));
3. `U = HbyA + rAU reconstruct((phig - p_rghEqn.flux())/rAUf)` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L102)).

The same `snGrad` operator, on the same faces, with the same `rAUf` weight, acts on `alpha` in the force and on `p_rgh` in the pressure flux. When `kappa_f` is one constant, `p_rgh = sigma kappa alpha + C` solves the Poisson equation face by face, the flux difference is zero, and the velocity stays zero ([SL deck 4/9](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/9)). What survives the projection is the face residual `R_f = phig - p_rghEqn.flux()`, which `fvc::reconstruct` turns into the parasitic velocity ([SL article `eq:fluxresidual`](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2213-L2217)).

## Why it matters

The pressure is a scalar cell field, so the projection removes only the part of `phig` in the range of the discrete gradient. On a uniform mesh with linear face interpolation the product rule is exact: `sigma <kappa>_f snGrad(alpha) = sigma snGrad(kappa alpha) - sigma <alpha>_f snGrad(kappa)`. The first term is absorbed identically. The driver of every parasitic current is the second term, the alpha-weighted face-normal gradient of `kappa` ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2218-L2228)). Two consequences follow. The driver differentiates `kappa`, so a second-order `kappa_f` gives a first-order force error. Only the tangential part of `snGrad(kappa)` is physical, because the exact `kappa` is constant along interface normals inside the force support. In continuum form, `curl(f_sigma) = sigma grad_t(kappa) x n`, and no pressure field can balance a curl ([METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L263-L272)).

The balance is a cancellation problem: the capillary predictor velocity is about 1.7 m/s at N = 128, and the 0.02 to 0.15 percent that the projection does not absorb is 2.5e-3 m/s ([CLAUDE](https://github.com/leia-openfoam/leia/blob/8867581/CLAUDE.md#L731-L750), [PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L578-L583)). The linear solver is not the limit on orthogonal hexahedra; see [[concepts/pressure-projection-and-linear-solvers]] for the non-orthogonal caveat.

## Where in the code

- The flux assembly: `surfaceTensionForce::integratedCSFFlux` returns `sigma kappaFace snGrad(forceWeight) magSf` ([`surfaceTensionForce.C`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.C#L145), declared in [`surfaceTensionForce.H`](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/surfaceTensionForce/surfaceTensionForce.H#L95); the code is shown on [SL deck 4/13](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/13)).
- The momentum equation evaluates `G_sigma` once and keeps it out of the matrix ([`UEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/UEqn.H#L8)).
- The flux-space residual on the active faces (`|snGrad(alpha)| > 0`) is written per final corrector to `capillaryFluxResidual.csv` ([`pEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/pEqn.H#L170-L215)).
- The pressure-velocity coupling was diffed line by line against interFoam: same `UEqn`, same `phig`, same `p_rghEqn.flux()`, same velocity update ([STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L968-L977)).

## Evidence

| claim | number | where |
|---|---|---|
| a constant curvature on a frozen circle gives no current | max abs U about 3e-11 m/s | [METHOD 5](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L259-L262), MEASURED |
| the CSF product form and the potential form agree on uniform meshes | 3.770e-9 / 7.635e-9 / 1.039e-8 m/s at N = 32 / 64 / 128; field difference 2.02e-14 to 1.50e-13 | [RM frozen circle](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1046-L1057), MEASURED |
| a 10 percent perturbed mesh breaks the balance with exact curvature | 8.584e-6 / 7.328e-4 / 1.392e-3 m/s at N = 32 / 64 / 128 | [RM frozen circle](https://github.com/leia-openfoam/leia/blob/8867581/docs/capillary-level-set-research-roadmap.md#L1046-L1057), MEASURED |
| the exact curvature removes the step-1 kick | 2.15e-3 to 1.69e-9, factor 1.27e6 | [STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L97-L102), MEASURED |
| the projection absorbs almost the whole reconstructed flux | 99.85 to 99.98 percent; R_f/phig 1.49e-3 (arithmetic) and 1.80e-4 (foot point) at t = 0, N = 128 | [PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L561-L583), MEASURED |
| the projection is converged | 30x tighter solve (2.9e-10 to 8.6e-12) moves R_f, U_res and max abs U by at most 2.4e-4 relative | [PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L569-L577), MEASURED |
| fvc::reconstruct is proportional, not amplifying | velocity-residual fraction 1.4e-3 against flux-residual fraction 1.5e-3 | [PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L578-L583), MEASURED |
| hanging nodes do not break the balance | Laplace jump 145.470 Pa in every arm (exact 145.48); max over time L2 velocity 3.8e-10 to 8.2e-10 m/s | [STATUS 4](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L1207-L1219), [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L1388-L1396), MEASURED |
| the remainder is the alpha-weighted gradient of kappa | product rule exact on a uniform mesh with w = 1/2 | [SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2218-L2228), DERIVED |
| the non-gradient content converges only for the cell-centre inverse | order +2.01 against +0.09 for every other cell curvature; 3200x smaller at N = 512 | [`cellCentreInverseCurvature.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/cellCentreInverseCurvature.H#L36-L48), MEASURED |
| the residual does net work on the stationary droplet | P/E = +344 1/s against viscous -256 1/s at R/h = 25, filter off | [SL deck 4/9](https://leia-openfoam.github.io/leia/decks/quadratic-semi-lagrangian-level-set.html#/4/9), MEASURED |

## Decisions

- The force stays a face flux; no model returns a cell vector. The `directCell` path of the conormal model was removed ([[models/surface-tension-force]]).
- The flux-space residual, not `max|U|` alone, is the instrument for the force balance ([PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L593-L602)).
- Production GAMG at tolerance 1e-9 and relTol 0.01 is kept for orthogonal droplet studies ([PCS 8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L600-L602), [[concepts/pressure-projection-and-linear-solvers]]).

## Open questions

1. The skew-mesh residual with exact curvature (1.392e-3 m/s at N = 128 on a 10 percent perturbed mesh) is insensitive to every lever except the `fvc::reconstruct` of the velocity correction ([METHOD 9.7](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L794-L799)).
2. A formulation in which the jump enters the pressure operator directly, or a delivery with one curvature per interface element, is the structural route the article names ([SL article](https://github.com/leia-openfoam/leia/blob/8867581/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2265-L2271), [[concepts/variational-capillary-force]]).

## Related

[[hubs/surface-tension]], [[models/surface-tension-force]], [[concepts/curvature-from-the-fit]], [[concepts/cell-centre-inverse-curvature]], [[concepts/parasitic-current-mechanism]], [[concepts/pressure-projection-and-linear-solvers]], [[concepts/well-balanced-exact-curvature-gate]], [[concepts/variational-capillary-force]], [[concepts/static-local-refinement]], [[decisions/surface-tension-reconstructed-curvature]], [[cases/stationary-droplet]].

## Log

### 2026-09-28
Created from METHOD 5, PCS section 8, the SL article sections on surface tension and the flux residual, the frozen-circle and hanging-node gates.
