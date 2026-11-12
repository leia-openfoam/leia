---
title: "The Fourier symbol of the semi-Lagrangian quadratic fit"
description: "The interior update on a perfect cube lattice has rho = 1 exactly for Co <= 1, the measured 1.00441 sits on boundary cells, a boundary cell amplifies exactly when its departure point lies outside its data, a perturbed stencil amplifies, and on a one-ring stencil the quadratic is the only degree that works (2026-10-13 to 2026-11-12)."
aliases: [amplification symbol, growth factor of the fit, von Neumann analysis of the SL fit, boundary extrapolation law]
kind: concept
status: open
part: advection
tags: [concept, part/advection]
date: 2026-11-12
date_settled:
decided_by: [workflow/scripts/sl_amplification_symbol.py, docs/analysis-amplification-sl/ms.tex]
code: [workflow/scripts/sl_amplification_symbol.py, docs/analysis-amplification-sl/sympy/boundary_lambda.py, docs/analysis-amplification-sl/sympy/hexagonal_lattice.py, docs/analysis-amplification-sl/sympy/perturbed_stencil.py, docs/analysis-amplification-sl/sympy/corner_face_data.py, docs/analysis-amplification-sl/sympy/degree_on_irregular_stencil.py, src/leiaLevelSet/semiLagrangian/uncachedQuadraticWeightedLeastSquaresReconstruction.C]
sources: ["docs/analysis-amplification-sl/ms.tex", "docs/analysis-amplification-sl/data/tables/sl_amplification_symbol.csv", "docs/analysis-amplification-sl/review-log.md"]
---
# The Fourier symbol of the semi-Lagrangian quadratic fit

> Verdict (2026-11-12). On a uniform lattice the semi-Lagrangian update is a constant-coefficient operator, so its growth factor is available in closed form: `g_j(c) = b(c)^T M_1^{-1} b(s_j) / |s_j|^2` contains no `h`, which derives the measured `h`-invariance of `rho - 1` ([[concepts/polyhedral-fit-amplification]]). In the interior of a perfect cube lattice `rho = 1` exactly for `Co <= 1` along an axis, so the measured whole-mesh `rho = 1.00441` is NOT the interior scheme; the solver's own field confirms that all 484 344 interior cells carry one value of `Lambda` (1.01639260, nine significant figures against the symbol) and the maximum 1.0527251 sits on a boundary cell. A boundary cell amplifies exactly when its departure point lies outside its data: `max|G| = 1 + 2 Co` with no boundary datum, and `1 + (4/3) Co (1/2 - Co)` with the production `zeroGradient` face datum, which freezes the cell at `Co = 1/2`. A perturbed stencil amplifies in the interior too. The reference is `docs/analysis-amplification-sl/ms.tex`.

## What it is

The anchored quadratic least-squares fit, read at the departure point, makes the new cell value a fixed weighted sum of the old values of its stencil. On a lattice that repeats, every cell applies the same weights, so a Fourier wave `exp(i theta . m)` comes out as the same wave times `G(theta) = (1 - sum_j g_j) + sum_j g_j exp(i theta . s_j)`, and `rho = max_theta |G|`. The exact solution multiplies every wave by a pure phase, so `|G| = 1` is fidelity to the material derivative, `|G| < 1` is numerical dissipation and `|G| > 1` is growth. The lattice must repeat; it need not be Cartesian (the honeycomb works with its own integer labels).

## Why it matters

It replaces a 2000-step power iteration on half a million cells with algebra, and it says WHERE growth comes from, which a single whole-mesh number cannot. It moved the search for the amplifier from the interior scheme to the boundary treatment and to mesh irregularity.

## Where in the code

`workflow/scripts/sl_amplification_symbol.py` builds the symbol (exact rationals, a 1D gate that must reproduce Lax-Wendroff symbolically, and a 96-row Courant sweep). The scripts under `docs/analysis-amplification-sl/sympy/` cover the boundary cell, the honeycomb, the perturbed stencil, corner and face data, and the polynomial degree. The operator modelled is `uncachedQuadraticWeightedLeastSquaresReconstruction`: anchored on `psi_i - psi_c`, constant-free basis, weights `1/|x_i - x_c|`.

## Evidence

| claim | number | where |
|---|---|---|
| 1D fit is Lax-Wendroff | `Lambda = 1 + C - C^2`, `abs(G)^2 = 1 - 4C^2(1-C^2)sin^4(theta/2)`, symbolic | `sl_amplification_symbol.py validate-1d` |
| interior, perfect cube, does not amplify | `rho = 1.00000000` in 96 configurations, Co 0.0164 to 0.95 | `data/tables/sl_amplification_symbol.csv` |
| the limit is the one-ring reach | first `rho > 1` at Co 1.000000 along an axis, 1.071207 body diagonal | ms.tex section on the Courant limit |
| the symbol is the code | `Lambda` 1.01639260 symbol against 1.01639260 code, 484 344 cells, std 2.2e-16 | ms.tex, check against the running code |
| growth is at the boundary | global max 1.0527251 on a boundary-adjacent cell | same |
| boundary law, no datum | `max abs(G) = 1 + 2 Co`, exact to 4e-16 | `sympy/boundary_lambda.py` |
| boundary law, production datum | `1 + (4/3) Co (1/2 - Co)`, peak 13/12 at Co 1/4, exact to 4e-16 | ms.tex Courant curve |
| inflow freezes at Co 1/2 | neighbour weights sum to 0.000000; `psi_c^{n+1} = psi_c^n` | ms.tex frozen inflow |
| face datum: presence, not weight | `max abs(G)` 1.083333 for weight 1 down to 1e-8; jumps to 1.5 at 0 | ms.tex |
| honeycomb limit | Co 0.866025404 = sqrt(3)/2 between neighbours | `sympy/hexagonal_lattice.py` |
| perturbed stencil amplifies | 10 % displacement: 1.001059 at Co 0.0164, 1.026140 at Co 0.3 | `sympy/perturbed_stencil.py` |

## Why it failed, or why we think so

Routes judged with the symbol and found worse than the current fit; each has a line in [[negative-log]].

- Moving the weights to the foot point: no change in `rho` anywhere, plus a pole where the foot meets a stencil point.
- A second interior layer at the wall: extrapolation better determined, so worse (2.125 against 1.083 at Co 0.25).
- Optimising the weights at one displacement: reaches `max|G| = 1` there and is worse elsewhere (1.015 against 1.001 along y).
- Reconstructing from the cell's own corner and face values: the perfect interior then amplifies at every Courant number (1.2667 at Co 0.5) and loses quadratic exactness. Interpolated data add no information: low order reads its O(h^2) error as curvature, high order from the same fit returns the current scheme.
- Changing the degree on the one-ring stencil: linear is forward-time centred-space, unstable at every Courant number (`abs(G)^2 = 1 + C^2 sin^2 theta`); cubic is singular on the perfect lattice because `x^3 = x` on {-1,0,1}, and ill-conditioned (1.3e3) and amplifying (2 to 6) when perturbed.

## Decisions

None taken. The analysis proposes, for the supervisor: an eikonal-consistent boundary face value `psi_c + (x_f - x_c).n` together with `stencilBoundaryFaces inflowOnly`, and replacing `Lambda` as a screen by a per-cell growth factor.

## Open questions

- The boundary-cell figures are indicators, not proofs: the direct test is the power iteration with boundary rows replaced by a non-amplifying rule.
- Modelled boundary stencils reach 1.033 at Co 0.0164; the measured maximum is 1.0527 (inlet, outlet and wall patches).
- Whether snapping explains the measured 36 % `rho - 1` penalty of a snapped hex mesh over blockMesh: the perturbation sizes match, the snapped mesh itself is not yet measured.
- The minimax weight design per cell over the displacements it sees: well posed, not attempted.
- Whether the frozen inflow at Co 0.5 moves the convergence orders of the 2D hex advection rung.

## Related

[[concepts/polyhedral-fit-amplification]] measured `Lambda` and `rho` on real meshes; this note explains them. [[negative-log]] holds the failed routes.

## Log

### 2026-11-12
Created, covering the analysis of 2026-10-13 to 2026-11-12 (`docs/analysis-amplification-sl/ms.tex`, three review passes in `review-log.md`). Added the corner and face data and the degree scan on this date.

### 2026-11-12 (second entry)
CANDIDATE, from the supervisor: vertex and face-centre values as independent unknowns, set exactly at t = 0 and advected by the same semi-Lagrangian step, each departure point interpolated in its host cell. Matrix symbol, uniform cubes: with the complete tensor point set (3x3 in 2D; 3x3x3 in 3D, which needs the 12 edge midpoints, 8 unknowns per cell) the growth factor is 1.000000 at every Courant number tested, 0.0164 to 1.25, in four flow directions; the current scheme has 2.125 at Co 1.25. The one-step error of a long wave at Co 0.3 is 2.8e-5 against 3.6e-4. In 3D, centre + faces + vertices WITHOUT edge midpoints amplifies with both spaces tried: least-squares quadratic 1.278 at Co 0.0164, 15-term interpolant 1.908 at Co 1.0. Open: cost, polyhedral cells, the departure point search in parallel. Script `docs/analysis-amplification-sl/sympy/multipoint_unknowns.py`; ms.tex section "Vertex and face values as their own unknowns".

