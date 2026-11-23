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

### 2026-11-12 (third entry)
Supervisor's decision: Courant numbers above 1.0 are out of scope. The candidate of the second entry is therefore judged for Co <= 1 only. In that range the current scheme also has growth factor 1 in the interior of a perfect lattice, so the candidate's gain is not the interior growth factor. What remains: the 13 times smaller one-step error (2.8e-5 against 3.6e-4 at Co 0.3) and the complete point set of a wall cell. The Co 1.25 column is removed from ms.tex.

### 2026-11-12 (fourth entry)
The point-unknown candidate on the honeycomb, Co <= 1, worst of three flow directions (`sympy/multipoint_hexagon.py`). No combination of centre (C), edge midpoints (E) and vertices (V) has growth factor exactly 1. Best: C+V with the interpolating cubic, 3 unknowns per cell, 1.0002 to 1.0015, one-step error 2.8e-5 against 2.6e-4 for the current scheme. C+E: up to 1.092 (the edge midpoints do not enclose a departure point near a vertex). Interpolating degree 4 on V+E and C+V+E: up to 1.667 and 1.75. Least-squares quadratic on C+V+E: 1.000 to 1.032. The current scheme is exactly 1 up to Co 0.82 and 1.667 at Co 1.0. V alone and E alone: rank 5 of 6 for a quadratic.

### 2026-11-12 (fifth entry)
The point-unknown candidate against the original scheme on a CONSISTENT perturbed quadrilateral mesh (periodic block of 3 x 3 cells, vertices displaced by 10 % of h, worst of 2 seeds and 3 flow directions; `sympy/multipoint_perturbed.py`, `data/tables/multipoint_perturbed.csv`). Growth factor - 1 at Co 0.0164 / 0.3 / 0.5: original 3.9e-4 / 1.1e-4 / 1.3e-4; tensor point set 4e-8 / 6e-7 / 9e-7. So the earlier single-stencil number (1.026 at 10 % and Co 0.3, `sympy/perturbed_stencil.py`) overstates a consistent mesh by two orders: one perturbed stencil applied to every cell is not a mesh. A block of 2 cells per direction is degenerate (the cell centres stay a perfect lattice); 3 is the minimum. The 3D run of that script did not finish in 55 minutes on the laptop and gave no number; `sympy/scheme_matrix.py` replaces it (mesh `hexPert`, on Lichtenberg). That script also names the schemes (OWN-NEIGHBOUR/SPACE, the original scheme is `C-Cp/L2`) and has four validity gates, which pass.

### 2026-11-12 (sixth entry)
The constrained hanging-point variant of the scheme matrix, suffix `+c` (`sympy/scheme_matrix.py --variant constrained`). Every cell is read as an unrefined cell. A point that lies on the boundary of a cell and is not a regular point of it carries no unknown; its value is the interpolant of that cell, read at the point. This is the ghost value of the octree level-set codes ([[concepts/sl-point-schemes-literature]]). Gates G5 (no change on a mesh without hanging points) and G6 (`SC/I+c` does not amplify at the transition) pass, and the tables of the free variant are byte-identical before and after the change. MEASURED, 2D, laptop (`data/tables/scheme_matrix_quad21_constrained.csv`, `scheme_matrix_quad21patch_constrained.csv`): `BFC/I+c` (tensor set) and `FCC/I+c` have growth factor 1.000000 at every Co <= 1, in five flow directions, on the planar transition and on the refined patch. Without the constraint `BFC/I` has 1.0337 at Co 0.3 and 1.0501 at Co 0.5 on the planar transition. Prediction P2 holds in 2D. Prediction P3 holds in 2D: the least-squares schemes do not gain (`BFC-Cf/L2+c` 1.2e-3 against 1.2e-3 at Co 0.1; `BCC-Cf/L2+c` 3.6e-3 against 3.3e-3). The interpolating schemes with neighbour centres stay unstable (`BCC-Cf/I+c` 4.2, `BFC-Cp/I+c` 4.0 at Co 0.3). 3D: one check only, `T/I+c` and `SC/I+c` are 1.000000 at Co 0.3 on the planar hexahedral transition; the full 3D tables are Lichtenberg job array 56255724 and have not landed.

### 2026-11-12 (seventh entry)
The ghost-cell variant for the cell-centred schemes, suffix `+g` (`sympy/scheme_matrix.py --variant ghost`). Every cell reads a uniform stencil at its own cell size, so it applies the weights of the perfect lattice. Restriction is the mean of the 2^d fine values; prolongation is the reconstruction of the coarse cell, read at the ghost position. Both are explicit. The mean is not the centre value of a quadratic, so `+g` is not exact for a quadratic in the cells at the transition. Gate G7 (no change on equal cells, square and cube) passes. MEASURED, laptop, growth factor - 1, largest of five flow directions, at Co 0.0164 / 0.1 / 0.3 / 0.5 / 1.0 (`data/tables/scheme_matrix_{quad21,quad21patch,hex21,hex21patch}_ghost.csv`). `C-Cp/L2+g`, planar quadrilateral transition: 4.3e-4 / 1.9e-3 / 2.9e-3 / 3.2e-3 / 3.3e-3, against 3.3e-3 / 2.0e-2 / 5.8e-2 / 9.7e-2 / 2.1e-1 without ghost cells. Planar hexahedral transition: 4.6e-4 / 2.2e-3 / 3.9e-3 / 4.4e-3 / 3.7e-3. Prediction P4 holds in 2D and in 3D (the criterion is 4.41e-3; the 3D value at Co 0.5 is 4.38e-3, so the margin is small). The growth is not zero: the worst direction is tangential to the transition. On the refined patches it is at most 2.2e-4 (2D) and 1.6e-4 (3D, only 3 waves per direction of the block). `C-Cp/I+g`: 6.5e-4 / 1.7e-3 / 0 / 0 / 0 (2D) and 4.8e-4 / 1.5e-3 / 0 / 0 / 0 (3D). `C-Cf/I+g` amplifies (2.2e-2 in 2D and 1.1e-1 in 3D at Co 0.5); `C-Cf/L2+g` is rank-deficient on cubes. Also on this date: the scheme names are in the SL pre-print (section "Names for semi-Lagrangian point schemes") and in the SL deck.

### 2026-11-12 (eighth entry)
The scheme matrix of the free variant landed from Lichtenberg (job array 56254576, ten meshes; seven tables also computed on the laptop, equal to round-off). Tables `data/tables/scheme_matrix_<mesh>.csv`; verdicts by the pre-registered criteria in `scheme_matrix_summary.tex` (`sympy/scheme_matrix_tables.py`); report section "The scheme matrix: every named scheme on every mesh" in `ms.tex`. MEASURED: (1) only `C-Cp/*` and the tensor set have growth factor 1 on the square and the cube; (2) on perturbed hexahedra (5 to 20 %) only `T/I` has less growth than `C-Cp/L2` (at most 1.7e-5 against 6.0e-4 for Co <= 0.5), which closes the 3D perturbed question; (3) `C-Cp/L2` fails the planar 2:1 transition in 2D and 3D with nearly equal numbers; (4) with hanging points as unknowns only `T-Cf/L2` passes the planar hexahedral transition. Five lines in [[negative-log]]. The `+c` tables in 3D failed on the first submission (the batch Python had no numpy) and run again as job array 56258883. Candidates for the supervisor's decision: `C-Cp/L2+g`, the tensor set with `+c`, `T-Cf/L2`.

### 2026-11-12 (ninth entry)
The figure of the scheme names is a TikZ drawing of the data stencils (`docs/analysis-amplification-sl/figures/scheme_stencils.tex`, `make scheme-stencils-figure`), in `ms.tex`, in the SL pre-print and in the SL deck. It replaces a UML class diagram of the same day, which the supervisor rejected: the figure must show each data stencil. Supervisor's decision: try `C-Cp/L2+g` first ([[decision-log]]).

### 2026-11-12 (tenth entry)
The `+c` tables in 3D landed (Lichtenberg, job array 56258883; the two 2D tables are byte-identical to the laptop tables). MEASURED: `T/I+c` has growth factor 1.000000 at every Co <= 1 in five flow directions on the planar hexahedral 2:1 transition, and at most 1 on the refined patch. Prediction P2 holds in 3D. The sets without edge midpoints do not gain (`BFC/I+c` 7.0e-2 at Co 0.5); `T-Cf/L2+c` 2.7e-3 at most, equal to the free variant within 3e-4 (P3 holds). `data/tables/scheme_matrix_{hex21,hex21patch}_constrained.csv`.
- 2026-11-13 · the first solver test of the 2:1 transition result: on the refined-box translation the original scheme is destroyed at the transition and the ghost-cell class passes its pre-registered gate ([[concepts/sl-ghost-cell-reconstruction]], STATUS 11.35). The symbol and the solver agree in sense; the measured growth was not compared with the predicted 1.096 per step.

### 2026-11-23
Three cell sizes (two refinement levels above the base mesh) for the `+g` schemes. Script `docs/analysis-amplification-sl/sympy/scheme_matrix_three_levels.py` (gates G8 to G11, predictions P5 and P6, written and committed before the first three-level run); meshes `quad3lev_b<k>` and `hex3lev_b<k>` of `scheme_matrix.py`; report section "Two refinement levels" in `ms.tex`. The block: 3 cells of size 4, k cells of size 2, 4 cells of size 1, k cells of size 2, normal to the transitions; tangential period one cell of size 4; Courant number on the cell of size 1. Unknowns per block 19 + 4k (2D) and 67 + 8k (3D). `assemble_g` applies the ghost rules level by level: a prolongation reads only larger cells, so the dependency ends at the largest size (DERIVED from the code). All runs on the laptop.

Gates, all pass. G8: the changed code reproduces the `+g` tables of `quad21`, `quad21patch`, `hex21`, `hex21patch` in all 8 decimals. G9: the three-level constructor with equal sizes returns growth factor 1 for `C-Cp/L2+g` and `C-Cp/I+g`, and every `+g` scheme equals the scheme without the suffix. G10: for k = 2 and 3 no cell has a larger AND a smaller cell-point-cell neighbour, and the longest prolongation chain is 1; for k = 1 the cells of size 2 have both (4 cells in 2D, 8 in 3D) and the chain is 2. G11: a linear function is returned exactly on the six meshes (at most 3e-13).

MEASURED, `C-Cp/L2+g`, growth factor - 1, largest of five flow directions, at Co 0.0164 / 0.1 / 0.3 / 0.5 / 1.0, at the finest wave sampling (`data/tables/scheme_matrix_three_levels_waves_{2d,3d}.csv`; t = tangential, n = normal to the transitions).

| mesh | buffer | Co 0.0164 | Co 0.1 | Co 0.3 | Co 0.5 | Co 1.0 | worst direction |
|---|---|---|---|---|---|---|---|
| 2D | one transition | 4.3e-4 | 1.9e-3 | 3.0e-3 | 3.2e-3 | 3.3e-3 | t |
| 2D | k = 1 | 3.3e-4 | 5.3e-4 | 1.2e-4 | 3.3e-6 | 0 | t, (2 1), (1 2) |
| 2D | k = 2 | 3.9e-4 | 1.2e-3 | 1.3e-3 | 1.6e-3 | 1.7e-3 | t |
| 2D | k = 3 | 3.6e-4 | 1.6e-3 | 2.8e-3 | 3.3e-3 | 3.7e-3 | t |
| 3D | one transition | 4.6e-4 | 2.2e-3 | 3.9e-3 | 4.4e-3 | 3.7e-3 | t |
| 3D | k = 1 | 4.1e-4 | 1.4e-3 | 2.2e-3 | 1.2e-3 | 0 | n |
| 3D | k = 2 | 2.8e-4 | 9.0e-4 | 1.7e-3 | 2.1e-3 | 2.4e-3 | t |
| 3D | k = 3 | 3.6e-4 | 1.7e-3 | 3.6e-3 | 4.5e-3 | 5.2e-3 | t |

`C-Cp/I+g` on the three-level meshes: at most 1.4e-3 (2D, k = 2, Co 0.1) and 1.3e-3 (3D, k = 1, Co 0.1); 0 at Co 1.0. The tables `data/tables/scheme_matrix_{quad,hex}3lev_b{1,2,3}_ghost.csv` use the wave sampling of the other meshes (2 waves normal to the transitions on these long blocks). That sampling misses the worst wave in places: up to 33 % in 2D, and 0 against 2.4e-3 in 3D at Co 1.0, k = 2. The table values are lower bounds; the numbers above are the ones to quote.

Verdict on P5. P5a HOLDS: with k = 2 the growth factor - 1 is at most 1.7e-3 (2D) and 2.4e-3 (3D), below 4.41e-3 at every Courant number. P5b HOLDS: k = 2 is 0.41 to 0.91 times the one-transition value. P5c FALSIFIED: k = 1 is not worse than k = 2. In 2D it is smaller at every Courant number; in 3D it is larger at Co 0.0164 to 0.3 (ratio 1.27 to 1.57) and smaller at Co 0.5 and 1.0. P5d MISSED: k = 3 has 1.3 to 2.2 times the growth of k = 2 for Co >= 0.1; in 3D it is 4.5e-3 at Co 0.5 and 5.2e-3 at Co 1.0, above the pass threshold and below the fail threshold 1.03e-2, so the symbol does not decide k = 3 in 3D. By the pre-registered rule P5 is falsified (through P5c). No buffer width fails the criterion, so there is no line in [[negative-log]]. The design buffer k = 2 passes. The growth factor does not give a reason to prefer k = 2 over k = 1; the reason for k = 2 remains the stencil rule of G10.

MEASURED, not pre-registered: the growth factor depends on the widths of the strips of the periodic block, and not monotonically. 2D scan of 70 blocks (`data/tables/scheme_matrix_three_levels_scan_2d.csv`: one transition with 2 to 8 coarse and 4 to 16 fine cells; three levels with k = 1 to 10): the largest growth factor - 1 of a block is between 2.2e-4 and 3.7e-3; the largest is the block k = 3 above. The one-transition value 3.3e-3 of `quad21` (3 coarse, 4 fine cells) falls to 6.6e-4 with 6 fine cells. So every number of the 2:1 tables is the value of one block, not a bound over all widths. HYPOTHESIS: the interior of each strip is neutral, so a disturbance made at a transition does not decay inside a strip, and the growth of a block mode depends on how the two transitions of the block combine; this is not tested.

Power-boundedness (P6; Lax and Richtmyer; Ferretti and Mehrenberger). P_max = max over the Bloch waves and n <= 2000 of the 2-norm of A^n, columns `pmax_Co*`, `npmax_Co*`, `pmaxhalf_Co*` APPENDED to the `+g` tables of `quad21`, `hex21` and the six three-level meshes (no existing column changes). MEASURED: where the growth factor g is above 1, P_max is at n = 2000 and is g^2000 times 1.3 to 36 (largest: 3.5e4, 3D, k = 3, Co 0.5; k = 2: 99 in 2D and 307 in 3D). Where the table growth factor is at most 1, P_max at N = 2000 is at most 9.2 in 2D and at most 18.5 in 3D: P6a (P_max < 10) HOLDS in 2D and is MISSED in 3D. Nine such cases (4 in 2D, 5 in 3D) grew by more than 5 % between N = 1000 and N = 2000 or were above 10. Repeated with N = 16000 (`scheme_matrix_three_levels_powers_long.csv`), every one is constant from N = 8000 on; the largest is 14.9 (3D, k = 1, `C-Cp/I+g`, normal flow, Co 0.0164, reached at n = 6880) and 18.5 (3D, k = 1, `C-Cp/L2+g`, tangential flow, Co 0.1, reached at n = 1660). The norm is sampled at every n up to 200 and at every 10th n after that. P6 is not falsified in substance: no neutral case shows algebraic growth. The 2-norm is the plain norm over cells of different sizes; a volume-weighted norm was not tried.

Which cell's polynomial (a check of whether the tables represent the solver). DERIVED from the code: for a scheme with OWN = C the free variant (`assemble`: `j, s = (u, 0) if own == "C"`) and the `+g` variant (`assemble_g`: `stencil(j, x, x - c)`) evaluate the polynomial of the ARRIVAL cell at the departure point, also when the departure point leaves the arrival cell. This is what the solver does. Only the schemes with point unknowns locate the HOST cell. The header sentence of `scheme_matrix.py` already states this exception. No arrival-cell option was needed, and no second set of tables exists.

Open: the worst block over all strip widths in 3D (no 3D scan); the patch geometry (edges and corners) with three cell sizes; a comparison with the measured growth of [[concepts/sl-ghost-cell-reconstruction]] on a dynamic mesh.
