---
title: "Literature on the stability of semi-Lagrangian point schemes"
description: "Online survey of 2026-11-12: the 1D result is classical (Strang 1962), the node-based tensor quadratic is computed stable (Bradley 2024), octree level sets store phi at vertices and give hanging nodes a ghost value, and no analysis of a cell-centred quadratic scheme at a 2:1 transition was found; the one report is negative."
aliases: [literature survey of semi-Lagrangian stability, what is known about SL point schemes]
kind: concept
status: open
part: advection
tags: [concept, part/advection]
date: 2026-11-12
date_settled:
decided_by: [docs/analysis-amplification-sl/literature-survey.md]
code: [docs/analysis-amplification-sl/sympy/scheme_matrix.py]
sources: ["docs/analysis-amplification-sl/literature-survey.md"]
---
# Literature on the stability of semi-Lagrangian point schemes

> Verdict (2026-11-12). We are not the first on the uniform lattice, and we found nobody before us on the 2:1 transition. The full survey, with every claim marked VERIFIED (read in the source in that session) or UNVERIFIED, is `docs/analysis-amplification-sl/literature-survey.md`. "Not found" means not found in one session of searches.

## What it is

A survey of eight topics, made to place the results of [[concepts/sl-amplification-symbol]]: Lagrange-interpolation stability, least-squares reconstruction on unstructured meshes, schemes with point unknowns, body-centred cubic data, octree level sets, Basilisk and Gerris, inflow boundaries, and matrix analysis on perturbed meshes.

## Why it matters

It says which of our results to cite and which to prove ourselves, and it shows how the octree level-set codes avoid the problem that our 2:1 analysis measures.

## Evidence

| our result | verdict | source |
|---|---|---|
| quadratic fit, growth factor 1 for Co <= 1 | known in 1D only: centred even degree has sup of the symbol equal to 1 on [-1, 1] | Strang 1962; Besse and Mehrenberger, Math. Comp. 77 (2008), Corollary 1 |
| the same in 2D and 3D for a least-squares fit | not found | |
| linear fit unstable, non-centred reconstruction unstable | known | Besse and Mehrenberger 2008 |
| boundary cell, 1 + 2 Co | known only as a statement that an inflow boundary needs a datum | Trefethen, Math. Comp. 45 (1985) |
| honeycomb limit sqrt(3)/2 | not found | |
| weak growth on perturbed meshes | known for node-based tensor schemes: the excess is small, at isolated shifts, of order delta^4 | Bradley, J. Comput. Phys. (2024), DOI 10.1016/j.jcp.2024.113034 |
| tensor quadratic point set, growth factor 1 | known as computed eigenvalues, no proof; the natural element interpolant is unstable only from degree 3 | Bradley 2024, citing Baptista 1987 and Oliveira and Baptista 1995 |
| point sets without edge midpoints, centre + vertices | not found | |
| 2:1 transition, hybrid schemes | not found | |

## What the octree level-set codes do

- Min and Gibou (J. Comput. Phys. 225, 2007) and Mirzadeh et al. (2016) store the level set at the cell VERTICES of a non-graded tree.
- A hanging node is not an unknown. It gets a ghost value: linear interpolation along the coarse edge with a transverse second-difference correction.
- The interpolation at the departure point is multilinear in the leaf plus a quadratic correction with the second derivative of minimum modulus, so it is limited.
- Min and Gibou wrote the full tensor quadratic on the 3 x 3 (x 3) nodes of the parent cell and rejected it as sensitive to kinks. This is a practical objection, not a linear instability.
- The finest cells lie on the interface, so the transition is away from it, and reinitialisation runs every few steps.
- Both papers state unconditional stability without an analysis.
- Basilisk gives every cell a uniform stencil at its own level through ghost cells (bilinear prolongation, averaging restriction).

## Why it failed, or why we think so

The one report on quadratic semi-Lagrangian interpolation at a graded octree is negative: Olshanskii, Terekhov and Vassilevski (Comput. Fluids 84, 2013) found it unstable for face-centred data and report spurious modes on fine-to-coarse interfaces. They do not localise the instability. Bradley gives a 1D example with cubic interpolation on nonuniform nodes whose largest eigenvalue exceeds 1 + 1e-3.

## Open questions

- Besse and Mehrenberger call the stability region on triangulations an open question. A proof for quadratic elements is announced in Ferretti and Mehrenberger (Math. Comp. 89, 2020) as forthcoming; no published proof was found.
- Berger (Math. Comp. 45, 1985) is summarised as showing that 1D Lax-Wendroff stays stable at a refinement interface. Unverified: the abstract could not be read.
- Ferretti and Perrone 2007 and Strang 1962 were not obtained in full.
- A Lax-Wendroff-type scheme from a least-squares polynomial on 5-, 9- and 13-point stencils is cited in arXiv:2210.05352; not checked.

## Related

[[concepts/sl-amplification-symbol]] holds our own results. [[hubs/advection]] is the hub.

## Log

### 2026-11-12
Created from the survey of the same day (one agent session, 93 searches and fetches; nine papers read in full or in the relevant sections).
