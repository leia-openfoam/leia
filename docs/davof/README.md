# davof: the dual area/volume-of-fluid method

The theme of the DAVOF studies (`config/davof/*.yaml`, `cases/davof/*`,
`studies/davof/`). Its results land in `davof-article/data/{figures,tables}`
through `workflow/scripts/make_davof_normal_table.py` and
`make_davof_advection_table.py`; there is no deck yet.

The port of the swept-face advection into leia (`src/davof`, `libleiaDavof`; the apps
`davofAdvectorFoam`, `davofInitFields`; `PORT-DESIGN.md`, `REPORT-swept.md` "The port into
leia") is built against OpenFOAM v2606's geometricVoF. Its one arithmetic difference from
the TwoPhaseFlow fork is the grazing-cut threshold of `cutCellPLIC::calcSubCell`
(`ROOTVSMALL` in v2606, `10*SMALL` in the fork), a round-off-level choice: the port
reproduces the fork bit for bit with the opt-in `grazingCut fork` (the `fvSolution` key of
the `davof` reconstruction scheme, the token `TPF_GRAZING_CUT`) in every test including the
Enright deformation on 64 ranks, and with the default the two runs separate at the first
grazing cut and the loop amplifies the one-ulp difference: the end norms agree to five
digits on every ladder rung, and on the chaotic Enright case the abort step becomes a
property of the realisation (the fork's N 32 run ends at T = 3, the default's aborts at
step 443 of 480; the numbers in PORT-DESIGN 1.2). Policy (Tomislav, 2026-11-07): the
default stays v2606's; a DAVOF row runs with `grazingCut fork` wherever a number is
reproduced against the fork's record, otherwise with the default; the key is in every case
stamp (`stampGrazingCut`) and in every table caption.

**The positioning rule (Tomislav, 2026-11-08): "We should never count on the closed form to
position PLIC, only always as input for CCS."** The explicit closed-form position
pPlane_c = (sum_f alpha_f (x_f - x_c).S_f^out - 3 alpha_c V_c)/|m_c| is the START of the
Consecutive Cubic Spline positioning and a diagnostic; no plane that enters a flux, a
normal model, a function object or a written field is the closed-form plane. The audit of
libleiaDavof and the apps (2026-11-08, port branch, no change needed), path by path:

| path | the plane's position comes from | closed form as the final plane? |
|---|---|---|
| `davof::positionPlane`, `positioning ccs` (the default) | CCS (`davofCcsPositioning.H`) on the cut volume along the model normal, from alpha, started at pPlane_c; a failed positioning leaves the cell WITHOUT a plane (normal zero), counted "end-of-interval cells (n without a plane)" | no |
| `davof::positionPlane`, `positioning vofCutCell` (opt-in) | the former `surfaceIteratorPLIC::vofCutCell` iteration (the cubic fit and the secant to isoFaceTol) | no |
| the explicit plane cut in `davof::reconstruct` (`cutCell_.calcSubCell(celli, pPlane_c, n)`) | the closed form | diagnostic only: `pPlane`, `planeVolDefect`, `planeAreaDefect` (written fields, the "explicit plane defects" line) |
| the swept scheme's cell surfaces (`buildCellSurfaces`, `surface plane`, `planeSource anchor`) | the frame from the used normal, the origin the CCS plane's cut centroid (the cell centre when the reconstruction gave no plane), the offset by the volume anchor (`fitAnchor volume`, `anchorMethod ccs`: CCS on the volume rule's V(c) to anchorTolRel); `anchorSkipConverged` keeps the CCS plane when its residual is within the tolerance; a failed anchor keeps the better bracket end (counted) | no |
| `planeSource exact` | the tangent plane of the exact sphere (a control) | no |
| the paraboloid (`surface paraboloid`) | the face fit (`fitParaboloidFull`) and the same volume anchor | no |
| the wisp surfaces (`sweptWispSurface matched`, the default) | the neighbour's anchored plane where its volume in the wisp is alpha V within the tolerance, else the neighbour's normal positioned by the wisp's alpha (`davof::solveOffset`, bisection on the volume rule); `anchored`: the latter always; `neighbour`: the neighbour's plane | no |
| `planeOf` (the swept clips, the face update, the consistency probe, the recompute of decision 2) | the anchored surface's plane (e3, o, c) | no |
| the normal models `exact`, `ringFit` | read the CCS-positioned planes of the first pass (the centroid as the evaluation point, the ring's centroids and normals for the fit); `ringArea`, `cellLocal`, `youngs`, `plicRDF`: normals only | no plane positioned by a model |
| `areaFractionError`, `reconstructionError` | the reconstruction's CCS planes (`surf.normal()`, `surf.centre()`); the exact reference from `davofInit` (the exact ball or cylinder: closed-form AREAS of the exact body, no PLIC) | no |
| `volumeFractionError` | the exact state, else `cutCellIso` (an alpha iso-value of the interpolated field, not a PLIC position) | no |
| `davofAdvectorFoam`, `davofInitFields` | the scheme above; the states (Detrixhe-Aslam, quadraticFaces, exactSphere, exactCylinder: no PLIC plane) | no |
| `davofStaticFaceFractions` (the static test) | the exact normal positioned by `vofCutCell` to positionTol | no (the former iteration, not the closed form) |
| `davofPlaneVolumeConsistency` (the diagnostic) | measures CCS from the explicit start and `vofCutCell`; forms no plane for anything | no |
| the diagnostics that write a plane | `pPlane` (the closed-form offset), `positioningResidual` (the CCS residual), `sweptPieceDump` (the anchored plane of a piece) | the closed form appears as the offset field only |

Outside the four apps: `leiaTestDavofNormal` (libleiaDavofState, the static gates) carries the
option `regenOrigin volume | explicit` (default volume); `explicit` would regenerate from the
closed-form position and is kept only as the gate's control.

Two document folders live here (2026-10-26):

- `davof-article/`: the preprint (`ms.tex`, its table bodies and figures,
  `make_table_advection.py`) and, in `data/`, the exported tables and figures
  of every study, the only tracked part. The preprint holds results only.
- `davof-discussion/`: the discussion document for the decisions on the
  method's next steps (`discussion.tex`, built with `latexmk -pdf`): the state
  of the implementation, every measured result with pgfplots charts of the
  study data, the predictions ledger, the open items and the options with
  their costs. Its data files under `data/` are written by
  `make_discussion_data.py` from `../davof-article/data/tables` (and, where
  no table exists, from the STATUS appendices named in that script); never
  edited by hand.

DAVOF transports, next to the cell volume fraction `alpha_c`, the liquid area
fraction `alpha_f` of every face, and recovers the interface area normal of a
cell locally from the Gauss identity on the liquid sub-volume,

    m_c = - sum_f alpha_f S_f^out ,

with `|m_c|` the interface area in the cell and `m_c/|m_c|` the normal out of
the liquid. Library: `src/leiaLevelSet/davof` (`libleiaDavofState`). Gates:

- `cases/davof/planeNormal3D/Allrun`: exact-solution gate (a plane; the
  normal, the explicit plane position, the tet/face consistency, the alpha
  against the exact plane cut and the realizability of the state by the
  recovered plane at round-off for both leia states, or the app exits
  non-zero).
- `config/davof/sphereNormal3D.yaml`: static gate 1, the normal, the total
  area and the explicit plane position of a sphere on the uniform ladder
  N = 20/40/80/160 from four states (`alphaSource exactSphere`: exact face and
  cell fractions; `quadraticFaces`: the quadratic interpolant of the exact
  signed distance cut exactly on every face triangle; `detrixheAslam` (alias
  `linearInterpolant`): one piecewise-linear surface, the identity exact for
  it; `planePhaseIndicator`: leia's plane-based production state), with
  OpenFOAM's geometricVoF plicRDF /
  gradAlpha / isoAlpha run on the identical alpha of each state as the
  cross-check. The pre-registered predictions and what was measured are in the
  config header.
- `config/davof/ellipsoidNormal3D.yaml`: static gate 2, the normal, the
  position and the CURVATURE of the DAVOF state on a triaxial ellipsoid
  (half-axes 1.0/0.8/0.6 mm, `signedDistanceEllipsoid`) on the same ladder,
  from the states `quadraticFaces` (q = 2), `detrixheAslam` and
  `planePhaseIndicator` (q = 1), with the curvature models `normalsOnly` (the
  quadric fit to the ring of centroids and normals in which the normals decide
  the second derivatives, the headline) and `hermite` (the naive fit), delivered
  to cells and faces with the parallel-surface closed form
  (`src/leiaLevelSet/davof/curvature`). Pre-registered prediction 5 and the
  measured orders in the config header. Nine GB at N = 160: run one case at a
  time (`--jobs 1`).
- `config/davof/planeNormal3DLadder.yaml`: the plane at N = 8/16/32 for the
  PLIC surface renderings (`figures/davof_plic_*.pdf`, DAVOF next to plicRDF
  and gradAlpha, written by `render_davof_interfaces.py` from the VTK
  surfaces of `libleiaDavofInterface`). `make davof-proposal DEST=...` copies
  the compact table and these figures into the DAVOF proposal. The matching TwoPhaseFlow benchmark is
  `run/benchmark/reconstruction/sphereNormal3D` in the sibling TwoPhaseFlow
  clone (branch `davof`).
- `config/davof/{plane,sphere}Translation3D_hex.yaml` (2026-10-13): the
  TRANSPORTED face fractions, driven through the fork's `advectorVoF`
  (TwoPhaseFlow branch `davof`: isoAdvection storing the end-of-step face
  areas, the `davof` reconstruction scheme with the normal models cellLocal /
  ringArea / ringFit / exact / plicRDF, `initDavofFields`, the abort guard).
  The plane of the static test translated along its normal and along x (the
  round-off test, N = 20..160, Co 0.05..0.5) and the sphere R = 0.15 carried
  over 3.3 diameters (N = 32/64/128, Co 0.1/0.2/0.5), seven rows each
  (`workflow/scripts/tpf_rows.py`), the predictions P-A1..3 and P-B1..3 in the
  config headers; `make_davof_advection_table.py` writes
  `tables/davof_advection_*` and `figures/davof_advection_*`
  (`workflow/README.md`, section "DAVOF: the TwoPhaseFlow advection studies");
  `make davof-preprint` copies them into `davof-article/`. Laptop smoke
  (`*_hexSmoke.yaml`, `tables/davof_advection_*Smoke*`, 2026-10-13, the same
  binary and dictionaries as the fork's laptop runs of that day): the plane
  at N = 20, Co 0.2, 20 steps: `davof-exact` and `davof-ringArea` at
  round-off on both arms (max L1(alphaf) below 1.2e-14, max L1(n) below
  1.3e-14), `davof-cellLocal` at round-off on the normal arm and at max
  L1(n) 1.1e-12 on the oblique arm (10 % above the 1e-12 gate over 20
  steps); the sphere at N = 32, Co 0.2, 160 steps (end values L1(alphaf) /
  L1(n) / L1(x) / Eshape, min alpha over the run): `davof-ringArea` 1.78e-2
  / 4.06e-2 / 1.35e-3 / 3.85e-4, -3.4e-21, 424 interface cells;
  `davof-exact` 1.04e-2 / 3.64e-3 / 9.83e-4 / 2.76e-4, -8.7e-15;
  `davof-plicRDFnormal` 2.07e-2 / 4.03e-2 / 1.58e-3 / 4.31e-4;
  `davof-ringFit` 7.93e-2 / 2.05e-1 / 1.32e-2 / 1.15e-3 (672 cells at the
  end, the fit fallbacks rising); `plicRDF` 3.89e-2 / 1.53e-3 / 4.32e-4 and
  `isoAlpha` 5.66e-2 / 1.60e-3 / 4.67e-4 (their initial alpha is the DAVOF
  state's, so not digit-identical to the fork's setAlphaField runs);
  `davof-cellLocal` ABORTED at step 111 (interface cells 4171 of initial
  416; L1(alphaf) 2.0e-1 there), recorded as `status diverged`. The ringArea
  and exact values reproduce the fork's laptop table to every printed digit.

## The ring normal's stencil: the opt-in keys of `normalModel ringArea` (2026-11-08)

The production normal of the swept-face advection is the one-ring area normal,
`d_k = sum_{j in R(k)} m_j`, `n_k = d_k/|d_k|`, the sum over the point-neighbour
ring of the cell (`src/davof/normalModels/ringAreaNormal.H`). Its stencil
variants are keys of the `ringArea` model in the davof dictionary of
`fvSolution`; every key is opt-in, and at the defaults the model runs the
production path bit for bit (the production variant of every study before this
date renders no key at all). In the studies the one word `TPF_NORMAL_STENCIL`
(`workflow/scripts/tpf_rows.py`, `NORMAL_STENCILS`) renders the keys as the
block `TPF_NORMAL_STENCIL_KEYS` after the `noNormalFallback` line of the four
swept-face templates (plane and sphere translation, rotation, deformation);
`production` renders an empty block.

| key | values | default | variant word |
|---|---|---|---|
| `ringStencil` | `one`, `two` (the union of the rings of the ring's cells, each cell once; the rings of ring cells on other processors are requested from their owners) | `one` | `twoRings` |
| `ringKernel`, `ringKernelSupport` | `none`, `gaussian` (`exp(-q^2)`, cut at the support), `wendland` (the C2 function `(1 - q/s)^4 (1 + 4 q/s)`); `q = |x_j - x_k|/h_k` with the cell centres, `h_k = V_k^(1/3)`, `s` the support in `h_k` | `none`, `2` | `kernelGaussian`, `kernelWendland` |
| `ringAngleFilter` | the angle in degrees beyond which a ring cell's own normal is dropped from the sum (0: off) | `0` | `angle30`, `angle60` |
| `ringPlaneFit`, `ringPlaneFitTol` | `none`, `mean`, `own`: the normal of the weighted least-squares plane through the ring's patch centroids (the positioned planes' centroids, two positioning passes), the eigenvector of the smallest eigenvalue of the weighted covariance with the area weights, the sign from `d_k`; `mean` about the weighted centroid, `own` about the cell's own centroid with the weights corrected so that `sum_j w'_j (x_j - x_k) = 0`; the tolerance of the eigenvalue tests (collinear or indefinite: `d_k`, a fit fallback) | `none`, `1e-10` | `planeFitMean`, `planeFitOwn` |
| `normalModel ringFit` | the ring quadric of the preprint (`weight area; alignTol 0.2; minCells 6;`) | - | `ringQuadric` |

Citations: the Gaussian kernel R. A. Gingold, J. J. Monaghan, Mon. Not. R.
Astron. Soc. 181 (1977) 375-389, doi:10.1093/mnras/181.3.375; the Wendland
function H. Wendland, Adv. Comput. Math. 4 (1995) 389-396,
doi:10.1007/BF02123482; the principal-component plane K. Pearson, Philos. Mag.
2 (1901) 559-572, doi:10.1080/14786440109462720, and H. Hoppe, T. DeRose, T.
Duchamp, J. McDonald, W. Stuetzle, Computer Graphics 26 (1992) 71-78,
doi:10.1145/133994.134011; the moment-corrected weights W. K. Liu, S. Jun, Y. F.
Zhang, Int. J. Numer. Methods Fluids 20 (1995) 1081-1106,
doi:10.1002/fld.1650200824. The studies: `config/davof/stencil_*.yaml`
(`make_stencil_configs.py`, the predictions in the headers), the tables
`davof_stencil_tables.py`, the perturbation of the face fractions
`davof_perturb_alphaf.py`; the record in `REPORT-swept.md`, "The ring normal's
stencil (2026-11-08)".
