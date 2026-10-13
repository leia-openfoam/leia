# davof: the dual area/volume-of-fluid method

The theme of the DAVOF studies (`config/davof/*.yaml`, `cases/davof/*`,
`studies/davof/`). Its results land in `davof-article/data/{figures,tables}`
through `workflow/scripts/make_davof_normal_table.py`; there is no deck yet.

DAVOF transports, next to the cell volume fraction `alpha_c`, the liquid area
fraction `alpha_f` of every face, and recovers the interface area normal of a
cell locally from the Gauss identity on the liquid sub-volume,

    m_c = - sum_f alpha_f S_f^out ,

with `|m_c|` the interface area in the cell and `m_c/|m_c|` the normal out of
the liquid. Library: `src/leiaLevelSet/davof` (`libleiaDavof`). Gates:

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
