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
