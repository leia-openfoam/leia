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
- `config/davof/sphereNormal3D.yaml`: static gate 1, the normal and the total
  area of a sphere on the uniform ladder N = 20/40/80/160 from three states
  (`alphaSource exactSphere`: exact face and cell fractions; `linearInterpolant`:
  one piecewise-linear surface, the identity exact for it; `planePhaseIndicator`:
  leia's plane-based production state), with OpenFOAM's geometricVoF plicRDF /
  gradAlpha / isoAlpha run on the identical alpha of each state as the
  cross-check. The pre-registered predictions and what was measured are in the
  config header. The matching TwoPhaseFlow benchmark is
  `run/benchmark/reconstruction/sphereNormal3D` in the sibling TwoPhaseFlow
  clone (branch `davof`).
