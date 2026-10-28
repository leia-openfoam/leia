#!/usr/bin/env python3
"""The rows and velocity arms of the TwoPhaseFlow advection studies
(cases/davof/{plane,sphere}Translation3D, solver advectorVoF of the fork
TwoPhaseFlow branch davof, run/benchmark/advection there).

A ROW names what the alpha solver dictionary of the fork's benchmark selects:
the reconstruction scheme and, for the davof scheme, its normal model
(davofNormalModel.H: normalModel, applyTo all | check, sliverTol, maxAngle,
noNormalFallback keep | ringArea). The study sweeps the single token
TPF_ROW; materialize._with_derived_tokens expands it through row_tokens()
into the six tokens the templates name (TPF_RECON_SCHEME, TPF_NORMAL_MODEL,
TPF_APPLY_TO, TPF_SLIVER_TOL, TPF_MAX_ANGLE, TPF_NO_NORMAL_FALLBACK), so
that a row is one word in a config and one column in case_params.json. The
comparator rows (plicRDF, isoAlpha) carry the davof keys inertly: the other
schemes ignore them, and every row starts from the same initDavofFields
state (the comparators ignore alphaf.water).

The row davof-plicRDFnormal-ring (2026-10-14) is the plicRDF-normal control
with noNormalFallback ringArea: the cells the davof scheme marks by their
face fractions but the helper plicRDF does not (alpha outside
(surfCellTol, 1 - surfCellTol); "kept m" in the fork's Info line, 16 to 24 of
1800 cells at N = 64 and 144 to 182 of 7100 at N = 128 in the hex study)
take the ring normal instead of the own m_c. Every other row keeps the
default keep, so its fvSolution gains one line and nothing else changes.

The seventh element (2026-10-15) is the fork's face-clipping model
(davofFaceClipping.H, keyword faceClipping of the alpha solver dictionary:
plane, today's plane clipping of the moved PLIC plane, or planeSagitta, the
plane clipping plus the second-order sagitta correction from the ring
curvature tensor; the sub-keys sagittaMinCos / sagittaMinChord /
sagittaMaxFraction / sagittaForm / tensorSource are written at their
defaults in the templates). The four sagitta rows: davof-sagitta-exact (the
exact normal: isolates the fractions' order), davof-sagitta-ringArea (the
ring normal everywhere), davof-sagitta-cellSliver (the OWN normal
everywhere except the slivers, |m_c| < 0.1 h^2, which take the ring normal:
normalModel ringArea, applyTo check, sliverTol 0.1, maxAngle 180) and
davof-sagitta-cellLocal (the own normal, undamped). Every existing row
carries plane, so its fvSolution gains the keyword block and nothing else
changes.

A VELOCITY ARM (token TPF_VELOCITY_ARM) names the uniform velocity, expanded
by velocity_tokens() into TPF_UX/UY/UZ: the plane test moves the plane along
its normal (`normal`, U = n, the sliver-amplifying arm) or obliquely
(`oblique`, U = (1 0 0)); the sphere is carried along x (`axial`).

The names follow the fork's advectionRows.py (davof-<normalModel>, the
comparators by scheme), without the state construction in the name: the
state is the study's TPF_DAVOF_METHOD (detrixheAslam on the plane,
quadraticFaces on the sphere), one per study.
"""
from collections import OrderedDict

# The eighth element (2026-10-23) is the fork's curved clipping surface of
# the dual advection (isoAdvectionCurved.C, keyword clippingSurface of the
# alpha solver dictionary): none (the plane sweep, today's loop),
# exactSphere (the exact sphere as the clipping surface for the volume flux
# AND the face update: the machinery gate, exact to the quadrature tolerance
# under uniform translation by construction) or paraboloid (the
# six-coefficient paraboloid of every interface cell fitted to its own face
# fractions with the ring tensor as prior, its own volume flux: the month-3
# milestone of the proposal). The exactSurface rows carry the normal model
# only for the recovered normal's error (the surface does not depend on it);
# the paraboloid rows' normal defines the frame of the fit (the slopes take
# the tilt): cellLocal (the own normal everywhere), cellSliver01 / cellSliver001
# (the own normal, the ring normal below |m_c| < 0.1 h^2 / 0.01 h^2), ringArea
# (the ring normal everywhere), exact (the control). Every existing row
# carries none, so its fvSolution gains the keyword block and nothing else
# changes.

# name -> (reconstructionScheme, normalModel, applyTo, sliverTol, maxAngle,
#          noNormalFallback, faceClipping, clippingSurface)
ROWS = OrderedDict([
    ("davof-cellLocal",          ("davof",    "cellLocal", "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("davof-ringArea",           ("davof",    "ringArea",  "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("davof-ringFit",            ("davof",    "ringFit",   "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("davof-exact",              ("davof",    "exact",     "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("davof-plicRDFnormal",      ("davof",    "plicRDF",   "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("davof-plicRDFnormal-ring", ("davof",    "plicRDF",   "all",   "0",   "15",  "ringArea", "plane",        "none")),
    ("davof-sagitta-exact",      ("davof",    "exact",     "all",   "0",   "15",  "keep",     "planeSagitta", "none")),
    ("davof-sagitta-ringArea",   ("davof",    "ringArea",  "all",   "0",   "15",  "keep",     "planeSagitta", "none")),
    ("davof-sagitta-cellSliver", ("davof",    "ringArea",  "check", "0.1", "180", "keep",     "planeSagitta", "none")),
    ("davof-sagitta-cellLocal",  ("davof",    "cellLocal", "all",   "0",   "15",  "keep",     "planeSagitta", "none")),
    ("davof-exactSurface-exact",     ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "exactSphere")),
    ("davof-exactSurface-cellLocal", ("davof", "cellLocal", "all",   "0",    "15",  "keep", "plane", "exactSphere")),
    ("davof-exactSurface-ringArea",  ("davof", "ringArea",  "all",   "0",    "15",  "keep", "plane", "exactSphere")),
    ("davof-paraboloid-cellLocal",   ("davof", "cellLocal", "all",   "0",    "15",  "keep", "plane", "paraboloid")),
    ("davof-paraboloid-cellSliver01",  ("davof", "ringArea", "check", "0.1",  "180", "keep", "plane", "paraboloid")),
    ("davof-paraboloid-cellSliver001", ("davof", "ringArea", "check", "0.01", "180", "keep", "plane", "paraboloid")),
    ("davof-paraboloid-ringArea",    ("davof", "ringArea",  "all",   "0",    "15",  "keep", "plane", "paraboloid")),
    ("davof-paraboloid-exact",       ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "paraboloid", "faces")),
    ("davof-paraboloidVol-exact",     ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume")),
    ("davof-paraboloidVol-cellLocal", ("davof", "cellLocal", "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume")),
    ("davof-paraboloidVol-cellSliver01",  ("davof", "ringArea", "check", "0.1",  "180", "keep", "plane", "paraboloid", "volume")),
    ("davof-paraboloidVol-cellSliver001", ("davof", "ringArea", "check", "0.01", "180", "keep", "plane", "paraboloid", "volume")),
    ("davof-paraboloidVol-ringArea",  ("davof", "ringArea",  "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume")),
    # the tenth element: the degrees of freedom the face fractions set in the
    # fit (six | one | none); the paraboloidRing rows take the frame of the
    # normal, the ring tensor and the offset from the volume alone (no face
    # information beyond the normal: the damping of classical VOF with a
    # second-order surface)
    ("davof-paraboloidRing-exact",     ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume", "none")),
    ("davof-paraboloidRing-cellLocal", ("davof", "cellLocal", "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume", "none")),
    ("davof-paraboloidRing-cellSliver01", ("davof", "ringArea", "check", "0.1", "180", "keep", "plane", "paraboloid", "volume", "none")),
    ("davof-paraboloidRing-ringArea",  ("davof", "ringArea",  "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume", "none")),
    # the frame from the ring quadric's normal at the cell (the Monge fit that
    # supplies the tensor), everything else as the ringArea row (2026-10-25)
    ("davof-paraboloidRing-ringFit",   ("davof", "ringFit",   "all",   "0",    "15",  "keep", "plane", "paraboloid", "volume", "none")),
    ("davof-paraboloidOne-exact",     ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "paraboloid", "faces",  "one")),
    ("davof-paraboloidOne-cellLocal",  ("davof", "cellLocal", "all",   "0",    "15",  "keep", "plane", "paraboloid", "faces",  "one")),
    # the plane clipping surface of the curved loop (fork branch davof-swept,
    # 2026-10-28, clippingSurface plane): the frame of the normal, a zero
    # tensor, the offset from the volume anchor; the surface the swept-face
    # flux (TPF_FLUX_INTEGRATION sweptFace) needs on a curved interface
    ("davof-planeSurf-exact",          ("davof", "exact",     "all",   "0",    "15",  "keep", "plane", "plane",      "volume", "none")),
    ("davof-planeSurf-ringArea",       ("davof", "ringArea",  "all",   "0",    "15",  "keep", "plane", "plane",      "volume", "none")),
    ("plicRDF",                  ("plicRDF",  "cellLocal", "all",   "0",   "15",  "keep",     "plane",        "none")),
    ("isoAlpha",                 ("isoAlpha", "cellLocal", "all",   "0",   "15",  "keep",     "plane",        "none")),
])
# the ninth element (2026-10-23, after the first smoke): the anchor of the
# paraboloid fit, faces (the six coefficients interpolate the cell's face
# fractions; the volume fraction never enters the surface: unstable under
# transport, growth 2.4 per step at N = 32, Co 0.2) or volume (the offset
# shifted after the fit so that the volume under the surface is alpha_k V,
# the feedback of the conserved volume fraction into the surface; the flux
# is the same surface's volume by the same quadrature). Rows without the
# element carry faces.
_FULL = {}
for _name, _vals in list(ROWS.items()):
    _v = _vals if len(_vals) >= 9 else (_vals + ("faces",))
    _v = _v if len(_v) >= 10 else (_v + ("six",))
    _FULL[_name] = _v
ROWS = OrderedDict(_FULL)

ROW_TOKENS = ("TPF_RECON_SCHEME", "TPF_NORMAL_MODEL", "TPF_APPLY_TO",
              "TPF_SLIVER_TOL", "TPF_MAX_ANGLE", "TPF_NO_NORMAL_FALLBACK",
              "TPF_FACE_CLIPPING", "TPF_CLIPPING_SURFACE", "TPF_FIT_ANCHOR",
              "TPF_FIT_DOF")

# arm -> (UX, UY, UZ); |U| = 1 for every arm (REFERENCE_SPEED 1 in the
# .parameter files, so dt = Co h). normal = the unit normal of the plane
# test (0.3 0.5 0.8124), |n| = 1 to 2e-5.
VELOCITY_ARMS = OrderedDict([
    ("normal",  ("0.3", "0.5", "0.8124")),
    ("oblique", ("1", "0", "0")),
    ("axial",   ("1", "0", "0")),
    # the solid-body rotation of cases/davof/sphereRotation3D (2026-10-28):
    # no uniform velocity; the field is omega x (x - axisPoint) from the
    # case's own tokens (TPF_OMEGA, TPF_AXIS_X/Y), the tokens here are the
    # U0 of the exact body, zero
    ("rotation", ("0", "0", "0")),
])

VELOCITY_TOKENS = ("TPF_UX", "TPF_UY", "TPF_UZ")

# Figure labels and ONE fixed hue per row, never cycled: the models that
# exist in make_davof_normal_table.COLORS keep their hue there (davof blue,
# plicRDF orange, isoAlpha yellow); the davof variants take distinct slots.
LABELS = OrderedDict([
    ("davof-cellLocal",     "DAVOF, own normal"),
    ("davof-ringArea",      "DAVOF, ring-area normal"),
    ("davof-ringFit",       "DAVOF, ring-fit normal"),
    ("davof-exact",         "DAVOF, exact normal (control)"),
    ("davof-plicRDFnormal", "DAVOF, plicRDF normal (control)"),
    ("davof-plicRDFnormal-ring", "DAVOF, plicRDF normal, ring fallback"),
    ("davof-sagitta-exact",      "DAVOF, sagitta clipping, exact normal"),
    ("davof-sagitta-ringArea",   "DAVOF, sagitta clipping, ring-area normal"),
    ("davof-sagitta-cellSliver", "DAVOF, sagitta clipping, own normal with sliver fallback"),
    ("davof-sagitta-cellLocal",  "DAVOF, sagitta clipping, own normal"),
    ("davof-exactSurface-exact",     "DAVOF, exact surface, exact normal (gate)"),
    ("davof-exactSurface-cellLocal", "DAVOF, exact surface, own normal"),
    ("davof-exactSurface-ringArea",  "DAVOF, exact surface, ring-area normal"),
    ("davof-paraboloid-cellLocal",   "DAVOF, paraboloid, own normal"),
    ("davof-paraboloid-cellSliver01",  "DAVOF, paraboloid, own normal, ring below 0.1 h^2"),
    ("davof-paraboloid-cellSliver001", "DAVOF, paraboloid, own normal, ring below 0.01 h^2"),
    ("davof-paraboloid-ringArea",    "DAVOF, paraboloid, ring-area normal"),
    ("davof-paraboloid-exact",       "DAVOF, paraboloid, exact normal (control)"),
    ("davof-paraboloidVol-exact",     "DAVOF, paraboloid, volume anchor, exact normal (control)"),
    ("davof-paraboloidVol-cellLocal", "DAVOF, paraboloid, volume anchor, own normal"),
    ("davof-paraboloidVol-cellSliver01",  "DAVOF, paraboloid, volume anchor, own normal, ring below 0.1 h^2"),
    ("davof-paraboloidVol-cellSliver001", "DAVOF, paraboloid, volume anchor, own normal, ring below 0.01 h^2"),
    ("davof-paraboloidVol-ringArea",  "DAVOF, paraboloid, volume anchor, ring-area normal"),
    ("davof-paraboloidRing-exact",     "DAVOF, ring-tensor paraboloid by volume, exact normal (control)"),
    ("davof-paraboloidRing-cellLocal", "DAVOF, ring-tensor paraboloid by volume, own normal"),
    ("davof-paraboloidRing-cellSliver01", "DAVOF, ring-tensor paraboloid by volume, own normal, ring below 0.1 h^2"),
    ("davof-paraboloidRing-ringArea",  "DAVOF, ring-tensor paraboloid by volume, ring-area normal"),
    ("davof-paraboloidRing-ringFit",   "DAVOF, ring-tensor paraboloid by volume, ring-quadric normal"),
    ("davof-paraboloidOne-exact",     "DAVOF, face-anchored paraboloid (one scalar), exact normal"),
    ("davof-paraboloidOne-cellLocal",  "DAVOF, face-anchored paraboloid (one scalar), own normal"),
    ("davof-planeSurf-exact",          "DAVOF, plane surface by volume, exact normal (control)"),
    ("davof-planeSurf-ringArea",       "DAVOF, plane surface by volume, ring-area normal"),
    ("plicRDF",             "plicRDF"),
    ("isoAlpha",            "isoAlpha"),
])

COLORS = {
    "davof-cellLocal":     "#2a78d6",   # DAVOF's blue (make_davof_normal_table.COLORS["davof"])
    "davof-ringArea":      "#1baf7a",   # aqua
    "davof-ringFit":       "#7a3bb5",   # purple
    "davof-exact":         "#444444",   # grey: the control
    "davof-plicRDFnormal": "#e87ba4",   # pink: the plicRDF-normal control
    "davof-plicRDFnormal-ring": "#9c2f5a",   # dark pink: the control with the ring fallback
    "davof-sagitta-exact":      "#8c8c8c",   # light grey: the sagitta control
    "davof-sagitta-ringArea":   "#0b6b48",   # dark aqua: ringArea with the sagitta clipping
    "davof-sagitta-cellSliver": "#1a4f8f",   # dark blue: the own normal with the sliver fallback
    "davof-sagitta-cellLocal":  "#7fb2ea",   # light blue: the own normal with the sagitta clipping
    "davof-exactSurface-exact":     "#222222",   # black: the gate
    "davof-exactSurface-cellLocal": "#5a5a5a",   # dark grey
    "davof-exactSurface-ringArea":  "#9a9a9a",   # grey
    "davof-paraboloid-cellLocal":   "#b5179e",   # magenta: the paraboloid with the own normal
    "davof-paraboloid-cellSliver01":  "#7209b7", # violet
    "davof-paraboloid-cellSliver001": "#3f37c9", # indigo
    "davof-paraboloid-ringArea":    "#4895ef",   # sky blue
    "davof-paraboloid-exact":       "#4cc9f0",   # cyan: the control
    "davof-paraboloidVol-exact":     "#8d99ae",   # slate: the anchored control
    "davof-paraboloidVol-cellLocal": "#d00000",   # red: the anchored own normal
    "davof-paraboloidVol-cellSliver01":  "#e85d04",  # orange-red
    "davof-paraboloidVol-cellSliver001": "#9d0208",  # dark red
    "davof-paraboloidVol-ringArea":  "#ffba08",   # amber
    "davof-paraboloidRing-exact":     "#006d77",  # teal: the ring-tensor paraboloid control
    "davof-paraboloidRing-cellLocal": "#83c5be",  # light teal
    "davof-paraboloidRing-cellSliver01": "#2a9d8f",  # sea green
    "davof-paraboloidRing-ringArea":  "#e29578",  # salmon
    "davof-paraboloidRing-ringFit":   "#7a3bb5",  # purple, as davof-ringFit
    "davof-paraboloidOne-exact":     "#6c757d",  # grey
    "davof-paraboloidOne-cellLocal":  "#adb5bd",  # light grey
    "davof-planeSurf-exact":          "#3a0ca3",  # dark violet: the plane-surface control
    "davof-planeSurf-ringArea":       "#f72585",  # pink: the plane surface with the ring normal
    "plicRDF":             "#eb6834",   # plicRDF's orange (make_davof_normal_table.COLORS)
    "isoAlpha":            "#eda100",   # isoAlpha's yellow (make_davof_normal_table.COLORS)
}

MARKERS = {
    "davof-cellLocal": "o", "davof-ringArea": "s", "davof-ringFit": "D",
    "davof-exact": "x", "davof-plicRDFnormal": "+", "davof-plicRDFnormal-ring": "P",
    "davof-sagitta-exact": "X", "davof-sagitta-ringArea": "p",
    "davof-sagitta-cellSliver": "h", "davof-sagitta-cellLocal": "8",
    "davof-exactSurface-exact": "*", "davof-exactSurface-cellLocal": "1", "davof-exactSurface-ringArea": "2",
    "davof-paraboloid-cellLocal": "<", "davof-paraboloid-cellSliver01": ">",
    "davof-paraboloid-cellSliver001": "d", "davof-paraboloid-ringArea": "H", "davof-paraboloid-exact": "|",
    "davof-paraboloidVol-exact": "_", "davof-paraboloidVol-cellLocal": "o",
    "davof-paraboloidVol-cellSliver01": "s", "davof-paraboloidVol-cellSliver001": "D",
    "davof-paraboloidVol-ringArea": "p",
    "davof-paraboloidRing-exact": "x", "davof-paraboloidRing-cellLocal": "o",
    "davof-paraboloidRing-cellSliver01": "s", "davof-paraboloidRing-ringArea": "D",
    "davof-paraboloidRing-ringFit": "p",
    "davof-paraboloidOne-exact": "+", "davof-paraboloidOne-cellLocal": "1",
    "davof-planeSurf-exact": "x", "davof-planeSurf-ringArea": "D",
    "plicRDF": "^", "isoAlpha": "v",
}


def is_davof(row):
    """True for the rows whose scheme stores the face fractions."""
    return ROWS[row][0] == "davof"


def row_tokens(row):
    """The seven template tokens of a row; RuntimeError on an unknown name
    (materialize runs in-process under snakemake: never SystemExit)."""
    if row not in ROWS:
        raise RuntimeError(f"[tpf_rows] unknown TPF_ROW {row!r}; known: {list(ROWS)}")
    return dict(zip(ROW_TOKENS, ROWS[row]))


def velocity_tokens(arm):
    """TPF_UX/UY/UZ of a velocity arm; RuntimeError on an unknown name."""
    if arm not in VELOCITY_ARMS:
        raise RuntimeError(f"[tpf_rows] unknown TPF_VELOCITY_ARM {arm!r}; "
                           f"known: {list(VELOCITY_ARMS)}")
    return dict(zip(VELOCITY_TOKENS, VELOCITY_ARMS[arm]))


if __name__ == "__main__":
    for name in ROWS:
        print(f"{name:22s} {row_tokens(name)}")
    for arm in VELOCITY_ARMS:
        print(f"{arm:22s} {velocity_tokens(arm)}")
