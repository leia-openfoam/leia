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

# name -> (reconstructionScheme, normalModel, applyTo, sliverTol, maxAngle,
#          noNormalFallback)
ROWS = OrderedDict([
    ("davof-cellLocal",          ("davof",    "cellLocal", "all", "0", "15", "keep")),
    ("davof-ringArea",           ("davof",    "ringArea",  "all", "0", "15", "keep")),
    ("davof-ringFit",            ("davof",    "ringFit",   "all", "0", "15", "keep")),
    ("davof-exact",              ("davof",    "exact",     "all", "0", "15", "keep")),
    ("davof-plicRDFnormal",      ("davof",    "plicRDF",   "all", "0", "15", "keep")),
    ("davof-plicRDFnormal-ring", ("davof",    "plicRDF",   "all", "0", "15", "ringArea")),
    ("plicRDF",                  ("plicRDF",  "cellLocal", "all", "0", "15", "keep")),
    ("isoAlpha",                 ("isoAlpha", "cellLocal", "all", "0", "15", "keep")),
])

ROW_TOKENS = ("TPF_RECON_SCHEME", "TPF_NORMAL_MODEL", "TPF_APPLY_TO",
              "TPF_SLIVER_TOL", "TPF_MAX_ANGLE", "TPF_NO_NORMAL_FALLBACK")

# arm -> (UX, UY, UZ); |U| = 1 for every arm (REFERENCE_SPEED 1 in the
# .parameter files, so dt = Co h). normal = the unit normal of the plane
# test (0.3 0.5 0.8124), |n| = 1 to 2e-5.
VELOCITY_ARMS = OrderedDict([
    ("normal",  ("0.3", "0.5", "0.8124")),
    ("oblique", ("1", "0", "0")),
    ("axial",   ("1", "0", "0")),
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
    "plicRDF":             "#eb6834",   # plicRDF's orange (make_davof_normal_table.COLORS)
    "isoAlpha":            "#eda100",   # isoAlpha's yellow (make_davof_normal_table.COLORS)
}

MARKERS = {
    "davof-cellLocal": "o", "davof-ringArea": "s", "davof-ringFit": "D",
    "davof-exact": "x", "davof-plicRDFnormal": "+", "davof-plicRDFnormal-ring": "P",
    "plicRDF": "^", "isoAlpha": "v",
}


def is_davof(row):
    """True for the rows whose scheme stores the face fractions."""
    return ROWS[row][0] == "davof"


def row_tokens(row):
    """The five template tokens of a row; RuntimeError on an unknown name
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
