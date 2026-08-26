#!/usr/bin/env python3
r"""Measure the 1D uniaxial-stretch case against its CLOSED-FORM solution.

WHY THIS CASE AND WHY THIS SCRIPT
---------------------------------
Every other verification case in this repository measures a SURROGATE: the
2D/3D benchmarks compare against a marker-traced reference, or against the
initial condition after a reversed flow, and both carry their own error. The 1D
uniaxial stretch has an exact solution in closed form for the field, its
gradient, the interface position, the phase measure AND the closest-point
velocity extension, so every number this script prints is a TRUE error.

It is also the cleanest possible gate for the SDPLS source, because the strain
the source must supply is an exact constant in space and time (below), so any
measured deviation is discretization and never a modelling error.

THE ANALYTIC SETUP
------------------
Velocity (OpenFOAM dimensionSet 0 0 -1 0 0 0 0, i.e. 1/s):

    v(x) = a (x - x0) e_x,        d(v_x)/dx = a  identically,   v_y = v_z = 0

with a the STRAIN RATE (default 1 1/s) and x0 the stagnation point. NOTE that
this flow is NOT divergence-free: div v = a. In one dimension it CANNOT be --
a 1D velocity with a non-zero strain has a non-zero divergence by definition.
That is intrinsic to the test, not a defect of the case, and it is restated
here because every reader of these numbers has to know it: the solver's psi
equation is assembled as

    ddt(psi) + div(phi,psi) - Sp(divPhi,psi) == S

and div(v psi) - psi div v = v . grad psi identically, so what is actually
solved is the advective form

    d(psi)/dt + v . grad(psi) = S

which is exact for a non-solenoidal v. No divergence-free assumption is used
anywhere in this script.

Level set: interface at x_i(0) = xi0, psi0(x) = x - xi0 (the exact 1D signed
distance). Interface motion, with and without the source alike, from the
characteristics dx/dt = a(x - x0):

    x_i(t) = x0 + (xi0 - x0) e^{a t}

The strain the source reads is, with nhat = grad(psi)/|grad psi| = e_x and
grad(v) = a e_x (x) e_x,

    a_strain = nhat . grad(v) . nhat = a         CONSTANT in space and time.

CLOSED FORMS (both derived along the characteristics x(t) = x0 + (x_s - x0)e^{at})

  (1) NO SOURCE, S = 0. psi is constant along characteristics:

        psi(x,t) = x0 + (x - x0) e^{-a t} - xi0,      d(psi)/dx = e^{-a t}

      The signed-distance property DECAYS EXPONENTIALLY. At a = 1 1/s, t = 1 s
      the gradient is e^{-1} = 0.3678794412 and |grad psi| - 1 = -0.6321205588.
      The interface is still at the exact x_i(t): losing the distance property
      does not move the zero contour, which is why the interface-position and
      phase-measure columns below are exact for BOTH arms.

  (2) SDPLS SOURCE R, S = a_strain psi = a psi. Then d(psi)/dt = a psi along
      the characteristics:

        psi(x,t) = x - x_i(t),                        d(psi)/dx = 1  EXACTLY

      The source holds the signed-distance property exactly, for all t. This is
      the verification target: the measured band |grad psi| must stay at 1 to
      discretization error while the sourceless arm decays as e^{-a t}, and
      BOTH arms must track x_i(t).

  (3) Rdiv / RdivStrictSp, the conservative form: A PREDICTION UNDER TEST, NOT
      an established result. Rdiv's derivation rests on a = -div_s(v), an
      identity its own header states holds for INCOMPRESSIBLE flow. In 1D the
      interface is a POINT: there are no tangential directions, div_s(v) = 0
      identically while the true strain is a = 1 1/s, so the identity fails.
      Concretely w = (nhat . v) nhat = v exactly here because v is entirely
      normal, hence phiW = phi, the two fvm::div terms cancel DISCRETELY as
      well as at the continuum, and div(w) - a = div(v) - a = 0: Rdiv assembles
      the ZERO MATRIX. PREDICTION: Rdiv behaves identically to noSource, with
      |grad psi| decaying as e^{-a t}. This script therefore scores an Rdiv arm
      against closed form (1) -- the reported errors ARE the test of that
      prediction, and are labelled as such in the `exactForm` column. Nothing
      here should be written up as established before the measurement exists.

SLAB VARIANT (two interfaces at x_L < x_R, phase inside). Both crossings ride
their own characteristic, so from x(t) = x0 + (x_s - x0)e^{at}

    x_L(t) = x0 + (x_L0 - x0) e^{a t},   x_R(t) = x0 + (x_R0 - x0) e^{a t}
    x_m(t) = (x_L + x_R)/2  (the medial axis),  W(t) = (x_R - x_L)/2 = W0 e^{a t}

and the exact phase measure is 2 W(t) = (x_R0 - x_L0) e^{a t} -- the slab is
advected AND stretched, so its measure is NOT constant; that expression is
derived here from the characteristics, it is not an assumption. The exact field
is psi = |x - x_m(t)| - W(t) for the source arm and e^{-a t} times the same for
the sourceless arm (substitute the foot point into psi0 and collect e^{-a t}),
so both closed forms collapse to

    psi_exact(x,t) = g(t) d_exact(x,t),      d(psi)/dx|_exact = g(t) s(x,t)

with g = 1 (source R) or e^{-a t} (sourceless), d_exact the exact signed
distance, and s = +-1 its sign. That single expression covers both geometries.

VELOCITY EXTENSION, the closestPoint diagnosis. Uext(x) = v(foot point of x):

  single interface: every band cell's foot point is x_i, so
        Uext(x) = a (x_i(t) - x0)      UNIFORM IN SPACE
      Then nhat . grad(Uext) = 0 exactly (not asymptotically), the advection is
      a rigid translation, |grad psi| = 1 is preserved with no source at all,
      and dx_i/dt = a(x_i - x0) is obeyed. In 1D with ONE interface,
      closestPoint is EXACT -- which is the regime the existing static gate
      measures, and why that gate reports O(h^1.98) while advection elsewhere
      is poor.

  slab: the closest-point map JUMPS at the medial axis,
        Uext(x) = a (x_L - x0) for x < x_m,   a (x_R - x0) for x > x_m
      a discontinuity of magnitude a (x_R - x_L). The controlling parameter is
      W / (nLayers h): the medial axis is inside the extension band exactly when
      nLayers h > W. A circle of radius 0.15 with a 3-cell band never has its
      medial axis (the centre) in the band, which is why the static gate never
      saw this; a filament thinning to ~1 cell does.

WHAT IS MEASURED, per written time
----------------------------------
 1. GRADIENT / signed-distance property. d(psi)/dx by second-order central
    differences on the cell centres (numpy.gradient's three-point non-uniform
    formula, which reduces to (psi_{i+1} - psi_{i-1})/(2h) on the uniform mesh
    the case runs; the two END cells get numpy's second-order one-sided
    formula and therefore inherit whatever the boundary condition did -- read
    the *Band columns, not the *All columns, when the ends matter). L1, L2 and
    Linf of (d(psi)/dx - exact) over the whole domain and over TWO bands:
      * ...Band     : |psi_measured| <= nLayers h -- the band the solver
                      actually has. On the sourceless arm psi shrinks like
                      e^{-a t}, so this band INFLATES in x; nBandCells makes
                      that visible and is part of the finding, not a nuisance.
      * ...GeomBand : |x - x_i(t)| <= nLayers h -- a mesh-fixed geometric band
                      that does not inflate, so the two arms stay comparable.
    The raw mean measured gradient is emitted too (signed and absolute), beside
    exp(-a t), so the exponential decay is readable directly off the CSV.
    Slab only: cells whose difference stencil STRADDLES the medial axis are
    excluded from the gradient norms (nKinkExcluded counts them). The exact
    d(psi)/dx jumps from -g to +g there, so a stencil across it carries an O(1)
    error that is a property of the exact solution's kink, not of the scheme.
 2. FIELD error. L1/L2/Linf of psi(x,t) minus the exact psi of whichever arm
    ran -- closed form (1) for noSource/Rdiv, closed form (2) for R (see
    --source). No kink exclusion: psi_exact is continuous there.
 3. INTERFACE POSITION -- the shape error in 1D. The zero crossing of psi by
    linear interpolation between the two cells that bracket it, signed against
    x_i(t), reported in metres AND in cell widths. Both crossings for a slab.
 4. PHASE VOLUME: the measure of {psi < 0}, from the written alpha field if
    present (alpha = 1 where psi < 0, detrixheAslamPhaseIndicator.C), else by
    sub-cell linear reconstruction of the zero crossing. Exact value: the
    distance from the domain end to x_i(t) for a single interface,
    (x_R0 - x_L0) e^{a t} for the slab. Absolute and relative error. Both
    estimates are emitted so alpha and the reconstruction can be compared.
 5. VELOCITY EXTENSION, the closestPoint diagnosis, when a Uext field is
    written. |Uext - Uext_exact| (the FULL vector norm, so a spurious y/z
    component is caught) with L2/Linf over the band, and -- the point of the
    whole diagnosis -- SPLIT for the slab into cells whose foot point the
    solver picked on their own side of the medial axis and cells where it did
    not. The side the solver used is inferred from Uext itself: the exact
    extension takes only the two values a(x_L - x0) and a(x_R - x0), so the
    nearer of the two names the interface the extension actually reached for.
    A jump at the medial axis shows up as an Linf that does not converge (and
    a non-zero nUextCrossSide) while the L2 away from it does. The location of
    the largest pointwise error is reported relative to the medial axis, in
    cell widths. uextSplitReliable is 0 when Linf exceeds half the jump
    a(x_R - x_L), because then the nearest-value classification can itself be
    flipped by the error and the split must not be read.

REPORTED TOGETHER, ALWAYS: gradient, shape (interface position) and volume
(phase measure) appear in every row and in the stdout summary. A single-metric
view of this case would misread it -- a source can hold the gradient while
displacing the interface, and a volume correction can restore the measure while
moving the zero set.

WHERE THE PARAMETERS COME FROM
------------------------------
Read from the case itself (in this order; the first that yields a value wins,
and the provenance of every one is printed):

  strain rate a   system/fvSolution `velocityModel` -- keys strainRate | a |
                  rate | strain; else a least-squares fit of the WRITTEN U
                  field, U_x = a(x - x0), which is also always computed and
                  cross-checked against the dictionary value (a disagreement is
                  reported: the fit is what the solver actually advected with).
  x0              `velocityModel` stagnationPoint | point | origin | x0 (a
                  vector's x-component, or a scalar); else the U fit.
  xi0, geometry   levelSet/implicitSurface (implicitPlane `position` with
                  `normal` along x, or hesseNormalPlane `n`/`d`); else, and
                  this is the robust route, the zero crossings of 0/psi -- one
                  crossing means a single interface, two mean a slab, and the
                  slab half-width W0 = (x_R0 - x_L0)/2 follows from them.
  nLayers         levelSet/velocityExtension/nLayers, else the sdplsSource
                  mollifier's nLayers, else case_params.json, else 3.
  source arm      levelSet/sdplsSource/type -> which closed form is exact.
  cell centres    a written C field (postProcess -func writeCellCentres) if
                  present -- no ordering assumption at all -- else the unique x
                  planes of constant/polyMesh/points, whose midpoints are the
                  cell centres; the assumed blockMesh ordering (x fastest, one
                  cell in y and z, so cells ascend in x) is then VERIFIED
                  against 0/psi, which must be monotone (single interface) or
                  V-shaped (slab). cellOrderVerified records the outcome.

Explicit arguments exist only as overrides and as the fallback for what cannot
be read: --a, --x0, --xi0, --half-width, --nlayers, --geometry, --source,
--domain.

Nothing here is allowed to crash on a missing optional field: a time directory
without alpha, without Uext, or with an unparsable field yields EMPTY columns
for exactly those quantities and a note on stdout. A reader that dies halfway
through has previously made a complete study look incomplete in this project.

OUTPUT
------
  <case>/analyse1Dstretch.csv    one row per written time, stable header
  stdout                         a short human-readable summary, with the
                                 measured mean |d(psi)/dx| beside e^{-a t} and
                                 beside 1, so the verification reads at a glance

Usage (from the repo root):
    python3 workflow/scripts/analyse_1Dstretch.py <case-dir>
    python3 workflow/scripts/analyse_1Dstretch.py <case-dir> --source R
    python3 workflow/scripts/analyse_1Dstretch.py <case-dir> --geometry slab \
        --a 1 --x0 0.5 --nlayers 3 --out /tmp/row.csv
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

# Reuse the repository's existing readers rather than adding a third one.
#   * fvschemes carries the tolerant regex dictionary parser that
#     check_discretization.py already trusts to read what a case ACTUALLY ran.
#     foamlib's dictionary parser is NOT usable here: it rejects the rendered
#     fvSolution outright ("Uext.*" { $noneUext; } -- regex keys plus macro
#     substitution), which is exactly the file every parameter below lives in.
#   * sdpls_source_error.read_internal is the field reader the SDPLS source
#     gate uses, kept as the fallback for environments without foamlib.
import fvschemes                      # noqa: E402
import sdpls_source_error             # noqa: E402

try:                                  # declared in requirements.txt, used by plots.py
    from foamlib import FoamFieldFile
except Exception:                     # noqa: BLE001 -- absence is not an error here
    FoamFieldFile = None


# --------------------------------------------------------------------------- #
# CSV schema. STABLE: every column is written on every row, empty when the
# quantity does not exist for that row (no alpha field, single interface so no
# medial axis, ...). Downstream readers may index by name and must never have
# to cope with a header that changes between arms.
# --------------------------------------------------------------------------- #
COLUMNS = [
    # -- what this row is, and the parameters it was scored with ------------ #
    "time", "nCells", "h", "strainRate", "x0", "geometry", "orientation",
    "sourceType", "exactForm", "velocityExtension", "nLayers", "bandWidth",
    "cellOrderVerified", "interfaceInDomain",
    # -- 1. gradient / signed-distance property ----------------------------- #
    "exactDpsiDx", "decayExpMinusAT",
    "meanDpsiDxAll", "meanAbsDpsiDxAll", "meanDpsiDxBand", "meanAbsDpsiDxBand",
    "gradErrL1All", "gradErrL2All", "gradErrLinfAll",
    "gradErrL1Band", "gradErrL2Band", "gradErrLinfBand", "nBandCells",
    "gradErrL1GeomBand", "gradErrL2GeomBand", "gradErrLinfGeomBand",
    "nGeomBandCells", "nKinkExcluded",
    # -- 2. field error ----------------------------------------------------- #
    "psiErrL1All", "psiErrL2All", "psiErrLinfAll",
    "psiErrL1Band", "psiErrL2Band", "psiErrLinfBand",
    # -- 3. interface position = the 1D shape error ------------------------- #
    "nZeroCrossings",
    "xiExact", "xiMeasured", "xiError", "xiErrorCells",
    "xLExact", "xLMeasured", "xLError", "xLErrorCells",
    "xRExact", "xRMeasured", "xRError", "xRErrorCells",
    # -- 4. phase volume: the measure of {psi < 0} -------------------------- #
    "phaseMeasure", "phaseMeasureSource", "phaseMeasureAlpha",
    "phaseMeasureRecon", "phaseMeasureExact", "phaseMeasureAbsError",
    "phaseMeasureRelError",
    # -- 5. velocity extension: the closestPoint diagnosis ------------------ #
    "uextPresent", "uextExactKind", "xMedial",
    "uextErrL2Band", "uextErrLinfBand", "nUextBandCells",
    "uextErrL2SameSide", "uextErrLinfSameSide", "nUextSameSide",
    "uextErrL2CrossSide", "uextErrLinfCrossSide", "nUextCrossSide",
    "uextSplitReliable",
    "uextArgmaxX", "uextArgmaxCellsFromMedial", "uextArgmaxCellsFromInterface",
]

# sdplsSource type -> which closed form is exact for that arm.
#   "none" = closed form (1), psi constant along characteristics
#   "R"    = closed form (2), psi = x - x_i(t)
# Rdiv/RdivStrictSp are scored against (1) because that is the PREDICTION under
# test (see the module docstring, item 3), never because it is established.
_SOURCE_EXACT_FORM = {
    "noSource": "none",
    "none": "none",
    "R": "R",
    "Rdiv": "none",
    "RdivStrictSp": "none",
}


# --------------------------------------------------------------------------- #
# Field reading
# --------------------------------------------------------------------------- #
def _is_uniform(path):
    """True/False if the file states a uniform/nonuniform internalField, else None.

    Needed because an OpenFOAM `internalField uniform (0 0 0)` parses to a bare
    3-vector, which would silently be read as a three-cell mesh. Uext is written
    for every model (velocityExtension.C constructs it AUTO_WRITE), and with
    `type none` before the first correct() it is exactly that uniform zero.
    """
    try:
        with open(path, "r", errors="replace") as fh:
            head = fh.read(8192)
    except OSError:
        return None
    m = re.search(r"internalField\s+(uniform|nonuniform)", head)
    return None if not m else (m.group(1) == "uniform")


def read_field(path, ncells=None, what=""):
    """internalField of an OpenFOAM volField as a numpy array, or None.

    Never raises: an absent, truncated or binary-unreadable field must leave
    empty columns, not kill the run.
    """
    if not os.path.isfile(path):
        return None
    arr = None
    if FoamFieldFile is not None:
        try:
            arr = np.asarray(FoamFieldFile(path).internal_field, dtype=float)
        except Exception as exc:                       # noqa: BLE001
            print(f"  [read] foamlib could not read {path}: "
                  f"{type(exc).__name__}: {exc} -- falling back to the regex reader")
            arr = None
    if arr is None:
        try:
            raw = sdpls_source_error.read_internal(path)
        except Exception as exc:                       # noqa: BLE001
            print(f"  [read] {path}: {type(exc).__name__}: {exc} -- column left empty")
            return None
        if raw is None:
            return None
        arr = np.asarray(raw, dtype=float)

    if _is_uniform(path) and ncells:
        flat = arr.reshape(-1)
        if flat.size == 1:
            arr = np.full(ncells, float(flat[0]))
        elif flat.size == 3:
            arr = np.tile(flat, (ncells, 1))
    if ncells and arr.shape[0] != ncells:
        print(f"  [read] {what or os.path.basename(path)}: {arr.shape[0]} entries "
              f"but {ncells} cells -- column left empty")
        return None
    return arr


def time_dirs(case):
    """[(t, path)] for every written time directory that HAS a psi field."""
    out, skipped = [], []
    for name in os.listdir(case):
        d = os.path.join(case, name)
        if not os.path.isdir(d):
            continue
        try:
            t = float(name)
        except ValueError:
            continue                              # 0.org, constant, system, ...
        if os.path.isfile(os.path.join(d, "psi")):
            out.append((t, d))
        else:
            skipped.append(name)
    if skipped:
        print(f"  [times] {len(skipped)} time director(ies) without a psi field, "
              f"skipped: {', '.join(sorted(skipped)[:8])}")
    return sorted(out)


# --------------------------------------------------------------------------- #
# Dictionaries
# --------------------------------------------------------------------------- #
def _shallow_entries(body):
    """Entries of THIS dictionary level only, sub-dictionaries removed.

    fvschemes._entries scans everything it is handed, so calling it on a block
    that CONTAINS sub-dictionaries merges their keys in and the LAST one wins:
    `sdplsSource { type R; mollifier { type band; } }` would be reported as
    type = band, i.e. the arm would be misidentified and the wrong closed form
    scored against. Strip nested blocks (and the sub-dictionary name in front
    of them) before parsing.
    """
    keep, depth = [], 0
    for ch in body:
        if ch == "{":
            if depth == 0:
                while keep and keep[-1] not in ";\n":
                    keep.pop()                    # drop the sub-dictionary's name
            depth += 1
        elif ch == "}":
            depth = max(0, depth - 1)
        elif depth == 0:
            keep.append(ch)
    return fvschemes._entries("".join(keep))


def read_dicts(case):
    """velocityModel / implicitSurface / velocityExtension / sdplsSource entries.

    Reads the RENDERED system/fvSolution, which is what the solver parsed --
    never the study tokens, for the reasons fvschemes.py documents at length
    (a token may be absent, hard-coded past, or an unexpanded alias).
    """
    out = {"velocityModel": {}, "implicitSurface": {},
           "velocityExtension": {}, "sdplsSource": {}, "levelSet": ""}
    path = os.path.join(case, "system", "fvSolution")
    if not os.path.isfile(path):
        print(f"  [dict] no {path} -- every parameter falls back to arguments")
        return out
    try:
        text = fvschemes._strip_comments(open(path).read())
    except OSError as exc:
        print(f"  [dict] cannot read {path}: {exc}")
        return out
    level_set = fvschemes._block(text, "levelSet")
    out["levelSet"] = level_set
    out["velocityModel"] = _shallow_entries(fvschemes._block(text, "velocityModel"))
    for key in ("implicitSurface", "velocityExtension", "sdplsSource"):
        out[key] = _shallow_entries(fvschemes._block(level_set, key))
    return out


def _tokens(case):
    path = os.path.join(case, "case_params.json")
    if not os.path.isfile(path):
        return {}
    try:
        return json.load(open(path)).get("tokens", {})
    except (OSError, ValueError) as exc:
        print(f"  [dict] cannot read {path}: {exc}")
        return {}


def _as_float(value):
    """First float in an OpenFOAM value string, or None. `(0.5 0 0)` -> 0.5."""
    if value is None:
        return None
    m = re.search(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", str(value))
    if not m:
        return None
    try:
        return float(m.group(0).replace("d", "e").replace("D", "e"))
    except ValueError:
        return None


def _as_vector(value):
    """`(a b c)` -> [a,b,c], else None."""
    if value is None:
        return None
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    if len(nums) < 3:
        return None
    return [float(v) for v in nums[:3]]


def _first(dct, keys):
    """(value, key) for the first present key, else (None, None)."""
    for k in keys:
        if k in dct and str(dct[k]).strip():
            return dct[k], k
    return None, None


# --------------------------------------------------------------------------- #
# Mesh
# --------------------------------------------------------------------------- #
def _points_x_planes(case):
    """Sorted unique x of constant/polyMesh/points, or None.

    For a 1D column mesh (one cell in y and z) these are the cell FACES normal
    to x, so consecutive midpoints are the cell centres -- exact for a graded
    mesh as well as a uniform one, and needing nothing from blockMeshDict
    (which may carry macros, `scale`, or `convertToMeters`).
    """
    path = os.path.join(case, "constant", "polyMesh", "points")
    if not os.path.isfile(path):
        return None
    try:
        text = open(path, "r", errors="replace").read()
    except OSError:
        return None
    xs = [float(m.group(1)) for m in
          re.finditer(r"\(\s*([-+0-9.eE]+)\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\)", text)]
    if not xs:
        return None
    xs = np.asarray(xs)
    # Collapse to distinct planes with a tolerance relative to the extent: the
    # ascii points are written at finite precision, so exact uniqueness would
    # split one plane into several.
    span = xs.max() - xs.min()
    if span <= 0:
        return None
    xs = np.sort(xs)
    keep = [xs[0]]
    for v in xs[1:]:
        if v - keep[-1] > 1e-9 * span:
            keep.append(v)
    return np.asarray(keep)


def cell_centres_x(case, tdirs, ncells, domain_arg):
    """(x centres ascending, cell widths, permutation, provenance string).

    `permutation` maps a field array in OpenFOAM cell order onto the ascending-x
    order used by everything downstream; it is the identity for the blockMesh
    ordering this case uses, and is only ever non-trivial when a written C field
    says so.
    """
    # 1. A written cell-centre field: no ordering assumption at all.
    for _t, d in [(None, os.path.join(case, "0"))] + [(t, d) for t, d in tdirs]:
        cfile = os.path.join(d, "C")
        if os.path.isfile(cfile):
            C = read_field(cfile, ncells, "C")
            if C is not None and C.ndim == 2 and C.shape[1] == 3:
                perm = np.argsort(C[:, 0], kind="stable")
                xs = C[perm, 0]
                planes = _points_x_planes(case)
                if planes is not None and planes.size == xs.size + 1:
                    dx = np.diff(planes)
                else:
                    dx = _widths_from_centres(xs)
                return xs, dx, perm, f"written C field ({os.path.relpath(cfile, case)})"

    # 2. The mesh points: cell centres are the midpoints of the x planes.
    planes = _points_x_planes(case)
    if planes is not None and planes.size == ncells + 1:
        xs = 0.5 * (planes[:-1] + planes[1:])
        return (xs, np.diff(planes), np.arange(ncells),
                "constant/polyMesh/points (x planes; cell order assumed "
                "blockMesh x-fastest, VERIFIED against 0/psi)")
    if planes is not None:
        print(f"  [mesh] constant/polyMesh/points has {planes.size} distinct x "
              f"planes but the field has {ncells} cells: this is not a 1D column "
              f"mesh (or it is decomposed). Falling back.")

    # 3. Explicit domain argument.
    if domain_arg:
        lo, hi = [float(v) for v in domain_arg.split(",")]
        planes = np.linspace(lo, hi, ncells + 1)
        xs = 0.5 * (planes[:-1] + planes[1:])
        return (xs, np.diff(planes), np.arange(ncells),
                f"--domain {lo},{hi} (uniform, {ncells} cells)")
    return None, None, None, None


def _widths_from_centres(xs):
    """Cell widths implied by centres alone (used only when the faces are unknown)."""
    dx = np.empty_like(xs)
    if xs.size == 1:
        return np.ones_like(xs)
    dx[1:-1] = 0.5 * (xs[2:] - xs[:-2])
    dx[0] = xs[1] - xs[0]
    dx[-1] = xs[-1] - xs[-2]
    return dx


def verify_cell_order(psi0):
    """Is psi at t = 0 monotone (single interface) or V-shaped (slab)?

    psi0 is the exact 1D signed distance, so in ascending-x order it is either
    strictly monotone or decreasing-then-increasing. Any other pattern means the
    assumed blockMesh cell ordering is wrong and every number below would be
    scrambled -- so it is checked rather than assumed.
    """
    d = np.diff(psi0)
    if d.size == 0:
        return True, "single cell"
    s = np.sign(d)
    s = s[s != 0]
    if s.size == 0:
        return False, "psi is constant at t = 0"
    flips = int(np.count_nonzero(np.diff(s) != 0))
    if flips == 0:
        return True, "monotone (single interface)"
    if flips == 1 and s[0] < 0:
        return True, "V-shaped (slab)"
    return False, f"{flips} sign changes in d(psi)/dx at t = 0"


# --------------------------------------------------------------------------- #
# Small numerics
# --------------------------------------------------------------------------- #
def norms(err, w):
    """Cell-width-weighted L1, L2 and Linf of an error array.

    The weighting is identical to plain averaging on the uniform mesh this case
    runs; it is written this way so a graded 1D mesh does not silently weight
    small cells as heavily as large ones.
    """
    err = np.asarray(err, dtype=float)
    if err.size == 0:
        return None, None, None
    w = np.asarray(w, dtype=float)
    tot = float(w.sum())
    if tot <= 0:
        return None, None, None
    return (float(np.abs(err) @ w / tot),
            float(math.sqrt(float((err ** 2) @ w / tot))),
            float(np.abs(err).max()))


def zero_crossings(x, f):
    """[(i, x_cross)] of f, by linear interpolation between bracketing cells."""
    out = []
    n = f.size
    for i in range(n - 1):
        a, b = f[i], f[i + 1]
        if a == 0.0:
            out.append((i, float(x[i])))
            continue
        if a * b < 0.0:
            out.append((i, float(x[i] + (x[i + 1] - x[i]) * (-a) / (b - a))))
    if n and f[-1] == 0.0:
        out.append((n - 1, float(x[-1])))
    return out


def negative_measure(x, psi, xlo, xhi, crossings):
    """Measure of {psi < 0} by sub-cell linear reconstruction of the crossings.

    The piecewise-linear interpolant of psi through the cell centres changes
    sign only at the interpolated crossings, so the breakpoints
    [xlo, crossings..., xhi] cut the domain into intervals of constant sign and
    the measure is the total length of the negative ones. Exact for the 1D
    configurations here and general in the number of crossings.
    """
    bps = [xlo] + [xc for _i, xc in crossings] + [xhi]
    bps = sorted(set(bps))
    total = 0.0
    for k in range(len(bps) - 1):
        mid = 0.5 * (bps[k] + bps[k + 1])
        if float(np.interp(mid, x, psi)) < 0.0:
            total += bps[k + 1] - bps[k]
    return total


def _fit_strain(x, ux):
    """(a, x0) from a least-squares fit of U_x = a (x - x0), or (None, None).

    This is what the solver ACTUALLY advected with, so it both supplies the
    parameters when the dictionary keys cannot be found and cross-checks them
    when they can.
    """
    if x.size < 2 or not np.all(np.isfinite(ux)):
        return None, None
    a, b = np.polyfit(x, ux, 1)                   # ux = a x + b
    if abs(a) < 1e-300:
        return float(a), None
    return float(a), float(-b / a)


# --------------------------------------------------------------------------- #
# Parameter resolution
# --------------------------------------------------------------------------- #
def resolve_parameters(case, args, dicts, tokens, x, psi0, u0):
    """Every geometry/flow parameter, with the provenance of each.

    Priority: an explicit command-line override, then the case's own rendered
    dictionaries, then a fit of the written U field, then case_params.json, then
    the documented default. Each is reported so a number can never be traced
    back to the wrong source.
    """
    p, why = {}, {}
    vm, isf, ve, sd = (dicts["velocityModel"], dicts["implicitSurface"],
                       dicts["velocityExtension"], dicts["sdplsSource"])

    # ---- strain rate a and stagnation point x0 ---------------------------- #
    a_fit = x0_fit = None
    if u0 is not None and u0.ndim == 2 and u0.shape[1] == 3:
        a_fit, x0_fit = _fit_strain(x, u0[:, 0])

    val, key = _first(vm, ("strainRate", "a", "rate", "strain"))
    if args.a is not None:
        p["a"], why["a"] = args.a, "--a (command line override)"
    elif _as_float(val) is not None:
        p["a"], why["a"] = _as_float(val), f"fvSolution velocityModel/{key}"
    elif a_fit is not None:
        p["a"], why["a"] = a_fit, "least-squares fit of the written U field"
    elif _as_float(tokens.get("STRAIN_RATE")) is not None:
        p["a"], why["a"] = _as_float(tokens["STRAIN_RATE"]), "case_params.json STRAIN_RATE"
    else:
        p["a"], why["a"] = 1.0, "DEFAULT 1 1/s (nothing in the case said otherwise)"

    val, key = _first(vm, ("stagnationPoint", "point", "origin", "x0", "centre", "center"))
    vec = _as_vector(val)
    if args.x0 is not None:
        p["x0"], why["x0"] = args.x0, "--x0 (command line override)"
    elif vec is not None:
        p["x0"], why["x0"] = vec[0], f"fvSolution velocityModel/{key} (x component)"
    elif _as_float(val) is not None:
        p["x0"], why["x0"] = _as_float(val), f"fvSolution velocityModel/{key}"
    elif x0_fit is not None:
        p["x0"], why["x0"] = x0_fit, "least-squares fit of the written U field"
    elif _as_float(tokens.get("STAGNATION_POINT") or tokens.get("X0")) is not None:
        p["x0"], why["x0"] = _as_float(tokens.get("STAGNATION_POINT") or tokens["X0"]), \
            "case_params.json STAGNATION_POINT/X0"
    else:
        p["x0"], why["x0"] = 0.0, "DEFAULT 0 m (nothing in the case said otherwise)"

    p["aFit"], p["x0Fit"] = a_fit, x0_fit

    # ---- geometry: single interface or slab, and the initial positions ---- #
    cr0 = zero_crossings(x, psi0)
    p["orientation"] = 1.0 if (psi0[-1] - psi0[0]) >= 0 else -1.0

    geom = args.geometry
    if geom == "auto":
        if len(cr0) >= 2:
            geom, why["geometry"] = "slab", f"{len(cr0)} zero crossings in 0/psi"
        elif len(cr0) == 1:
            geom, why["geometry"] = "single", "one zero crossing in 0/psi"
        else:
            geom, why["geometry"] = "single", \
                "NO zero crossing in 0/psi -- assumed single interface; " \
                "the interface columns will be empty"
    else:
        why["geometry"] = f"--geometry {geom} (command line override)"
    p["geometry"] = geom

    if geom == "slab":
        if args.xi0 is not None and args.half_width is not None:
            p["xL0"], p["xR0"] = args.xi0 - args.half_width, args.xi0 + args.half_width
            why["slab"] = "--xi0 (centre) and --half-width (command line override)"
        elif len(cr0) >= 2:
            p["xL0"], p["xR0"] = cr0[0][1], cr0[-1][1]
            why["slab"] = "the two zero crossings of 0/psi"
        else:
            hw = args.half_width if args.half_width is not None \
                else _as_float(tokens.get("SLAB_HALF_WIDTH"))
            c = args.xi0 if args.xi0 is not None else _as_float(tokens.get("XI0"))
            if hw is None or c is None:
                raise SystemExit(
                    "[analyse_1Dstretch] --geometry slab but 0/psi has fewer than two "
                    "zero crossings and no --xi0/--half-width were given.")
            p["xL0"], p["xR0"] = c - hw, c + hw
            why["slab"] = "--xi0/--half-width or case_params.json XI0/SLAB_HALF_WIDTH"
        p["W0"] = 0.5 * (p["xR0"] - p["xL0"])
        p["xi0"] = None
    else:
        val, key = _first(isf, ("position",))
        vec = _as_vector(val)
        nrm = _as_vector(isf.get("normal"))
        if args.xi0 is not None:
            p["xi0"], why["xi0"] = args.xi0, "--xi0 (command line override)"
        elif cr0:
            p["xi0"], why["xi0"] = cr0[0][1], "the zero crossing of 0/psi"
        elif vec is not None and (nrm is None or abs(nrm[0]) > 0.9):
            p["xi0"], why["xi0"] = vec[0], \
                f"fvSolution levelSet/implicitSurface/{key} (x component)"
        elif isf.get("type", "").strip() == "hesseNormalPlane":
            n = _as_vector(isf.get("n"))
            d = _as_float(isf.get("d"))
            if n and d is not None and abs(n[0]) > 1e-12:
                p["xi0"], why["xi0"] = d / n[0], "levelSet/implicitSurface d/n_x"
            else:
                p["xi0"], why["xi0"] = None, "hesseNormalPlane not resolvable"
        elif _as_float(tokens.get("XI0")) is not None:
            p["xi0"], why["xi0"] = _as_float(tokens["XI0"]), "case_params.json XI0"
        else:
            p["xi0"], why["xi0"] = None, \
                "NOT RESOLVED -- interface/volume columns will be empty"
        p["xL0"] = p["xR0"] = p["W0"] = None

    # ---- narrow band width ------------------------------------------------ #
    if args.nlayers is not None:
        p["nLayers"], why["nLayers"] = args.nlayers, "--nlayers (command line override)"
    elif _as_float(ve.get("nLayers")) is not None:
        p["nLayers"], why["nLayers"] = _as_float(ve["nLayers"]), \
            "fvSolution levelSet/velocityExtension/nLayers"
    else:
        moll = _shallow_entries(fvschemes._block(dicts["levelSet"], "mollifier"))
        if _as_float(moll.get("nLayers")) is not None:
            p["nLayers"], why["nLayers"] = _as_float(moll["nLayers"]), \
                "fvSolution levelSet/sdplsSource/mollifier/nLayers"
        elif _as_float(tokens.get("SDPLS_NLAYERS") or tokens.get("N_LAYERS")) is not None:
            p["nLayers"], why["nLayers"] = _as_float(
                tokens.get("SDPLS_NLAYERS") or tokens["N_LAYERS"]), "case_params.json"
        else:
            p["nLayers"], why["nLayers"] = 3.0, "DEFAULT 3 layers"

    # ---- which closed form is exact for this arm -------------------------- #
    p["sourceType"] = (sd.get("type") or tokens.get("SDPLS_SOURCE") or "").strip()
    if args.source != "auto":
        p["exactForm"], why["exactForm"] = args.source, \
            f"--source {args.source} (command line override)"
    elif p["sourceType"] in _SOURCE_EXACT_FORM:
        p["exactForm"] = _SOURCE_EXACT_FORM[p["sourceType"]]
        why["exactForm"] = f"levelSet/sdplsSource/type = {p['sourceType']}"
        if p["sourceType"] in ("Rdiv", "RdivStrictSp"):
            why["exactForm"] += (" -- scored against closed form (1), which is the "
                                 "PREDICTION under test (Rdiv assembles the zero "
                                 "matrix in 1D), NOT an established result")
    else:
        # `beta` (f_nl = beta - |grad psi|, so S = psi(beta - |grad psi|)) is
        # nonlinear and has no closed form derived in this script; the same goes
        # for any source type added later. Score against the sourceless form and
        # SAY SO in the column, so nobody reads the number as a true error.
        p["exactForm"] = "none"
        why["exactForm"] = (
            f"levelSet/sdplsSource/type = {p['sourceType']!r} has no closed form "
            "derived here; scored against the SOURCELESS form (1), so the field and "
            "gradient errors are DEVIATIONS FROM (1) and must be read as such")

    p["velocityExtension"] = (ve.get("type") or tokens.get("VELOCITY_EXTENSION") or "").strip()
    return p, why


# --------------------------------------------------------------------------- #
# Per-time evaluation
# --------------------------------------------------------------------------- #
def exact_state(t, p):
    """Everything the closed forms need at time t.

    x_i(t) = x0 + (xi0 - x0) e^{a t} for the single interface; the same map
    applied to BOTH crossings for the slab, whence x_m(t) and W(t) = W0 e^{a t}.
    g(t) = 1 for the R arm (closed form 2), e^{-a t} for the sourceless arm
    (closed form 1), and psi_exact = g d_exact, d(psi)/dx|_exact = g s in both
    geometries (see the module docstring).
    """
    a, x0 = p["a"], p["x0"]
    E = math.exp(a * t)
    st = {"E": E, "decay": math.exp(-a * t)}
    st["g"] = 1.0 if p["exactForm"] == "R" else st["decay"]
    if p["geometry"] == "slab":
        st["xL"] = x0 + (p["xL0"] - x0) * E
        st["xR"] = x0 + (p["xR0"] - x0) * E
        st["xm"] = 0.5 * (st["xL"] + st["xR"])
        st["W"] = 0.5 * (st["xR"] - st["xL"])
        st["xi"] = None
    else:
        st["xi"] = (x0 + (p["xi0"] - x0) * E) if p["xi0"] is not None else None
        st["xL"] = st["xR"] = st["xm"] = st["W"] = None
    return st


def evaluate(t, tdir, ctx):
    """One CSV row: gradient, field, interface, phase measure and Uext at time t."""
    x, dx, perm, p = ctx["x"], ctx["dx"], ctx["perm"], ctx["p"]
    n, h = x.size, ctx["h"]
    xlo, xhi = ctx["xlo"], ctx["xhi"]
    st = exact_state(t, p)

    row = {c: "" for c in COLUMNS}
    row.update(time=t, nCells=n, h=h, strainRate=p["a"], x0=p["x0"],
               geometry=p["geometry"], orientation=p["orientation"],
               sourceType=p["sourceType"], exactForm=p["exactForm"],
               velocityExtension=p["velocityExtension"], nLayers=p["nLayers"],
               bandWidth=p["nLayers"] * h, cellOrderVerified=ctx["orderOK"],
               decayExpMinusAT=st["decay"])

    psi = read_field(os.path.join(tdir, "psi"), n, f"psi at t = {t:g}")
    if psi is None:
        return row
    psi = psi[perm]

    # ---- exact signed distance, its sign, and the exact field ------------- #
    if p["geometry"] == "slab":
        s_exact = np.where(x >= st["xm"], 1.0, -1.0)
        d_exact = np.abs(x - st["xm"]) - st["W"]
        nearest = np.minimum(np.abs(x - st["xL"]), np.abs(x - st["xR"]))
        in_domain = (xlo < st["xL"] < xhi) and (xlo < st["xR"] < xhi)
    elif st["xi"] is not None:
        s_exact = np.full(n, p["orientation"])
        d_exact = p["orientation"] * (x - st["xi"])
        nearest = np.abs(x - st["xi"])
        in_domain = xlo < st["xi"] < xhi
    else:
        s_exact = d_exact = nearest = None
        in_domain = None
    row["interfaceInDomain"] = "" if in_domain is None else int(in_domain)

    # ---- 1. gradient / signed-distance property --------------------------- #
    # numpy.gradient: the three-point non-uniform central formula in the
    # interior (identical to (psi_{i+1} - psi_{i-1})/(2h) on this uniform mesh)
    # and second-order one-sided at the two end cells.
    dpsi = np.gradient(psi, x, edge_order=2) if n >= 3 else None
    band = np.abs(psi) <= p["nLayers"] * dx
    geo_band = (nearest <= p["nLayers"] * dx) if nearest is not None \
        else np.zeros(n, dtype=bool)

    if dpsi is not None:
        row["meanDpsiDxAll"] = float(dpsi @ dx / dx.sum())
        row["meanAbsDpsiDxAll"] = float(np.abs(dpsi) @ dx / dx.sum())
        if band.any():
            row["meanDpsiDxBand"] = float(dpsi[band] @ dx[band] / dx[band].sum())
            row["meanAbsDpsiDxBand"] = float(
                np.abs(dpsi[band]) @ dx[band] / dx[band].sum())
    row["exactDpsiDx"] = st["g"]        # magnitude; the sign is s_exact
    row["nBandCells"] = int(band.sum())
    row["nGeomBandCells"] = int(geo_band.sum())

    if dpsi is not None and s_exact is not None:
        gerr = dpsi - st["g"] * s_exact
        # Slab only: drop the cells whose difference stencil straddles the
        # medial axis. The EXACT derivative jumps from -g to +g there, so such a
        # stencil carries an O(1) error belonging to the kink in the exact
        # solution, not to the scheme, and it would otherwise pin Linf at ~g
        # forever and hide whatever the scheme is really doing.
        ok = np.ones(n, dtype=bool)
        if p["geometry"] == "slab" and n >= 3:
            lo, hi = np.empty(n), np.empty(n)
            lo[1:-1], hi[1:-1] = x[:-2], x[2:]
            lo[0], hi[0] = x[0], x[2]
            lo[-1], hi[-1] = x[-3], x[-1]
            ok = ~((lo <= st["xm"]) & (st["xm"] <= hi))
        row["nKinkExcluded"] = int((~ok).sum())
        for tag, mask in (("All", ok),
                          ("Band", band & ok),
                          ("GeomBand", geo_band & ok)):
            l1, l2, li = norms(gerr[mask], dx[mask])
            row[f"gradErrL1{tag}"] = "" if l1 is None else l1
            row[f"gradErrL2{tag}"] = "" if l2 is None else l2
            row[f"gradErrLinf{tag}"] = "" if li is None else li

    # ---- 2. field error --------------------------------------------------- #
    if d_exact is not None:
        perr = psi - st["g"] * d_exact
        for tag, mask in (("All", np.ones(n, dtype=bool)), ("Band", band)):
            l1, l2, li = norms(perr[mask], dx[mask])
            row[f"psiErrL1{tag}"] = "" if l1 is None else l1
            row[f"psiErrL2{tag}"] = "" if l2 is None else l2
            row[f"psiErrLinf{tag}"] = "" if li is None else li

    # ---- 3. interface position: the 1D shape error ------------------------ #
    cr = zero_crossings(x, psi)
    row["nZeroCrossings"] = len(cr)
    if p["geometry"] == "slab":
        targets = (("xL", st["xL"]), ("xR", st["xR"]))
    else:
        targets = (("xi", st["xi"]),)
    used = set()
    for name, exact in targets:
        if exact is None:
            continue
        row[f"{name}Exact"] = exact
        # Pick the crossing NEAREST the exact position, never simply the first:
        # a spurious extra crossing from an overshoot elsewhere in the domain
        # must not be reported as the interface.
        cand = [(abs(xc - exact), k, xc) for k, (_i, xc) in enumerate(cr)
                if k not in used]
        if not cand:
            continue
        _d, k, xc = min(cand)
        used.add(k)
        row[f"{name}Measured"] = xc
        row[f"{name}Error"] = xc - exact
        row[f"{name}ErrorCells"] = (xc - exact) / h

    # ---- 4. phase volume: the measure of {psi < 0} ------------------------ #
    recon = negative_measure(x, psi, xlo, xhi, cr)
    row["phaseMeasureRecon"] = recon
    alpha = read_field(os.path.join(tdir, "alpha"), n, f"alpha at t = {t:g}")
    meas, src = recon, "reconstruction"
    if alpha is not None and alpha.ndim == 1:
        alpha = alpha[perm]
        # alpha = 1 where psi < 0 (detrixheAslamPhaseIndicator.C). Verified per
        # run rather than assumed, because a differently-signed indicator would
        # turn the phase measure into its complement without any other symptom.
        neg, pos = psi < 0, psi > 0
        if neg.any() and pos.any() and alpha[pos].mean() > alpha[neg].mean():
            if not ctx["alphaWarned"]:
                print("  [alpha] the written alpha is LARGER where psi > 0: the "
                      "indicator marks the psi > 0 phase, so 1 - alpha is "
                      "integrated for the measure of {psi < 0}.")
                ctx["alphaWarned"] = True
            alpha = 1.0 - alpha
        row["phaseMeasureAlpha"] = float(alpha @ dx)
        meas, src = row["phaseMeasureAlpha"], "alpha"
    row["phaseMeasure"], row["phaseMeasureSource"] = meas, src

    exact_measure = None
    if p["geometry"] == "slab" and st["xL"] is not None:
        # 2 W(t) = (x_R0 - x_L0) e^{a t}, from the characteristics -- the slab is
        # advected AND stretched. Clipped to the domain so a slab that has partly
        # left the box is not compared against a length the mesh cannot hold.
        exact_measure = max(0.0, min(st["xR"], xhi) - max(st["xL"], xlo))
    elif st["xi"] is not None:
        # Single interface: {psi < 0} is the part of the domain on the negative
        # side of x_i(t), i.e. the distance from the domain end to the interface.
        exact_measure = (min(max(st["xi"], xlo), xhi) - xlo) if p["orientation"] > 0 \
            else (xhi - min(max(st["xi"], xlo), xhi))
    if exact_measure is not None:
        row["phaseMeasureExact"] = exact_measure
        row["phaseMeasureAbsError"] = meas - exact_measure
        row["phaseMeasureRelError"] = ((meas - exact_measure) / exact_measure
                                       if exact_measure > 0 else "")

    # ---- 5. velocity extension: the closestPoint diagnosis ---------------- #
    uext = read_field(os.path.join(tdir, "Uext"), n, f"Uext at t = {t:g}")
    row["uextPresent"] = int(uext is not None)
    if uext is None or uext.ndim != 2 or uext.shape[1] != 3 or s_exact is None:
        return row
    ue = uext[perm]
    exact_u = np.zeros((n, 3))
    if p["geometry"] == "slab":
        uL, uR = p["a"] * (st["xL"] - p["x0"]), p["a"] * (st["xR"] - p["x0"])
        exact_u[:, 0] = np.where(x < st["xm"], uL, uR)
        row["uextExactKind"], row["xMedial"] = "twoSided", st["xm"]
    else:
        uL = uR = None
        exact_u[:, 0] = p["a"] * (st["xi"] - p["x0"])
        row["uextExactKind"] = "uniform"
    # FULL vector norm: a spurious y/z component of the extension is a defect
    # the x-component alone would not show.
    uerr = np.linalg.norm(ue - exact_u, axis=1)

    l2, li = None, None
    if band.any():
        _l1, l2, li = norms(uerr[band], dx[band])
        row["uextErrL2Band"], row["uextErrLinfBand"] = l2, li
    row["nUextBandCells"] = int(band.sum())

    if band.any():
        idx = np.flatnonzero(band)[int(np.argmax(uerr[band]))]
        row["uextArgmaxX"] = float(x[idx])
        if st["xm"] is not None:
            row["uextArgmaxCellsFromMedial"] = float((x[idx] - st["xm"]) / h)
        if st["xi"] is not None:
            row["uextArgmaxCellsFromInterface"] = float((x[idx] - st["xi"]) / h)
        elif st["xL"] is not None:
            near = st["xL"] if abs(x[idx] - st["xL"]) <= abs(x[idx] - st["xR"]) \
                else st["xR"]
            row["uextArgmaxCellsFromInterface"] = float((x[idx] - near) / h)

    if p["geometry"] == "slab" and band.any():
        # Which interface did the extension ACTUALLY reach for? The exact
        # closest-point extension takes only the two values uL and uR, so the
        # nearer of the two names the foot point the solver used. Cells whose
        # inferred foot point is on the other side of the medial axis from
        # themselves are the misassignment this case exists to expose.
        inferred = np.where(np.abs(ue[:, 0] - uL) <= np.abs(ue[:, 0] - uR), -1.0, 1.0)
        geom_side = np.where(x < st["xm"], -1.0, 1.0)
        same = band & (inferred == geom_side)
        cross = band & (inferred != geom_side)
        for tag, mask in (("SameSide", same), ("CrossSide", cross)):
            _l1, sl2, sli = norms(uerr[mask], dx[mask])
            row[f"uextErrL2{tag}"] = "" if sl2 is None else sl2
            row[f"uextErrLinf{tag}"] = "" if sli is None else sli
        row["nUextSameSide"], row["nUextCrossSide"] = int(same.sum()), int(cross.sum())
        # The nearest-value classification is only meaningful while the error is
        # smaller than half the jump it is classifying against.
        jump = abs(uR - uL)
        row["uextSplitReliable"] = int(li is not None and jump > 0 and li < 0.5 * jump)
    return row


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _fmt(v):
    if v is None or v == "":
        return ""
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def write_csv(path, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: _fmt(r.get(c, "")) for c in COLUMNS})


def _col(row, key, fmt="%11.4e", width=11):
    v = row.get(key, "")
    if v == "" or v is None:
        return " " * (width - 3) + "n/a"
    try:
        return fmt % float(v)
    except (TypeError, ValueError):
        return str(v).rjust(width)


def summarize(rows, p, ctx):
    """The verification, readable at a glance: measured gradient vs e^{-a t} vs 1."""
    print()
    print("=" * 118)
    print("1D UNIAXIAL STRETCH -- measured against the closed-form solution")
    print("=" * 118)
    form = {"R": "(2)  psi = x - x_i(t),  d(psi)/dx = 1 exactly",
            "none": "(1)  psi = x0 + (x - x0)e^{-a t} - xi0,  d(psi)/dx = e^{-a t}"}
    print(f"  exact form scored : {form.get(p['exactForm'], p['exactForm'])}")
    print(f"  v(x) = a (x - x0) e_x with a = {p['a']:g} 1/s, x0 = {p['x0']:g} m; "
          f"div v = a = {p['a']:g} 1/s (NOT divergence-free -- intrinsic to a 1D "
          f"stretch, and the solver's advective form is exact for it)")
    print()
    hdr = (f"{'t [s]':>9} {'mean|dpsi/dx|':>14} {'exact':>9} {'e^{-a t}':>9} "
           f"{'grad L2 band':>13} {'psi L2 band':>12} {'x_i err [h]':>12} "
           f"{'vol rel err':>12} {'Uext Linf':>11} {'band':>5}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        xerr = r.get("xiErrorCells", "")
        if xerr == "":
            xerr = r.get("xLErrorCells", "")
        print(f"{float(r['time']):>9.4g} "
              f"{_col(r, 'meanAbsDpsiDxBand', '%14.8f', 14)} "
              f"{_col(r, 'exactDpsiDx', '%9.6f', 9)} "
              f"{_col(r, 'decayExpMinusAT', '%9.6f', 9)} "
              f"{_col(r, 'gradErrL2Band', '%13.4e', 13)} "
              f"{_col(r, 'psiErrL2Band', '%12.4e', 12)} "
              f"{(('%12.4e' % float(xerr)) if xerr not in ('', None) else '         n/a')} "
              f"{_col(r, 'phaseMeasureRelError', '%12.4e', 12)} "
              f"{_col(r, 'uextErrLinfBand', '%11.4e', 11)} "
              f"{_col(r, 'nBandCells', '%5d', 5)}")
    print("-" * len(hdr))
    print("  GRADIENT, SHAPE and VOLUME are all above and are reported together by "
          "construction:")
    print("    * mean|dpsi/dx| against BOTH references makes the verification direct --")
    print("      the R arm must hold 1 to discretization error, the sourceless arm must")
    print("      track e^{-a t} (0.36787944 at a = 1 1/s, t = 1 s).")
    print("    * x_i err [h] is the shape error in cell widths; the interface must track")
    print("      x_i(t) = x0 + (xi0 - x0) e^{a t} on BOTH arms -- losing the distance")
    print("      property does not move the zero contour.")
    print("    * vol rel err is the measure of {psi < 0} against its exact value.")
    if p["geometry"] == "slab":
        W0, nl, h = p["W0"], p["nLayers"], ctx["h"]
        print(f"  SLAB: W0 = {W0:g} m, nLayers h = {nl * h:g} m, "
              f"W0/(nLayers h) = {W0 / (nl * h):.3f} -- the medial axis is inside the "
              f"extension band when this is < 1")
        cross = [int(r["nUextCrossSide"]) for r in rows
                 if str(r.get("nUextCrossSide", "")) != ""]
        if cross:
            print(f"  closestPoint foot-point misassignment: max {max(cross)} band "
                  f"cell(s) whose extension velocity matches the interface on the "
                  f"FAR side of the medial axis")
    else:
        print("  SINGLE INTERFACE: the exact closest-point extension is UNIFORM, "
              "a (x_i(t) - x0), so n.grad(Uext) = 0 exactly and there is no medial "
              "axis in the domain -- the Uext split columns are empty by construction.")
    if p["exactForm"] == "none" and p["sourceType"] in ("Rdiv", "RdivStrictSp"):
        print(f"  NOTE: sourceType = {p['sourceType']} is scored against closed form (1). "
              "That is the PREDICTION under test (in 1D w = v, phiW = phi, the two "
              "fvm::div terms cancel discretely and the source assembles the zero "
              "matrix), not an established result. Read these rows as the measurement "
              "of that prediction.")


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", help="1D stretch case directory (holds 0/, system/, ...)")
    ap.add_argument("--source", choices=("auto", "none", "R"), default="auto",
                    help="which closed form is exact: (1) sourceless = none, "
                         "(2) SDPLS R = R. Default auto: read levelSet/sdplsSource/type.")
    ap.add_argument("--geometry", choices=("auto", "single", "slab"), default="auto",
                    help="default auto: one zero crossing of 0/psi = single "
                         "interface, two = slab")
    ap.add_argument("--a", type=float, default=None,
                    help="strain rate [1/s]; override, normally read from the case")
    ap.add_argument("--x0", type=float, default=None,
                    help="stagnation point x [m]; override")
    ap.add_argument("--xi0", type=float, default=None,
                    help="initial interface x [m] (slab: its CENTRE); override")
    ap.add_argument("--half-width", type=float, default=None,
                    help="slab half-width W0 [m]; override")
    ap.add_argument("--nlayers", type=float, default=None,
                    help="narrow-band width in cell layers; override")
    ap.add_argument("--domain", default=None,
                    help="xmin,xmax fallback when constant/polyMesh/points is "
                         "unreadable (assumes a uniform mesh)")
    ap.add_argument("--out", default=None,
                    help="CSV path (default <case>/analyse1Dstretch.csv)")
    args = ap.parse_args(argv)

    case = os.path.abspath(args.case)
    if not os.path.isdir(case):
        raise SystemExit(f"[analyse_1Dstretch] no such case directory: {case}")
    print(f"[analyse_1Dstretch] {case}")

    tdirs = time_dirs(case)
    if not tdirs:
        if glob.glob(os.path.join(case, "processor*")):
            raise SystemExit(
                "[analyse_1Dstretch] only processor*/ directories hold fields here. "
                "This script reads the RECONSTRUCTED case: run `reconstructPar` "
                "first (a 1D stretch case is meant to run serial anyway).")
        raise SystemExit(f"[analyse_1Dstretch] no written time directory with a psi "
                         f"field under {case}")

    psi0_raw = read_field(os.path.join(tdirs[0][1], "psi"), None, "psi")
    if psi0_raw is None:
        raise SystemExit(f"[analyse_1Dstretch] cannot read {tdirs[0][1]}/psi")
    ncells = int(psi0_raw.shape[0])

    x, dx, perm, mesh_why = cell_centres_x(case, tdirs, ncells, args.domain)
    if x is None:
        raise SystemExit(
            "[analyse_1Dstretch] could not establish the cell centres: no written C "
            "field, and constant/polyMesh/points does not give exactly "
            f"{ncells + 1} distinct x planes. Pass --domain xmin,xmax if the mesh is "
            "uniform, or write cell centres with "
            "`postProcess -func writeCellCentres`.")
    psi0 = psi0_raw[perm]
    order_ok, order_why = verify_cell_order(psi0)
    if not order_ok:
        print(f"  [mesh] *** CELL ORDER CHECK FAILED: {order_why}. 0/psi should be "
              f"monotone (single interface) or V-shaped (slab) in ascending x. "
              f"Every number below is suspect; check that this is a 1D column mesh "
              f"and a serial run. ***")

    dicts, tokens = read_dicts(case), _tokens(case)
    u0 = read_field(os.path.join(tdirs[0][1], "U"), ncells, "U")
    if u0 is not None:
        u0 = u0[perm]
    p, why = resolve_parameters(case, args, dicts, tokens, x, psi0, u0)

    h = float(np.mean(dx))
    xlo, xhi = float(x[0] - 0.5 * dx[0]), float(x[-1] + 0.5 * dx[-1])
    ctx = dict(x=x, dx=dx, perm=perm, p=p, h=h, xlo=xlo, xhi=xhi,
               orderOK=int(order_ok), alphaWarned=False)

    print(f"  mesh              : {ncells} cells, x in [{xlo:g}, {xhi:g}] m, "
          f"h = {h:g} m   [{mesh_why}]")
    print(f"  cell order        : {order_why}")
    print(f"  strain rate a     : {p['a']:g} 1/s   [{why['a']}]")
    print(f"  stagnation x0     : {p['x0']:g} m   [{why['x0']}]")
    if p["aFit"] is not None:
        # The fit is what the solver ACTUALLY advected with; a disagreement with
        # the dictionary means the case is not the case its dictionary describes.
        da = abs(p["aFit"] - p["a"]) / max(abs(p["a"]), 1e-300)
        dx0 = (abs(p["x0Fit"] - p["x0"]) / max(abs(xhi - xlo), 1e-300)
               if p["x0Fit"] is not None else 0.0)
        tag = "consistent" if (da < 1e-6 and dx0 < 1e-6) else "*** DISAGREES ***"
        print(f"  U-field fit       : a = {p['aFit']:.10g} 1/s, "
              f"x0 = {'n/a' if p['x0Fit'] is None else format(p['x0Fit'], '.10g')} m "
              f"-> {tag} (relative: a {da:.2e}, x0 {dx0:.2e} of the domain)")
    else:
        print("  U-field fit       : no readable U field -- the strain rate and "
              "stagnation point could not be cross-checked against what ran")
    if p["geometry"] == "slab":
        print(f"  geometry          : slab, x_L0 = {p['xL0']:g} m, "
              f"x_R0 = {p['xR0']:g} m, W0 = {p['W0']:g} m   [{why['geometry']}; "
              f"{why['slab']}]")
    else:
        xi0 = p["xi0"]
        print(f"  geometry          : single interface, xi0 = "
              f"{'NOT RESOLVED' if xi0 is None else format(xi0, 'g') + ' m'}   "
              f"[{why['geometry']}; {why['xi0']}]")
    print(f"  narrow band       : nLayers = {p['nLayers']:g} -> "
          f"|psi| <= {p['nLayers'] * h:g} m   [{why['nLayers']}]")
    print(f"  sdpls source      : {p['sourceType'] or 'not recorded'} -> exact form "
          f"'{p['exactForm']}'   [{why['exactForm']}]")
    print(f"  velocityExtension : {p['velocityExtension'] or 'not recorded'}")
    print(f"  written times     : {len(tdirs)}")

    rows = [evaluate(t, d, ctx) for t, d in tdirs]
    out = args.out or os.path.join(case, "analyse1Dstretch.csv")
    write_csv(out, rows)
    summarize(rows, p, ctx)
    print(f"\n  wrote {out}  ({len(rows)} row(s), {len(COLUMNS)} columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
