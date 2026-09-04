#!/usr/bin/env python3
"""leia_refine -- shared helpers for the static interface-refinement drivers.

Used by leiaRefineHexMesh.py and leiaRefinePolyMesh.py (and by
check_refined_band.py to re-check an existing case). Standard library only: the
drivers run inside the Snakemake `mesh` rule, which under profiles/slurm is a
serial sbatch job on a compute node with no conda environment.

WHAT THE DRIVERS PRODUCE. A statically refined mesh in constant/polyMesh AND the
initial fields on it, RE-INITIALISED by leiaSetFields on the final mesh: 0/ is
reset from 0.org before every leiaSetFields call and once more after the last
refinement pass, so no field that was mapped through a refinement
(refineHexMesh maps every field it finds in 0/) is ever used as a criterion or
kept as an initial condition. A mapped volume fraction is a smeared one.

THE BAND CHECK (band_check) reads the FINAL fields and the FINAL mesh and asks
the two questions that decide whether the refined mesh is fit for the solver:
  * is every interface cell (0 < alpha < 1) inside the fine region?  The
    capillary force, its curvature and the level-set transport must never see a
    resolution transition (see the layer table in the plan / the article).
  * how many fine cells separate the interface from the first coarse cell
    centre?  min|psi| over coarse cells, in units of the fine spacing. FAIL
    below 4 (the stencil minimum), WARN below the requested band width.
  * hex only: does N_CELLS (the capillary-dt handle) equal the FINE spacing?
    controlDict has adjustTimeStep no, so a wrong N_CELLS is a wrong dt for the
    whole run.  poly: the same question as a suggested N_CELLS pin.

Nothing here reads a solver log with an ad-hoc grep: cell counts come from the
owner-file header, set sizes from the set file, and checkMesh's own summary
lines are the only log lines parsed.
"""
import csv
import glob
import json
import math
import os
import re
import shutil
import subprocess

TAG = "[leia_refine]"

# The two bounds that make a cell an interface cell. 1e-9 rather than 0 so a
# round-off alpha of 1e-16 in a bulk cell is not "interface".
ALPHA_LO, ALPHA_HI = 1e-9, 1.0 - 1e-9

# The stencil minimum: complete fine layers between the interface and the first
# coarse cell, at the worst point. The first coarse cell centre lies one fine
# cell inside that cell, so layers = min|psi|/h_fine - 1 over coarse cells.
MIN_FINE_LAYERS = 4.0


def say(msg):
    print(f"{TAG} {msg}", flush=True)


def tail(path, n=25):
    try:
        with open(path, errors="replace") as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return ""


def run(cmd, logfile, cwd=".", dry_run=False):
    """Run a shell command in cwd, stdout+stderr into cwd/logfile, raise on rc != 0."""
    say(f"$ {cmd}   > {logfile}")
    if dry_run:
        return
    path = os.path.join(cwd, logfile)
    with open(path, "w") as fh:
        rc = subprocess.call(cmd, shell=True, cwd=cwd, stdout=fh,
                             stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"{cmd!r} exited {rc}; tail of {logfile}:\n{tail(path)}")


def reset_zero(case, dry_run=False):
    """0/ := 0.org with the *.template files dropped.

    Called before EVERY leiaSetFields and after the last pass, so that a field
    mapped by refineHexMesh (it maps all of 0/) or left by a previous pass is
    never used as a criterion nor kept as an initial condition."""
    say("0/ := 0.org (drop every mapped or stale field)")
    if dry_run:
        return
    zero, org = os.path.join(case, "0"), os.path.join(case, "0.org")
    if not os.path.isdir(org):
        raise RuntimeError(f"{org} is missing; the drivers need a 0.org to reset from")
    shutil.rmtree(zero, ignore_errors=True)
    shutil.copytree(org, zero)
    for root, _dirs, files in os.walk(zero):
        for f in files:
            if f.endswith(".template"):
                os.remove(os.path.join(root, f))


def tokens(case):
    """The materialized tokens of a rendered case ({} when run by hand elsewhere)."""
    path = os.path.join(case, "case_params.json")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return json.load(fh).get("tokens", {})


def tok_float(toks, key, default=None):
    try:
        return float(toks[key])
    except (KeyError, TypeError, ValueError):
        return default


def tok_int(toks, key, default=None):
    v = tok_float(toks, key)
    return default if v is None else int(round(v))


def alpha_name_from(setfields_cmd, default="alpha.water"):
    m = re.search(r"-alphaName\s+(\S+)", setfields_cmd or "")
    return m.group(1) if m else default


# --------------------------------------------------------------------------
# Mesh and field readers (ASCII OpenFOAM files; the cases write ascii).
# --------------------------------------------------------------------------
_NCELLS = re.compile(r"nCells:\s*(\d+)")
_FORMAT = re.compile(r"\bformat\s+(\w+)\s*;")
_UNIFORM_LIST = re.compile(r"(\d+)\s*\{\s*(-?[\d.eE+-]+)\s*\}")


def n_cells(case):
    """Cell count of the CURRENT mesh, never from a log.

    blockMesh and refineHexMesh write `note "nCells:..."` into the owner header;
    cfMesh's pMesh does not. The fallback reads the owner list itself and returns
    max(owner) + 1, which is exact for any mesh in which every cell owns a face (all
    OpenFOAM meshes) and cannot be confused by a stale field of another mesh."""
    path = os.path.join(case, "constant", "polyMesh", "owner")
    with open(path, "rb") as fh:
        head = fh.read(4096).decode(errors="replace")
    m = _NCELLS.search(head)
    if m:
        return int(m.group(1))
    with open(path) as fh:
        text = fh.read()
    _assert_ascii(text, path)
    body = _strip_comments(text[text.index("}") + 1:])
    m = re.search(r"(\d+)\s*\(", body)
    if not m:
        raise RuntimeError(f"cannot read the owner list in {path}")
    return max(int(v) for v in _list_after(body, m.end() - 1)) + 1


def _strip_comments(text):
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def boundary_types(case):
    """{patchName: type} from constant/polyMesh/boundary (the authority on patches)."""
    with open(os.path.join(case, "constant", "polyMesh", "boundary")) as fh:
        text = _strip_comments(fh.read())
    return dict(re.findall(r"(\w+)\s*\{\s*type\s+(\w+)\s*;", text))


def _assert_ascii(text, path):
    m = _FORMAT.search(text[:2048])
    if m and m.group(1) != "ascii":
        raise RuntimeError(f"{path} is {m.group(1)}; the readers handle ascii only "
                           "(controlDict writeFormat ascii)")


def _list_after(text, start):
    """Numbers of the parenthesised list that starts at text[start] == '('."""
    end = text.index(")", start)
    return text[start + 1:end].split()


def read_scalar_field(path, ncells):
    """internalField of an ascii volScalarField as a list of ncells floats."""
    with open(path) as fh:
        text = fh.read()
    _assert_ascii(text, path)
    m = re.search(r"internalField\s+uniform\s+(-?[\d.eE+-]+)\s*;", text)
    if m:
        return [float(m.group(1))] * ncells
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", text)
    if not m:
        raise RuntimeError(f"no scalar internalField in {path}")
    n = int(m.group(1))
    vals = [float(v) for v in _list_after(text, m.end() - 1)]
    if n != ncells or len(vals) != ncells:
        raise RuntimeError(f"{path}: {len(vals)} values for {ncells} cells")
    return vals


def read_label_list(path, ncells):
    """An ascii labelIOList (e.g. constant/polyMesh/cellLevel), uniform form included."""
    with open(path) as fh:
        text = fh.read()
    _assert_ascii(text, path)
    body = text[text.index("}") + 1:] if "FoamFile" in text else text
    body = _strip_comments(body)
    m = _UNIFORM_LIST.search(body)
    if m and int(m.group(1)) == ncells:
        return [int(float(m.group(2)))] * ncells
    m = re.search(r"(\d+)\s*\(", body)
    if not m:
        raise RuntimeError(f"no label list in {path}")
    vals = [int(v) for v in _list_after(body, m.end() - 1)]
    if int(m.group(1)) != ncells or len(vals) != ncells:
        raise RuntimeError(f"{path}: {len(vals)} labels for {ncells} cells")
    return vals


def read_set(case, name):
    """Labels of a written cellSet/faceSet (constant/polyMesh/sets/<name>)."""
    path = os.path.join(case, "constant", "polyMesh", "sets", name)
    with open(path) as fh:
        text = fh.read()
    body = _strip_comments(text[text.index("}") + 1:])
    m = re.search(r"(\d+)\s*\(", body)
    if not m:
        return []
    labels = [int(v) for v in _list_after(body, m.end() - 1)]
    if len(labels) != int(m.group(1)):
        raise RuntimeError(f"{path}: {len(labels)} labels announced as {m.group(1)}")
    return labels


def worst_case_layers_after_split(case, selected, ncells):
    """Complete fine layers the SELECTED cells will leave between the interface and
    the first unselected cell once every selected cell is split in two per
    direction: min over unselected cells c of 2|psi_c|/h_c - 1, with h_c = V_c^(1/3)
    of that (coarse) cell. Needs 0/psi and 0/V on the CURRENT mesh. This is the
    measurement that decides how many face dilations a pass needs: face dilation
    grows a Manhattan diamond, so along a sphere's diagonals each dilation buys only
    ~1.15 fine cells instead of 2, and no closed-form count survives a change of
    shape or mesh."""
    psi = read_scalar_field(os.path.join(case, "0", "psi"), ncells)
    vol = read_scalar_field(os.path.join(case, "0", "V"), ncells)
    sel = bytearray(ncells)
    for l in selected:
        sel[l] = 1
    worst = float("inf")
    for c in range(ncells):
        if not sel[c]:
            d = 2.0 * abs(psi[c]) / vol[c] ** (1.0 / 3.0) - 1.0
            if d < worst:
                worst = d
    return worst


def set_size(case, name, kind="cellSet"):
    """Size of a written topoSet from its file, not from the topoSet log."""
    path = os.path.join(case, "constant", "polyMesh", "sets", name)
    with open(path) as fh:
        text = fh.read()
    body = _strip_comments(text[text.index("}") + 1:])
    m = re.search(r"(\d+)\s*[\(\{]", body)
    if not m:
        raise RuntimeError(f"cannot read the size of {kind} {name} from {path}")
    return int(m.group(1))


# --------------------------------------------------------------------------
# checkMesh summary (the only log lines any of this parses).
# --------------------------------------------------------------------------
def checkmesh_stats(path):
    out = {}
    try:
        with open(path, errors="replace") as fh:
            text = fh.read()
    except OSError:
        return out
    for key, rx in (("nCellsCheckMesh", r"^\s*cells:\s+(\d+)"),
                    ("nHexahedra", r"^\s*hexahedra:\s+(\d+)"),
                    ("nPolyhedra", r"^\s*polyhedra:\s+(\d+)"),
                    ("maxNonOrtho", r"(?:Mesh non-orthogonality Max:|Max non-orthogonality =)\s*([\d.eE+-]+)"),
                    ("maxSkewness", r"Max skewness = ([\d.eE+-]+)")):
        m = re.search(rx, text, flags=re.M)
        if m:
            out[key] = m.group(1)
    out["checkMeshFailed"] = "Failed" if re.search(r"^Failed \d+ mesh checks", text, re.M) else ""
    return out


# --------------------------------------------------------------------------
# The band check.
# --------------------------------------------------------------------------
def _median(sorted_vals):
    n = len(sorted_vals)
    return sorted_vals[n // 2] if n % 2 else 0.5 * (sorted_vals[n // 2 - 1] + sorted_vals[n // 2])


def band_check(case, mode, band_cells, alpha_name="alpha.water", toks=None,
               level_ratio=2.0, pin_tolerance=None, allow_pin_mismatch=False):
    """Classify the final mesh into fine/coarse and measure the band around the
    interface. Needs 0/V (postProcess -func writeCellVolumes -time 0), 0/psi and
    0/<alpha> written by leiaSetFields on the FINAL mesh, and for mode "hex"
    constant/polyMesh/cellLevel (written by refineHexMesh).

    Returns a flat dict (CSV row) whose "status" is PASS, WARN or FAIL."""
    toks = toks or {}
    n = n_cells(case)
    vol = read_scalar_field(os.path.join(case, "0", "V"), n)
    psi = read_scalar_field(os.path.join(case, "0", "psi"), n)
    alpha = read_scalar_field(os.path.join(case, "0", alpha_name), n)
    size = [v ** (1.0 / 3.0) for v in vol]
    band = [ALPHA_LO < a < ALPHA_HI for a in alpha]
    n_band = sum(band)
    if n_band == 0:
        raise RuntimeError(f"no interface cells ({ALPHA_LO} < {alpha_name} < {ALPHA_HI})")
    band_sizes = sorted(size[i] for i in range(n) if band[i])
    h_band_med, h_band_min, h_band_max = _median(band_sizes), band_sizes[0], band_sizes[-1]

    res = {"mode": mode, "nCells": n, "nBand": n_band, "bandCellsRequested": band_cells,
           "hBandMedian": h_band_med, "hBandMin": h_band_min, "hBandMax": h_band_max}
    reasons = []

    if mode == "hex":
        lev = read_label_list(os.path.join(case, "constant", "polyMesh", "cellLevel"), n)
        top = max(lev)
        fine = [l == top for l in lev]
        top_sizes = sorted(size[i] for i in range(n) if fine[i])
        h_fine = _median(top_sizes)
        counts = {}
        for l in lev:
            counts[l] = counts.get(l, 0) + 1
        res["topLevel"] = top
        res["fracPerLevel"] = " ".join(f"L{l}:{counts[l] / n:.4f}" for l in sorted(counts))
        # cross-check: a top-level cell whose size is off by > 5 % is not a hex
        if top_sizes and (top_sizes[-1] / top_sizes[0] > 1.05):
            reasons.append(f"top-level cell sizes spread {top_sizes[0]:.3e}..{top_sizes[-1]:.3e}")
    elif mode == "poly":
        # cfMesh writes no levels: the fine region is the set of cells within
        # sqrt(level_ratio) of the median interface-cell size -- midway, in log
        # scale, between one octree level and the next.
        h_fine = h_band_med
        thr = h_fine * math.sqrt(level_ratio)
        fine = [s <= thr for s in size]
        res["fineThreshold"] = thr
    else:
        raise ValueError(f"mode must be hex or poly, got {mode!r}")

    n_fine = sum(fine)
    n_band_outside = sum(1 for i in range(n) if band[i] and not fine[i])
    coarse_psi = [abs(psi[i]) for i in range(n) if not fine[i]]
    min_psi_over_h = (min(coarse_psi) / h_fine) if coarse_psi else float("inf")
    fine_layers = min_psi_over_h - 1.0     # complete fine layers at the worst point
    # how far the fine region reaches from the interface (its widest point), in fine
    # cells -- the cost side of the band: hex ~ requested + 1, cfMesh grades further
    fine_psi = [abs(psi[i]) for i in range(n) if fine[i]]
    res.update({"hFine": h_fine, "nFine": n_fine, "fracFine": n_fine / n,
                "nBandOutsideFine": n_band_outside, "minPsiOverHfine": min_psi_over_h,
                "fineLayersWorst": fine_layers,
                "fineExtentOverHfine": (max(fine_psi) / h_fine) if fine_psi else 0.0})
    # Where the grading ends, without a size threshold: median cell size in bands of
    # distance from the interface (in fine cells). On a cfMesh mesh the dual cells of one
    # octree level scatter in size, so a single fine/coarse split over-counts "fine"
    # cells far away; this profile is the honest picture for both mesh kinds.
    edges = [0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 1e9]
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        vals = sorted(size[i] for i in range(n) if lo <= abs(psi[i]) / h_fine < hi)
        if vals:
            tag = f"{lo:g}-{hi:g}" if hi < 1e8 else f"{lo:g}+"
            prof.append(f"{tag}:{_median(vals) / h_fine:.2f}")
    res["sizeProfileOverHfine"] = " ".join(prof)

    # The dt handle: N_CELLS must encode the FINE spacing (hex) / be pinned to
    # the measured fine spacing (poly).
    length = tok_float(toks, "DOMAIN_LENGTH")
    n_tok = tok_float(toks, "N_CELLS")
    if length and length > 0:
        res["N_CELLS_suggested"] = int(round(length / h_fine))
        res["N_CELLS_suggested_fromMin"] = int(round(length / h_band_min))
        if n_tok and n_tok > 0:
            h_tok = length / n_tok
            rel = abs(h_tok - h_fine) / h_fine
            res["N_CELLS"] = int(n_tok)
            res["hFromN_CELLS"] = h_tok
            res["pinRelError"] = rel
            tol = pin_tolerance if pin_tolerance is not None else (0.02 if mode == "hex" else 0.05)
            if rel > tol and not allow_pin_mismatch:
                reasons.append(f"N_CELLS={int(n_tok)} gives h={h_tok:.4e} but the fine "
                               f"spacing is {h_fine:.4e} ({rel * 100:.1f} % off): the "
                               f"capillary dt is wrong; pin N_CELLS={res['N_CELLS_suggested']}")

    if n_band_outside:
        reasons.append(f"{n_band_outside} interface cells outside the fine region")
    if fine_layers < MIN_FINE_LAYERS:
        reasons.append(f"only {fine_layers:.2f} complete fine layers beyond the interface at "
                       f"the worst point < the stencil minimum {MIN_FINE_LAYERS:g}")
    status = "FAIL" if reasons else ("WARN" if fine_layers < band_cells else "PASS")
    if status == "WARN":
        reasons.append(f"{fine_layers:.2f} complete fine layers at the worst point "
                       f"< requested {band_cells}")
    res["status"] = status
    res["reason"] = "; ".join(reasons)
    return res


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def report(res):
    say(f"band check: {res['status']}  nCells={res['nCells']}  nBand={res['nBand']}  "
        f"nBandOutsideFine={res['nBandOutsideFine']}  fineLayersWorst={res['fineLayersWorst']:.2f}  "
        f"(first coarse centre at {res['minPsiOverHfine']:.2f} h_fine)  hFine={res['hFine']:.4e}  "
        f"fracFine={res['fracFine']:.4f}")
    if res.get("N_CELLS_suggested") is not None:
        say(f"N_CELLS pin: suggested {res['N_CELLS_suggested']} (from the median fine spacing; "
            f"{res['N_CELLS_suggested_fromMin']} from the smallest interface cell)"
            + (f"; case has N_CELLS={res['N_CELLS']} ({res['pinRelError'] * 100:.2f} % off)"
               if "N_CELLS" in res else ""))
    if res.get("reason"):
        say(f"reason: {res['reason']}")


def ensure_cell_volumes(case, dry_run=False):
    if dry_run or not os.path.isfile(os.path.join(case, "0", "V")):
        run("postProcess -func writeCellVolumes -time 0", "log.writeCellVolumes", case, dry_run)


def find_stl(case, func="isoInterface"):
    return sorted(glob.glob(os.path.join(case, "postProcessing", func, "**", "*.stl"),
                            recursive=True))
