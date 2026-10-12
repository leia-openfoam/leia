#!/usr/bin/env python3
"""DAVOF static gate: order table + log-log figures of the Gauss-identity normal.

Reads, for every case of a leiaTestDavofNormal study,

    <case>/case_params.json               mesh, tokens.N_CELLS, DOMAIN_LENGTH, DAVOF_ALPHA_SOURCE
    <case>/leiaTestDavofNormalModels.csv  tidy: one row per MODEL (davof, plicRDF, ...)
    <case>/leiaTestDavofNormal.csv        the DAVOF diagnostics (consistency, closure)
    <case>/leiaTestDavofCurvature.csv     tidy: one row per curvature model (optional)

and writes into the theme's data source (paths.tables_dir/figs_dir):

    tables/davof_normal_<study>.csv           the long table, every (alphaSource, model, N)
    tables/davof_normal_orders_<study>.csv    least-squares + pairwise orders per series
    tables/davof_normal_orders_<study>.tex    booktabs body of the same (\\input-able)
    tables/davof_normal_proposal_<study>.tex  the compact table of the DAVOF proposal
    figures/davof_normal_<metric>_<study>.png log-log error vs h, one per metric

and, when the study is one of the mesh-class ladders of a surface
(<surface>Normal3D[_hex|_perturbed|_poly], 2026-10-12), the COMBINED table of
the three mesh classes from the sibling study directories that exist (a
missing one is skipped with a note):

    tables/davof_normal_proposal_<surface>Normal3D_meshes.tex  mesh column first,
                                              \\midrule between the mesh blocks
    tables/davof_normal_orders_<surface>Normal3D_meshes.csv    the orders, mesh column

Series = (mesh, alphaSource, model); the mesh class (hex | perturbed | poly) is
the "mesh" entry of every case's case_params.json (the study config's `mesh`).
A study has one mesh class, so its own outputs are unchanged by the key (the
hex tables the proposal \\input{}s stay identical). Metrics = E_L1_N, E_L2_N, E_LINF_N (normal error
|n - n_exact| over the interface cells: mean, rms, max), E_L1_N_AW, E_L2_N_AW
(area-weighted), E_AREA_REL (|sum |m_c| - A_exact| / A_exact), E_POS_* (the
plane-centroid distance from the surface), and the curvature of the DAVOF
state: E_KAPPA_* at the interface centroids (|kappa - kappa_exact|, kappa =
kappa_1 + kappa_2), E_KAPPA_REL_L2 = E_KAPPA_L2/KAPPA_REF_L2, E_K_L2 (Gaussian),
E_KAPPA_CELL_L2 (the contour-referenced cell field on the force band),
E_KAPPA_FACE_L2 (interpolate + parallel-surface inverse at the active faces),
E_KAPPA_FACE_FOOT_L2 (the models' foot value). The models CSV carries the
HEADLINE curvature model's columns on the davof row; every curvature model is
a series of its own, model "davof:curv:<CURV_MODEL>", from the curvature CSV.
h = L/N on every mesh class (the polyhedral meshDict asks cfMesh for the
interface cell size h, README "DAVOF"; the mean over all cells is DELTA_X).

--tpf-results <csv> merges the rows of the TwoPhaseFlow benchmark
run/benchmark/reconstruction/sphereNormal3D (results_sphereNormal3D.csv) as
model "tpf:<scheme>/<setAlphaMethod>", alphaSource "tpf", with its LNormalDiff1/2/Inf
and E_area_rel columns mapped onto E_L1_N / E_L2_N / E_LINF_N / E_AREA_REL.

Usage (from the repo root):
    python3 workflow/scripts/make_davof_normal_table.py studies/davof/sphereNormal3D \\
        [--theme davof] [--tpf-results /path/results_sphereNormal3D.csv]
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

CSV_MODELS = "leiaTestDavofNormalModels.csv"
CSV_WIDE = "leiaTestDavofNormal.csv"
CSV_CURV = "leiaTestDavofCurvature.csv"
CURV_METRICS = ["E_KAPPA_L1", "E_KAPPA_L2", "E_KAPPA_LINF", "E_KAPPA_REL_L2", "E_K_L2",
                "E_KAPPA_CELL_L2", "E_KAPPA_FACE_L2", "E_KAPPA_FACE_FOOT_L2"]
METRICS = ["E_L1_N", "E_L2_N", "E_LINF_N", "E_L1_N_AW", "E_L2_N_AW", "E_AREA_REL",
           "E_POS_L1", "E_POS_L2", "E_POS_LINF"] + CURV_METRICS
FIG_METRICS = ["E_L2_N", "E_L1_N", "E_LINF_N", "E_AREA_REL", "E_POS_L2", "E_POS_LINF",
               "E_KAPPA_REL_L2", "E_KAPPA_FACE_L2", "E_KAPPA_CELL_L2", "E_K_L2"]
DIAG = ["MAX_CONSISTENCY", "SUM_M_REL", "MAX_ALPHA_DIFF_DA", "N_WISP"]

# The mesh classes of the workflow (config `mesh:`), in the order of the
# combined table's blocks, and their table labels.
MESH_ORDER = ("hex", "perturbed", "poly")
MESH_LABEL = {"hex": "hexahedral", "perturbed": "perturbed hexahedral (0.1 h)",
              "poly": "polyhedral (cfMesh)"}
# The combined table is built for a study named <surface>Normal3D with one of
# these suffixes (or none: the original hex studies sphereNormal3D and
# ellipsoidNormal3D carry no suffix).
MESH_SUFFIXES = ("_hex", "_perturbed", "_poly")

# Fixed hue per model (categorical slots in the validated order of the dataviz
# reference palette); the line style carries the alpha source. Never cycled.
COLORS = {
    "davof":     "#2a78d6",   # blue
    "plicRDF":   "#eb6834",   # orange
    "gradAlpha": "#1baf7a",   # aqua
    "isoAlpha":  "#eda100",   # yellow
    "davof:curv:normalsOnly": "#2a78d6",   # the headline curvature model, DAVOF's hue
    "davof:curv:hermite":     "#7a3bb5",   # the naive Hermite fit, purple
}
TPF_COLORS = {"plicRDF": "#e87ba4", "gradAlpha": "#008300", "isoAlpha": "#4a3aa7"}
MARKERS = {"davof": "o", "plicRDF": "s", "gradAlpha": "D", "isoAlpha": "^",
           "davof:curv:normalsOnly": "o", "davof:curv:hermite": "v"}
STYLES = {"exactSphere": "-", "quadraticFaces": "--", "detrixheAslam": "-.",
          "linearInterpolant": "-.", "planePhaseIndicator": ":", "tpf": ":"}
SOURCE_LABEL = {"exactSphere": "exact sphere",
                "quadraticFaces": "DA + quadratic faces",
                "detrixheAslam": "Detrixhe-Aslam",
                "linearInterpolant": "Detrixhe-Aslam",   # the name of the first runs
                "planePhaseIndicator": "plane indicator",
                "tpf": "TwoPhaseFlow"}
METRIC_LABEL = {
    "E_L1_N": r"$L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L2_N": r"$L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_LINF_N": r"$L_\infty$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L1_N_AW": r"area-weighted $L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L2_N_AW": r"area-weighted $L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_AREA_REL": r"$|\sum_c |\mathbf{m}_c| - 4\pi R^2| / 4\pi R^2$",
    "E_POS_L1": r"$L_1$ of the plane-centroid distance to the surface [m]",
    "E_POS_L2": r"$L_2$ of the plane-centroid distance to the surface [m]",
    "E_POS_LINF": r"$L_\infty$ of the plane-centroid distance to the surface [m]",
    "E_KAPPA_L1": r"$L_1$ of $|\kappa_c - \kappa_{\mathrm{exact}}|$ at the centroids [1/m]",
    "E_KAPPA_L2": r"$L_2$ of $|\kappa_c - \kappa_{\mathrm{exact}}|$ at the centroids [1/m]",
    "E_KAPPA_LINF": r"$L_\infty$ of $|\kappa_c - \kappa_{\mathrm{exact}}|$ at the centroids [1/m]",
    "E_KAPPA_REL_L2": r"$L_2(\kappa_c - \kappa_{\mathrm{exact}}) / L_2(\kappa_{\mathrm{exact}})$",
    "E_K_L2": r"$L_2$ of $|K_c - K_{\mathrm{exact}}|$ (Gaussian) [1/m$^2$]",
    "E_KAPPA_CELL_L2": r"$L_2$ of the cell field error on the force band [1/m]",
    "E_KAPPA_FACE_L2": r"$L_2$ of $\kappa_f$ (interpolate + inverse) on the active faces [1/m]",
    "E_KAPPA_FACE_FOOT_L2": r"$L_2$ of $\kappa_f$ (foot value) on the active faces [1/m]",
}


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _fit(h, e):
    """Least-squares order p of e ~ h^p over the points with h, e > 0; (p, R2, n)."""
    m = [(a, b) for a, b in zip(h, e) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None, len(m)
    lh, le = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, le, 1)
    r = le - np.polyval(p, lh)
    ss, st = float((r**2).sum()), float(((le - le.mean())**2).sum())
    return float(p[0]), (1.0 - ss/st if st > 0 else 1.0), len(m)


def _pairwise(h, e):
    out = []
    for k in range(len(h) - 1):
        a, b = e[k], e[k + 1]
        if a and b and a > 0 and b > 0 and h[k] > 0 and h[k + 1] > 0 and h[k] != h[k + 1]:
            out.append(np.log(a/b)/np.log(h[k]/h[k + 1]))
        else:
            out.append(None)
    return out


def _read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def collect(study_dir):
    """{(mesh, alphaSource, model): [record, ...]} with record = {N, h, metrics..., diag...}.

    mesh is the "mesh" entry of case_params.json (hex | perturbed | poly; "hex"
    when absent, the studies before 2026-10-12 all carry it)."""
    series = {}
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        cdir = os.path.dirname(meta)
        mpath = os.path.join(cdir, CSV_MODELS)
        if not os.path.isfile(mpath) or os.path.getsize(mpath) == 0:
            continue
        with open(meta) as fh:
            params = json.load(fh)
        tok = params.get("tokens", {})
        mesh = params.get("mesh") or "hex"
        n, L = _f(tok.get("N_CELLS")), _f(tok.get("DOMAIN_LENGTH"))
        if not n or not L:
            continue
        wide = {}
        wpath = os.path.join(cdir, CSV_WIDE)
        if os.path.isfile(wpath) and os.path.getsize(wpath) > 0:
            rows = _read_rows(wpath)
            if rows:
                wide = rows[-1]
        for row in _read_rows(mpath):
            model, src = row.get("MODEL"), row.get("ALPHA_SOURCE") or tok.get("DAVOF_ALPHA_SOURCE")
            if not model or not src:
                continue
            rec = {"N": int(n), "h": L/n, "R_OVER_H": _f(row.get("R_OVER_H")),
                   "DELTA_X": _f(row.get("DELTA_X")),
                   "N_INTERFACE": row.get("N_INTERFACE"), "N_WISP": row.get("N_WISP"),
                   "CPU_SECONDS": _f(row.get("CPU_SECONDS")), "case": os.path.basename(cdir)}
            for m in METRICS:
                rec[m] = _f(row.get(m))
            _relative_curvature(rec, row)
            for d in DIAG:
                rec[d] = (wide.get(d) if model == "davof" else "")
            series.setdefault((mesh, src, model), []).append(rec)
        # Every curvature model as a series of its own (the davof row above
        # carries the headline model's columns only).
        cpath = os.path.join(cdir, CSV_CURV)
        if os.path.isfile(cpath) and os.path.getsize(cpath) > 0:
            for row in _read_rows(cpath):
                cm, src = row.get("CURV_MODEL"), row.get("ALPHA_SOURCE") or tok.get("DAVOF_ALPHA_SOURCE")
                if not cm or not src:
                    continue
                rec = {"N": int(n), "h": L/n, "R_OVER_H": _f(row.get("R_OVER_H")),
                       "DELTA_X": _f(row.get("DELTA_X")),
                       "N_INTERFACE": row.get("N_INTERFACE"), "N_WISP": "",
                       "CPU_SECONDS": _f(row.get("CPU_SECONDS")), "case": os.path.basename(cdir)}
                rec.update({m: None for m in METRICS})
                for m in CURV_METRICS:
                    rec[m] = _f(row.get(m))
                _relative_curvature(rec, row)
                for d in DIAG:
                    rec[d] = ""
                series.setdefault((mesh, src, f"davof:curv:{cm}"), []).append(rec)
    for recs in series.values():
        recs.sort(key=lambda r: -r["h"])          # coarse -> fine
    return series


def _relative_curvature(rec, row):
    """E_KAPPA_REL_L2 = E_KAPPA_L2 / KAPPA_REF_L2 (rms of the exact curvature)."""
    e, ref = _f(row.get("E_KAPPA_L2")), _f(row.get("KAPPA_REF_L2"))
    rec["E_KAPPA_REL_L2"] = (e/ref) if (e is not None and ref) else None
    rec["N_CURV_FALLBACK"] = row.get("N_CURV_FALLBACK", "")


def merge_tpf(series, path, mesh="hex"):
    """TwoPhaseFlow results_sphereNormal3D.csv rows -> series tpf:<scheme>/<setAlpha>."""
    for row in _read_rows(path):
        scheme, setA = row.get("scheme"), row.get("setAlphaMethod")
        n, h = _f(row.get("N")), _f(row.get("h"))
        if not scheme or not n or not h:
            continue
        rec = {"N": int(n), "h": h, "R_OVER_H": _f(row.get("R_over_h")), "DELTA_X": None,
               "N_INTERFACE": row.get("nInterface"), "N_WISP": "",
               "CPU_SECONDS": _f(row.get("k1")), "case": row.get("case", "")}
        rec.update({m: None for m in METRICS})
        rec["E_L1_N"] = _f(row.get("LNormalDiff1"))
        rec["E_L2_N"] = _f(row.get("LNormalDiff2"))
        rec["E_LINF_N"] = _f(row.get("LNormalDiffInf"))
        # reconstructionError's original centre columns: mean and max of the
        # distance of the PLIC centre from the surface (no rms).
        rec["E_POS_L1"] = _f(row.get("LCentre1"))
        rec["E_POS_LINF"] = _f(row.get("LCentreInf"))
        rec["E_AREA_REL"] = _f(row.get("E_area_rel"))
        rec["N_CURV_FALLBACK"] = ""
        for d in DIAG:
            rec[d] = ""
        series.setdefault((mesh, "tpf", f"tpf:{scheme}/{setA}"), []).append(rec)
    for recs in series.values():
        recs.sort(key=lambda r: -r["h"])
    return series


def _style(src, model):
    base = model.split(":")[-1].split("/")[0] if model.startswith("tpf:") else model
    color = TPF_COLORS.get(base, "#e34948") if model.startswith("tpf:") else COLORS.get(base, "#e34948")
    return color, MARKERS.get(base, "x"), STYLES.get(src, "-")


def _mesh_rank(mesh):
    return MESH_ORDER.index(mesh) if mesh in MESH_ORDER else len(MESH_ORDER)


def _sorted_series(series):
    """Series in table order: mesh block (hex, perturbed, poly), then (alphaSource, model)."""
    return sorted(series.items(), key=lambda kv: (_mesh_rank(kv[0][0]), kv[0][0], kv[0][1], kv[0][2]))


def long_rows_of(series):
    """One row per (series, N); the "mesh" key is dropped for a single-mesh study's CSV."""
    long_rows = []
    for (mesh, src, model), recs in _sorted_series(series):
        for r in recs:
            row = {"mesh": mesh, "alphaSource": src, "model": model, "N": r["N"], "h": r["h"],
                   "R_over_h": r["R_OVER_H"], "N_interface": r["N_INTERFACE"],
                   "N_wisp": r["N_WISP"]}
            row.update({m: r[m] for m in METRICS})
            row.update({d: r[d] for d in DIAG})
            row["N_curv_fallback"] = r.get("N_CURV_FALLBACK", "")
            row["cpu_s"] = r["CPU_SECONDS"]
            row["case"] = r["case"]
            long_rows.append(row)
    return long_rows


def order_rows_of(series):
    """Least-squares and pairwise orders per series (one row each, "mesh" key first)."""
    order_rows = []
    for (mesh, src, model), recs in _sorted_series(series):
        h = [r["h"] for r in recs]
        row = {"mesh": mesh, "alphaSource": src, "model": model, "N_coarsest": recs[0]["N"],
               "N_finest": recs[-1]["N"], "n_levels": len(recs)}
        for m in METRICS:
            e = [r[m] for r in recs]
            p, r2, n = _fit(h, e)
            row[f"{m}_finest"] = e[-1]
            row[f"p_{m}"] = p
            row[f"R2_{m}"] = r2
            row[f"pairwise_{m}"] = " ".join(
                "--" if q is None else f"{q:.2f}" for q in _pairwise(h, e))
        # The orders over the three finest levels (the proposal's success
        # criterion counts "the three finest of at least four refinements").
        for m in ("E_L1_N", "E_L2_N", "E_POS_L2", "E_KAPPA_L2"):
            e = [r[m] for r in recs]
            p3 = _fit(h[-3:], e[-3:])[0] if len(recs) >= 3 else None
            row[f"p3_{m}"] = p3
        order_rows.append(row)
    return order_rows


def _write_csv(path, rows, drop=()):
    fields = [k for k in rows[0] if k not in drop]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _num(x, fmt="{:.2f}"):
    return "--" if x is None else fmt.format(x)


def _p(v):   # an order: two decimals, no "-0.00"
    return "--" if v is None else f"{(0.0 if abs(v) < 0.005 else v):.2f}"


# The proposal table: the informative states and the three comparators
# (isoAlpha as a row since 2026-10-06: the fourth review asked for the
# comparator the text quotes; the plane-indicator state stays in the text).
PROPOSAL_MODELS = ("davof", "plicRDF", "gradAlpha", "isoAlpha")
PROPOSAL_SOURCES = ("exactSphere", "quadraticFaces", "detrixheAslam", "linearInterpolant")
_SRC_ORDER = ["exactSphere", "quadraticFaces", "detrixheAslam", "linearInterpolant",
              "planePhaseIndicator"]


def _same(a, b):
    # Equal to the digits the table prints (the CSVs carry 12 significant
    # digits; on a polyhedral mesh the comparators of the two states differ in
    # the last one, 2026-10-12), so 1e-9 relative.
    for k in ("E_L2_N_finest", "E_LINF_N_finest", "p_E_L2_N", "p_E_L1_N",
              "E_POS_L2_finest"):
        x, y = a.get(k), b.get(k)
        if (x is None) != (y is None):
            return False
        if x is not None and abs(x - y) > 1e-9*max(1.0, abs(x)):
            return False
    return True


def write_proposal_table(path, order_rows, with_mesh=False):
    """The compact table for the proposal: one row per (state, method).

    p(n): least-squares order of E_L2_N, p(x) of E_POS_L2 (E_POS_L1 where no
    rms exists, i.e. the TwoPhaseFlow rows), p(kappa) of E_KAPPA_L2 (the
    headline curvature model on the davof row, "--" for the comparators and
    for a study without curvature data); the errors at the finest mesh (the
    normal's rms and its maximum, the position's rms), the curvature one
    relative to the rms of the exact curvature. \\input by the DAVOF proposal
    (figures/davof/). The area order (2.00 for every state) is quoted in the
    proposal's text, not tabulated.

    with_mesh: a "mesh" column first (the label on the first row of each mesh
    block) and a \\midrule between the blocks (the combined three-mesh table).
    Without it the output of a single-mesh study is the table of before
    2026-10-12, byte for byte. Returns the list of collapsed blocks."""
    prop_rows = [r for r in order_rows
                 if r["alphaSource"] in PROPOSAL_SOURCES and r["model"] in PROPOSAL_MODELS]
    prop_rows.sort(key=lambda r: (_mesh_rank(r["mesh"]), r["mesh"],
                                  _SRC_ORDER.index(r["alphaSource"])
                                  if r["alphaSource"] in _SRC_ORDER else 99,
                                  PROPOSAL_MODELS.index(r["model"])))
    meshes = []
    for r in prop_rows:
        if r["mesh"] not in meshes:
            meshes.append(r["mesh"])
    collapsed_blocks = []
    with open(path, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_davof_normal_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{" + ("l" if with_mesh else "") + "llcccccccc}\n\\toprule\n")
        fh.write(("mesh & " if with_mesh else "")
                 + "$(\\alpha_k,\\alpha_f)$ & method & $p_1(\\mathbf{n})$ & $p_2(\\mathbf{n})$ & "
                 "$L_2(\\mathbf{n})$ & $L_\\infty(\\mathbf{n})$ & $p(\\mathbf{x})$ & $L_2(d)$ [m] & $p(\\kappa)$ & "
                 "$L_2(\\kappa)/L_2(\\kappa_\\Sigma)$ \\\\\n\\midrule\n")

        def _pp(r, m):   # the order over all levels (over the three finest)
            p3 = r.get(f"p3_{m}")
            return _p(r.get(f"p_{m}")) + ("" if p3 is None else f" ({_p(p3)})")

        pending = [None]   # the mesh label still to print in the block

        def _write(label, r):
            model = "DAVOF" if r["model"] == "davof" else r["model"]
            px = r.get("p_E_POS_L2") if r.get("p_E_POS_L2") is not None else r.get("p_E_POS_L1")
            pk = r.get("p_E_KAPPA_L2") if r["model"] == "davof" else None
            ek = r.get("E_KAPPA_REL_L2_finest") if r["model"] == "davof" else None
            lead = ""
            if with_mesh:
                lead = (pending[0] or "") + " & "
                pending[0] = None
            fh.write(f"{lead}{label} & {model} & "
                     f"{_pp(r, 'E_L1_N')} & {_pp(r, 'E_L2_N')} & "
                     f"{_num(r['E_L2_N_finest'], '{:.2e}')} & "
                     f"{_num(r['E_LINF_N_finest'], '{:.2e}')} & "
                     f"{_p(px)} & {_num(r.get('E_POS_L2_finest'), '{:.2e}')} & "
                     f"{_p(pk)} & {_num(ek, '{:.2e}')} \\\\\n")

        for i, mesh in enumerate(meshes):
            if i:
                fh.write("\\midrule\n")
            pending[0] = MESH_LABEL.get(mesh, mesh)
            block = [r for r in prop_rows if r["mesh"] == mesh]
            # The comparators depend on alpha_k only. When every state of the
            # block shares alpha_k (the ellipsoid ladder: quadratic faces keep
            # the Detrixhe-Aslam cell fractions), their rows are identical to
            # the digit and are printed once, after the DAVOF rows, under one
            # label.
            davof_rows = [r for r in block if r["model"] == "davof"]
            comp = {}
            for r in block:
                if r["model"] != "davof":
                    comp.setdefault(r["model"], []).append(r)
            collapsed = (len(davof_rows) > 1 and comp
                         and all(len(v) == len(davof_rows) and all(_same(v[0], x) for x in v[1:])
                                 for v in comp.values()))
            if collapsed:
                collapsed_blocks.append(mesh)
                for r in davof_rows:
                    _write(SOURCE_LABEL.get(r["alphaSource"], r["alphaSource"]), r)
                fh.write("\\midrule\n")
                first = True
                for m in PROPOSAL_MODELS:
                    if m in comp:
                        _write("both states ($\\alpha_k$ shared)" if first else "", comp[m][0])
                        first = False
            else:
                last_src = None
                for r in block:
                    if last_src is not None and r["alphaSource"] != last_src:
                        fh.write("\\midrule\n")
                    last_src = r["alphaSource"]
                    _write(SOURCE_LABEL.get(r["alphaSource"], r["alphaSource"]), r)
        fh.write("\\bottomrule\n\\end{tabular}\n")
    return collapsed_blocks


def mesh_class_siblings(study_dir):
    """(surface base name, [study dirs]) of the mesh-class ladders this study belongs to.

    A study <surface>Normal3D[_hex|_perturbed|_poly] has the siblings
    <surface>Normal3D (the original hex study), _hex, _perturbed and _poly next
    to it; the ones that exist are returned in mesh order. (None, []) for any
    other study."""
    study = os.path.basename(os.path.normpath(study_dir))
    parent = os.path.dirname(os.path.normpath(study_dir))
    base = study
    for suf in MESH_SUFFIXES:
        if study.endswith(suf):
            base = study[:-len(suf)]
            break
    if not re.fullmatch(r"[A-Za-z]+Normal3D", base):
        return None, []
    dirs = []
    for name in (base, base + "_hex", base + "_perturbed", base + "_poly"):
        d = os.path.join(parent, name)
        if os.path.isdir(d) and glob.glob(os.path.join(d, "*", CSV_MODELS)):
            dirs.append(d)
    return base, dirs


def write_combined(study_dir, tables):
    """The three-mesh table and orders CSV of the surface this study belongs to."""
    base, dirs = mesh_class_siblings(study_dir)
    if base is None:
        return
    expected = {"hex": [base, base + "_hex"], "perturbed": [base + "_perturbed"],
                "poly": [base + "_poly"]}
    present = {os.path.basename(d) for d in dirs}
    for mesh, names in expected.items():
        if not any(n in present for n in names):
            print(f"[davof] combined {base}: no completed {mesh} study ({' or '.join(names)}), skipped")
    series = {}
    for d in dirs:
        for key, recs in collect(d).items():
            series.setdefault(key, []).extend(recs)
    for recs in series.values():
        recs.sort(key=lambda r: -r["h"])
    if not series:
        return
    order_rows = order_rows_of(series)
    opath = os.path.join(tables, f"davof_normal_orders_{base}_meshes.csv")
    _write_csv(opath, order_rows)
    ppath = os.path.join(tables, f"davof_normal_proposal_{base}_meshes.tex")
    collapsed = write_proposal_table(ppath, order_rows, with_mesh=True)
    meshes = sorted({k[0] for k in series}, key=_mesh_rank)
    print(f"[davof] wrote {opath} and {ppath} (mesh classes: {', '.join(meshes)}; "
          f"from {', '.join(os.path.basename(d) for d in dirs)}"
          + (f"; comparator rows collapsed on {', '.join(collapsed)}" if collapsed else "") + ")")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="davof")
    ap.add_argument("--tpf-results", default=None,
                    help="results_sphereNormal3D.csv of the TwoPhaseFlow benchmark")
    ap.add_argument("--no-combined", action="store_true",
                    help="skip the three-mesh table of the sibling studies")
    a = ap.parse_args(argv)

    series = collect(a.study_dir)
    if a.tpf_results and os.path.isfile(a.tpf_results):
        merge_tpf(series, a.tpf_results)
    if not series:
        print(f"[davof] no completed cases under {a.study_dir}")
        return 1

    study = os.path.basename(os.path.normpath(a.study_dir))
    tables, figs = paths.tables_dir(a.theme), paths.figs_dir(a.theme)
    meshes = sorted({k[0] for k in series}, key=_mesh_rank)
    if len(meshes) > 1:
        print(f"[davof] WARNING: {study} mixes mesh classes {meshes}; its own tables omit the mesh column")

    # --- long table (the mesh column only in the combined outputs) ---------
    long_rows = long_rows_of(series)
    lpath = os.path.join(tables, f"davof_normal_{study}.csv")
    _write_csv(lpath, long_rows, drop=("mesh",))
    print(f"[davof] wrote {lpath} ({len(long_rows)} rows; mesh class {', '.join(meshes)})")

    # --- orders ---------------------------------------------------------------
    order_rows = order_rows_of(series)
    opath = os.path.join(tables, f"davof_normal_orders_{study}.csv")
    _write_csv(opath, order_rows, drop=("mesh",))
    print(f"[davof] wrote {opath} ({len(order_rows)} rows)")

    tpath = os.path.join(tables, f"davof_normal_orders_{study}.tex")
    with open(tpath, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_davof_normal_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{llrcccccccc}\n\\toprule\n")
        fh.write("$\\alpha$ source & model & $N$ & $L_2(\\mathbf{n})$ & $p(L_2)$ & $R^2$ & "
                 "$p(L_1)$ & $p(L_\\infty)$ & $p(A)$ & $L_2(\\mathbf{x})$ [m] & $p(\\mathbf{x})$ \\\\\n\\midrule\n")
        for r in order_rows:
            fh.write(f"{SOURCE_LABEL.get(r['alphaSource'], r['alphaSource'])} & "
                     f"{r['model'].replace('_', chr(92) + '_')} & {r['N_finest']} & "
                     f"{_num(r['E_L2_N_finest'], '{:.3e}')} & {_num(r['p_E_L2_N'])} & "
                     f"{_num(r['R2_E_L2_N'], '{:.3f}')} & {_num(r['p_E_L1_N'])} & "
                     f"{_num(r['p_E_LINF_N'])} & {_num(r['p_E_AREA_REL'])} & "
                     f"{_num(r.get('E_POS_L2_finest'), '{:.3e}')} & "
                     f"{_num(r.get('p_E_POS_L2') if r.get('p_E_POS_L2') is not None else r.get('p_E_POS_L1'))} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[davof] wrote {tpath}")

    # --- the compact table for the proposal: one row per (state, method) -------
    ppath = os.path.join(tables, f"davof_normal_proposal_{study}.tex")
    collapsed = write_proposal_table(ppath, order_rows, with_mesh=False)
    print(f"[davof] wrote {ppath}" + (" (comparator rows collapsed)" if collapsed else ""))

    # --- the three-mesh table of the surface (sibling studies) -------------------
    if not a.no_combined:
        write_combined(a.study_dir, tables)

    # --- figures: one log-log panel per metric --------------------------------
    hs = sorted({r["h"] for recs in series.values() for r in recs})
    for metric in FIG_METRICS:
        # No figure for a metric no series carries (a study without curvature
        # models has no curvature figures).
        if not any(any(r.get(metric) for r in recs) for recs in series.values()):
            continue
        fig, ax = plt.subplots(figsize=(8.4, 5.0))
        ymin, ymax = np.inf, -np.inf
        for (mesh, src, model), recs in _sorted_series(series):
            h = [r["h"] for r in recs]
            e = [r[metric] for r in recs]
            pts = [(a_, b) for a_, b in zip(h, e) if a_ and b and b > 0]
            if not pts:
                continue
            p, _, _ = _fit(h, e)
            color, marker, ls = _style(src, model)
            lab = f"{model}, {SOURCE_LABEL.get(src, src)}"
            if len(meshes) > 1:
                lab += f", {mesh}"
            if p is not None:
                lab += rf"  ($\propto h^{{{p:.2f}}}$)"
            ax.loglog([q[0] for q in pts], [q[1] for q in pts], ls, color=color,
                      marker=marker, ms=5, lw=1.6, label=lab)
            ymin = min(ymin, min(q[1] for q in pts))
            ymax = max(ymax, max(q[1] for q in pts))
        if np.isfinite(ymin) and len(hs) >= 2:
            h0, h1 = hs[0], hs[-1]
            y1 = 0.6*ymax
            ax.loglog([h0, h1], [y1*h0/h1, y1], color="#7f7f7f", ls="-", lw=0.8, alpha=0.7)
            ax.annotate(r"$h^1$", xy=(h0*1.05, 1.15*y1*h0/h1), color="#7f7f7f", fontsize=9)
            y2 = 0.3*ymax
            ax.loglog([h0, h1], [y2*(h0/h1)**2, y2], color="#7f7f7f", ls="-", lw=0.8, alpha=0.7)
            ax.annotate(r"$h^2$", xy=(h0*1.05, 1.15*y2*(h0/h1)**2), color="#7f7f7f", fontsize=9)
        ax.set_xlabel(r"$h$ [m]")
        ax.set_ylabel(METRIC_LABEL.get(metric, metric))
        ax.set_title(f"DAVOF static gate {study}: " + metric.replace("_", " "),
                     fontsize=11)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fpath = os.path.join(figs, f"davof_normal_{metric.lower()}_{study}.png")
        fig.savefig(fpath, dpi=200); plt.close(fig)
        print(f"[davof] wrote {fpath}")

    for r in order_rows:
        print(f"[davof] {r['mesh']:10s} {r['alphaSource']:20s} {r['model']:28s} N={r['N_coarsest']}..{r['N_finest']}"
              f"  p(L1)={_num(r['p_E_L1_N'])}  p(L2)={_num(r['p_E_L2_N'])}"
              f"  p(Linf)={_num(r['p_E_LINF_N'])}  p(A)={_num(r['p_E_AREA_REL'])}"
              f"  L2@finest={_num(r['E_L2_N_finest'], '{:.3e}')}"
              f"  pos: p(L2)={_num(r.get('p_E_POS_L2'))} p(L1)={_num(r.get('p_E_POS_L1'))}"
              f" p(Linf)={_num(r.get('p_E_POS_LINF'))}"
              f" L2@finest={_num(r.get('E_POS_L2_finest'), '{:.3e}')}")
        if r.get("p_E_KAPPA_L2") is not None:
            print(f"[davof] {'':10s} {'':20s} {'':28s} curvature: p(L2)={_num(r.get('p_E_KAPPA_L2'))}"
                  f" p(L1)={_num(r.get('p_E_KAPPA_L1'))} p(Linf)={_num(r.get('p_E_KAPPA_LINF'))}"
                  f" p(K)={_num(r.get('p_E_K_L2'))} p(cell)={_num(r.get('p_E_KAPPA_CELL_L2'))}"
                  f" p(face)={_num(r.get('p_E_KAPPA_FACE_L2'))}"
                  f" p(foot)={_num(r.get('p_E_KAPPA_FACE_FOOT_L2'))}"
                  f" rel L2@finest={_num(r.get('E_KAPPA_REL_L2_finest'), '{:.3e}')}"
                  f" pairwise={r.get('pairwise_E_KAPPA_L2')}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
