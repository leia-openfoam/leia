#!/usr/bin/env python3
"""Post-process the exported table bodies for the manuscript (never edit the bodies by hand).

For every table_*.tex written by the leia workflow this writes table_*_pp.tex with
  1. the state labels of the review of 2026-10-27: "DA + quadratic faces" -> "DA-Q",
     "Detrixhe-Aslam" -> "DA" (defined in the captions);
  2. the better DAVOF results in bold: in the static tables (plane, sphere, ellipsoid) a
     DAVOF cell is bold when it beats every comparator row (plicRDF, gradAlpha, isoAlpha) of
     the same mesh class, higher orders and lower errors; in the transport table of the
     plane-clipped loop a DAVOF row's cell is bold when it beats plicRDF and isoAlpha at the
     same Courant number. Tables without comparator rows are relabelled only.
Run before latexmk: python3 prepare_tables.py
"""
import glob
import re

LABELS = {"DA + quadratic faces": "DA-Q", "Detrixhe-Aslam": "DA"}
COMPARATORS = {"plicRDF", "gradAlpha", "isoAlpha"}
# column index -> True if higher is better (orders), False if lower is better (errors)
STATIC_COLS = {3: True, 4: True, 5: False, 6: False, 7: True, 8: False, 9: True, 10: False}
TRANSPORT_COLS = {4: False, 5: True, 6: False, 7: True, 8: False}


def num(cell):
    m = re.match(r"\s*(-?\d+(?:\.\d+)?(?:e[-+]?\d+)?)", cell.strip())
    return float(m.group(1)) if m else None


def split_row(line):
    body = line.rstrip()
    if not body.endswith(r"\\"):
        return None
    return [c.strip() for c in body[:-2].split("&")]


def join_row(cells):
    return " & ".join(cells) + r" \\"


def bold(cell):
    return r"\textbf{" + cell + "}"


def static_table(lines):
    rows = []  # (index, cells, meshclass)
    mesh = None
    for i, line in enumerate(lines):
        cells = split_row(line)
        if cells is None or len(cells) < 11 or cells[0].startswith("mesh"):
            continue
        if cells[0]:
            mesh = cells[0]
        rows.append((i, cells, mesh))
    for i, cells, mesh in rows:
        if cells[2] != "DAVOF":
            continue
        comps = [c for _, c, m in rows if m == mesh and c[2] in COMPARATORS]
        for col, higher in STATIC_COLS.items():
            v = num(cells[col])
            if v is None:
                continue
            cv = [num(c[col]) for c in comps]
            cv = [x for x in cv if x is not None]
            if not cv:
                continue
            if (higher and all(v > x for x in cv)) or (not higher and all(v < x for x in cv)):
                cells[col] = bold(cells[col])
        lines[i] = join_row(cells)


def transport_table(lines):
    rows = []
    name = None
    for i, line in enumerate(lines):
        cells = split_row(line)
        if cells is None or len(cells) < 9 or cells[0] == "normal":
            continue
        if cells[0]:
            name = cells[0]
        rows.append((i, cells, name))
    for i, cells, name in rows:
        if not name.startswith("DAVOF"):
            continue
        comps = [c for _, c, n in rows if n in COMPARATORS and c[1] == cells[1]]
        for col, higher in TRANSPORT_COLS.items():
            v = num(cells[col])
            if v is None:
                continue
            cv = [x for x in (num(c[col]) for c in comps) if x is not None]
            if not cv:
                continue
            if (higher and all(v > x for x in cv)) or (not higher and all(v < x for x in cv)):
                cells[col] = bold(cells[col])
        lines[i] = join_row(cells)


for src in sorted(glob.glob("table_*.tex")):
    if src.endswith("_pp.tex"):
        continue
    lines = open(src).read().splitlines()
    for k, v in LABELS.items():
        lines = [l.replace(k, v) for l in lines]
    if src in ("table_plane.tex", "table_sphere.tex", "table_ellipsoid.tex"):
        static_table(lines)
    elif src == "table_advection_sphere.tex":
        transport_table(lines)
    dst = src[:-4] + "_pp.tex"
    with open(dst, "w") as fh:
        fh.write("%% written by prepare_tables.py from %s -- do not edit\n" % src)
        fh.write("\n".join(lines) + "\n")
    print("wrote", dst)
