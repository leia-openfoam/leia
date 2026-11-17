#!/usr/bin/env python3
"""Cell-by-cell comparison of psi between a serial run and a reconstructed parallel run.

    ghost_locate_compare.py <serial case> <parallel case> [--out PREFIX] [--field psi]
                            [--rel 1e-10] [--band 3]

Both cases have the SAME mesh (the parallel case is reconstructed, so the cell order is the
order of the serial mesh) and psi written at the same steps, in ASCII. The parallel case
holds constant/cellDecomposition (decomposePar -cellDist). The serial case holds 0/C
(postProcess -func writeCellCentres) and constant/polyMesh/cellLevel.

PRE-REGISTERED READ-OUT (coordinator, 2026-11-14, before the runs; the full text is in
config/integration/ghostDecompLocateN64.yaml).
  A cell DIFFERS when |psi_serial - psi_parallel| > rel * max|psi_serial| (rel = 1e-10).
  H1 (wall tie) HOLDS if, after step 1, every differing cell is a level-1 cell that touches
     the transition and lies within 2 fine cells of a boundary patch. FALSIFIED if any
     differing cell lies more than 2 fine cells from every boundary patch.
  H2 (processor boundary) HOLDS if the differing cells after step 1 lie within 2 cells of a
     processor boundary and away from the boundary patches.
  SECOND QUESTION: the first step with a differing cell in the interface band
     (|psi| < band * h, h the size of the cell), and where these cells are.
  Counts and L1 only. The largest difference is a locator, not a verdict.

Definitions.
  level        cellLevel of refineHexMesh (0 coarse, 1 fine).
  transition   layer 1 = the cell has a face neighbour of another level; layer k = k - 1
               face steps from such a cell.
  wall         the distance from the cell centre to the nearest face of the domain box, in
               FINE cell sizes (0.5 = the cell is at the boundary, for a fine cell). Every
               physical patch counts; the script prints the patches and their types.
               "Within 2 fine cells" means a distance below 2.
  processor    layer 1 = the cell has a face neighbour on another rank; layer k as above.
  band         |psi_serial| < band * h_cell.

Output: a table with one row per written step (stdout and <PREFIX>_steps.csv), histograms at
selected steps, and the read-out.
"""
import argparse
import os
import re
import sys

import numpy as np


def _text(path):
    with open(path) as fh:
        return fh.read()


def read_scalar_field(path):
    t = _text(path)
    i = t.index("internalField")
    m = re.match(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", t[i:])
    if not m:
        raise SystemExit("%s: internalField is not a nonuniform scalar list" % path)
    n = int(m.group(1))
    a = i + m.end()
    b = t.index(")", a)
    v = np.array(t[a:b].split(), dtype=float)
    if v.size != n:
        raise SystemExit("%s: %d values, header says %d" % (path, v.size, n))
    return v


def read_vector_field(path):
    t = _text(path)
    i = t.index("internalField")
    m = re.match(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", t[i:])
    n = int(m.group(1))
    a = i + m.end()
    b = t.index("boundaryField", a)
    v = np.array(t[a:b].replace("(", " ").replace(")", " ").replace(";", " ").split(), dtype=float)
    return v[: 3 * n].reshape(n, 3)


def read_label_list(path):
    t = _text(path)
    t = t[t.index("}", t.index("FoamFile")) + 1:]
    t = re.sub(r"//.*", "", t)
    m = re.search(r"(\d+)\s*\{\s*(-?\d+)\s*\}", t)
    m2 = re.search(r"(\d+)\s*\(", t)
    if m and (not m2 or m.start() < m2.start()):
        return np.full(int(m.group(1)), int(m.group(2)), dtype=np.int64)
    n = int(m2.group(1))
    a = m2.end()
    b = t.index(")", a)
    v = np.array(t[a:b].split(), dtype=np.int64)
    if v.size != n:
        raise SystemExit("%s: %d labels, header says %d" % (path, v.size, n))
    return v


def read_boundary(path):
    t = re.sub(r"//.*", "", _text(path))
    out = []
    for m in re.finditer(r"(\w+)\s*\{([^{}]*)\}", t):
        body = m.group(2)
        ty = re.search(r"\btype\s+(\w+);", body)
        nf = re.search(r"\bnFaces\s+(\d+);", body)
        sf = re.search(r"\bstartFace\s+(\d+);", body)
        if ty and nf and sf:
            out.append((m.group(1), ty.group(1), int(nf.group(1)), int(sf.group(1))))
    return out


def layers(seed, own, nei, n_cells, max_layer=12):
    """Layer 1 = the seed cells; layer k = k - 1 face steps away; 99 = further."""
    d = np.full(n_cells, 99, dtype=np.int64)
    d[seed] = 1
    for k in range(1, max_layer):
        f = d == k
        if not f.any():
            break
        new = np.zeros(n_cells, dtype=bool)
        new[nei[f[own]]] = True
        new[own[f[nei]]] = True
        new &= d == 99
        d[new] = k + 1
    return d


def times(case, field):
    out = []
    for name in os.listdir(case):
        try:
            v = float(name)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(case, name, field)):
            out.append((v, name))
    return sorted(out)


def hist(values, bins, labels):
    return ", ".join("%s: %d" % (lab, int(((values >= lo) & (values < hi)).sum()))
                     for (lo, hi), lab in zip(bins, labels))


LAYER_BINS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 9), (9, 1000)]
LAYER_LAB = ["1", "2", "3", "4", "5-8", ">8"]
WALL_BINS = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 8), (8, 1e9)]
WALL_LAB = ["<1", "1-2", "2-3", "3-4", "4-8", ">=8"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("serial")
    ap.add_argument("parallel")
    ap.add_argument("--field", default="psi")
    ap.add_argument("--rel", type=float, default=1e-10)
    ap.add_argument("--band", type=float, default=3.0)
    ap.add_argument("--out", default="locate")
    a = ap.parse_args()

    pm = os.path.join(a.serial, "constant", "polyMesh")
    own_all = read_label_list(os.path.join(pm, "owner"))
    nei = read_label_list(os.path.join(pm, "neighbour"))
    own = own_all[: nei.size]
    n_cells = int(max(own_all.max(), nei.max())) + 1
    level = read_label_list(os.path.join(pm, "cellLevel"))
    proc = read_label_list(os.path.join(a.parallel, "constant", "cellDecomposition"))
    C = read_vector_field(os.path.join(a.serial, "0", "C"))
    patches = read_boundary(os.path.join(pm, "boundary"))
    if not (level.size == proc.size == C.shape[0] == n_cells):
        raise SystemExit("sizes differ: cells %d, level %d, proc %d, C %d" % (n_cells, level.size, proc.size, C.shape[0]))

    # cell size from the level: the largest cells have level 0
    wall_seed = np.zeros(n_cells, dtype=bool)
    print("## Mesh")
    print("cells: %d; level counts: %s; ranks: %s" % (
        n_cells, dict(zip(*np.unique(level, return_counts=True))), dict(zip(*np.unique(proc, return_counts=True)))))
    print()
    print("| patch | type | faces | x range of its cells | y range | z range | fine cells at the patch | coarse cells |")
    print("|---|---|---|---|---|---|---|---|")
    for name, ty, nf, sf in patches:
        c = own_all[sf: sf + nf]
        if ty == "processor" or nf == 0:
            continue
        wall_seed[c] = True
        cc = C[c]
        print("| %s | %s | %d | %.4g to %.4g | %.4g to %.4g | %.4g to %.4g | %d | %d |" % (
            name, ty, nf, cc[:, 0].min(), cc[:, 0].max(), cc[:, 1].min(), cc[:, 1].max(),
            cc[:, 2].min(), cc[:, 2].max(), int((level[c] > 0).sum()), int((level[c] == 0).sum())))
    # the domain box and the cell sizes: h of a cell from the layer-1 distance to its patch
    lo = C.min(axis=0)
    hi = C.max(axis=0)
    # coarse h: twice the smallest distance of a coarse cell centre to the box of centres is
    # zero, so use the x spacing of sorted unique coarse x centres
    xs = np.unique(np.round(C[level == level.min(), 0], 12))
    h0 = float(np.min(np.diff(xs))) if xs.size > 1 else 1.0
    h = h0 / (2.0 ** (level - level.min()))
    hf = float(h.min())
    box_lo = lo - 0.5 * h[np.argmin(C, axis=0)]
    box_hi = hi + 0.5 * h[np.argmax(C, axis=0)]
    d_wall = np.minimum(C - box_lo, box_hi - C).min(axis=1) / hf
    which = np.argmin(np.concatenate([C - box_lo, box_hi - C], axis=1), axis=1)
    side = np.array(["xmin", "ymin", "zmin", "xmax", "ymax", "zmax"])[which]
    print()
    print("domain box: %s to %s; coarse h = %.6g, fine h = %.6g" % (np.round(box_lo, 9), np.round(box_hi, 9), h0, hf))
    fine = level > level.min()
    if fine.any():
        print("fine cells: x from %.6g to %.6g, y from %.6g to %.6g, z from %.6g to %.6g (cell centres)" % (
            C[fine, 0].min(), C[fine, 0].max(), C[fine, 1].min(), C[fine, 1].max(), C[fine, 2].min(), C[fine, 2].max()))

    tr_seed = np.zeros(n_cells, dtype=bool)
    f = level[own] != level[nei]
    tr_seed[own[f]] = True
    tr_seed[nei[f]] = True
    pr_seed = np.zeros(n_cells, dtype=bool)
    f = proc[own] != proc[nei]
    pr_seed[own[f]] = True
    pr_seed[nei[f]] = True
    L_tr = layers(tr_seed, own, nei, n_cells)
    L_pr = layers(pr_seed, own, nei, n_cells)
    L_wall = layers(wall_seed, own, nei, n_cells)
    print("cells that touch the transition: %d (fine %d, coarse %d); cells at a processor boundary: %d; cells at a boundary patch: %d" % (
        int(tr_seed.sum()), int((tr_seed & fine).sum()), int((tr_seed & ~fine).sum()), int(pr_seed.sum()), int(wall_seed.sum())))
    print("fine cells that touch the transition and lie within 2 fine cells of a boundary patch: %d" % int((tr_seed & fine & (d_wall < 2)).sum()))

    ts, tp = times(a.serial, a.field), times(a.parallel, a.field)
    n = min(len(ts), len(tp))
    print()
    print("written times: serial %d, parallel %d; compared %d" % (len(ts), len(tp), n))
    mism = [(s[1], p[1]) for s, p in zip(ts[:n], tp[:n]) if s[1] != p[1]]
    if mism:
        print("time names differ in %d steps, first: %s against %s (paired by order)" % (len(mism), mism[0][0], mism[0][1]))

    def describe(mask, title, psi_s, d):
        k = int(mask.sum())
        print()
        print("### %s: %d cells" % (title, k))
        if k == 0:
            return
        print("- level: " + ", ".join("%d: %d" % (lv, int((level[mask] == lv).sum())) for lv in np.unique(level)))
        print("- transition layer: " + hist(L_tr[mask], LAYER_BINS, LAYER_LAB))
        print("- distance to the domain boundary in fine cells: " + hist(d_wall[mask], WALL_BINS, WALL_LAB))
        print("- boundary-patch layer (in cells): " + hist(L_wall[mask], LAYER_BINS, LAYER_LAB))
        print("- nearest side of the box: " + ", ".join("%s: %d" % (s, int((side[mask] == s).sum())) for s in np.unique(side[mask])))
        print("- processor-boundary layer: " + hist(L_pr[mask], LAYER_BINS, LAYER_LAB))
        print("- rank: " + ", ".join("%d: %d" % (r, int((proc[mask] == r).sum())) for r in np.unique(proc[mask])))
        band = np.abs(psi_s) < a.band * h
        print("- inside the interface band: yes %d, no %d" % (int((mask & band).sum()), int((mask & ~band).sum())))
        print("- x range %.6g to %.6g; y range %.6g to %.6g; z range %.6g to %.6g" % (
            C[mask, 0].min(), C[mask, 0].max(), C[mask, 1].min(), C[mask, 1].max(), C[mask, 2].min(), C[mask, 2].max()))
        print("- combinations (level, touches the transition, within 2 fine cells of the boundary, within 2 cells of a processor boundary): count")
        key = np.stack([level[mask], (L_tr[mask] == 1).astype(int), (d_wall[mask] < 2).astype(int), (L_pr[mask] <= 2).astype(int)], axis=1)
        u, cnt = np.unique(key, axis=0, return_counts=True)
        for row, c in zip(u, cnt):
            print("  - (%d, %s, %s, %s): %d" % (row[0], "yes" if row[1] else "no", "yes" if row[2] else "no", "yes" if row[3] else "no", c))
        idx = np.where(mask)[0]
        top = idx[np.argsort(-np.abs(d[idx]))[:8]]
        print("- the 8 largest differences (locators): cell, x, y, z, level, transition layer, wall distance, processor layer, psi serial, difference")
        for c in top:
            print("  - %d, %.6g, %.6g, %.6g, %d, %d, %.2f, %d, %.6e, %.3e" % (
                c, C[c, 0], C[c, 1], C[c, 2], level[c], L_tr[c], d_wall[c], L_pr[c], psi_s[c], d[c]))

    rows = []
    first_band = None
    detail_steps = {1, 2, 3, 5, 10, 20, 30, 60, 100, n - 1}
    details = []
    h1 = h2 = None
    for i in range(n):
        ps = read_scalar_field(os.path.join(a.serial, ts[i][1], a.field))
        pp = read_scalar_field(os.path.join(a.parallel, tp[i][1], a.field))
        d = pp - ps
        thr = a.rel * np.abs(ps).max()
        m = np.abs(d) > thr
        band = np.abs(ps) < a.band * h
        rows.append([i, ts[i][1], int(m.sum()), int((d != 0).sum()), float(np.abs(d).max()), float(np.abs(d).mean()),
                     int((m & (level == level.min())).sum()), int((m & fine).sum()), int((m & (L_tr == 1)).sum()),
                     int((m & (d_wall < 2)).sum()), int((m & (d_wall >= 2)).sum()), int((m & (L_pr <= 2)).sum()),
                     int((m & band).sum()), float(np.abs(d[band]).max()) if band.any() else 0.0,
                     float(np.abs(d[band]).mean()) if band.any() else 0.0, int(band.sum())])
        if i == 1:
            k = int(m.sum())
            if k:
                h1 = bool(np.all(fine[m] & (L_tr[m] == 1) & (d_wall[m] < 2)))
                h1_fals = bool(np.any(d_wall[m] >= 2))
                h2 = bool(np.all((L_pr[m] <= 2) & (d_wall[m] >= 2)))
            else:
                h1 = h1_fals = h2 = None
        if first_band is None and (m & band).any():
            first_band = i
            details.append((i, "Differing cells INSIDE the interface band at the first step where they exist (step %d, time %s)" % (i, ts[i][1]), m & band, ps, d))
        if i in detail_steps:
            details.append((i, "Differing cells at step %d (time %s)" % (i, ts[i][1]), m, ps, d))

    head = ["step", "time", "differing", "any difference", "largest abs diff (locator)", "mean abs diff, all cells",
            "coarse", "fine", "touch transition", "wall < 2 fine cells", "wall >= 2", "processor layer <= 2",
            "differing in band", "largest abs diff in band (locator)", "mean abs diff in band", "cells in band"]
    print()
    print("## Table per written step (threshold %.0e * max|psi|)" % a.rel)
    print()
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for r in rows:
        print("| " + " | ".join(("%.3e" % x) if isinstance(x, float) else str(x) for x in r) + " |")
    with open(a.out + "_steps.csv", "w") as fh:
        fh.write(",".join(hh.replace(",", ";") for hh in head) + "\n")
        for r in rows:
            fh.write(",".join(("%.17g" % x) if isinstance(x, float) else str(x) for x in r) + "\n")

    print()
    print("## Histograms")
    for i, title, mask, ps, d in sorted(details, key=lambda x: x[0]):
        describe(mask, title, ps, d)

    print()
    print("## Read-out")
    if h1 is None:
        print("- Step 1: no cell differs above the threshold. H1 and H2 cannot be tested at step 1.")
        nz = next((r for r in rows if r[2] > 0), None)
        print("- First step with a differing cell: %s" % ("step %d, time %s, %d cells" % (nz[0], nz[1], nz[2]) if nz else "none"))
    else:
        print("- H1 (wall tie): %s. Falsifying cells (more than 2 fine cells from every boundary patch): %s." % (
            "HOLDS" if h1 else "does NOT hold", "present" if h1_fals else "none"))
        print("- H2 (processor boundary): %s." % ("HOLDS" if h2 else "does NOT hold"))
    if first_band is None:
        print("- Second question: no cell of the interface band differs above the threshold at any written step.")
    else:
        print("- Second question: the first step with a differing cell in the interface band is step %d, time %s (see the histogram above)." % (
            first_band, ts[first_band][1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
