#!/usr/bin/env python3
"""Per-cell view of the one-step sagitta diagnostic (hexSagittaDiag, 2026-10-16) on the
blockMesh hexahedral mesh of sphereTranslation3D (2 N x N x N cells of size h = 1/N in the
2 x 1 x 1 box, i fastest): the cell's volume fraction against the exact sphere's at t (the
chord state's deficit), and the plane row's face errors grouped by the face's orientation
(x faces carry the flux, y and z faces have none) and by the sign of the cell's deficit.

    python3 davof_sagitta_cell_check.py <plane case> <sagitta case> <time> --N 32 --t 0.00625
"""
import argparse
import os
import re
import sys

import numpy as np


def read_field(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
    if not m:
        sys.exit(f"{path}: no nonuniform scalar internalField")
    n = int(m.group(1))
    body = txt[m.end():]
    end = body.index("\n)\n") if "\n)\n" in body else body.index(")\n;")
    vals = np.array([float(v) for v in body[:end].split()])
    assert len(vals) == n, (path, len(vals), n)
    return vals


def read_labels(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"\n(\d+)\s*\n\(", txt)
    n = int(m.group(1))
    body = txt[m.end():]
    end = body.index("\n)\n") if "\n)\n" in body else body.index(")\n;")
    vals = np.array([int(v) for v in body[:end].split()])
    assert len(vals) == n, (path, len(vals), n)
    return vals


def sphere_fraction(c, h, R, ncol=160):
    """alpha of the ball |x| < R in the cube centred at c (origin at the sphere centre)."""
    x = c[0] - h/2 + (np.arange(ncol) + 0.5)*h/ncol
    y = c[1] - h/2 + (np.arange(ncol) + 0.5)*h/ncol
    X, Y = np.meshgrid(x, y, indexing="ij")
    r2 = R*R - X*X - Y*Y
    zext = np.sqrt(np.maximum(r2, 0.0))
    zlo, zhi = c[2] - h/2, c[2] + h/2
    length = np.clip(np.minimum(zext, zhi) - np.maximum(-zext, zlo), 0.0, None)
    length[r2 <= 0] = 0.0
    return length.sum()*(h/ncol)**2/h**3


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("plane_case"); ap.add_argument("sagitta_case"); ap.add_argument("time")
    ap.add_argument("--N", type=int, default=32); ap.add_argument("--t", type=float, default=0.00625)
    ap.add_argument("--R", type=float, default=0.15); ap.add_argument("--group", default="water")
    a = ap.parse_args(argv)
    N = a.N; h = 1.0/N
    alpha = read_field(os.path.join(a.plane_case, a.time, f"alpha.{a.group}"))
    ep = read_field(os.path.join(a.plane_case, a.time, f"alphafError.{a.group}"))
    es = read_field(os.path.join(a.sagitta_case, a.time, f"alphafError.{a.group}"))
    dA = read_field(os.path.join(a.sagitta_case, a.time, f"sagittaCorrection.{a.group}"))
    own = read_labels(os.path.join(a.plane_case, "constant/polyMesh/owner"))
    nei = read_labels(os.path.join(a.plane_case, "constant/polyMesh/neighbour"))
    nInt = len(nei)
    ncell = len(alpha)
    assert ncell == 2*N*N*N, ncell
    idx = np.arange(ncell)
    ix = idx % (2*N); iy = (idx // (2*N)) % N; iz = idx // (2*N*N)
    C = np.stack([(ix + 0.5)*h, (iy + 0.5)*h, (iz + 0.5)*h], axis=1)
    centre = np.array([0.5 + a.t, 0.5, 0.5])
    # the cells whose faces carry an error entry: interface cells of either state
    face_cells = set(own[:nInt][ep != 0].tolist()) | set(nei[ep != 0].tolist())
    cells = sorted(face_cells)
    print(f"cells adjacent to a counted face: {len(cells)}")
    dalpha = {}
    for c in cells:
        aex = sphere_fraction(C[c] - centre, h, a.R)
        dalpha[c] = alpha[c] - aex
    d = np.array([dalpha[c] for c in cells])
    inner = np.array([0 < alpha[c] < 1 for c in cells])
    print(f"alpha - alpha_exact over these cells: mean {d.mean():+.3e}, mean|.| {np.abs(d).mean():.3e}, "
          f"min {d.min():+.3e}, max {d.max():+.3e}; cells with 0 < alpha < 1: {inner.sum()}")
    print(f"  over the cut cells only: mean {d[inner].mean():+.3e}, negative in {(d[inner] < 0).sum()} of {inner.sum()}")
    # the face errors by orientation: faces between cells i and i+1 (x), i and i+2N (y), else z
    diff = nei - own[:nInt]
    orient = np.where(diff == 1, 0, np.where(diff == 2*N, 1, 2))
    for o, name in enumerate("xyz"):
        sel = (ep != 0) & (orient == o)
        selc = sel & (dA != 0)
        print(f"{name} faces: counted {sel.sum()}, corrected {selc.sum()}; e_plane mean {ep[sel].mean():+.3e} mean|.| {np.abs(ep[sel]).mean():.3e}; "
              f"e_sag mean|.| {np.abs(es[sel]).mean():.3e}; dA mean {dA[selc].mean():+.3e}; e_plane+dA mean|.| {np.abs(ep[selc] + dA[selc]).mean():.3e}")
    # the x faces: the upwind (owner, lower x) cell's plane sets the face; its deficit
    sel = (ep != 0) & (orient == 0)
    up = own[:nInt][sel]
    dd = np.array([dalpha.get(c, 0.0) for c in up])
    e = ep[sel]
    r = np.corrcoef(dd, e)[0, 1] if len(e) > 2 else float("nan")
    print(f"x faces: correlation of e_plane with the upwind cell's alpha deficit: {r:+.3f}; "
          f"faces with a negative-deficit upwind cell: e_plane mean {e[dd < 0].mean():+.3e} (n {int((dd < 0).sum())}), "
          f"positive: {e[dd > 0].mean() if (dd > 0).any() else float('nan'):+.3e} (n {int((dd > 0).sum())})")
    # the y and z faces: which cell set them? compare the error with the owner's and the neighbour's deficit
    for o, name in ((1, "y"), (2, "z")):
        sel = (ep != 0) & (orient == o)
        do = np.array([dalpha.get(c, 0.0) for c in own[:nInt][sel]])
        dn = np.array([dalpha.get(c, 0.0) for c in nei[sel]])
        e = ep[sel]
        print(f"{name} faces: corr(e, owner deficit) {np.corrcoef(do, e)[0,1]:+.3f}, corr(e, neighbour deficit) {np.corrcoef(dn, e)[0,1]:+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
