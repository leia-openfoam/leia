#!/usr/bin/env python3
r"""Measure the SDPLS source's own discretization order, in isolation.

WHAT THIS ISOLATES. The band gradient error converges at ~+1.2 (2D) and ~+0.67
(3D shear) although every discrete operator is second order. The error estimate
of Section "Where the order is lost" attributes that to

    delta_g = (grad psi_h)_h - grad psi = C2 h^2 grad^3 psi  +  (grad e)_h
                                          \____ O(h^2) ____/   \_ O(h^(p-1)) _/

with e = psi_h - psi the SOLUTION error. e is O(h^p) in magnitude but is not
smooth on the scale h, so differentiating it costs an order and the second term
dominates for p = 2.

The estimate makes a sharp, falsifiable prediction: at t = 0 the discrete field
IS the sampled exact signed distance, so e = 0 identically, the error-field term
vanishes, and the order must be exactly 2. This script tests that directly on
the strain a = n.grad(v).n itself -- not on |grad psi|, but on the quantity the
source actually applies -- by comparing the solver's own registered R field
against the analytic value.

It reads t = 0 of an existing study, so it costs no new simulation.

    a_exact = n . grad(v) . n,   n = (x - c)/|x - c|      (exact sphere normal)

with v the prescribed Rider-Kothe field the case ran with. Errors are reported
over the narrow band, where the source is what matters.

Usage (from repo root):
    python3 workflow/scripts/sdpls_source_error.py <study> [--case 2Dvortex]
"""
import argparse
import glob
import json
import math
import os
import re
import sys

import numpy as np


def read_internal(path):
    """internalField of an OpenFOAM scalar/vector field, as a numpy array."""
    with open(path) as fh:
        s = fh.read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*\n(\d+)\s*\n\(",
                  s)
    if m:
        kind, n = m.group(1), int(m.group(2))
        start = s.index("(", m.end() - 1) + 1
        end = s.index("\n)", start)
        body = s[start:end]
        if kind == "scalar":
            return np.fromstring(body, sep=" ")
        vals = np.fromstring(body.replace("(", " ").replace(")", " "), sep=" ")
        return vals.reshape(-1, 3)
    m = re.search(r"internalField\s+uniform\s+\(([^)]*)\)\s*;", s)
    if m:
        return np.array([float(x) for x in m.group(1).split()])
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", s)
    if m:
        return float(m.group(1))
    return None


def cell_centres_uniform(case, N, dims):
    """Cell centres of the uniform unit-cube/square block mesh, in OpenFOAM order.

    blockMesh numbers hexes x fastest, then y, then z, which is what this
    reproduces. Cheaper and more robust than parsing 8M points back out of
    polyMesh, and exact for the uniform meshes these studies use.
    """
    h = 1.0 / N
    xs = (np.arange(N) + 0.5) * h
    if dims == 2:
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        Z = np.full(X.size, 0.5 * h)          # one cell thick
        return np.column_stack([X.ravel(), Y.ravel(), Z])
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="xy")
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])


def grad_v_riderkothe_2d(p):
    r"""grad(v) for v = (sin(2 pi y) sin^2(pi x), -sin(2 pi x) sin^2(pi y), 0).

    du/dx =  pi sin(2 pi x) sin(2 pi y)      dv/dy = -du/dx   (divergence free)
    du/dy =  2 pi cos(2 pi y) sin^2(pi x)
    dv/dx = -2 pi cos(2 pi x) sin^2(pi y)
    """
    x, y = p[:, 0], p[:, 1]
    pi = math.pi
    dudx = pi * np.sin(2 * pi * x) * np.sin(2 * pi * y)
    dudy = 2 * pi * np.cos(2 * pi * y) * np.sin(pi * x) ** 2
    dvdx = -2 * pi * np.cos(2 * pi * x) * np.sin(pi * y) ** 2
    dvdy = -dudx
    G = np.zeros((p.shape[0], 3, 3))
    G[:, 0, 0] = dudx; G[:, 0, 1] = dudy
    G[:, 1, 0] = dvdx; G[:, 1, 1] = dvdy
    return G


def order(hs, es):
    pts = [(h, abs(e)) for h, e in zip(hs, es) if h and e and abs(e) > 0]
    if len(pts) < 2:
        return None
    h, e = zip(*pts)
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--case", default="2Dvortex")
    ap.add_argument("--centre", default="0.5,0.75,0")
    ap.add_argument("--radius", type=float, default=0.15)
    a = ap.parse_args(argv)

    c = np.array([float(x) for x in a.centre.split(",")])
    rows = []
    for cp in sorted(glob.glob(os.path.join("studies", a.study,
                                            f"{a.case}_[0-9]*", "case_params.json"))):
        cd = os.path.dirname(cp)
        t = json.load(open(cp))["tokens"]
        if t.get("SDPLS_SOURCE") != "R":
            continue
        # NO linearization filter: R = n.grad(v).n at t=0 is identical for every
        # linearization -- the split into Sc/Sp changes how the source enters the
        # matrix, not the strain field itself. Filtering would only shrink the
        # ladder for no reason.
        N = int(t["N_CELLS"])
        Rf = os.path.join(cd, "0", "R")
        nb = os.path.join(cd, "0", "NarrowBand")
        if not os.path.isfile(Rf):
            continue
        Rh = read_internal(Rf)
        if not isinstance(Rh, np.ndarray) or Rh.size < 2:
            print(f"  N={N}: 0/R is uniform (the t=0 strain-ordering bug) -- skipped")
            continue
        p = cell_centres_uniform(a.case, N, 2 if a.case.startswith("2") else 3)
        if p.shape[0] != Rh.size:
            print(f"  N={N}: {Rh.size} cells but {p.shape[0]} centres -- mesh mismatch, skipped")
            continue
        d = p - c
        if a.case.startswith("2"):
            # 2D: the mesh is one cell thick and its centres sit at z = h/2,
            # while the circle centre is on z = 0. Leaving that in gives the
            # "exact" normal a spurious z-component of order h/r -- 10% at
            # N=32 -- which decays like h and would be misread as an O(h)
            # source error. The interface normal is in-plane by construction.
            d[:, 2] = 0.0
        r = np.linalg.norm(d, axis=1)
        n = d / np.maximum(r, 1e-300)[:, None]
        G = grad_v_riderkothe_2d(p)
        a_ex = np.einsum("ki,kij,kj->k", n, G, n)
        err = np.abs(Rh - a_ex)
        if os.path.isfile(nb):
            band = read_internal(nb)
            mask = np.asarray(band) > 0.5 if isinstance(band, np.ndarray) else None
        else:
            mask = np.abs(r - a.radius) < 2.0 / N
        if mask is None or mask.sum() < 4:
            mask = np.abs(r - a.radius) < 2.0 / N
        e = err[mask]
        rows.append((N, 1.0 / N, e.mean(), math.sqrt((e ** 2).mean()), e.max(), int(mask.sum())))

    if not rows:
        print("no usable t=0 R fields found")
        return 1
    print(f"\nSDPLS source error at t=0:  || a_h - a_exact ||  over the narrow band")
    print(f"study={a.study}  case={a.case}  sphere c={tuple(c)} r={a.radius}\n")
    print(f"{'N':>6} {'h':>10} {'L1':>12} {'L2':>12} {'Linf':>12} {'band cells':>11}")
    for N, h, l1, l2, li, nb_ in rows:
        print(f"{N:>6} {h:>10.5f} {l1:>12.4e} {l2:>12.4e} {li:>12.4e} {nb_:>11d}")
    hs = [r[1] for r in rows]
    print()
    for j, name in ((2, "L1"), (3, "L2"), (4, "Linf")):
        o = order(hs, [r[j] for r in rows])
        print(f"  order({name:4s}) = {o:+.3f}" if o is not None else f"  order({name}) = --")
    return 0


if __name__ == "__main__":
    sys.exit(main())
