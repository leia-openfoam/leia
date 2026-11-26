#!/usr/bin/env python3
"""Boundary-face values as advected unknowns: growth factor and exactness of a wall slab.

ROUTE (supervisor's decision 2026-11-26; the candidate of 2026-11-12, "face values as their
own unknowns", used on boundary faces only). Every boundary face centre carries a psi value.
On an outflow or tangential wall it is advected by the same semi-Lagrangian step, from its
own departure point, evaluated in the fit of its wall cell. On an inflow face the boundary
condition gives the value (exact data: zero perturbation). The wall cell fits its quadratic
from its 26-point stencil (cells that exist) plus the boundary faces of its stencil cells.

MODEL. A slab of K cell layers normal to x (h = 1), walls at x = -1/2 and x = K - 1/2,
periodic in y and z (Bloch waves exp(i (ky y + kz z))). Unknowns per wave: K cell values
and one face value per wall. A constant velocity with displacement d = Co * e per step; the
departure point of a cell is its centre - d, evaluated in the fit of the ARRIVAL cell (the
solver's rule); the departure point of a face is the face centre - d, evaluated in the fit
of its wall cell. The fit is the anchored weighted least-squares quadratic of the code
(9 coefficients, weights 1/|s|^2, the rank ladder quadratic -> linear).

VARIANTS.
  inflowOnly    the production setting: outflow and tangential faces are not data.
  zeroGradient  every face carries the value of its wall cell (stencilBoundaryFaces include).
  faceUnknowns  outflow and tangential faces are advected unknowns; inflow faces are exact.

PRE-REGISTERED READ-OUT (written before the first run, 2026-11-26).
  P1  faceUnknowns reproduces a quadratic field exactly (|error| <= 1e-12) in every cell and
      face update, for normal, tangential and oblique flow at every Co.
  P2  faceUnknowns: the largest |eigenvalue| over the Bloch waves is <= 1 + 1e-10 for every
      Co in (0, 1] and every direction, on both walls.
  P3  faceUnknowns: max over the waves of ||A^n||_2 <= 10 for n <= 2000 (power-bounded).
  FALSIFIED if P2 fails for any direction with Co <= 1. Then the fallback is the complete
  tensor point set at the walls (the expensive form of the same candidate).
  Expected for comparison (not a criterion): inflowOnly falls back to LINEAR in the wall
  cell of an outflow or tangential wall, so it fails P1 there.

Usage: python3 wall_face_unknowns.py [--layers K] [--waves N] [--powers]
"""
import argparse, itertools, sys
import numpy as np

def basis(s):
    s = np.atleast_2d(np.asarray(s, float))
    return np.stack([s[:, 0], s[:, 1], s[:, 2],
                     .5*s[:, 0]**2, .5*s[:, 1]**2, .5*s[:, 2]**2,
                     s[:, 0]*s[:, 1], s[:, 0]*s[:, 2], s[:, 1]*s[:, 2]], axis=1)

def fit_weights(slots, x):
    """Weights g_s with value(x) = psi_c + sum_s g_s (psi_s - psi_c); anchored fit at 0.
    Returns g and the rank used (9 quadratic, 3 linear)."""
    slots = np.asarray(slots, float)
    w = 1.0/np.sum(slots**2, axis=1)
    for nc in (9, 3):
        B = basis(slots)[:, :nc]
        M = (B*w[:, None]).T @ B
        if np.linalg.cond(M) < 1e12:
            z = np.linalg.solve(M, basis(x)[0, :nc])
            return w*(B @ z), nc
    raise RuntimeError("singular at every rank")

def wall_kind(d, side):
    """side -1: wall at x = -1/2 (outward normal -x); +1: wall at x = K - 1/2."""
    un = side*d[0]            # u.n_out with the displacement d = u dt
    if un > 1e-14: return "outflow"
    if un < -1e-14: return "inflow"
    return "tangential"

def slots_of_cell(i, K, variant, d):
    """Stencil of the cell in layer i: list of (offset, kind, layer or wall side, dy, dz)."""
    out = []
    for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
        if (dx, dy, dz) == (0, 0, 0): continue
        j = i + dx
        if 0 <= j < K:
            out.append(((dx, dy, dz), "cell", j, dy, dz))
    # boundary faces of the stencil cells: the cells of layer 0 (or K-1) in the 3x3 columns
    for side, layer, xf in ((-1, 0, -0.5), (+1, K - 1, K - 0.5)):
        if abs(layer - i) > 1: continue
        kind = wall_kind(d, side)
        if variant == "inflowOnly" and kind != "inflow": continue
        for dy, dz in itertools.product((-1, 0, 1), repeat=2):
            out.append(((xf - i, dy, dz), "face", side, dy, dz))
    return out

def face_value_rule(variant, kind):
    if variant == "zeroGradient": return "wallcell"
    if kind == "inflow": return "exact"
    if variant == "faceUnknowns": return "unknown"
    return "none"              # inflowOnly: the face is not data

def assemble(K, d, variant, ky, kz):
    """Step matrix A (cells 0..K-1, then face of side -1, then face of side +1)."""
    n = K + 2
    A = np.zeros((n, n), complex)
    ph = lambda dy, dz: np.exp(1j*(ky*dy + kz*dz))
    fidx = {-1: K, +1: K + 1}
    info = {}
    def add_slot_value(row, coef, kind, ref, dy, dz):
        """Add coef * (value of the slot) to row."""
        if kind == "cell":
            A[row, ref] += coef*ph(dy, dz)
        else:
            rule = face_value_rule(variant, wall_kind(d, ref))
            wall_layer = 0 if ref == -1 else K - 1
            if rule == "wallcell": A[row, wall_layer] += coef*ph(dy, dz)
            elif rule == "unknown": A[row, fidx[ref]] += coef*ph(dy, dz)
            # exact: zero perturbation
    for i in range(K):
        sl = slots_of_cell(i, K, variant, d)
        g, nc = fit_weights([s[0] for s in sl], -np.asarray(d))
        info[i] = nc
        A[i, i] += 1 - g.sum()
        for gs, (off, kind, ref, dy, dz) in zip(g, sl):
            add_slot_value(i, gs, kind, ref, dy, dz)
    for side in (-1, +1):
        rule = face_value_rule(variant, wall_kind(d, side))
        if rule != "unknown": continue          # not a state variable in this variant
        layer = 0 if side == -1 else K - 1
        xf = -0.5 if side == -1 else K - 0.5
        sl = slots_of_cell(layer, K, variant, d)
        x = np.array([xf - layer, 0, 0]) - np.asarray(d)
        g, nc = fit_weights([s[0] for s in sl], x)
        info[("face", side)] = nc
        r = fidx[side]
        A[r, layer] += 1 - g.sum()
        for gs, (off, kind, ref, dy, dz) in zip(g, sl):
            add_slot_value(r, gs, kind, ref, dy, dz)
    keep = [k for k in range(n) if k < K or face_value_rule(variant, wall_kind(d, -1 if k == K else 1)) == "unknown"]
    return A[np.ix_(keep, keep)], info

def exactness(K, d, variant):
    """Apply one step to a quadratic field on a 3x3 block of columns; largest error."""
    rng = np.random.default_rng(1)
    a = rng.normal(size=10)
    q = lambda p: (a[0] + a[1]*p[0] + a[2]*p[1] + a[3]*p[2] + a[4]*p[0]**2 + a[5]*p[1]**2
                   + a[6]*p[2]**2 + a[7]*p[0]*p[1] + a[8]*p[0]*p[2] + a[9]*p[1]*p[2])
    d = np.asarray(d, float); err = 0.0
    def val(off, kind, ref, centre):
        if kind == "cell": return q(centre + np.asarray(off, float))
        rule = face_value_rule(variant, wall_kind(d, ref))
        if rule == "wallcell":
            wl = 0 if ref == -1 else K - 1
            return q(np.array([wl, centre[1] + off[1], centre[2] + off[2]], float))
        return q(centre + np.asarray(off, float))      # unknown or exact: the exact value
    for i in range(K):
        c = np.array([i, 0., 0.]); sl = slots_of_cell(i, K, variant, d)
        g, _ = fit_weights([s[0] for s in sl], -d)
        v = q(c) + sum(gs*(val(*s[:3], c) - q(c)) for gs, s in zip(g, sl))
        err = max(err, abs(v - q(c - d)))
    for side in (-1, +1):
        if face_value_rule(variant, wall_kind(d, side)) != "unknown": continue
        layer = 0 if side == -1 else K - 1; c = np.array([layer, 0., 0.])
        xf = np.array([-0.5 if side == -1 else K - 0.5, 0., 0.])
        sl = slots_of_cell(layer, K, variant, d)
        g, _ = fit_weights([s[0] for s in sl], xf - c - d)
        v = q(c) + sum(gs*(val(*s[:3], c) - q(c)) for gs, s in zip(g, sl))
        err = max(err, abs(v - q(xf - d)))
    return err

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--layers", type=int, default=8)
    ap.add_argument("--waves", type=int, default=24)
    ap.add_argument("--powers", action="store_true")
    a = ap.parse_args()
    K = a.layers
    dirs = {"normal +x": (1, 0, 0), "normal -x": (-1, 0, 0), "tangential y": (0, 1, 0),
            "oblique (1 1 0)": (1, 1, 0), "oblique (-1 1 0)": (-1, 1, 0),
            "oblique (1 1 1)": (1, 1, 1), "oblique (-1 2 1)": (-1, 2, 1),
            "grazing (0.05 1 0)": (0.05, 1, 0)}
    Cos = (0.0164, 0.1, 0.25, 0.3, 0.5, 0.75, 1.0)
    ks = np.linspace(-np.pi, np.pi, a.waves, endpoint=False)
    print(f"slab K = {K} layers, {a.waves}^2 tangential Bloch waves; growth = max|eig| - 1")
    worst = {}
    for variant in ("inflowOnly", "zeroGradient", "faceUnknowns"):
        print(f"\n=== {variant}")
        for name, e in dirs.items():
            e = np.asarray(e, float); e /= np.linalg.norm(e)
            row = []
            for Co in Cos:
                d = Co*e
                rho = 0.0; pmax = 0.0
                for ky in ks:
                    for kz in ks:
                        A, info = assemble(K, d, variant, ky, kz)
                        r = np.abs(np.linalg.eigvals(A)).max()
                        rho = max(rho, r)
                        if a.powers and r > 1 - 1e-9:
                            P = np.eye(len(A)); best = 1.0
                            for n in range(1, 2001):
                                P = A @ P
                                if n % 10 == 0 or n < 50: best = max(best, np.linalg.norm(P, 2))
                            pmax = max(pmax, best)
                ex = exactness(K, d, variant)
                ranks = sorted(set(str(v) for v in info.values()))
                row.append((Co, rho - 1, ex, ranks, pmax))
                worst[(variant, name, Co)] = (rho - 1, ex, pmax)
            print(f"  {name:20s} " + "  ".join(f"Co {Co:g}: {g:+.2e} ex {ex:.0e} r{'/'.join(rk)}"
                                               + (f" P {pm:.2f}" if a.powers else "")
                                               for Co, g, ex, rk, pm in row))
    fu = [(k, v) for k, v in worst.items() if k[0] == "faceUnknowns"]
    p1 = max(v[1] for k, v in fu); p2 = max(v[0] for k, v in fu)
    print(f"\nfaceUnknowns: largest exactness error {p1:.2e} (P1: <= 1e-12), "
          f"largest growth {p2:+.2e} (P2: <= 1e-10)")
    if a.powers:
        print(f"faceUnknowns: largest power norm {max(v[2] for k, v in fu):.3f} (P3: <= 10)")

if __name__ == "__main__":
    main()
