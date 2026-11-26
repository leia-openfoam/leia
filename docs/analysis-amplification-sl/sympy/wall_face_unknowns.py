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

RESULT OF THE FIRST RUN (2026-11-26, laptop, K = 6, 12^2 waves; MEASURED): P1 holds
(exactness error 2e-14). P2 FALSIFIED: the outflow wall amplifies, growth 3.6e-2 at Co 0.1,
4.0e-2 at Co 0.25, in the mode (ky, kz) = (pi, pi) of the wall cell and its face
(eigenvector |v| 1.0 in the wall cell, 0.6 to 0.8 in the face, 0.1 one layer in). For
comparison inflowOnly: 4.8e-3 at Co 0.1, in the second layer of the INFLOW wall.

SECOND CANDIDATE, pre-registered after the first run and before its own run (2026-11-26):
  faceUnknownsTensorI  the same unknowns; the wall cell and its face are updated by the
      TRIQUADRATIC INTERPOLANT on the 27-point tensor set x in {wall face, layer 0, layer 1}
      x {-1, 0, 1}^2 (the 9 face values of the 3x3 columns, the 9 centres of layer 0, the 9
      centres of layer 1). Every other cell keeps the least-squares quadratic. The
      departure point of the wall cell and of its face lies inside this set for outflow and
      tangential flow (Co <= 1/2 normal to the wall).
  P4  faceUnknownsTensorI: P1, and growth <= 1e-10 for every direction and Co <= 1.
  FALSIFIED if any direction amplifies by more than 1e-10.

RESULT faceUnknownsTensorI (2026-11-26, laptop, K = 6, 12^2 waves; MEASURED): exactness
3.6e-14. Tangential flow: growth <= 4.4e-16 at every Co (the outflow/tangential defect of the
first candidate is gone). P4 FALSIFIED by the INFLOW wall: growth 4.9e-3 at Co 0.1 (normal)
and 2.3e-3 (oblique), mode in layer 1 next to the inflow wall; the production inflowOnly has
the same growth there (4.8e-3). At Co 1 (normal) growth 1.04 in layer 0: the departure point
of the inflow wall cell lies outside the domain (beyond the face at -1/2).

THIRD CANDIDATE, pre-registered after that run and before its own run (2026-11-26):
  faceUnknownsTensorI2  faceUnknownsTensorI, and (a) a cell that is NOT a wall cell takes no
      boundary face data (its stencil is cells only, the interior C-Cp/L2 fit), and (b) a
      departure point outside the domain takes the boundary value of the inflow face that
      its characteristic crosses (exact data: zero perturbation).
  P5  faceUnknownsTensorI2: P1, and growth <= 1e-10 for every direction and Co <= 1.
  FALSIFIED if any direction amplifies by more than 1e-10.

RESULT faceUnknownsTensorI2 (2026-11-26, laptop; MEASURED): exactness 3.6e-14. K = 6,
12^2 waves: tangential flow <= 2.2e-16; every Co >= 0.1 and every direction <= 0 (decays);
the largest growth is 2.6e-4 at Co 0.0164 (normal) and 1.7e-4 (oblique). P5 FALSIFIED by
the pre-registered threshold. NOT pre-registered, measured after the verdict: the growing
mode spans the whole slab (eigenvector |v| 0.3 to 1.0 in every layer, k = pi in one
tangential direction), and its growth falls with the slab thickness: Co 0.0164: 2.6e-4 /
1.3e-4 / 5.1e-5 at K = 6 / 12 / 24; Co 0.05: 4.6e-4 / 8.2e-5 / -1.5e-11. The production
inflowOnly has 1.2e-3 (Co 0.0164) and 3.1e-3 (Co 0.05) at every K: a growth that sits at
the wall and does not fall with K.

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
        if variant == "faceUnknownsTensorI2" and layer != i: continue
        kind = wall_kind(d, side)
        if variant == "inflowOnly" and kind != "inflow": continue
        for dy, dz in itertools.product((-1, 0, 1), repeat=2):
            out.append(((xf - i, dy, dz), "face", side, dy, dz))
    return out

def face_value_rule(variant, kind):
    if variant == "zeroGradient": return "wallcell"
    if kind == "inflow": return "exact"
    if variant.startswith("faceUnknowns"): return "unknown"
    return "none"              # inflowOnly: the face is not data

def lagrange(nodes, x):
    w = []
    for j, xj in enumerate(nodes):
        v = 1.0
        for m, xm in enumerate(nodes):
            if m != j: v *= (x - xm)/(xj - xm)
        w.append(v)
    return w

def tensor_slots(i, K, x):
    """Triquadratic interpolation weights for the wall cell of layer i (0 or K-1) at the
    displacement x from its centre: list of (weight, kind, ref, dy, dz), kind cell/face/centre."""
    side = -1 if i == 0 else +1
    xn = (-0.5, 0.0, 1.0) if side == -1 else (0.5, 0.0, -1.0)
    wx = lagrange(xn, x[0]); wy = lagrange((-1, 0, 1), x[1]); wz = lagrange((-1, 0, 1), x[2])
    out = []
    for a, xa in enumerate(xn):
        for b, yb in enumerate((-1, 0, 1)):
            for c, zc in enumerate((-1, 0, 1)):
                w = wx[a]*wy[b]*wz[c]
                if a == 0: out.append((w, "face", side, yb, zc))
                elif a == 1: out.append((w, "cell", i, yb, zc))
                else: out.append((w, "cell", i - side, yb, zc))
    return out

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
    tens = variant.startswith("faceUnknownsTensorI")
    for i in range(K):
        if variant == "faceUnknownsTensorI2" and (i - d[0] < -0.5 or i - d[0] > K - 0.5):
            info[i] = "BC"; continue          # departure point outside: boundary value, exact
        if tens and i in (0, K - 1):
            for w, kind, ref, dy, dz in tensor_slots(i, K, -np.asarray(d)):
                add_slot_value(i, w, kind, ref, dy, dz)
            info[i] = "T"; continue
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
        x = np.array([xf - layer, 0, 0]) - np.asarray(d)
        if tens:
            r = fidx[side]
            for w, kind, ref, dy, dz in tensor_slots(layer, K, x):
                add_slot_value(r, w, kind, ref, dy, dz)
            info[("face", side)] = "T"; continue
        sl = slots_of_cell(layer, K, variant, d)
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
    tens = variant.startswith("faceUnknownsTensorI")
    def tval(w, kind, ref, dy, dz, layer):
        side = -1 if layer == 0 else 1
        if kind == "face":
            if face_value_rule(variant, wall_kind(d, ref)) == "exact" or True:
                return q(np.array([-0.5 if side == -1 else K - 0.5, dy, dz], float))
        return q(np.array([ref, dy, dz], float))
    for i in range(K):
        if variant == "faceUnknownsTensorI2" and (i - d[0] < -0.5 or i - d[0] > K - 0.5):
            continue
        if tens and i in (0, K - 1):
            c = np.array([i, 0., 0.])
            v = sum(w*tval(w, kind, ref, dy, dz, i) for w, kind, ref, dy, dz in tensor_slots(i, K, -d))
            err = max(err, abs(v - q(c - d)))
            xf = np.array([-0.5 if i == 0 else K - 0.5, 0., 0.])
            v = sum(w*tval(w, kind, ref, dy, dz, i) for w, kind, ref, dy, dz in tensor_slots(i, K, xf - c - d))
            err = max(err, abs(v - q(xf - d)))
            continue
        c = np.array([i, 0., 0.]); sl = slots_of_cell(i, K, variant, d)
        g, _ = fit_weights([s[0] for s in sl], -d)
        v = q(c) + sum(gs*(val(*s[:3], c) - q(c)) for gs, s in zip(g, sl))
        err = max(err, abs(v - q(c - d)))
    for side in (-1, +1):
        if tens or face_value_rule(variant, wall_kind(d, side)) != "unknown": continue
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
    ap.add_argument("--variants", nargs="+", default=["inflowOnly", "zeroGradient", "faceUnknowns", "faceUnknownsTensorI", "faceUnknownsTensorI2"])
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
    for variant in a.variants:
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
    for vn in a.variants:
        fu = [(k, v) for k, v in worst.items() if k[0] == vn]
        p1 = max(v[1] for k, v in fu); p2 = max(v[0] for k, v in fu)
        print(f"{vn}: largest exactness error {p1:.2e}, largest growth {p2:+.2e}"
              + (f", largest power norm {max(v[2] for k, v in fu):.3f}" if a.powers else ""))

if __name__ == "__main__":
    main()
