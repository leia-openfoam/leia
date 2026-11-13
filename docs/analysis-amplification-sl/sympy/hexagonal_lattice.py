#!/usr/bin/env python3
"""The amplification analysis on a regular HEXAGONAL tiling (a honeycomb), in 2D.

WHY THIS EXISTS. The rest of the analysis was done on meshes of boxes. The wave argument
needs the same rule to repeat from cell to cell -- it does NOT need the mesh to be square --
so it should carry over to any repeating tiling. A honeycomb is the sharpest available test
of that, and it turns out to have a tighter stability limit than boxes do, for a reason that
generalises the boundary result.

Note on words: HEXAHEDRAL means a box (six faces, a solid); HEXAGONAL means a six-sided flat
tile. This script is about the latter.

WHAT THIS MEASURES.
  1. The structure of the stencil. In a regular hexagonal tiling every vertex is shared by
     exactly THREE cells, so "shares a vertex" and "shares an edge" select the SAME six
     cells. The stencil has six members, against eight for a Cartesian 2D mesh.
  2. Six data for five coefficients: one spare, against three on Cartesian 2D. And all six
     neighbours sit at the same distance, so every inverse-distance weight is equal and the
     weighting has no effect whatsoever on a regular honeycomb.
  3. The growth factor against Courant number and flow direction, and the limit at which it
     first exceeds one.
  4. The same for a cell against a wall.

THE TWO SETS OF OFFSETS, which a Cartesian mesh hides. The PHYSICAL positions build the
basis and the weights (the fit); the WHOLE-NUMBER lattice labels build the wave phases (the
symbol). On a box mesh these are the same list; on a honeycomb they are not.

VALIDITY GATE. workflow/scripts/sl_amplification_symbol.py validate-1d must pass first: it
checks the shared machinery symbolically against the known 1D answer.

Usage:
  hexagonal_lattice.py            all four sections
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "workflow", "scripts"))
import sl_amplification_symbol as M                                  # noqa: E402

A1 = np.array([1.0, 0.0])
A2 = np.array([0.5, math.sqrt(3.0) / 2.0])
INTS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
PHYS = [m1 * A1 + m2 * A2 for m1, m2 in INTS]


def basis(p, ncoeff):
    """The code's constant-free 2D quadratic basis, truncated to the fitted rank."""
    P = np.atleast_2d(np.asarray(p, dtype=float))
    full = np.stack([P[:, 0], P[:, 1],
                     0.5 * P[:, 0] ** 2, 0.5 * P[:, 1] ** 2,
                     P[:, 0] * P[:, 1]], axis=1)
    return full[:, :ncoeff]


def growth(phys, ints, c, faces=(), ngrid=320):
    """max|G| for a stencil, with the rank ladder of the code (quadratic, else linear)."""
    slots = np.array(list(phys) + [f for f in faces], dtype=float)
    ncell = len(phys)
    w2 = np.sum(slots ** 2, axis=1) ** -1.0
    for ncoeff in (5, 2):
        B = basis(slots, ncoeff)
        Mn = (B * w2[:, None]).T @ B
        if len(slots) >= ncoeff and np.linalg.cond(Mn) < 1e12:
            break
    z = np.linalg.solve(Mn, basis(c, ncoeff)[0])
    g = w2 * (B @ z)
    gcell, gface = g[:ncell], g[ncell:]
    centre = 1.0 - g.sum() + gface.sum()          # a zeroGradient face carries psi_c
    lam = abs(centre) + np.abs(gcell).sum()
    sv = np.array(list(ints), dtype=float)
    ax = [np.linspace(-np.pi, np.pi, ngrid, endpoint=False)] * 2
    G = np.meshgrid(*ax, indexing="ij")
    ph = sum(G[a][..., None] * sv[None, ..., a] for a in range(2))
    mag = np.abs(centre + (gcell * np.exp(1j * ph)).sum(axis=-1)).max()
    return lam, float(mag), ("quadratic" if ncoeff == 5 else "LINEAR")


def foot(Co, angle_deg):
    t = math.radians(angle_deg)
    return np.array([[-Co * math.cos(t), -Co * math.sin(t)]])


def limit(phys, ints, angle_deg, hi=1.4):
    lo = 0.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if growth(phys, ints, foot(mid, angle_deg))[1] > 1 + 1e-9:
            hi = mid
        else:
            lo = mid
    return lo


def main():
    print("=" * 78)
    print("1. THE STENCIL")
    print("=" * 78)
    print(f"  primitive vectors  a1 = {A1},  a2 = ({A2[0]}, sqrt3/2 = {A2[1]:.6f})")
    for (m1, m2), p in zip(INTS, PHYS):
        ang = math.degrees(math.atan2(p[1], p[0])) % 360
        print(f"    lattice ({m1:+d},{m2:+d})  physical ({p[0]:+.6f},{p[1]:+.6f})  "
              f"|s| = {np.linalg.norm(p):.6f}  angle = {ang:5.1f} deg")
    print("  every vertex of the tiling meets 3 cells, so vertex- and edge-neighbours")
    print("  are the SAME six cells (on a box mesh they are 8 and 4).")
    B = basis(np.array(PHYS), 5)
    Mn = B.T @ B
    print(f"\n  6 neighbours, 5 coefficients -> 1 spare (Cartesian 2D: 8 and 5, 3 spare)")
    print(f"  normal matrix: condition {np.linalg.cond(Mn):.4f}, determinant "
          f"{np.linalg.det(Mn):.6f} -> well posed")
    print("  all six at distance 1 -> every weight equal -> the weighting cannot matter")

    print("\n" + "=" * 78)
    print("2. GROWTH FACTOR, INTERIOR CELL")
    print("=" * 78)
    print(f"  {'Co':>7} {'direction':>26} {'Lambda':>10} {'max|G|':>12}")
    for Co in (0.0164, 0.1, 0.25, 0.5, 0.85, 0.9, 1.0):
        for ang, lab in ((0, "straight at a neighbour"), (30, "between two neighbours")):
            lam, mag, _ = growth(PHYS, INTS, foot(Co, ang))
            print(f"  {Co:7.4f} {lab:>26} {lam:10.6f} {mag:12.8f}")

    print("\n" + "=" * 78)
    print("3. THE LIMIT, AND ITS GEOMETRY")
    print("=" * 78)
    for ang, lab in ((0, "straight at a neighbour"), (30, "between two neighbours")):
        print(f"  {lab:>26}: Co limit = {limit(PHYS, INTS, ang):.9f}")
    print(f"  {'sqrt(3)/2':>26} = {math.sqrt(3)/2:.9f}"
          "   <- centre to the nearest side of the hexagon of neighbours")
    print("  Cartesian 2D for comparison: 1.000000 along an axis, 1.041283 on the diagonal.")

    print("\n" + "=" * 78)
    print("4. A CELL AGAINST A WALL (wall on the -x side)")
    print("=" * 78)
    keep = [i for i, p in enumerate(PHYS) if p[0] > 1e-9]
    kp = [PHYS[i] for i in keep]
    ki = [INTS[i] for i in keep]
    print(f"  {len(kp)} neighbours remain for 5 coefficients -> the quadratic CANNOT be")
    print("  formed; the fit drops to a straight line along the whole boundary.")
    print("  (A box mesh keeps 17 neighbours against 9 coefficients and holds its quadratic.)")
    face = (-math.sqrt(3) / 4, 0.0)
    for Co in (0.0164, 0.25):
        print(f"\n  Co = {Co}")
        for lab, ph_, it_, c, fs in (
            ("interior, all six neighbours", PHYS, INTS, foot(Co, 0), ()),
            ("at the wall, foot INTO the domain", kp, ki, -foot(Co, 0), ()),
            ("at the wall, foot TOWARDS the wall", kp, ki, foot(Co, 0), ()),
            ("the same, with one face value added", kp, ki, foot(Co, 0), (face,)),
        ):
            lam, mag, fit = growth(ph_, it_, c, fs)
            print(f"    {lab:38s} {fit:9s} max|G| = {mag:.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
