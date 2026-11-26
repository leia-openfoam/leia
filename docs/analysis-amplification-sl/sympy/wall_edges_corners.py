#!/usr/bin/env python3
"""The boundary treatment `faceUnknownsTensorI2` at EDGES and CORNERS of a cubic mesh.

DECISION (supervisor, 2026-11-26): boundary-face values as advected unknowns for the AMR line
on cubic cells (docs/plan-sl-amr-preprint.md section 6 item 8; the slab analysis is
wall_face_unknowns.py). This script extends the construction to cells that touch two walls
(edge) or three walls (corner).

CONSTRUCTION (the general form of faceUnknownsTensorI2).
  - Unknowns: the cell centres, and every boundary point of the lattice: a point whose
    coordinates are cell-centre coordinates or wall coordinates, with at least one wall
    coordinate. One wall coordinate: a face centre; two: a point on a boundary edge (at the
    cell-centre position along the edge); three: a corner.
  - A cell that touches a wall is updated by the tensor-product quadratic interpolant on the
    27 points with the nodes {wall, 0, 1} (cell coordinates relative to the cell; mirrored
    at the high wall) in every bounded direction where the cell touches a wall, and
    {-1, 0, 1} in the other directions. A boundary point is updated by the interpolant of
    its host cell (the boundary cell next to it), at its own departure point.
  - Every other cell: the interior least-squares quadratic on its 26 neighbour cells, no
    boundary data.
  - A departure point outside the domain takes the boundary value where its characteristic
    enters; the data are exact, so the perturbation is zero (a zero row).
  The step matrix is built for the perturbation; Bloch waves exp(i k z) in a periodic
  direction.

DOMAINS.  edge:   K x K cells in x and y (four walls), periodic in z (Bloch kz).
          corner: K x K x K cells, six walls.

PRE-REGISTERED READ-OUT (2026-11-26, before the first run of this script).
  P6a every update reproduces a quadratic field exactly (|error| <= 1e-12; for the edge domain
      a quadratic without z, which is periodic).
  P6b growth (largest |eigenvalue| - 1) <= 1e-10 for every direction and 0.1 <= Co <= 1.
  P6c for Co < 0.1: growth not above the slab remainder at the same K
      (wall_face_unknowns.py, K = 6: 2.6e-4 at Co 0.0164, 4.6e-4 at Co 0.05).
  P6d the largest ||A^n||_2, n <= 2000, is <= 10.
  FALSIFIED if P6a or P6b fails at an edge or a corner.

RESULT (2026-11-26, laptop; MEASURED). Edge, K = 6, 12 waves along the edge: exactness
2.1e-14 (P6a PASS); growth for Co >= 0.1 at most 2.2e-16 (P6b PASS); Co < 0.1 at most
2.2e-4 (P6c PASS); largest ||A^n||_2 for n <= 2000: 4.76 (P6d PASS). Corner, K = 5:
exactness 4.3e-14 (P6a PASS); every arm decays (growth -4.4e-3 or less for Co >= 0.1,
-2.5e-4 or less for Co < 0.1; P6b, P6c PASS); P6d not run for the corner (343 unknowns,
2000 powers per arm).

Usage: python3 wall_edges_corners.py [--K 6] [--waves 12] [--domain edge|corner|both]
"""
import argparse, itertools
import numpy as np

def basis(s):
    s = np.atleast_2d(np.asarray(s, float))
    return np.stack([s[:, 0], s[:, 1], s[:, 2], .5*s[:, 0]**2, .5*s[:, 1]**2, .5*s[:, 2]**2,
                     s[:, 0]*s[:, 1], s[:, 0]*s[:, 2], s[:, 1]*s[:, 2]], axis=1)

def ls_weights(offsets, x):
    offsets = np.asarray(offsets, float); w = 1.0/np.sum(offsets**2, axis=1)
    B = basis(offsets); M = (B*w[:, None]).T @ B
    z = np.linalg.solve(M, basis(x)[0]); return w*(B @ z)

def lagrange(nodes, x):
    out = []
    for j, xj in enumerate(nodes):
        v = 1.0
        for m, xm in enumerate(nodes):
            if m != j: v *= (x - xm)/(xj - xm)
        out.append(v)
    return out

class Box:
    def __init__(self, K, bounded):
        self.K = K; self.bounded = bounded          # bounded[a]: walls on axis a
        self.index = {}; self.pos = []
        ranges = [list(range(K)) + [-0.5, K - 0.5] if bounded[a] else [0] for a in range(3)]
        for p in itertools.product(*ranges):
            self.index[p] = len(self.pos); self.pos.append(np.array(p, float))
        self.ncell = sum(1 for p in self.pos if not self.is_boundary_point(p))
    def is_wall(self, a, v): return self.bounded[a] and v in (-0.5, self.K - 0.5)
    def is_boundary_point(self, p): return any(self.is_wall(a, p[a]) for a in range(3))
    def inside(self, x):
        return all((not self.bounded[a]) or (-0.5 - 1e-12 <= x[a] <= self.K - 0.5 + 1e-12) for a in range(3))
    def key(self, p):
        """Index of an absolute lattice point and the periodic shift (for the Bloch phase)."""
        q = list(p); shift = [0, 0, 0]
        for a in range(3):
            if not self.bounded[a]: shift[a] = q[a]; q[a] = 0
        return self.index[tuple(float(v) if v != int(v) else int(v) for v in q)], shift
    def touches(self, c):
        return [a for a in range(3) if self.bounded[a] and c[a] in (0, self.K - 1)]
    def tensor(self, c, x):
        """Weights of the tensor interpolant of boundary cell c at displacement x."""
        nodes = []
        for a in range(3):
            if self.bounded[a] and c[a] == 0: nodes.append((-0.5, 0.0, 1.0))
            elif self.bounded[a] and c[a] == self.K - 1: nodes.append((0.5, 0.0, -1.0))
            else: nodes.append((-1.0, 0.0, 1.0))
        w = [lagrange(nodes[a], x[a]) for a in range(3)]
        out = []
        for i, j, k in itertools.product(range(3), repeat=3):
            p = c + np.array([nodes[0][i], nodes[1][j], nodes[2][k]])
            out.append((w[0][i]*w[1][j]*w[2][k], p))
        return out
    def interior(self, c, x):
        offs = [np.array(o, float) for o in itertools.product((-1, 0, 1), repeat=3) if any(o)]
        g = ls_weights(offs, x)
        return [(1 - g.sum(), c.copy())] + [(gs, c + o) for gs, o in zip(g, offs)]
    def host(self, p):
        return np.array([min(max(p[a], 0), self.K - 1) if self.bounded[a] else p[a] for a in range(3)], float)
    def rule(self, p, d):
        """Weights of the update of unknown p: list of (weight, absolute point); None = BC."""
        if not self.inside(p - d): return None
        if self.is_boundary_point(p):
            c = self.host(p); return self.tensor(c, p - d - c)
        if self.touches(p): return self.tensor(p, -d)
        return self.interior(p, -d)
    def matrix(self, d, kz=0.0):
        n = len(self.pos); A = np.zeros((n, n), complex)
        for r, p in enumerate(self.pos):
            ru = self.rule(p, d)
            if ru is None: continue
            for w, q in ru:
                idx, shift = self.key(tuple(q))
                A[r, idx] += w*np.exp(1j*kz*shift[2])
        return A
    def exactness(self, d):
        rng = np.random.default_rng(3); a = rng.normal(size=10)
        if not self.bounded[2]: a[[3, 6, 8, 9]] = 0          # no z: periodic
        q = lambda x: (a[0] + a[1]*x[0] + a[2]*x[1] + a[3]*x[2] + a[4]*x[0]**2 + a[5]*x[1]**2
                       + a[6]*x[2]**2 + a[7]*x[0]*x[1] + a[8]*x[0]*x[2] + a[9]*x[1]*x[2])
        err = 0.0
        for p in self.pos:
            ru = self.rule(p, d)
            if ru is None: continue
            err = max(err, abs(sum(w*q(qq) for w, qq in ru) - q(p - d)))
        return err

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--K", type=int, default=6); ap.add_argument("--waves", type=int, default=12)
    ap.add_argument("--domain", default="both"); ap.add_argument("--powers", action="store_true")
    a = ap.parse_args()
    dirs = {"x": (1, 0, 0), "z (along the edge)": (0, 0, 1), "(1 1 0)": (1, 1, 0), "(1 -1 0)": (1, -1, 0),
            "(1 1 1)": (1, 1, 1), "(-1 2 1)": (-1, 2, 1), "(0.05 0.05 1)": (0.05, 0.05, 1),
            "(1 0.05 0)": (1, 0.05, 0)}
    Cos = (0.0164, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
    doms = {"edge": (True, True, False), "corner": (True, True, True)}
    for name, bounded in doms.items():
        if a.domain not in ("both", name): continue
        K = a.K if name == "edge" else min(a.K, 5)
        box = Box(K, bounded); print(f"\n=== {name}: K = {K}, unknowns {len(box.pos)} "
                                     f"({box.ncell} cells, {len(box.pos) - box.ncell} boundary points)")
        ks = np.linspace(-np.pi, np.pi, a.waves, endpoint=False) if not bounded[2] else [0.0]
        worst_hi = -1; worst_lo = -1; worst_ex = 0; worst_p = 0
        for dn, e in dirs.items():
            e = np.asarray(e, float); e /= np.linalg.norm(e); row = []
            for Co in Cos:
                d = Co*e; g = -1.0; pm = 0.0
                for kz in ks:
                    A = box.matrix(d, kz); ev = np.abs(np.linalg.eigvals(A)).max(); g = max(g, ev - 1)
                    if a.powers:
                        P = np.eye(len(A)); best = 1.0
                        for n in range(1, 2001):
                            P = A @ P
                            if n < 50 or n % 20 == 0: best = max(best, np.linalg.norm(P, 2))
                        pm = max(pm, best)
                ex = box.exactness(d); worst_ex = max(worst_ex, ex); worst_p = max(worst_p, pm)
                if Co >= 0.1: worst_hi = max(worst_hi, g)
                else: worst_lo = max(worst_lo, g)
                row.append(f"Co {Co:g}: {g:+.2e}" + (f" P {pm:.2f}" if a.powers else ""))
            print(f"  {dn:16s} " + "  ".join(row))
        print(f"  {name}: exactness {worst_ex:.1e} (P6a <= 1e-12); growth Co >= 0.1: {worst_hi:+.2e} "
              f"(P6b <= 1e-10); Co < 0.1: {worst_lo:+.2e} (P6c <= slab 2.6e-4 / 4.6e-4)"
              + (f"; power norm {worst_p:.2f} (P6d <= 10)" if a.powers else ""))

if __name__ == "__main__":
    main()
