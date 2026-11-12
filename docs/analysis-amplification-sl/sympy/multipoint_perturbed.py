#!/usr/bin/env python3
"""Point unknowns against the current scheme on a PERTURBED hexahedral mesh, Co <= 1.

WHY. On a perfect lattice both schemes have growth factor 1 for Co <= 1, so the perfect
lattice cannot separate them. The current scheme amplifies on a perturbed stencil
(sympy/perturbed_stencil.py). The question is whether the point-unknown scheme does better
on the meshes where the current one fails.

THE MESH. A block of P^d cells with randomly displaced vertices, repeated periodically. This
is a consistent mesh (every cell is a multilinear image of the unit cell), and it repeats, so
the wave analysis holds with one matrix per wave over the unknowns of the block (Bloch form).
The growth factor is the largest eigenvalue magnitude over all waves.

THE MINIMAL SUFFICIENT BLOCK IS P = 3. MEASURED 2026-11-12: with P = 2 every cell contains
ALL the distinct vertices of the block, so every cell centre is shifted by the same mean
displacement and the centres form a perfect lattice. The current scheme, which reads only
centres, then sees no perturbation and returns 1.000000 at every amplitude. P = 1 is a rigid
shift of the whole mesh. So P = 3 is the default, and `--block 2` reproduces the degenerate
result.

THE ARMS.
  current   one unknown per cell centre; anchored quadratic least-squares fit through the
            3^d - 1 neighbour centres, inverse-distance weights (the production fit).
  iso       point unknowns at vertices, edge midpoints, face centres and the cell centre
            (3^d points per cell). The departure point is located in its host cell, mapped to
            the reference cell by Newton iteration, and the tensor quadratic is interpolated
            there (the isoparametric Q2 element).
  phys      the same unknowns; the tensor quadratic x^a y^b z^c, a,b,c <= 2, is interpolated
            in PHYSICAL coordinates through the host's 3^d points. No Newton iteration.

PRE-REGISTERED (2026-11-12, before the run).
  P1. amplitude 0: all three arms give growth factor 1.000000 at every Co <= 1. This is the
      validity gate of the script: it must reproduce the perfect-lattice results.
  P2. amplitude 0.10: `iso` stays within 1e-3 of 1 at every Co, because the departure point
      is always enclosed by its host's data. `current` exceeds 1.005 at Co 0.3.
  FALSIFIED IF the growth factor of `iso` at amplitude 0.10 is equal to or larger than that
  of `current` at Co 0.3 and Co 0.5. Then the point unknowns do not help where the current
  scheme fails, and the candidate has no advantage in growth for Co <= 1.

Also reported: the error of one step for a long wave, RMS over the unknowns of the block.

Usage: multipoint_perturbed.py [--dim 2|3]
"""
import argparse, itertools, sys
import numpy as np

class Mesh:
    def __init__(self, d, P, amp, seed):
        self.d, self.P = d, P
        rng = np.random.default_rng(seed)
        self.pert = (rng.random((P,)*d + (d,)) - 0.5)*2*amp

    def vertex(self, m):
        m = np.asarray(m)
        return m + self.pert[tuple(np.mod(m, self.P))]

    def pos(self, n):
        """Position of node n of the doubled grid: the mean of the vertices of its entity."""
        n = np.asarray(n)
        opts = [[ni//2] if ni % 2 == 0 else [(ni-1)//2, (ni+1)//2] for ni in n]
        return np.mean([self.vertex(m) for m in itertools.product(*opts)], axis=0)

    def cell_vertices(self, k):
        return [(np.array(j), self.vertex(np.asarray(k)+np.array(j)))
                for j in itertools.product((0, 1), repeat=self.d)]

    def to_reference(self, k, p):
        """Newton inverse of the multilinear map of cell k. Returns xi in [-1/2,1/2]^d or None."""
        d = self.d; V = self.cell_vertices(k); xi = np.zeros(d)
        for _ in range(30):
            x = np.zeros(d); J = np.zeros((d, d))
            for j, X in V:
                f = [(0.5 + xi[i]) if j[i] else (0.5 - xi[i]) for i in range(d)]
                x += np.prod(f)*X
                for a in range(d):
                    g = np.prod([f[i] for i in range(d) if i != a])*(1.0 if j[a] else -1.0)
                    J[:, a] += g*X
            dx = np.linalg.solve(J, p - x); xi += dx
            if np.abs(dx).max() < 1e-13: break
        return xi if np.abs(xi).max() <= 0.5 + 1e-9 else None

    def host(self, p):
        k0 = np.floor(p).astype(int)
        for off in sorted(itertools.product((-1, 0, 1), repeat=self.d), key=lambda o: sum(map(abs, o))):
            k = k0 + np.array(off); xi = self.to_reference(k, p)
            if xi is not None: return k, xi
        raise RuntimeError("departure point not located")

NODES1 = (-0.5, 0.0, 0.5)
def lagrange(x, i):
    others = [NODES1[j] for j in range(3) if j != i]
    return np.prod([(x - o)/(NODES1[i] - o) for o in others])

def rows_point(mesh, c, arm):
    d, M = mesh.d, 2*mesh.P; rows = []
    J = list(itertools.product(range(3), repeat=d))
    for n in itertools.product(range(M), repeat=d):
        foot = mesh.pos(n) - c; k, xi = mesh.host(foot)
        nodes = [2*k + np.array(j) for j in J]
        if arm == "iso":
            w = [np.prod([lagrange(xi[i], j[i]) for i in range(d)]) for j in J]
        else:
            xc = mesh.pos(2*k + 1)
            mono = lambda X: np.array([np.prod(X**np.array(e)) for e in J])
            Vm = np.array([mono(mesh.pos(nd) - xc) for nd in nodes])
            w = np.linalg.solve(Vm.T, mono(foot - xc))
        rows.append((np.array(n), list(zip(nodes, w))))
    return rows, M

def rows_current(mesh, c):
    d, P = mesh.d, mesh.P; rows = []
    S = [np.array(s) for s in itertools.product((-1, 0, 1), repeat=d) if any(s)]
    mons = [e for e in itertools.product(range(3), repeat=d) if 0 < sum(e) <= 2]
    b = lambda X: np.array([np.prod(X**np.array(e)) for e in mons])
    for k in itertools.product(range(P), repeat=d):
        k = np.array(k); xc = mesh.pos(2*k + 1)
        D = np.array([mesh.pos(2*(k + s) + 1) - xc for s in S])
        w2 = np.sum(D**2, axis=1)**-1.0; B = np.array([b(x) for x in D])
        g = w2*(B@np.linalg.solve((B*w2[:, None]).T@B, b(-np.asarray(c))))
        rows.append((k, [(k, 1 - g.sum())] + [(k + s, gs) for s, gs in zip(S, g)]))
    return rows, P

def growth(rows, M, d, ntheta):
    N = M**d; flat = lambda loc: int(np.ravel_multi_index(tuple(loc), (M,)*d))
    ent = [(flat(r), flat(np.mod(n, M)), np.floor_divide(n, M), w) for r, lst in rows for n, w in lst]
    I = np.array([e[0] for e in ent]); Jc = np.array([e[1] for e in ent])
    SH = np.array([e[2] for e in ent], float); W = np.array([e[3] for e in ent])
    ax = np.linspace(-np.pi, np.pi, ntheta, endpoint=False)
    TH = np.array(list(itertools.product(ax, repeat=d))); worst = 0.0
    for chunk in np.array_split(TH, max(1, len(TH)*N*N//4_000_000)):
        A = np.zeros((len(chunk), N, N), complex)
        ph = W[None, :]*np.exp(1j*(chunk@SH.T))
        for t in range(len(chunk)): np.add.at(A[t], (I, Jc), ph[t])
        worst = max(worst, np.abs(np.linalg.eigvals(A)).max())
    return worst

def step_error(mesh, rows, c, position, kappa=0.2):
    kv = kappa*np.asarray(c)/np.linalg.norm(c); f = lambda x: np.sin(x@kv + 0.3)
    e = [sum(w*f(position(n)) for n, w in lst) - f(position(r) - c) for r, lst in rows]
    return float(np.sqrt(np.mean(np.square(e))))

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dim", type=int, default=2, choices=(2, 3))
    ap.add_argument("--block", type=int, default=3,
                    help="cells per direction of the periodic block. 3 is the minimum that "
                         "perturbs the cell centres; 2 is degenerate (see the header)")
    ap.add_argument("--ntheta", type=int, default=None, help="waves per direction of the block")
    ap.add_argument("--seeds", type=int, nargs="+", default=None)
    ap.add_argument("--cos", type=float, nargs="+", default=[0.0164, 0.1, 0.3, 0.5, 1.0])
    ap.add_argument("--amps", type=float, nargs="+", default=[0.0, 0.05, 0.10, 0.20])
    ap.add_argument("--csv", default=None, help="append the rows to this CSV")
    a = ap.parse_args(); d, P = a.dim, a.block
    seeds = tuple(a.seeds) if a.seeds else {2: (1, 2), 3: (1,)}[d]   # fixed: reproducible
    nth = a.ntheta or {2: 16, 3: 6}[d]   # waves per direction of the block
    dirs = [(1, 0), (1, 1), (2, 1)] if d == 2 else [(1, 0, 0), (1, 1, 1), (2, 1, 0)]
    COS = tuple(a.cos)
    print(f"{d}D, periodic block of {P}^{d} cells, worst over {len(seeds)} meshes "
          f"(seeds {seeds}) and {len(dirs)} flow directions\n")
    print(f"{'amplitude':>9} {'arm':>8} " + " ".join(f"{'Co '+str(c):>10}" for c in COS)
          + "   one-step error, Co 0.3")
    out = []
    for amp in a.amps:
        res = {arm: np.zeros(len(COS)) for arm in ("current", "iso", "phys")}
        err = {arm: 0.0 for arm in res}
        for seed in (seeds if amp else seeds[:1]):
            mesh = Mesh(d, P, amp, seed)
            for ic, Co in enumerate(COS):
                for v in dirs:
                    c = Co*np.array(v, float)/np.linalg.norm(v)
                    rc, Mc = rows_current(mesh, c)
                    res["current"][ic] = max(res["current"][ic], growth(rc, Mc, d, nth))
                    if Co == 0.3:
                        err["current"] = max(err["current"],
                                             step_error(mesh, rc, c, lambda k: mesh.pos(2*np.asarray(k)+1)))
                    for arm in ("iso", "phys"):
                        rp, Mp = rows_point(mesh, c, arm)
                        res[arm][ic] = max(res[arm][ic], growth(rp, Mp, d, nth))
                        if Co == 0.3:
                            err[arm] = max(err[arm], step_error(mesh, rp, c, mesh.pos))
        for arm in ("current", "iso", "phys"):
            print(f"{amp:9.2f} {arm:>8} " + " ".join(f"{x:10.6f}" for x in res[arm])
                  + f"   {err[arm]:.2e}", flush=True)
            out.append([d, P, amp, arm] + [f"{x:.8f}" for x in res[arm]] + [f"{err[arm]:.6e}"])
        print()
    if a.csv:
        import csv, os
        new = not os.path.exists(a.csv)
        with open(a.csv, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["dim", "block", "amplitude", "arm"] + [f"growth_Co{c}" for c in COS]
                           + ["step_error_Co0.3"])
            w.writerows(out)
        print(f"appended {len(out)} rows to {a.csv}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
