#!/usr/bin/env python3
"""Independent point unknowns on a regular HEXAGONAL tiling, for Courant numbers up to 1.

WHY. sympy/multipoint_unknowns.py showed that on cubes a complete point set gives growth
factor 1. A hexagon has no tensor point set, so the question is which combination of its
natural points -- centre (C), 6 edge midpoints (E), 6 vertices (V) -- works, and with which
interpolation space.

WHAT THIS MEASURES. Each unknown is moved along its own trajectory. Its departure point is
located in the host hexagon, and the host's own point values are interpolated there. Per
unit cell there are 1 centre, 3 edge midpoints and 2 vertices. One step multiplies a wave by
a matrix; the growth factor is its largest eigenvalue magnitude over all waves.

Two spaces per combination:
  interp   the lowest-degree polynomial space that returns every point value exactly. It is
           built from the functions r^(k-m) cos(m t), r^(k-m) sin(m t) in order of degree k,
           keeping each one that raises the rank. The space is then symmetric by construction.
  lsqP2    the least-squares fit of a complete quadratic (6 coefficients), where the
           combination has more than 6 points.
The baseline is the current scheme: centres only, quadratic fit through the 6 neighbours.

Usage: multipoint_hexagon.py
"""
import itertools, sys
import numpy as np

A1 = np.array([1.0, 0.0]); A2 = np.array([0.5, np.sqrt(3)/2]); LAT = np.array([A1, A2]).T
ang = lambda deg, r: r*np.array([np.cos(np.radians(deg)), np.sin(np.radians(deg))])
PTS = {"C": [np.zeros(2)],
       "E": [ang(60*k, 0.5) for k in range(6)],
       "V": [ang(30+60*k, 1/np.sqrt(3)) for k in range(6)]}
TYPES = {"C": [np.zeros(2)], "E": [A1/2, A2/2, (A2-A1)/2], "V": [(A1+A2)/3, (2*A2-A1)/3]}

def harmonics(kmax=8):
    fs = []
    for k in range(kmax+1):
        for m in range(k % 2, k+1, 2):
            def re(X, k=k, m=m):
                z = X[:,0]+1j*X[:,1]; return (np.abs(z)**(k-m))*(z**m).real
            def im(X, k=k, m=m):
                z = X[:,0]+1j*X[:,1]; return (np.abs(z)**(k-m))*(z**m).imag
            fs.append((k, re));
            if m: fs.append((k, im))
    return fs

def space(Q, kind):
    Q = np.array(Q)
    if kind == "lsqP2":
        fs = [f for k, f in harmonics(2)]
        B = np.stack([f(Q) for f in fs], axis=1)
        return fs, np.linalg.pinv(B), 2, np.linalg.matrix_rank(B)
    chosen, cols, deg = [], [], 0
    for k, f in harmonics():
        col = f(Q)
        if np.linalg.matrix_rank(np.stack(cols+[col], axis=1), tol=1e-9) > len(cols):
            cols.append(col); chosen.append(f); deg = k
        if len(cols) == len(Q): break
    B = np.stack(cols, axis=1)
    return chosen, np.linalg.inv(B), deg, np.linalg.cond(B)

def host_of(p):
    m = np.round(np.linalg.solve(LAT, p)); best = None
    for d in itertools.product((-1,0,1), repeat=2):
        k = m+np.array(d); c = LAT@k
        if best is None or np.linalg.norm(p-c) < best[0]-1e-12: best = (np.linalg.norm(p-c), k, c)
    return best[1], best[2]

def classify(p, tlist):
    for i, (name, pos) in enumerate(tlist):
        L = np.linalg.solve(LAT, p-pos)
        if np.allclose(L, np.round(L), atol=1e-8): return i, np.round(L)
    raise RuntimeError("point not on the lattice of unknowns")

def symbol_rows(combo, kind, c):
    tlist = [(n, pos) for n in combo for pos in TYPES[n]]
    Q = [q for n in combo for q in PTS[n]]
    fs, Binv, deg, info = space(Q, kind)
    rows = []
    for name, pos in tlist:
        foot = pos - c; k, cen = host_of(foot)
        w = np.array([f((foot-cen).reshape(1,2))[0] for f in fs])@Binv
        rows.append([(*classify(cen+q, tlist), wq) for q, wq in zip(Q, w)])
    return len(tlist), rows, deg, info, len(Q)

def growth(combo, kind, c, n=120):
    nt, rows, *_ = symbol_rows(combo, kind, c)
    ax = np.linspace(-np.pi, np.pi, n, endpoint=False)
    TH = np.stack(np.meshgrid(ax, ax, indexing="ij"), axis=-1).reshape(-1,2)
    A = np.zeros((len(TH), nt, nt), complex)
    for i, ent in enumerate(rows):
        for j, L, w in ent: A[:,i,j] += w*np.exp(1j*(TH@L))
    return np.abs(np.linalg.eigvals(A)).max()

def step_error(combo, kind, c, eps):
    nt, rows, *_ = symbol_rows(combo, kind, c)
    kv = eps*np.array(c)/np.linalg.norm(c); th = LAT.T@kv
    A = np.zeros((nt, nt), complex)
    for i, ent in enumerate(rows):
        for j, L, w in ent: A[i,j] += w*np.exp(1j*(th@L))
    ev = np.linalg.eigvals(A); ex = np.exp(-1j*(kv@c))
    return abs(ev[np.argmin(abs(ev-ex))]-ex)

def baseline(c, n=240, eps=None):
    S = np.array([ang(60*k, 1.0) for k in range(6)]); I = np.round(np.linalg.solve(LAT, S.T).T)
    b = lambda X: np.stack([X[:,0], X[:,1], .5*X[:,0]**2, .5*X[:,1]**2, X[:,0]*X[:,1]], axis=1)
    B = b(S); g = B@np.linalg.solve(B.T@B, b(-np.array(c).reshape(1,2))[0])
    if eps is not None:
        kv = eps*np.array(c)/np.linalg.norm(c)
        return abs((1-g.sum())+(g*np.exp(1j*(I@(LAT.T@kv)))).sum()-np.exp(-1j*(kv@c)))
    ax = np.linspace(-np.pi, np.pi, n, endpoint=False)
    TH = np.stack(np.meshgrid(ax, ax, indexing="ij"), axis=-1).reshape(-1,2)
    return np.abs((1-g.sum())+(np.exp(1j*TH@I.T)*g).sum(axis=1)).max()

def main():
    COS = (0.0164, 0.1, 0.25, 0.5, 0.75, 1.0); DIRS = (0, 15, 30)
    worst = lambda f: [max(f(ang(d, Co)) for d in DIRS) for Co in COS]
    print("largest growth factor over flow directions 0, 15, 30 degrees; Courant number <= 1\n")
    print(f"{'points':8s} {'space':7s} {'unk/cell':>8} {'pts':>4} {'deg':>4} "
          + " ".join(f"{'Co '+str(c):>9}" for c in COS) + "   one-step error (Co 0.3)")
    r = worst(baseline); e = [baseline(ang(0, 0.3), eps=x) for x in (0.2, 0.1)]
    print(f"{'current':8s} {'P2 nbrs':7s} {1:8d} {6:4d} {2:4d} " + " ".join(f"{x:9.5f}" for x in r)
          + f"   {e[0]:.1e} {e[1]:.1e} (x{e[0]/e[1]:.1f})")
    for combo in ("CV", "CE", "VE", "CVE"):
        for kind in ("interp", "lsqP2"):
            nt, _, deg, info, npts = symbol_rows(combo, kind, ang(0, 0.3))
            if kind == "lsqP2" and npts <= 6: continue
            r = worst(lambda c: growth(combo, kind, c))
            e = [step_error(combo, kind, ang(0, 0.3), x) for x in (0.2, 0.1)]
            print(f"{'+'.join(combo):8s} {kind:7s} {nt:8d} {npts:4d} {deg:4d} "
                  + " ".join(f"{x:9.5f}" for x in r) + f"   {e[0]:.1e} {e[1]:.1e} (x{e[0]/e[1]:.1f})")
    for combo in ("V", "E"):
        Q = np.array(PTS[combo]); B = np.stack([f(Q) for k, f in harmonics(2)], axis=1)
        print(f"{combo:8s} complete quadratic on its 6 points: rank {np.linalg.matrix_rank(B)} of 6"
              " (x^2 + y^2 is constant on a circle)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
