#!/usr/bin/env python3
"""Growth factor of every named semi-Lagrangian point scheme on every mesh type, Co <= 1.

THE NAMES.  OWN-NEIGHBOUR/SPACE.
  OWN        the points of a cell that carry advected unknowns
               C    cell centre                      SC   vertices (simple cubic)
               BCC  centre + vertices                FCC  vertices + face centres
               BFC  centre + vertices + face centres
               T    centre + vertices + face centres + edge midpoints (3D only; in 2D a
                    face is an edge, so BFC is the tensor set there)
  NEIGHBOUR  cell-centre values read from neighbour cells (needs C in OWN)
               -Cf  face neighbours (OpenFOAM cell-face-cell)
               -Cp  point neighbours (OpenFOAM cell-point-cell)
  SPACE      /I   the lowest-degree polynomial space that returns every datum exactly
             /L2  least-squares complete quadratic, weights 1/distance to the cell centre,
                  anchored on the centre value when the centre is a datum
  The original scheme is C-Cp/L2.

THE STEP.  Every unknown moves along its own trajectory. Its departure point is located in
its host cell, and the data of the host cell are interpolated there. Exception: a scheme
with OWN = C uses the cell of the unknown itself, as the original scheme does.

THE MESH.  A periodic block of polygonal or polyhedral cells. A coarse cell next to fine
cells carries the hanging nodes and the split faces, as OpenFOAM stores it. One step
multiplies every wave of the block by a matrix over the unknowns of the block (Bloch form).
The growth factor is the largest eigenvalue magnitude over all waves. The cell size of the
smallest cell is 1, and the Courant number is defined with it.

THE SPACE /I.  Degree by degree, the homogeneous polynomials of degree k are projected on
the complement of what the lower degrees already span on the data; the directions that
remain are added. The monomials carry the Bombieri scale sqrt(k!/(a! b! c!)), which makes
the choice invariant under rotation, so the space has the symmetry of the points. On a
perturbed mesh the space is selected on the unperturbed cell and used in physical
coordinates on the perturbed one.

PRE-REGISTERED (2026-11-12, before any run on a refined mesh or of a hybrid scheme).
  Validity gates (`--smoke`):
    G1. C-Cp/L2 and the tensor set (BFC/I in 2D, T/I in 3D) return 1.000000 on the perfect
        square and cube at every Co <= 1.
    G2. A block of 3 x 3 equal cells returns the value of the single-cell lattice.
    G3. With `--centre mean`, C-Cp/L2 and BFC/I on the perturbed quadrilateral reproduce
        the rows `current` and `phys` of data/tables/multipoint_perturbed.csv.
    G4. BCC/I on the honeycomb reproduces multipoint_hexagon.py (C+V interp) within 5e-5:
        the reference has five decimals and 120 waves per direction, this script uses 48.
  Predictions:
    P1. BCC/L2 has rank 8 of 10 on the cube. BCC-Cf has full rank.
  A scheme PASSES a mesh type if (a) the complete quadratic has full rank, (b) on the
  perfect lattice |growth - 1| <= 1e-6 for every Co <= 1, (c) on the perturbed block its
  growth is not larger than that of C-Cp/L2.
  A scheme PASSES the 2:1 transition if growth - 1 <= 4.41e-3 at every Co. It FAILS if
  growth - 1 >= 1.03e-2 at any Co <= 0.5. Between the two the symbol does not decide.

Usage: scheme_matrix.py --mesh NAME [--csv FILE] [--jobs N] | --smoke | --list
"""
import argparse, itertools, math, os, sys
# The matrices are small. BLAS threads on small matrices cost more than they give, and they
# multiply with --jobs (MEASURED 2026-11-12: 4 jobs x 16 threads ran 200 times slower).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
import numpy as np

OWN = {"C": "C", "SC": "V", "BCC": "CV", "FCC": "VF", "BFC": "CVF", "T": "CVFE"}
NUDGE = 1e-7*np.array([0.37, 0.61, 0.83])      # fixed tie-break for a point on a cell boundary

# ---------------------------------------------------------------------------- the mesh
class PMesh:
    """Periodic block. cells: list of cells, a cell is a list of faces, a face is an (n, d)
    array of reference coordinates (2D: n = 2; 3D: the polygon in order)."""
    def __init__(self, d, box, cells, name, lattice=None, centre="centroid"):
        self.d, self.box, self.ref, self.name = d, np.asarray(box, float), cells, name
        self.lattice, self.centre_mode = lattice, centre
        self.disp = {}
        self.build()

    def key(self, x):
        f = np.asarray(x)/self.box; s = np.floor(f + 1e-8).astype(int)
        return tuple(np.round((f - s)*1e6).astype(int)), s

    def perturb(self, amp, seed):
        keys = sorted({self.key(v)[0] for c in self.ref for f in c for v in f})
        rng = np.random.default_rng(seed)
        if self.lattice:                      # the same numbers as multipoint_perturbed.py
            P = self.lattice; pert = (rng.random((P,)*self.d + (self.d,)) - 0.5)*2*amp
            self.disp = {k: pert[tuple(np.round(np.array(k)/1e6*self.box).astype(int) % P)] for k in keys}
        else:
            pert = (rng.random((len(keys), self.d)) - 0.5)*2*amp
            self.disp = dict(zip(keys, pert))
        self.build()

    def act(self, P):
        return np.array([v + self.disp.get(self.key(v)[0], 0.0) for v in P])

    def fcentre(self, P):
        if self.d == 2 or self.centre_mode == "mean": return P.mean(0)
        m = P.mean(0); a = P - m; b = np.roll(P, -1, axis=0) - m
        w = np.linalg.norm(np.cross(a, b), axis=1)
        return ((m + P + np.roll(P, -1, axis=0))/3*w[:, None]).sum(0)/w.sum()

    def simplices(self, faces, c):
        out = []
        for P in faces:
            if self.d == 2: out.append(np.array([c, P[0], P[1]]))
            else:
                f = self.fcentre(P)
                out += [np.array([c, f, P[i], P[(i+1) % len(P)]]) for i in range(len(P))]
        return out

    def ccentre(self, faces):
        V = np.unique(np.round(np.vstack(faces), 12), axis=0); m = V.mean(0)
        if self.centre_mode == "mean": return m
        S = self.simplices(faces, m)
        w = np.array([abs(np.linalg.det(s[1:] - s[0])) for s in S])
        return (np.array([s.mean(0) for s in S])*w[:, None]).sum(0)/w.sum()

    def build(self):
        d = self.d; self.cells = []; self.uid = {t: {} for t in "CVFE"}; self.upos = {t: [] for t in "CVFE"}
        def reg(t, xref, xact):
            k, s = self.key(xref)
            if k not in self.uid[t]:
                self.uid[t][k] = len(self.upos[t]); self.upos[t].append(xact - s*self.box)
            return (self.uid[t][k], s, xact)
        vmap, fmap = {}, {}
        for j, fr in enumerate(self.ref):
            fa = [self.act(P) for P in fr]
            cr, ca = self.ccentre(fr), self.ccentre(fa)
            pts = {"C": [reg("C", cr, ca)], "V": [], "F": [], "E": []}; seen = set()
            def add(t, xr, xa):
                k, s = self.key(xr)
                if (t, k, tuple(s)) in seen: return
                seen.add((t, k, tuple(s))); pts[t].append(reg(t, xr, xa))
                if t == "V": vmap.setdefault(k, []).append((j, s))
                if t == "F": fmap.setdefault(k, []).append((j, s))
            for Pr, Pa in zip(fr, fa):
                for xr, xa in zip(Pr, Pa): add("V", xr, xa)
                add("F", self.fcentre(Pr), self.fcentre(Pa))
                if d == 3:
                    for i in range(len(Pr)):
                        i2 = (i+1) % len(Pr); add("E", (Pr[i]+Pr[i2])/2, (Pa[i]+Pa[i2])/2)
            self.cells.append({"pts": pts, "cref": cr, "c": ca, "faces": fa,
                               "href": max(np.linalg.norm(v - cr) for P in fr for v in P)})
        # neighbours, with their periodic shift
        inv = {t: {i: k for k, i in self.uid[t].items()} for t in "VF"}
        for j, c in enumerate(self.cells):
            for t, mp, nm in (("V", vmap, "Cp"), ("F", fmap, "Cf")):
                nb = set()
                for (u, s, _) in c["pts"][t]:
                    for (j2, s2) in mp[inv[t][u]]:
                        sh = tuple(int(x) for x in (s - s2))
                        if j2 != j or any(sh): nb.add((j2, sh))
                c[nm] = sorted(nb)
        S, O, cid = [], [], []
        for j, c in enumerate(self.cells):
            for s in self.simplices(c["faces"], c["c"]):
                S.append(np.linalg.inv((s[1:] - s[0]).T)); O.append(s[0]); cid.append(j)
        self.Tinv, self.O, self.scell = np.array(S), np.array(O), np.array(cid)
        self.shifts = np.array(list(itertools.product((-1, 0, 1), repeat=d)))

    def locate(self, p):
        q = p + NUDGE[:self.d]
        X = q[None, None, :] - (self.shifts*self.box)[:, None, :] - self.O[None]
        lam = np.einsum("nij,mnj->mni", self.Tinv, X)
        mb = np.minimum(lam.min(-1), 1 - lam.sum(-1))
        m, n = np.unravel_index(np.argmax(mb), mb.shape)
        if mb[m, n] < -1e-6: raise RuntimeError(f"departure point {p} not located ({mb[m, n]:.2e})")
        return int(self.scell[n]), self.shifts[m]

def box_mesh(d, boxes, L, name, **kw):
    """Axis-aligned cubes of size 1 or 2 -> cells with hanging nodes and split faces."""
    L = np.asarray(L, float)
    def key(x):
        f = np.asarray(x)/L; s = np.floor(f + 1e-8); return tuple(np.round((f - s)*1e6).astype(int))
    allv = {key(np.array(o) + sz*np.array(j)) for o, sz in boxes for j in itertools.product((0, 1), repeat=d)}
    def with_mid(P):
        out = []
        for i in range(len(P)):
            a, b = P[i], P[(i+1) % len(P)]; out.append(a)
            if key((a+b)/2) in allv: out.append((a+b)/2)
        return out
    cells = []
    for o, sz in boxes:
        o = np.array(o, float); e = np.eye(d)*sz
        if d == 2:
            poly = with_mid([o, o+e[0], o+e[0]+e[1], o+e[1]])
            cells.append([np.array([poly[i], poly[(i+1) % len(poly)]]) for i in range(len(poly))]); continue
        faces = []
        for ax in range(3):
            u, v = [a for a in range(3) if a != ax]
            for side in (0, 1):
                b = o + side*e[ax]; c = [b, b+e[u], b+e[u]+e[v], b+e[v]]; ctr = b + (e[u]+e[v])/2
                quads = [c] if key(ctr) not in allv else \
                        [[c[i], (c[i]+c[(i+1) % 4])/2, ctr, (c[i-1]+c[i])/2] for i in range(4)]
                faces += [np.array(with_mid(q)) for q in quads]
        cells.append(faces)
    return PMesh(d, L, cells, name, **kw)

def lattice_mesh(d, P, name, **kw):
    return box_mesh(d, [(o, 1) for o in itertools.product(range(P), repeat=d)], (P,)*d, name, lattice=P, **kw)

def strip_mesh(d, name, **kw):
    """x: 3 coarse cells (size 2), then 4 fine cells (size 1). Tangential period: one coarse cell."""
    t = d - 1; boxes = [((2*i,) + (0,)*t, 2) for i in range(3)]
    boxes += [((6+i,) + j, 1) for i in range(4) for j in itertools.product((0, 1), repeat=t)]
    return box_mesh(d, boxes, (10,) + (2,)*t, name, **kw)

def patch_mesh(d, name, **kw):
    """4^d coarse cells (size 2), the central 2^d of them refined: transition edges and corners."""
    boxes = []
    for o in itertools.product(range(4), repeat=d):
        if all(1 <= x <= 2 for x in o):
            boxes += [(tuple(2*np.array(o) + np.array(j)), 1) for j in itertools.product((0, 1), repeat=d)]
        else: boxes.append((tuple(2*np.array(o)), 2))
    return box_mesh(d, boxes, (8,)*d, name, **kw)

def honeycomb_mesh(nx, ny, name, **kw):
    r = 1/np.sqrt(3); h = np.sqrt(3)
    ring = [r*np.array([np.cos(np.radians(30+60*k)), np.sin(np.radians(30+60*k))]) for k in range(6)]
    cells = []
    for i in range(nx):
        for j in range(ny):
            for c in (np.array([i, j*h]), np.array([i+0.5, (j+0.5)*h])):
                cells.append([np.array([c+ring[k], c+ring[(k+1) % 6]]) for k in range(6)])
    return PMesh(2, (nx, ny*h), cells, name, **kw)

AX2 = [(1, 0), (1, 1), (2, 1)]; AX3 = [(1, 0, 0), (1, 1, 1), (2, 1, 0)]
TR2 = [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2)]
TR3 = [(1, 0, 0), (0, 1, 0), (1, 1, 1), (2, 1, 0), (1, 2, 0)]
HEX = [(np.cos(np.radians(a)), np.sin(np.radians(a))) for a in (0, 15, 30)]
PERT = [(0.05, 1), (0.05, 2), (0.10, 1), (0.10, 2), (0.20, 1), (0.20, 2)]
MESHES = {   # name: (builder, flow directions, [(amplitude, seed)])
    "square":        (lambda **k: lattice_mesh(2, 1, "square", **k), AX2, [(0, 0)]),
    "squareBlock":   (lambda **k: lattice_mesh(2, 3, "squareBlock", **k), AX2, [(0, 0)]),
    "quadPert":      (lambda **k: lattice_mesh(2, 3, "quadPert", **k), AX2, PERT),
    "quad21":        (lambda **k: strip_mesh(2, "quad21", **k), TR2, [(0, 0)]),
    "quad21patch":   (lambda **k: patch_mesh(2, "quad21patch", **k), TR2, [(0, 0)]),
    "cube":          (lambda **k: lattice_mesh(3, 1, "cube", **k), AX3, [(0, 0)]),
    "hexPert":       (lambda **k: lattice_mesh(3, 3, "hexPert", **k), AX3, [(0.05, 1), (0.10, 1), (0.20, 1)]),
    "hex21":         (lambda **k: strip_mesh(3, "hex21", **k), TR3, [(0, 0)]),
    "hex21patch":    (lambda **k: patch_mesh(3, "hex21patch", **k), TR3, [(0, 0)]),
    "honeycomb":     (lambda **k: honeycomb_mesh(1, 1, "honeycomb", **k), HEX, [(0, 0)]),
    "honeycombPert": (lambda **k: honeycomb_mesh(3, 2, "honeycombPert", **k), HEX, PERT),
}

# ---------------------------------------------------------------------------- the spaces
def exps(d, k):
    return [e for e in itertools.product(range(k+1), repeat=d) if sum(e) == k]

def mono(X, E, bombieri=True):
    X = np.atleast_2d(X); k = sum(E[0])
    sc = [math.sqrt(math.factorial(k)/np.prod([math.factorial(a) for a in e])) if bombieri else 1.0 for e in E]
    return np.stack([s*np.prod(X**np.array(e), axis=1) for e, s in zip(E, sc)], axis=1)

def least_space(X, tol=1e-8, kmax=16):
    """Lowest-degree space that interpolates on the points X. Returns ([(E, coef)], degree, unique)."""
    N, d = X.shape; Q = np.zeros((N, 0)); fs = []; unique = True
    for k in range(kmax+1):
        E = exps(d, k); B = mono(X, E); R = B - Q@(Q.T@B)
        U, s, Vt = np.linalg.svd(R, full_matrices=False)
        r = int((s > tol).sum()); need = N - Q.shape[1]
        if r > need:
            if s[need] > 1e-3*s[need-1] and abs(s[need]-s[need-1]) < 1e-6: unique = False
            r = need
        fs += [(E, Vt[i]) for i in range(r)]
        Q = np.linalg.qr(np.hstack([Q, U[:, :r]]))[0] if r else Q
        if Q.shape[1] == N: return fs, k, unique
    raise RuntimeError("no interpolating space found")

def feval(fs, X):
    return np.stack([mono(X, E)@c for E, c in fs], axis=1)

class CellRule:
    """The data of one host cell for one scheme, and the weights of a point in it."""
    def __init__(self, mesh, j, own, nbr, space):
        d = mesh.d; c = mesh.cells[j]; self.xc = c["c"]; self.h = c["href"]; self.space = space
        data = []                                   # (type, uid, shift, x actual, x reference-local)
        for t in own:
            for (u, s, xa) in c["pts"][t]:
                data.append((t, u, s, xa))
        if nbr:
            for (j2, sh) in c[nbr]:
                data.append(("C", j2, np.array(sh), mesh.cells[j2]["c"] + np.array(sh)*mesh.box))
        self.data = data; X = (np.array([x for *_, x in data]) - self.xc)/self.h
        self.anch = "C" in own
        ncoef = len(exps(d, 1)) + len(exps(d, 2)) + 1
        full = np.hstack([mono(X, exps(d, k), False) for k in range(3)])
        self.rank = int(np.linalg.matrix_rank(full, tol=1e-9)); self.full = self.rank == ncoef
        self.deg, self.unique, self.ok = 2, True, True
        if space == "I":
            # the space is chosen on the unperturbed cell
            Xr = self._ref(mesh, j, own, nbr)
            self.fs, self.deg, self.unique = least_space(Xr)
            self.B = feval(self.fs, X); self.ok = np.linalg.cond(self.B) < 1e12
        else:
            self.ok = self.full
            if self.anch:
                D = X[1:]; self.W = 1/np.sum(D**2, axis=1)
                self.Bq = np.hstack([mono(D, exps(d, k), False) for k in (1, 2)])
            else:
                self.W = 1/np.sum(X**2, axis=1); self.Bq = full
            if self.ok: self.Minv = np.linalg.inv((self.Bq*self.W[:, None]).T@self.Bq)

    def _ref(self, mesh, j, own, nbr):
        if not mesh.disp:
            return (np.array([x for *_, x in self.data]) - self.xc)/self.h
        m0 = mesh._ref_mesh; c0 = m0.cells[j]
        X = [xa for t in own for (_, _, xa) in c0["pts"][t]]
        if nbr: X += [m0.cells[j2]["c"] + np.array(sh)*m0.box for (j2, sh) in c0[nbr]]
        return (np.array(X) - c0["c"])/c0["href"]

    def weights(self, p):
        d = len(p); x = ((p - self.xc)/self.h)[None, :]
        if self.space == "I":
            return np.linalg.solve(self.B.T, feval(self.fs, x)[0])
        if self.anch:
            b = np.hstack([mono(x, exps(d, k), False) for k in (1, 2)])[0]
            g = self.W*(self.Bq@(self.Minv@b)); return np.concatenate([[1 - g.sum()], g])
        b = np.hstack([mono(x, exps(d, k), False) for k in range(3)])[0]
        return self.W*(self.Bq@(self.Minv@b))

# ---------------------------------------------------------------------------- the symbol
def parse(name):
    sch, space = name.split("/"); own, _, nbr = sch.partition("-")
    return own, nbr, space

def scheme_list(d, honeycomb=False):
    owns = ["C", "SC", "BCC", "FCC", "BFC"] + (["T"] if d == 3 else [])
    out = []
    for o in owns:
        for n in ("", "Cf", "Cp"):
            if n and "C" not in OWN[o]: continue
            if o == "C" and not n: continue
            if honeycomb and n == "Cp": continue      # the same six cells as Cf
            out += [f"{o}{'-'+n if n else ''}/{s}" for s in ("I", "L2")]
    return out

def assemble(mesh, name, c, rules):
    own, nbr, space = parse(name); types = OWN[own]
    off = {}; n = 0
    for t in types:
        off[t] = n; n += len(mesh.upos[t])
    ent = []; pos = []; npts = []
    for t in types:
        for u, x in enumerate(mesh.upos[t]):
            p = x - c
            j, s = (u, np.zeros(mesh.d, int)) if own == "C" else mesh.locate(p)
            if j not in rules: rules[j] = CellRule(mesh, j, types, nbr, space)
            R = rules[j]
            if not R.ok: return None
            w = R.weights(p - s*mesh.box); npts.append(len(w))
            lst = [(off[tt] + uu, ss + s, xx + s*mesh.box, ww) for (tt, uu, ss, xx), ww in zip(R.data, w)]
            ent.append((off[t] + u, lst)); pos.append(x)
    return n, ent, np.array(pos)

def growth(mesh, n, ent, target):
    I = np.array([r for r, lst in ent for _ in lst]); J = np.array([e[0] for _, lst in ent for e in lst])
    SH = np.array([e[1] for _, lst in ent for e in lst], float); W = np.array([e[3] for _, lst in ent for e in lst])
    nth = [max(2, int(round(target/L))) for L in mesh.box]
    TH = np.array(list(itertools.product(*[np.linspace(-np.pi, np.pi, m, endpoint=False) for m in nth])))
    TH = TH[[i for i, t in enumerate(TH) if tuple(t) >= tuple(-t) or np.any(np.isclose(np.abs(t), np.pi))]]
    flat = I*n + J; worst = 0.0
    for ch in np.array_split(TH, max(1, len(TH)*n*n*16//200_000_000)):
        ph = W[None, :]*np.exp(1j*(ch@SH.T)); A = np.zeros((len(ch), n*n), complex)
        for t in range(len(ch)): A[t] = np.bincount(flat, ph[t].real, n*n) + 1j*np.bincount(flat, ph[t].imag, n*n)
        worst = max(worst, np.abs(np.linalg.eigvals(A.reshape(len(ch), n, n))).max())
    return worst

def step_error(ent, pos, c, kappa):
    kv = kappa*c/np.linalg.norm(c); f = lambda x: np.sin(x@kv + 0.3)
    e = [sum(w*f(x) for _, _, x, w in lst) - f(pos[i] - c) for i, (_, lst) in enumerate(ent)]
    return float(np.sqrt(np.mean(np.square(e))))

MESH = None
def run_one(task):
    name, amp, seed, COS, dirs, target = task; mesh = MESH[(amp, seed)]
    own, nbr, space = parse(name); rows = []
    for v in dirs:
        v = np.array(v, float); g = []; info = None; rules = {}
        for Co in COS + (0.3,):
            c = Co*v/np.linalg.norm(v); a = assemble(mesh, name, c, rules)
            if info is None:
                R = list(rules.values())
                info = (min(r.rank for r in R), max(r.deg for r in R), int(all(r.unique for r in R)))
            if a is None: g.append(float("nan")); continue
            n, ent, pos = a
            if len(g) < len(COS): g.append(growth(mesh, n, ent, target))
            else: err = [step_error(ent, pos, c, k) for k in (0.2, 0.1)]
        if a is None: err = [float("nan")]*2; n = 0; npt = (0, 0)
        else: npt = (min(len(l) for _, l in ent), max(len(l) for _, l in ent))
        rows.append([mesh.name, amp, seed, name, "(%s)" % " ".join(f"{x:.4f}" for x in v), n, len(mesh.cells),
                     npt[0], npt[1], info[0], info[1], info[2]] + [f"{x:.8f}" for x in g[:len(COS)]]
                    + [f"{err[0]:.6e}", f"{err[1]:.6e}"])
    return rows

def run(mesh_name, schemes=None, COS=(0.0164, 0.1, 0.3, 0.5, 1.0), centre="centroid", jobs=1, csv_path=None,
        cases=None, quiet=False):
    global MESH
    build, dirs, pert = MESHES[mesh_name]; cases = cases or pert; MESH = {}
    for amp, seed in cases:
        m = build(centre=centre)
        if amp: m._ref_mesh = build(centre=centre); m.perturb(amp, seed)
        MESH[(amp, seed)] = m
    m0 = next(iter(MESH.values())); d = m0.d
    schemes = schemes or scheme_list(d, mesh_name.startswith("honeycomb"))
    target = {2: 48, 3: 24}[d]
    tasks = [(s, amp, seed, tuple(COS), dirs, target) for amp, seed in cases for s in schemes]
    if jobs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(jobs) as pool: res = pool.map(run_one, tasks, chunksize=1)
    else: res = [run_one(t) for t in tasks]
    rows = [r for rr in res for r in rr]
    head = ["mesh", "amplitude", "seed", "scheme", "direction", "unknowns", "cells", "pts_min", "pts_max",
            "rank_P2", "degree", "unique"] + [f"growth_Co{c}" for c in COS] + ["err_k0.2_Co0.3", "err_k0.1_Co0.3"]
    if csv_path:
        import csv
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(head); w.writerows(rows); fh.write("# done\n")
    if not quiet: summary(rows, COS, d)
    return rows

def summary(rows, COS, d):
    full = {2: 6, 3: 10}[d]; nC = len(COS)
    print(f"{'mesh':14s} {'amp':>5} {'scheme':10s} {'unk/cell':>8} {'pts':>7} {'rank':>5} {'deg':>3} "
          + " ".join(f"{'Co '+str(c):>10}" for c in COS) + "   one-step error k=0.2, order")
    keys = []
    for r in rows:
        k = (r[0], r[1], r[3])
        if k not in keys: keys.append(k)
    for k in keys:
        R = [r for r in rows if (r[0], r[1], r[3]) == k]
        g = np.nanmax(np.array([[float(x) for x in r[12:12+nC]] for r in R]), axis=0) if R[0][5] else [float("nan")]*nC
        e2 = max(float(r[12+nC]) for r in R); e1 = max(float(r[13+nC]) for r in R)
        order = np.log2(e2/e1) if e1 > 0 else float("nan")
        upc = R[0][5]/R[0][6] if R[0][5] else 0
        print(f"{k[0]:14s} {k[1]:5.2f} {k[2]:10s} {upc:8.2f} {str(R[0][7])+'-'+str(R[0][8]):>7} "
              f"{R[0][9]:>2}/{full:<2} {R[0][10]:>3} " + " ".join(f"{x:10.6f}" for x in g)
              + f"   {e2:.2e} {order:5.2f}" + ("" if R[0][11] else "  (space not unique)"), flush=True)

def smoke():
    ok = True
    def worst(rows, scheme, n=5):
        return np.nanmax(np.array([[float(x) for x in r[12:12+n]] for r in rows if r[3] == scheme]), axis=0)
    def check(label, cond):
        nonlocal ok; ok &= bool(cond); print(f"  {'PASS' if cond else 'FAIL'}  {label}")
    r = run("square", ["C-Cp/L2", "BFC/I"], quiet=True)
    check("G1 square: C-Cp/L2 and BFC/I are 1.000000", all(abs(worst(r, s) - 1).max() < 1e-6 for s in ("C-Cp/L2", "BFC/I")))
    r = run("cube", ["C-Cp/L2", "T/I", "BCC/L2", "BCC-Cf/L2"], quiet=True)
    check("G1 cube: C-Cp/L2 and T/I are 1.000000", all(abs(worst(r, s) - 1).max() < 1e-6 for s in ("C-Cp/L2", "T/I")))
    rk = {x[3]: x[9] for x in r}
    check(f"P1 cube: rank of BCC/L2 is {rk['BCC/L2']} (predicted 8), of BCC-Cf/L2 is {rk['BCC-Cf/L2']} (predicted 10)",
          rk["BCC/L2"] == 8 and rk["BCC-Cf/L2"] == 10)
    sl = ["C-Cp/L2", "BFC/I", "BCC-Cf/L2", "FCC/I"]
    a = run("square", sl, quiet=True); b = run("squareBlock", sl, quiet=True)
    check("G2 block of equal cells equals the single-cell lattice",
          all(abs(worst(a, s) - worst(b, s)).max() < 1e-7 for s in sl))
    here = os.path.dirname(os.path.abspath(__file__))
    ref = os.path.join(here, "..", "data", "tables", "multipoint_perturbed.csv")
    if os.path.exists(ref):
        import csv
        T = {(x["dim"], float(x["amplitude"]), x["arm"]): [float(x[k]) for k in x if k.startswith("growth")]
             for x in csv.DictReader(open(ref))}
        r = run("quadPert", ["C-Cp/L2", "BFC/I"], centre="mean", cases=[(0.10, 1), (0.10, 2)], quiet=True)
        dev = max(abs(worst(r, s) - np.array(T[("2", 0.1, arm)])).max() for s, arm in (("C-Cp/L2", "current"), ("BFC/I", "phys")))
        check(f"G3 perturbed quadrilateral reproduces multipoint_perturbed.csv (largest deviation {dev:.1e})", dev < 2e-6)
    r = run("honeycomb", ["BCC/I"], COS=(0.0164, 0.1, 0.25, 0.5, 0.75, 1.0), quiet=True)
    refv = np.array([1.00030, 1.00127, 1.00148, 1.00035, 1.00020, 1.00024]); dev = abs(worst(r, "BCC/I", 6) - refv).max()
    check(f"G4 honeycomb BCC/I reproduces multipoint_hexagon.py C+V interp (largest deviation {dev:.1e})", dev < 5e-5)
    return 0 if ok else 1

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mesh", choices=sorted(MESHES)); ap.add_argument("--schemes", nargs="+")
    ap.add_argument("--cos", type=float, nargs="+", default=[0.0164, 0.1, 0.3, 0.5, 1.0])
    ap.add_argument("--centre", choices=("centroid", "mean"), default="centroid")
    ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--csv")
    ap.add_argument("--smoke", action="store_true"); ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.smoke: return smoke()
    if a.list or not a.mesh:
        for m in MESHES: print(m)
        return 0
    run(a.mesh, a.schemes, tuple(a.cos), a.centre, a.jobs, a.csv)
    return 0

if __name__ == "__main__":
    sys.exit(main())
