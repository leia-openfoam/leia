#!/usr/bin/env python3
"""The +g schemes on meshes with THREE cell sizes (two refinement levels above the base mesh).

THE MESHES.  `quad3lev_b<k>` (2D) and `hex3lev_b<k>` (3D), k = 1, 2, 3, built by
`three_level_mesh` of scheme_matrix.py. The periodic block holds, in the direction normal
to the transitions:
    3 cells of size 4 | k cells of size 2 | 4 cells of size 1 | k cells of size 2
The tangential period is one cell of size 4. This is the convention of `strip_mesh` (3
coarse cells, 4 fine cells, tangential period one coarse cell). The middle cell of size 4
has only neighbours of its own size. The two middle columns of size 1 have only neighbours
of their own size. A cell of size 2 has only neighbours of its own size for k = 3 only.
The Courant number is defined with the smallest cell (size 1), as on every mesh.
    unknowns per block: 19 + 4 k (2D), 67 + 8 k (3D)
    waves: as on the other meshes, `target` points per unit length of the reciprocal block
The buffer of the solver design is k = 2.

HOW THE GHOST VALUES DEPEND ON EACH OTHER.  See the section THREE CELL SIZES in the header
of scheme_matrix.py. A prolongation reads only larger cells, so the dependency ends at the
largest size. `assemble_g` records the longest chain in `mesh.hcache["depth"]`.

PRE-REGISTERED (2026-11-23, before any run on a mesh with three cell sizes).
  Validity gates (`--gates`). No three-level row is read before the gates pass.
    G8.  The generalised `assemble_g` reproduces the committed tables
         data/tables/scheme_matrix_{quad21,hex21,quad21patch}_ghost.csv within 1e-8 in every
         growth column of every row. (`--gates-hex21patch` adds the hexahedral patch.)
    G9.  On the blocks `quad3lev_u` and `hex3lev_u` (the constructor of the three-level
         mesh with the three sizes equal, k = 2) C-Cp/L2+g and C-Cp/I+g return
         |growth - 1| <= 1e-6 at every Co <= 1, and every +g scheme returns the value of
         the scheme without the suffix on the same block within 1e-9 (the form of G7).
         (The first form of this gate compared C-Cf/I+g with the single-cell lattice. It
         failed by 1.8e-4 and 5.2e-3: C-Cf/I amplifies, and a block of 11 cells samples
         other waves than a block of 1 cell. The comparison on the same block replaces
         it. No mesh with three cell sizes had run at that time.)
    G10. For k = 2 and k = 3 no cell has both a larger and a smaller cell among its
         cell-point-cell neighbours, and the longest prolongation chain is 1. The gate
         asserts it. For k = 1 the script reports the count of such cells and the chain
         length; that is a result, not a gate.
    G11. On every three-level mesh, C-Cp/L2+g and C-Cp/I+g return a LINEAR function
         exactly: for every cell, the weighted sum of the linear function over the cell
         values of its row equals the function at the departure point, within 1e-11, at
         Co 0.3 and 1.0 in the direction (1, 2[, 0]). Both ghost rules are exact for a
         linear function, so this checks the positions, the periodic shifts and the
         weights of the chained ghost values. It reads no growth factor.
  PREDICTION P5, for C-Cp/L2+g, growth factor - 1, the largest of the five flow directions,
  at Co 0.0164 / 0.1 / 0.3 / 0.5 / 1.0. The one-transition values (measured 2026-11-12):
      2D (quad21)  4.3231e-4 / 1.88607e-3 / 2.94209e-3 / 3.15372e-3 / 3.27854e-3
      3D (hex21)   4.6401e-4 / 2.17126e-3 / 3.91455e-3 / 4.38432e-3 / 3.68499e-3
    P5a. Buffer k = 2: growth - 1 <= 4.41e-3 at every one of the five Co, in 2D and in 3D.
    P5b. Buffer k = 2: growth - 1 <= 1.2 times the one-transition value of the same
         dimension, at every one of the five Co.
    P5c. Buffer k = 1: growth - 1 is larger than with k = 2 at Co 0.3 and at Co 0.5, in 2D
         and in 3D.
    P5d. Buffer k = 3: growth - 1 is within 10 % of the k = 2 value at every one of the five
         Co (|g3 - g2| <= 0.1 g2).
    Reason: with k >= 2 every cell at a transition sees one transition only, so the block
    is two one-transition problems that are 2 cells apart. With k = 1 the cell of size 2
    reads a restriction on one side and a prolongation on the other side.
    FALSIFIED IF k = 2 exceeds 4.41e-3 at any of the five Co in 2D or in 3D (P5a), or if
    k = 1 is not larger than k = 2 at Co 0.3 or at Co 0.5 (P5c). P5b and P5d are reported
    as HOLDS or MISSED; a miss does not falsify P5.
    The criterion of the 2:1 transition stays: PASS if growth - 1 <= 4.41e-3 at every Co,
    FAIL if growth - 1 >= 1.03e-2 at any Co <= 0.5, between the two the symbol does not
    decide. The script never reports an order.
  PREDICTION P6, power-boundedness (added 2026-11-23 on the coordinator's request, before
  any three-level run and before any run of `--powers`). The spectral radius does not
  bound the powers of a non-normal matrix (Lax and Richtmyer; Ferretti and Mehrenberger
  give a 2 x 2 matrix with its eigenvalues in the unit disk and unbounded powers).
  `--powers` therefore computes, for C-Cp/L2+g and C-Cp/I+g, per flow direction and Co:
      P_max = max over the Bloch waves of the block and over n = 1..N of ||A_theta^n||_2,
  N = 2000, the 2-norm = the largest singular value. The norm is taken at every n <= 200
  and at every 10th n after that (the cost of 2000 singular value decompositions per wave
  does not fit the laptop limit in 3D). The script records P_max, the n of the maximum,
  and the maximum over n <= 1000. Meshes: quad21, hex21, and the six three-level meshes.
    P6a. Where growth factor - 1 <= 1e-12 (neutral or damped): P_max < 10.
    P6b. Where the growth factor g is above 1: ||A^n|| grows like g^n, so the maximum is at
         n = N; the script reports P_max with its n and the ratio P_max / g^N.
    (The tables hold the growth factor with 8 decimals, so "<= 1e-12" reads as "the table
    value is at most 1.00000000".)
    FALSIFIED IF a case with growth factor - 1 <= 1e-12 has a P_max that grows with N:
    P_max(N = 2000) > 1.05 P_max(N = 1000) (algebraic growth from a defective eigenvalue).
    The columns pmax_Co*, npmax_Co*, pmaxhalf_Co* are APPENDED to the tables; no existing
    column changes. A row of another scheme gets empty fields.
  WHICH CELL'S POLYNOMIAL (read from the code, 2026-11-23). The free variant of a scheme
  with OWN = C (`assemble`: `j, s = (u, 0) if own == "C"`) and the +g variant
  (`assemble_g`: `stencil(j, x, x - c)`) evaluate the polynomial of the ARRIVAL cell at
  the departure point, also when the departure point leaves the arrival cell. This is what
  the solver does. Only the schemes with point unknowns locate the HOST cell
  (`mesh.locate`). So the +g tables need no arrival-cell option.
  A CHECK OF THE MEASUREMENT (`--waves`): the block is 20 to 28 long, so the wave sampling
  of the other meshes gives 2 waves in the normal direction. The script repeats C-Cp/L2+g
  with 4 times the waves per direction and reports the difference. A difference above 5 %
  of growth - 1 means that the table values are a lower bound, and the report says so.

Usage: scheme_matrix_three_levels.py --gates | --run [--dim 2 3] [--jobs N] | --powers [--dim 2 3] | --verdict | --waves [--dim 2 3] [--mesh NAME ...] | --wide | --scan
"""
import argparse, csv, os, sys
import scheme_matrix as sm          # first: it limits the BLAS threads before numpy loads
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(HERE, "..", "data", "tables")
COS = (0.0164, 0.1, 0.3, 0.5, 1.0)
PREFIX = {2: "quad", 3: "hex"}
ONE = {2: "quad21", 3: "hex21"}
PASS, FAIL = 4.41e-3, 1.03e-2

def read(path):
    rows = list(csv.reader(l for l in open(path) if not l.startswith("#")))
    head, rows = rows[0], rows[1:]; gi = [i for i, h in enumerate(head) if h.startswith("growth_Co")]
    return {(r[3], r[4]): np.array([float(r[i]) for i in gi]) for r in rows}

def table(mesh):
    return os.path.join(TABLES, f"scheme_matrix_{mesh}_ghost.csv")

def worst(T, scheme):
    """growth - 1, the largest over the directions, per Co; and the direction of the largest."""
    D = [(k[1], v) for k, v in T.items() if k[0] == scheme]
    G = np.array([v for _, v in D]) - 1
    return G.max(0), [D[i][0] for i in G.argmax(0)]

def gates(with_hex21patch=False):
    ok = True
    def check(label, cond):
        nonlocal ok; ok &= bool(cond); print(f"  {'PASS' if cond else 'FAIL'}  {label}", flush=True)
    for m in ["quad21", "quad21patch", "hex21"] + (["hex21patch"] if with_hex21patch else []):
        ref = read(table(m)); rows = sm.run(m, variant="ghost", quiet=True)
        new = {(r[3], r[4]): np.array([float(x) for x in r[12:17]]) for r in rows}
        same = set(new) == set(ref) and all(np.array_equal(np.isnan(new[k]), np.isnan(ref[k])) for k in ref)
        dev = max(np.nanmax(np.abs(new[k] - ref[k]), initial=0.0) for k in ref) if same else float("nan")
        check(f"G8 {m}: the committed +g table is reproduced (largest deviation {dev:.1e}, {len(ref)} rows)",
              same and dev <= 1e-8)
    for u in ("quad3lev_u", "hex3lev_u"):
        sl = ["C-Cp/L2", "C-Cp/I", "C-Cf/I"]
        r = sm.run(u, [x + "+g" for x in sl], quiet=True); a = sm.run(u, sl, quiet=True)
        w = lambda rows, s: np.nanmax(np.array([[float(x) for x in q[12:17]] for q in rows if q[3] == s]), axis=0)
        dev = max(abs(w(r, x + "+g") - 1).max() for x in sl[:2]); dfr = max(abs(w(r, x + "+g") - w(a, x)).max() for x in sl)
        check(f"G9 {u} ({r[0][6]} equal cells): C-Cp/L2+g and C-Cp/I+g are 1 (largest deviation {dev:.1e}); "
              f"every +g scheme equals the scheme without the suffix ({dfr:.1e})", dev <= 1e-6 and dfr <= 1e-9)
    for d in (2, 3):
        for k in (1, 2, 3):
            name = f"{PREFIX[d]}3lev_b{k}"; mesh = sm.MESHES[name][0]()
            size = [float(np.ptp(np.vstack(P)[:, 0])) for P in mesh.plain]
            both = [j for j, c in enumerate(mesh.cells)
                    if any(size[j2] > size[j] for j2, _ in c["Cp"]) and any(size[j2] < size[j] for j2, _ in c["Cp"])]
            sm.assemble_g(mesh, "C-Cp/L2+g", 0.3*np.ones(d)/np.sqrt(d), {}); depth = mesh.hcache["depth"]
            n = {s: size.count(s) for s in sorted(set(size), reverse=True)}
            msg = (f"{name}: box {tuple(int(x) for x in mesh.box)}, cells per size {n}, {len(both)} cell(s) with a larger "
                   f"AND a smaller cell-point-cell neighbour (sizes {sorted({size[j] for j in both})}), "
                   f"longest prolongation chain {depth}")
            lin = 0.0; a = np.array([0.7, -0.4, 0.3])[:d]
            for sch in ("C-Cp/L2+g", "C-Cp/I+g"):
                for Co in (0.3, 1.0):
                    c = Co*np.array(sm.MESHES[name][1][-1], float); c /= np.linalg.norm(c)/Co
                    _, ent, pos = sm.assemble_g(mesh, sch, c, {})
                    lin = max(lin, max(abs(sum(w*(x@a) for _, _, x, w in lst) - (pos[i] - c)@a) for i, (_, lst) in enumerate(ent)))
            check(f"G11 {name}: a linear function is returned exactly (largest deviation {lin:.1e})", lin <= 1e-11)
            if k == 1: print(f"  NOTE  G10 {msg}")
            else: check(f"G10 {msg}", not both and depth == 1)
    return ok

def run(dims, jobs):
    for d in dims:
        for k in (1, 2, 3):
            name = f"{PREFIX[d]}3lev_b{k}"; print(f"==== {name}", flush=True)
            sm.run(name, variant="ghost", jobs=jobs, csv_path=table(name))

def finest(d, name, scheme="C-Cp/L2+g"):
    """The row of the finest wave sampling of `--waves`, in the form of worst()."""
    R = [r for r in csv.reader(open(os.path.join(TABLES, f"scheme_matrix_three_levels_waves_{d}d.csv"))) if r[:2] == [name, scheme]]
    r = max(R, key=lambda r: np.prod([int(x) for x in r[2].split("x")]))
    return np.array([float(x) for x in r[3:8]]), r[8:13]

def verdict(fine=False):
    falsified = False
    for d in (2, 3):
        if not all(os.path.exists(table(f"{PREFIX[d]}3lev_b{k}")) for k in (1, 2, 3)):
            print(f"{d}D: the tables have not landed"); continue
        if fine:
            one = finest(d, ONE[d])[0]; g = {k: finest(d, f"{PREFIX[d]}3lev_b{k}") for k in (1, 2, 3)}
        else:
            one = worst(read(table(ONE[d])), "C-Cp/L2+g")[0]
            g = {k: worst(read(table(f"{PREFIX[d]}3lev_b{k}")), "C-Cp/L2+g") for k in (1, 2, 3)}
        print(f"{d}D, C-Cp/L2+g, growth - 1, the largest of the flow directions, {'the finest wave sampling' if fine else 'the sampling of the tables'}, Co " + " / ".join(map(str, COS)))
        print(f"  one transition : " + " / ".join(f"{x:.3e}" for x in one))
        for k in (1, 2, 3):
            crit = "PASS" if g[k][0].max() <= PASS else ("FAIL" if g[k][0][:4].max() >= FAIL else "NOT DECIDED")
            print(f"  buffer k = {k}   : " + " / ".join(f"{x:.3e}" for x in g[k][0]) + f"   criterion: {crit}")
            print(f"     worst direction: " + " / ".join(g[k][1]))
        a = g[2][0].max() <= PASS; b = np.all(g[2][0] <= 1.2*one)
        c = g[1][0][2] > g[2][0][2] and g[1][0][3] > g[2][0][3]
        e = np.all(np.abs(g[3][0] - g[2][0]) <= 0.1*g[2][0])
        print(f"  P5a k = 2 at most 4.41e-3 at every Co: {'HOLDS' if a else 'FALSIFIED'} (largest {g[2][0].max():.3e})")
        print(f"  P5b k = 2 at most 1.2 times one transition: {'HOLDS' if b else 'MISSED'} "
              f"(ratios " + " / ".join(f"{x:.3f}" for x in g[2][0]/one) + ")")
        print(f"  P5c k = 1 larger than k = 2 at Co 0.3 and 0.5: {'HOLDS' if c else 'FALSIFIED'} "
              f"(ratios k1/k2 " + " / ".join(f"{x:.3f}" for x in g[1][0]/g[2][0]) + ")")
        print(f"  P5d k = 3 within 10 % of k = 2: {'HOLDS' if e else 'MISSED'} "
              f"(ratios k3/k2 " + " / ".join(f"{x:.3f}" for x in g[3][0]/g[2][0]) + ")")
        falsified |= not (a and c)
    print("P5: " + ("FALSIFIED" if falsified else "not falsified (see the lines above for each part)"))

SAMPLING = {2: [(2, 12), (8, 12), (2, 48), (8, 48), (8, 96), (8, 192)],        # waves normal, tangential (period 4)
            3: [(2, 6, 6), (4, 6, 6), (2, 12, 12), (2, 24, 24), (4, 24, 24)]}
def waves(dims, wide=False, only=None):
    """C-Cp/L2+g and C-Cp/I+g with more waves per direction. The first sampling is the one
    of the tables (for the three-level meshes). The tangential count is per period 4; a mesh
    with the tangential period 2 gets half the count, which gives the same wave numbers.
    wide: the buffers k = 4 and 6 (NOT pre-registered, exploratory)."""
    for d in dims:
        names = [ONE[d]] + [f"{PREFIX[d]}3lev_b{k}" for k in (1, 2, 3)]
        if wide: names = [f"{PREFIX[d]}3lev_b{k}" for k in (4, 6)]
        if only: names = [n for n in names if n in only]
        path = os.path.join(TABLES, f"scheme_matrix_three_levels_waves{'_wide' if wide else ''}_{d}d.csv")
        new = not (only and os.path.exists(path))              # --mesh appends to the table (one run per mesh)
        with open(path, "w" if new else "a") as fh:
            if new: fh.write("mesh,scheme,waves," + ",".join(f"growth_minus_1_Co{c}" for c in COS) + "," + ",".join(f"direction_Co{c}" for c in COS) + "\n")
            for name in names:
                if name in sm.MESHES: build, dirs, _ = sm.MESHES[name]; mesh = build()
                else: mesh = sm.three_level_mesh(d, int(name.split("_b")[1]), name); dirs = sm.TR2 if d == 2 else sm.TR3
                for sch in PSCH:
                    A = {}
                    for v in dirs:
                        v = np.array(v, float); rules = {}
                        for i, Co in enumerate(COS): A[(tuple(v), i)] = sm.assemble_g(mesh, sch, Co*v/np.linalg.norm(v), rules)
                    for nth in SAMPLING[d][-2:] if wide else SAMPLING[d]:
                        nt = (nth[0],) + tuple(max(2, int(round(x*mesh.box[1]/4))) for x in nth[1:])
                        g = np.full(len(COS), -1.0); wd = [""]*len(COS)
                        for (v, i), (n, ent, pos) in A.items():
                            x = sm.growth(mesh, n, ent, nt) - 1
                            if x > g[i]: g[i] = x; wd[i] = "(%s)" % " ".join(f"{y:g}" for y in v)
                        fh.write(f"{name},{sch},{'x'.join(map(str, nt))}," + ",".join(f"{x:.6e}" for x in g) + "," + ",".join(wd) + "\n"); fh.flush()
                        print(f"{name:12s} {sch:10s} {'x'.join(map(str, nt)):9s} " + " / ".join(f"{x:.3e}" for x in g) + "   " + " ".join(wd), flush=True)

def scan():
    """NOT pre-registered, exploratory, 2D, C-Cp/L2+g: how the growth factor depends on the
    widths of the strips of the periodic block. One transition (sizes 2 and 1) and three
    levels. Waves: 4 normal, 48 per tangential period 2 (96 per period 4)."""
    path = os.path.join(TABLES, "scheme_matrix_three_levels_scan_2d.csv")
    cases = [(0, nc, nf) for nc in (2, 3, 4, 5, 6, 8) for nf in (4, 6, 8, 12, 16)]
    cases += [(k, nc, nf) for k in (1, 2, 3, 4, 5, 6, 8, 10) for nc, nf in ((3, 4), (3, 8), (4, 4), (5, 8), (4, 12))]
    with open(path, "w") as fh:
        fh.write("sizes,buffer_k,cells_coarsest,cells_finest,unknowns," + ",".join(f"growth_minus_1_Co{c}" for c in COS) + "\n")
        for k, nc, nf in cases:
            sizes = (4, 2, 1) if k else (2, 2, 1); mesh = sm.three_level_mesh(2, k, "scan", sizes, (nc, nf)); g = np.full(len(COS), -1.0)
            for v in sm.TR2:
                v = np.array(v, float); rules = {}
                for i, Co in enumerate(COS):
                    n, ent, pos = sm.assemble_g(mesh, "C-Cp/L2+g", Co*v/np.linalg.norm(v), rules)
                    g[i] = max(g[i], sm.growth(mesh, n, ent, (4, 96 if k else 48)) - 1)
            line = f"{'4-2-1' if k else '2-1'},{k},{nc},{nf},{n}," + ",".join(f"{x:.6e}" for x in g)
            fh.write(line + "\n"); fh.flush(); print(line, flush=True)

PSCH = ("C-Cp/L2+g", "C-Cp/I+g")
def powers(meshes, N=2000):
    """Append the power-bound columns to the +g tables. The existing text of every row stays."""
    for name in meshes:
        raw = open(table(name), newline="").read().splitlines(keepends=True)      # keep the line ends of the table
        lines = [l.rstrip("\r\n") for l in raw]; end = [l[len(l.rstrip("\r\n")):] for l in raw]; head = lines[0].split(",")
        if any(h.startswith("pmax_") for h in head): print(f"{name}: the columns exist"); continue
        build, dirs, _ = sm.MESHES[name]; mesh = build(); d = mesh.d; out = []
        gi = [i for i, h in enumerate(head) if h.startswith("growth_Co")]
        out.append(lines[0] + "," + ",".join(f"{c}_Co{x}" for c in ("pmax", "npmax", "pmaxhalf") for x in COS))
        for l in lines[1:]:
            if l.startswith("#"): out.append(l); continue
            r = l.split(",")
            if r[3] not in PSCH: out.append(l + ","*(3*len(COS))); continue
            v = np.array([float(x) for x in r[4].strip("()").split()]); P = []; rules = {}
            for Co in COS:
                n, ent, pos = sm.assemble_g(mesh, r[3], Co*v/np.linalg.norm(v), rules)
                P.append(sm.power_bound(mesh, n, ent, {2: 48, 3: 24}[d], N))
            out.append(l + "," + ",".join([f"{p[0]:.6e}" for p in P] + [str(p[1]) for p in P] + [f"{p[2]:.6e}" for p in P]))
            g = np.array([float(r[i]) for i in gi])
            print(f"{name} {r[3]} {r[4]}: P_max " + " / ".join(f"{p[0]:.3e}@{p[1]}" for p in P), flush=True)
        open(table(name), "w", newline="").write("".join(o + e for o, e in zip(out, end)))

def powers_long(meshes, N=16000):
    """The neutral cases that miss P6 at N = 2000 (P_max >= 10, or more than 5 % above the
    value at N/2), again with N = 16000. A defective eigenvalue of magnitude 1 gives a
    P_max that doubles with N. Writes scheme_matrix_three_levels_powers_long.csv."""
    path = os.path.join(TABLES, "scheme_matrix_three_levels_powers_long.csv"); new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new: fh.write("mesh,scheme,direction,Co,growth," + ",".join(f"pmax_N{N >> i}" for i in (3, 2, 1, 0)) + ",npmax\n")
        for name in meshes:
            rows = list(csv.DictReader(l for l in open(table(name)) if not l.startswith("#")))
            mesh = sm.MESHES[name][0](); d = mesh.d
            for q in rows:
                if q["scheme"] not in PSCH: continue
                for c in COS:
                    pq, ph = float(q[f"pmax_Co{c}"]), float(q[f"pmaxhalf_Co{c}"])
                    if float(q[f"growth_Co{c}"]) - 1 > 1e-12 or not (pq >= 10 or pq > 1.05*ph): continue
                    v = np.array([float(x) for x in q["direction"].strip("()").split()])
                    n, ent, pos = sm.assemble_g(mesh, q["scheme"], c*v/np.linalg.norm(v), {})
                    P = [sm.power_bound(mesh, n, ent, {2: 48, 3: 24}[d], N >> i) for i in (3, 2, 1, 0)]
                    line = f"{name},{q['scheme']},{q['direction']},{c},{q[f'growth_Co{c}']}," + ",".join(f"{p[0]:.6e}" for p in P) + f",{P[-1][1]}"
                    fh.write(line + "\n"); fh.flush(); print(line, flush=True)

def powers_verdict(meshes, N=2000):
    bad = False
    for name in meshes:
        rows = list(csv.DictReader(l for l in open(table(name)) if not l.startswith("#")))
        if "pmax_Co0.1" not in rows[0]: print(f"{name}: no power columns"); continue
        for s in PSCH:
            R = [r for r in rows if r["scheme"] == s]
            print(f"{name} {s}: Co, growth - 1, P_max (n of the maximum), P_max/growth^N, direction; the direction of the largest P_max")
            for c in COS:
                r = max(R, key=lambda r: float(r[f"pmax_Co{c}"])); g = float(r[f"growth_Co{c}"]); p = float(r[f"pmax_Co{c}"])
                print(f"   Co {c}: {g-1: .3e}  {p:.4e} ({r[f'npmax_Co{c}']})  {p/max(g, 1)**N:.3f}  {r['direction']}")
                for q in R:
                    if float(q[f"growth_Co{c}"]) - 1 <= 1e-12:
                        pq, ph = float(q[f"pmax_Co{c}"]), float(q[f"pmaxhalf_Co{c}"])
                        if pq >= 10 or pq > 1.05*ph:
                            bad = True; print(f"   P6 MISS: neutral case {q['direction']} Co {c}: P_max {pq:.4e}, to N/2 {ph:.4e}")
    print("P6: " + ("a neutral case misses, see above" if bad else "no neutral case has P_max >= 10 or growth with N"))

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gates", action="store_true"); ap.add_argument("--gates-hex21patch", action="store_true")
    ap.add_argument("--run", action="store_true"); ap.add_argument("--verdict", action="store_true")
    ap.add_argument("--waves", action="store_true"); ap.add_argument("--powers", action="store_true"); ap.add_argument("--wide", action="store_true")
    ap.add_argument("--scan", action="store_true"); ap.add_argument("--powers-long", action="store_true"); ap.add_argument("--mesh", nargs="+")
    ap.add_argument("--dim", type=int, nargs="+", default=[2, 3]); ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args()
    if a.gates or a.gates_hex21patch:
        if not gates(a.gates_hex21patch): return 1
    if a.run: run(a.dim, a.jobs)
    if a.waves or a.wide: waves(a.dim, a.wide, a.mesh)
    if a.scan: scan()
    pm = [m for m in ["quad21", "hex21"] + [f"{PREFIX[d]}3lev_b{k}" for d in (2, 3) for k in (1, 2, 3)]
          if (3 if m.startswith("hex") else 2) in a.dim and os.path.exists(table(m))]
    if a.powers: powers([m for m in pm if not a.mesh or m in a.mesh])
    if a.powers_long: powers_long([m for m in pm if not a.mesh or m in a.mesh])
    if a.verdict: verdict(); verdict(True); powers_verdict(pm)
    return 0

if __name__ == "__main__":
    sys.exit(main())
