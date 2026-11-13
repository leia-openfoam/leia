#!/usr/bin/env python3
"""LaTeX tables of the scheme matrix, from data/tables/scheme_matrix_*.csv.

Reads the tables of scheme_matrix.py (the free variant, `_constrained` and `_ghost`) and
writes, into data/tables/:
  scheme_matrix_<mesh>.tex     one row per scheme: growth factor - 1 per Courant number,
                               the largest over the flow directions (and, on a perturbed
                               mesh, over the amplitudes and seeds)
  scheme_matrix_summary.tex    scheme x mesh, one verdict per cell

THE VERDICTS are the criteria that scheme_matrix.py pre-registered. Nothing is tuned here.
  r   the data of a cell do not carry a complete quadratic (rank-deficient)
  lattice (square, cube, honeycomb):   P if |growth - 1| <= 1e-6 at every Co <= 1, else F
  perturbed mesh:                      P if the scheme passes its lattice and its growth is
                                       not larger than that of the original scheme at any Co
                                       (C-Cp/L2; on the honeycomb C-Cf/L2), else F
  2:1 transition:                      P if growth - 1 <= 4.41e-3 at every Co;
                                       F if growth - 1 >= 1.03e-2 at any Co <= 0.5;
                                       ? between the two (the symbol does not decide)
  --  not run on this mesh (the scheme does not exist there, or the table has not landed)

Usage: scheme_matrix_tables.py [--tables DIR]
"""
import argparse, collections, csv, glob, os, sys
import numpy as np

MESHES = ["square", "quadPert", "quad21", "quad21patch", "cube", "hexPert", "hex21", "hex21patch",
          "honeycomb", "honeycombPert"]
TITLE = {"square": "square lattice", "quadPert": "perturbed quadrilaterals (5, 10 and 20 \\%)",
         "quad21": "quadrilaterals, planar 2:1 transition", "quad21patch": "quadrilaterals, refined patch",
         "cube": "cube lattice", "hexPert": "perturbed hexahedra (5, 10 and 20 \\%)",
         "hex21": "hexahedra, planar 2:1 transition", "hex21patch": "hexahedra, refined patch",
         "honeycomb": "honeycomb", "honeycombPert": "perturbed honeycomb (5, 10 and 20 \\%)"}
SHORT = {"square": "sq", "quadPert": "sq-p", "quad21": "sq-2:1", "quad21patch": "sq-patch", "cube": "cube",
         "hexPert": "cube-p", "hex21": "cube-2:1", "hex21patch": "cube-patch", "honeycomb": "hc", "honeycombPert": "hc-p"}
LATTICE = {"quadPert": "square", "hexPert": "cube", "honeycombPert": "honeycomb"}
PASS, FAIL = 4.41e-3, 1.03e-2

def load(tables):
    """{mesh: {scheme: dict(g=array of growth - 1 per Co, rank, full, deg, upc, err)}} and the Co list."""
    out = collections.defaultdict(dict); cos = None
    for f in sorted(glob.glob(os.path.join(tables, "scheme_matrix_*.csv"))):
        rows = list(csv.reader(l for l in open(f) if not l.startswith("#")))
        head, rows = rows[0], rows[1:]
        gi = [i for i, h in enumerate(head) if h.startswith("growth_Co")]
        cos = [float(head[i][len("growth_Co"):]) for i in gi]
        by = collections.defaultdict(list)
        for r in rows: by[(r[0], r[3])].append(r)
        for (mesh, scheme), R in by.items():
            g = np.array([[float(r[i]) for i in gi] for r in R])
            d = 3 if mesh in ("cube", "hexPert", "hex21", "hex21patch") else 2
            rank = min(int(r[9]) for r in R); full = {2: 6, 3: 10}[d]
            out[mesh][scheme] = dict(g=(np.max(g, axis=0) - 1) if not np.isnan(g).all() else np.full(len(gi), np.nan),
                                     rank=rank, full=full, deg=max(int(r[10]) for r in R),
                                     upc=max(int(r[5]) for r in R)/max(1, int(R[0][6])),
                                     err=max(float(r[gi[-1]+1]) for r in R))
    return out, cos

def verdict(D, mesh, scheme, cos):
    if scheme not in D.get(mesh, {}): return "--"
    e = D[mesh][scheme]; g = e["g"]; c = np.array(cos)
    if e["rank"] < e["full"] or np.isnan(g).any(): return "r"
    if mesh in ("square", "cube", "honeycomb"):
        return "P" if np.abs(g[c <= 1]).max() <= 1e-6 else "F"
    if mesh in LATTICE:
        if verdict(D, LATTICE[mesh], scheme, cos) != "P": return "F"
        ref = D[mesh]["C-Cf/L2" if mesh.startswith("honeycomb") else "C-Cp/L2"]["g"]
        return "P" if (g <= ref + 1e-9).all() else "F"
    if (g <= PASS).all(): return "P"
    if (g[c <= 0.5] >= FAIL).any(): return "F"
    return "?"

def num(x):
    if np.isnan(x): return "--"
    if abs(x) < 5e-7: return "0"
    m, e = f"{x:.1e}".split("e"); return f"{m}e{int(e)}"

def key(s):
    base = s.replace("+c", "").replace("+g", ""); sch, space = base.split("/"); own, _, nbr = sch.partition("-")
    return (["C", "SC", "BCC", "FCC", "BFC", "T"].index(own), ["", "Cf", "Cp"].index(nbr), space, s)

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--tables", default=os.path.join(here, "..", "data", "tables"))
    a = ap.parse_args(); D, cos = load(a.tables)
    if not D: sys.exit("no scheme_matrix_*.csv found")
    tt = lambda s: "\\texttt{" + s + "}"
    for mesh in MESHES:
        if mesh not in D: continue
        L = ["% generated by sympy/scheme_matrix_tables.py -- do not edit", "\\begin{tabular}{lrrr" + "r"*len(cos) + "rc}", "\\toprule",
             "scheme & unknowns & rank & degree & " + " & ".join(f"Co {c:g}" for c in cos) + " & one-step & verdict\\\\",
             " & per cell & & & " + " & "*len(cos) + "error & \\\\", "\\midrule"]
        for s in sorted(D[mesh], key=key):
            e = D[mesh][s]
            L.append(f"{tt(s)} & {e['upc']:.2f} & {e['rank']}/{e['full']} & {e['deg']} & "
                     + " & ".join(num(x) for x in e["g"]) + f" & {num(e['err'])} & {verdict(D, mesh, s, cos)}\\\\")
        L += ["\\bottomrule", "\\end{tabular}", ""]
        open(os.path.join(a.tables, f"scheme_matrix_{mesh}.tex"), "w").write("\n".join(L))
    schemes = sorted({s for m in D for s in D[m]}, key=key)
    L = ["% generated by sympy/scheme_matrix_tables.py -- do not edit", "\\begin{tabular}{l" + "c"*len(MESHES) + "}", "\\toprule",
         "scheme & " + " & ".join(SHORT[m] for m in MESHES) + "\\\\", "\\midrule"]
    last = None
    for s in schemes:
        if last is not None and key(s)[0] != last: L.append("\\addlinespace")
        last = key(s)[0]
        L.append(tt(s) + " & " + " & ".join(verdict(D, m, s, cos) for m in MESHES) + "\\\\")
    L += ["\\bottomrule", "\\end{tabular}", ""]
    open(os.path.join(a.tables, "scheme_matrix_summary.tex"), "w").write("\n".join(L))
    w = max(len(s) for s in schemes)
    print(f"{'scheme':{w}s} " + " ".join(f"{SHORT[m]:>10s}" for m in MESHES))
    for s in schemes: print(f"{s:{w}s} " + " ".join(f"{verdict(D, m, s, cos):>10s}" for m in MESHES))
    return 0

if __name__ == "__main__":
    sys.exit(main())
