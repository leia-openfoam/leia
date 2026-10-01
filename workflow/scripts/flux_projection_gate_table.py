"""The flux-projection gate (STATUS 11.25): fluxProjection none (A) against helmholtz (B) per rung,
the read-out at T and T/2 with consecutive orders. Usage: flux_projection_gate_table.py <study dir>"""
import csv, glob, math, os, re, sys
d = sys.argv[1]
def rows(p):
    with open(p) as f:
        return list(csv.DictReader(f))
def ncells(case):
    t = open(os.path.join(case, "constant", "polyMesh", "owner"), errors="replace").read()
    m = re.search(r"\n(\d+)\s*\n\(", t); s = m.end(); e = t.index(")", s)
    return max(int(x) for x in t[s:e].split()) + 1
def at(rs, t):
    return min(rs, key=lambda r: abs(float(r["TIME"]) - t))
res = []
for a in sorted(glob.glob(os.path.join(d, "*_none"))):
    b = a[:-5] + "_helmholtz"
    ra, rb = rows(os.path.join(a, "leiaSemiLagrangeLevelSetFoam.csv")), rows(os.path.join(b, "leiaSemiLagrangeLevelSetFoam.csv"))
    ga, gb = rows(os.path.join(a, "gradPsiError.csv")), rows(os.path.join(b, "gradPsiError.csv"))
    T = float(ra[-1]["TIME"]); Tb = float(rb[-1]["TIME"])
    n = ncells(a)
    r = dict(name=os.path.basename(a)[:-5], n=n, h=(1.0/n)**(1/3), steps=(len(ra)-1, len(rb)-1), T=(T, Tb))
    for key, col, src, tt in (("geomT", "E_GEOM_ALPHA_REL", "m", 1.0), ("volH", "E_VOL_ALPHA_REL", "m", 0.5),
                              ("volT", "E_VOL_ALPHA_REL", "m", 1.0), ("gradBandH", "E_NARROW_L2_GRAD_PSI", "g", 0.5),
                              ("boundT", "E_BOUND_ALPHA", "m", 1.0)):
        A = ra if src == "m" else ga; B = rb if src == "m" else gb
        r[key] = (float(at(A, tt*T)[col]), float(at(B, tt*Tb)[col]))
    res.append(r)
res.sort(key=lambda r: -r["h"])
def order(e1, e2, h1, h2):
    return math.log(e1/e2)/math.log(h1/h2) if e1 > 0 and e2 > 0 else float("nan")
for r in res:
    print(f"{r['name']}: cells {r['n']}, h_eff {r['h']:.5f}, steps none/helmholtz {r['steps'][0]}/{r['steps'][1]}, T {r['T'][0]:g}/{r['T'][1]:g}")
for key, lab in (("geomT", "E_GEOM_ALPHA_REL at T"), ("volH", "E_VOL_ALPHA_REL at T/2"), ("volT", "E_VOL_ALPHA_REL at T"),
                 ("gradBandH", "E_NARROW_L2_GRAD_PSI at T/2"), ("boundT", "E_BOUND_ALPHA at T")):
    print(lab)
    for i, r in enumerate(res):
        x, y = r[key]
        ch = (y - x)/x*100 if x else float("nan")
        line = f"  h={r['h']:.5f}  none {x:.4e}  helmholtz {y:.4e}  change {ch:+.1f} %"
        if i:
            q = res[i-1]
            line += f"   order none {order(q[key][0], x, q['h'], r['h']):.3f}  helmholtz {order(q[key][1], y, q['h'], r['h']):.3f}"
        print(line)
