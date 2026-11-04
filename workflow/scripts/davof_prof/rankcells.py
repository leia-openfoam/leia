"""per-rank cell and interface-cell counts from a decomposed case (ascii or binary alpha fields)"""
import sys, re, os, glob, struct
import numpy as np
def read_field(path):
    b = open(path, "rb").read()
    h = b[:2000].decode(errors="replace")
    fmt = re.search(r"format\s+(\w+);", h).group(1)
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", b)
    if not m:
        m2 = re.search(rb"internalField\s+uniform\s+([^;]+);", b)
        return None, float(m2.group(1)) if m2 else None
    n = int(m.group(1)); o = m.end()
    if fmt == "binary":
        return np.frombuffer(b, dtype="<f8", count=n, offset=o), None
    s = b[o:].split(b")", 1)[0].decode()
    return np.array([float(x) for x in s.split()]), None
case = sys.argv[1]; t = sys.argv[2] if len(sys.argv) > 2 else "0"
procs = sorted(glob.glob(case + "/processor*"), key=lambda p: int(p.split("processor")[-1]))
rows = []
for p in procs:
    a, u = read_field(p + "/" + t + "/alpha.water")
    if a is None:
        nc = None; ni = 0
    else:
        nc = len(a); ni = int(((a > 1e-8) & (a < 1 - 1e-8)).sum()); nw = int(((a > 1e-8) & (a < 1 - 1e-8) & ((a < 1e-3) | (a > 1 - 1e-3))).sum())
    rows.append((os.path.basename(p), nc, ni))
ni = np.array([r[2] for r in rows]); nc = np.array([r[1] or 0 for r in rows])
print("%s t=%s: ranks %d, cells/rank %d..%d, partial cells (1e-8<a<1-1e-8) per rank: total %d, mean %.1f, max %d, min %d, ranks with 0: %d, max/mean %.2f" % (case.split("/")[-2] + "/" + case.split("/")[-1], t, len(rows), nc.min(), nc.max(), ni.sum(), ni.mean(), ni.max(), ni.min(), int((ni == 0).sum()), ni.max()/max(ni.mean(), 1e-9)))
print("  per rank:", " ".join(str(x) for x in ni))
