"""cellweight_par.py <case> <time> <W> [<layers>]: the cell weights (cellweight.py) of a DECOMPOSED
case, written per processor at <time> from each processor's alpha.water and polyMesh (the
face-neighbour layers are grown inside each processor; a band cell's layer across a processor
boundary is missed, at most LAYERS cells deep). Claude, 2026-11-04."""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
case, t, W = sys.argv[1], sys.argv[2], float(sys.argv[3])
layers = int(sys.argv[4]) if len(sys.argv) > 4 else 2

def read_list(path, dtype):
    b = open(path, "rb").read()
    fmt = re.search(r"format\s+(\w+);", b[:1500].decode(errors="replace")).group(1)
    m = re.search(rb"\n(\d+)\s*\n\(", b)
    n = int(m.group(1)); o = m.end()
    if fmt == "binary": return np.frombuffer(b, dtype=dtype, count=n, offset=o)
    return np.array(b[o:].split(b")", 1)[0].split(), dtype=dtype)

def read_field(path):
    b = open(path, "rb").read()
    fmt = re.search(r"format\s+(\w+);", b[:1500].decode(errors="replace")).group(1)
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", b)
    if not m:
        m2 = re.search(rb"internalField\s+uniform\s+([^;]+);", b)
        return None, float(m2.group(1))
    n = int(m.group(1)); o = m.end()
    if fmt == "binary": return np.frombuffer(b, dtype="<f8", count=n, offset=o), None
    return np.array(b[o:].split(b")", 1)[0].split(), dtype="<f8"), None

tot = 0; totb = 0; totw = 0; per = []
for p in sorted(glob.glob(case + "/processor*"), key=lambda p: int(p.split("processor")[-1])):
    own = read_list(p + "/constant/polyMesh/owner", "<i4"); nei = read_list(p + "/constant/polyMesh/neighbour", "<i4")
    nc = int(own.max()) + 1 if len(own) else 0
    a, u = read_field(p + "/" + t + "/alpha.water")
    if a is None: a = np.full(nc, u)
    band = (a > 1e-8) & (a < 1 - 1e-8)
    mark = band.copy(); ni = len(nei)
    for _ in range(layers):
        new = mark.copy(); o, n = own[:ni], nei
        sel = mark[o] | mark[n]; new[o[sel]] = True; new[n[sel]] = True; mark = new
    w = np.where(mark, W, 1.0)
    with open(p + "/" + t + "/cellWeight", "w") as f:
        f.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       volScalarField;\n    location    \"%s\";\n    object      cellWeight;\n}\ndimensions      [0 0 0 0 0 0 0];\ninternalField   nonuniform List<scalar>\n%d\n(\n" % (t, len(w)))
        f.write("\n".join("%g" % x for x in w))
        f.write("\n)\n;\nboundaryField\n{\n    \".*\"\n    {\n        type            calculated;\n        value           uniform 1;\n    }\n}\n")
    tot += len(w); totb += int(band.sum()); totw += int(mark.sum()); per.append(int(band.sum()))
print("cellWeight per processor at t=%s: %d cells, %d partial, %d weighted (W %g, %d layers); partial per rank max %d mean %.1f" % (t, tot, totb, totw, W, layers, max(per), sum(per)/len(per)))
