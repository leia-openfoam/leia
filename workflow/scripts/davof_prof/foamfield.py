"""foamfield.py: OpenFOAM list and field readers without numpy (the compute nodes' python has none);
Claude, 2026-11-04. read_labels(path) -> list of int; read_scalars(path) -> (list of float, uniform value)."""
import re, struct, array

def _fmt(b):
    return re.search(r"format\s+(\w+);", b[:2000].decode(errors="replace")).group(1)

def read_labels(path):
    b = open(path, "rb").read()
    m = re.search(rb"\n(\d+)\s*\n\(", b)
    n = int(m.group(1)); o = m.end()
    if _fmt(b) == "binary":
        a = array.array("i"); a.frombytes(b[o:o + 4*n]); return a.tolist()
    return [int(x) for x in b[o:].split(b")", 1)[0].split()]

def read_scalars(path):
    b = open(path, "rb").read()
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", b)
    if not m:
        m2 = re.search(rb"internalField\s+uniform\s+([^;]+);", b)
        return None, float(m2.group(1))
    n = int(m.group(1)); o = m.end()
    if _fmt(b) == "binary":
        a = array.array("d"); a.frombytes(b[o:o + 8*n]); return a.tolist(), None
    return [float(x) for x in b[o:].split(b")", 1)[0].split()], None

def band_weights(own, nei, a, W, layers):
    """the cells within <layers> face layers of a partial cell (1e-8 < a < 1 - 1e-8) get W, the rest 1"""
    nc = len(a)
    mark = [1e-8 < x < 1 - 1e-8 for x in a]
    nb = sum(mark)
    ni = len(nei)
    for _ in range(layers):
        new = mark[:]
        for f in range(ni):
            o, n = own[f], nei[f]
            if mark[o] or mark[n]:
                new[o] = True; new[n] = True
        mark = new
    return [W if m else 1.0 for m in mark], nb, sum(mark)

def write_weights(path, t, w):
    with open(path, "w") as f:
        f.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       volScalarField;\n    location    \"%s\";\n    object      cellWeight;\n}\ndimensions      [0 0 0 0 0 0 0];\ninternalField   nonuniform List<scalar>\n%d\n(\n" % (t, len(w)))
        f.write("\n".join("%g" % x for x in w))
        f.write("\n)\n;\nboundaryField\n{\n    \".*\"\n    {\n        type            calculated;\n        value           uniform 1;\n    }\n}\n")
