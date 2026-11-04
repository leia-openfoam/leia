"""cellweight.py <case> <time> <W> [<layers>]: the cell weights of an interface-weighted
decomposition (Claude, 2026-11-04): a cell within <layers> face-neighbour layers of a partial
cell (1e-8 < alpha < 1 - 1e-8) at <time> gets the weight W, every other cell 1; written as the
volScalarField cellWeight in <case>/<time>/ (ascii; the mesh read from constant/polyMesh, ascii
or binary owner/neighbour). Serial (undecomposed) cases."""
import sys, re, os
import numpy as np
case, t, W = sys.argv[1], sys.argv[2], float(sys.argv[3])
layers = int(sys.argv[4]) if len(sys.argv) > 4 else 2

def read_list(path, dtype):
    b = open(path, "rb").read()
    h = b[:1500].decode(errors="replace")
    fmt = re.search(r"format\s+(\w+);", h).group(1)
    m = re.search(rb"\n(\d+)\s*\n\(", b)
    n = int(m.group(1)); o = m.end()
    if fmt == "binary":
        return np.frombuffer(b, dtype=dtype, count=n, offset=o)
    return np.array(b[o:].split(b")", 1)[0].split(), dtype=dtype)

def read_field(path):
    b = open(path, "rb").read()
    h = b[:1500].decode(errors="replace")
    fmt = re.search(r"format\s+(\w+);", h).group(1)
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", b)
    n = int(m.group(1)); o = m.end()
    if fmt == "binary":
        return np.frombuffer(b, dtype="<f8", count=n, offset=o)
    return np.array(b[o:].split(b")", 1)[0].split(), dtype="<f8")

own = read_list(os.path.join(case, "constant/polyMesh/owner"), "<i4")
nei = read_list(os.path.join(case, "constant/polyMesh/neighbour"), "<i4")
a = read_field(os.path.join(case, t, "alpha.water"))
nc = len(a)
band = (a > 1e-8) & (a < 1 - 1e-8)
mark = band.copy()
ni = len(nei)
for _ in range(layers):
    new = mark.copy()
    o, n = own[:ni], nei
    sel = mark[o] | mark[n]
    new[o[sel]] = True; new[n[sel]] = True
    mark = new
w = np.where(mark, W, 1.0)
hdr = """FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "%s";
    object      cellWeight;
}
dimensions      [0 0 0 0 0 0 0];
internalField   nonuniform List<scalar>
%d
(
""" % (t, nc)
with open(os.path.join(case, t, "cellWeight"), "w") as f:
    f.write(hdr)
    f.write("\n".join("%g" % x for x in w))
    f.write("\n)\n;\nboundaryField\n{\n    \".*\"\n    {\n        type            calculated;\n        value           uniform 1;\n    }\n}\n")
print("cellWeight: %d cells, %d partial, %d weighted (W = %g, %d layers), total weight %g, interface share of the weight %.3f" % (nc, int(band.sum()), int(mark.sum()), W, layers, w.sum(), (w.sum() - (nc - mark.sum()))/w.sum()))
