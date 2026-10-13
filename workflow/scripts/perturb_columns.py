"""Perturb an OpenFOAM blockMesh hex mesh column-wise: every point gets a
displacement (dx, dy) that depends on its (x, y) only, so all points above
each other move together; z is not changed and points on the x/y boundary
planes stay fixed. Then every face stays PLANAR (lateral faces are vertical
quads with two vertical edges, top/bottom faces lie in z = const planes),
which separates the effect of warped faces from that of an irregular mesh.

usage: python3 perturb_columns.py <case> <amplitude in cell sizes> [seed]
"""
import random
import re
import sys

case, amp = sys.argv[1], float(sys.argv[2])
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
path = case + "/constant/polyMesh/points"
txt = open(path).read()
head, rest = txt.split("(\n", 1)
body, tail = rest.rsplit(")\n", 1)
pts = []
for line in body.strip().splitlines():
    m = re.match(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", line.strip())
    pts.append(tuple(float(v) for v in m.groups()))
xs = sorted({round(p[0], 12) for p in pts})
ys = sorted({round(p[1], 12) for p in pts})
h = xs[1] - xs[0]
rng = random.Random(seed)
disp = {}
out = []
for x, y, z in pts:
    kx, ky = round(x, 12), round(y, 12)
    if kx in (xs[0], xs[-1]) or ky in (ys[0], ys[-1]):
        out.append((x, y, z))
        continue
    if (kx, ky) not in disp:
        disp[(kx, ky)] = (rng.uniform(-amp, amp)*h, rng.uniform(-amp, amp)*h)
    dx, dy = disp[(kx, ky)]
    out.append((x + dx, y + dy, z))
with open(path, "w") as fh:
    fh.write(head + "(\n")
    for p in out:
        fh.write("(%.15g %.15g %.15g)\n" % p)
    fh.write(")\n" + tail)
print("perturbed %d of %d points column-wise, amplitude %.2f h, h = %g" % (len(disp), len(pts), amp, h))
