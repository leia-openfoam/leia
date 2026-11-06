#!/usr/bin/env python3
"""The alpha = 0.5 surface of a 3D deformation run (2026-11-06) without ParaView
(the cluster's ParaView modules carry no pvpython): the reconstructed
alpha.water of the N^3 blockMesh cube is reshaped (blockMesh orders the cells
x fastest, then y, then z), the 0.5 surface is extracted with the marching
cubes of scikit-image (W. E. Lorensen, H. E. Cline, Marching cubes: a high
resolution 3D surface construction algorithm, ACM SIGGRAPH Computer Graphics
21 (1987) 163-169, doi:10.1145/37402.37422; the implementation of T. Lewiner,
H. Lopes, A. W. Vieira, G. Tavares, Efficient implementation of marching
cubes' cases with topological guarantees, J. Graphics Tools 8 (2003) 1-15,
doi:10.1080/10867651.2003.10487582) and drawn with matplotlib.
    python3 render_alpha_mc.py <case dir> <time> <N> <out.png> [<title>]
The case: cases/davof/sphereDeformation3D (LeVeque 1996, doi:10.1137/0733033;
Enright et al. 2002, doi:10.1006/jcph.2002.7166).
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

case, t, N, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
title = sys.argv[5] if len(sys.argv) > 5 else ""
fn = os.path.join(case, t, "alpha.water")
txt = open(fn).read()
m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
if not m:
    sys.exit("no nonuniform internalField in " + fn)
n = int(m.group(1)); start = m.end(); end = txt.index("\n)", start)
a = np.array(txt[start:end].split(), dtype=float)
assert a.size == n == N**3, (a.size, n, N**3)
A = a.reshape((N, N, N)).transpose(2, 1, 0)   # [ix, iy, iz]: x fastest in the file
h = 1.0/N
# pad with zeros (the box is dry outside) so that the surface closes at the walls
Ap = np.pad(A, 1, constant_values=0.0)
if Ap.max() < 0.5 or Ap.min() > 0.5:
    sys.exit("no 0.5 level at %s (min %g max %g)" % (t, Ap.min(), Ap.max()))
verts, faces, normals, _ = measure.marching_cubes(Ap, level=0.5, spacing=(h, h, h))
verts = verts - 0.5*h   # cell centres at (i + 0.5) h after the pad's shift of one cell
fig = plt.figure(figsize=(7, 6.2))
ax = fig.add_subplot(111, projection="3d")
mesh = Poly3DCollection(verts[faces], alpha=1.0, linewidths=0.0)
mesh.set_facecolor((0.30, 0.50, 0.85))
mesh.set_edgecolor("none")
ax.add_collection3d(mesh)
# the unit cube's edges
for s, e in [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1)), ((1, 1, 1), (0, 1, 1)), ((1, 1, 1), (1, 0, 1)), ((1, 1, 1), (1, 1, 0)),
             ((1, 0, 0), (1, 1, 0)), ((1, 0, 0), (1, 0, 1)), ((0, 1, 0), (1, 1, 0)), ((0, 1, 0), (0, 1, 1)), ((0, 0, 1), (1, 0, 1)), ((0, 0, 1), (0, 1, 1))]:
    ax.plot(*zip(s, e), color="0.4", linewidth=0.6)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
ax.set_box_aspect((1, 1, 1))
ax.view_init(elev=22, azim=-60)
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
if title: ax.set_title(title, fontsize=10)
fig.text(0.02, 0.02, "alpha = 0.5 (marching cubes), %d triangles, %d cells with 1e-8 < alpha < 1 - 1e-8" % (len(faces), int(((a > 1e-8) & (a < 1 - 1e-8)).sum())), fontsize=7)
fig.tight_layout()
fig.savefig(out, dpi=150)
print("wrote", out, "triangles", len(faces))
