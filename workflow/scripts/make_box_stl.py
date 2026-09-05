#!/usr/bin/env python3
"""Write an ASCII STL of an axis-aligned box with NAMED solids for cfMesh.

cfMesh turns every `solid <name>` of the STL into a patch of that name, so the box
comes out with exactly the patches the case's 0.org fields declare: `inlet` (x = 0),
`outlet` (x = xlen) and `walls` (the four sides), outward normals, two triangles per
face. The stationary-droplet box (box6mm-edges.stl) had six VTK-named solids that
meshDict had to rename; here the names are chosen at the source.

    python3 workflow/scripts/make_box_stl.py --xlen 2 --ylen 1 --zlen 1 \\
        --out cases/popinetTranslating3D/box2x1x1.stl
"""
import argparse


def tri(fh, n, a, b, c):
    fh.write(f" facet normal {n[0]:g} {n[1]:g} {n[2]:g}\n  outer loop\n")
    for p in (a, b, c):
        fh.write(f"   vertex {p[0]:.10g} {p[1]:.10g} {p[2]:.10g}\n")
    fh.write("  endloop\n endfacet\n")


def quad(fh, n, p0, p1, p2, p3):
    """Quad p0-p1-p2-p3 counter-clockwise seen from outside (normal n)."""
    tri(fh, n, p0, p1, p2)
    tri(fh, n, p0, p2, p3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlen", type=float, required=True)
    ap.add_argument("--ylen", type=float, required=True)
    ap.add_argument("--zlen", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    X, Y, Z = a.xlen, a.ylen, a.zlen
    P = lambda x, y, z: (x, y, z)
    with open(a.out, "w") as fh:
        fh.write("solid inlet\n")
        quad(fh, (-1, 0, 0), P(0, 0, 0), P(0, 0, Z), P(0, Y, Z), P(0, Y, 0))
        fh.write("endsolid inlet\nsolid outlet\n")
        quad(fh, (1, 0, 0), P(X, 0, 0), P(X, Y, 0), P(X, Y, Z), P(X, 0, Z))
        fh.write("endsolid outlet\nsolid walls\n")
        quad(fh, (0, -1, 0), P(0, 0, 0), P(X, 0, 0), P(X, 0, Z), P(0, 0, Z))   # y = 0
        quad(fh, (0, 1, 0), P(0, Y, 0), P(0, Y, Z), P(X, Y, Z), P(X, Y, 0))    # y = Y
        quad(fh, (0, 0, -1), P(0, 0, 0), P(0, Y, 0), P(X, Y, 0), P(X, 0, 0))   # z = 0
        quad(fh, (0, 0, 1), P(0, 0, Z), P(X, 0, Z), P(X, Y, Z), P(0, Y, Z))    # z = Z
        fh.write("endsolid walls\n")
    print(f"wrote {a.out}: box {X} x {Y} x {Z}, solids inlet/outlet/walls, 12 facets")


if __name__ == "__main__":
    main()
