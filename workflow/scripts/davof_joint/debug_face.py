#!/usr/bin/env python3
"""Per-tetrahedron trace of the evaluator on one face of a dump (the defect hunt)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from joint_prototype import (read_dump, triangle_list, sub, add, mul, dot, cross, norm, newell,
                             tet_faces, clip_polyhedron, polyhedron_volume, clip_polygon)


def closed_form(T, n, x0):
    total = (0.0, 0.0, 0.0); moment = 0.0
    for (a, b, c) in T:
        poly, _ = clip_polygon([a, b, c], n, x0)
        if not poly:
            continue
        aw = newell(poly)
        total = add(total, aw); moment += dot(a, aw)
    cap = -dot(n, total)
    return (moment + dot(x0, n)*cap)/3.0


def main():
    path, facei = sys.argv[1], int(sys.argv[2])
    hdr, faces = read_dump(path)
    f = [g for g in faces if g.facei == facei][0]
    T = triangle_list(f)
    print("face", f.facei, "vTotal", f.vTotal, "vExact", f.vExact, "vPieces", f.vPieces, "preSign", f.preSign)
    planes = [(c.n, c.x0) for c in f.cells if c.wetPlane and c.cov > 0]
    for pid, (n, x0) in enumerate(planes):
        print("plane", pid, n, x0, "closed form", closed_form(T, n, x0))
        tot = 0.0
        for k, (a, b, c) in enumerate(T):
            tp = dot(cross(sub(a, f.xf), sub(b, f.xf)), sub(c, f.xf))
            if abs(tp) <= 1e-300:
                print("  tet", k, "degenerate"); continue
            sign = 1.0 if tp > 0 else -1.0
            faces_ = tet_faces(f.xf, a, b, c) if tp > 0 else tet_faces(f.xf, a, c, b)
            v0 = polyhedron_volume(faces_)
            cl = clip_polyhedron(faces_, n, x0, pid)
            v1 = polyhedron_volume(cl) if cl else 0.0
            s = [dot(sub(p, x0), n) for p in (f.xf, a, b, c)]
            print("  tet %2d sign %+d vol %.6e clipped %.6e  s = %s  faces %d" % (k, int(sign), v0, v1, ["%.2e" % x for x in s], len(cl)))
            if cl:
                for poly, lab in cl:
                    av = newell(poly)
                    print("      face lab %2d n %d area %.3e  pts %s" % (lab, len(poly), norm(av), ["(%.6f,%.6f,%.6f)" % p for p in poly]))
            tot += sign*v1
        print("  sum of signed clipped volumes", tot)


if __name__ == "__main__":
    main()
