#!/usr/bin/env python3
"""Generate a round-channel STL, for meshing a geometry that does NOT align with the mesh.

WHY THIS EXISTS. The four-mesher amplification comparison of 2026-09-09 used a BOX, and a box
is the easy case: it is axis-aligned with both hexahedral meshers, so cartesianMesh cut
nothing (8 cells below 0.8 h, none below 0.6 h) and snappyHexMesh's layers sat on flat walls.
That measured "the hexahedral family is clean WHERE THE GEOMETRY ALIGNS WITH THE MESH", which
is much weaker than "the hexahedral family is clean". A complex microfluidic geometry forces
both meshers to CUT and SNAP cells at a curved, non-aligned wall, and cut cells are small and
one-sided -- the amplifier signature.

A round channel is the cheapest geometry that breaks alignment: one curved wall, no corners, no
junctions. If the amplification bound rises here, no further geometry is needed to answer the
question. If it does not, escalate to corners and junctions.

The cylinder is written OPEN and OVERLONG in x, so it extends past the background mesh at both
ends. The meshed region is then closed laterally by this surface and in x by the background
blockMesh's own inlet and outlet patches, which survive as the cut circles. That is how a pipe
is normally meshed, and it means NO patch renaming: the case's fields declare inlet, outlet and
walls, and all three exist. `constant/polyMesh/boundary` is the authority -- check it after
meshing, because OpenFOAM silently ignores a field entry that matches no mesh patch.

The surface is triangulated finer than the target cell size, or snapping resolves the mesh
rather than the geometry. Default: 240 divisions around a 0.9 mm radius gives 23.6 um facets
against a 39 um cell.

Usage:
  python3 workflow/scripts/make_channel_stl.py out.stl [--radius 9e-4] [--x0 -5e-4]
      [--x1 5.5e-3] [--yc 1.25e-3] [--zc 1.25e-3] [--ntheta 240] [--nx 80] [--name walls]
"""
import argparse
import math


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--radius", type=float, default=9.0e-4)
    ap.add_argument("--x0", type=float, default=-5.0e-4,
                    help="start of the cylinder, BEFORE the domain inlet")
    ap.add_argument("--x1", type=float, default=5.5e-3,
                    help="end of the cylinder, BEYOND the domain outlet")
    ap.add_argument("--yc", type=float, default=1.25e-3)
    ap.add_argument("--zc", type=float, default=1.25e-3)
    ap.add_argument("--ntheta", type=int, default=240)
    ap.add_argument("--nx", type=int, default=80)
    ap.add_argument("--name", default="walls",
                    help="STL solid name; snappyHexMesh names the patch after it")
    args = ap.parse_args()

    R, yc, zc = args.radius, args.yc, args.zc
    xs = [args.x0 + (args.x1 - args.x0) * i / args.nx for i in range(args.nx + 1)]
    th = [2.0 * math.pi * j / args.ntheta for j in range(args.ntheta + 1)]

    def p(i, j):
        return (xs[i], yc + R * math.cos(th[j]), zc + R * math.sin(th[j]))

    def facet(fh, a, b, c):
        # Outward normal of the cylinder side; sign is irrelevant to snappyHexMesh's
        # inside test (locationInMesh decides), but a consistent orientation keeps
        # surfaceCheck quiet.
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        m = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
        fh.write(f"  facet normal {nx/m:.6e} {ny/m:.6e} {nz/m:.6e}\n    outer loop\n")
        for q in (a, b, c):
            fh.write(f"      vertex {q[0]:.9e} {q[1]:.9e} {q[2]:.9e}\n")
        fh.write("    endloop\n  endfacet\n")

    n = 0
    with open(args.out, "w") as fh:
        fh.write(f"solid {args.name}\n")
        for i in range(args.nx):
            for j in range(args.ntheta):
                a, b, c, d = p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)
                facet(fh, a, b, c)
                facet(fh, a, c, d)
                n += 2
        fh.write(f"endsolid {args.name}\n")

    circ = 2.0 * math.pi * R
    print(f"wrote {args.out}: {n} facets, solid '{args.name}'")
    print(f"  cylinder R = {R:.4e} m, axis (y,z) = ({yc:.4e}, {zc:.4e}), "
          f"x from {args.x0:.4e} to {args.x1:.4e}")
    print(f"  facet size around the circumference = {circ/args.ntheta:.4e} m "
          f"(keep it below the target cell size)")


if __name__ == "__main__":
    main()
