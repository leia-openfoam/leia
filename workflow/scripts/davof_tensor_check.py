#!/usr/bin/env python3
"""Check of the ring curvature tensor and vertex offset the sagitta correction uses
(study sphereTranslation3D_hexSagittaDiag, 2026-10-16): reads the ascii fields
curvatureTensorDavof.<group>, curvatureOffsetDavof.<group> and interfaceNormal.<group> of
one time directory of a case and prints, over the cells with a non-zero tensor, the
eigenvalues of B against 1/R (two at 1/R, one at zero on a sphere), the sign and size of
z_0 against -h^2/(12 R), and the alignment of B's null vector with the normal.

    python3 davof_tensor_check.py <case> <time> [--radius 0.15] [--group water]
"""
import argparse
import math
import os
import re
import sys


def read_list(path):
    """The internalField of an ascii volField as a list of tuples of floats."""
    with open(path, encoding="utf-8", errors="replace") as f:
        txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*\n?(\d+)\s*\n?\(", txt)
    if not m:
        m2 = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
        if m2:
            sys.exit(f"{path}: uniform internalField {m2.group(1).strip()} (nothing written per cell)")
        sys.exit(f"{path}: no internalField list found")
    n = int(m.group(2))
    body = txt[m.end():]
    end = body.index("\n)\n") if "\n)\n" in body else body.index(")\n;")
    body = body[:end]
    if m.group(1) == "scalar":
        vals = [(float(v),) for v in body.split()]
    else:
        vals = [tuple(float(v) for v in t.split()) for t in re.findall(r"\(([^()]*)\)", body)]
    if len(vals) != n:
        sys.exit(f"{path}: {len(vals)} entries read, {n} declared")
    return vals


def sym_eig(s):
    """Eigenvalues of a symmetric 3x3 (xx xy xz yy yz zz), by the closed form."""
    xx, xy, xz, yy, yz, zz = s
    p1 = xy*xy + xz*xz + yz*yz
    q = (xx + yy + zz)/3.0
    p2 = (xx - q)**2 + (yy - q)**2 + (zz - q)**2 + 2*p1
    p = math.sqrt(p2/6.0)
    if p < 1e-300:
        return [xx, yy, zz]
    b = [(xx - q)/p, xy/p, xz/p, (yy - q)/p, yz/p, (zz - q)/p]
    bxx, bxy, bxz, byy, byz, bzz = b
    r = (bxx*(byy*bzz - byz*byz) - bxy*(bxy*bzz - byz*bxz) + bxz*(bxy*byz - byy*bxz))/2.0
    r = max(-1.0, min(1.0, r))
    phi = math.acos(r)/3.0
    e1 = q + 2*p*math.cos(phi)
    e3 = q + 2*p*math.cos(phi + 2*math.pi/3)
    e2 = 3*q - e1 - e3
    return sorted([e1, e2, e3])


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("time")
    ap.add_argument("--radius", type=float, default=0.15)
    ap.add_argument("--group", default="water")
    ap.add_argument("--h", type=float, default=1.0/32)
    a = ap.parse_args(argv)
    td = os.path.join(a.case, a.time)
    B = read_list(os.path.join(td, f"curvatureTensorDavof.{a.group}"))
    z0 = read_list(os.path.join(td, f"curvatureOffsetDavof.{a.group}"))
    nrm = read_list(os.path.join(td, f"interfaceNormal.{a.group}"))
    kappa = 1.0/a.radius
    cells = [i for i, s in enumerate(B) if any(abs(v) > 0 for v in s)]
    print(f"cells with a tensor: {len(cells)} of {len(B)}; 1/R = {kappa:.4f} 1/m; "
          f"-h^2/(12 R) = {-a.h*a.h/(12*a.radius):.3e} m")
    n_ok = n_neg = n_null = 0
    ev_mid = []
    z0s = []
    rel = []
    for i in cells:
        ev = sym_eig(B[i])             # ascending
        lam_small = min(abs(e) for e in ev)
        others = sorted([e for e in ev if abs(e) != lam_small] or ev[1:])
        if len(others) < 2:
            others = ev[1:]
        r1, r2 = others[0]/kappa, others[1]/kappa
        rel.append((r1, r2))
        if abs(r1 - 1) < 0.25 and abs(r2 - 1) < 0.25:
            n_ok += 1
        if z0[i][0] < 0:
            n_neg += 1
        z0s.append(z0[i][0])
        ev_mid.append(0.5*(others[0] + others[1]))
        # alignment: B n should vanish on a sphere (n the surface normal)
        nx, ny, nz = nrm[i]
        nn = math.sqrt(nx*nx + ny*ny + nz*nz)
        if nn > 0:
            xx, xy, xz, yy, yz, zz = B[i]
            bn = ((xx*nx + xy*ny + xz*nz)/nn, (xy*nx + yy*ny + yz*nz)/nn, (xz*nx + yz*ny + zz*nz)/nn)
            if math.sqrt(sum(v*v for v in bn)) < 0.25*kappa:
                n_null += 1
    nc = max(1, len(cells))
    z0s.sort()
    ev_mid.sort()
    rel1 = sorted(r[0] for r in rel)
    rel2 = sorted(r[1] for r in rel)
    med = lambda L: L[len(L)//2] if L else float("nan")
    print(f"eigenvalues within 25 % of 1/R (both tangential): {n_ok} ({100*n_ok/nc:.1f} %)")
    print(f"B n small (|B n| < 0.25/R): {n_null} ({100*n_null/nc:.1f} %)")
    print(f"z0 < 0: {n_neg} ({100*n_neg/nc:.1f} %); z0 median {med(z0s):.3e}, "
          f"10 % {z0s[len(z0s)//10]:.3e}, 90 % {z0s[9*len(z0s)//10]:.3e}")
    print(f"tangential eigenvalues / (1/R): medians {med(rel1):.3f} {med(rel2):.3f}; "
          f"10 % {rel1[len(rel1)//10]:.3f} {rel2[len(rel2)//10]:.3f}; "
          f"90 % {rel1[9*len(rel1)//10]:.3f} {rel2[9*len(rel2)//10]:.3f}")
    print(f"mean curvature (B11 + B22)/2 relative to 1/R: median {med(ev_mid)/kappa:.3f}")
    print("first ten cells: eigenvalues/kappa, z0")
    for i in cells[:10]:
        ev = sym_eig(B[i])
        print("  ", i, [round(e/kappa, 3) for e in ev], f"{z0[i][0]:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
