#!/usr/bin/env python3
"""compare_psi_fields.py: WHERE do two runs of the same case differ? (PDE level set, step 0.2)

    compare_psi_fields.py <caseA> <caseB> <time> [--field psi] [--celldist <file>] [--n N]

Reads the internalField of <field> at <time> in two case directories (ascii, a serial case or a
reconstructed parallel case), and prints:
  * the L1, L2 and maximum of the cell-wise difference, absolute and relative to max|field|;
  * the twelve largest differences with the cell index, its (i, j) on an N x N uniform 2D hex
    mesh (blockMesh order, x fastest), and its processor if a cell-to-processor list is given;
  * if a cell-to-processor list is given: the share of the total |difference| that sits in cells
    with a face neighbour on another processor (the seam cells), against the seam cells' share of
    the mesh. A seam defect concentrates the difference there; solver round-off spreads it.
The cell-to-processor list is the `cellDist` field (decomposePar -cellDist) or the
constant/cellDecomposition labelList. Stdlib only. Exit 0 always: this is a diagnostic, not a gate.
"""
import argparse, math, os, re, sys


def read_internal(path):
    txt = open(path, errors="ignore").read()
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", txt)
    if m:
        return None, float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<(?:scalar|label)>\s*(\d+)\s*\(", txt)
    if not m:
        m = re.search(r"^(\d+)\s*\n\(", txt, re.M)          # a bare labelList (cellDecomposition)
        if not m:
            raise SystemExit(f"cannot parse {path}")
    n = int(m.group(1))
    body = txt[m.end():]
    vals = body.split(")", 1)[0].split()
    if len(vals) != n:
        raise SystemExit(f"{path}: expected {n} values, read {len(vals)}")
    return [float(v) for v in vals], None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("time")
    ap.add_argument("--field", default="psi")
    ap.add_argument("--celldist", default=None)
    ap.add_argument("--n", type=int, default=None, help="cells per side of a uniform 2D hex mesh")
    a = ap.parse_args()
    A, ua = read_internal(os.path.join(a.a, a.time, a.field))
    B, ub = read_internal(os.path.join(a.b, a.time, a.field))
    if A is None or B is None:
        print("uniform field in one case; nothing to compare"); return 0
    if len(A) != len(B):
        raise SystemExit(f"cell counts differ: {len(A)} against {len(B)}")
    nc = len(A)
    N = a.n or int(round(math.sqrt(nc)))
    if N * N != nc:
        N = None
    d = [x - y for x, y in zip(A, B)]
    scale = max(abs(x) for x in A) or 1.0
    l1 = sum(abs(x) for x in d) / nc
    l2 = math.sqrt(sum(x * x for x in d) / nc)
    mx = max(abs(x) for x in d)
    print(f"[compare_psi_fields] {a.field} at t = {a.time}: {nc} cells; max|{a.field}| = {scale:.6g}")
    print(f"  L1 {l1:.3e}  L2 {l2:.3e}  max {mx:.3e}  (relative to max|field|: {l1/scale:.2e} / {l2/scale:.2e} / {mx/scale:.2e})")
    print(f"  cells with |diff| > 1e-14 * scale: {sum(1 for x in d if abs(x) > 1e-14*scale)}; "
          f"> 1e-10: {sum(1 for x in d if abs(x) > 1e-10*scale)}; > 1e-6: {sum(1 for x in d if abs(x) > 1e-6*scale)}")
    proc = None
    if a.celldist:
        p, u = read_internal(a.celldist)
        proc = [int(round(x)) for x in p] if p is not None else None
    order = sorted(range(nc), key=lambda c: -abs(d[c]))[:12]
    for c in order:
        ij = f"(i {c % N:4d}, j {c // N:4d})" if N else ""
        pr = f" proc {proc[c]}" if proc else ""
        print(f"    cell {c:7d} {ij}{pr}  diff {d[c]: .3e}  {a.field} {A[c]: .4e}")
    if proc and N:
        seam = set()
        for c in range(nc):
            i, j = c % N, c // N
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < N and 0 <= jj < N and proc[jj * N + ii] != proc[c]:
                    seam.add(c); break
        tot = sum(abs(x) for x in d) or 1.0
        share = sum(abs(d[c]) for c in seam) / tot
        print(f"  seam cells: {len(seam)} of {nc} ({len(seam)/nc:.1%} of the mesh) carry {share:.1%} of the total |diff|")
    return 0


if __name__ == "__main__":
    sys.exit(main())
