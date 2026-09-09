#!/usr/bin/env python3
"""Where does the quasi-monotone clip fire? Classify the cells of slClipFired.

The clip's counter says how often the fit created a new extremum. It does not say
WHERE, and where is what decides the region rule: a firing at a genuine extremum of
the level set needs an exemption, a firing in a small one-sided cell needs the bound,
and a firing at the interface is a method change. This script reads slClipFired (and
slClipEligible, and psi) from every write time of a probe arm and classifies each
firing.

The classification needs the reconstruction's cell-point-cell stencil, which exists
only inside OpenFOAM. On a UNIFORM STRUCTURED HEX mesh built by blockMesh the stencil
is exactly the surrounding cells in index space, and this script rebuilds it from
--nx/--ny/--nz: blockMesh numbers cells with i fastest, then j, then k. That covers
the hexahedral probe. For a polyhedral mesh use the solver's own slStencilPosition
field instead, which carries (psi_c - lo)/(hi - lo) per cell and needs no connectivity.

Per fired cell it reports:
  - the index position (i, j, k) and the cell centre;
  - BOUNDARY or interior (a cell on a mesh face of the block);
  - whether psi_c is the stencil MINIMUM, MAXIMUM or neither;
  - the stencil position (psi_c - lo)/(hi - lo), so a NEAR-extremum (0.02, say) is
    visible as such -- the exact test psi_c == lo cannot catch a cell the fit has
    already pushed just off its extremum;
  - the distance from the interface in cell widths, |psi_c|/h, which separates a
    far-field firing from one next to the zero set.

Usage:
  python3 workflow/scripts/locate_clip_firings.py <arm-dir> --nx 128 --ny 64 [--nz 1]
      [--lx 0.005 --ly 0.0025 --lz 0]  [--max-list 40]
"""
import argparse
import os
import re
import sys


def read_internal(path, ncells):
    """The internalField of an OpenFOAM volScalarField, uniform or nonuniform."""
    with open(path) as fh:
        s = fh.read()
    m = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\((.*?)\n\)\s*;",
        s, re.S)
    if m:
        vals = [float(v) for v in m.group(2).split()]
        if len(vals) != ncells:
            raise ValueError(f"{path}: {len(vals)} values, expected {ncells}")
        return vals
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", s)
    if m:
        return [float(m.group(1))] * ncells
    raise ValueError(f"{path}: no internalField found")


def time_dirs(arm):
    out = []
    for d in os.listdir(arm):
        if d == "0.org":
            continue
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return [d for _, d in sorted(out)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", help="study arm directory holding the time directories")
    ap.add_argument("--nx", type=int, required=True)
    ap.add_argument("--ny", type=int, required=True)
    ap.add_argument("--nz", type=int, default=1)
    ap.add_argument("--lx", type=float, default=None)
    ap.add_argument("--ly", type=float, default=None)
    ap.add_argument("--lz", type=float, default=None)
    ap.add_argument("--max-list", type=int, default=40,
                    help="print at most this many fired cells per time")
    args = ap.parse_args()

    NX, NY, NZ = args.nx, args.ny, args.nz
    n = NX * NY * NZ
    hx = args.lx / NX if args.lx else None

    def stencil(c):
        """cell-point-cell neighbours in blockMesh index order (i fastest)."""
        i = c % NX
        j = (c // NX) % NY
        k = c // (NX * NY)
        out = []
        for dk in ((-1, 0, 1) if NZ > 1 else (0,)):
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    if di == dj == dk == 0:
                        continue
                    ii, jj, kk = i + di, j + dj, k + dk
                    if 0 <= ii < NX and 0 <= jj < NY and 0 <= kk < NZ:
                        out.append((kk * NY + jj) * NX + ii)
        return out

    def on_boundary(c):
        i = c % NX
        j = (c // NX) % NY
        k = c // (NX * NY)
        return (i in (0, NX - 1) or j in (0, NY - 1)
                or (NZ > 1 and k in (0, NZ - 1)))

    times = time_dirs(args.arm)
    if not times:
        sys.exit(f"no time directories under {args.arm}")

    print(f"clip firing locations: {args.arm}")
    print(f"  mesh {NX} x {NY} x {NZ} = {n} cells, {len(times)} write times")
    print()

    seen = {}
    for t in times:
        fp = os.path.join(args.arm, t, "slClipFired")
        if not os.path.isfile(fp):
            continue
        fired = read_internal(fp, n)
        psi = read_internal(os.path.join(args.arm, t, "psi"), n)
        ep = os.path.join(args.arm, t, "slClipEligible")
        elig = read_internal(ep, n) if os.path.isfile(ep) else [1.0] * n
        idx = [c for c, v in enumerate(fired) if v > 0.5]

        nb_ext = nmin = nmax = nnear = nbnd = 0
        rows = []
        for c in idx:
            nb = [psi[q] for q in stencil(c)]
            lo, hi = min(nb + [psi[c]]), max(nb + [psi[c]])
            pos = (psi[c] - lo) / (hi - lo) if hi > lo else 0.0
            is_min, is_max = psi[c] <= min(nb), psi[c] >= max(nb)
            if is_min:
                nmin += 1
            if is_max:
                nmax += 1
            if is_min or is_max:
                nb_ext += 1
            elif pos < 0.05 or pos > 0.95:
                nnear += 1
            if on_boundary(c):
                nbnd += 1
            rows.append((c, psi[c], pos, is_min, is_max, on_boundary(c)))
            seen.setdefault(c, []).append(t)

        print(f"  t = {t}:  {len(idx)} fired  ({int(sum(elig))} eligible)"
              f"  -- {nb_ext} exact extrema ({nmin} min, {nmax} max),"
              f" {nnear} near-extrema, {nbnd} boundary cells")
        for c, p, pos, is_min, is_max, bnd in rows[:args.max_list]:
            i = c % NX
            j = (c // NX) % NY
            k = c // (NX * NY)
            kind = ("MIN" if is_min else "MAX" if is_max else
                    "near-min" if pos < 0.05 else
                    "near-max" if pos > 0.95 else "interior of range")
            hcell = f" |psi|/h={abs(p)/hx:6.2f}" if hx else ""
            print(f"     cell {c:>6} (i,j,k)=({i:>3},{j:>3},{k:>3}) psi={p:+.4e}"
                  f" pos={pos:5.3f} {kind:<17}"
                  f" {'BOUNDARY' if bnd else 'interior':<8}{hcell}")
        if len(rows) > args.max_list:
            print(f"     ... {len(rows) - args.max_list} more not listed")
        print()

    print(f"  distinct cells that ever fired: {len(seen)}")
    persistent = [c for c, ts in seen.items() if len(ts) == len(times)]
    print(f"  cells that fired at EVERY write time: {len(persistent)}"
          f" -- a fixed set points to the mesh, a changing set to the field")


if __name__ == "__main__":
    main()
