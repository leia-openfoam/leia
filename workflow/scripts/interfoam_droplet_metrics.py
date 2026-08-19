#!/usr/bin/env python3
"""Convert interFoam function-object output into the droplet-metrics CSV.

WHY THIS EXISTS. interFoam writes no metrics CSV of its own, and the comparison
against the level-set solver is only worth anything if BOTH sides are scored by
identically defined quantities. So this script reads the function-object output
that the interFoam case's controlDict writes and emits `interFoam.csv` with the
same column names the level-set solver writes, letting one analysis path serve
both:

  TIME                  time
  maxMagU               max|U| over the domain           <- fieldMinMax(U)
  phaseVolume           integral of alpha                <- volFieldValue
  phaseVolumeRelError   |V(t) - V(0)|/V(0), a CONSERVATION error relative to the
                        INITIAL DISCRETE volume -- the same definition the
                        level-set solver uses, and domain-size independent
  pLaplace              p_rgh(centre) - p_rgh(far field) <- probes
  zeroSetRadialL2       RMS of |x - centre| - R over the alpha = 0.5 isosurface
                        points, at the write instants only  <- surfaces(raw)

TWO HONEST CAVEATS, because they decide what may be compared.

1. max|U| and phaseVolumeRelError ARE directly comparable between the solvers:
   both are max|U| over the same mesh, and both volume errors are relative to
   each run's own initial discrete volume.

2. zeroSetRadialL2 is NOT directly comparable to the level-set solver's column of
   the same name. Here it is the radial spread of the alpha = 0.5 isosurface;
   there it is computed from the psi zero set. Each solver's shape TREND under
   refinement is meaningful, the absolute values are not interchangeable, and the
   column is named identically only so one analysis path can read both files.
   Getting a common shape metric needs the same isosurface function object added
   to the level-set cases, which would mean re-running them.

Usage: interfoam_droplet_metrics.py [--case .] [--out interFoam.csv]
                                    [--radius 1e-3] [--centre x,y,z]
"""
import argparse
import csv
import glob
import math
import os
import re
import sys


def _read_dat(path, ncols_min=2):
    """OpenFOAM function-object .dat: '#' comments, whitespace columns."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            if len(parts) < ncols_min:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return rows


def _latest_dir(pattern):
    """Function objects restart into a new time-stamped dir; take them all and
    let the caller sort by time so a restarted run is not silently truncated."""
    return sorted(glob.glob(pattern))


def collect(case, name, fname, ncols_min=2):
    out = []
    for d in _latest_dir(os.path.join(case, "postProcessing", name, "*")):
        f = os.path.join(d, fname)
        if os.path.isfile(f):
            out += _read_dat(f, ncols_min)
    # de-duplicate on time, keeping the last occurrence (a restart overwrites)
    seen = {}
    for r in out:
        seen[round(r[0], 12)] = r
    return [seen[k] for k in sorted(seen)]


def find_file(case, name, stem):
    for d in _latest_dir(os.path.join(case, "postProcessing", name, "*")):
        for f in sorted(glob.glob(os.path.join(d, stem))):
            return os.path.basename(f)
    return None


def shape_error(case, centre, radius):
    """RMS radial error over the raw alpha=0.5 isosurface points, per write time."""
    out = {}
    for d in _latest_dir(os.path.join(case, "postProcessing", "interface", "*")):
        for raw in sorted(glob.glob(os.path.join(d, "*.raw"))):
            t = os.path.basename(os.path.dirname(raw))
            try:
                tv = float(t)
            except ValueError:
                continue
            s2 = n = 0
            with open(raw) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    p = line.split()
                    if len(p) < 3:
                        continue
                    try:
                        x, y, z = float(p[0]), float(p[1]), float(p[2])
                    except ValueError:
                        continue
                    r = math.sqrt((x-centre[0])**2 + (y-centre[1])**2
                                  + (z-centre[2])**2)
                    s2 += (r - radius)**2
                    n += 1
            if n:
                out[round(tv, 12)] = math.sqrt(s2/n)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--out", default="interFoam.csv")
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--centre", default=None,
                    help="x,y,z of the droplet centre; default = mesh box centre "
                         "read from the setAlphaFieldDict origin")
    args = ap.parse_args()
    case = args.case

    centre = None
    if args.centre:
        centre = [float(v) for v in args.centre.split(",")]
    else:
        # Read it from the dict the case was actually initialised with: the
        # interFoam arm shares leiaSetFields with the level-set cases, so the
        # droplet centre lives in levelSet.implicitSurface.center of fvSolution.
        # setAlphaFieldDict is also accepted, for a case initialised that way.
        for rel, key in (("system/fvSolution", "center"),
                         ("system/setAlphaFieldDict", "origin")):
            d = os.path.join(case, rel)
            if not os.path.isfile(d):
                continue
            for line in open(d):
                t = line.strip()
                if t.startswith("//") or not t.startswith(key):
                    continue
                v = t.split(key, 1)[1]
                v = v.replace("(", " ").replace(")", " ").replace(";", " ")
                try:
                    nums = [float(x) for x in v.split()]
                except ValueError:
                    continue
                if len(nums) == 3:
                    centre = nums
                    break
            if centre is not None:
                break
    if centre is None:
        sys.exit("[interfoam-metrics] could not determine the droplet centre")

    # fieldMinMax, mode magnitude, writes one line per FIELD per time:
    #
    #   # Time   field   min   location(min)   max   location(max)
    #   2.12e-05   mag(U)   0.0   (0 7.8e-05 0)   3.47e-02   (4.6e-03 6.0e-03 0)
    #
    # Two traps, both hit on the first attempt: the field is named `mag(U)` and
    # not `U` under mode magnitude, and the location columns contain SPACES inside
    # their parentheses, so a whitespace split does not line up with the header.
    # Stripping every parenthesised group first collapses `mag(U)` -> `mag` and
    # each location -> nothing, leaving exactly [time, field, min, max].
    umax = {}
    for d in _latest_dir(os.path.join(case, "postProcessing", "parasiticU", "*")):
        for f in sorted(glob.glob(os.path.join(d, "*.dat"))):
            for line in open(f):
                if line.startswith("#") or not line.strip():
                    continue
                raw = line.split()
                if len(raw) < 2:
                    continue
                field = raw[1]
                if field not in ("mag(U)", "U"):
                    continue
                flat = re.sub(r"\([^)]*\)", " ", line).split()
                if len(flat) < 4:
                    continue
                try:
                    umax[round(float(flat[0]), 12)] = abs(float(flat[3]))
                except ValueError:
                    continue
    vol = {round(r[0], 12): r[1] for r in collect(case, "phaseVolume",
                                                 "volFieldValue.dat")}
    if not vol:
        f = find_file(case, "phaseVolume", "*.dat")
        if f:
            vol = {round(r[0], 12): r[1]
                   for r in collect(case, "phaseVolume", f)}
    probes = collect(case, "pressureProbes", "p_rgh", 3)
    pjump = {round(r[0], 12): r[1] - r[2] for r in probes if len(r) >= 3}
    shape = shape_error(case, centre, args.radius)

    times = sorted(set(umax) | set(vol))
    if not times:
        sys.exit("[interfoam-metrics] no function-object output found under "
                 f"{os.path.join(case, 'postProcessing')}")
    v0 = vol.get(times[0])
    cols = ["TIME", "maxMagU", "phaseVolume", "phaseVolumeRelError", "pLaplace",
            "zeroSetRadialL2"]
    with open(os.path.join(case, args.out), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for t in times:
            v = vol.get(t)
            w.writerow([
                f"{t:.10g}",
                f"{umax.get(t, float('nan')):.10g}",
                f"{v:.10g}" if v is not None else "",
                f"{abs(v-v0)/v0:.10g}" if (v is not None and v0) else "",
                f"{pjump.get(t, float('nan')):.10g}" if t in pjump else "",
                f"{shape.get(t, float('nan')):.10g}" if t in shape else "",
            ])
    print(f"[interfoam-metrics] wrote {args.out}: {len(times)} rows, "
          f"centre={centre}, R={args.radius}, "
          f"{len(shape)} isosurface instants")
    return 0


if __name__ == "__main__":
    sys.exit(main())
