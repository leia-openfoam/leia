#!/usr/bin/env python3
"""The clipRegion gate: is the band-aware quasi-monotone clip inert at the interface?

Reads a clipRegion gate study whose arms are the cross product of SL_CLIP
[false, true] and SL_CLIP_REGION [all, outsideBand], and answers the three
pre-registered questions in the config header:

  1. TOKEN INERTNESS. (false, outsideBand) must be BIT-IDENTICAL to (false, all):
     the region is read only when the clip is on. A difference means the new token
     is not inert, and no other row may be believed.
  2. THE CONTROL SEES THE DAMAGE. (true, all) is the GLOBAL clip, already measured
     to cost +30 % volume error on Popinet's 2D hexahedral translating droplet at
     N = 64. It must move the interface metrics, or the matrix cannot show that
     (true, outsideBand) avoided anything.
  3. THE CANDIDATE IS INERT. (true, outsideBand) must match (false, all) to within
     the tolerance, on every interface metric.

Reports L1 and L2 of the parasitic velocity, the volume error, the shape error, the
centroid error, the Laplace pressure and the band curvature error -- together, never
one alone. Lengths are divided by DROPLET_RADIUS, because zeroSetRadialL2 and
centroidError are absolute lengths in metres and a fixed threshold on them reads every
SI run as clean. L_inf is never reported: it does not converge.

Also reports the clip's own activity counter from the solver log,

    slCorrector: quasi-monotone clip (clipRegion = ...): N eligible, M bounded

because M = 0 in the candidate arm makes the gate vacuous -- the clip never acted, so
it proves nothing.

Arms that did not reach the same step count as the baseline are reported as incomplete
and excluded: a number nobody has seen land is not reported.

Usage:  python3 workflow/scripts/make_clip_region_gate_table.py <study> [--root .]
            [--tol 0.01] [--out <dir>]
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import re
import sys

SOLVER_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
SOLVER_LOG = "log.leiaSemiLagrangianLevelSetTwoPhaseFoam"

# (column, label, normalise by the droplet radius?)
METRICS = (
    ("meanMagUPrime", "L1|u'|", False),
    ("l2MagUPrime", "L2|u'|", False),
    ("phaseVolumeRelError", "volume err", False),
    ("zeroSetRadialL2", "shape L2/R", True),
    ("centroidError", "centroid/R", True),
    ("pLaplace", "p Laplace", False),
    ("kErrL2Band", "kappa L2 band", False),
)

CLIP_LINE = re.compile(
    r"quasi-monotone clip \(clipRegion = (\w+)\):\s+(\d+) eligible cells,\s+(\d+) bounded"
)


def _f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def load_arms(root, study):
    """{(clip, region, keepExtrema): {...}} keyed by the three clip tokens."""
    arms = {}
    pattern = os.path.join(root, "studies", study, "*_[0-9]*")
    for d in sorted(glob.glob(pattern)):
        cp = os.path.join(d, "case_params.json")
        csv_path = os.path.join(d, SOLVER_CSV)
        if not (os.path.isfile(cp) and os.path.isfile(csv_path)):
            continue
        with open(cp) as fh:
            tokens = json.load(fh).get("tokens", {})
        clip = str(tokens.get("SL_CLIP", "?")).lower()
        region = str(tokens.get("SL_CLIP_REGION", "all"))
        keep = str(tokens.get("SL_CLIP_KEEP_EXTREMA", "false")).lower()
        with open(csv_path) as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            continue
        arms[(clip, region, keep)] = {
            "dir": d,
            "csv": csv_path,
            "rows": rows,
            "nsteps": len(rows),
            "radius": _f(tokens.get("DROPLET_RADIUS"), 1.0),
            "mesh": mesh_digest(d),
        }
    return arms


def mesh_digest(arm_dir):
    """Hash of constant/polyMesh/points, or None when it is absent.

    A bit-identity criterion is only meaningful when the arms SHARE a mesh. blockMesh
    is deterministic, so hexahedral arms do. cfMesh is NOT bit-reproducible and the
    workflow rebuilds the mesh per arm, so polyhedral arms each get their own -- and
    then a byte difference in the metrics measures the mesher, not the tokens. This
    hash is what separates the two readings.
    """
    path = os.path.join(arm_dir, "constant", "polyMesh", "points")
    if not os.path.isfile(path):
        return None
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clip_activity(arm_dir):
    """The last 'N eligible, M bounded' line of the solver log, or None."""
    log = os.path.join(arm_dir, SOLVER_LOG)
    if not os.path.isfile(log):
        return None
    last = None
    with open(log, errors="replace") as fh:
        for line in fh:
            m = CLIP_LINE.search(line)
            if m:
                last = (m.group(1), int(m.group(2)), int(m.group(3)))
    return last


def rel(a, b):
    """Relative difference of a against the reference b."""
    if b is None or a is None:
        return None
    if b == 0.0:
        return 0.0 if a == 0.0 else float("inf")
    return (a - b) / abs(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--root", default=".")
    ap.add_argument("--tol", type=float, default=0.01,
                    help="pass tolerance on every interface metric (default 1 %%)")
    ap.add_argument("--out", default=None,
                    help="directory for clip_region_gate.csv (default: print only)")
    args = ap.parse_args()

    arms = load_arms(args.root, args.study)
    if not arms:
        sys.exit(f"no arms with a metrics CSV under studies/{args.study}")

    # The baseline is ANY clip-off arm: the region and the extremum exemption are read
    # only when the clip is on, so every clip-off arm must be bit-identical (which is
    # criterion 1 below). Prefer region "all" for a stable column order; a study that
    # holds one of the two tokens fixed is then read without editing this script.
    off_keys = sorted(k for k in arms if k[0] == "false")
    if not off_keys:
        sys.exit("no clip-off arm: there is nothing to compare the clip against")
    base_key = next((k for k in off_keys if k[1] == "all"), off_keys[0])
    base = arms[base_key]
    radius = base["radius"]

    print(f"clipRegion gate: {args.study}")
    print(f"  arms found: {len(arms)}   baseline steps: {base['nsteps'] - 1}"
          f"   droplet radius: {radius} m   tolerance: {args.tol:.1%}")
    print()

    # ---- arm inventory, completeness and clip activity ------------------
    print("  arm                                        steps   clip activity (last write)")
    for key in sorted(arms):
        a = arms[key]
        act = clip_activity(a["dir"])
        act_s = (f"{act[1]} eligible, {act[2]} bounded" if act
                 else "-- (clip off, no counter)")
        flag = "" if a["nsteps"] == base["nsteps"] else "  INCOMPLETE"
        print(f"  clip={key[0]:<5} region={key[1]:<12} keepExtrema={key[2]:<5}"
              f" {a['nsteps'] - 1:>6}   {act_s}{flag}")
    print()

    complete = {k: v for k, v in arms.items() if v["nsteps"] == base["nsteps"]}

    # ---- 1. token inertness: every clip-off arm must agree bitwise -------
    verdicts = []
    off = off_keys
    with open(base["csv"], "rb") as fh:
        ref = fh.read()
    diffs = []
    for k in sorted(off):
        if k == base_key:
            continue
        with open(arms[k]["csv"], "rb") as fh:
            if fh.read() != ref:
                diffs.append(k)
    shared_mesh = len({arms[k]["mesh"] for k in off}) == 1 and arms[base_key]["mesh"]
    if len(off) < 2:
        print("  1. token inertness   only one clip-off arm -- NOT TESTED")
        verdicts.append(("token inertness", None))
    elif not diffs:
        print(f"  1. token inertness   {len(off)} clip-off arms: ALL BIT-IDENTICAL -- PASS")
        verdicts.append(("token inertness", True))
    elif not shared_mesh:
        # Each arm built its own mesh, so a byte difference measures the MESHER, not the
        # tokens. cfMesh is not bit-reproducible. Report the largest relative difference
        # on the primary metrics instead, and leave the criterion untested.
        worst, worst_metric = 0.0, "-"
        for col, label, norm in METRICS:
            x = _f(arms[base_key]["rows"][-1].get(col))
            for k in diffs:
                y = _f(arms[k]["rows"][-1].get(col))
                r = rel(y, x)
                if r is not None and abs(r) > worst:
                    worst, worst_metric = abs(r), label
        print(f"  1. token inertness   NOT TESTABLE -- the clip-off arms have DIFFERENT"
              f" MESHES (cfMesh is not bit-reproducible and the workflow rebuilds per arm)."
              f" Largest relative difference at the last step: {worst:.3e} on {worst_metric}."
              f" Test bit-identity on a blockMesh study instead.")
        verdicts.append(("token inertness", None))
    else:
        print(f"  1. token inertness   {len(off)} clip-off arms share a mesh but DIFFER"
              f" -- FAIL {diffs}")
        verdicts.append(("token inertness", False))

    # ---- 2 & 3. the metric table at the last common step ----------------
    print()
    others = [k for k in sorted(complete) if k != base_key and k[0] == "true"]
    header = f"  {'metric':<16}{'baseline':>13}"
    for key in others:
        tag = ("bnd" if key[1] == "outsideBand" else "all") + \
              ("+ext" if key[2] == "true" else "")
        header += f"{tag:>22}"
    print(header)

    table = []
    for col, label, norm in METRICS:
        b = _f(base["rows"][-1].get(col))
        if b is not None and norm:
            b /= radius
        line = f"  {label:<16}{b:>13.4e}" if b is not None else f"  {label:<16}{'--':>13}"
        row = {"metric": label, "baseline": b}
        for key in others:
            v = _f(complete[key]["rows"][-1].get(col))
            if v is not None and norm:
                v /= radius
            r = rel(v, b)
            tag = f"{key[0]}/{key[1]}/{key[2]}"
            row[tag] = v
            row[tag + "_rel"] = r
            line += (f"{v:>13.4e} {r:>+7.1%}" if v is not None and r is not None
                     else f"{'--':>22}")
        print(line)
        table.append(row)

    # ---- the pre-registered verdicts ------------------------------------
    print()
    # Arms with keepExtrema false are CONTROLS and must MOVE the metrics (they flatten the
    # level set's own extrema); arms with it true are candidates and must NOT move them.
    checks = []
    for key in sorted(k for k in complete if k[0] == "true"):
        region = "bnd" if key[1] == "outsideBand" else "all"
        if key[2] == "false":
            checks.append((key, f"control, must move   ({region}, keep=false)", True))
        else:
            checks.append((key, f"candidate, must not  ({region}, keep=true) ", False))
    for i, (key, name, want_move) in enumerate(checks, start=2):
        name = f"{i}. {name}"
        worst, worst_metric = 0.0, "-"
        for row in table:
            r = row.get(f"{key[0]}/{key[1]}/{key[2]}_rel")
            if r is not None and abs(r) > worst:
                worst, worst_metric = abs(r), row["metric"]
        if want_move:
            ok = worst > args.tol
            print(f"  {name}   largest move {worst:.1%} on {worst_metric} -- "
                  f"{'PASS (the metrics respond to a clip at the interface)' if ok else 'FAIL (the matrix cannot see the damage)'}")
        else:
            ok = worst <= args.tol
            print(f"  {name}   largest move {worst:.1%} on {worst_metric} -- "
                  f"{'PASS' if ok else 'FAIL'}")
        verdicts.append((name, ok))

    print()
    decided = [v for _, v in verdicts if v is not None]
    if len(decided) == len(verdicts) and all(decided):
        print("  GATE: PASS -- every pre-registered criterion is met.")
    elif any(v is False for _, v in verdicts):
        print("  GATE: FAIL -- see the failing criterion above.")
    else:
        untested = [n for n, v in verdicts if v is None]
        print(f"  GATE: PARTIAL -- every criterion that could be tested PASSED;"
              f" not tested here: {', '.join(untested)}."
              f" Do not read the untested criterion as a pass.")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "clip_region_gate.csv")
        keys = ["metric", "baseline"] + [k for k in table[0] if k not in ("metric", "baseline")]
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(table)
        print(f"\n  wrote {path}")


if __name__ == "__main__":
    main()
