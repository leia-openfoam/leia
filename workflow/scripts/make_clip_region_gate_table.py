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
        }
    return arms


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

    base_key = ("false", "all", "false")
    if base_key not in arms:
        sys.exit("the baseline arm (SL_CLIP false, SL_CLIP_REGION all) is missing")
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
    off = [k for k in arms if k[0] == "false"]
    with open(base["csv"], "rb") as fh:
        ref = fh.read()
    diffs = []
    for k in sorted(off):
        if k == base_key:
            continue
        with open(arms[k]["csv"], "rb") as fh:
            if fh.read() != ref:
                diffs.append(k)
    if len(off) > 1:
        ok = not diffs
        print(f"  1. token inertness   {len(off)} clip-off arms: "
              f"{'ALL BIT-IDENTICAL -- PASS' if ok else 'DIFFER -- FAIL ' + str(diffs)}")
        verdicts.append(("token inertness", ok))
    else:
        print("  1. token inertness   only one clip-off arm -- NOT TESTED")
        verdicts.append(("token inertness", None))

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
    for key, name, want_move in (
        (("true", "all", "false"), "2. control sees the damage    (all, keep=false)", True),
        (("true", "outsideBand", "false"), "3. control reproduces it (bnd, keep=false)", True),
        (("true", "all", "true"), "4. extremum exemption alone   (all, keep=true) ", False),
        (("true", "outsideBand", "true"), "5. THE CANDIDATE          (bnd, keep=true) ", False),
    ):
        if key not in complete:
            print(f"  {name}   arm {key} incomplete or absent -- NOT TESTED")
            verdicts.append((name, None))
            continue
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
        print("  GATE: UNDECIDED -- an arm has not landed. Do not read a verdict.")

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
