#!/usr/bin/env python3
"""Compare the cases of two study directories for bit-identity of their metric CSVs.

    integration_identity_check.py <study dir A> <study dir B> [--match tokens|name]
        [--only-a 00002,00003] [--skip COLS] [--out table.csv] [--label TEXT]

A and B hold the same study, run by two builds (for example the merged tree and the
pre-merge tree). The script pairs the cases, then compares every CSV in the top level of
each case directory (the solver CSV, gradPsiError.csv, leiaSetFields.csv and any other
per-step CSV).

Pairing. `--match tokens` (default) pairs two cases when the tokens of their
case_params.json are equal on the keys that both studies have. A token that only one tree
has (a new token with an inert default) is listed and is not used. `--match name` pairs the
directories with the same name (the polyhedral arms of advect_bound_arm.sh carry the
case_params.json of their mesh case, so their tokens do not separate them).

Comparison of one CSV pair, both parts must pass:
  1. workflow/scripts/compare_metrics_csv.py A B --tol 0 --skip <wall-clock columns>
  2. an exact text comparison of every cell outside the skipped columns. The comparator of
     item 1 does not see a difference between two non-numeric cells or two NaN cells; the
     text comparison does.
For a difference the script reports the column, the first row and the largest relative
difference. The log state comes from workflow/scripts/foam_log_state.sh. A case that
diverged in both trees with equal CSVs is identical.

Exit code: 0 all pairs identical, 2 a difference, 1 a structural problem (a case without a
partner, a missing CSV, no case found).
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_DEFAULT = "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"


def cases(study, match):
    out = {}
    for d in sorted(glob.glob(os.path.join(study, "*"))):
        if not os.path.isdir(d):
            continue
        cp = os.path.join(d, "case_params.json")
        if match == "tokens" and not os.path.isfile(cp):
            continue
        if match == "name" and not glob.glob(os.path.join(d, "*.csv")) and not os.path.isfile(cp):
            continue
        tokens = json.load(open(cp)).get("tokens", {}) if os.path.isfile(cp) else {}
        out[os.path.basename(d)] = (d, tokens)
    return out


def log_state(d):
    logs = [p for p in glob.glob(os.path.join(d, "log.solve")) + glob.glob(os.path.join(d, "log.leia*Foam"))]
    if not logs:
        return "NO_LOG"
    r = subprocess.run(["bash", os.path.join(HERE, "foam_log_state.sh"), logs[0], "--stall", "1"],
                       capture_output=True, text=True)
    words = r.stdout.split()
    keep = [w for w in words if w.startswith(("steps=", "nprocs="))]
    return " ".join(words[:1] + keep) if words else "NO_STATE"


def md5(path):
    for p in (path, path + ".gz"):
        if os.path.isfile(p):
            h = hashlib.md5()
            with open(p, "rb") as fh:
                for blk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(blk)
            return h.hexdigest()[:12]
    return "none"


def load(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    return (rows[0], rows[1:]) if rows else ([], [])


def exact(pa, pb, skip):
    """Return (identical, rows_a, rows_b, message)."""
    ha, ra = load(pa)
    hb, rb = load(pb)
    if ha != hb:
        return False, len(ra), len(rb), "headers differ: only A %s, only B %s" % (
            sorted(set(ha) - set(hb)), sorted(set(hb) - set(ha)))
    if len(ra) != len(rb):
        return False, len(ra), len(rb), "row counts differ: %d against %d" % (len(ra), len(rb))
    first, worst = None, {}
    tcol = ha.index("TIME") if "TIME" in ha else None
    for i, (a, b) in enumerate(zip(ra, rb)):
        if a == b:
            continue
        for j, name in enumerate(ha):
            if name in skip or j >= len(a) or j >= len(b) or a[j] == b[j]:
                continue
            try:
                va, vb = float(a[j]), float(b[j])
                rel = abs(va - vb) / max(abs(va), abs(vb), 1e-300)
            except ValueError:
                rel = float("inf")
            if first is None:
                t = a[tcol] if tcol is not None else "?"
                first = "first difference: column %s, row %d (TIME %s): %s against %s" % (name, i, t, a[j], b[j])
            if rel > worst.get(name, (0.0, 0))[0]:
                worst[name] = (rel, i)
    if first is None:
        return True, len(ra), len(rb), ""
    top = sorted(worst.items(), key=lambda kv: -kv[1][0])[:4]
    return False, len(ra), len(rb), first + "; largest relative differences: " + ", ".join(
        "%s %.3e (row %d)" % (n, m, i) for n, (m, i) in top) + "; %d column(s) differ" % len(worst)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--match", choices=("tokens", "name"), default="tokens")
    ap.add_argument("--only-a", default="", help="comma-separated suffixes of the case directories of A to use")
    ap.add_argument("--only-b", default="", help="the same for B")
    ap.add_argument("--skip", default=SKIP_DEFAULT)
    ap.add_argument("--csvs", default="", help="comma-separated CSV names; default: every top-level CSV")
    ap.add_argument("--out", default="")
    ap.add_argument("--label", default="")
    ap.add_argument("--show", default="", help="comma-separated tokens to print per case; default: the tokens that vary")
    a = ap.parse_args()
    skip = {c for c in a.skip.split(",") if c}
    A, B = cases(a.a, a.match), cases(a.b, a.match)
    for only, D in ((a.only_a, A), (a.only_b, B)):
        if only:
            suf = tuple(s for s in only.split(",") if s)
            for k in [k for k in D if not k.endswith(suf)]:
                del D[k]
    if not A or not B:
        print("[identity] no case found in %s or in %s" % (a.a, a.b))
        return 1
    label = a.label or os.path.basename(os.path.normpath(a.a))
    structural = False
    if a.match == "tokens":
        ka = set.intersection(*(set(t) for _, t in A.values()))
        kb = set.intersection(*(set(t) for _, t in B.values()))
        common = sorted(ka & kb)
        if ka ^ kb:
            print("[identity] %s: tokens only in A: %s; only in B: %s (not used for the pairing)" % (
                label, sorted(ka - kb), sorted(kb - ka)))
        key = lambda t: tuple((k, t[k]) for k in common)
        ib = {}
        for name, (d, t) in B.items():
            if key(t) in ib:
                print("[identity] %s: two cases of B have equal tokens: %s, %s" % (label, ib[key(t)], name))
                structural = True
            ib[key(t)] = name
        pairs = []
        for name, (d, t) in A.items():
            if key(t) in ib:
                pairs.append((name, ib.pop(key(t))))
            else:
                print("[identity] %s: case %s of A has no partner in B" % (label, name))
                structural = True
        for name in ib.values():
            print("[identity] %s: case %s of B has no partner in A" % (label, name))
            structural = True
        vary = [k for k in common if len({t[k] for _, t in A.values()}) > 1]
    else:
        pairs = [(n, n) for n in A if n in B]
        for n in sorted(set(A) ^ set(B)):
            print("[identity] %s: case %s is in one study only" % (label, n))
            structural = True
        vary = []
    show = [s for s in a.show.split(",") if s] or vary
    rows, differ = [], False
    for na, nb in pairs:
        da, db = A[na][0], B[nb][0]
        par = " ".join("%s=%s" % (k, A[na][1].get(k, "?")) for k in show) or na
        names = ([c for c in a.csvs.split(",") if c] or
                 sorted({os.path.basename(p) for p in glob.glob(os.path.join(da, "*.csv")) +
                         glob.glob(os.path.join(db, "*.csv"))}))
        sa, sb = log_state(da), log_state(db)
        mesh = (md5(os.path.join(da, "constant/polyMesh/points")), md5(os.path.join(db, "constant/polyMesh/points")))
        if not names:
            print("[identity] %s: %s has no CSV" % (label, na))
            structural = True
        for n in names:
            pa, pb = os.path.join(da, n), os.path.join(db, n)
            if not (os.path.isfile(pa) and os.path.isfile(pb)):
                rows.append([label, na, nb, par, n, "", "", sa, sb, mesh[0], mesh[1], "MISSING", "", "the CSV is missing in one tree"])
                structural = True
                continue
            if os.path.getsize(pa) == 0 and os.path.getsize(pb) == 0:
                rows.append([label, na, nb, par, n, 0, 0, sa, sb, mesh[0], mesh[1], "yes", "both empty", ""])
                continue
            r = subprocess.run([sys.executable, os.path.join(HERE, "compare_metrics_csv.py"), pa, pb,
                                "--tol", "0", "--skip", ",".join(sorted(skip))], capture_output=True, text=True)
            same, ra, rb, msg = exact(pa, pb, skip)
            ok = same and r.returncode == 0
            if not ok:
                differ = True
                if not msg:
                    msg = "comparator: " + (r.stdout.strip().splitlines() or ["rc %d" % r.returncode])[-1]
            rows.append([label, na, nb, par, n, ra, rb, sa, sb, mesh[0], mesh[1], "yes" if ok else "NO",
                         "comparator rc %d" % r.returncode, msg])
    head = ["study", "caseA", "caseB", "parameters", "csv", "rowsA", "rowsB", "logStateA", "logStateB",
            "meshA", "meshB", "identical", "comparator", "difference"]
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for r in rows:
        print("| " + " | ".join(str(x) for x in r) + " |")
    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(head)
            w.writerows(rows)
    n_ok = sum(1 for r in rows if r[11] == "yes")
    print("[identity] %s: %d of %d CSV pairs identical in %d case pairs" % (label, n_ok, len(rows), len(pairs)))
    return 1 if structural else (2 if differ else 0)


if __name__ == "__main__":
    sys.exit(main())
