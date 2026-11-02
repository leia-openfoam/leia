#!/usr/bin/env python3
"""The re-scoring of the preprint's shape errors (2026-11-03; arXiv v1, tag
davof-v1.0). The old E_shape is the volumeFractionError of the runs, whose
alpha_exact was the cutCellIso cut of the implicit function's vertex values;
the new one is sum |alpha - alpha_exact| V against the exact sphere
(davofInit::exactSphere), the Eshape of the fixed volumeFractionError
(fork davof-swept 1f3f238), computed from the written fields by the app
staticFaceFractions -entry volumeFractionError -davofMethod exactSphere
-evolOnly at t = 0 and the end time, without re-solving.

Inputs: --article DIR (the v1 data tables davof_advection_<study>.csv: the
old Eshape_end exactly as tabulated), --rescore DIR (<study>/<case>.dat of the
app). Per (study, row, Co): the end values per N, the least-squares order over
N in {32, 64, 128} (make_table_advection._order), and the pairwise orders.

    python3 workflow/scripts/davof_rescore_eshape.py --article A --rescore R [--csv out.csv]
"""
import argparse
import csv
import glob
import math
import os


def order(ns, es):
    pts = [(1.0 / n, e) for n, e in zip(ns, es) if e and e > 0]
    if len(pts) < 2:
        return None
    lh = [math.log(p[0]) for p in pts]
    le = [math.log(p[1]) for p in pts]
    mh, me = sum(lh) / len(lh), sum(le) / len(le)
    sxx = sum((a - mh) ** 2 for a in lh)
    return sum((a - mh) * (b - me) for a, b in zip(lh, le)) / sxx if sxx > 0 else None


def read_app(path):
    rows = []
    names = None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            q = s[1:].split()
            if q and q[0] == "Time":
                names = q
            continue
        rows.append(dict(zip(names, map(float, s.split()))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--article", required=True)
    ap.add_argument("--rescore", required=True)
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    out = []
    for p in sorted(glob.glob(os.path.join(a.article, "davof_advection_*.csv"))):
        study = os.path.basename(p)[len("davof_advection_"):-4]
        for r in csv.DictReader(open(p, newline="")):
            app = os.path.join(a.rescore, study, r["case"] + ".dat")
            rec = dict(study=study, case=r["case"], row=r["row"], N=int(r["N"]), Co=r["Co"], status=r["status"],
                       old=float(r["Eshape_end"]) if r.get("Eshape_end") else None, new=None, new_t0=None, tEnd=None)
            if os.path.isfile(app):
                rows = read_app(app)
                if rows:
                    rec["new_t0"] = rows[0]["Evol"] if abs(rows[0]["Time"]) < 1e-12 else None
                    rec["new"] = rows[-1]["Evol"]
                    rec["tEnd"] = rows[-1]["Time"]
                    rec["Ec_end"] = rows[-1]["Ec"]
            out.append(rec)

    # orders per (study, row, Co)
    groups = {}
    for r in out:
        groups.setdefault((r["study"], r["row"], r["Co"]), []).append(r)
    print("study\trow\tCo\tN\tstatus\tE_old\tE_new\tE_new(t=0)\tnew/old")
    summary = []
    for key in sorted(groups):
        g = sorted(groups[key], key=lambda r: r["N"])
        for r in g:
            f = lambda v: "%.3e" % v if isinstance(v, float) else "-"
            ratio = r["new"] / r["old"] if r["new"] is not None and r["old"] else None
            print("\t".join([key[0], key[1], key[2], str(r["N"]), r["status"], f(r["old"]), f(r["new"]), f(r["new_t0"]), f(ratio)]))
        done = [r for r in g if r["status"] == "completed" and r["new"] is not None and r["N"] in (32, 64, 128)]
        if len(done) >= 2:
            ns = [r["N"] for r in done]
            po = order(ns, [r["old"] for r in done])
            pn = order(ns, [r["new"] for r in done])
            pair_old = [order([x["N"], y["N"]], [x["old"], y["old"]]) for x, y in zip(done, done[1:])]
            pair_new = [order([x["N"], y["N"]], [x["new"], y["new"]]) for x, y in zip(done, done[1:])]
            fo = lambda v: "%.2f" % v if v is not None else "-"
            print("  -> order over N %s: old %s, new %s; pairwise old %s, new %s"
                  % (ns, fo(po), fo(pn), "/".join(map(fo, pair_old)), "/".join(map(fo, pair_new))))
            for r in g:
                r["order_old"], r["order_new"] = po, pn
        summary.append(key)

    if a.csv:
        keys = ["study", "case", "row", "N", "Co", "status", "old", "new", "new_t0", "tEnd", "Ec_end", "order_old", "order_new"]
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
