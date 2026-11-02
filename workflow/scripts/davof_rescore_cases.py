#!/usr/bin/env python3
"""The fixed shape error of a list of cases against their own former one
(2026-11-02). For every case of --list (one case directory per line) the
app's output <rescore>/<study>/<case>.dat (staticFaceFractions -entry
volumeFractionError -davofMethod exactSphere -evolOnly over the written
times) is matched with the case's volumeFractionError.dat (the former Eshape,
cutCellIso, unless its header names exactState exactSphere) at the same
times. One row per case and written time:

    study case row arm N flux Co face t Eshape_old construction_old
    Eshape_new Ec Eout

plus the summary per case at the end time and, per (study, row, arm, flux,
Co), the order of the new end value over N.

    python3 workflow/scripts/davof_rescore_cases.py --list L --rescore R --csv out.csv [--summary sum.csv]
"""
import argparse
import csv
import json
import math
import os


def table(path):
    names, rows = None, []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            q = s[1:].split()
            if q and q[0] == "Time":
                names = q
            continue
        rows.append(dict(zip(names, map(float, s.split()[:len(names)]))))
    return rows


def order(ns, es):
    pts = [(1.0 / n, e) for n, e in zip(ns, es) if e and e > 0]
    if len(pts) < 2:
        return None
    lh = [math.log(p[0]) for p in pts]
    le = [math.log(p[1]) for p in pts]
    mh, me = sum(lh) / len(lh), sum(le) / len(le)
    sxx = sum((a - mh) ** 2 for a in lh)
    return sum((a - mh) * (b - me) for a, b in zip(lh, le)) / sxx if sxx > 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True)
    ap.add_argument("--rescore", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--summary", default="")
    a = ap.parse_args()
    out, ends = [], []
    for c in [l.strip() for l in open(a.list) if l.strip()]:
        study, case = os.path.basename(os.path.dirname(c)), os.path.basename(c)
        app = os.path.join(a.rescore, study, case + ".dat")
        if not os.path.isfile(app):
            print("missing", app)
            continue
        tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        dat = os.path.join(c, "postProcessing", "volumeFractionError", "0", "volumeFractionError.dat")
        old, cons = {}, "-"
        if os.path.isfile(dat):
            cons = "exactSphere" if "exactState exactSphere" in open(dat).readline() else "cutCellIso"
            for r in table(dat):
                old[round(r["Time"], 9)] = r["Eshape"]
        status = "-"
        conv = os.path.join(c, "advectorVoF.csv")
        if os.path.isfile(conv):
            rr = list(csv.DictReader(open(conv)))
            if rr:
                status = rr[-1].get("status", "-")
        base = dict(study=study, case=case, row=tok.get("TPF_ROW"), arm=tok.get("TPF_VELOCITY_ARM"),
                    N=int(float(tok["N_CELLS"])), flux=tok.get("TPF_FLUX_INTEGRATION", "-"),
                    Co=tok.get("TARGET_ADVECTIVE_CO"), face=tok.get("TPF_SWEPT_FACE_UPDATE", "-"),
                    status=status, construction_old=cons)
        rows = table(app)
        for r in rows:
            t = r["Time"]
            out.append(dict(base, t=t, Eshape_old=old.get(round(t, 9)), Eshape_new=r["Evol"], Ec=r["Ec"], Eout=r["Eout"]))
        if rows:
            ends.append(dict(base, t=rows[-1]["Time"], Eshape_old=old.get(round(rows[-1]["Time"], 9)),
                             Eshape_new=rows[-1]["Evol"], Ec=rows[-1]["Ec"], Eshape_new_t0=rows[0]["Evol"]))
    keys = ["study", "case", "row", "arm", "N", "flux", "face", "Co", "status", "t", "Eshape_old", "construction_old", "Eshape_new", "Ec", "Eout"]
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    # the orders of the new end value over N
    groups = {}
    for e in ends:
        groups.setdefault((e["study"], e["row"], e["arm"], e["flux"], e["face"], e["Co"]), []).append(e)
    for g in groups.values():
        g.sort(key=lambda e: e["N"])
        done = [e for e in g if e["status"] in ("completed", "-")]
        po = order([e["N"] for e in done], [e["Eshape_old"] for e in done]) if all(e["Eshape_old"] for e in done) else None
        pn = order([e["N"] for e in done], [e["Eshape_new"] for e in done])
        for e in g:
            e["order_old"], e["order_new"] = po, pn
    f = lambda v: ("%.4e" % v) if isinstance(v, float) else str(v)
    skeys = ["study", "row", "arm", "flux", "face", "Co", "N", "status", "t", "Eshape_old", "construction_old", "Eshape_new", "Eshape_new_t0", "Ec", "order_old", "order_new"]
    print("\t".join(skeys))
    for e in sorted(ends, key=lambda e: (e["study"], str(e["row"]), str(e["arm"]), e["N"])):
        print("\t".join(f(e.get(k)) for k in skeys))
    if a.summary:
        with open(a.summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=skeys, extrasaction="ignore")
            w.writeheader()
            w.writerows(ends)


if __name__ == "__main__":
    main()
