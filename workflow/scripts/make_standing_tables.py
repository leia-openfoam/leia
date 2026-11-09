#!/usr/bin/env python3
"""The tables of the standing study matrix (Tomislav, 2026-11-07; docs/davof/PORT-DESIGN.md
5b; the configs config/davof/standing/<family>_N<N>.yaml from make_standing_configs.py).

    python3 workflow/scripts/make_standing_tables.py <family> [--studies-dir DIR] [--out DIR]

Reads every case of studies/davof/standing_<family>_N*/ and writes, under
docs/davof/data/standing/ (the davof theme's data folder):

  standing_<family>.csv          one row per (row, N, OF Co): the end-of-run errors, the
                                 bounds and the mass over the run, the counts, the wall
                                 clock (the median, min and max over the repetitions) and the
                                 per-step time, the status, the stamps, the repetitions'
                                 identity verdict
  standing_<family>.tex          the same as a booktabs tabular body (\\input it)
  standing_<family>_orders.csv   the orders over N per (row, Co) of E_shape, L1(x), L1(n)
  standing_<family>_orders.tex   (least squares over the completed rungs; pairwise too)
  standing_<family>_caption.tex  the caption's provenance: the stamps (the ids per rung,
                                 the solver binary, leia's describe), the grazingCut key,
                                 dt and the steps per rung, the ranks

The measures (the converter tpf_advection_metrics.py's columns): E_shape (volumeFractionError
against the exact end state), L1(x) = LCentreDiff1, L1(n) = LNormalDiff1, L1(alphaf) (the
DAVOF rows; -1 for the comparators), the mass deviation devPercent, the bounds min alpha and
the largest excess over 1 from the log (every step, after the conservative bounding), the
interface-cell count at the end, the OpenFOAM Courant number (the log's maximum), the status
(ended | aborted: the abort guard's step and reason), the wall clock of the solve step from
.leia_launch and the median per-step time (timeStepTime). The repetitions must agree in
every result column; the table says so ("identical") or names the differing one.
"""
import argparse, csv, glob, json, math, os, re, statistics, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RESULT_COLS = ("Eshape", "L1x", "L1n", "L1af", "massDevPct", "minAlpha", "excess", "nInterfaceEnd", "steps", "status", "abortStep")
LABELS = {"davof-planeSurf-ringArea": "DAVOF (swept-face, ring normal)", "plicRDF": "isoAdvector + plicRDF (v2606)",
          "davof-planeSurf-ringArea [recommended]": "DAVOF, signedTets + scale + normalAndAlpha",
          "davof-planeSurf-ringArea [signedTets]": "DAVOF, signedTets",
          "davof-planeSurf-ringArea [closureThenTets]": "DAVOF, closureThenTets",
          "davof-planeSurf-ringArea [onePlane]": "DAVOF, one plane per face (KISS)"}


def f(x, nd=4):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "-"
    if v != v:
        return "-"
    return "%.*e" % (nd, v)


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def read_case(c):
    p = os.path.join(c, "case_params.json")
    if not os.path.isfile(p):
        return None
    t = json.load(open(p))["tokens"]
    # the row of the table = the TPF_ROW plus the Enright setting (the variant rows)
    setting = t.get("TPF_ENRIGHT_SETTING", "reference")
    row = t.get("TPF_ROW") + ("" if setting in ("reference", None) else " [" + setting + "]")
    rec = dict(case=c, row=row, setting=setting, N=int(float(t.get("N_CELLS", 0))), co=float(t.get("TARGET_ADVECTIVE_CO", 0)),
               rep=int(float(t.get("TPF_REPETITION", 1))), grazingCut=t.get("TPF_GRAZING_CUT", ""), dt=t.get("MAX_DELTA_T", ""), T=t.get("END_TIME", ""))
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) and os.path.getsize(csvp) else []
    last = rows[-1] if rows else {}
    rec.update(Eshape=num(last.get("Eshape")), L1x=num(last.get("LCentreDiff1")), L1n=num(last.get("LNormalDiff1")),
               L1af=num(last.get("L1_alphaf")) if rec["row"].startswith("davof") else None, massDevPct=num(last.get("devPercent")),
               nInterfaceEnd=num(last.get("nInterface")), tLast=num(last.get("TIME")), status=last.get("status", "noCSV"), abortStep=last.get("abortStep", ""))
    st = [num(r.get("timeStepTime")) for r in rows]
    st = [x for x in st if x is not None and x > 0]
    rec["stepMedian"] = statistics.median(st) if st else None
    logp = next((os.path.join(c, n) for n in ("log.davofAdvectorFoam", "log.advectorVoF") if os.path.isfile(os.path.join(c, n))), None)
    rec.update(minAlpha=None, excess=None, ofCo=None, steps=None)
    if logp:
        log = open(logp, errors="replace").read()
        mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
        exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)]
        cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
        rec.update(minAlpha=min(mins) if mins else None, excess=max(exc) if exc else 0.0, ofCo=max(cos) if cos else None,
                   steps=len(re.findall(r"\nTime = ", log)))
        if "ABORT at step" in log and rec["status"] == "noCSV":
            rec["status"] = "aborted"
    launch = {}
    lp = os.path.join(c, ".leia_launch")
    if os.path.isfile(lp):
        for l in open(lp):
            if "=" in l:
                k, v = l.strip().split("=", 1); launch[k] = v
    try:
        rec["wall"] = float(launch["end"]) - float(launch["start"])
    except (KeyError, ValueError):
        rec["wall"] = None
    rec["np"] = int(launch.get("np", 0) or 0)
    rec["rc"] = launch.get("rc", "")
    sp = os.path.join(c, ".leia_stamp.json")
    s = json.load(open(sp)) if os.path.isfile(sp) else {}
    bins = s.get("trigger", {}).get("binaries", {})
    rec["stampId"] = s.get("id", "")
    rec["stampBinary"] = (next(iter(bins.values()), {}).get("sha256") or "")[:12] if bins else ""
    rec["stampLeia"] = s.get("record", {}).get("leiaAtRun") or s.get("record", {}).get("leia", "")
    return rec


def collect(studies_dir, family):
    recs = []
    for st in sorted(glob.glob(os.path.join(studies_dir, f"standing_{family}_N*"))):
        for c in sorted(glob.glob(os.path.join(st, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = read_case(c)
            if r:
                r["study"] = os.path.basename(st); recs.append(r)
    return recs


def cells(recs):
    groups = {}
    for r in recs:
        groups.setdefault((r["row"], r["N"], r["co"]), []).append(r)
    out = []
    for k in sorted(groups, key=lambda k: (ROW_ORDER.get(k[0], 9), k[1], k[2])):
        rs = sorted(groups[k], key=lambda r: r["rep"]); r0 = rs[0]
        same = all(all(r[c] == r0[c] for c in RESULT_COLS) for r in rs[1:])
        differing = sorted({c for r in rs[1:] for c in RESULT_COLS if r[c] != r0[c]})
        walls = [r["wall"] for r in rs if r["wall"] is not None]
        steps = [r["stepMedian"] for r in rs if r["stepMedian"] is not None]
        cell = dict(row=k[0], N=k[1], coTarget=k[2], ofCo=r0["ofCo"], np=r0["np"], dt=r0["dt"], T=r0["T"], reps=len(rs),
                    repsIdentical="identical" if same else "DIFFER: " + ",".join(differing),
                    **{c: r0[c] for c in RESULT_COLS}, tLast=r0["tLast"],
                    wallMedian=statistics.median(walls) if walls else None, wallMin=min(walls) if walls else None, wallMax=max(walls) if walls else None,
                    stepMedian=statistics.median(steps) if steps else None,
                    grazingCut=r0["grazingCut"], stampIds=";".join(sorted({r["stampId"] for r in rs if r["stampId"]})),
                    stampBinary=";".join(sorted({r["stampBinary"] for r in rs if r["stampBinary"]})), stampLeia=";".join(sorted({r["stampLeia"] for r in rs if r["stampLeia"]})),
                    study=r0["study"])
        out.append(cell)
    # the wall-clock and error ratios DAVOF / plicRDF at the same (N, Co)
    ref = {(c["N"], c["coTarget"]): c for c in out if c["row"] == "plicRDF"}
    for c in out:
        o = ref.get((c["N"], c["coTarget"]))
        for key, name in (("wallMedian", "wallRatioToRDF"), ("stepMedian", "stepRatioToRDF"), ("Eshape", "EshapeRatioToRDF"), ("L1n", "L1nRatioToRDF")):
            c[name] = (c[key] / o[key]) if (o and c[key] is not None and o.get(key)) else None
    return out


ROW_ORDER = {"davof-planeSurf-ringArea": 0, "davof-planeSurf-ringArea [recommended]": 1,
             "davof-planeSurf-ringArea [signedTets]": 2, "davof-planeSurf-ringArea [closureThenTets]": 3,
             "davof-planeSurf-ringArea [onePlane]": 4, "plicRDF": 9}


def orders(cells_):
    """The least-squares order over the completed rungs per (row, Co) and the pairwise ones."""
    out = []
    by = {}
    for c in cells_:
        by.setdefault((c["row"], c["coTarget"]), []).append(c)
    for (row, co), cs in sorted(by.items(), key=lambda kv: (ROW_ORDER.get(kv[0][0], 9), kv[0][1])):
        cs = sorted(cs, key=lambda c: c["N"])
        rec = dict(row=row, coTarget=co, rungs=",".join(str(c["N"]) for c in cs))
        for m in ("Eshape", "L1x", "L1n"):
            pts = [(c["N"], c[m]) for c in cs if c["status"] == "ended" and c[m] is not None and c[m] > 0]
            if len(pts) >= 2:
                xs = [math.log(1.0 / n) for n, _ in pts]; ys = [math.log(e) for _, e in pts]
                xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
                den = sum((x - xm) ** 2 for x in xs)
                rec[m + "_order"] = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / den if den else None
                rec[m + "_pairwise"] = ";".join("%.2f" % (math.log(pts[i][1] / pts[i + 1][1]) / math.log(pts[i + 1][0] / pts[i][0])) for i in range(len(pts) - 1))
            else:
                rec[m + "_order"] = None; rec[m + "_pairwise"] = ""
            rec[m + "_completedRungs"] = ",".join(str(n) for n, _ in pts)
        out.append(rec)
    return out


def tex_escape(s):
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def write_tables(family, cells_, ords, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, f"standing_{family}")
    if cells_:
        with open(base + ".csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cells_[0].keys())); w.writeheader(); w.writerows(cells_)
    lines = [r"% " + f"standing_{family}: make_standing_tables.py; one row per (row, N, OF Co); the three repetitions' median wall clock",
             r"\begin{tabular}{llrrrrrrrrrrrrl}", r"\toprule",
             r"row & $N$ & np & Co & steps & $E_\mathrm{shape}$ & $L_1(x)$ & $L_1(n)$ & $L_1(\alpha_f)$ & mass \% & $\min\alpha$ & excess & wall s & s/step & status \\", r"\midrule"]
    last_row = None
    for c in cells_:
        if last_row is not None and c["row"] != last_row:
            lines.append(r"\midrule")
        last_row = c["row"]
        status = c["status"] + (f" ({c['abortStep']})" if c["status"] == "aborted" and c["abortStep"] else "")
        lines.append(" & ".join([tex_escape(LABELS.get(c["row"], c["row"])), str(c["N"]), str(c["np"]), f(c["ofCo"], 2), str(c["steps"] or "-"),
                                 f(c["Eshape"], 3), f(c["L1x"], 3), f(c["L1n"], 3), f(c["L1af"], 3) if c["L1af"] is not None and c["L1af"] >= 0 else "-",
                                 f(c["massDevPct"], 1), f(c["minAlpha"], 1), f(c["excess"], 1),
                                 ("%.0f" % c["wallMedian"]) if c["wallMedian"] is not None else "-",
                                 f(c["stepMedian"], 2), tex_escape(status) + ("" if c["repsIdentical"] == "identical" else r" \textbf{reps differ}")]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(base + ".tex", "w").write("\n".join(lines) + "\n")
    if ords:
        with open(base + "_orders.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ords[0].keys())); w.writeheader(); w.writerows(ords)
    ol = [r"\begin{tabular}{llllll}", r"\toprule", r"row & Co & rungs (completed) & order $E_\mathrm{shape}$ & order $L_1(x)$ & order $L_1(n)$ \\", r"\midrule"]
    for o in ords:
        ol.append(" & ".join([tex_escape(LABELS.get(o["row"], o["row"])), str(o["coTarget"]), tex_escape(o["Eshape_completedRungs"] or "-")]
                             + [("%.2f (%s)" % (o[m + "_order"], o[m + "_pairwise"])) if o[m + "_order"] is not None else "-" for m in ("Eshape", "L1x", "L1n")]) + r" \\")
    ol += [r"\bottomrule", r"\end{tabular}"]
    open(base + "_orders.tex", "w").write("\n".join(ol) + "\n")
    # the caption's provenance
    rungs = {}
    for c in cells_:
        rungs.setdefault(c["N"], c)
    prov = ["The standing matrix, family " + tex_escape(family) + ": leia " + tex_escape(";".join(sorted({c["stampLeia"] for c in cells_ if c["stampLeia"]})) or "-")
            + ", the solver binary " + tex_escape(";".join(sorted({c["stampBinary"] for c in cells_ if c["stampBinary"]})) or "-")
            + "; the DAVOF rows with grazingCut " + tex_escape(";".join(sorted({c["grazingCut"] for c in cells_ if c["row"].startswith("davof")})) or "-") + "."]
    for N, c in sorted(rungs.items()):
        prov.append(f"$N = {N}$: {c['np']} rank(s), $T = {c['T']}$, " + "; ".join(f"Co {x['coTarget']}: dt {x['dt']}, {x['steps'] or '-'} steps" for x in cells_ if x["N"] == N and x["row"] == "plicRDF") + ".")
    prov.append("Stamp ids: " + tex_escape(", ".join(sorted({sid for c in cells_ for sid in c["stampIds"].split(";") if sid}))) + ".")
    open(base + "_caption.tex", "w").write("\n".join(prov) + "\n")
    print("written:", base + ".{csv,tex}", base + "_orders.{csv,tex}", base + "_caption.tex", f"({len(cells_)} cells)")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("family", choices=("shear2D", "shear3D", "deformation3D"))
    ap.add_argument("--studies-dir", default=os.path.join(REPO, "studies", "davof"))
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "davof", "data", "standing"))
    a = ap.parse_args(argv)
    recs = collect(a.studies_dir, a.family)
    if not recs:
        print("no cases under", a.studies_dir, "for", a.family); return 1
    cs = cells(recs)
    write_tables(a.family, cs, orders(cs), a.out)
    print("| row | N | np | OF Co | steps | Eshape | L1x | L1n | L1af | mass % | min alpha | excess | wall med [min,max] | s/step | wall/rdf | status | reps |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in cs:
        print(f"| {c['row']} | {c['N']} | {c['np']} | {f(c['ofCo'], 3)} | {c['steps']} | {f(c['Eshape'])} | {f(c['L1x'])} | {f(c['L1n'])} | {f(c['L1af'])} | {f(c['massDevPct'], 2)} | {f(c['minAlpha'], 2)} | {f(c['excess'], 2)} | "
              f"{f(c['wallMedian'], 1)} [{f(c['wallMin'], 1)}, {f(c['wallMax'], 1)}] | {f(c['stepMedian'], 3)} | {f(c['wallRatioToRDF'], 2)} | {c['status']} {c['abortStep']} | {c['reps']} {c['repsIdentical']} |")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
