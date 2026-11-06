#!/usr/bin/env python3
"""The scoring of the timed verification of the isoAdvector package (2026-11-06; the
pre-registration F-1 .. F-5 in leia config/davof/iso*3D_hexN32.yaml).

Reads <W>/time/<FLOW>_N<N>_co<CO>/<method>_r<rep>/ (log.advectorVoF, the converter's
advectorVoF.csv, timing*.csv of iso-time.sbatch; the profiled run <method>_rprof with
uniform/profiling of the last time) and writes under <W>/csv/:
  iso_runs.csv      every run: the end errors (E_shape, L1(x), L1(n), L1(alphaf)), the
                    bounds and the mass over the whole run, the interface-cell count at the
                    writes, the status (completed | diverged | crashed), the timing
  iso_table.csv     per (flow, N, Co, method): the errors of repetition 1 (the repetitions
                    are checked to be identical), the median / min / max of the wall clock
                    and of the per-step time, the ratios to RDF, the profiling shares
  iso_orders.csv    the orders over N of the end errors per flow, Co and method
  iso_verdicts.txt  F-1 .. F-5 against the measurements
The methods: RDF (plicRDF, vofCutCell), RDFccs (plicRDF, ccs), RING (davof ring normal,
ccs), RINGvcc (davof ring normal, vofCutCell), ISO (isoAlpha).
"""
import csv, glob, os, re, statistics, sys, math

W = os.environ.get("ISO_W", "/work/scratch/tm83tomy/davof/agents/swept/iso")
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
METHODS = ("RDF", "RDFccs", "RING", "RINGvcc", "ISO")
DAVOF = ("RING", "RINGvcc")
T_END = {"Rotation": 1.0, "Shear": 3.0, "Deformation": 3.0}
CO_LABEL = {"0": 0.2, "1": 0.5, "2": 1.0}

_RE_T = re.compile(r"^t = ([^,]+),\s*sum\(alpha\*V\) = ([^,]+),\s*dev = ([^%]+)%,\s*1-max\(alpha1\) = (\S+) at cell \S+,\s*min\(alpha1\) = (\S+) at cell", re.M)
_RE_BOUND = re.compile(r"After conservative bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)")
_RE_CO = re.compile(r"Courant Number mean: \S+ max: (\S+)")
_RE_ABORT = re.compile(r"advectorVoF: ABORT at step (\d+), t = (\S+): (.*)")
_RE_STEP = re.compile(r"timeStepTime = (\S+) s")
_RE_DAVOF = re.compile(r"^davof: (\d+) interface cells", re.M)
_RE_AF = re.compile(r"alphaf set by a plane on (\d+) faces \((\d+) in the second ring; the plane's flux on (\d+) faces beyond[^)]*\), (\d+) fallbacks, (\d+) clips(?:; content bound on (\d+) faces of (\d+) cells)?")
_RE_CCS = re.compile(r"positioning ccs: (\d+) positionings, cuts per cell (\S+) \(max (\d+)[^)]*\), residual max (\S+) mean (\S+)")

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

def read_csv(p):
    with open(p) as f:
        return list(csv.DictReader(f))

def last_row(rows, key, t_end):
    """The row at the end time (the last non-blank value of key)."""
    best = None
    for r in rows:
        if r.get(key, "") not in ("", None) and r.get("TIME", "") != "":
            best = r
    return best

def profile(c):
    """OpenFOAM profiling of the last time: {name: (calls, total_s, child_s)}."""
    ts = sorted((d for d in glob.glob(c + "/[0-9]*") if os.path.isdir(d) and os.path.exists(d + "/uniform/profiling")),
                key=lambda d: float(os.path.basename(d)))
    if not ts:
        ts = sorted((d for d in glob.glob(c + "/processor0/[0-9]*") if os.path.exists(d + "/uniform/profiling")),
                    key=lambda d: float(os.path.basename(d)))
    if not ts: return {}
    s = open(ts[-1] + "/uniform/profiling").read()
    out = {}
    for m in re.finditer(r"\{\s*id\s+\d+;\s*description\s+([^;]+);\s*calls\s+(\d+);\s*totalTime\s+(\S+);\s*childTime\s+(\S+);", s):
        out.setdefault(m.group(1).strip().strip('"'), [0, 0.0, 0.0])
        o = out[m.group(1).strip().strip('"')]
        o[0] += int(m.group(2)); o[1] += float(m.group(3)); o[2] += float(m.group(4))
    return out

def scan_run(c, flow, N, co, m, rep):
    rec = dict(flow=flow, N=N, co_index=co, co=CO_LABEL.get(co, co), method=m, rep=rep, case=c)
    logp = c + "/log.advectorVoF"
    if not os.path.exists(logp):
        rec["status"] = "missing"; return rec
    log = open(logp).read()
    tl = _RE_T.findall(log)
    rec["steps"] = len(re.findall(r"\nTime = ", log))
    rec["co_max_measured"] = max((fnum(x) for x in _RE_CO.findall(log)), default=float("nan"))
    if tl:
        s0 = fnum(tl[0][1])
        rec["mass_dev_max_pct"] = max(abs(fnum(x[2])) for x in tl)
        rec["min_alpha_run"] = min(fnum(x[4]) for x in tl)
        rec["max_excess_run"] = max(-fnum(x[3]) for x in tl)   # 1 - max(alpha) negative = excess
    b = _RE_BOUND.findall(log)
    if b:
        rec["min_alpha_bounding"] = min(fnum(x[0]) for x in b)
        rec["max_excess_bounding"] = max(fnum(x[1]) for x in b)
    ab = _RE_ABORT.search(log)
    rec["status"] = "diverged" if ab else ("completed" if re.search(r"^End\s*$", log, re.M) else "incomplete")
    if ab:
        rec["abort_step"] = int(ab.group(1)); rec["abort_t"] = fnum(ab.group(2)); rec["abort_reason"] = ab.group(3).strip()
    st = [fnum(x) for x in _RE_STEP.findall(log)][1:-1]
    if st:
        rec["step_s_median"] = statistics.median(st); rec["step_s_mean"] = statistics.mean(st)
    ic = _RE_DAVOF.findall(log)
    if ic:
        rec["interface_cells_first"] = int(ic[0]); rec["interface_cells_last"] = int(ic[-1]); rec["interface_cells_max"] = max(int(x) for x in ic)
    af = _RE_AF.findall(log)
    if af:
        rec["af_plane_faces_last"] = int(af[-1][0]); rec["af_second_ring_last"] = int(af[-1][1]); rec["af_fallbacks_max"] = max(int(x[3]) for x in af)
        rec["content_bound_faces_max"] = max(int(x[5] or 0) for x in af)
    cc = _RE_CCS.findall(log)
    if cc:
        rec["ccs_cuts_per_cell"] = statistics.mean(fnum(x[1]) for x in cc); rec["ccs_residual_max"] = max(fnum(x[3]) for x in cc)
    # the converter's end values
    csvp = c + "/advectorVoF.csv"
    if os.path.exists(csvp):
        rows = read_csv(csvp)
        for key, name in (("Eshape", "Eshape_end"), ("LNormalDiff1", "L1n_end"), ("LCentreDiff1", "L1x_end"),
                          ("L1_alphaf", "L1af_end"), ("Emass", "Emass_end"), ("Ebound", "Ebound_end"), ("nInterface", "nInterface_end")):
            r = last_row(rows, key, T_END[flow])
            rec[name] = fnum(r[key]) if r else float("nan")
            if r: rec[name + "_t"] = fnum(r["TIME"])
    pr = profile(c)
    if pr:
        tot = sum(v[1] for k, v in pr.items() if k.startswith("application"))  # the run's total
        for k, v in pr.items():
            kk = k.replace("::", "_").replace(" ", "_")
            if any(s in k for s in ("reconstruct", "advect", "application", "isoAdvection", "davof", "plicRDF", "isoAlpha")):
                rec["prof_" + kk + "_s"] = v[1]
                rec["prof_" + kk + "_calls"] = v[0]
    return rec

def main():
    runs = []
    for d in sorted(glob.glob(W + "/time/*_N*_co*")):
        mm = re.match(r"(\w+?)_(hex|perturbed|poly|tet)_N(\d+)_co(\d)$", os.path.basename(d))
        if mm:
            flow, mc, N, co = mm.group(1), mm.group(2), int(mm.group(3)), mm.group(4)
        else:
            mm = re.match(r"(\w+)_N(\d+)_co(\d)$", os.path.basename(d))
            if not mm: continue
            flow, mc, N, co = mm.group(1), "hex", int(mm.group(2)), mm.group(3)
        tim = {}
        for tp in glob.glob(d + "/timing*.csv"):
            for r in read_csv(tp):
                tim[(r["method"], r["rep"])] = r
        for c in sorted(glob.glob(d + "/*_r*")):
            m, rep = os.path.basename(c).rsplit("_r", 1)
            if m not in METHODS: continue
            rec = scan_run(c, flow, N, co, m, rep); rec["mesh"] = mc
            t = tim.get((m, rep))
            if t:
                rec["wall_s"] = fnum(t["wall_s"]); rec["clock_s"] = fnum(t["clock_s"]); rec["np"] = int(t["np"]); rec["node"] = t["node"]
            runs.append(rec)
    if not runs:
        print("no runs under", W + "/time"); return
    keys = sorted({k for r in runs for k in r})
    with open(OUT + "/iso_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(runs)
    # the table per cell and method
    table = []
    cells = sorted({(r["flow"], r["mesh"], r["N"], r["co_index"]) for r in runs})
    for (flow, mc, N, co) in cells:
        for m in METHODS:
            rs = [r for r in runs if (r["flow"], r["mesh"], r["N"], r["co_index"], r["method"]) == (flow, mc, N, co, m) and r["rep"] != "prof"]
            if not rs: continue
            pr = [r for r in runs if (r["flow"], r["mesh"], r["N"], r["co_index"], r["method"]) == (flow, mc, N, co, m) and r["rep"] == "prof"]
            r1 = sorted(rs, key=lambda r: r["rep"])[0]
            row = dict(flow=flow, mesh=mc, N=N, co=CO_LABEL.get(co, co), method=m, reps=len(rs), status=r1.get("status"),
                       steps=r1.get("steps"), np=r1.get("np"), co_max_measured=r1.get("co_max_measured"))
            for k in ("Eshape_end", "L1x_end", "L1n_end", "L1af_end", "Emass_end", "min_alpha_run", "max_excess_run", "mass_dev_max_pct",
                      "interface_cells_first", "interface_cells_max", "interface_cells_last", "abort_step", "abort_reason", "ccs_cuts_per_cell", "ccs_residual_max"):
                row[k] = r1.get(k)
            # identical repetitions?
            e = [r.get("Eshape_end") for r in rs]
            row["reps_identical"] = all(x == e[0] for x in e)
            for q in ("wall_s", "step_s_median"):
                v = [r[q] for r in rs if q in r and not math.isnan(r[q])]
                if v:
                    row[q + "_median"] = statistics.median(v); row[q + "_min"] = min(v); row[q + "_max"] = max(v)
                    row[q + "_spread_pct"] = 100*(max(v) - min(v))/statistics.median(v)
            if pr:
                p = pr[0]
                for k, v in p.items():
                    if k.startswith("prof_"): row[k] = v
            table.append(row)
    # ratios to RDF
    for row in table:
        ref = next((t for t in table if (t["flow"], t["mesh"], t["N"], t["co"]) == (row["flow"], row["mesh"], row["N"], row["co"]) and t["method"] == "RDF"), None)
        if ref:
            for q in ("Eshape_end", "L1x_end", "L1n_end", "wall_s_median", "step_s_median_median"):
                a, b = row.get(q), ref.get(q)
                row[q + "_over_RDF"] = (a/b if (a not in (None, "") and b not in (None, "", 0) and not (isinstance(a, float) and math.isnan(a))) else "")
    keys = []
    for row in table:
        for k in row:
            if k not in keys: keys.append(k)
    with open(OUT + "/iso_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(table)
    # orders over N
    orders = []
    for (flow, mc) in sorted({(t["flow"], t["mesh"]) for t in table}):
        for co in sorted({t["co"] for t in table if t["flow"] == flow}):
            for m in METHODS:
                ts = sorted([t for t in table if (t["flow"], t["mesh"], t["co"], t["method"]) == (flow, mc, co, m) and t["status"] == "completed"], key=lambda t: t["N"])
                for a, b in zip(ts, ts[1:]):
                    o = dict(flow=flow, mesh=mc, co=co, method=m, N_from=a["N"], N_to=b["N"])
                    for q in ("Eshape_end", "L1x_end", "L1n_end", "L1af_end"):
                        try:
                            o["order_" + q] = math.log2(a[q]/b[q])
                        except (TypeError, ValueError, ZeroDivisionError):
                            o["order_" + q] = ""
                    orders.append(o)
    if orders:
        with open(OUT + "/iso_orders.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(orders[0].keys())); w.writeheader(); w.writerows(orders)
    # a readable summary
    lines = []
    for row in table:
        lines.append("%-12s %-9s N %3d Co %-4s %-8s %-10s steps %5s wall %8s s (%s reps, spread %s%%) step %s s E_shape %s L1n %s L1x %s minA %s excess %s mass %s%% ic %s/%s %s" % (
            row["flow"], row["mesh"], row["N"], row["co"], row["method"], row["status"], row["steps"], fmt(row.get("wall_s_median")), row["reps"], fmt(row.get("wall_s_spread_pct"), 1),
            fmt(row.get("step_s_median_median"), 4), fmt(row.get("Eshape_end"), 3), fmt(row.get("L1n_end"), 3), fmt(row.get("L1x_end"), 3),
            fmt(row.get("min_alpha_run"), 2), fmt(row.get("max_excess_run"), 2), fmt(row.get("mass_dev_max_pct"), 2), row.get("interface_cells_first"), row.get("interface_cells_max"),
            ("ABORT step %s: %s" % (row.get("abort_step"), row.get("abort_reason"))) if row.get("abort_step") else ""))
    open(OUT + "/iso_summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote", OUT + "/iso_runs.csv", OUT + "/iso_table.csv", OUT + "/iso_orders.csv")

def fmt(x, d=3):
    if x in (None, "") or (isinstance(x, float) and math.isnan(x)): return "-"
    if isinstance(x, float):
        return ("%." + str(d) + "e") % x if (abs(x) < 1e-2 or abs(x) >= 1e4) else ("%." + str(d) + "f") % x
    return str(x)

if __name__ == "__main__":
    main()
