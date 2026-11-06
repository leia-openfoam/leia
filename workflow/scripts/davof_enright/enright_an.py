#!/usr/bin/env python3
"""The scoring of the timed 3D deformation comparison (2026-11-06; the
pre-registration E-1 .. E-5 in config/davof/sphereDeformation3D_hexSweptN32.yaml).

Reads /work/scratch/tm83tomy/davof/agents/enright/time/N<N>_co<CO>/{ring,rdf}_r<rep>/
(log.advectorVoF, the converter's advectorVoF.csv, the written alpha fields, the
timing*.csv and sacct_steps*.txt of enright-time.sbatch) and writes under
.../enright/csv/:
  enright_runs.csv      every run: the end errors, the bounds and the mass over the
                        whole run, the counts at T/2 and T, the swept counters, the
                        timings
  enright_table.csv     per (N, Co, method): the errors of repetition 1 (the
                        repetitions are checked to be identical), the median / min /
                        max of the timings, the ratios to plicRDF
  enright_orders.csv    the orders over N of the end errors per Co and method
  enright_counters.csv  ours: the swept counters at the four write times and their
                        maxima over the run
The case: cases/davof/sphereDeformation3D (LeVeque 1996, doi:10.1137/0733033;
Enright et al. 2002, doi:10.1006/jcph.2002.7166).
"""
import csv, glob, os, re, statistics, sys, math

W = os.environ.get("ENRIGHT_W", "/work/scratch/tm83tomy/davof/agents/enright")
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
V_SPHERE = 4.0/3.0*math.pi*0.15**3
T_END, T_HALF = 3.0, 1.5
# ring = the swept-face advection with the ring normal (ours), rdf = plicRDF + isoAdvector,
# prism = the plane-clipped loop with the ring normal and ccs (the third row, 2026-11-06)
METHODS = ("ring", "rdf", "prism")
DAVOF = ("ring", "prism")

_RE_T = re.compile(r"^t = ([^,]+),\s*sum\(alpha\*V\) = ([^,]+),\s*dev = ([^%]+)%,\s*1-max\(alpha1\) = (\S+) at cell \S+,\s*min\(alpha1\) = (\S+) at cell", re.M)
_RE_BOUND = re.compile(r"After conservative bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)")
_RE_CO = re.compile(r"Courant Number mean: \S+ max: (\S+)")
_RE_CONT = re.compile(r"Continuity error: max\(mag\(sum\(phi\)\)\) = (\S+)")
_RE_DAVOF = re.compile(r"^davof: (\d+) interface cells", re.M)
_RE_ABORT = re.compile(r"advectorVoF: ABORT at step (\d+), t = (\S+): (.*)")
# the swept counters of the curved-clipping line (fork c5c00c5)
_RE_SWEPT = {
    "rankFallbacks": r"rank fallbacks (\d+)",
    "beyondLayer": r"beyond layer (\d+) \(max (\S+)\)",
    "pyramidDeltaMaxOverH": r"max \|delta\|/h ([0-9.eE+-]+)",
    "pyramidVolDevMax": r"max \|V_f - V_g\|/V_f ([0-9.eE+-]+)",
    "twistedFaces": r"twisted faces (\d+)",
    "pyramidSkipped": r"pyramid skipped (\d+)",
    "wispPiecesNeighbour": r"pieces clipped by a neighbour's plane (\d+)",
    "wispPiecesNoSurface": r"wisps without a neighbour surface (\d+)",
    "wispMatched": r"wisp neighbour plane matched (\d+)",
    "wispAnchorNotBracketed": r"wisp anchor not bracketed (\d+)",
    "facesNoSurface": r"faces without a surface (\d+)",
    "pieces": r"cells: pieces (\d+)",
    "partialCells": r"partial cells (\d+)",
    "degenerate": r"degenerate (\d+)",
    "pieceClosureMax": r"piece closure max ([0-9.eE+-]+)",
    "sweptFaces": r"swept faces (\d+)",
    "outOfDonor": r"out of donor (\d+)",
    "anchorFailed": r"volume anchor: failed (\d+)",
    "ruleFallbacks": r"volume rule \w+ fallbacks (\d+)",
    "boundingSecondLayer": r"bounding extended by a second layer: (\d+) cells",
    "remoteCells": r"remote cells (\d+)",
    "eulerTraces": r"Euler traces (\d+)",
}
_RE_SWEPT = {k: re.compile(v) for k, v in _RE_SWEPT.items()}


def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


def read_log(path):
    """Per-step blocks of the log: the time, the mass line, the bounds, the counters."""
    txt = open(path, errors="replace").read()
    blocks = re.split(r"\nTime = ", txt)
    steps = []
    for b in blocks[1:]:
        t = f(b.split("\n", 1)[0].strip())
        d = {"t": t}
        m = _RE_T.search(b)
        if m:
            d["sumAlphaV"] = f(m.group(2)); d["devPct"] = f(m.group(3)); d["excess"] = -f(m.group(4)); d["minAlpha"] = f(m.group(5))
        m = _RE_BOUND.search(b)
        if m: d["minAlphaBound"] = f(m.group(1)); d["excessBound"] = f(m.group(2))
        m = _RE_CO.search(b)
        if m: d["ofCo"] = f(m.group(1))
        m = _RE_CONT.search(b)
        if m: d["contErr"] = f(m.group(1))
        m = _RE_DAVOF.search(b)
        if m: d["nInterfaceAdv"] = int(m.group(1))
        for k, r in _RE_SWEPT.items():
            m = r.search(b)
            if m:
                d[k] = f(m.group(1)) if k != "beyondLayer" else f(m.group(1))
                if k == "beyondLayer": d["beyondLayerMax"] = f(m.group(2))
        steps.append(d)
    abort = _RE_ABORT.search(txt)
    ended = "\nEnd\n" in txt or txt.rstrip().endswith("End")
    return steps, (abort.groups() if abort else None), ended


def read_alpha(case, tname):
    """The internal alpha.water of a time directory (serial or decomposed): the
    partial-cell count (1e-8 < alpha < 1 - 1e-8), min, max and sum over the cells."""
    files = sorted(glob.glob("%s/processor*/%s/alpha.water" % (case, tname))) or [case + "/%s/alpha.water" % tname]
    n_part = 0; amin = float("inf"); amax = -float("inf"); n = 0
    for fn in files:
        if not os.path.isfile(fn): return None
        txt = open(fn).read()
        m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
        if not m:
            m2 = re.search(r"internalField\s+uniform\s+(\S+);", txt)
            if m2: v = float(m2.group(1)); amin = min(amin, v); amax = max(amax, v); continue
            return None
        start = m.end(); end = txt.index("\n)", start)
        vals = [float(x) for x in txt[start:end].split()]
        n += len(vals)
        n_part += sum(1 for v in vals if 1e-8 < v < 1 - 1e-8)
        amin = min(amin, min(vals)); amax = max(amax, max(vals))
    return dict(nPartial=n_part, minAlpha=amin, maxAlpha=amax, nCells=n)


def time_dirs(case):
    ds = glob.glob(case + "/processor0/[0-9]*") or glob.glob(case + "/[0-9]*")
    out = {}
    for d in ds:
        b = os.path.basename(d)
        try: out[float(b)] = b
        except ValueError: pass
    return out


def nearest(d, t, tol=1e-6):
    for k, v in d.items():
        if abs(k - t) <= tol: return v
    return None


def conv_rows(case):
    fn = case + "/advectorVoF.csv"
    return list(csv.DictReader(open(fn))) if os.path.isfile(fn) else []


def conv_at(rows, t, key, tol=1e-6):
    best = None
    for r in rows:
        if r.get(key) in (None, ""): continue
        tt = f(r.get("TIME"))
        if abs(tt - t) <= tol: best = r[key]
    return f(best)


def main():
    runs = []
    for d in sorted(glob.glob(W + "/time/N*_co*")):
        tfs = sorted(glob.glob(d + "/timing*.csv"))
        if not tfs: continue
        sac = {}
        for sf in glob.glob(d + "/sacct_steps*.txt"):
            for l in open(sf):
                p = l.strip().split("|")
                if len(p) >= 4 and "_r" in p[1]: sac[p[1]] = int(p[3])
        for t in (row for tf in tfs for row in csv.DictReader(open(tf))):
            c = "%s/%s_r%s" % (d, t["method"], t["rep"])
            if not os.path.isfile(c + "/log.advectorVoF"): continue
            steps, abort, ended = read_log(c + "/log.advectorVoF")
            R = conv_rows(c)
            tds = time_dirs(c)
            aT = read_alpha(c, nearest(tds, T_END)) if nearest(tds, T_END) else None
            aH = read_alpha(c, nearest(tds, T_HALF)) if nearest(tds, T_HALF) else None
            g = lambda k: [s[k] for s in steps if k in s]
            mins = g("minAlpha") + g("minAlphaBound"); exc = g("excess") + g("excessBound")
            sH = min(steps, key=lambda s: abs(s["t"] - T_HALF)) if steps else {}
            sT = steps[-1] if steps else {}
            run = dict(N=int(t["N"]), np=int(t["np"]), co_index=int(t["co_index"]), method=t["method"], rep=int(t["rep"]),
                       status="aborted" if abort else ("ended" if ended else "incomplete"), abortStep=abort[0] if abort else "", abortReason=abort[2] if abort else "",
                       steps=len(steps), dt=f(t["dt"]), end_time=f(t["end_time"]), tLast=sT.get("t", float("nan")),
                       ofCoMax=max(g("ofCo")) if g("ofCo") else float("nan"), contErrMax=max(g("contErr")) if g("contErr") else float("nan"),
                       minAlphaRun=min(mins) if mins else float("nan"), excessRun=max(exc) if exc else float("nan"),
                       massDevPctMax=max(abs(x) for x in g("devPct")) if g("devPct") else float("nan"),
                       Eshape=conv_at(R, T_END, "Eshape"), EshapeRel=conv_at(R, T_END, "Eshape")/V_SPHERE,
                       L1x=conv_at(R, T_END, "LCentreDiff1"), L1n=conv_at(R, T_END, "LNormalDiff1"),
                       L1af=conv_at(R, T_END, "L1_alphaf") if t["method"] in DAVOF else float("nan"),
                       EshapeHalf=conv_at(R, T_HALF, "Eshape"),
                       nInterfaceT=conv_at(R, T_END, "nInterface"), nInterfaceHalf=conv_at(R, T_HALF, "nInterface"),
                       nInterfaceAdvT=sT.get("nInterfaceAdv", float("nan")), nInterfaceAdvHalf=sH.get("nInterfaceAdv", float("nan")),
                       nInterfaceAdvMax=max(g("nInterfaceAdv")) if g("nInterfaceAdv") else float("nan"),
                       nPartialT=aT["nPartial"] if aT else float("nan"), nPartialHalf=aH["nPartial"] if aH else float("nan"),
                       minAlphaFieldT=aT["minAlpha"] if aT else float("nan"), maxAlphaFieldT=aT["maxAlpha"] if aT else float("nan"),
                       rankFallbacksSum=sum(g("rankFallbacks")), beyondLayerSum=sum(g("beyondLayer")), beyondLayerMax=max(g("beyondLayerMax")) if g("beyondLayerMax") else float("nan"),
                       pyramidDeltaMaxOverH=max(g("pyramidDeltaMaxOverH")) if g("pyramidDeltaMaxOverH") else float("nan"),
                       twistedFacesSum=sum(g("twistedFaces")), wispPiecesNoSurfaceMax=max(g("wispPiecesNoSurface")) if g("wispPiecesNoSurface") else float("nan"),
                       wispPiecesNeighbourMax=max(g("wispPiecesNeighbour")) if g("wispPiecesNeighbour") else float("nan"),
                       boundingSecondLayerSum=sum(g("boundingSecondLayer")), boundingSecondLayerMax=max(g("boundingSecondLayer")) if g("boundingSecondLayer") else float("nan"),
                       anchorFailedSum=sum(g("anchorFailed")), ruleFallbacksSum=sum(g("ruleFallbacks")), degenerateMax=max(g("degenerate")) if g("degenerate") else float("nan"),
                       clock_s=f(t["clock_s"]), exec_s=f(t["exec_s"]), wall_s=f(t["wall_s"]), sacct_s=sac.get("%s_r%s" % (t["method"], t["rep"]), float("nan")),
                       median_step_s=f(t["median_step_s"]), node=t["node"], case=c)
            # the wisps: the partial cells (1e-8 < alpha < 1 - 1e-8 in the written field)
            # without a reconstructed surface (reconstructionError's nInterface at the
            # same write; for the davof rows the cells with |m_c| > wispTol V^(2/3), for
            # plicRDF the cells inside (surfCellTol, 1 - surfCellTol), so plicRDF's count
            # is zero by definition); davof's own per-step count nInterfaceAdv is kept
            # next to it (it counts by the face fractions, alpha unrestricted)
            run["nWispsT"] = run["nPartialT"] - run["nInterfaceT"]
            run["nWispsHalf"] = run["nPartialHalf"] - run["nInterfaceHalf"]
            runs.append(run)
            if t["method"] in DAVOF and int(t["rep"]) == 1:
                with open(OUT + "/enright_counters_N%s_co%s_%s.csv" % (t["N"], t["co_index"], t["method"]), "w", newline="") as fh:
                    keys = ["t", "nInterfaceAdv", "minAlpha", "excess", "devPct", "contErr"] + list(_RE_SWEPT.keys()) + ["beyondLayerMax"]
                    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(steps)
    if not runs:
        print("no runs"); return
    with open(OUT + "/enright_runs.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)
    groups = {}
    for r in runs: groups.setdefault((r["N"], r["co_index"], r["method"]), []).append(r)
    table = []
    for k in sorted(groups):
        rs = sorted(groups[k], key=lambda r: r["rep"]); r0 = rs[0]
        same = all(abs(r["Eshape"] - r0["Eshape"]) <= 1e-12*abs(r0["Eshape"]) + 1e-300 for r in rs if r["Eshape"] == r["Eshape"])
        med = lambda x: statistics.median([r[x] for r in rs])
        row = dict(N=k[0], np=r0["np"], co_index=k[1], ofCoMax=r0["ofCoMax"], method=k[2], reps=len(rs), repsIdentical=same, status=",".join(sorted(set(r["status"] for r in rs))),
                   steps=r0["steps"], dt=r0["dt"], tLast=r0["tLast"], contErrMax=r0["contErrMax"],
                   Eshape=r0["Eshape"], EshapeRel=r0["EshapeRel"], L1x=r0["L1x"], L1n=r0["L1n"], L1af=r0["L1af"], EshapeHalf=r0["EshapeHalf"],
                   minAlphaRun=r0["minAlphaRun"], excessRun=r0["excessRun"], massDevPctMax=r0["massDevPctMax"],
                   nInterfaceHalf=r0["nInterfaceHalf"], nInterfaceT=r0["nInterfaceT"], nInterfaceAdvHalf=r0["nInterfaceAdvHalf"], nInterfaceAdvT=r0["nInterfaceAdvT"],
                   nPartialHalf=r0["nPartialHalf"], nPartialT=r0["nPartialT"], nWispsHalf=r0["nWispsHalf"], nWispsT=r0["nWispsT"],
                   rankFallbacksSum=r0["rankFallbacksSum"], beyondLayerSum=r0["beyondLayerSum"], pyramidDeltaMaxOverH=r0["pyramidDeltaMaxOverH"], twistedFacesSum=r0["twistedFacesSum"],
                   boundingSecondLayerSum=r0["boundingSecondLayerSum"], anchorFailedSum=r0["anchorFailedSum"], ruleFallbacksSum=r0["ruleFallbacksSum"],
                   wall_med=med("wall_s"), wall_min=min(r["wall_s"] for r in rs), wall_max=max(r["wall_s"] for r in rs),
                   clock_med=med("clock_s"), exec_med=med("exec_s"), sacct_med=med("sacct_s"), step_med=med("median_step_s"),
                   node=",".join(sorted(set(r["node"] for r in rs))))
        table.append(row)
    for row in table:
        rr = [x for x in table if x["N"] == row["N"] and x["co_index"] == row["co_index"] and x["method"] == "rdf"]
        for key in ("wall_med", "Eshape", "L1x", "L1n"):
            row[key + "_ratio_to_rdf"] = row[key]/rr[0][key] if rr and rr[0][key] else float("nan")
    with open(OUT + "/enright_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    orders = []
    for co in sorted(set(r["co_index"] for r in table)):
        for m in METHODS:
            rows = sorted([r for r in table if r["co_index"] == co and r["method"] == m], key=lambda r: r["N"])
            for a, b in zip(rows, rows[1:]):
                o = dict(co_index=co, method=m, N1=a["N"], N2=b["N"])
                for key in ("Eshape", "L1x", "L1n", "L1af"):
                    o["order_" + key] = math.log2(a[key]/b[key]) if a[key] > 0 and b[key] > 0 else float("nan")
                orders.append(o)
    with open(OUT + "/enright_orders.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(orders[0].keys()) if orders else ["co_index"]); w.writeheader(); w.writerows(orders)
    print("N np co OFCo method status steps | Eshape Eshape/Vs L1x L1n L1af | minAlpha excess mass%max | nInt(T/2) nInt(T) nPart(T/2) nPart(T) wisps(T/2) wisps(T) | wall med [min,max] s, step med | wall/rdf E/rdf L1x/rdf | reps id | rankFb beyond pyrDelta twisted bound2")
    for r in table:
        print("%d %d %d %.3f %s %s %d | %.3e %.4f %.3e %.3e %s | %.2e %.2e %.1e | %s %s %s %s %s %s | %.1f [%.1f, %.1f] %.3f | %.2f %.3f %.3f | %d %s | %d %d %s %d %d" % (
            r["N"], r["np"], r["co_index"], r["ofCoMax"], r["method"], r["status"], r["steps"], r["Eshape"], r["EshapeRel"], r["L1x"], r["L1n"],
            "%.3e" % r["L1af"] if r["L1af"] == r["L1af"] else "-", r["minAlphaRun"], r["excessRun"], r["massDevPctMax"],
            r["nInterfaceHalf"], r["nInterfaceT"], r["nPartialHalf"], r["nPartialT"], r["nWispsHalf"], r["nWispsT"],
            r["wall_med"], r["wall_min"], r["wall_max"], r["step_med"], r["wall_med_ratio_to_rdf"], r["Eshape_ratio_to_rdf"], r["L1x_ratio_to_rdf"], r["reps"], r["repsIdentical"],
            int(r["rankFallbacksSum"]), int(r["beyondLayerSum"]), "%.2e" % r["pyramidDeltaMaxOverH"] if r["pyramidDeltaMaxOverH"] == r["pyramidDeltaMaxOverH"] else "-", int(r["twistedFacesSum"]), int(r["boundingSecondLayerSum"])))
    print("orders:")
    for o in orders:
        print("  co %d %s %d -> %d: Eshape %.2f L1x %.2f L1n %.2f L1af %s" % (o["co_index"], o["method"], o["N1"], o["N2"], o["order_Eshape"], o["order_L1x"], o["order_L1n"],
              "%.2f" % o["order_L1af"] if o["order_L1af"] == o["order_L1af"] else "-"))


if __name__ == "__main__":
    main()
