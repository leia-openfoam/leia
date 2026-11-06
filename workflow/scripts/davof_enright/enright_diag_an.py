#!/usr/bin/env python3
"""The scoring of the fragmentation diagnosis on the 3D deformation (2026-11-06;
the predictions D-1 .. D-6 in REPORT-swept "The 3D Enright deformation"). Reads
<W>/<tag>/N<N>_co<i>/<run>/ (log.advectorVoF, advectorVoF.csv, the written alpha
fields, postProcessing/sweptFragmentLog/cells_*.csv) and writes <W>/csv/
  diag_runs.csv      per run: min alpha / excess over the run, E_shape at T, the
                     face-fraction interface count and the partial cells at the
                     four writes, the rank fallbacks and the second-layer bounding
                     at step 1 / over the run, the fragment-log class counts
  diag_frag_<run>.csv the first 200 fragmenting cells with their classes
  diag_history_<run>.csv per step: t, min alpha, excess, interface cells, rank
                     fallbacks, second-layer bounding, wisp pieces, new fragments
    python3 enright_diag_an.py [tag]      (default diag3)
"""
import csv, glob, os, re, sys, math

W = os.environ.get("ENRIGHT_W", "/work/scratch/tm83tomy/davof/agents/enright")
TAG = sys.argv[1] if len(sys.argv) > 1 else "diag3"
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
V_SPHERE = 4.0/3.0*math.pi*0.15**3
WRITES = (0.75, 1.5, 2.25, 3.0)


def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


def read_alpha(case, tname):
    files = sorted(glob.glob("%s/processor*/%s/alpha.water" % (case, tname))) or [case + "/%s/alpha.water" % tname]
    n_part = 0; amin = float("inf"); amax = -float("inf")
    for fn in files:
        if not os.path.isfile(fn): return None
        txt = open(fn).read()
        m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
        if not m: return None
        start = m.end(); end = txt.index("\n)", start)
        vals = [float(x) for x in txt[start:end].split()]
        n_part += sum(1 for v in vals if 1e-8 < v < 1 - 1e-8)
        amin = min(amin, min(vals)); amax = max(amax, max(vals))
    return dict(nPartial=n_part, minAlpha=amin, maxAlpha=amax)


def time_dirs(case):
    out = {}
    for d in glob.glob(case + "/processor0/[0-9]*") or glob.glob(case + "/[0-9]*"):
        try: out[float(os.path.basename(d))] = os.path.basename(d)
        except ValueError: pass
    return out


def nearest(d, t, tol=1e-6):
    for k, v in d.items():
        if abs(k - t) <= tol: return v
    return None


_RE = {
    "minAlpha": r"min\(alpha1\) = (\S+) at cell", "oneMinusMax": r"1-max\(alpha1\) = (\S+) at cell",
    "ifc": r"\ndavof: (\d+) interface cells", "rankFb": r"rank fallbacks (\d+)", "swept": r"swept faces (\d+)",
    "bound2": r"second layer: (\d+) cells", "wispNoSurf": r"wisps without a neighbour surface (\d+)",
    "newFrag": r"sweptFragmentLog: new [^:]*: (\d+) \(fallback-adjacent (\d+), second-layer bounded (\d+), wisp without a surface (\d+), none (\d+)\)",
    "reinit": r"re-initialised from the planes[^:]*: (\d+) faces between two planes[^,]*, (\d+) by one plane, (\d+) from alpha",
}


def main():
    rows = []
    for d in sorted(glob.glob("%s/%s/N*_co*" % (W, TAG))):
        cell = os.path.basename(d)
        for c in sorted(glob.glob(d + "/*")):
            L = c + "/log.advectorVoF"
            if not os.path.isfile(L): continue
            run = os.path.basename(c)
            txt = open(L).read()
            blocks = re.split(r"\nTime = ", txt)[1:]
            hist = []
            for b in blocks:
                h = {"t": f(b.split("\n", 1)[0])}
                for k, r in _RE.items():
                    m = re.search(r, b)
                    if not m: continue
                    if k == "newFrag":
                        h["newFrag"], h["newFb"], h["newB2"], h["newW"], h["newNone"] = [int(x) for x in m.groups()]
                    elif k == "reinit":
                        h["reinitTwo"], h["reinitOne"], h["reinitAlpha"] = [int(x) for x in m.groups()]
                    else:
                        h[k] = f(m.group(1))
                hist.append(h)
            with open("%s/diag_history_%s_%s.csv" % (OUT, cell, run), "w", newline="") as fh:
                keys = ["t", "minAlpha", "oneMinusMax", "ifc", "rankFb", "swept", "bound2", "wispNoSurf", "newFrag", "newFb", "newB2", "newW", "newNone", "reinitTwo", "reinitOne", "reinitAlpha"]
                w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(hist)
            g = lambda k: [h[k] for h in hist if k in h]
            ab = re.search(r"ABORT at step (\d+), t = (\S+): ([^(]*)", txt)
            R = list(csv.DictReader(open(c + "/advectorVoF.csv"))) if os.path.isfile(c + "/advectorVoF.csv") else []
            def conv(t, key):
                v = [r[key] for r in R if r.get(key) and abs(f(r["TIME"]) - t) < 1e-6]
                return f(v[-1]) if v else float("nan")
            tds = time_dirs(c)
            row = dict(cell=cell, run=run, steps=len(hist), status="aborted" if ab else ("ended" if "\nEnd" in txt else "incomplete"),
                       abort=("step %s t %s %s" % ab.groups()).strip() if ab else "",
                       minAlphaRun=min(g("minAlpha")) if g("minAlpha") else float("nan"), excessRun=max(-x for x in g("oneMinusMax")) if g("oneMinusMax") else float("nan"),
                       minAlphaStep20=min(g("minAlpha")[:20]) if g("minAlpha") else float("nan"), excessStep20=max(-x for x in g("oneMinusMax")[:20]) if g("oneMinusMax") else float("nan"),
                       rankFbStep1=g("rankFb")[0] if g("rankFb") else float("nan"), rankFbSum=sum(g("rankFb")), rankFbMax=max(g("rankFb")) if g("rankFb") else float("nan"),
                       bound2Step1=g("bound2")[0] if g("bound2") else float("nan"), bound2Sum=sum(g("bound2")), bound2Max=max(g("bound2")) if g("bound2") else float("nan"),
                       wispNoSurfMax=max(g("wispNoSurf")) if g("wispNoSurf") else float("nan"),
                       ifcStep1=g("ifc")[0] if g("ifc") else float("nan"), ifcMax=max(g("ifc")) if g("ifc") else float("nan"),
                       fragTotal=sum(g("newFrag")), fragFb=sum(g("newFb")), fragB2=sum(g("newB2")), fragW=sum(g("newW")), fragNone=sum(g("newNone")),
                       fragFirstStep=next((i + 1 for i, h in enumerate(hist) if h.get("newFrag", 0) > 0), float("nan")),
                       EshapeT=conv(3.0, "Eshape"), EshapeTRel=conv(3.0, "Eshape")/V_SPHERE, L1xT=conv(3.0, "LCentreDiff1"), L1nT=conv(3.0, "LNormalDiff1"), L1afT=conv(3.0, "L1_alphaf"))
            for t in WRITES:
                tn = nearest(tds, t)
                a = read_alpha(c, tn) if tn else None
                row["nPartial_%g" % t] = a["nPartial"] if a else float("nan")
                row["nInterface_%g" % t] = conv(t, "nInterface")
                ht = [h for h in hist if abs(h["t"] - t) < 1e-6]
                row["ifcAdv_%g" % t] = ht[0].get("ifc", float("nan")) if ht else float("nan")
                row["ratioIfcPartial_%g" % t] = (row["ifcAdv_%g" % t]/row["nPartial_%g" % t]) if a and a["nPartial"] else float("nan")
            # the fragment log
            frag = []
            for fn in sorted(glob.glob(c + "/postProcessing/sweptFragmentLog/cells_*.csv")):
                frag += list(csv.DictReader(open(fn)))
            if frag:
                with open("%s/diag_frag_%s_%s.csv" % (OUT, cell, run), "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(frag[0].keys())); w.writeheader(); w.writerows(frag)
                n = len(frag)
                row["frag200_n"] = n
                row["frag200_firstStep"] = min(int(r["step"]) for r in frag); row["frag200_lastStep"] = max(int(r["step"]) for r in frag)
                row["frag200_fallbackAdj"] = sum(1 for r in frag if r["fallbackPrev1"] == "1" or r["fallbackPrev2"] == "1")
                row["frag200_bound2"] = sum(1 for r in frag if r["bound2Prev1"] == "1" or r["bound2Prev2"] == "1")
                row["frag200_wisp"] = sum(1 for r in frag if r["wispPrev1"] == "1" or r["wispPrev2"] == "1")
                row["frag200_none"] = sum(1 for r in frag if not any(r[k] == "1" for k in ("fallbackPrev1", "fallbackPrev2", "bound2Prev1", "bound2Prev2", "wispPrev1", "wispPrev2")))
                row["frag200_dry"] = sum(1 for r in frag if f(r["alpha"]) < 1e-6)
                row["frag200_full"] = n - row["frag200_dry"]
                mm = sorted(f(r["mOverV23"]) for r in frag)
                row["frag200_mOverV23_med"] = mm[n//2]; row["frag200_mOverV23_max"] = mm[-1]
                afmax = sorted(max(f(r["af%d" % i]) for i in range(6)) for r in frag)
                row["frag200_afMax_med"] = afmax[n//2]
                row["frag200_nWetFaces_med"] = sorted(sum(1 for i in range(6) if f(r["af%d" % i]) > 1e-12 and f(r["af%d" % i]) < 1 - 1e-12) for r in frag)[n//2]
            rows.append(row)
    if not rows:
        print("no runs"); return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(OUT + "/diag_runs_%s.csv" % TAG, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print("cell run status steps | min excess (run) | min excess (first 20) | rankFb step1/max | bound2 step1/max | ifc step1/max | frag total: fb b2 wisp none, first step | frag200: n firstStep dry fbAdj b2 wisp none |m|med afMax_med wet_med | ifc/partial at T/4 T/2 3T/4 T | partial T/2, T | Eshape T (rel) L1x L1af")
    def s(x, fmt="%g"):
        try:
            return "-" if x != x else (fmt % x)
        except TypeError:
            return str(x)
    for r in rows:
        gg = lambda k: r.get(k, float("nan"))
        print("%s %-12s %-8s %4d | %s %s | %s %s | %s/%s | %s/%s | %s/%s | %s: %s %s %s %s, %s | %s %s %s %s %s %s %s %s %s %s | %s %s %s %s | %s %s | %s (%s) %s %s %s" % (
            r["cell"], r["run"], r["status"], r["steps"], s(r["minAlphaRun"], "%.2e"), s(r["excessRun"], "%.2e"), s(r["minAlphaStep20"], "%.2e"), s(r["excessStep20"], "%.2e"),
            s(r["rankFbStep1"], "%.0f"), s(r["rankFbMax"], "%.0f"), s(r["bound2Step1"], "%.0f"), s(r["bound2Max"], "%.0f"), s(r["ifcStep1"], "%.0f"), s(r["ifcMax"], "%.0f"),
            s(r["fragTotal"], "%.0f"), s(r["fragFb"], "%.0f"), s(r["fragB2"], "%.0f"), s(r["fragW"], "%.0f"), s(r["fragNone"], "%.0f"), s(r["fragFirstStep"], "%.0f"),
            s(gg("frag200_n"), "%.0f"), s(gg("frag200_firstStep"), "%.0f"), s(gg("frag200_dry"), "%.0f"), s(gg("frag200_fallbackAdj"), "%.0f"), s(gg("frag200_bound2"), "%.0f"), s(gg("frag200_wisp"), "%.0f"), s(gg("frag200_none"), "%.0f"),
            s(gg("frag200_mOverV23_med"), "%.2e"), s(gg("frag200_afMax_med"), "%.3f"), s(gg("frag200_nWetFaces_med"), "%.0f"),
            s(r["ratioIfcPartial_0.75"], "%.2f"), s(r["ratioIfcPartial_1.5"], "%.2f"), s(r["ratioIfcPartial_2.25"], "%.2f"), s(r["ratioIfcPartial_3"], "%.2f"),
            s(r["nPartial_1.5"], "%.0f"), s(r["nPartial_3"], "%.0f"), s(r["EshapeT"], "%.3e"), s(r["EshapeTRel"], "%.3f"), s(r["L1xT"], "%.3f"), s(r["L1afT"], "%.3e"), r["abort"]))


if __name__ == "__main__":
    main()
