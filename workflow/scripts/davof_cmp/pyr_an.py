"""Item 1 (2026-11-03): the pyramid against the scaling. Writes CSVs for the
writing agent to /work/scratch/tm83tomy/davof/agents/swept/csv/."""
import json, glob, os, re, csv, sys
S = "/work/scratch/tm83tomy/davof/agents/swept/leia/studies/davof"
OUT = "/work/scratch/tm83tomy/davof/agents/swept/csv"
os.makedirs(OUT, exist_ok=True)

def cases(st):
    out = {}
    for c in sorted(glob.glob(os.path.join(S, st, "*_0*"))):
        if "." in os.path.basename(c) or not os.path.isfile(c + "/advectorVoF.csv"): continue
        t = json.load(open(c + "/case_params.json"))["tokens"]
        key = (t.get("TPF_VELOCITY_ARM"), t["TPF_ROW"].replace("davof-planeSurf-", ""), float(t["TARGET_ADVECTIVE_CO"]), int(float(t["N_CELLS"])))
        out[key] = c
    return out

def info(c):
    R = list(csv.DictReader(open(c + "/advectorVoF.csv")))
    log = open(c + "/log.advectorVoF").read()
    f = lambda k: [float(r[k]) for r in R if r.get(k) not in ("", None)]
    vc = re.findall(r"volume correction (\w+): max \|V_f - V_g\|/V_f (\S+), max \|delta\|/h (\S+), twisted faces (\d+), pyramid skipped (\d+), max axis displacement/h ([0-9.eE+-]+)", log)
    mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
    exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)] + [0.0]
    cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
    beyond = [float(x) for x in re.findall(r"beyond layer \d+ \(max (\S+)\)", log)] + [0.0]
    rf = sum(int(x) for x in re.findall(r"rank fallbacks (\d+)", log))
    d = dict(L1af_max=max(f("L1_alphaf")), L1n_max=max(f("LNormalDiff1")), L1af_end=f("L1_alphaf")[-1], L1x_end=f("LCentreDiff1")[-1],
             Eshape_end=f("Eshape")[-1], L1n_end=f("LNormalDiff1")[-1], minAlpha=min(mins) if mins else float("nan"), excess=max(exc),
             massDevPct=float(R[-1]["devPercent"]), EshapeRel=f("Eshape")[-1] / float(R[-1]["sumAlphaV"]), ofCo=max(cos) if cos else float("nan"),
             defectMax=max(float(v[1]) for v in vc) if vc else float("nan"), deltaMax=max(float(v[2]) for v in vc) if vc else float("nan"),
             twisted=sum(int(v[3]) for v in vc), skipped=sum(int(v[4]) for v in vc), dispMax=max(float(v[5]) for v in vc) if vc else float("nan"),
             beyondMax=max(beyond), rankFb=rf, status=R[-1]["status"], clock=float(R[-1]["clockTime"]) if R[-1].get("clockTime") else float("nan"))
    return d

rows = []
# the plane: pyramid against M9 (scale)
P, M = cases("planeTranslation3D_hexSweptPyr"), cases("planeTranslation3D_hexSweptWispM9")
print("== plane: pyramid against M9 (scale)")
gate = lambda d: d["L1af_max"] < 1e-12 and d["L1n_max"] < 1e-12 and d["minAlpha"] >= -1e-12 and d["excess"] <= 1e-11 and d["EshapeRel"] < 1e-12
np_ = nm = 0; dmax = 0; defmax = 0; delmax = 0
for k in sorted(P):
    a, b = info(P[k]), info(M[k])
    np_ += gate(a); nm += gate(b)
    if a["ofCo"] <= 1: dmax = max(dmax, abs(a["L1af_max"] - b["L1af_max"]), abs(a["L1n_max"] - b["L1n_max"]))
    defmax = max(defmax, a["defectMax"]); delmax = max(delmax, a["deltaMax"])
    rows.append(dict(study="plane", arm=k[0], row=k[1], Co=k[2], N=k[3], corr="pyramid", **a)); rows.append(dict(study="plane", arm=k[0], row=k[1], Co=k[2], N=k[3], corr="scale(M9)", **b))
print("  cases %d; gate passes pyramid %d, M9 %d; max |dL1| (OpenFOAM Co <= 1) %.2e; max |V_f - V_g|/V_f %.2e, max |delta|/h %.2e, twisted %d" % (
      len(P), np_, nm, dmax, defmax, delmax, sum(info(P[k])["twisted"] for k in P)))
# curved pairs
for name in ("hexLineXCo", "hexRotCo", "hexSpinCo", "hexRotN128Co"):
    A, B = cases("sphereRotation3D_" + name + "Scale"), cases("sphereRotation3D_" + name + "Pyr")
    if not A or not B: print("==", name, "missing"); continue
    print("== %s: pyramid against scale" % name)
    for k in sorted(B):
        if k not in A: continue
        a, b = info(A[k]), info(B[k])
        rel = lambda x: (b[x] - a[x]) / a[x] if a[x] else float("nan")
        print("  %s %s Co %.2f N %d (OF Co %.3f): L1af %.4e/%.4e (%+.2e) L1x %.4e/%.4e (%+.2e) E %.4e/%.4e (%+.2e) | min %.1e/%.1e exc %.1e/%.1e mass %.1e/%.1e%% | defect %.2e delta/h %.2e twisted %d/%d disp %.3f beyond %.2e rankFb %d/%d | %s/%s" % (
            k[0], k[1], k[2], k[3], b["ofCo"], a["L1af_end"], b["L1af_end"], rel("L1af_end"), a["L1x_end"], b["L1x_end"], rel("L1x_end"), a["Eshape_end"], b["Eshape_end"], rel("Eshape_end"),
            a["minAlpha"], b["minAlpha"], a["excess"], b["excess"], a["massDevPct"], b["massDevPct"], b["defectMax"], b["deltaMax"], a["twisted"], b["twisted"], b["dispMax"], b["beyondMax"], a["rankFb"], b["rankFb"], a["status"], b["status"]))
        for corr, d in (("scale", a), ("pyramid", b)):
            rows.append(dict(study=name, arm=k[0], row=k[1], Co=k[2], N=k[3], corr=corr, **d))
with open(OUT + "/pyramid_item1.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote", OUT + "/pyramid_item1.csv", len(rows))
