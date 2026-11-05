"""The rerun of the reference setting (fork 407f969, no keys in the dictionaries), scored against the
opt-in runs (Claude, 2026-11-04): the plane gate against planeTranslation3D_hexSweptWispM9Ccs2, the
timed cells against the comparison (csv/cmp_table.csv), the spin against sphereRotation3D_hexSpinCoPyr,
the parallel identity against np 1. Writes csv/ref_setting.csv (one row per case and method) and
prints the markdown tables and the verdicts D-1 .. D-4.
    python3 ref_an.py"""
import csv, glob, json, os, re, statistics, math
W = "/work/scratch/tm83tomy/davof/agents/swept"
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
nan = float("nan")

def logs(c):
    log = open(c + "/log.advectorVoF").read()
    mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
    exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)] + [0.0]
    cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
    rf = sum(int(x) for x in re.findall(r"rank fallbacks (\d+)", log))
    st = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)]
    keys = re.search(r"sweptPreTest (\w+), sweptStateShortcut (\d), sweptTraceCache (\d)", log)
    return dict(minAlpha=min(mins) if mins else nan, excess=max(exc), ofCo=max(cos) if cos else nan, rankFb=rf, steps=len(re.findall(r"\nTime = ", log)),
                step_med=statistics.median(st) if st else nan, keys=keys.group(0) if keys else "", clock=float(re.findall(r"ClockTime = (\S+) s", log)[-1]) if "ClockTime" in log else nan,
                wisp=re.search(r"sweptWispSurface \w+, sweptWispTol \S+, sweptWispMatchTol \S+", log).group(0) if re.search(r"sweptWispSurface \w+, sweptWispTol", log) else "")

def metrics(c):
    p = c + "/advectorVoF.csv"
    if not os.path.isfile(p): return {}
    R = list(csv.DictReader(open(p)))
    if not R: return {}
    f = lambda k: [float(r[k]) for r in R if r.get(k) not in ("", None)]
    last = R[-1]
    g = lambda k: float(last[k]) if last.get(k) not in (None, "") else nan
    d = dict(L1af_end=g("L1_alphaf"), L1x_end=g("LCentreDiff1"), Eshape_end=g("Eshape"), L1n_end=g("LNormalDiff1"), massDevPct=g("devPercent"), status=last.get("status", ""))
    d["L1af_max"] = max(f("L1_alphaf")) if f("L1_alphaf") else nan
    d["L1n_max"] = max(f("LNormalDiff1")) if f("LNormalDiff1") else nan
    d["EshapeRel"] = d["Eshape_end"] / float(last["sumAlphaV"]) if last.get("sumAlphaV") else nan
    return d

def tokens(c):
    src = open(c + "/.source").read().strip() if os.path.isfile(c + "/.source") else c
    return json.load(open(src + "/case_params.json"))["tokens"], src

rel = lambda a, b: (a - b) / b if b not in (0, None) and b == b else nan
rows = []
# ---- the plane gate
gate = lambda d: d["L1af_max"] < 1e-12 and d["L1n_max"] < 1e-12 and d["minAlpha"] >= -1e-12 and d["excess"] <= 1e-11 and d["EshapeRel"] < 1e-12
print("### The plane gate (32 cases, no keys) against M9Ccs2\n")
print("| N | Co | arm | row | gate ref / new | L1(alphaf) max ref / new | L1(n) max ref / new | min alpha new | excess new | mass % new | steps |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
nsame = ntot = npass = 0
planes = {}   # the split submissions (ref/plane, ref/plane_small): one case directory per case, the newest scored copy
for c in glob.glob(W + "/ref/plane*/planeTranslation3D_0*"):
    if not os.path.isfile(c + "/advectorVoF.csv"): continue
    k = os.path.basename(c)
    if k not in planes or os.path.getmtime(c + "/advectorVoF.csv") > os.path.getmtime(planes[k] + "/advectorVoF.csv"): planes[k] = c
for c in sorted(planes.values()):
    t, src = tokens(c)
    a = {**metrics(c), **logs(c)}; b = {**metrics(src), **logs(src)}
    if not a.get("status"): print("| %s | | | | NO CSV | | | | | | |" % os.path.basename(c)); continue
    ga, gb = gate(a), gate(b); ntot += 1; nsame += (ga == gb); npass += ga
    print("| %s | %s | %s | %s | %s / %s | %.2e / %.2e | %.2e / %.2e | %.1e | %.1e | %.1e | %d |" % (int(float(t["N_CELLS"])), t["TARGET_ADVECTIVE_CO"], t["TPF_VELOCITY_ARM"], t["TPF_ROW"].replace("davof-planeSurf-", ""),
          gb, ga, b["L1af_max"], a["L1af_max"], b["L1n_max"], a["L1n_max"], a["minAlpha"], a["excess"], a["massDevPct"], a["steps"]))
    rows.append(dict(cell="plane", flow="plane", N=int(float(t["N_CELLS"])), Co=t["TARGET_ADVECTIVE_CO"], arm=t["TPF_VELOCITY_ARM"], row=t["TPF_ROW"].replace("davof-planeSurf-", ""), method="ref", np=1,
                     gate=ga, gate_ref=gb, L1af_max=a["L1af_max"], L1af_max_ref=b["L1af_max"], L1n_max=a["L1n_max"], L1n_max_ref=b["L1n_max"], minAlpha=a["minAlpha"], excess=a["excess"], massDevPct=a["massDevPct"],
                     steps=a["steps"], keys=a["keys"], wisp=a["wisp"], status=a["status"]))
print("\nplane: %d cases, %d verdicts equal to M9Ccs2's, %d pass" % (ntot, nsame, npass))
# ---- the timed cells
cmp = {}
for r in csv.DictReader(open(W + "/csv/cmp_table.csv")):
    cmp[(r["flow"], int(r["N"]), int(round(float(r["ofCo"]))) if float(r["ofCo"]) > 0.5 else 0, r["method"])] = r
print("\n### The timed cells (ring 3 repetitions, exact 1) against the comparison of 2026-11-03\n")
print("| cell | OF Co | row | E_shape new (rel. to cmp) | L1(alphaf) (rel.) | L1(x) (rel.) | L1(n) (rel.) | min / excess | mass % | wall new med [min, max] s | wall cmp med | cmp / new | step new / cmp s |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for d in sorted(glob.glob(W + "/ref/cmp/*_N*_co*")):
    tf = d + "/timing.csv"
    if not os.path.isfile(tf): continue
    m = re.match(r".*/(\w+)_N(\d+)_co(\d)$", d); flow, N, co = m.group(1), int(m.group(2)), int(m.group(3))
    T = list(csv.DictReader(open(tf)))
    for meth in ("ring", "exact"):
        reps = [t for t in T if t["method"] == meth]
        if not reps: continue
        c = "%s/%s_r1" % (d, meth)
        a = {**metrics(c), **logs(c)}
        walls = [float(t["wall_s"]) for t in reps]
        ref = cmp.get((flow, N, co, meth))
        g = lambda k: float(ref[k]) if ref and ref.get(k) not in ("", None) else nan
        wr = g("wall_med") / statistics.median(walls) if ref else nan
        print("| %s N %d, %d ranks | %.3f | %s | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.1e / %.1e | %.1e | %.1f [%.1f, %.1f] | %.1f | %.2f | %.3f / %.3f |" % (
            flow, N, int(reps[0]["np"]), a["ofCo"], meth, a["Eshape_end"], rel(a["Eshape_end"], g("Eshape")), a["L1af_end"], rel(a["L1af_end"], g("L1af")), a["L1x_end"], rel(a["L1x_end"], g("L1x")),
            a["L1n_end"], rel(a["L1n_end"], g("L1n")), a["minAlpha"], a["excess"], a["massDevPct"], statistics.median(walls), min(walls), max(walls), g("wall_med"), wr, a["step_med"], g("step_med")))
        rows.append(dict(cell="cmp", flow=flow, N=N, Co=co, arm="", row=meth, method="ref", np=int(reps[0]["np"]), ofCo=a["ofCo"], steps=a["steps"],
                         Eshape=a["Eshape_end"], Eshape_rel=rel(a["Eshape_end"], g("Eshape")), L1af=a["L1af_end"], L1af_rel=rel(a["L1af_end"], g("L1af")), L1x=a["L1x_end"], L1x_rel=rel(a["L1x_end"], g("L1x")),
                         L1n=a["L1n_end"], L1n_rel=rel(a["L1n_end"], g("L1n")), minAlpha=a["minAlpha"], excess=a["excess"], massDevPct=a["massDevPct"], rankFb=a["rankFb"],
                         wall_med=statistics.median(walls), wall_min=min(walls), wall_max=max(walls), wall_cmp_med=g("wall_med"), wall_ratio_cmp_over_new=wr, step_med=a["step_med"], step_cmp=g("step_med"),
                         keys=a["keys"], wisp=a["wisp"], status=a["status"]))
# ---- the spin
print("\n### The pure spin against sphereRotation3D_hexSpinCoPyr\n")
print("| N | target Co | OF Co | E_shape new (rel.) | L1(alphaf) (rel.) | L1(x) (rel.) | L1(n) (rel.) | min / excess | mass % | step s |")
print("|---|---|---|---|---|---|---|---|---|---|")
for c in sorted(glob.glob(W + "/ref/spin/sphereRotation3D_0*")):
    t, src = tokens(c)
    a = {**metrics(c), **logs(c)}; b = {**metrics(src), **logs(src)}
    if not a.get("status"): continue
    print("| %s | %s | %.3f | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.1e / %.1e | %.1e | %.3f |" % (int(float(t["N_CELLS"])), t["TARGET_ADVECTIVE_CO"], a["ofCo"],
          a["Eshape_end"], rel(a["Eshape_end"], b["Eshape_end"]), a["L1af_end"], rel(a["L1af_end"], b["L1af_end"]), a["L1x_end"], rel(a["L1x_end"], b["L1x_end"]), a["L1n_end"], rel(a["L1n_end"], b["L1n_end"]),
          a["minAlpha"], a["excess"], a["massDevPct"], a["step_med"]))
    rows.append(dict(cell="spin", flow="spin", N=int(float(t["N_CELLS"])), Co=t["TARGET_ADVECTIVE_CO"], arm="rotation", row="exact", method="ref", np=1, ofCo=a["ofCo"], steps=a["steps"],
                     Eshape=a["Eshape_end"], Eshape_rel=rel(a["Eshape_end"], b["Eshape_end"]), L1af=a["L1af_end"], L1af_rel=rel(a["L1af_end"], b["L1af_end"]), L1x=a["L1x_end"], L1x_rel=rel(a["L1x_end"], b["L1x_end"]),
                     L1n=a["L1n_end"], L1n_rel=rel(a["L1n_end"], b["L1n_end"]), minAlpha=a["minAlpha"], excess=a["excess"], massDevPct=a["massDevPct"], step_med=a["step_med"], keys=a["keys"], wisp=a["wisp"], status=a["status"]))
# ---- the parallel identity
print("\n### The parallel identity (rotation N 32, 20 steps) against np 1\n")
print("| row | np | step 1 max dalpha / dalphaf | max over 20 steps dalpha / dalphaf | step of the max | nMissing / nDup | rank fallbacks |")
print("|---|---|---|---|---|---|---|")
for i, rowname in ((0, "exact"), (1, "ring")):
    for k in (2, 4, 8, 16):
        f = "%s/ref/par/cmp_np%d_%d.csv" % (W, k, i)
        if not os.path.isfile(f): continue
        R = list(csv.DictReader(open(f)))
        if not R: continue
        g = lambda r, key: float(r[key]) if r.get(key) not in (None, "") else nan
        imax = max(range(len(R)), key=lambda j: g(R[j], "maxAlphaDiff"))
        lg = logs("%s/ref/par/np%d_%d" % (W, k, i))
        print("| %s | %d | %.1e / %.1e | %.1e / %.1e | %s | %d / %d | %d |" % (rowname, k, g(R[0], "maxAlphaDiff"), g(R[0], "maxAlphafDiff"), g(R[imax], "maxAlphaDiff"), max(g(r, "maxAlphafDiff") for r in R),
              R[imax].get("time", R[imax].get("step", imax + 1)), sum(int(float(r.get("nMissing", 0) or 0)) for r in R), sum(int(float(r.get("nDup", 0) or 0)) for r in R), lg["rankFb"]))
        rows.append(dict(cell="par", flow="Rot", N=32, Co=0.2, arm="rotation", row=rowname, method="ref", np=k, steps=len(R), dalpha_step1=g(R[0], "maxAlphaDiff"), dalphaf_step1=g(R[0], "maxAlphafDiff"),
                         dalpha_max=g(R[imax], "maxAlphaDiff"), dalphaf_max=max(g(r, "maxAlphafDiff") for r in R), nMissing=sum(int(float(r.get("nMissing", 0) or 0)) for r in R), nDup=sum(int(float(r.get("nDup", 0) or 0)) for r in R), rankFb=lg["rankFb"]))
keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open(OUT + "/ref_setting.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
print("\nwrote", OUT + "/ref_setting.csv", len(rows), "rows")
