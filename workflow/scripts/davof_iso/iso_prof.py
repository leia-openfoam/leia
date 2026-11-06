#!/usr/bin/env python3
"""The reconstruction's share of the step from OpenFOAM's profiling of the profiled
runs (<W>/time/<cell>/<method>_rprof, uniform/profiling of the last write, rank 0 in
parallel): per run the wall time of time.run(), the steps, and the totals of
reconstruct() (plicRDF / davof / isoAlpha), the davof normal model's select(),
isoAdvection::advect(), timeIntegratedFlux(), limitFluxes(), sweepSecondRing(),
boundCellContent(), defaultAreaFractions() in seconds and per step; writes
<W>/csv/iso_prof.csv and prints a table. The profiled run is one run (not timed
for the table); its time.run() per step is also listed."""
import csv, glob, os, re, sys
W = os.environ.get("ISO_W", "/work/scratch/tm83tomy/davof/agents/swept/iso")
KEYS = [("reconstruct", r"::reconstruct\(bool\)"), ("select", r"davofNormalModel::select"), ("advect", r"isoAdvection::advect\("),
        ("flux", r"timeIntegratedFlux"), ("limitFluxes", r"::limitFluxes"), ("secondRing", r"sweepSecondRing"),
        ("contentBound", r"boundCellContent"), ("defaults", r"defaultAreaFractions"), ("timerun", r"^time\.run\(\)"), ("main", r"^application::main")]
def parse(p):
    s = open(p).read()
    out = {}
    for blk in re.finditer(r"\{(.*?)\n\s*\}", s, re.S):
        b = blk.group(1)
        md = re.search(r'description\s+"?([^";]+)"?;', b); mc = re.search(r"calls\s+(\d+);", b); mt = re.search(r"totalTime\s+(\S+);", b)
        if not (md and mc and mt): continue
        d = md.group(1).strip(); calls = int(mc.group(1)); tot = float(mt.group(1))
        for k, pat in KEYS:
            if re.search(pat, d) and calls >= 10 or (re.search(pat, d) and k in ("main",)):
                out[k] = out.get(k, 0.0) + tot
                out[k + "_calls"] = out.get(k + "_calls", 0) + calls
    return out
rows = []
for c in sorted(glob.glob(W + "/time/*/*_rprof")):
    cell = os.path.basename(os.path.dirname(c)); m = os.path.basename(c).replace("_rprof", "")
    fs = sorted(glob.glob(c + "/processor0/[0-9]*/uniform/profiling")) or sorted(glob.glob(c + "/[0-9]*/uniform/profiling"))
    if not fs: continue
    f = sorted(fs, key=lambda p: float(p.split("/")[-3]))[-1]
    pr = parse(f)
    log = open(c + "/log.advectorVoF").read(); steps = len(re.findall(r"\nTime = ", log))
    mm = re.match(r"(\w+?)_(hex|perturbed|poly|tet)_N(\d+)_co(\d)$", cell) or re.match(r"(\w+)_N(\d+)_co(\d)$", cell)
    g = mm.groups(); flow = g[0]; mesh = g[1] if len(g) == 4 else "hex"; N = int(g[-2]); co = {"0": 0.2, "1": 0.5, "2": 1.0}[g[-1]]
    r = dict(flow=flow, mesh=mesh, N=N, co=co, method=m, steps=steps, timerun_s=pr.get("timerun", float("nan")))
    for k, _ in KEYS:
        if k in ("timerun", "main"): continue
        r[k + "_s"] = pr.get(k, 0.0); r[k + "_per_step_ms"] = 1e3*pr.get(k, 0.0)/max(steps, 1)
        r[k + "_share_pct"] = 100*pr.get(k, 0.0)/r["timerun_s"] if r["timerun_s"] else float("nan")
    r["step_ms"] = 1e3*r["timerun_s"]/max(steps, 1)
    rows.append(r)
os.makedirs(W + "/csv", exist_ok=True)
with open(W + "/csv/iso_prof.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("%-12s %-9s %4s %4s %-8s %6s %8s | %9s %6s | %8s %8s | %8s %8s %8s %8s" % ("flow", "mesh", "N", "Co", "method", "steps", "step ms", "recon ms", "share%", "advect ms", "flux ms", "ring ms", "2ndRing", "content", "defaults"))
for r in rows:
    print("%-12s %-9s %4d %4.1f %-8s %6d %8.2f | %9.3f %6.1f | %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f" % (
        r["flow"], r["mesh"], r["N"], r["co"], r["method"], r["steps"], r["step_ms"], r["reconstruct_per_step_ms"], r["reconstruct_share_pct"],
        r["advect_per_step_ms"], r["flux_per_step_ms"], r["select_per_step_ms"], r["secondRing_per_step_ms"], r["contentBound_per_step_ms"], r["defaults_per_step_ms"]))
print("wrote", W + "/csv/iso_prof.csv")
