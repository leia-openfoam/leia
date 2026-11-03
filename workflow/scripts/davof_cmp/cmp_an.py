"""The timed comparison (2026-11-03): error norms with the total wall clock at
the same ranks. Reads /work/.../cmp-time/<FLOW>_N<N>_co<CO>/, writes
csv/cmp_runs.csv (every run) and csv/cmp_table.csv (per case and method: the
metrics, median / min / max of the timings, the ratio to plicRDF)."""
import csv, glob, os, re, statistics, json
W = "/work/scratch/tm83tomy/davof/agents/swept"
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
runs = []
for d in sorted(glob.glob(W + "/cmp-time/*_N*_co*")):
    tfs = sorted(glob.glob(d + "/timing*.csv"))
    if not tfs: continue
    sac = {}
    for sf in glob.glob(d + "/sacct_steps*.txt"):
        for l in open(sf):
            p = l.strip().split("|")
            if len(p) >= 4 and "_r" in p[1]: sac[p[1]] = int(p[3])
    for t in (row for tf in tfs for row in csv.DictReader(open(tf))):
        c = "%s/%s_r%s" % (d, t["method"], t["rep"])
        log = open(c + "/log.advectorVoF").read()
        R = list(csv.DictReader(open(c + "/advectorVoF.csv"))) if os.path.isfile(c + "/advectorVoF.csv") else []
        mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
        exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)] + [0.0]
        cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
        beyond = [float(x) for x in re.findall(r"beyond layer \d+ \(max (\S+)\)", log)] + [0.0]
        rf = sum(int(x) for x in re.findall(r"rank fallbacks (\d+)", log))
        last = R[-1] if R else {}
        g = lambda k: float(last[k]) if last.get(k) not in (None, "") else float("nan")
        runs.append(dict(flow=t["flow"], N=int(t["N"]), np=int(t["np"]), co_index=int(t["co_index"]), method=t["method"], rep=int(t["rep"]), corr=t["corr"],
                         ofCo=max(cos) if cos else float("nan"), steps=int(t["steps"]), dt=float(t["dt"]), end_time=float(t["end_time"]),
                         clock_s=float(t["clock_s"] or "nan"), exec_s=float(t["exec_s"] or "nan"), wall_s=float(t["wall_s"]), sacct_s=sac.get("%s_r%s" % (t["method"], t["rep"]), float("nan")),
                         median_step_s=float(t["median_step_s"] or "nan"), node=t["node"], status=last.get("status", "noCSV"),
                         L1x=g("LCentreDiff1"), Eshape=g("Eshape"), L1n=g("LNormalDiff1"), L1af=g("L1_alphaf") if t["method"] != "rdf" else float("nan"),
                         minAlpha=min(mins) if mins else float("nan"), excess=max(exc), massDevPct=g("devPercent"), beyondMax=max(beyond), rankFb=rf))
with open(OUT + "/cmp_runs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)
groups = {}
for r in runs: groups.setdefault((r["flow"], r["N"], r["co_index"], r["method"]), []).append(r)
table = []
for k in sorted(groups):
    rs = groups[k]; r0 = rs[0]
    same = all(abs(r["L1x"] - r0["L1x"]) <= 1e-12 * abs(r0["L1x"]) + 1e-300 for r in rs if r["L1x"] == r["L1x"])
    med = lambda x: statistics.median([r[x] for r in rs])
    row = dict(flow=k[0], N=k[1], np=r0["np"], ofCo=r0["ofCo"], method=k[3], reps=len(rs), repsIdentical=same, end_time=r0["end_time"], steps=r0["steps"],
               L1x=r0["L1x"], Eshape=r0["Eshape"], L1n=r0["L1n"], L1af=r0["L1af"], minAlpha=r0["minAlpha"], excess=r0["excess"], massDevPct=r0["massDevPct"],
               beyondMax=r0["beyondMax"], rankFb=r0["rankFb"],
               clock_med=med("clock_s"), clock_min=min(r["clock_s"] for r in rs), clock_max=max(r["clock_s"] for r in rs),
               exec_med=med("exec_s"), wall_med=med("wall_s"), wall_min=min(r["wall_s"] for r in rs), wall_max=max(r["wall_s"] for r in rs),
               sacct_med=med("sacct_s"), sacct_min=min(r["sacct_s"] for r in rs), sacct_max=max(r["sacct_s"] for r in rs),
               step_med=med("median_step_s"), node=",".join(sorted(set(r["node"] for r in rs))), status=",".join(sorted(set(r["status"] for r in rs))))
    table.append(row)
for row in table:
    ref = [x for x in table if (x["flow"], x["N"], x["ofCo"] // 0.5, x["method"]) == (row["flow"], row["N"], row["ofCo"] // 0.5, "rdf")]
    rr = [x for x in table if x["flow"] == row["flow"] and x["N"] == row["N"] and x["method"] == "rdf" and abs(x["ofCo"] - row["ofCo"]) < 0.05]
    row["wall_ratio_to_rdf"] = row["wall_med"] / rr[0]["wall_med"] if rr else float("nan")
    row["L1x_ratio_to_rdf"] = row["L1x"] / rr[0]["L1x"] if rr else float("nan")
    row["Eshape_ratio_to_rdf"] = row["Eshape"] / rr[0]["Eshape"] if rr else float("nan")
with open(OUT + "/cmp_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
print("flow N np OFCo method | L1x Eshape L1n L1af | min excess mass% | wall med [min,max] s, sacct med, clock med, step med | wall/rdf L1x/rdf E/rdf | reps id status")
for r in table:
    print("%-8s %3d %2d %.3f %-5s | %.3e %.3e %.3e %.3e | %.1e %.1e %.1e | %.1f [%.1f,%.1f] %.0f %.0f %.4f | %.2f %.3f %.3f | %d %s %s" % (
        r["flow"], r["N"], r["np"], r["ofCo"], r["method"], r["L1x"], r["Eshape"], r["L1n"], r["L1af"], r["minAlpha"], r["excess"], r["massDevPct"],
        r["wall_med"], r["wall_min"], r["wall_max"], r["sacct_med"], r["clock_med"], r["step_med"], r["wall_ratio_to_rdf"], r["L1x_ratio_to_rdf"], r["Eshape_ratio_to_rdf"],
        r["reps"], r["repsIdentical"], r["status"]))
print("wrote", OUT + "/cmp_runs.csv", OUT + "/cmp_table.csv")
