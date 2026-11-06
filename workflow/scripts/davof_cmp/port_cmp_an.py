#!/usr/bin/env python3
"""The port's timed comparison, scored (docs/davof/PORT-DESIGN.md 4.4, 2026-11-06). Reads
<workspace>/cmp-time/<FLOW>_N<N>_co<CO>/ written by port-cmp-time.sbatch (the methods ring,
exact, rdf on both routes, leia and fork) and writes <workspace>/csv/port_cmp_runs.csv (every
run) and port_cmp_table.csv (per cell, method and route: the end-of-run norms, the
boundedness, the wall-clock medians, the ratio leia/fork of the wall clock, and for the
DAVOF arm the identity verdict of the two routes by port_identity_check.py on the first
repetition's case directories).

    python3 port_cmp_an.py [workspace]    (default /work/scratch/tm83tomy/davof/agents/port)

The metrics follow the fork's cmp_an.py (2026-11-03) line by line: L1x = LCentreDiff1,
Eshape, L1n = LNormalDiff1, L1af = L1_alphaf (the DAVOF arm), minAlpha and the largest
excess over 1 after the conservative bounding from the log, massDevPct = devPercent.
"""
import csv, glob, os, re, statistics, subprocess, sys
W = sys.argv[1] if len(sys.argv) > 1 else "/work/scratch/tm83tomy/davof/agents/port"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
nan = float("nan")
runs = []
for d in sorted(glob.glob(W + "/cmp-time/*_N*_co*")):
    tf = d + "/timing.csv"
    if not os.path.isfile(tf): continue
    for t in csv.DictReader(open(tf)):
        c = "%s/%s_%s_r%s" % (d, t["method"], t["route"], t["rep"])
        logn = "log.davofAdvectorFoam" if t["route"] == "leia" else "log.advectorVoF"
        if not os.path.isfile(c + "/" + logn): continue
        log = open(c + "/" + logn).read()
        R = list(csv.DictReader(open(c + "/advectorVoF.csv"))) if os.path.isfile(c + "/advectorVoF.csv") else []
        mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
        exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)] + [0.0]
        cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
        rf = sum(int(x) for x in re.findall(r"rank fallbacks (\d+)", log))
        last = R[-1] if R else {}
        g = lambda k: float(last[k]) if last.get(k) not in (None, "") else nan
        runs.append(dict(flow=t["flow"], N=int(t["N"]), np=int(t["np"]), co_index=int(t["co_index"]), method=t["method"], route=t["route"], rep=int(t["rep"]),
                         ofCo=max(cos) if cos else nan, steps=int(t["steps"]), dt=float(t["dt"]), end_time=float(t["end_time"]),
                         clock_s=float(t["clock_s"] or "nan"), exec_s=float(t["exec_s"] or "nan"), wall_s=float(t["wall_s"]),
                         median_step_s=float(t["median_step_s"] or "nan"), node=t["node"], status=last.get("status", "noCSV"),
                         L1x=g("LCentreDiff1"), Eshape=g("Eshape"), L1n=g("LNormalDiff1"), L1af=g("L1_alphaf") if t["method"] != "rdf" else nan,
                         minAlpha=min(mins) if mins else nan, excess=max(exc), massDevPct=g("devPercent"), rankFb=rf, case=c))
if not runs:
    sys.exit("no runs under " + W + "/cmp-time")
with open(OUT + "/port_cmp_runs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)
groups = {}
for r in runs: groups.setdefault((r["flow"], r["N"], r["co_index"], r["method"], r["route"]), []).append(r)
table = []
for k in sorted(groups):
    rs = sorted(groups[k], key=lambda r: r["rep"]); r0 = rs[0]
    med = lambda x: statistics.median([r[x] for r in rs])
    row = dict(flow=k[0], N=k[1], co_index=k[2], method=k[3], route=k[4], np=r0["np"], ofCo=r0["ofCo"], reps=len(rs), end_time=r0["end_time"], steps=r0["steps"],
               L1x=r0["L1x"], Eshape=r0["Eshape"], L1n=r0["L1n"], L1af=r0["L1af"], minAlpha=r0["minAlpha"], excess=r0["excess"], massDevPct=r0["massDevPct"], rankFb=r0["rankFb"],
               wall_med=med("wall_s"), wall_min=min(r["wall_s"] for r in rs), wall_max=max(r["wall_s"] for r in rs), clock_med=med("clock_s"), step_med=med("median_step_s"),
               node=",".join(sorted(set(r["node"] for r in rs))), status=",".join(sorted(set(r["status"] for r in rs))), case=r0["case"])
    table.append(row)
# the two routes side by side: the wall-clock ratio and the identity of the DAVOF arm
byk = {(r["flow"], r["N"], r["co_index"], r["method"], r["route"]): r for r in table}
for r in table:
    o = byk.get((r["flow"], r["N"], r["co_index"], r["method"], "fork" if r["route"] == "leia" else "leia"))
    r["wall_ratio_leia_over_fork"] = (r["wall_med"] / o["wall_med"]) if (o and r["route"] == "leia" and o["wall_med"]) else nan
    r["identity"] = ""
    if o and r["route"] == "leia":
        same = all(abs(r[x] - o[x]) <= 1e-300 or r[x] == o[x] for x in ("L1x", "Eshape", "L1n", "minAlpha", "excess", "massDevPct") if r[x] == r[x])
        if r["method"] != "rdf":
            p = subprocess.run([sys.executable, os.path.join(HERE, "..", "port_identity_check.py"), o["case"], r["case"]], capture_output=True, text=True)
            r["identity"] = ("bit for bit" if p.stdout.strip().endswith("ALL IDENTICAL") else "DIFFERENT: " + " | ".join(x.strip() for x in p.stdout.split("\n") if x.startswith("  ") and "identical" not in x)[:300])
        else:
            r["identity"] = "same digits" if same else "v2606 vs TwoPhaseFlow: " + ", ".join(f"{x} {r[x]:.6e} / {o[x]:.6e}" for x in ("L1x", "Eshape", "L1n", "minAlpha", "excess") if r[x] == r[x])
    # the comparison of the DAVOF arm with the plicRDF arm on the same route
    rdf = byk.get((r["flow"], r["N"], r["co_index"], "rdf", r["route"]))
    r["wall_ratio_to_rdf"] = (r["wall_med"] / rdf["wall_med"]) if (rdf and rdf["wall_med"]) else nan
    r["L1x_ratio_to_rdf"] = (r["L1x"] / rdf["L1x"]) if (rdf and rdf["L1x"]) else nan
    r["Eshape_ratio_to_rdf"] = (r["Eshape"] / rdf["Eshape"]) if (rdf and rdf["Eshape"]) else nan
with open(OUT + "/port_cmp_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
print("| flow | N | np | OF Co | method | route | L1x | Eshape | L1n | min alpha | excess | mass % | wall med [min, max] s | wall leia/fork | wall/rdf | identity |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in table:
    f = lambda x: ("%.4e" % r[x]) if r[x] == r[x] else "-"
    print(f"| {r['flow']} | {r['N']} | {r['np']} | {r['ofCo']:.3g} | {r['method']} | {r['route']} | {f('L1x')} | {f('Eshape')} | {f('L1n')} | {f('minAlpha')} | {f('excess')} | {f('massDevPct')} | {r['wall_med']:.1f} [{r['wall_min']:.1f}, {r['wall_max']:.1f}] | {f('wall_ratio_leia_over_fork')} | {f('wall_ratio_to_rdf')} | {r['identity']} |")
print("written:", OUT + "/port_cmp_runs.csv", OUT + "/port_cmp_table.csv")
