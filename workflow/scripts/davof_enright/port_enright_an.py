#!/usr/bin/env python3
"""Phase 3 of the port (docs/davof/PORT-DESIGN.md, 2026-11-06): the scoring of the
two-route Enright comparison written by port-enright-time.sbatch under
<workspace>/enright-time/N<N>_co<CO>/<method>_<route>_r<rep>/ (ring and rdf on the
routes leia and fork, prism on the fork route). The measures are enright_an.py's own
(its read_log, read_alpha, conv_at and counter regexes are imported unchanged); this
file adds the route, the identity verdict of the DAVOF arm between the routes
(port_identity_check.py on the first repetition's case directories) and the
wall-clock ratio leia/fork. Writes <workspace>/csv/port_enright_{runs,table}.csv.

    python3 port_enright_an.py [workspace]   (default /work/scratch/tm83tomy/davof/agents/port)
"""
import csv, glob, os, statistics, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import enright_an as E   # noqa: E402  (its module-level W/OUT are not used here)

W = sys.argv[1] if len(sys.argv) > 1 else "/work/scratch/tm83tomy/davof/agents/port"
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
nan = float("nan")
f = E.f

def main():
    runs = []
    for d in sorted(glob.glob(W + "/enright-time/N*_co*")):
        tf = d + "/timing.csv"
        if not os.path.isfile(tf): continue
        sac = {}
        for sf in glob.glob(d + "/sacct_steps*.txt"):
            for l in open(sf):
                p = l.strip().split("|")
                if len(p) >= 4 and "_r" in p[1]: sac[p[1]] = int(p[3])
        for t in csv.DictReader(open(tf)):
            c = "%s/%s_%s_r%s" % (d, t["method"], t["route"], t["rep"])
            logn = "log.davofAdvectorFoam" if t["route"] == "leia" else "log.advectorVoF"
            if not os.path.isfile(c + "/" + logn): continue
            steps, abort, ended = E.read_log(c + "/" + logn)
            R = E.conv_rows(c)
            tds = E.time_dirs(c)
            aT = E.read_alpha(c, E.nearest(tds, E.T_END)) if E.nearest(tds, E.T_END) else None
            aH = E.read_alpha(c, E.nearest(tds, E.T_HALF)) if E.nearest(tds, E.T_HALF) else None
            g = lambda k: [s[k] for s in steps if k in s]
            mins = g("minAlpha") + g("minAlphaBound"); exc = g("excess") + g("excessBound")
            sH = min(steps, key=lambda s: abs(s["t"] - E.T_HALF)) if steps else {}
            sT = steps[-1] if steps else {}
            davof = t["method"] in E.DAVOF
            run = dict(N=int(t["N"]), np=int(t["np"]), co_index=int(t["co_index"]), method=t["method"], route=t["route"], rep=int(t["rep"]),
                       status="aborted" if abort else ("ended" if ended else "incomplete"), abortStep=abort[0] if abort else "", abortReason=abort[2] if abort else "",
                       steps=len(steps), dt=f(t["dt"]), end_time=f(t["end_time"]), tLast=sT.get("t", nan),
                       ofCoMax=max(g("ofCo")) if g("ofCo") else nan, contErrMax=max(g("contErr")) if g("contErr") else nan,
                       minAlphaRun=min(mins) if mins else nan, excessRun=max(exc) if exc else nan,
                       massDevPctMax=max(abs(x) for x in g("devPct")) if g("devPct") else nan,
                       Eshape=E.conv_at(R, E.T_END, "Eshape"), EshapeRel=E.conv_at(R, E.T_END, "Eshape")/E.V_SPHERE,
                       L1x=E.conv_at(R, E.T_END, "LCentreDiff1"), L1n=E.conv_at(R, E.T_END, "LNormalDiff1"),
                       L1af=E.conv_at(R, E.T_END, "L1_alphaf") if davof else nan, EshapeHalf=E.conv_at(R, E.T_HALF, "Eshape"),
                       nInterfaceT=E.conv_at(R, E.T_END, "nInterface"), nInterfaceHalf=E.conv_at(R, E.T_HALF, "nInterface"),
                       nInterfaceAdvT=sT.get("nInterfaceAdv", nan), nInterfaceAdvHalf=sH.get("nInterfaceAdv", nan),
                       nPartialT=aT["nPartial"] if aT else nan, nPartialHalf=aH["nPartial"] if aH else nan,
                       rankFallbacksSum=sum(g("rankFallbacks")), beyondLayerSum=sum(g("beyondLayer")),
                       pyramidDeltaMaxOverH=max(g("pyramidDeltaMaxOverH")) if g("pyramidDeltaMaxOverH") else nan,
                       twistedFacesSum=sum(g("twistedFaces")), boundingSecondLayerSum=sum(g("boundingSecondLayer")),
                       anchorFailedSum=sum(g("anchorFailed")), ruleFallbacksSum=sum(g("ruleFallbacks")),
                       clock_s=f(t["clock_s"]), exec_s=f(t["exec_s"]), wall_s=f(t["wall_s"]),
                       sacct_s=sac.get("%s_%s_r%s" % (t["method"], t["route"], t["rep"]), nan),
                       median_step_s=f(t["median_step_s"]), node=t["node"], case=c)
            run["nWispsT"] = run["nPartialT"] - run["nInterfaceT"]; run["nWispsHalf"] = run["nPartialHalf"] - run["nInterfaceHalf"]
            runs.append(run)
    if not runs:
        print("no runs under", W + "/enright-time"); return
    with open(OUT + "/port_enright_runs.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)
    groups = {}
    for r in runs: groups.setdefault((r["N"], r["co_index"], r["method"], r["route"]), []).append(r)
    table = []
    for k in sorted(groups):
        rs = sorted(groups[k], key=lambda r: r["rep"]); r0 = rs[0]
        same = all(abs(r["Eshape"] - r0["Eshape"]) <= 1e-12*abs(r0["Eshape"]) + 1e-300 for r in rs if r["Eshape"] == r["Eshape"])
        med = lambda x: statistics.median([r[x] for r in rs])
        row = dict(N=k[0], np=r0["np"], co_index=k[1], ofCoMax=r0["ofCoMax"], method=k[2], route=k[3], reps=len(rs), repsIdentical=same,
                   status=",".join(sorted(set(r["status"] for r in rs))), steps=r0["steps"], dt=r0["dt"], tLast=r0["tLast"],
                   Eshape=r0["Eshape"], EshapeRel=r0["EshapeRel"], L1x=r0["L1x"], L1n=r0["L1n"], L1af=r0["L1af"], EshapeHalf=r0["EshapeHalf"],
                   minAlphaRun=r0["minAlphaRun"], excessRun=r0["excessRun"], massDevPctMax=r0["massDevPctMax"],
                   nInterfaceHalf=r0["nInterfaceHalf"], nInterfaceT=r0["nInterfaceT"], nPartialHalf=r0["nPartialHalf"], nPartialT=r0["nPartialT"],
                   nWispsHalf=r0["nWispsHalf"], nWispsT=r0["nWispsT"], rankFallbacksSum=r0["rankFallbacksSum"], beyondLayerSum=r0["beyondLayerSum"],
                   pyramidDeltaMaxOverH=r0["pyramidDeltaMaxOverH"], twistedFacesSum=r0["twistedFacesSum"], boundingSecondLayerSum=r0["boundingSecondLayerSum"],
                   anchorFailedSum=r0["anchorFailedSum"], ruleFallbacksSum=r0["ruleFallbacksSum"],
                   wall_med=med("wall_s"), wall_min=min(r["wall_s"] for r in rs), wall_max=max(r["wall_s"] for r in rs),
                   clock_med=med("clock_s"), step_med=med("median_step_s"), node=",".join(sorted(set(r["node"] for r in rs))), case=r0["case"])
        table.append(row)
    byk = {(r["N"], r["co_index"], r["method"], r["route"]): r for r in table}
    for r in table:
        o = byk.get((r["N"], r["co_index"], r["method"], "fork" if r["route"] == "leia" else "leia"))
        r["wall_ratio_leia_over_fork"] = (r["wall_med"]/o["wall_med"]) if (o and r["route"] == "leia" and o["wall_med"]) else nan
        r["identity"] = ""
        if o and r["route"] == "leia":
            if r["method"] in E.DAVOF:
                p = subprocess.run([sys.executable, os.path.join(HERE, "..", "port_identity_check.py"), o["case"], r["case"]], capture_output=True, text=True)
                r["identity"] = "bit for bit" if p.stdout.strip().endswith("ALL IDENTICAL") else "DIFFERENT (see port_identity_check)"
            else:
                r["identity"] = "same digits" if all(r[x] == o[x] for x in ("Eshape", "L1x", "L1n") if r[x] == r[x]) else \
                    "v2606 vs TwoPhaseFlow: " + ", ".join(f"{x} {r[x]:.6e} / {o[x]:.6e}" for x in ("Eshape", "L1x", "L1n", "minAlphaRun", "excessRun") if r[x] == r[x])
        rr = byk.get((r["N"], r["co_index"], "rdf", r["route"]))
        for key in ("wall_med", "Eshape", "L1x", "L1n"):
            r[key + "_ratio_to_rdf"] = (r[key]/rr[key]) if (rr and rr[key]) else nan
    with open(OUT + "/port_enright_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    print("| N | np | OF Co | method | route | status | steps | Eshape | Eshape/Vs | L1x | L1n | L1af | min alpha | excess | mass % max | nInt T/2, T | wisps T/2, T | wall med [min, max] s | wall leia/fork | wall/rdf | E/rdf | identity |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in table:
        fm = lambda x: ("%.3e" % r[x]) if r[x] == r[x] else "-"
        print(f"| {r['N']} | {r['np']} | {r['ofCoMax']:.3g} | {r['method']} | {r['route']} | {r['status']} | {r['steps']} | {fm('Eshape')} | {r['EshapeRel']:.4f} | {fm('L1x')} | {fm('L1n')} | {fm('L1af')} | {fm('minAlphaRun')} | {fm('excessRun')} | {fm('massDevPctMax')} | {r['nInterfaceHalf']:g}, {r['nInterfaceT']:g} | {r['nWispsHalf']:g}, {r['nWispsT']:g} | {r['wall_med']:.1f} [{r['wall_min']:.1f}, {r['wall_max']:.1f}] | {fm('wall_ratio_leia_over_fork')} | {fm('wall_med_ratio_to_rdf')} | {fm('Eshape_ratio_to_rdf')} | {r['identity']} |")
    print("written:", OUT + "/port_enright_runs.csv", OUT + "/port_enright_table.csv")

if __name__ == "__main__":
    main()
