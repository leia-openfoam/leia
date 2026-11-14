#!/usr/bin/env python3
"""diag2 Part 2 analysis (markdown to stdout). Usage: an2.py <leia>/studies/davof

Pieces cases: the davofProf lines of the log (steps >= 2). plicRDF / isoAlpha cases: OpenFOAM's profiling
(uniform/profiling of the last write time, mean and max over the ranks) and the log."""
import csv, glob, os, re, sys, importlib.util

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("prof", here + "/prof.py")
prof = importlib.util.module_from_spec(spec); spec.loader.exec_module(prof)

root = sys.argv[1]
rx_co = re.compile(r"^Courant Number mean: \S+ max: (\S+)")
rx_exec = re.compile(r"^ExecutionTime = (\S+) s\s+ClockTime = (\S+) s")
rx_int = re.compile(r"^davof: (\d+) interface cells")

def db(study):
    f = glob.glob("%s/%s/*_database.csv" % (root, study))
    return {r["index"]: r for r in csv.DictReader(open(f[0]))} if f else {}

def _isnum(x):
    try:
        float(x); return True
    except ValueError:
        return False

def parse_prof_line(line):
    t = line.split()
    d = {"T": {}, "TM": {}, "Tn": {}, "N": {}, "NM": {}, "X": []}
    d["step"] = int(t[2]); d["ranks"] = int(t[4])
    i = 5
    while i < len(t):
        k = t[i]
        if k.startswith("T_") or k.startswith("N_") or k == "X_swept":
            j = i + 1; vals = []
            while j < len(t) and _isnum(t[j]):
                vals.append(float(t[j])); j += 1
            if k.startswith("T_"):
                d["T"][k[2:]] = vals[0]; d["TM"][k[2:]] = vals[1]
                if len(vals) > 2: d["Tn"][k[2:]] = vals[2]
            elif k.startswith("N_"):
                d["N"][k[2:]] = vals[0]
                if len(vals) > 1: d["NM"][k[2:]] = vals[1]
            else:
                d["X"] = vals[:6]
            i = j
        else:
            i += 1
    return d

def read_log(path):
    steps = 0; co = 0.0; ex = (0, 0); ni = []
    prof_lines = []
    for line in open(path, errors="replace"):
        if line.startswith("Time = "): steps += 1
        m = rx_co.match(line)
        if m: co = max(co, float(m.group(1)))
        m = rx_exec.match(line)
        if m: ex = (float(m.group(1)), float(m.group(2)))
        m = rx_int.match(line)
        if m: ni.append(int(m.group(1)))
        if line.startswith("davofProf:"): prof_lines.append(line)
    return steps, co, ex, ni, prof_lines

def mean(v): return sum(v)/len(v) if v else 0.0

LEAVES = ["between", "reinit", "rcAreaNormal", "rcCells", "rcExplicit", "rcSelect", "rcPositioning", "rcTail",
          "cfInit", "cfSurfaces", "cfVelSetup", "cfExch1", "cfWisp", "cfCand", "cfExch2", "cfSweptExch",
          "fTrace", "fPre", "fVol", "fArea", "fAssembly", "fOverhead", "cfProcSync", "cfTail",
          "alphaUpdate", "limit", "brute", "defaultAf", "report"]
DESC = {
    "between": "outside advect(): velocity update, Courant number, function objects, IO",
    "reinit": "reinitFaceFractions", "rcAreaNormal": "reconstruction: m_k from alpha_f, hasValue",
    "rcCells": "reconstruction: interface-cell rule", "rcExplicit": "reconstruction: face moment + explicit-plane cut",
    "rcSelect": "reconstruction: ring normal (normal model)", "rcPositioning": "reconstruction: CCS positioning",
    "rcTail": "reconstruction: boundary updates, reductions, report", "cfInit": "flux: counters, cell boxes",
    "cfSurfaces": "flux: buildCellSurfaces", "cfVelSetup": "flux: point velocities, tet cache, octree",
    "cfExch1": "flux: exchange 1 (cell surfaces)", "cfWisp": "flux: wisp surfaces", "cfCand": "flux: candidate cells",
    "cfExch2": "flux: exchange 2 (second ring)", "cfSweptExch": "flux: swept exchange (remote ring cells)",
    "fTrace": "face: trace to the pre-image (tets)", "fPre": "face: ring, planes, pre-test",
    "fVol": "face: signed-tet clips against the cell planes (volume)", "fArea": "face: clips of the pre-image (area)",
    "fAssembly": "face: flux scaling, piece sums, alpha_f update (remainder of sweptFaceCellsFlux)",
    "fOverhead": "face loop: surfaceOf, velocity interpolation, face set (remainder of the loop)",
    "cfProcSync": "flux: processor-patch sync", "cfTail": "flux: surface defect diagnostic",
    "alphaUpdate": "alpha update, mass sums", "limit": "limiter (limitFluxes, boundFlux, scale factors)",
    "brute": "brute-force bounding", "defaultAf": "default face fractions, wet faces", "report": "reductions of the reports"}

def pieces_case(case, label, dbrow):
    steps, co, ex, ni, pl = read_log(case + "/log.davofAdvectorFoam")
    recs = [parse_prof_line(l) for l in pl]
    if len(recs) < 3: return None
    use = recs[1:]
    R = recs[0]["ranks"]
    T = {k: mean([r["T"][k] for r in use]) for k in use[0]["T"]}
    TM = {k: mean([r["TM"][k] for r in use]) for k in use[0]["TM"]}
    Tn = {k: mean([r["Tn"][k] for r in use if k in r["Tn"]]) for k in use[0]["Tn"]}
    N = {k: mean([r["N"][k] for r in use]) for k in use[0]["N"]}
    NM = {k: mean([r["NM"][k] for r in use if k in r["NM"]]) for k in use[0]["NM"]}
    # derived leaves
    for D in (T, TM, Tn):
        if not D: continue
        D["fAssembly"] = max(D["fFlux"] - D["fTrace"] - D["fPre"] - D["fVol"] - D["fArea"], 0.0)
        D["fOverhead"] = max(D["cfFaceLoop"] - D["fFlux"], 0.0)
    step = T["between"] + T["advect"]
    stepM = TM["between"] + TM["advect"]
    # between at write steps vs others
    btw = sorted(r["T"]["between"] for r in use)
    return dict(label=label, case=case, steps=steps, co=co, ranks=R, T=T, TM=TM, Tn=Tn, N=N, NM=NM, step=step, stepM=stepM,
                exec=ex, nint=mean(ni), between_median=btw[len(btw)//2], first=recs[0], db=dbrow, nrecs=len(use))

def fmt(x, p=3):
    return ("%." + str(p) + "g") % x

def pieces_tables(c):
    T, TM, N = c["T"], c["TM"], c["N"]
    R = c["ranks"]; step = c["step"]
    out = []
    out.append("#### %s: %s, %d ranks, steps %d (profiled %d), max OF Co %.3f, mean interface cells %.0f, wall step (between + advect) %.4g s mean, %.4g s max over ranks"
               % (c["label"], os.path.basename(os.path.dirname(c["case"])) + "/" + os.path.basename(c["case"]), R, c["steps"], c["nrecs"], c["co"], N["interface"], step, c["stepM"]))
    out.append("")
    Tn = c["Tn"]
    out.append("| phase | what | s/step mean over ranks | share of step by the mean (%) | s/step max over ranks | share of step by the max (%) | s/step min over ranks | core-us per interface cell per step (mean) |")
    out.append("|---|---|---|---|---|---|---|---|")
    for k in LEAVES:
        if T[k] <= 0 and TM[k] <= 0: continue
        out.append("| %s | %s | %s | %.1f | %s | %.1f | %s | %.2f |" % (k, DESC[k], fmt(T[k]), 100*T[k]/step, fmt(TM[k]), 100*TM[k]/c["stepM"],
                                                                 fmt(Tn[k]) if Tn else "-", 1e6*R*T[k]/max(N["interface"], 1)))
    if Tn:
        out.append("")
        out.append("Imbalance and wait (s/step): load imbalance of the face loop max/mean = %.2f (flux time) and %.2f (faces through the flux per rank, counts); "
                   "wait of the average rank in the processor-patch sync (mean - min of the sync, the minimum being the exchange itself on the most loaded rank) = %s s (%.1f %% of the step); exchange itself (min over ranks of the sync) = %s s."
                   % (TM["fFlux"]/max(T["fFlux"], 1e-30), (c["NM"].get("facesFlux", 0)*R/max(N["facesFlux"], 1)),
                      fmt(T["cfProcSync"] - Tn["cfProcSync"]), 100*(T["cfProcSync"] - Tn["cfProcSync"])/step, fmt(Tn["cfProcSync"])))
    sub = lambda ks: sum(T[k] for k in ks)
    out.append("")
    out.append("Aggregates (mean over ranks, s/step): advect() %s (%.1f %%), reconstruction %s (%.1f %%), swept flux curvedFlux %s (%.1f %%), "
               "face loop %s (%.1f %%), parallel exchanges (exch1 + exch2 + swept exchange + proc sync) %s (%.1f %%; max over ranks %s), between %s (%.1f %%; median %s)."
               % (fmt(T["advect"]), 100*T["advect"]/step, fmt(T["recon"]), 100*T["recon"]/step, fmt(T["curved"]), 100*T["curved"]/step,
                  fmt(T["cfFaceLoop"]), 100*T["cfFaceLoop"]/step,
                  fmt(sub(["cfExch1", "cfExch2", "cfSweptExch", "cfProcSync"])), 100*sub(["cfExch1", "cfExch2", "cfSweptExch", "cfProcSync"])/step,
                  fmt(TM["cfExch1"] + TM["cfExch2"] + TM["cfSweptExch"] + TM["cfProcSync"]),
                  fmt(T["between"]), 100*T["between"]/step, fmt(c["between_median"])))
    out.append("")
    f = max(N["facesFlux"], 1); sw = max(N["facesSweep"], 1)
    out.append("Counts per step (sums over ranks): interface cells %.0f; CCS cuts %.0f (%.2f per interface cell) + explicit-plane cuts %.0f; faces swept (processFace) %.0f, through the flux %.0f, with a swept volume %.0f; "
               "swept-volume triangles %.0f (%.2f per face); ring cells %.0f (%.2f per face); cells kept by the pre-test (planes the swept volume is clipped against) %.0f (%.2f per face); "
               "signed tets built and clipped %.0f (%.2f per face); tet-plane clips %.0f (%.2f per face, %.2f per tet); polygon clips of the tet faces %.0f (%.2f per face); "
               "coverage-polygon clips %.0f (%.2f per face); pre-image (area) clips %.0f (%.2f per face); pieces produced %.0f (%.2f per face); limiter cells %.0f; exchange rounds %.1f"
               % (N["interface"], N["ccsCuts"], N["ccsCuts"]/max(N["interface"], 1), N["explicitCuts"], N["faces"], N["facesFlux"], N["facesSweep"],
                  N["tets"], N["tets"]/f, N["ring"], N["ring"]/f, N["touched"], N["touched"]/f, N["tetsBuilt"], N["tetsBuilt"]/f,
                  N["tetClips"], N["tetClips"]/f, N["tetClips"]/max(N["tetsBuilt"], 1), N["polyClips"], N["polyClips"]/f,
                  N["covClips"], N["covClips"]/f, N["areaClips"], N["areaClips"]/f, N["pieces"], N["pieces"]/f, N["limitCells"], N["exch"]))
    out.append("")
    return "\n".join(out)

def of_case(case, label, dbrow):
    steps, co, ex, ni, _ = read_log(case + "/log.davofAdvectorFoam")
    P = prof.case_profile(case)
    if not P: return None
    def get(sub):
        for d, v in P.items():
            if sub in d: return v
        return (0, 0.0, 0.0, "-")
    run = [v for d, v in P.items() if d.startswith("time.run()")]
    run = run[0] if run else (0, 0, 0, "-")
    return dict(label=label, case=case, steps=steps, co=co, P=P, run=run, get=get, exec=ex, nint=mean(ni), db=dbrow)

def of_row(c):
    g = c["get"]; run = c["run"]; n = max(c["steps"], 1)
    adv = g("isoAdvection::advect"); rec = g("::reconstruct(bool)")
    tif = g("timeIntegratedFlux"); lim = g("isoAdvection::limitFluxes"); fo_e = g("functionObjects.execute()")
    out = dict(step=run[1]/n, adv=adv[1]/n, rec=rec[1]/n, flux=tif[1]/n, lim=lim[1]/n, fo=fo_e[1]/n, runmax=run[2]/n, advmax=adv[2]/n,
               recmax=rec[2]/n)
    out["other"] = out["step"] - out["adv"] - out["fo"]
    return out

def main():
    cases = []
    for st in (sorted(glob.glob(root + "/diag2_prof_*")) + sorted(glob.glob(root + "/diag2_prof2_*"))
               + sorted(glob.glob(root + "/diag2_scal_*")) + sorted(glob.glob(root + "/diag2_scal2_*"))):
        name = os.path.basename(st)
        rows = db(name)
        for cdir in sorted(glob.glob(st + "/iso*_0*")):
            if not os.path.exists(cdir + "/log.davofAdvectorFoam"): continue
            idx = cdir.split("_")[-1]
            r = rows.get(idx, {})
            row = r.get("TPF_ROW", "?")
            label = name.replace("diag2_", "")
            if "_davof" in name or name.startswith("diag2_scal"):
                c = pieces_case(cdir, label, r)
                kind = "pieces"
            else:
                c = of_case(cdir, label, r)
                kind = "of"
            if c: cases.append((kind, name, idx, r, c))
    # --- per-case phase tables of the pieces rows
    print("## Pieces rows: time per step per phase\n")
    for kind, name, idx, r, c in cases:
        if kind == "pieces" and ("prof2_" in name or "scal2_" in name):
            print(pieces_tables(c))
    # --- plicRDF and isoAlpha rows (OpenFOAM profiling)
    print("## plicRDF and isoAlpha rows: OpenFOAM profiling (per step = totalTime / steps; mean over ranks)\n")
    print("| case | row | ranks | steps | max OF Co | s/step main loop (time.run()) | advect() | reconstruct() | timeIntegratedFlux | limitFluxes | function objects (execute) | outside advect and FOs | main loop max over ranks | wall (s) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for kind, name, idx, r, c in cases:
        if kind == "of":
            o = of_row(c)
            print("| %s %s | %s | %s | %d | %.3f | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                name.replace("diag2_prof_", ""), idx, r.get("TPF_ROW", "?"), r.get("np", "?"), c["steps"], c["co"],
                fmt(o["step"]), fmt(o["adv"]), fmt(o["rec"]), fmt(o["flux"]), fmt(o["lim"]), fmt(o["fo"]), fmt(o["other"]), fmt(o["runmax"]), r.get("wallClockSolve", "?")))
    # --- summary across rows
    print("\n## Summary: seconds per step (mean over ranks), the three rows side by side\n")
    print("| family N (ranks) | OF Co | window | steps | interface cells | pieces: step | pieces: reconstruction | pieces: advect minus reconstruction | pieces: between | plicRDF: step | plicRDF: reconstruct | plicRDF: advect minus reconstruct | isoAlpha: step | isoAlpha: reconstruct | isoAlpha: advect minus reconstruct | step ratio pieces/plicRDF | pieces/isoAlpha | reconstruction ratio pieces/plicRDF | advection ratio pieces/plicRDF | pieces core-us per interface cell per step | plicRDF | isoAlpha |")
    print("|" + "---|"*22)
    groups = {}
    for kind, name, idx, r, c in cases:
        if name.startswith("diag2_scal"): continue
        fam = name.replace("diag2_prof2_", "").replace("diag2_prof_", "").rsplit("_", 1)[0]
        key = (fam, str(r.get("TARGET_ADVECTIVE_CO")))
        sub = kind + ("_" + r.get("TPF_ROW", "") if kind == "of" else ("2" if "prof2_" in name else "1"))
        groups.setdefault(key, {})[sub] = (c, r)
    for (fam, tgt), g in sorted(groups.items()):
        pc = g.get("pieces2") or g.get("pieces1"); pl = g.get("of_plicRDF"); isoc = g.get("of_isoAlpha")
        if not pc: continue
        c, r = pc
        T = c["T"]; N = c["N"]
        ps = c["step"]; prec = T["recon"]; padv = T["advect"] - T["recon"]
        def ofv(x):
            if not x: return (float("nan"),)*3
            o = of_row(x[0]); return (o["step"], o["rec"], o["adv"] - o["rec"])
        pls = ofv(pl); iss = ofv(isoc)
        R = c["ranks"]
        cu = lambda s_, n_: 1e6*R*s_/max(n_, 1)
        ni = N["interface"]
        print("| %s (%d) | %s | END_TIME %s | %d | %.0f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f |" % (
            fam, R, c["co"] and ("%.3f" % c["co"]), r.get("END_TIME"), c["steps"], ni,
            fmt(ps), fmt(prec), fmt(padv), fmt(T["between"]),
            fmt(pls[0]), fmt(pls[1]), fmt(pls[2]), fmt(iss[0]), fmt(iss[1]), fmt(iss[2]),
            ps/pls[0] if pls[0] == pls[0] else float("nan"), ps/iss[0] if iss[0] == iss[0] else float("nan"),
            prec/pls[1] if pls[1] == pls[1] and pls[1] > 0 else float("nan"), padv/pls[2] if pls[2] == pls[2] and pls[2] > 0 else float("nan"),
            cu(ps, ni), cu(pls[0], ni) if pls[0] == pls[0] else float("nan"), cu(iss[0], ni) if iss[0] == iss[0] else float("nan")))
    print("\n## Repetition: the pieces step (between + advect, mean over ranks, s) of the first and the second batch (identical cases, separate jobs)\n")
    print("| family N (ranks) | OF Co target | first batch | second batch | ratio second/first |")
    print("|---|---|---|---|---|")
    for (fam, tgt), g in sorted(groups.items()):
        a = g.get("pieces1"); b = g.get("pieces2")
        if a and b:
            print("| %s | %s | %s | %s | %.3f |" % (fam, tgt, fmt(a[0]["step"], 4), fmt(b[0]["step"], 4), b[0]["step"]/a[0]["step"]))
    # --- strong scaling
    print("\n## Strong scaling: deformation3D N 64, pieces, first 200 steps (steps 2..200 averaged)\n")
    print("| ranks | steps | step s (mean over ranks) | step s (max over ranks) | speedup vs 1 rank | parallel efficiency | advect | between | reconstruction | flux face loop | exchanges (exch1 + exch2 + swept + proc sync, mean) | exchanges max over ranks | exchange rounds per step | wait before the swept exchange (max) |")
    print("|" + "---|"*14)
    base = None
    sc = [(c["ranks"], c) for kind, name, idx, r, c in cases if name.startswith("diag2_scal2") and kind == "pieces"]
    for R, c in sorted(sc, key=lambda x: x[0]):
        T, TM, N = c["T"], c["TM"], c["N"]
        wall = c["step"]
        if base is None: base = wall
        ex = T["cfExch1"] + T["cfExch2"] + T["cfSweptExch"] + T["cfProcSync"]
        exm = TM["cfExch1"] + TM["cfExch2"] + TM["cfSweptExch"] + TM["cfProcSync"]
        print("| %d | %d | %s | %s | %.2f | %.2f | %s | %s | %s | %s | %s | %s | %.1f | %s |" % (
            R, c["steps"], fmt(wall), fmt(c["stepM"]), base/wall, base/(R*wall), fmt(T["advect"]), fmt(T["between"]), fmt(T["recon"]), fmt(T["cfFaceLoop"]),
            fmt(ex), fmt(exm), N["exch"], fmt(c["first"]["X"][0] if c["first"]["X"] else 0)))
    sc1 = sorted([(c["ranks"], c["step"]) for kind, name, idx, r, c in cases if name.startswith("diag2_scal_") and kind == "pieces"])
    print("\nFirst batch (identical cases, library without the rank minima): step s " + ", ".join("%d ranks %s" % (R, fmt(t, 4)) for R, t in sc1))
    # --- three largest phases
    print("\n## Three largest leaf phases per pieces case (second batch; by the maximum over the ranks, the critical path)\n")
    print("| case | ranks | 1st | 2nd | 3rd |")
    print("|---|---|---|---|---|")
    for kind, name, idx, r, c in cases:
        if kind != "pieces" or not ("prof2_" in name or "scal2_" in name): continue
        T = c["T"]; TM = c["TM"]; step = c["step"]; stepM = c["stepM"]
        topm = sorted([(TM[k], k) for k in LEAVES if k not in ("cfProcSync",)], reverse=True)[:3]
        print("| %s %s (Co %.3f) | %d | %s |" % (name.replace("diag2_", ""), idx, c["co"], c["ranks"],
              " | ".join("%s %s s max over ranks (%.1f %% of the step), mean %s s" % (k, fmt(v), 100*v/stepM, fmt(T[k])) for v, k in topm)))

if __name__ == "__main__":
    main()
