#!/usr/bin/env python3
"""Tables of the switch experiment step 2 (2026-11-16; docs/davof/REPORT-swept.md, "One plane per face: the switch
experiment"). usage: sw_an.py <study> [<study> ...]   (studies under studies/davof). Reading and printing only."""
import csv, glob, os, re, sys
S = "/work/scratch/tm83tomy/davof/agents/diag1/leia/studies/davof/"
STAGEC = "/work/scratch/tm83tomy/davof/agents/diag2/stageC_tables.md"
VAR = ["base", "S1", "S2low", "S2all", "S1+S2low", "S1+S2all"]
F = r"([-+0-9.eE]+|nan|inf)"
RE_AFTER = re.compile(r"After conservative bounding: min\(alpha\) = " + F + r", max\(alpha\) = 1 \+ " + F)
RE_SW = re.compile(r"switches acted on faces: S1 (\d+), S2 (\d+) \(upwind alpha < 0.5: (\d+), >= 0.5: (\d+)\)")
RE_LIM = re.compile(r"one plane: limiter faces changed (\d+), cells (\d+), volume moved " + F)
RE_DSTEP = re.compile(r"DIAG step (\d+) t \S+ : unlimited update: cells outside \[-1e-10, 1\+1e-10\] (\d+), non-finite (\d+), min alpha " + F + r", max alpha - 1 " + F + r", summed violation \(cell volumes\) " + F + r", cells with violation > 1e-3 (\d+)")
RE_DLIM = re.compile(r"DIAG step (\d+) : limiter passes (\d+), faces changed (\d+), cells touched (\d+), volume moved \(cell volumes\) " + F)


def stored(fam, N):
    out = {}
    if not os.path.isfile(STAGEC):
        return out
    sec = None
    for l in open(STAGEC):
        m = re.match(r"### (\w+), OpenFOAM Co 0.2", l)
        if m:
            sec = m.group(1); continue
        if l.startswith("### "):
            sec = None
        if sec == fam and l.startswith("| ") :
            c = [x.strip() for x in l.strip().strip("|").split("|")]
            if len(c) > 17 and c[0] in ("pieces", "plicRDF") and c[1] == str(N):
                out[c[0]] = dict(E=c[6], wall=c[16], ranks=c[2], status=c[5])
    return out


def main():
    for study in sys.argv[1:]:
        fam = study.split("_")[1]
        N = int(re.search(r"N(\d+)", study).group(1))
        ef = f"{S}{study}/{study}_errors.csv"
        errs = {r["caseDir"]: r for r in csv.DictReader(open(ef))} if os.path.isfile(ef) else {}
        dbs = {r["index"]: r for r in csv.DictReader(open(f"{S}{study}/{study}_database.csv"))}
        print(f"\n### {study} (OF Co 0.2; diagnostics on; 6 cases concurrently on 6 cores)\n")
        print("| variant | status | steps | first step with violation > 1e-3 V_k | run max violation before limiting (V_k) | run max violation after limiting (V_k) | S1 faces/step mean (max) | S2 faces/step mean (max) [alpha_up<0.5 / >=0.5 mean] | limiter faces/step mean (max) | limiter cells/step mean (max) | limiter volume total (cell volumes) | limiter passes max | E_shape | Emass (csv) | wall clock (s) | ranks |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for k, v in enumerate(VAR):
            cd = sorted(d for d in os.listdir(S + study) if d.endswith("_%05d" % k))
            if not cd:
                print(f"| {v} | no case | | | | | | | | | | | | | | |"); continue
            cdir = S + study + "/" + cd[0]
            lp = cdir + "/log.davofAdvectorFoam"
            log = open(lp, errors="replace").read() if os.path.isfile(lp) else ""
            dsteps = [m for m in RE_DSTEP.finditer(log)]
            steps = len(dsteps)
            first3 = next((int(m.group(1)) for m in dsteps if int(m.group(7)) > 0), None)
            vb = max([max(-float(m.group(4)), float(m.group(5))) for m in dsteps] or [0])
            va = max([max(-float(m.group(1)), float(m.group(2))) for m in RE_AFTER.finditer(log)] or [0])
            sw = [tuple(int(x) for x in m.groups()) for m in RE_SW.finditer(log)]
            s1 = [a for a, b, c, d in sw]; s2 = [b for a, b, c, d in sw]; s2l = [c for a, b, c, d in sw]; s2h = [d for a, b, c, d in sw]
            mean = lambda x: (sum(x)/len(x)) if x else 0
            lim = [(int(m.group(1)), int(m.group(2))) for m in RE_LIM.finditer(log)]
            dl = [(int(m.group(2)), float(m.group(5))) for m in RE_DLIM.finditer(log)]
            lf = [a for a, b in lim]; lc = [b for a, b in lim]
            status = "completed"
            mm = re.search(r"ABORT at step (\d+)", log)
            if mm: status = "abort step " + mm.group(1)
            elif "sigFpe" in log or "Floating point exception" in log.replace("trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).", ""):
                status = "SIGFPE (last DIAG step %d)" % (steps)
            elif not re.search(r"\nEnd\n", log): status = "not finished"
            er = errs.get(cd[0], {})
            E = er.get("shapeError", "n/a"); wall = er.get("wallClockSolve", "n/a"); nr = er.get("nRanks", "n/a")
            emass = "n/a"
            cp = cdir + "/advectorVoF.csv"
            if os.path.isfile(cp) and os.path.getsize(cp):
                rows = list(csv.DictReader(open(cp)))
                if rows:
                    emass = rows[-1].get("Emass", "n/a")
                    if status == "completed": E = rows[-1].get("Eshape", "n/a")
            def fe(x):
                try: return "%.3g" % float(x)
                except Exception: return x
            sdisp = lambda x: f"{mean(x):.3g} ({max(x) if x else 0})"
            print(f"| {v} | {status} | {steps} | {first3 if first3 else 'none'} | {vb:.3g} | {va:.3g} | {sdisp(s1)} | {sdisp(s2)} [{mean(s2l):.3g} / {mean(s2h):.3g}] | {sdisp(lf)} | {sdisp(lc)} | "
                  f"{fe(sum(b for a, b in dl)) if dl else 'n/a'} | {max([a for a, b in dl] or [0])} | {fe(E)} | {fe(emass)} | {fe(wall)} | {nr} |")
        st = stored(fam, N)
        for m in ("pieces", "plicRDF"):
            if m in st:
                print(f"\n{m} stored (Co 0.2, N {N}, ranks {st[m]['ranks']}): status {st[m]['status']}, E_shape {st[m]['E']}, wall clock {st[m]['wall']} s")


if __name__ == "__main__":
    main()
