#!/usr/bin/env python3
"""Stage B one-plane matrix of the best switch variant (2026-11-17): tables as stageB_tables.md with orders over N and the
ratio to plicRDF and to the pieces row (E_shape, wall clock). The pieces / plicRDF / isoAlpha rows are read from the
stored table /work/scratch/tm83tomy/davof/agents/diag2/stageC_tables.md. usage: sm_an.py <variantLabel>"""
import csv, math, os, re, sys, glob
S = "/work/scratch/tm83tomy/davof/agents/diag1/leia/studies/davof/"
STAGEC = "/work/scratch/tm83tomy/davof/agents/diag2/stageC_tables.md"
FAMS = {"shear2D": (64, 128, 256), "shear3D": (32, 64, 128), "deformation3D": (32, 64, 128)}
label = sys.argv[1] if len(sys.argv) > 1 else "onePlane S1+S2all"


def num(x):
    try:
        v = float(x); return v if v == v else None
    except Exception:
        return None


def stored():
    out = {}
    sec = None; hdr = None
    for l in open(STAGEC):
        m = re.match(r"### (\w+), OpenFOAM Co ([0-9.]+)", l)
        if m:
            sec = (m.group(1), float(m.group(2))); hdr = None; continue
        if l.startswith("### "):
            sec = None; continue
        if sec and l.startswith("| method"):
            hdr = [x.strip() for x in l.strip().strip("|").split("|")]
            if "N" not in hdr: hdr = None
            continue
        if sec and hdr and l.startswith("| ") and not l.startswith("|---"):
            c = [x.strip() for x in l.strip().strip("|").split("|")]
            if len(c) == len(hdr):
                r = dict(zip(hdr, c))
                out[(sec[0], sec[1], r["method"], int(r["N"]))] = r
    return out


def order(pts):
    pts = [(math.log(h), math.log(e)) for h, e in pts if e and e > 0]
    if len(pts) < 2:
        return None
    n = len(pts); mx = sum(p[0] for p in pts)/n; my = sum(p[1] for p in pts)/n
    sxx = sum((p[0]-mx)**2 for p in pts)
    return sum((p[0]-mx)*(p[1]-my) for p in pts)/sxx if sxx else None


def case_rows():
    rows = []
    for fam, Ns in FAMS.items():
        for N in Ns:
            st = f"swm_{fam}_N{N}"
            dbf = f"{S}{st}/{st}_database.csv"
            if not os.path.isfile(dbf):
                continue
            db = list(csv.DictReader(open(dbf)))
            tg = sorted({float(r["TARGET_ADVECTIVE_CO"]) for r in db})
            for k, r in enumerate(db):
                cd = sorted(glob.glob(f"{S}{st}/*_{r['index']}"))
                if not cd:
                    continue
                c = cd[0]
                lp = c + "/log.davofAdvectorFoam"
                log = open(lp, errors="replace").read() if os.path.isfile(lp) else ""
                status = "completed" if re.search(r"\nEnd\n", log[-3000:]) else "unfinished"
                mm = re.search(r"ABORT at step (\d+)", log)
                if mm: status = "aborted step " + mm.group(1)
                elif "sigFpe::sigHandler" in log: status = "SIGFPE"
                csvp = c + "/advectorVoF.csv"
                last = {}
                if os.path.isfile(csvp) and os.path.getsize(csvp):
                    rr = list(csv.DictReader(open(csvp)))
                    last = rr[-1] if rr else {}
                cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: ([0-9.eE+-]+)", log) if re.match(r"^[0-9.eE+-]+$", x)]
                ef = f"{S}{st}/{st}_errors.csv"
                wall = None; nr = None
                if os.path.isfile(ef):
                    for e in csv.DictReader(open(ef)):
                        if e.get("caseDir", "") == os.path.basename(c):
                            wall = num(e.get("wallClockSolve")); nr = e.get("nRanks")
                lim = re.findall(r"limiter faces changed (\d+), cells (\d+), volume moved (\S+);", log)
                sw = re.findall(r"switches acted on faces: S1 (\d+), S2 (\d+)", log)
                nst = len(re.findall(r"After conservative bounding", log))
                rows.append(dict(fam=fam, N=N, target=float(r["TARGET_ADVECTIVE_CO"]), ofCo=(max(cos) if cos else None),
                                 np=r["np"], status=status, steps=nst, E=(num(last.get("Eshape")) if status == "completed" else None),
                                 Emass=num(last.get("Emass")), minA=num(last.get("minAlpha")), wall=wall, nr=nr,
                                 limFaces=sum(int(a) for a, b, c_ in lim), limVol=sum(float(c_) for a, b, c_ in lim if num(c_) is not None),
                                 s2=(sum(int(b) for a, b in sw)/len(sw) if sw else 0)))
    return rows


def main():
    st = stored()
    rows = case_rows()
    for fam in FAMS:
        cos_sorted = sorted({r["target"] for r in rows if r["fam"] == fam})
        for k, tgt in enumerate(cos_sorted):
            rs = sorted([r for r in rows if r["fam"] == fam and r["target"] == tgt], key=lambda r: r["N"])
            if not rs:
                continue
            ofco = max([r["ofCo"] for r in rs if r["ofCo"]] or [0])
            tag = 0.2 if k == 0 else 0.5 if k == 1 else 0.95
            print(f"\n### {fam}, OpenFOAM Co {tag} (target {tgt}; measured max Co {ofco:.3f})\n")
            print("| method | N | ranks | steps | status | E_shape | E_shape / plicRDF | E_shape / pieces | Emass (csv) | min alpha | limiter faces (run) | limiter volume (run, V) | S2 faces/step mean | wall clock (s) | wall / plicRDF | wall / pieces |")
            print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for r in rs:
                pl = st.get((fam, tag, "plicRDF", r["N"])); pc = st.get((fam, tag, "pieces", r["N"]))
                f = lambda x: ("%.3g" % x) if x is not None else "n/a"
                plE = num(pl["E_shape"]) if pl else None; pcE = num(pc["E_shape"]) if pc else None
                plW = num(pl["wall clock (s)"]) if pl else None; pcW = num(pc["wall clock (s)"]) if pc else None
                print(f"| {label} | {r['N']} | {r['nr'] or r['np']} | {r['steps']} | {r['status']} | {f(r['E'])} | {f(r['E']/plE) if r['E'] and plE else 'n/a'} | {f(r['E']/pcE) if r['E'] and pcE else 'n/a'} | {f(r['Emass'])} | {f(r['minA'])} | {r['limFaces']} | {f(r['limVol'])} | {f(r['s2'])} | {f(r['wall'])} | {f(r['wall']/plW) if r['wall'] and plW else 'n/a'} | {f(r['wall']/pcW) if r['wall'] and pcW else 'n/a'} |")
                for name, d in (("pieces", pc), ("plicRDF", pl)):
                    if d:
                        print(f"| {name} (stored) | {r['N']} | {d['ranks']} | {d['steps']} | {d['status']} | {d['E_shape']} | | | | | | | | {d['wall clock (s)']} | | |")
            h = lambda N: 1.0/N
            for name, get in (("S1+S2all", lambda r: r["E"]), ("pieces", None), ("plicRDF", None)):
                if get:
                    pts = [(h(r["N"]), r["E"]) for r in rs if r["E"]]
                else:
                    pts = [(h(N), num(st[(fam, tag, name, N)]["E_shape"])) for N in FAMS[fam] if (fam, tag, name, N) in st and st[(fam, tag, name, N)]["status"] == "completed"]
                o = order(pts)
                print(f"\norder of E_shape over N, {name}: {('%.2f' % o) if o is not None else 'n/a'} ({len(pts)} completed rungs)")


if __name__ == "__main__":
    main()
