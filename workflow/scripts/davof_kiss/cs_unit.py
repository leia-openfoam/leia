import csv, glob, hashlib, os, re, sys
S = "/work/scratch/tm83tomy/davof/agents/diag1/leia/studies/davof/"
R = "/work/scratch/tm83tomy/davof/agents/ccs/leia/studies/davof/"
RE = re.compile(r"contentScale upwind: faces examined (\d+), scaled (\d+) \(negative factor (\d+)\), min factor (\S+), removed volume / V_k (\S+); unit check, faces with L_U and V_D\(U\) above 1e-12 V_U: (\d+), raw ratio L_U/V_D\(U\) min (\S+) max (\S+), of them scaled by more than 1e-12: (\d+)")

def lines(p):
    return [l.rstrip() for l in open(p, errors="replace") if "After conservative bounding" in l or "mass change relative" in l or "limiter faces changed" in l]

for arm in ("normal", "oblique"):
    st = f"cs_unit_plane_{arm}"
    rows = list(csv.DictReader(open(S + st + "/" + st + "_database.csv")))
    old = list(csv.DictReader(open(R + f"ccs_plane_{arm}_np1/ccs_plane_{arm}_np1_database.csv")))
    print(f"\n### {arm} arm, np1")
    print("| N | target Co | setting | steps | faces examined (sum) | scaled (sum) | negative | min factor | big class faces | raw ratio min | raw ratio max | big-class faces scaled by > 1e-12 | limiter faces | removed volume/V_k (sum) | L1(alpha_f) | key-off vs Stage B |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        cd = S + st + "/" + [d for d in os.listdir(S + st) if d.endswith("_" + r["index"])][0]
        txt = open(cd + "/log.davofAdvectorFoam", errors="replace").read()
        cs = [m.groups() for m in RE.finditer(txt)]
        isCS = "sweptContentScale upwind" in " ".join(str(v) for v in r.values())
        steps = len(re.findall("After conservative bounding", txt))
        lim = sum(int(a) for a, b, c in re.findall(r"limiter faces changed (\d+), cells (\d+), volume moved (\S+);", txt))
        l1 = ""
        cp = cd + "/advectorVoF.csv"
        if os.path.isfile(cp) and os.path.getsize(cp):
            rr = list(csv.DictReader(open(cp))); l1 = rr[-1].get("L1_alphaf", "") if rr else ""
            l1 = "%.2e" % max(float(x["L1_alphaf"]) for x in rr)
        ident = ""
        if not isCS:
            o = [x for x in old if x["N_CELLS"] == r["N_CELLS"] and abs(float(x["TARGET_ADVECTIVE_CO"]) - float(r["TARGET_ADVECTIVE_CO"])) < 1e-9 and "sweptClip onePlane" in " ".join(x.values())]
            if o:
                od = R + f"ccs_plane_{arm}_np1/" + [d for d in os.listdir(R + f"ccs_plane_{arm}_np1") if d.endswith("_" + o[0]["index"])][0]
                a = hashlib.md5("\n".join(lines(cd + "/log.davofAdvectorFoam")).encode()).hexdigest()
                b = hashlib.md5("\n".join(lines(od + "/log.davofAdvectorFoam")).encode()).hexdigest()
                ident = "IDENTICAL" if a == b else "DIFFERENT"
        if cs:
            ne = sum(int(c[0]) for c in cs); ns = sum(int(c[1]) for c in cs); nn = sum(int(c[2]) for c in cs)
            mf = min(float(c[3]) for c in cs); rem = sum(float(c[4]) for c in cs); nb = sum(int(c[5]) for c in cs)
            rmin = min(float(c[6]) for c in cs); rmax = max(float(c[7]) for c in cs); nbo = sum(int(c[8]) for c in cs)
            print(f"| {r['N_CELLS']} | {r['TARGET_ADVECTIVE_CO']} | S1+contentScale | {steps} | {ne} | {ns} | {nn} | {mf:.6g} | {nb} | {rmin:.15g} | {rmax:.15g} | {nbo} | {lim} | {rem:.3g} | {l1} | |")
        else:
            print(f"| {r['N_CELLS']} | {r['TARGET_ADVECTIVE_CO']} | onePlane (key off) | {steps} | | | | | | | | | {lim} | | {l1} | {ident} |")
