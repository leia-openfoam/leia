#!/usr/bin/env python3
"""Tables and verdicts of the pointIsoAdvector study (docs/davof/REPORT-swept.md "pointIsoAdvector").
Reads studies/davof/pointIso_<family>_N<N> (row pointIso) and the stored rows of the Stage B table
(ccs_tables/stageB_cells.csv: plicRDF, isoAlpha, pieces, onePlane), writes one markdown file.
    python3 workflow/scripts/davof_pointiso/pointiso_an.py [--out FILE]
"""
import argparse, csv, glob, math, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "davof_kiss"))
import kiss_an as ka
import ccs_stageB_an as sb
REPO = ka.REPO
FAMS = {"shear2D": (64, 128, 256), "deformation3D": (32, 64, 128)}
COS = (0.2, 0.5, 0.95)
STAGEB = "/work/scratch/tm83tomy/davof/agents/ccs/leia/studies/davof/ccs_tables/stageB_cells.csv"
PL_RE = re.compile(r"pointIsoAdvection: planes (\d+), fit residual \(rms/V\^\(1/3\)\) mean (\S+) max (\S+), max angle\(n', n_k\) (\S+) deg, "
                   r"cells with boundary points (\d+), with more than three non-collinear boundary points (\d+)")
def f(v, nd=2):
    return "-" if v is None else "%.*e" % (nd, v)
def nearest_co(x):
    return min(COS, key=lambda c: abs(c - x))
def planes(c):
    p = os.path.join(c, "log.davofAdvectorFoam")
    if not os.path.isfile(p):
        return {}
    m = PL_RE.findall(open(p, errors="replace").read())
    if not m:
        return {}
    n = [int(x[0]) for x in m]
    tot = sum(n)
    return dict(steps=len(m), planesMean=tot/len(m), planesMax=max(n),
                resMean=sum(int(x[0])*float(x[1]) for x in m)/max(tot, 1), resMax=max(float(x[2]) for x in m),
                angMax=max(float(x[3]) for x in m), bndCells=max(int(x[4]) for x in m), over3=max(int(x[5]) for x in m))
def mine():
    out = []
    for fam, Ns in FAMS.items():
        for N in Ns:
            for c in sorted(glob.glob(os.path.join(REPO, "studies", "davof", "pointIso_%s_N%d" % (fam, N), "*_[0-9][0-9][0-9][0-9][0-9]"))):
                r = ka.read_case(c)
                if not r:
                    continue
                sb.extra(c, r)
                r["fam"] = fam; r["method"] = "pointIso"
                r["coN"] = nearest_co(r["ofCo"] if r["ofCo"] else r["co"])
                r["pl"] = planes(c)
                out.append(r)
    return out
def stored():
    out = []
    for r in csv.DictReader(open(STAGEB)):
        if r["family"] not in FAMS:
            continue
        g = lambda k: sb.num(r.get(k))
        out.append(dict(fam=r["family"], method=r["method"], N=int(r["N"]), coN=float(r["co"]), np=int(float(r["np"])), ofCo=g("ofCo"),
                        steps=int(float(r["steps"])) if r["steps"] else None, status=r["status"], Eshape=g("Eshape"), Emass=g("Emass_pct"),
                        minAlpha=g("minAlpha"), excess=g("excess"), wall=g("wall"), src=r["src"]))
    return out
def status_of(r):
    s = r.get("status", "")
    return "completed" if s in ("ended", "completed") else s
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "studies", "davof", "pointiso_tables.md"))
    a = ap.parse_args()
    M = mine(); S = stored()
    by = {}
    for r in S:
        by[(r["fam"], r["method"], r["N"], r["coN"])] = r
    for r in M:
        by[(r["fam"], "pointIso", r["N"], r["coN"])] = r
    L = []
    verd = dict(done=0, n=0, pp2=[], pp3=[], pp5=[])
    orders = {}
    for fam, Ns in FAMS.items():
        for co in COS:
            L.append("\n### %s, OpenFOAM Co %s\n" % (fam, co))
            L.append("| method | N | ranks | OF Co | steps | status | E_shape | E_mass (%) | min alpha | max alpha - 1 | wall clock (s) | E_shape / plicRDF | wall / plicRDF | source |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for meth in ("pointIso", "plicRDF", "isoAlpha", "pieces"):
                for N in Ns:
                    r = by.get((fam, meth, N, co))
                    if r is None:
                        continue
                    pl = by.get((fam, "plicRDF", N, co))
                    rat = (r["Eshape"]/pl["Eshape"]) if (r.get("Eshape") and pl and pl.get("Eshape")) else None
                    wr = (r["wall"]/pl["wall"]) if (r.get("wall") and pl and pl.get("wall")) else None
                    npr = r.get("np")
                    L.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                        meth, N, npr, "%.3f" % r["ofCo"] if r.get("ofCo") else "-", r.get("steps"), status_of(r)[:60],
                        f(r.get("Eshape")), f(r.get("Emass") if "Emass" in r and meth != "pointIso" else r.get("Emass")),
                        f(r.get("minAlpha")), f(r.get("excess")), "%.0f" % r["wall"] if r.get("wall") else "-",
                        "%.3f" % rat if (meth == "pointIso" and rat) else "-", "%.3f" % wr if (meth == "pointIso" and wr) else "-",
                        r.get("src", r.get("study", ""))))
                    if meth == "pointIso":
                        verd["n"] += 1
                        if status_of(r) == "completed":
                            verd["done"] += 1
                        if pl and status_of(r) == "completed":
                            verd["pp3"].append((fam, N, co, rat))
                            verd["pp5"].append((fam, N, co, wr))
                            for k in ("minAlpha", "excess"):
                                x, y = r.get(k), pl.get(k)
                                if x is None or y is None:
                                    continue
                                ok = (abs(x) <= 1e-12 and abs(y) <= 1e-12) or (abs(y) > 0 and 0.5 <= abs(x)/abs(y) <= 2.0)
                                verd["pp2"].append((fam, N, co, k, x, y, ok))
            L.append("\nOrders of E_shape over N (least squares; pairwise in brackets), completed rungs only:\n")
            L.append("| method | E_shape order |")
            L.append("|---|---|")
            for meth in ("pointIso", "plicRDF", "isoAlpha", "pieces"):
                pts = [(N, by[(fam, meth, N, co)].get("Eshape")) for N in Ns if (fam, meth, N, co) in by and status_of(by[(fam, meth, N, co)]) == "completed"]
                ls, pw = sb.fit(pts)
                orders[(fam, co, meth)] = ls
                L.append("| %s | %s |" % (meth, "-" if ls is None else "%.2f (%s)" % (ls, ", ".join("%.2f" % x for x in pw))))
    L.append("\n### Fitted planes per case (pointIso)\n")
    L.append("| family | N | OF Co | steps with the line | planes/step mean (max) | fit residual mean / max (rms/V^(1/3)) | max angle(n', n_k) (deg) | max cells with boundary points | max cells with >3 non-collinear |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for r in sorted(M, key=lambda r: (r["fam"], r["coN"], r["N"])):
        p = r["pl"]
        if p:
            L.append("| %s | %d | %.3f | %d | %.1f (%d) | %s / %s | %s | %d | %d |" % (r["fam"], r["N"], r["ofCo"] or 0, p["steps"], p["planesMean"], p["planesMax"], f(p["resMean"]), f(p["resMax"]), f(p["angMax"]), p["bndCells"], p["over3"]))
    L.append("\n### Verdicts\n")
    L.append("- PP1: %d of %d completed (abort or crash is the falsifier)." % (verd["done"], verd["n"]))
    bad2 = [x for x in verd["pp2"] if not x[6]]
    L.append("- PP2: bounds within a factor 2 of plicRDF's: %d of %d comparisons within; outside: %s" % (len(verd["pp2"]) - len(bad2), len(verd["pp2"]), "; ".join("%s N%d Co%s %s %.2e vs %.2e" % (x[0], x[1], x[2], x[3], x[4], x[5]) for x in bad2) or "none"))
    below = [x for x in verd["pp3"] if x[3] is not None and x[3] < 1]
    high = [x for x in verd["pp3"] if x[2] >= 0.5]
    highbad = [x for x in high if x[3] is None or x[3] >= 1]
    L.append("- PP3: E_shape below plicRDF's in %d of %d compared cases (>= 12 of 18 required); at Co 0.5 and 0.95 above plicRDF's in: %s" % (len(below), len(verd["pp3"]), "; ".join("%s N%d Co%s ratio %.3f" % (x[0], x[1], x[2], x[3] or float("nan")) for x in highbad) or "none"))
    d4 = []
    for fam in FAMS:
        for co in COS:
            a_, b_ = orders.get((fam, co, "pointIso")), orders.get((fam, co, "plicRDF"))
            d4.append((fam, co, a_, b_, (a_ - b_) if (a_ is not None and b_ is not None) else None))
    L.append("- PP4: order(pointIso) - order(plicRDF) per family and Co: " + "; ".join("%s Co%s: %s (%s vs %s)" % (x[0], x[1], "-" if x[4] is None else "%.2f" % x[4], "-" if x[2] is None else "%.2f" % x[2], "-" if x[3] is None else "%.2f" % x[3]) for x in d4) + "; falsified where the difference is below -0.1")
    w5 = [x for x in verd["pp5"] if x[3] is not None]
    L.append("- PP5: wall clock ratio to plicRDF at the same ranks: " + "; ".join("%s N%d Co%s %.2f" % (x[0], x[1], x[2], x[3]) for x in w5) + "; above 1.3 in %d of %d" % (len([x for x in w5 if x[3] > 1.3]), len(w5)))
    open(a.out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
if __name__ == "__main__":
    main()
