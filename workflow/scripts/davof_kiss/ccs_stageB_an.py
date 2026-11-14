#!/usr/bin/env python3
"""Stage B tables of the CCS alignment (2026-11-14; docs/davof/REPORT-swept.md "CCS aligned with the
published method"). Reads the cases of studies/davof/ccsB_<family>_N<N>_{davof,plic,iso} and the stored
plicRDF cases (Co 0.2, 0.5) of the standing studies of the KISS workspace, writes one table per family
and Co (rows: the four methods, one block per N) with the orders over N, and a CSV.

    python3 workflow/scripts/davof_kiss/ccs_stageB_an.py [--out DIR] [--stored DIR]
"""
import argparse, csv, glob, json, math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kiss_an as ka

REPO = ka.REPO
FAMS = {"shear2D": (64, 128, 256), "shear3D": (32, 64, 128), "deformation3D": (32, 64, 128)}
COS = {"shear2D": {64: 0.95175, 128: 0.95096, 256: 0.95077}}  # unused; the third target is read from the case
CUTS_RE = re.compile(r"positioning ccs: (\d+) positionings, cuts per cell (\S+) \(max (\d+)")
METHODS = ("onePlane", "pieces", "plicRDF", "isoAlpha")


def num(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


def extra(c, rec):
    """The columns that kiss_an.read_case does not hold: Emass, L2 and Linf of the normal, interface cells, cuts, abort reason."""
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) and os.path.getsize(csvp) else []
    last = rows[-1] if rows else {}
    rec["Emass"] = num(last.get("Emass"))
    rec["L2n"] = num(last.get("LNormalDiff2"))
    rec["Linfn"] = num(last.get("LNormalDiffInf"))
    rec["nIf"] = num(last.get("nInterface"))
    rec["L1n2"] = (rec["L2n"]**2/2) if rec["L2n"] is not None else None        # 1 - n.n_e = |n - n_e|^2 / 2 on unit vectors
    rec["Linfn2"] = (rec["Linfn"]**2/2) if rec["Linfn"] is not None else None
    logp = next((os.path.join(c, n) for n in ("log.davofAdvectorFoam", "log.advectorVoF") if os.path.isfile(os.path.join(c, n))), None)
    rec["cuts"] = None; rec["cutsMax"] = None; rec["reason"] = ""
    if logp:
        log = open(logp, errors="replace").read()
        m = CUTS_RE.findall(log)
        if m:
            tot = sum(int(a) for a, b, c_ in m)
            rec["cuts"] = sum(int(a)*float(b) for a, b, c_ in m)/max(tot, 1)
            rec["cutsMax"] = max(int(c_) for a, b, c_ in m)
        mm = re.search(r"ABORT at step (\d+)[^\n]*", log)
        sg = re.search(r"#1\s+Foam::(sig\w+)::sigHandler", log)
        fr = re.search(r"#3\s+(?:void |bool |Foam::\w+ )?(Foam::[A-Za-z_:]+)", log)
        if mm:
            rec["reason"] = mm.group(0)[:150].replace("|", "/")
        elif sg:
            rec["reason"] = ("%s in %s" % (sg.group(1), fr.group(1) if fr else "?")).replace("|", "/")
            st_ = re.search(r"\nTime = (\S+)[^\n]*\n(?:(?!\nTime = ).)*\[stack trace\]", log, re.S)
            rec["reason"] += " (last Time = %s)" % (st_.group(1) if st_ else "?")
        elif rec["status"] not in ("ended", "completed"):
            tail = [l for l in log.strip().split("\n") if l.strip()][-3:]
            rec["reason"] = " / ".join(t.strip()[:80] for t in tail).replace("|", "/").replace("#", "")
    return rec


def load(dirs, rowmatch):
    out = []
    for d in dirs:
        for c in sorted(glob.glob(os.path.join(d, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = ka.read_case(c)
            if not r or not rowmatch(r):
                continue
            extra(c, r)
            out.append(r)
    return out


def method_of(r):
    if r["row"] == "plicRDF": return "plicRDF"
    if r["row"] == "isoAlpha": return "isoAlpha"
    return "onePlane" if r["setting"] == "onePlane" else "pieces"


def fmt(v, nd=2):
    return "-" if v is None else "%.*e" % (nd, v)


def fit(pts):
    """least squares slope of log E over log h (h = 1/N) and the pairwise slopes; pts = [(N, E)]"""
    pts = [(n, e) for n, e in pts if e is not None and e > 0]
    pts.sort()
    if len(pts) < 2:
        return None, []
    xs = [math.log(1.0/n) for n, e in pts]; ys = [math.log(e) for n, e in pts]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    ls = sum((x - mx)*(y - my) for x, y in zip(xs, ys))/sum((x - mx)**2 for x in xs)
    pw = [(math.log(pts[i+1][1]) - math.log(pts[i][1]))/(math.log(1.0/pts[i+1][0]) - math.log(1.0/pts[i][0])) for i in range(len(pts)-1)]
    return ls, pw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies-dir", default=os.path.join(REPO, "studies", "davof"))
    ap.add_argument("--stored", default="/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof")
    ap.add_argument("--out", default=os.path.join(REPO, "studies", "davof", "ccs_tables"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    allrows = []
    md = []
    for fam, Ns in FAMS.items():
        recs = []
        for N in Ns:
            for kind in ("davof", "plic", "iso"):
                recs += [dict(r, src="ccsB_%s_N%d_%s" % (fam, N, kind)) for r in load([os.path.join(a.studies_dir, "ccsB_%s_N%d_%s" % (fam, N, kind))], lambda r: True)]
            # the stored plicRDF cases at Co 0.2 and 0.5, repetition 1 (the standing studies)
            recs += [dict(r, src="standing_%s_N%d (stored)" % (fam, N)) for r in
                     load([os.path.join(a.stored, "standing_%s_N%d" % (fam, N))], lambda r: r["row"] == "plicRDF" and r["rep"] == 1 and r["co"] in (0.2, 0.5))]
        for r in recs:
            r["method"] = method_of(r)
            r["family"] = fam
        # the Co index: 0.2, 0.5, then the third target
        for coi, colabel in ((0, "0.2"), (1, "0.5"), (2, "0.95")):
            def sel(r):
                if coi == 0: return abs(r["co"] - 0.2) < 1e-9
                if coi == 1: return abs(r["co"] - 0.5) < 1e-9
                return r["co"] > 0.8 and r["co"] < 0.99
            rows = [r for r in recs if sel(r)]
            md.append("\n### %s, OpenFOAM Co %s (target %s)\n" % (fam, colabel, "/".join(sorted({str(r['co']) for r in rows}))))
            md.append("| method | N | ranks | OF Co | steps | status | E_shape | E_mass (%) | min alpha | max alpha - 1 | L1 abs(n-n_e) | Linf abs(n-n_e) | L1 (1-n.n_e) | Linf (1-n.n_e) | interface cells | limiter cells | limiter volume (cells) | CCS cuts/positioning (max) | wall clock (s) | source |")
            md.append("|" + "---|"*20)
            for mth in METHODS:
                for N in Ns:
                    rs = [r for r in rows if r["method"] == mth and r["N"] == N]
                    if not rs:
                        md.append("| %s | %d | | | | not run | | | | | | | | | | | | | | |" % (mth, N)); continue
                    r = rs[0]
                    ok = r["status"] in ("ended", "completed") and r["Eshape"] is not None
                    div = "div." if not ok else None
                    h = 1.0/N
                    lim = r["rule2"]
                    limc = int(lim["limCellsSum"]) if lim else ("n/a" if mth == "pieces" else "-")
                    limv = ("%.2e" % (lim["limVolumeSum"]/h**3)) if lim else ("n/a" if mth == "pieces" else "-")
                    cuts = ("%.2f (%d)" % (r["cuts"], r["cutsMax"])) if r["cuts"] is not None else "-"
                    st = r["status"] if ok else ("%s: %s" % (r["status"], r["reason"][:60]))
                    md.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                        mth, N, r["np"], ("%.3f" % r["ofCo"]) if r["ofCo"] else "-", r["steps"], st,
                        div or fmt(r["Eshape"]), div or fmt(r["massDevPct"]), div or fmt(r["minAlpha"]), div or fmt(r["excess"]),
                        div or fmt(r["L1n"]), div or fmt(r["Linfn"]), div or fmt(r["L1n2"]), div or fmt(r["Linfn2"]),
                        div or ("%d" % r["nIf"] if r["nIf"] is not None else "-"), limc, limv, cuts,
                        ("%.0f" % r["wall"]) if r["wall"] is not None else "-", r["src"]))
                    allrows.append(dict(family=fam, co=colabel, method=mth, N=N, np=r["np"], ofCo=r["ofCo"], steps=r["steps"], status=r["status"],
                                        reason=r["reason"], Eshape=r["Eshape"], Emass_pct=r["massDevPct"], minAlpha=r["minAlpha"], excess=r["excess"],
                                        L1n=r["L1n"], Linfn=r["Linfn"], L1n2=r["L1n2"], Linfn2=r["Linfn2"], nInterface=r["nIf"],
                                        limCells=lim.get("limCellsSum") if lim else None, limVolumeCells=(lim["limVolumeSum"]/h**3) if lim else None,
                                        cuts=r["cuts"], cutsMax=r["cutsMax"], wall=r["wall"], src=r["src"]))
            # the orders over N (least squares and pairwise), div. rungs excluded
            md.append("\nOrders over N (least squares; pairwise in brackets), diverged rungs excluded:\n")
            md.append("| method | E_shape | L1 abs(n-n_e) | L1 (1-n.n_e) |")
            md.append("|---|---|---|---|")
            for mth in METHODS:
                cells = []
                for key in ("Eshape", "L1n", "L1n2"):
                    pts = []
                    for N in Ns:
                        rs = [r for r in rows if r["method"] == mth and r["N"] == N]
                        if rs and rs[0]["status"] in ("ended", "completed"):
                            pts.append((N, rs[0][key]))
                    ls, pw = fit(pts)
                    cells.append("-" if ls is None else "%.2f (%s)" % (ls, ", ".join("%.2f" % p for p in pw)))
                md.append("| %s | %s |" % (mth, " | ".join(cells)))
    open(os.path.join(a.out, "stageB_tables.md"), "w").write("\n".join(md) + "\n")
    if allrows:
        with open(os.path.join(a.out, "stageB_cells.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)
    print("\n".join(md))


if __name__ == "__main__":
    main()
