#!/usr/bin/env python3
"""diag2 task C tables: the pieces rows of the Stage B matrix recomputed (studies diag2_c_<fam>_N<N>), the plicRDF and
isoAlpha rows of Stage B reused (stored standing plicRDF at Co 0.2 / 0.5, ccsB_*_plic for the third target, ccsB_*_iso).
Columns as stageB_tables.md with the extremes of alpha before AND after the conservative bounding, and the ratios of the
pieces row to plicRDF in E_shape and in wall clock.

    python3 stageC_an.py --studies <leia>/studies/davof [--stageB <ccs leia>/studies/davof] [--stored <kiss standing dir>]
"""
import argparse, glob, os, sys, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(os.path.dirname(here), "davof_kiss"))
import ccs_stageB_an as sb
import kiss_an as ka
spec = importlib.util.spec_from_file_location("an1", here + "/an1.py")
an1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(an1)

FAMS = sb.FAMS


def load(dirs, rowmatch):
    out = []
    for d in dirs:
        for c in sorted(glob.glob(os.path.join(d, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = ka.read_case(c)
            if not r or not rowmatch(r):
                continue
            sb.extra(c, r)
            r["dir"] = c
            lg = os.path.join(c, "log.davofAdvectorFoam")
            r["before"] = None
            if os.path.isfile(lg):
                recs, abort = an1.parse(lg)
                r["abortStep"] = abort[0] if abort else None
                bmin = [x["bmin"] for x in recs if x["bmin"] is not None]
                bmax = [x["bmax"] for x in recs if x["bmax"] is not None]
                amin = [x["amin"] for x in recs if x["amin"] is not None]
                amax = [x["amax"] for x in recs if x["amax"] is not None]
                if bmin:
                    r["before"] = (min(bmin), max(bmax))
                if amin:
                    r["after"] = (min(amin), max(amax))
            out.append(r)
    return out


def f(v, nd=2):
    return "-" if v is None else "%.*e" % (nd, v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", required=True)
    ap.add_argument("--stageB", default="/work/scratch/tm83tomy/davof/agents/ccs/leia/studies/davof")
    ap.add_argument("--stored", default="/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    md = []
    rowsAll = []
    for fam, Ns in FAMS.items():
        recs = []
        for N in Ns:
            recs += [dict(r, src="diag2_c_%s_N%d" % (fam, N)) for r in load([os.path.join(a.studies, "diag2_c_%s_N%d" % (fam, N))], lambda r: True)]
            for kind in ("plic", "iso"):
                recs += [dict(r, src="ccsB_%s_N%d_%s" % (fam, N, kind)) for r in load([os.path.join(a.stageB, "ccsB_%s_N%d_%s" % (fam, N, kind))], lambda r: True)]
            recs += [dict(r, src="standing_%s_N%d (stored)" % (fam, N)) for r in
                     load([os.path.join(a.stored, "standing_%s_N%d" % (fam, N))], lambda r: r["row"] == "plicRDF" and r["rep"] == 1 and r["co"] in (0.2, 0.5))]
        for r in recs:
            r["method"] = sb.method_of(r); r["family"] = fam
        for coi, colabel in ((0, "0.2"), (1, "0.5"), (2, "0.95")):
            def sel(r):
                if coi == 0: return abs(r["co"] - 0.2) < 1e-9
                if coi == 1: return abs(r["co"] - 0.5) < 1e-9
                return 0.8 < r["co"] < 0.99
            rows = [r for r in recs if sel(r)]
            md.append("\n### %s, OpenFOAM Co %s\n" % (fam, colabel))
            md.append("| method | N | ranks | OF Co | steps | status | E_shape | E_shape / plicRDF | E_mass (%) | min alpha before bounding | min alpha after bounding | max alpha - 1 before bounding | max alpha - 1 after bounding | L1 abs(n-n_e) | L1 (1-n.n_e) | Linf (1-n.n_e) | wall clock (s) | wall / plicRDF | source |")
            md.append("|" + "---|"*19)
            for mth in ("pieces", "plicRDF", "isoAlpha"):
                for N in Ns:
                    rs = [r for r in rows if r["method"] == mth and r["N"] == N]
                    if not rs:
                        md.append("| %s | %d | | | | not run |" % (mth, N) + " |"*13); continue
                    r = rs[0]
                    pl = [x for x in rows if x["method"] == "plicRDF" and x["N"] == N]
                    ok = r["status"] in ("ended", "completed") and r["Eshape"] is not None
                    plE = pl[0]["Eshape"] if pl and pl[0]["Eshape"] is not None else None
                    plW = pl[0]["wall"] if pl and pl[0]["wall"] else None
                    rE = ("%.2f" % (r["Eshape"]/plE)) if (ok and plE) else "-"
                    rW = ("%.1f" % (r["wall"]/plW)) if (r["wall"] and plW) else "-"
                    b = r.get("before"); af = r.get("after")
                    if mth == "pieces":
                        bmin = f(b[0]) if b else "-"; bmax = f(b[1]) if b else "-"
                        amin = f(af[0]) if af else f(r["minAlpha"]); amax = f(af[1]) if af else f(r["excess"])
                    else:
                        bmin = bmax = "n/a"; amin = f(r["minAlpha"]); amax = f(r["excess"])
                    st = r["status"] if ok else ("%s: %s" % (r["status"], (r.get("reason") or "")[:60]))
                    md.append("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                        mth, N, r["np"], ("%.3f" % r["ofCo"]) if r["ofCo"] else "-", r["steps"], st,
                        f(r["Eshape"]) if ok else "div.", rE, f(r["massDevPct"]) if ok else "div.", bmin, amin, bmax, amax,
                        f(r["L1n"]) if ok else "div.", f(r["L1n2"]) if ok else "div.", f(r["Linfn2"]) if ok else "div.",
                        ("%.0f" % r["wall"]) if r["wall"] is not None else "-", rW, r["src"]))
                    rowsAll.append((fam, colabel, mth, N, r["status"], r["Eshape"], plE, b))
            md.append("\nOrders over N (least squares; pairwise in brackets), non-completed rungs excluded:\n")
            md.append("| method | E_shape | L1 abs(n-n_e) | L1 (1-n.n_e) |")
            md.append("|---|---|---|---|")
            for mth in ("pieces", "plicRDF", "isoAlpha"):
                cells = []
                for key in ("Eshape", "L1n", "L1n2"):
                    pts = [(N, [r for r in rows if r["method"] == mth and r["N"] == N][0][key]) for N in Ns
                           if [r for r in rows if r["method"] == mth and r["N"] == N and r["status"] in ("ended", "completed")]]
                    ls, pw = sb.fit(pts)
                    cells.append("-" if ls is None else "%.2f (%s)" % (ls, ", ".join("%.2f" % p for p in pw)))
                md.append("| %s | %s |" % (mth, " | ".join(cells)))
    # verdicts
    pieces = [x for x in rowsAll if x[2] == "pieces"]
    done = [x for x in pieces if x[4] in ("ended", "completed")]
    md.append("\n### Verdicts\n")
    md.append("PC1: %d of %d pieces cases complete." % (len(done), 27))
    exc = [x for x in pieces if x[7] and (x[7][0] < -1e-3 or x[7][1] > 1e-3)]
    md.append("PC2: %d cases with the extremes of alpha before the bounding outside [-1e-3, 1 + 1e-3]: %s" % (len(exc), "; ".join("%s N%d Co %s: min %.3e max-1 %.3e" % (x[0], x[3], x[1], x[7][0], x[7][1]) for x in exc)))
    over = [x for x in done if x[6] and x[5] > 1.5*x[6]]
    md.append("PC3: %d cases with E_shape above 1.5 x plicRDF: %s" % (len(over), "; ".join("%s N%d Co %s: %.2f" % (x[0], x[3], x[1], x[5]/x[6]) for x in over)))
    txt = "\n".join(md) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
