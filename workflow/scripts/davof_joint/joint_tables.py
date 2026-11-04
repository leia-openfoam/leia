#!/usr/bin/env python3
"""
The tables and figure data of the joint-function note from the prototype's
output (joint_prototype.py: <an>/<run>_s<step>_{summary,faces,donors,worst}.csv).

    joint_tables.py --an /work/.../an --out /work/.../tex

Writes LaTeX table bodies (booktabs rows only; the note supplies the
tabular environment) and pgfplots coordinate files:

  table_l1.tex        per run and step: the L1 over faces of the pieces, the
                      joint min, max, the state-guarded band variant against
                      the exact sphere, in units of the summed exact liquid,
                      and the ratio joint / pieces (the pre-registered 1.2)
  table_jump.tex      the jump volume (max - min) over the exact liquid, its
                      largest share of a swept volume, and the faces with k
                      planes
  table_k.tex         the distribution of k (planes touching) and k_caps
  table_af.tex        the face fractions: the mean L1 over faces of the pieces,
                      the donor plane, the joint min and the band variant
  table_bound.tex     the per-donor overdraw of the pieces, min, max, band
  table_cost.tex      the clip counts per face: the pieces now, after the
                      pre-test, the joint, the band variant
  hist_<run>_s<step>.dat  the histogram of (V - V_exact)/V_total for the
                      pieces and the joint min (pgfplots table)
  verdict.tex         the pre-registered rule evaluated per run and overall
"""
import argparse
import csv
import glob
import math
import os


def read_summary(path):
    with open(path) as fh:
        r = list(csv.reader(fh))
    return dict(zip(r[0], r[1]))


def fnum(x):
    v = float(x)
    if v == 0:
        return "0"
    if abs(v) >= 1e-2 and abs(v) < 1e3:
        return "%.3g" % v
    m, e = ("%.2e" % v).split("e")
    return "%s\\,\\cdot\\,10^{%d}" % (m, int(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--an", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    runs = []
    for p in sorted(glob.glob(os.path.join(args.an, "*_summary.csv"))):
        s = read_summary(p)
        name = os.path.basename(p)[:-len("_summary.csv")]
        run, step = name.rsplit("_s", 1)
        N = int(run.split("_")[0][1:]); co = run.split("_")[1]; row = run.split("_")[2]
        runs.append((N, {"co02": 0.2, "co1": 1.0}[co], row, int(step), s, name))
    runs.sort(key=lambda t: (t[0], t[1], t[2], t[3]))

    def label(N, co, row, step):
        return "%d & %.1f & %s & %d" % (N, co, row, step)

    verdict_rows = []
    with open(os.path.join(args.out, "table_l1.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            E = float(s["sumVexact"])
            l1p = float(s["L1_vPieces"])/E; l1m = float(s["L1_vMin"])/E; l1x = float(s["L1_vMax"])/E; l1b = float(s["L1_vBand"])/E
            l1c = float(s["L1_vCurv"])/E
            ratio = l1m/l1p if l1p > 0 else float("inf")
            ratiob = l1b/l1p if l1p > 0 else float("inf")
            fh.write("%s & %s & $%s$ & $%s$ & $%s$ & $%s$ & %.2f & %.2f \\\\\n" % (
                label(N, co, row, step), s["faces"], fnum(l1p), fnum(l1m), fnum(l1x), fnum(l1b), ratio, ratiob))
            over = float(s["overdrawMax_min"]); overb = float(s["overdrawMax_band"])
            share = float(s["kCapsLe3Share"])
            ok_l1 = ratio <= 1.2; ok_over = over <= 1e-12; ok_k = share >= 0.99
            verdict_rows.append((label(N, co, row, step), ratio, ok_l1, over, ok_over, share, ok_k, ratiob, overb))
    with open(os.path.join(args.out, "table_net.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            E = float(s["sumVexact"])
            fh.write("%s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\\n" % (
                label(N, co, row, step), fnum(float(s["net_vPieces"])/E), fnum(float(s["net_vMin"])/E),
                fnum(float(s["net_vMax"])/E), fnum(float(s["net_vBand"])/E), fnum(float(s["max_vMin"])/float(s["sumVtotal"])*float(s["faces"]))))
    with open(os.path.join(args.out, "table_jump.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            fh.write("%s & $%s$ & $%s$ & %s & %s & %s & %s \\\\\n" % (
                label(N, co, row, step), fnum(s["sumJumpRel"]), fnum(s["maxJumpRel"]),
                s["facesWithFull"], s["facesWithEmpty"], s["facesWithWisp"], s["fallbacks"]))
    with open(os.path.join(args.out, "table_k.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            n = int(s["faces"])
            ks = " & ".join(s["k%d" % k] for k in range(0, 6))
            fh.write("%s & %s & %s & %s & %.4f & %s / %s / %s \\\\\n" % (
                label(N, co, row, step), ks, s["kGt3"], s["kCapsGt3"], float(s["kCapsLe3Share"]),
                s["curv_convex"], s["curv_concave"], s["curv_saddle"]))
    with open(os.path.join(args.out, "table_af.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            fh.write("%s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\\n" % (
                label(N, co, row, step), fnum(s["L1_afPieces"]), fnum(s["L1_afDonor"]), fnum(s["L1_afMin"]),
                fnum(s["L1_afMax"]), fnum(s["L1_afBand"]), fnum(s["max_afMin"])))
    with open(os.path.join(args.out, "table_bound.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            # the overdraw of the joint min beyond the pieces' own, per donor:
            # max over donors of (out_min - out_pieces)/V, and the donors where
            # the joint min overdraws (> 1e-12 V) while the pieces do not
            drows = list(csv.DictReader(open(os.path.join(args.an, name + "_donors.csv"))))
            dmax = max((float(r["overMin"]) - float(r["overPieces"]) for r in drows), default=0.0)
            dnew = sum(1 for r in drows if float(r["overMin"]) > 1e-12 and float(r["overPieces"]) <= 1e-12)
            fh.write("%s & %s & $%s$ & $%s$ & $%s$ & $%s$ & %s / %s / %s & $%s$ & %d \\\\\n" % (
                label(N, co, row, step), s["donors"], fnum(s["overdrawMax_pieces"]), fnum(s["overdrawMax_min"]),
                fnum(s["overdrawMax_max"]), fnum(s["overdrawMax_band"]),
                s["overdrawCount_pieces"], s["overdrawCount_min"], s["overdrawCount_band"], fnum(dmax), dnew))
            s["overdrawBeyondPieces"] = dmax; s["overdrawNewDonors"] = dnew
    with open(os.path.join(args.out, "table_cost.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            n = int(s["faces"])
            fh.write("%s & %.1f & %.2f & %.0f & %.0f & %.0f & %.0f & %.1f \\\\\n" % (
                label(N, co, row, step), float(s["ringCellsPerFace"]), float(s["piecesPerFace"]),
                float(s["clipsPiecesNow"])/n, float(s["clipsPiecesPretest"])/n, float(s["clipsJoint"])/n,
                float(s["clipsBand"])/n, float(s["clipsPiecesPretest"])/max(1.0, float(s["clipsJoint"]))))
    with open(os.path.join(args.out, "table_check.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            fh.write("%s & %s & $%s$ & $%s$ & %.0f \\\\\n" % (
                label(N, co, row, step), s["piecesCheckFaces"], fnum(s["piecesCheckMaxRel"]), fnum(s["listCheckMaxRel"]), float(s["seconds"])))
    with open(os.path.join(args.out, "verdict.tex"), "w") as fh:
        allok = True; allok_rel = True
        for (lab, ratio, ok1, over, ok2, share, ok3, ratiob, overb), (N, co, row, step, s, name) in zip(verdict_rows, runs):
            overp = float(s["overdrawMax_pieces"]); beyond = float(s.get("overdrawBeyondPieces", 0.0))
            ok2rel = beyond <= 1e-12
            fh.write("%s & %.2f & %s & $%s$ & $%s$ & %s & $%s$ & %s & %.4f & %s \\\\\n" % (
                lab, ratio, "yes" if ok1 else "no", fnum(over), fnum(overp), "yes" if ok2 else "no",
                fnum(beyond), "yes" if ok2rel else "no", share, "yes" if ok3 else "no"))
            allok = allok and ok1 and ok2 and ok3
            allok_rel = allok_rel and ok1 and ok2rel and ok3
        fh.write("%% overall, the literal rule: %s; with the overdraw clause read against the pieces: %s\n" % (
            "POSITIVE" if allok else "NEGATIVE", "POSITIVE" if allok_rel else "NEGATIVE"))
    # the breakdown by the class of the face: the number of planes k (0: no
    # plane, identical by construction; 1: one plane, the extension into the
    # cells without a plane is the only departure; 2; >= 3) and whether a
    # wisp's matched plane is among them; per class the L1 of the pieces, the
    # joint min and the band variant against the truth, in units of the
    # class's exact liquid, and the ratios
    with open(os.path.join(args.out, "table_class.tex"), "w") as fh:
        for N, co, row, step, s, name in runs:
            if step == 1:
                continue
            rows = list(csv.DictReader(open(os.path.join(args.an, name + "_faces.csv"))))
            classes = [("$k = 0$", lambda r: int(r["k"]) == 0), ("$k = 1$", lambda r: int(r["k"]) == 1),
                       ("$k = 2$", lambda r: int(r["k"]) == 2), ("$k \\ge 3$", lambda r: int(r["k"]) >= 3),
                       ("with a wisp plane", lambda r: int(r["nWisp"]) == 0 and "wisp" in r.get("classes", "") ),
                       ("with a full cell", lambda r: int(r["nFull"]) > 0),
                       ("with an empty cell", lambda r: int(r["nEmpty"]) > 0)]
            first = True
            for lab, pred in classes:
                sel = [r for r in rows if pred(r) and float(r["vExact"]) >= 0]
                if not sel:
                    continue
                E = sum(float(r["vExact"]) for r in sel)
                if E <= 0:
                    continue
                l1p = sum(abs(float(r["vPieces"]) - float(r["vExact"])) for r in sel)/E
                l1m = sum(abs(float(r["vMin"]) - float(r["vExact"])) for r in sel)/E
                l1b = sum(abs(float(r["vBand"]) - float(r["vExact"])) for r in sel)/E
                netm = sum(float(r["vMin"]) - float(r["vExact"]) for r in sel)/E
                netp = sum(float(r["vPieces"]) - float(r["vExact"]) for r in sel)/E
                fh.write("%s & %s & %d & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & %.2f & %.2f \\\\\n" % (
                    label(N, co, row, step) if first else " & & & ", lab, len(sel), fnum(l1p), fnum(l1m), fnum(l1b),
                    fnum(netp), fnum(netm), l1m/l1p if l1p > 0 else float("inf"), l1b/l1p if l1p > 0 else float("inf")))
                first = False
            fh.write("\\midrule\n")
    # the histograms of the signed error over the swept volume
    for N, co, row, step, s, name in runs:
        rows = list(csv.DictReader(open(os.path.join(args.an, name + "_faces.csv"))))
        edges = [-0.2 + 0.01*i for i in range(41)]
        def hist(key):
            cnt = [0]*40
            for r in rows:
                vE = float(r["vExact"]); vT = float(r["vTotal"])
                if vE < 0 or vT <= 0:
                    continue
                x = (float(r[key]) - vE)/vT
                i = int(math.floor((x + 0.2)/0.01))
                i = max(0, min(39, i))
                cnt[i] += 1
            return cnt
        hp = hist("vPieces"); hm = hist("vMin"); hb = hist("vBand")
        with open(os.path.join(args.out, "hist_%s.dat" % name), "w") as fh:
            fh.write("x pieces min band\n")
            for i in range(40):
                fh.write("%.3f %d %d %d\n" % (0.5*(edges[i] + edges[i + 1]), hp[i], hm[i], hb[i]))
    # the worst faces of the joint min on the N = 64, Co 1, ring, mid-turn dump (or the last run)
    pick = [t for t in runs if t[0] == 64 and t[1] == 1.0 and t[2] == "ring" and t[3] > 1] or runs[-1:]
    if pick:
        name = pick[0][5]
        rows = list(csv.DictReader(open(os.path.join(args.an, name + "_worst.csv"))))[:8]
        with open(os.path.join(args.out, "table_worst.tex"), "w") as fh:
            fh.write("%% %s\n" % name)
            for r in rows:
                vT = float(r["vTotal"]); vE = float(r["vExact"])
                fh.write("%s & %s & $%s$ & $%s$ & $%s$ & $%s$ & %s & %s & %s & %s & %s \\\\\n" % (
                    r["facei"], r["donor"], fnum((float(r["vPieces"]) - vE)/vT), fnum((float(r["vMin"]) - vE)/vT),
                    fnum((float(r["vMax"]) - vE)/vT), fnum((float(r["vBand"]) - vE)/vT),
                    r["k"], r["kCaps"], r["nFull"], r["nEmpty"], r["curv"]))
    print("wrote tables for %d run-steps to %s" % (len(runs), args.out))


if __name__ == "__main__":
    main()
