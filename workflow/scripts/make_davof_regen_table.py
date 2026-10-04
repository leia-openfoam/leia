#!/usr/bin/env python3
"""DAVOF zero-step regeneration (davofRegeneration.H): order tables and figures.

Reads, for every case of a leiaTestDavofNormal study run with a regeneration
arm (config/davof/*Regeneration3D.yaml),

    <case>/case_params.json          tokens.N_CELLS, DOMAIN_LENGTH,
                                     DAVOF_ALPHA_SOURCE, DAVOF_REGEN_ARM
    <case>/leiaTestDavofRegen.csv    one row: the face-fraction errors, the
                                     owner/neighbour mismatch, the regenerated
                                     normal and position, the counts
    <case>/leiaTestDavofNormalModels.csv  the davof row: the ORIGINAL state's
                                     normal and position (for comparison)

and writes into the theme's data source (paths.tables_dir / figs_dir):

    tables/davof_regen_<study>.csv            the long table, every (source, arm, N)
    tables/davof_regen_orders_<study>.csv     least-squares, three-finest and pairwise
                                              orders per (source, arm) and metric
    tables/davof_regen_orders_<study>.tex     booktabs body of the same
    tables/davof_regen_proposal_<study>.tex   the compact table (one row per arm:
                                              the orders of the face fractions, the
                                              normal and the position, three finest
                                              rungs, all rungs in parentheses)
    figures/davof_regen_<metric>_<study>.png  log-log error against h per source

Series = (alphaSource, arm). h = L/N. Orders as in make_davof_normal_table.py:
least squares over the rungs with e > 0, the same over the three finest, and
the pairwise orders.

Usage (from the repo root):
    python3 workflow/scripts/make_davof_regen_table.py studies/davof/ellipsoidRegeneration3D [--theme davof]
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

CSV_REGEN = "leiaTestDavofRegen.csv"
CSV_MODELS = "leiaTestDavofNormalModels.csv"
METRICS = ["E_ALPHAF_L1", "E_ALPHAF_L2", "E_ALPHAF_LINF",
           "E_ALPHAF_MISMATCH_L1", "E_ALPHAF_MISMATCH_LINF",
           "E_L1_N", "E_L2_N", "E_LINF_N", "E_POS_L1", "E_POS_L2", "E_POS_LINF"]
ORIG = {"E_L1_N": "E0_L1_N", "E_L2_N": "E0_L2_N", "E_POS_L2": "E0_POS_L2"}
COUNTS = ["N_INTERFACE", "N_INTERFACE_REGEN", "N_REGEN_CELLS", "N_NO_MODEL",
          "N_NO_BRACKET", "N_TRI_FALLBACK", "N_REF_FALLBACK", "N_PAIRS",
          "N_MISMATCH_FACES", "ANCHOR_RESIDUAL_MAX", "MAX_VOL_DIFF_PLANE_REGEN",
          "CPU_SECONDS"]
FIG_METRICS = ["E_ALPHAF_L1", "E_ALPHAF_L2", "E_ALPHAF_MISMATCH_L1", "E_L1_N", "E_L2_N", "E_POS_L2"]
ARMS = ["plane", "paraboloidPlane", "paraboloidFaces", "paraboloidVolume", "paraboloidVolumeExact"]
ARM_LABEL = {"plane": "plane (2)",
             "paraboloidPlane": "paraboloid, vertex on the plane",
             "paraboloidFaces": "paraboloid anchored by the face fractions",
             "paraboloidVolume": "paraboloid, volume on the DA tets",
             "paraboloidVolumeExact": "paraboloid, own cut volume"}
COLORS = {"plane": "#2a78d6", "paraboloidPlane": "#eb6834", "paraboloidFaces": "#1baf7a",
          "paraboloidVolume": "#eda100", "paraboloidVolumeExact": "#7a3bb5"}
MARKERS = {"plane": "o", "paraboloidPlane": "s", "paraboloidFaces": "D",
           "paraboloidVolume": "^", "paraboloidVolumeExact": "v"}
SOURCE_LABEL = {"exactSphere": "exact sphere", "quadraticFaces": "DA + quadratic faces",
                "detrixheAslam": "Detrixhe-Aslam", "planePhaseIndicator": "plane indicator"}
METRIC_LABEL = {
    "E_ALPHAF_L1": r"$L_1$ of $|\alpha_f^{\mathrm{regen}} - \alpha_f^{\mathrm{ref}}|$",
    "E_ALPHAF_L2": r"$L_2$ of $|\alpha_f^{\mathrm{regen}} - \alpha_f^{\mathrm{ref}}|$",
    "E_ALPHAF_MISMATCH_L1": r"$L_1$ of the owner/neighbour mismatch of $\alpha_f$",
    "E_L1_N": r"$L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$ (regenerated)",
    "E_L2_N": r"$L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$ (regenerated)",
    "E_POS_L2": r"$L_2$ of the plane-centroid distance to the surface [m] (regenerated)",
}


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _fit(h, e):
    m = [(a, b) for a, b in zip(h, e) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, len(m)
    lh, le = np.log([a for a, _ in m]), np.log([b for _, b in m])
    return float(np.polyfit(lh, le, 1)[0]), len(m)


def _pairwise(h, e):
    out = []
    for k in range(len(h) - 1):
        a, b = e[k], e[k + 1]
        if a and b and a > 0 and b > 0 and h[k] > 0 and h[k + 1] > 0 and h[k] != h[k + 1]:
            out.append(np.log(a/b)/np.log(h[k]/h[k + 1]))
        else:
            out.append(None)
    return out


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def collect(study_dir):
    series = {}
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        cdir = os.path.dirname(meta)
        rpath = os.path.join(cdir, CSV_REGEN)
        if not os.path.isfile(rpath) or os.path.getsize(rpath) == 0:
            continue
        with open(meta) as fh:
            tok = json.load(fh).get("tokens", {})
        n, L = _f(tok.get("N_CELLS")), _f(tok.get("DOMAIN_LENGTH"))
        if not n or not L:
            continue
        row = _rows(rpath)[-1]
        src = row.get("ALPHA_SOURCE") or tok.get("DAVOF_ALPHA_SOURCE")
        arm = row.get("REGEN_ARM") or tok.get("DAVOF_REGEN_ARM")
        rec = {"N": int(n), "h": L/n, "R_OVER_H": _f(row.get("R_OVER_H")), "case": os.path.basename(cdir)}
        for m in METRICS:
            rec[m] = _f(row.get(m))
        for c in COUNTS:
            rec[c] = row.get(c, "")
        for m in ORIG.values():
            rec[m] = None
        mpath = os.path.join(cdir, CSV_MODELS)
        if os.path.isfile(mpath) and os.path.getsize(mpath) > 0:
            for mrow in _rows(mpath):
                if mrow.get("MODEL") == "davof":
                    for m, m0 in ORIG.items():
                        rec[m0] = _f(mrow.get(m))
        series.setdefault((src, arm), []).append(rec)
    for k in series:
        series[k].sort(key=lambda r: r["N"])
    return series


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="davof")
    a = ap.parse_args(argv)
    study = os.path.basename(os.path.normpath(a.study_dir))
    series = collect(a.study_dir)
    if not series:
        print(f"[regen] {study}: no {CSV_REGEN} found; nothing to do")
        return 0
    tables, figs = paths.tables_dir(a.theme), paths.figs_dir(a.theme)

    # --- the long table ------------------------------------------------------
    lpath = os.path.join(tables, f"davof_regen_{study}.csv")
    with open(lpath, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ALPHA_SOURCE", "REGEN_ARM", "N", "h", "R_OVER_H"] + METRICS + list(ORIG.values()) + COUNTS + ["case"])
        for (src, arm), recs in sorted(series.items(), key=lambda kv: (kv[0][0], ARMS.index(kv[0][1]) if kv[0][1] in ARMS else 9)):
            for r in recs:
                w.writerow([src, arm, r["N"], r["h"], r["R_OVER_H"]]
                           + [r[m] for m in METRICS] + [r[m] for m in ORIG.values()]
                           + [r[c] for c in COUNTS] + [r["case"]])
    print(f"[regen] wrote {lpath}")

    # --- the orders ------------------------------------------------------------
    orders = {}
    opath = os.path.join(tables, f"davof_regen_orders_{study}.csv")
    with open(opath, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ALPHA_SOURCE", "REGEN_ARM", "METRIC", "P_ALL", "N_ALL", "P_3FINEST",
                    "PAIRWISE", "E_FINEST", "N_FINEST"])
        for (src, arm), recs in series.items():
            h = [r["h"] for r in recs]
            for m in METRICS + list(ORIG.values()):
                e = [r[m] for r in recs]
                p, n = _fit(h, e)
                p3, _ = _fit(h[-3:], e[-3:]) if len(h) >= 3 else (None, 0)
                pw = _pairwise(h, e)
                orders[(src, arm, m)] = (p, p3, pw, e[-1] if e else None)
                w.writerow([src, arm, m,
                            "" if p is None else f"{p:.3f}", n,
                            "" if p3 is None else f"{p3:.3f}",
                            " ".join("" if q is None else f"{q:.2f}" for q in pw),
                            "" if e[-1] is None else f"{e[-1]:.6e}", recs[-1]["N"]])
    print(f"[regen] wrote {opath}")

    def _p(v):
        if v is None:
            return "--"
        s = f"{v:.2f}"
        return "0.00" if s == "-0.00" else s

    def _e(v):
        return "--" if v is None else f"{v:.2e}".replace("e-0", "e-").replace("e+0", "e+")

    # --- the full booktabs table -------------------------------------------
    tpath = os.path.join(tables, f"davof_regen_orders_{study}.tex")
    with open(tpath, "w") as fh:
        fh.write("% generated by make_davof_regen_table.py from " + study + "\n")
        fh.write("\\begin{tabular}{llcccccccc}\n\\toprule\n")
        fh.write("state & arm & $p(\\alpha_f, L_1)$ & $p(\\alpha_f, L_2)$ & $p(\\alpha_f, L_\\infty)$ & "
                 "$p(\\text{mismatch})$ & $p(\\mathbf{n}, L_1)$ & $p(\\mathbf{n}, L_2)$ & $p(\\mathbf{x})$ & "
                 "$E_{\\alpha_f, L_1}$ finest \\\\\n\\midrule\n")
        for (src, arm), recs in sorted(series.items(), key=lambda kv: (kv[0][0], ARMS.index(kv[0][1]) if kv[0][1] in ARMS else 9)):
            cells = []
            for m in ("E_ALPHAF_L1", "E_ALPHAF_L2", "E_ALPHAF_LINF", "E_ALPHAF_MISMATCH_L1", "E_L1_N", "E_L2_N", "E_POS_L2"):
                p, p3, _, _ = orders[(src, arm, m)]
                cells.append(f"{_p(p3)} ({_p(p)})")
            ef = orders[(src, arm, "E_ALPHAF_L1")][3]
            fh.write(f"{SOURCE_LABEL.get(src, src)} & {ARM_LABEL.get(arm, arm)} & " + " & ".join(cells) + f" & {_e(ef)} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[regen] wrote {tpath}")

    # --- the compact table for the proposal ---------------------------------
    ppath = os.path.join(tables, f"davof_regen_proposal_{study}.tex")
    with open(ppath, "w") as fh:
        fh.write("% generated by make_davof_regen_table.py from " + study + "; orders over the three finest rungs (all rungs)\n")
        fh.write("\\begin{tabular}{llcccccc}\n\\toprule\n")
        fh.write("state & clipping surface & $p_1(\\alpha_f)$ & $p_2(\\alpha_f)$ & $p_1(\\mathbf{n})$ & $p_2(\\mathbf{n})$ & "
                 "$p(\\mathbf{x})$ & $L_1(\\alpha_f)$ finest \\\\\n\\midrule\n")
        for (src, arm), recs in sorted(series.items(), key=lambda kv: (kv[0][0], ARMS.index(kv[0][1]) if kv[0][1] in ARMS else 9)):
            cells = []
            for m in ("E_ALPHAF_L1", "E_ALPHAF_L2", "E_L1_N", "E_L2_N", "E_POS_L2"):
                p, p3, _, _ = orders[(src, arm, m)]
                cells.append(f"{_p(p3)} ({_p(p)})")
            ef = orders[(src, arm, "E_ALPHAF_L1")][3]
            fh.write(f"{SOURCE_LABEL.get(src, src)} & {ARM_LABEL.get(arm, arm)} & " + " & ".join(cells) + f" & {_e(ef)} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[regen] wrote {ppath}")

    # --- figures -------------------------------------------------------------
    sources = sorted({src for src, _ in series})
    for metric in FIG_METRICS:
        fig, axes = plt.subplots(1, len(sources), figsize=(5.2*len(sources), 4.2), squeeze=False)
        for ax, src in zip(axes[0], sources):
            hmin, hmax, eref = None, None, None
            for arm in ARMS:
                recs = series.get((src, arm))
                if not recs:
                    continue
                h = [r["h"] for r in recs if r[metric]]
                e = [r[metric] for r in recs if r[metric]]
                if not h:
                    continue
                ax.loglog(h, e, marker=MARKERS[arm], color=COLORS[arm], lw=1.6, ms=6,
                          label=ARM_LABEL[arm])
                hmin = min(h) if hmin is None else min(hmin, min(h))
                hmax = max(h) if hmax is None else max(hmax, max(h))
                # The reference slopes hang from the largest coarsest-rung error.
                ecoarse = e[h.index(max(h))]
                eref = ecoarse if eref is None else max(eref, ecoarse)
            if hmin and eref:
                for p, ls in ((1, ":"), (2, "--")):
                    ax.loglog([hmax, hmin], [eref, eref*(hmin/hmax)**p], ls, color="0.5", lw=1,
                              label=f"$h^{p}$")
            ax.set_title(SOURCE_LABEL.get(src, src))
            ax.set_xlabel("h [m]")
            ax.set_ylabel(METRIC_LABEL.get(metric, metric))
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fpath = os.path.join(figs, f"davof_regen_{metric}_{study}.png")
        fig.savefig(fpath, dpi=150)
        plt.close(fig)
        print(f"[regen] wrote {fpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
