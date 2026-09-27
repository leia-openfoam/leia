#!/usr/bin/env python3
"""Result figures and tables of the pre-print, from the archive of ONE software version.

  make_result_figures.py --archive data/archive/<version> --figures data/figures --tables data/tables

Reads only the archive (make_archive.py writes it) and the gate and candidate definitions
(config/gates/methodGate2D.yaml, config/candidates/*.yaml). Writes PDF figures and LaTeX tables;
the gate tables themselves come from workflow/scripts/make_gate_tables.py on
<archive>/gate/summaries.
"""
import argparse
import csv
import math
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
CANDS = ["baseline", "S1", "HL0", "HL1q", "HL1z", "HL2", "FP0"]
COLOR = {"baseline": "black", "S1": "tab:orange", "HL0": "tab:blue", "HL1q": "tab:green",
         "HL1z": "tab:red", "HL2": "tab:purple", "FP0": "tab:brown"}
MARK = {"baseline": "o", "S1": "s", "HL0": "^", "HL1q": "v", "HL1z": "D", "HL2": "P", "FP0": "X"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "legend.fontsize": 7,
                     "lines.linewidth": 1.0, "lines.markersize": 3.5, "figure.dpi": 150})


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def fmt(v, digits=2):
    if v is None:
        return "--"
    if v == 0:
        return "0"
    return r"\num{" + f"{v:.{digits}e}" + "}"


def save(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"[figures] {os.path.relpath(path)}")


def write_tex(path, text):
    with open(path, "w") as f:
        f.write(text)
    print(f"[tables] {os.path.relpath(path)}")


# --------------------------------------------------------------------------- figures
def fig_shear(A, out):
    S = {c: [r for r in rows(os.path.join(A, "gate", "summaries", c, "summary.csv")) if r["arm"] == "shear"]
         for c in CANDS if os.path.isfile(os.path.join(A, "gate", "summaries", c, "summary.csv"))}
    H = rows(os.path.join(A, "gate", "histories_shear.csv"))
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for c, R in S.items():
        h = [num(r["h"]) for r in R]; e = [num(r["shapeError"]) for r in R]
        ax[0].loglog(h, e, marker=MARK[c], color=COLOR[c], label=c)
    h0 = 1 / 136.0
    for p, lab, y0 in ((1, r"$h^1$", 3e-4), (3, r"$h^3$", 2.2e-4)):
        hs = [h0, h0 * 2]
        ax[0].loglog(hs, [y0 * (x / h0) ** p for x in hs], "k:", lw=0.8)
        ax[0].text(hs[1] * 1.04, y0 * 2 ** p, lab, fontsize=7, va="center")
    ax[0].set_xticks([1 / 136, 1 / 96, 1 / 68])
    ax[0].set_xticklabels(["1/136", "1/96", "1/68"])
    ax[0].minorticks_off()
    ax[0].set_xlabel("$h$"); ax[0].set_ylabel(r"shape error $e_\mathrm{shape}(T)$")
    ax[0].set_title("(a) shape error at $T = 2$ against $h$")
    for c in S:
        R = [r for r in H if r["candidate"] == c and r["N"] == "136"]
        t = [num(r["TIME"]) for r in R]; g = [num(r["E_NARROW_L2_GRAD_PSI"]) for r in R]
        ax[1].semilogy(t, g, color=COLOR[c], label=c)
    ax[1].axvline(1.0, color="gray", lw=0.6, ls="--")
    ax[1].set_xlabel("$t$"); ax[1].set_ylabel(r"$\Vert q-1\Vert_{2,\mathrm{band}}$")
    ax[1].set_title("(b) band gradient error, $N = 136$ ($T/2$ dashed)")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=7, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, os.path.join(out, "gcls_shear.pdf"))


def fig_droplets(A, out):
    fig, ax = plt.subplots(1, 3, figsize=(6.8, 2.4))
    panels = [("stationary", "l2MagUPrime", r"$\Vert\mathbf{u}\Vert_2$ [m/s]", "(a) stationary, $N = 200$"),
              ("translating", "zeroSetRadialL2", r"$e_\mathrm{shape}\,R$ [m]", "(b) translating, $N = 200$"),
              ("oscillating", "gradPsiL2ErrorBand", r"$\Vert q-1\Vert_{2,\mathrm{band}}$", "(c) oscillating, $N = 200$")]
    for k, (arm, col, lab, title) in enumerate(panels):
        H = rows(os.path.join(A, "gate", f"histories_{arm}.csv"))
        for c in CANDS:
            R = [r for r in H if r["candidate"] == c and r["N"] == "200"]
            t = [num(r["TIME"]) for r in R]; v = [num(r[col]) for r in R]
            pts = [(x, y) for x, y in zip(t, v) if x is not None and y is not None and y > 0]
            if pts:
                ax[k].semilogy([p[0] for p in pts], [p[1] for p in pts], color=COLOR[c], label=c)
        ax[k].set_xlabel("$t$ [s]"); ax[k].set_ylabel(lab); ax[k].set_title(title)
    h, l = ax[2].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=7, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, os.path.join(out, "gcls_droplets.pdf"))


def fig_oscillating_drift(A, out):
    H = rows(os.path.join(A, "gate", "histories_oscillating.csv"))
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
    for N, ls in (("100", ":"), ("142", "--"), ("200", "-")):
        R = [r for r in H if r["candidate"] == "baseline" and r["N"] == N]
        t = [num(r["TIME"]) for r in R]
        ax[0].semilogy(t, [num(r["gradPsiL2ErrorBand"]) for r in R], "k" + ls, label=f"$N = {N}$")
        ax[1].semilogy(t, [num(r["meanMagUPrime"]) for r in R], "k" + ls, label=f"$N = {N}$")
    ax[0].set_ylabel(r"$\Vert q-1\Vert_{2,\mathrm{band}}$"); ax[1].set_ylabel(r"$\Vert\mathbf{u}\Vert_1$ [m/s]")
    for a in ax:
        a.set_xlabel("$t$ [s]")
    ax[0].set_title("(a) baseline: band gradient error"); ax[1].set_title("(b) baseline: velocity (the oscillation)")
    ax[0].legend(frameon=False)
    fig.tight_layout()
    save(fig, os.path.join(out, "gcls_oscillating_drift.pdf"))


def fig_translating_boxes(A, out):
    H = rows(os.path.join(A, "laptop", "translating_histories.csv"))
    runs = [("trans100_np4fix", "10 mm, $N = 100$", "tab:red", "-"),
            ("std142", "10 mm, $N = 142$", "tab:red", "--"),
            ("long100_T25", "20 mm, $N = 100$", "tab:blue", "-"),
            ("long142_T25", "20 mm, $N = 142$", "tab:blue", "--"),
            ("long200", "20 mm, $N = 200$", "tab:blue", ":"),
            ("box40", "40 mm, $N = 100$", "tab:green", "-"),
            ("box40_142", "40 mm, $N = 142$", "tab:green", "--")]
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.6))
    for run, lab, col, ls in runs:
        R = [r for r in H if r["run"] == run]
        t = [num(r["TIME"]) for r in R]
        L = float(R[0]["boxLength_mm"]) if R else 10.0
        # the volume-weighted L2 norm scales with 1/sqrt(domain volume): refer every box to 10 mm
        s = math.sqrt(L / 10.0)
        ax[0].semilogy(t, [num(r["l2MagUPrime"]) * s for r in R], color=col, ls=ls, label=lab)
        ax[1].plot(t, [num(r["centroidLead_mm"]) for r in R], color=col, ls=ls, label=lab)
    ax[0].set_ylabel(r"$\Vert\mathbf{u}-\mathbf{U}_0\Vert_2\,\sqrt{L/10\,\mathrm{mm}}$ [m/s]")
    ax[1].set_ylabel(r"centroid $x_c - (x_0 + U_0 t)$ [mm]")
    ax[1].set_ylim(-0.5, 4.5)
    for a in ax:
        a.set_xlabel("$t$ [s]")
    ax[0].set_title("(a) spurious current, three box lengths"); ax[1].set_title("(b) lead of the droplet over the stream")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save(fig, os.path.join(out, "gcls_translating_boxes.pdf"))


SEAM_COLS = [("zeroSetRadialL2", "shape"), ("gradPsiL2ErrorBand", "gradient band"),
             ("kErrL2Band", "curvature"), ("m2Amplitude", "mode-2 amplitude"),
             ("l2MagUPrime", r"$\Vert\mathbf{u}'\Vert_2$"), ("meanMagUPrime", r"$\Vert\mathbf{u}'\Vert_1$"),
             ("phaseVolumeRelError", "volume"), ("pLaplace", "pressure jump"), ("centroidX", "centroid")]


def fig_seam(A, out):
    T = {r["column"]: r for r in rows(os.path.join(A, "laptop", "seam_np4_vs_serial.csv"))}
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    x = range(len(SEAM_COLS)); w = 0.27
    for k, (key, lab, col) in enumerate((("beforeTheFixes", "before the fixes", "tab:red"),
                                         ("coupledFaceFix", "coupled-face fix", "tab:orange"),
                                         ("bothFixes", "both fixes", "tab:green"))):
        ax.bar([i + (k - 1) * w for i in x], [max(num(T[c][key]), 1e-12) for c, _ in SEAM_COLS],
               width=w, color=col, label=lab)
    ax.set_yscale("log"); ax.set_xticks(list(x)); ax.set_xticklabels([l for _, l in SEAM_COLS], rotation=25, ha="right")
    ax.set_ylabel("np 4 against serial")
    ax.set_ylim(1e-10, 30)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    save(fig, os.path.join(out, "gcls_seam.pdf"))


# --------------------------------------------------------------------------- tables
def tab_candidates(out):
    gate = yaml.safe_load(open(os.path.join(REPO, "config", "gates", "methodGate2D.yaml")))
    tref = {a: float(v["T_REF"]) for a, v in gate["arms"].items() if "T_REF" in v}
    lines = []
    for c in CANDS:
        d = yaml.safe_load(open(os.path.join(REPO, "config", "candidates", f"{c}.yaml")))
        t = d.get("tokens", {}) or {}
        ext = t.get("VELOCITY_EXTENSION", "none")
        if ext == "haloLimited":
            ext = rf"\texttt{{haloLimited}}, $R={t.get('HL_RADIUS_CELLS')}h$, $m={t.get('HL_M')}$, $\beta={t.get('HL_BETA')}$"
        else:
            ext = rf"\texttt{{{ext}}}"
        law = t.get("GC_LAW", "none") if t.get("SL_SOURCE", "none") != "none" else "none"
        if law == "softWall":
            src = (rf"\texttt{{softWall}}, $C_\kappa={t['SW_C_KAPPA']}$, $\delta_s={t['SW_DELTA_S']}$, "
                   rf"$p={t['SW_P']}$, $\gamma=\mathrm{{artanh}}\,0.9$")
        elif law in ("linearQ", "linearZ"):
            m = (d.get("rates") or {}).get("GC_M_MU", 1)
            src = rf"\texttt{{{law}}}, $\mu = {m}/T_\mathrm{{ref}}$"
        else:
            src = "none"
        band = rf"$3h$" if law != "none" else "--"
        lines.append(rf"\texttt{{{c}}} & {ext} & {src} & {band}\\")
    body = "\n".join(lines)
    tr = ", ".join(rf"{a} {v:g}\,s" for a, v in tref.items() if a != "exact1D")
    write_tex(os.path.join(out, "gcls_candidates.tex"), rf"""% Generated by figures/make_result_figures.py from config/candidates and config/gates/methodGate2D.yaml.
\begin{{tabular}}{{@{{}}llll@{{}}}}
  \toprule
  candidate & velocity extension (trace flux) & source law (strain weight \texttt{{none}}) & band\\
  \midrule
{body}
  \bottomrule
\end{{tabular}}
""")
    write_tex(os.path.join(out, "gcls_candidates_tref.tex"),
              rf"% Generated by figures/make_result_figures.py from config/gates/methodGate2D.yaml." "\n"
              rf"$T_\mathrm{{ref}}$ per arm: {tr}.")


def tab_baseline(A, out):
    S = rows(os.path.join(A, "gate", "summaries", "baseline", "summary.csv"))
    O = {(r["arm"], r["metric"]): r for r in rows(os.path.join(A, "gate", "summaries", "baseline", "orders.csv"))}
    spec = {"shear": ["shapeError", "volumeError", "gradientBandError"],
            "stationary": ["shapeError", "volumeError", "gradientBandError", "l2MagUPrime", "pressureJumpError", "curvatureError"],
            "translating": ["shapeError", "volumeError", "gradientBandError", "l2MagUPrime", "pressureJumpError", "curvatureError", "travelledFractionError"],
            "oscillating": ["volumeError", "gradientBandError", "oscPeriod", "oscDampingRate"]}
    label = {"shapeError": r"$e_\mathrm{shape}$", "volumeError": r"$e_V$", "gradientBandError": r"$\lVert q-1\rVert_{2,b}$",
             "l2MagUPrime": r"$\lVert\mathbf u'\rVert_2$", "pressureJumpError": r"$e_{\Delta p}$",
             "curvatureError": r"$e_\kappa$", "travelledFractionError": r"$e_\mathrm{travel}$",
             "oscPeriod": r"$T_\mathrm{osc}$ [s]", "oscDampingRate": r"$\gamma$ [1/s]"}
    lines = []
    for arm, mets in spec.items():
        R = sorted([r for r in S if r["arm"] == arm], key=lambda r: int(r["N"]))
        for m in mets:
            vals = " & ".join(fmt(num(r.get(m))) for r in R)
            o = O.get((arm, m), {})
            order = num(o.get("lsqOrder"))
            if order is None and num(o.get("celikOrder")) is not None:
                oc = f"{num(o['celikOrder']):.2f} ({o.get('convergence', '')[:3]}.)"
            else:
                oc = "--" if order is None else f"{order:.2f}"
            lines.append(rf"{arm} & {label[m]} & {vals} & {oc}\\")
        lines.append(r"\addlinespace")
    write_tex(os.path.join(out, "gcls_baseline.tex"), r"""% Generated by figures/make_result_figures.py from the archived baseline summary.
\begin{tabular}{@{}llrrrr@{}}
  \toprule
  arm & quantity & rung 1 & rung 2 & rung 3 & order\\
  \midrule
""" + "\n".join(lines[:-1]) + r"""
  \bottomrule
\end{tabular}
""")


def tab_translating(A, out):
    R = rows(os.path.join(A, "laptop", "translating_runs.csv"))
    lines = []
    for r in R:
        onset = num(r["growthOnsetTime"])
        lines.append(" & ".join([r["boxLength_mm"], r["N"], r["np"], r["startX_mm"], r["variant"].replace("_", r"\_"),
                                 r["binaries"], r["endTime"], r["state"], f"{num(r['lastTime']):.4f}",
                                 "--" if onset is None else f"{onset:.4f}"]) + r"\\")
    write_tex(os.path.join(out, "gcls_translating_runs.tex"), r"""% Generated by figures/make_result_figures.py from the archived laptop runs.
\begin{tabular}{@{}rrrrllrlrr@{}}
  \toprule
  $L$ [mm] & $N$ & np & $x_0$ [mm] & change & binaries & $T$ [s] & state & last $t$ [s] & onset $t$ [s]\\
  \midrule
""" + "\n".join(lines) + r"""
  \bottomrule
\end{tabular}
""")


def tab_eulerian(A, out):
    R = rows(os.path.join(A, "laptop", "eulerian_mass_flux.csv"))
    lines = [" & ".join([r["solver"], f"{num(r['densityRatio']):g}", r["np"], f"{num(r['travelledFraction']):.5f}",
                         fmt(num(r["L2_U_minus_U0"])), fmt(num(r["L1_U_minus_U0"]))]) + r"\\" for r in R]
    write_tex(os.path.join(out, "gcls_eulerian.tex"), r"""% Generated by figures/make_result_figures.py from the archived laptop runs.
\begin{tabular}{@{}lrrrrr@{}}
  \toprule
  Eulerian solver & $\rho_1/\rho_2$ & np & travelled fraction & $\lVert\mathbf u-\mathbf U_0\rVert_2$ [m/s] & $\lVert\mathbf u-\mathbf U_0\rVert_1$ [m/s]\\
  \midrule
""" + "\n".join(lines) + r"""
  \bottomrule
\end{tabular}
""")


def tab_seam(A, out):
    T = {r["column"]: r for r in rows(os.path.join(A, "laptop", "seam_np4_vs_serial.csv"))}
    lines = [" & ".join([lab, fmt(num(T[c]["beforeTheFixes"])), fmt(num(T[c]["coupledFaceFix"])),
                         fmt(num(T[c]["bothFixes"]))]) + r"\\" for c, lab in SEAM_COLS]
    F = rows(os.path.join(A, "laptop", "seam_face_density.csv"))
    fl = [" & ".join([r["code"], r["time"], r["field"], r["facesCompared"], r["facesDiffering"], fmt(num(r["worstRelative"]))]) + r"\\"
          for r in F if r["field"] in ("alphaf", "rhof")]
    write_tex(os.path.join(out, "gcls_seam.tex"), r"""% Generated by figures/make_result_figures.py from the archived laptop runs.
\begin{tabular}{@{}lrrr@{}}
  \toprule
  quantity & before the fixes & coupled-face fix & both fixes\\
  \midrule
""" + "\n".join(lines) + r"""
  \bottomrule
\end{tabular}
""")
    write_tex(os.path.join(out, "gcls_face_density.tex"), r"""% Generated by figures/make_result_figures.py from the archived laptop runs.
\begin{tabular}{@{}llllrr@{}}
  \toprule
  code & $t$ [s] & face field & faces & differing & worst relative\\
  \midrule
""" + "\n".join(fl) + r"""
  \bottomrule
\end{tabular}
""")


def tab_gate_seam(A, out):
    lines = []
    for c in CANDS:
        p = os.path.join(A, "gate", "summaries", c, "seam.csv")
        if not os.path.isfile(p):
            continue
        cells = {}
        for r in rows(p):
            key = r.get("decomposition", "")
            v = num(r.get("maxColumnScaledDiff"))
            prev = cells.get(key)
            if r.get("verdict") == "NOT_COMPARABLE":
                cells[key] = "n.c."
            elif v is not None and (prev is None or (isinstance(prev, float) and v > prev)):
                cells[key] = v
        def cell(k):
            v = cells.get(k)
            return v if isinstance(v, str) else fmt(v)
        lines.append(rf"\texttt{{{c}}} & {cell('seamNp1')} & {cell('seamNp8')} & {cell('translatingSeamNp1')}\\")
    write_tex(os.path.join(out, "gcls_gate_seam.tex"), r"""% Generated by figures/make_result_figures.py from the archived gate summaries.
\begin{tabular}{@{}lrrr@{}}
  \toprule
  candidate & shear, np 1 & shear, np 8 & translating, np 1\\
  \midrule
""" + "\n".join(lines) + r"""
  \bottomrule
\end{tabular}\\[3pt]
{\footnotesize The column-scaled maximum difference to the np 4 run over all rows; tolerance $10^{-10}$
(shear) and $10^{-5}$ (translating); n.c.: not comparable, the np 4 reference diverged.}
""")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", required=True)
    ap.add_argument("--figures", required=True)
    ap.add_argument("--tables", required=True)
    a = ap.parse_args()
    os.makedirs(a.figures, exist_ok=True); os.makedirs(a.tables, exist_ok=True)
    fig_shear(a.archive, a.figures)
    fig_droplets(a.archive, a.figures)
    fig_oscillating_drift(a.archive, a.figures)
    fig_translating_boxes(a.archive, a.figures)
    fig_seam(a.archive, a.figures)
    tab_candidates(a.tables)
    tab_baseline(a.archive, a.tables)
    tab_translating(a.archive, a.tables)
    tab_eulerian(a.archive, a.tables)
    tab_seam(a.archive, a.tables)
    tab_gate_seam(a.archive, a.tables)
    subprocess.run([sys.executable, os.path.join(REPO, "workflow", "scripts", "make_gate_tables.py"),
                    "--summary-root", os.path.join(a.archive, "gate", "summaries"),
                    "--gate", "methodGate2D", "--out", a.tables], check=True)


if __name__ == "__main__":
    main()
