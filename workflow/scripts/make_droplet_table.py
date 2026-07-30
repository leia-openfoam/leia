#!/usr/bin/env python3
"""Stationary-droplet parasitic-current table + log-log convergence figure.

Reads every case of the ``stationaryDroplet2D`` study (each case's
``leiaSemiLagrangianLevelSetTwoPhaseFoam.csv`` = per-step TIME,maxMagU,meanMagU,
pLaplace, and ``case_params.json`` for N_CELLS), reduces the time series to a
*settled* value (mean over the last 20% of the run, after the initial transient
has decayed), and emits into the semi-Lagrangian theme's data source:

    data/tables/droplet_parasitic.csv    one row per resolution
    data/tables/droplet_parasitic.tex    booktabs body (\\input by the article),
                                          with the log-log fit order in a note row
    data/figures/sl_droplet_parasitic.png  log-log max|U| & mean|U| vs h (fitted
                                          slope) + the Laplace jump vs the exact
                                          sigma/R

The parasitic-current convergence order is the least-squares slope of
log(settled max|U|) vs log(h); its R^2 is reported so the reader can judge the fit.

Usage (from repo root):
    python3 workflow/scripts/make_droplet_table.py studies/stationaryDroplet2D
"""
import csv
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths  # thematic docs layout (single source of truth for output dirs)

# Benchmark constants (arXiv:2212.02904, Sec. 3.3.1), fixed for this case.
BOX_L = 0.01        # domain edge length [m]  -> h = BOX_L / N
SIGMA = 0.07274     # surface tension [N/m]
RADIUS = 1.0e-3     # drop radius [m]
DP_EXACT = SIGMA / RADIUS   # Young-Laplace jump = 72.74 Pa

THEME = "semi-lagrangian-level-set"
CSV_NAME = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
SETTLE_FRAC = 0.2   # average over the last 20% of the run (settled currents)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _read_series(path):
    """Return (time, maxU, meanU, pLap) numpy arrays from a metrics CSV."""
    t, mx, mn, pl = [], [], [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            ti = _f(row.get("TIME"))
            if ti is None:
                continue
            t.append(ti)
            mx.append(_f(row.get("maxMagU")) or 0.0)
            mn.append(_f(row.get("meanMagU")) or 0.0)
            pl.append(_f(row.get("pLaplace")) or 0.0)
    return (np.array(t), np.array(mx), np.array(mn), np.array(pl))


def _settled(t, y):
    """Mean of y over the last SETTLE_FRAC of the time span (settled value)."""
    if t.size == 0:
        return None
    t_cut = t.max() - SETTLE_FRAC * (t.max() - t.min())
    sel = t >= t_cut
    return float(y[sel].mean()) if sel.any() else float(y[-1])


def _fit(h, y):
    """Least-squares slope of log y vs log h, plus R^2. None if < 2 valid points."""
    m = [(a, b) for a, b in zip(h, y) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None
    lh, ly = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, ly, 1)
    resid = ly - np.polyval(p, lh)
    ss_res = float((resid**2).sum())
    ss_tot = float(((ly - ly.mean())**2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(p[0]), r2


def _case_rows(study_dir):
    """One reduced record per case, sorted by N (coarse -> fine)."""
    recs = []
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        case_dir = os.path.dirname(meta_path)
        csv_path = os.path.join(case_dir, CSV_NAME)
        if not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
            print(f"[droplet] no {CSV_NAME} in {os.path.basename(case_dir)}; skip")
            continue
        with open(meta_path) as fh:
            tokens = json.load(fh).get("tokens", {})
        n = _f(tokens.get("N_CELLS"))
        t_end = _f(tokens.get("END_TIME"))
        if not n:
            continue
        t, mx, mn, pl = _read_series(csv_path)
        if t.size == 0:
            print(f"[droplet] empty series in {os.path.basename(case_dir)}; skip")
            continue
        settled_mx = _settled(t, mx)
        early = float(mx[:max(1, t.size // 10)].max())   # initial-transient peak
        # A resolution is "settled/convergent" if the parasitic current decayed
        # below its initial transient; "diverged" if it instead grew above it (the
        # fine-mesh transport-coupled instability). Only the settled points define a
        # convergence order; the diverged ones are shown on the figure as a distinct
        # marker (their last-window magnitude), never as a settled value.
        diverged = bool(settled_mx is not None and settled_mx > early)
        # Completeness gate: a still-SETTLING case must reach (near) the end time,
        # else its not-yet-settled value would be misleading. A diverged case is
        # exempt -- it has already blown up, so its state is decided regardless of
        # whether it reached t_end (and it never will settle).
        if t_end and t.max() < 0.95 * t_end and not diverged:
            print(f"[droplet] {os.path.basename(case_dir)} at t={t.max():.4g} "
                  f"< 0.95*{t_end} and still settling; skip")
            continue
        recs.append({
            "N": int(n),
            "h": BOX_L / n,
            "maxU": settled_mx,
            "meanU": _settled(t, mn),
            "maxUpeak": float(mx.max()),
            "pLaplace": _settled(t, pl),
            "diverged": diverged,
        })
    recs.sort(key=lambda r: r["N"])
    for r in recs:
        r["pRelErr"] = abs(r["pLaplace"] - DP_EXACT) / DP_EXACT
    return recs


def main(argv):
    if not argv:
        print("usage: make_droplet_table.py <study_dir>")
        return 1
    study_dir = argv[0]
    recs = _case_rows(study_dir)
    if not recs:
        print(f"[droplet] no completed cases under {study_dir}; nothing written")
        return 1

    conv = [r for r in recs if not r["diverged"]]   # settled, define the order
    div = [r for r in recs if r["diverged"]]         # grow in time (fine-mesh)
    hc = [r["h"] for r in conv]
    order_max, r2_max = _fit(hc, [r["maxU"] for r in conv])
    order_mean, r2_mean = _fit(hc, [r["meanU"] for r in conv])

    figs_dir = paths.figs_dir(THEME)
    tables_dir = paths.tables_dir(THEME)

    # --- CSV (every case, with a settled/diverged state column) --------------
    csv_path = os.path.join(tables_dir, "droplet_parasitic.csv")
    cols = ["N", "h", "maxU", "meanU", "maxUpeak", "pLaplace", "pRelErr", "state"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in recs:
            row = {c: r.get(c) for c in cols}
            row["state"] = "diverged" if r["diverged"] else "settled"
            w.writerow(row)
    print(f"[droplet] wrote {csv_path} "
          f"({len(conv)} settled, {len(div)} diverged)")

    # --- LaTeX booktabs body (\input by the article): SETTLED resolutions only
    # (a diverged case has no meaningful settled value; it appears on the figure).
    tex_path = os.path.join(tables_dir, "droplet_parasitic.tex")
    with open(tex_path, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_droplet_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{rccccc}\n\\toprule\n")
        fh.write("$N$ & $h$ [m] & $\\max|\\mathbf{u}|$ [m/s] & "
                 "$\\overline{|\\mathbf{u}|}$ [m/s] & $\\Delta p$ [Pa] & "
                 "rel.\\ err. \\\\\n\\midrule\n")
        for r in conv:
            fh.write(
                f"{r['N']} & {r['h']:.3e} & {r['maxU']:.3e} & {r['meanU']:.3e} & "
                f"{r['pLaplace']:.2f} & {100*r['pRelErr']:.1f}\\% \\\\\n")
        fh.write("\\midrule\n")
        if order_max is None:
            note = "rate of $\\max|\\mathbf{u}|$ in $h$: --"
        elif len(conv) < 3:
            note = (f"two-resolution rate of $\\max|\\mathbf{{u}}|$ in $h$: "
                    f"$p \\approx {order_max:.1f}$")
        else:
            note = (f"least-squares order of $\\max|\\mathbf{{u}}|$ in $h$: "
                    f"$p = {order_max:.2f}$ ($R^2 = {r2_max:.3f}$)")
        if div:
            dl = ", ".join(f"$N{{=}}{r['N']}$" for r in div)
            note += f"; {dl} diverge (Figure)"
        fh.write("\\multicolumn{6}{l}{\\footnotesize " + note +
                 f"; exact $\\sigma/R = {DP_EXACT:.2f}$ Pa.}} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[droplet] wrote {tex_path}")

    # --- figure: convergent trend (settled) + diverged markers ---------------
    fig, (axU, axP) = plt.subplots(1, 2, figsize=(9.5, 3.8))

    if conv:
        axU.loglog(hc, [r["maxU"] for r in conv], "o-", color="#0072B2",
                   label=r"$\max|\mathbf{u}|$ (settled)")
        axU.loglog(hc, [r["meanU"] for r in conv], "s--", color="#009E73",
                   label=r"$\overline{|\mathbf{u}|}$ (settled)")
    if order_max is not None and conv:
        hf = np.array(sorted(hc))
        yf = conv[-1]["maxU"] * (hf / conv[-1]["h"])**order_max
        axU.loglog(hf, yf, "k:", lw=1, label=rf"$\propto h^{{{order_max:.2f}}}$")
    if div:
        axU.loglog([r["h"] for r in div], [r["maxU"] for r in div], "X",
                   color="#D55E00", ms=11, ls="none",
                   label="diverged (grows in time)")
    axU.set_xlabel(r"$h$ [m]")
    axU.set_ylabel("spurious current [m/s]")
    axU.set_title("Parasitic currents vs resolution")
    axU.grid(True, which="both", ls=":", alpha=0.5)
    axU.legend(frameon=False, fontsize=8)

    if conv:
        axP.semilogx(hc, [r["pLaplace"] for r in conv], "o-", color="#0072B2",
                     label=r"measured $\Delta p$")
    axP.axhline(DP_EXACT, color="k", ls="--", lw=1,
                label=rf"exact $\sigma/R = {DP_EXACT:.2f}$ Pa")
    axP.set_xlabel(r"$h$ [m]")
    axP.set_ylabel(r"Laplace jump $\Delta p$ [Pa]")
    axP.set_title("Young--Laplace jump")
    axP.grid(True, which="both", ls=":", alpha=0.5)
    axP.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    png_path = os.path.join(figs_dir, "sl_droplet_parasitic.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    print(f"[droplet] wrote {png_path}")

    if order_max is not None:
        print(f"[droplet] max|U| order p={order_max:.3f} (R^2={r2_max:.3f}); "
              f"mean|U| order p={order_mean:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
