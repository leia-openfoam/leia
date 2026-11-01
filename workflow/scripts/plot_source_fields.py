#!/usr/bin/env python3
"""Field figures of the source-term studies for the gradient-control deck (STATUS 11.27, 11.28).

    plot_source_fields.py mechanism --out <dir>
    plot_source_fields.py static    --study <staticSourceGcBand2D dir> --out <dir>
    plot_source_fields.py rescale   --once <...RescaleOnce2D> --every <...RescaleEvery2D> --out <dir>
    plot_source_fields.py extension --study <projectionBlend2Dvortex dir> --out <dir>

Reads the reconstructed psi of the uniform N x N unit-square cases (cell (i, j) at
((i + 1/2)/N, (j + 1/2)/N), x fastest) and selects each case by the tokens of its
case_params.json. The figures show psi - d in units of h (static circle: d from the
case's circle) and |grad psi| by central differences (a DISPLAY estimate, not the
solver's). Dark deck surface; the dataviz reference palette (dark steps)."""
import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import SymLogNorm, TwoSlopeNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from static_source_metrics import read_scalar_field, time_dirs  # noqa: E402

SURFACE, TEXT, MUTED, GRID = "#191919", "#f0efec", "#c3c2b7", "#3a3a38"
SLOT = ["#3987e5", "#d95926", "#199e70"]
DIV = "RdBu_r"      # diverging: two hues, neutral midpoint
SEQ = "viridis"     # sequential magnitude


def dark(fig, axes):
    fig.patch.set_facecolor(SURFACE)
    for ax in np.atleast_1d(axes).ravel():
        ax.set_facecolor(SURFACE)
        for s in ax.spines.values():
            s.set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=10)
        ax.xaxis.label.set_color(TEXT)
        ax.yaxis.label.set_color(TEXT)
        ax.title.set_color(TEXT)


def cbar(fig, im, ax, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(colors=MUTED, labelsize=9)
    cb.outline.set_edgecolor(MUTED)
    cb.set_label(label, color=TEXT, fontsize=10)
    return cb


def find_case(study, **want):
    hits = []
    for c in sorted(glob.glob(os.path.join(study, "*_0*"))):
        p = os.path.join(c, "case_params.json")
        if not os.path.isfile(p):
            continue
        tok = json.load(open(p))["tokens"]
        if all(str(tok.get(k)) == str(v) for k, v in want.items()):
            hits.append(c)
    if len(hits) != 1:
        raise SystemExit(f"{study}: {len(hits)} cases match {want}")
    return hits[0]


def tokens(case):
    return json.load(open(os.path.join(case, "case_params.json")))["tokens"]


def field(case, tname, N, name="psi"):
    return read_scalar_field(os.path.join(case, tname, name), N * N).reshape(N, N)


def at_time(case, t):
    times = time_dirs(case)
    return min(times, key=lambda x: abs(x[0] - t))


def circle_d(case):
    tok = tokens(case)
    N = int(float(tok["N_CELLS"]))
    x = (np.arange(N) + 0.5) / N
    X, Y = np.meshgrid(x, x)
    d = np.hypot(X - float(tok["CIRCLE_X0"]), Y - float(tok["CIRCLE_Y0"])) - float(tok["CIRCLE_R"])
    return N, X, Y, d


def save(fig, out, name):
    fig.savefig(os.path.join(out, name), facecolor=SURFACE, dpi=150)
    plt.close(fig)
    print("wrote", os.path.join(out, name))


# --------------------------------------------------------------------------- mechanism (DERIVED)
def cmd_mechanism(a):
    """1D: psi_i = d_i (1 + E s_i), s_i = (-1)^i, the mode (-1)^i d. Centred and one-sided q per cell,
    and the sign of the linearQ factor F = mu (1 - q)."""
    h, E = 1.0, 0.25
    i = np.arange(-6, 7)
    d = (i + 0.4) * h                         # zero set between cells -1 and 0
    s = np.where(i % 2 == 0, 1.0, -1.0)
    psi = d * (1 + E * s)
    qc = np.full_like(psi, np.nan)
    qc[1:-1] = np.abs(psi[2:] - psi[:-2]) / (2 * h)
    qu = np.full_like(psi, np.nan)             # one-sided: the neighbour toward the zero set
    for k in range(1, len(i) - 1):
        up = k - 1 if d[k] > 0 else k + 1
        qu[k] = abs(psi[k] - psi[up]) / h
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 1]))
    dark(fig, [ax1, ax2])
    ax1.plot(d, d, "-", color=MUTED, lw=1.5, label=r"$d$")
    ax1.plot(d, psi, "o", color=SLOT[0], ms=9, label=r"$\psi_i=d_i(1+E\,(-1)^i)$, $E=0.25$")
    for k in range(1, len(i) - 1):
        F = 1 - qc[k]                          # linearQ with mu = 1, centred q
        if abs(F) > 1e-12 and abs(d[k]) > 0.5:
            ax1.annotate("", xy=(d[k], psi[k] + 0.9 * np.sign(F * psi[k])), xytext=(d[k], psi[k]),
                         arrowprops=dict(arrowstyle="->", color=SLOT[1], lw=1.6))
    ax1.axhline(0, color=GRID, lw=1)
    ax1.set_ylabel(r"$\psi$  (arrows: $\psi F$, linearQ, centred $q$)")
    leg = ax1.legend(frameon=False, fontsize=10, loc="upper left")
    [t.set_color(TEXT) for t in leg.get_texts()]
    ax2.axhline(1, color=MUTED, lw=1, ls="--")
    ax2.plot(d, qc, "o-", color=SLOT[1], ms=8, lw=1.5, label=r"centred $q_i=|\psi_{i+1}-\psi_{i-1}|/2h$")
    ax2.plot(d, qu, "s-", color=SLOT[2], ms=7, lw=1.5, label=r"one-sided $q_i=|\psi_i-\psi_{up}|/h$")
    ax2.set_xlabel(r"$d/h$")
    ax2.set_ylabel(r"$q$")
    leg = ax2.legend(frameon=False, fontsize=10, loc="upper left", ncol=1)
    [t.set_color(TEXT) for t in leg.get_texts()]
    fig.tight_layout()
    save(fig, a.out, "mech_centred_q.png")


# --------------------------------------------------------------------------- static test (MEASURED)
def cmd_static(a):
    c = find_case(a.study, SDPLS_SOURCE="gradientControl", GC_LAW="linearQ", INIT_FIELD="sdf",
                  N_CELLS=128, MOLLIFIER="band", GC_STRAIN_WEIGHT="none")
    N, X, Y, d = circle_d(c)
    h = 1.0 / N
    taus = [20, 40, 60, 80]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    dark(fig, axes)
    win = (X > 0.62) & (X < 0.82) & (Y > 0.40) & (Y < 0.60)
    for ax, t in zip(axes, taus):
        tv, tn = at_time(c, t)
        e = (field(c, tn, N) - d) / h
        lim = max(np.abs(e[win]).max(), 1e-12)
        im = ax.pcolormesh(X, Y, e, cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim), shading="nearest")
        ax.contour(X, Y, d, [0], colors=[TEXT], linewidths=0.8)
        ax.set_xlim(0.62, 0.82); ax.set_ylim(0.40, 0.60); ax.set_aspect("equal")
        ax.set_title(rf"$\tau={tv:g}$", fontsize=12)
        cbar(fig, im, ax, r"$(\psi-d)/h$")
    fig.suptitle(f"linearQ, band of 3 cells, start d, N = {N} (window near the right of the circle; "
                 "white: exact circle)", color=TEXT, fontsize=12)
    fig.tight_layout()
    save(fig, a.out, "static_linearQ_fields.png")

    # a row cut through the centre
    j = int(float(tokens(c)["CIRCLE_Y0"]) * N)
    x = X[j]
    sel = (x > 0.64) & (x < 0.86)
    fig, ax = plt.subplots(figsize=(9, 4.6))
    dark(fig, ax)
    for t, col in zip((40, 60, 80), SLOT):
        tv, tn = at_time(c, t)
        e = (field(c, tn, N)[j] - d[j]) / h
        ax.plot(x[sel], e[sel], "o-", color=col, ms=5, lw=1.2, label=rf"$\tau={tv:g}$")
    ax.axvline(float(tokens(c)["CIRCLE_X0"]) + float(tokens(c)["CIRCLE_R"]), color=MUTED, ls="--", lw=1)
    ax.text(0.753, ax.get_ylim()[1] * 0.9, " zero set", color=MUTED, fontsize=10)
    ax.set_xlabel("x (cell centres of the row through the circle centre)")
    ax.set_ylabel(r"$(\psi-d)/h$")
    ax.grid(True, color=GRID, lw=0.6)
    leg = ax.legend(frameon=False, fontsize=10)
    [t.set_color(TEXT) for t in leg.get_texts()]
    fig.tight_layout()
    save(fig, a.out, "static_linearQ_cut.png")


# --------------------------------------------------------------------------- rescaling (MEASURED)
def _resc(study, order, events, mag, init, N, hess="SdplsLSQ"):
    return find_case(study, REDISTANCER="rescaling", RESCALE_ORDER=order, RESCALE_MAX_EVENTS=events,
                     RESCALE_MAG_GRAD=mag, RESCALE_HESSIAN_SCHEME=hess, INIT_FIELD=init, N_CELLS=N)


def cmd_rescale(a):
    N = 128
    # once: the normalised algebraic start, then (0a) and (0b)
    c1 = _resc(a.once, "first", 1, "centred", "algNorm", N)
    c2 = _resc(a.once, "second", 1, "centred", "algNorm", N)
    _, X, Y, d = circle_d(c1)
    h = 1.0 / N
    panels = [(c1, 0, "start: $(r^2-R^2)/2R$"), (c1, 1, "after (0a), centred"), (c2, 1, "after (0b), centred")]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    dark(fig, axes)
    lim = 0.5
    for ax, (c, t, lab) in zip(axes, panels):
        tv, tn = at_time(c, t)
        e = (field(c, tn, N) - d) / h
        band = np.abs(d) <= 3 * h
        l2 = np.sqrt(np.mean(e[band] ** 2))
        im = ax.pcolormesh(X, Y, np.where(np.abs(d) <= 6 * h, e, np.nan), cmap=DIV,
                           norm=TwoSlopeNorm(0, -lim, lim), shading="nearest")
        ax.contour(X, Y, d, [0], colors=[TEXT], linewidths=0.6)
        ax.set_xlim(0.2, 0.8); ax.set_ylim(0.2, 0.8); ax.set_aspect("equal")
        ax.set_title(f"{lab}\nL2 in |d| <= 3h: {l2:.2e} h", fontsize=11)
        cbar(fig, im, ax, r"$(\psi-d)/h$ in $|d|\leq 6h$")
    fig.tight_layout()
    save(fig, a.out, "rescale_once_fields.png")

    # every step: centred second order at tau = 1, one-sided first order at tau = 1 and 10
    ce = _resc(a.every, "second", -1, "centred", "sdf", N)
    ge = _resc(a.every, "first", -1, "godunov", "sdf", N)
    panels = [(ce, 1, "centred (0b), every step"), (ge, 1, "one-sided (0a), every step"),
              (ge, 10, "one-sided (0a), every step")]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    dark(fig, axes)
    for k, (ax, (c, t, lab)) in enumerate(zip(axes, panels)):
        tv, tn = at_time(c, t)
        e = (field(c, tn, N) - d) / h
        show = np.where(np.abs(d) <= 6 * h, e, np.nan)
        if k == 0:
            norm = SymLogNorm(linthresh=1.0, vmin=-np.nanmax(np.abs(show)), vmax=np.nanmax(np.abs(show)))
        else:
            norm = TwoSlopeNorm(0, -3, 3)
        im = ax.pcolormesh(X, Y, show, cmap=DIV, norm=norm, shading="nearest")
        ax.contour(X, Y, d, [0], colors=[TEXT], linewidths=0.6)
        ax.set_xlim(0.62, 0.82); ax.set_ylim(0.40, 0.60); ax.set_aspect("equal")
        ax.set_title(rf"{lab}, $\tau={tv:g}$", fontsize=11)
        cbar(fig, im, ax, r"$(\psi-d)/h$" + (" (symlog)" if k == 0 else ""))
    fig.suptitle(f"start d, N = {N}, window near the right of the circle", color=TEXT, fontsize=12)
    fig.tight_layout()
    save(fig, a.out, "rescale_every_fields.png")


# --------------------------------------------------------------------------- extension (MEASURED)
def cmd_extension(a):
    N = 256
    base = find_case(a.study, VELOCITY_EXTENSION="none", PB_DISTANCE="second", SDPLS_SOURCE="noSource", N_CELLS=N)
    ext = find_case(a.study, VELOCITY_EXTENSION="projectionBlend", PB_DISTANCE="second", SDPLS_SOURCE="noSource", N_CELLS=N)
    T = float(tokens(base)["END_TIME"])
    h = 1.0 / N
    x = (np.arange(N) + 0.5) / N
    X, Y = np.meshgrid(x, x)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
    dark(fig, axes)
    for ax, (c, lab) in zip(axes, ((base, "no extension"), (ext, "projectionBlend, second distance"))):
        tv, tn = at_time(c, T / 2)
        p = field(c, tn, N)
        gy, gx = np.gradient(p, h)
        q = np.hypot(gx, gy)
        im = ax.pcolormesh(X, Y, np.where(np.abs(p) <= 6 * h, q, np.nan), cmap=SEQ, vmin=0, vmax=3,
                           shading="nearest")
        ax.contour(X, Y, p, [0], colors=[TEXT], linewidths=0.7)
        ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(rf"{lab}, $t={tv:.3f}$ (T/2)", fontsize=12)
        cbar(fig, im, ax, r"$|\nabla\psi|$ in $|\psi|\leq 6h$ (central differences)")
    fig.suptitle(f"2D vortex, N = {N}: the gradient at the largest deformation", color=TEXT, fontsize=13)
    fig.tight_layout()
    save(fig, a.out, "extension_vortex_gradient.png")

    # shape at T: the zero contours against the initial circle
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    dark(fig, ax)
    t0, n0 = at_time(base, 0)
    ax.contour(X, Y, field(base, n0, N), [0], colors=[MUTED], linewidths=1.2, linestyles="--")
    for c, col, lab in ((base, SLOT[0], "no extension"), (ext, SLOT[1], "projectionBlend")):
        tv, tn = at_time(c, T)
        ax.contour(X, Y, field(c, tn, N), [0], colors=[col], linewidths=1.4)
        ax.plot([], [], color=col, label=f"{lab}, t = {tv:g}")
    ax.plot([], [], color=MUTED, ls="--", label="initial circle")
    ax.set_aspect("equal"); ax.set_xlim(0.33, 0.67); ax.set_ylim(0.58, 0.92)
    ax.set_title(f"2D vortex, N = {N}: the zero level set at T", fontsize=12)
    leg = ax.legend(frameon=False, fontsize=10, loc="lower center")
    [t.set_color(TEXT) for t in leg.get_texts()]
    fig.tight_layout()
    save(fig, a.out, "extension_vortex_shape.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mechanism"); m.add_argument("--out", required=True)
    s = sub.add_parser("static"); s.add_argument("--study", required=True); s.add_argument("--out", required=True)
    r = sub.add_parser("rescale"); r.add_argument("--once", required=True); r.add_argument("--every", required=True)
    r.add_argument("--out", required=True)
    e = sub.add_parser("extension"); e.add_argument("--study", required=True); e.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    dict(mechanism=cmd_mechanism, static=cmd_static, rescale=cmd_rescale, extension=cmd_extension)[a.cmd](a)


if __name__ == "__main__":
    main()
