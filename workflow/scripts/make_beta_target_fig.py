#!/usr/bin/env python3
"""Is the sdplsBeta gradient error a RESIDUAL or a wrong TARGET?

sdplsBeta drives the level set with f_nl = beta - |grad psi|. That is a
relaxation, and it competes with the interfacial normal strain
a = n . grad(U) . n, so the interfacial gradient relaxes toward

    g* = beta - a

-- a target set by the FLOW, not by the mesh. If that is what is happening, no
amount of refinement can fix it, and the fix is not a finer grid but a different
source. The rival explanation is the ordinary one: the offset is a
discretization residual and vanishes as h -> 0.

The two are distinguishable, and this script runs the three tests that
distinguish them (`R`, which has no target at all, is the control arm):

  (i)   h-INDEPENDENCE. A residual shrinks with h; a target offset does not.
  (ii)  ONE-FOR-ONE SHIFT. d(band mean |grad psi|)/d(beta) = 1 if the target is
        beta - a. This is the decisive test and the reason the sweep exists.
  (iii) SPREAD INVARIANCE. The band spread (max - min) is set by the RANGE of a
        over the interface, so shifting beta must not widen it.

and one direct check that needs no inference at all:

  (iv)  RESIDUAL vs the PREDICTED TARGET. gradPsiErrorCSV now reports the band
        statistics of the strain field `R` itself, so `beta - meanStrainBand` is
        a number we can subtract from the measured `meanGradPsiBand`. If the
        mechanism is right, that residual is small and does not grow with beta.

and one consequence that is sharper than any of them:

  (v)   UNREACHABLE TARGET. Taking the gradient of the transport equation and
        evaluating at psi = 0 gives, for g = |grad psi| and a = n.grad(U).n,

            D g / Dt = g (f_nl - a)

        so `R` (f_nl = a) cancels stretching EXACTLY -- its error is purely
        numerical, which is why it converges -- while `beta` (f_nl = beta - g)
        has the fixed point g* = beta - a. Wherever a > beta that fixed point is
        NEGATIVE and a magnitude cannot reach it, so g is driven toward ZERO:
        the level set FLATTENS exactly where the flow stretches hardest. The
        band strain measured on this benchmark spans roughly +-2.9 against
        beta = 1, so this is not a corner case. maxStrainBand > beta is a direct
        test, and it predicts that the damage concentrates at t = T/2 and shows
        up in the band MINIMUM -- which is what the min collapsing to ~0.42
        while the mean sits near 0.84 already looks like.

Falsification is a real outcome and is reported as such: if the error is
minimised at some intermediate beta, or the slope in (ii) is not ~1, the
g* = beta - a story is wrong and the script says so.

Reads the curated CSVs (no re-running):
    studies/sdplsBetaSweep/sdplsBetaSweep_errors.csv     the sweep
    studies/benchVortexEulerT2/benchVortexEulerT2_errors.csv   R / noSource control

Writes into the sdpls theme data/:
    figures/sdpls_beta_target.png
    tables/sdpls_beta_target.csv, tables/sdpls_beta_target.tex

    python3 make_beta_target_fig.py
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import method_label
import paths

THEME = "sdpls-level-set"
SWEEP = "sdplsBetaSweep"
CONTROL = "benchVortexEulerT2"
# t = T/2 is the meaningful meter on a reversed flow: at t = T the integral of
# the strain over the period is zero, which flatters every arm including the
# sourceless one. See the SDPLS article's "Why T/2".
MEAN, MINC, MAXC, STRAIN = ("meanGradPsiBandHalf", "minGradPsiBandHalf",
                            "maxGradPsiBandHalf", "meanStrainBandHalf")
# The largest interfacial stretching anywhere in the band. Where a > beta the
# fixed point g* = beta - a is negative and a magnitude cannot reach it, so g is
# driven to zero instead -- test (v).
STRAIN_MAX = "maxStrainBandHalf"


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(v) else v


def _rows(study):
    p = os.path.join(paths.REPO, "studies", study, f"{study}_errors.csv")
    if not os.path.isfile(p):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def _beta_of(r):
    """The sweep's beta. Prefer the explicit column; fall back to parsing the
    method label so this still works on a CSV written before sdplsBeta existed."""
    b = _f(r.get("sdplsBeta"))
    if b is not None:
        return b
    m = r.get("method", "")
    if "SDPLS:beta(" in m:
        return _f(m.split("SDPLS:beta(", 1)[1].split(")", 1)[0])
    return 1.0 if "SDPLS:beta" in m else None


def collect():
    """{beta: [(h, mean, min, max, strain), ...]} for the sweep, sorted by h."""
    out = {}
    for r in _rows(SWEEP):
        if r.get("solverFailed"):
            continue
        b, h = _beta_of(r), _f(r.get("h"))
        mean, lo, hi = _f(r.get(MEAN)), _f(r.get(MINC)), _f(r.get(MAXC))
        if b is None or h is None or mean is None:
            continue
        out.setdefault(b, []).append((h, mean, lo, hi, _f(r.get(STRAIN)),
                                      _f(r.get(STRAIN_MAX))))
    for b in out:
        out[b].sort(key=lambda t: t[0])
    return out


def control():
    """The `R` arm of the control study: (h, mean) sorted by h."""
    pts = []
    for r in _rows(CONTROL):
        if r.get("solverFailed") or "SDPLS:R" not in r.get("method", ""):
            continue
        h, mean = _f(r.get("h")), _f(r.get(MEAN))
        if h is not None and mean is not None:
            pts.append((h, mean))
    return sorted(pts)


def main():
    data = collect()
    if not data:
        print(f"[beta_target] no usable rows in studies/{SWEEP}; "
              f"run `make studies-one STUDY={SWEEP}` first")
        return 1
    betas = sorted(data)
    figs, tabs = paths.figs_dir(THEME), paths.tables_dir(THEME)

    # ---- (ii) one-for-one shift, measured at the FINEST common h ------------
    finest = max(min(h for h, *_ in v) for v in data.values())
    at_fine = {}
    for b, v in data.items():
        row = min(v, key=lambda t: abs(t[0] - finest))
        at_fine[b] = row
    xs = np.array(betas)
    ys = np.array([at_fine[b][1] for b in betas])
    slope, icept = (np.polyfit(xs, ys, 1) if len(xs) > 1 else (float("nan"),)*2)

    # ---- (iii) spread, and (iv) the predicted-target residual ---------------
    rows = []
    for b in betas:
        h, mean, lo, hi, strain, smax = at_fine[b]
        spread = (hi - lo) if (hi is not None and lo is not None) else None
        pred = (b - strain) if strain is not None else None
        # (v) the target is unreachable wherever a > beta.
        unreachable = (smax is not None and smax > b)
        rows.append(dict(beta=b, h=h, mean=mean, min=lo, max=hi,
                         spread=spread, strain=strain, predicted=pred,
                         residual=(mean - pred) if pred is not None else None,
                         maxStrain=smax, targetUnreachable=int(unreachable)))

    # ---- (i) h-independence: how much does the mean move over the ladder? ---
    drift = {}
    for b, v in data.items():
        m = [t[1] for t in v]
        drift[b] = (max(m) - min(m)) if len(m) > 1 else float("nan")
    ctrl = control()
    ctrl_drift = (max(m for _h, m in ctrl) - min(m for _h, m in ctrl)) \
        if len(ctrl) > 1 else float("nan")

    # ------------------------------- report ---------------------------------
    print(f"[beta_target] sweep: {len(data)} beta value(s), "
          f"{sum(len(v) for v in data.values())} rows; finest common h = {finest:g}")
    print()
    print("  (ii) band mean |grad psi| at T/2 vs beta, at the finest common h")
    print(f"       {'beta':>6} {'mean':>9} {'spread':>9} {'strain a':>10} "
          f"{'beta - a':>10} {'mean-(beta-a)':>14}")
    for r in rows:
        def s(x, w, p=4):
            return f"{x:>{w}.{p}f}" if x is not None else f"{'--':>{w}}"
        print(f"       {r['beta']:>6.2f} {s(r['mean'],9)} {s(r['spread'],9)} "
              f"{s(r['strain'],10)} {s(r['predicted'],10)} {s(r['residual'],14)}")
    print()
    print(f"       d(mean)/d(beta) = {slope:.4f}   "
          f"(1.0 => the target IS beta - a; 0.0 => beta does not set the level)")
    print()
    print("  (i) drift of the band mean across the whole ladder "
          "(small => h-INDEPENDENT => a target, not a residual)")
    for b in betas:
        print(f"       beta={b:<5.2f} drift={drift[b]:.4f}")
    print(f"       control R (no target): drift={ctrl_drift:.4f} over "
          f"{len(ctrl)} resolutions")
    print()
    print("  (v) unreachable target: max band strain a vs beta. Where a > beta the")
    print("      fixed point beta - a is NEGATIVE, so |grad psi| is driven to 0")
    print("      (flattening) instead of to a finite level.")
    for r in rows:
        if r["maxStrain"] is None:
            print(f"       beta={r['beta']:<5.2f} max a = --   (strain not recorded)")
            continue
        mark = "UNREACHABLE over part of the band" if r["targetUnreachable"] \
            else "target reachable everywhere"
        print(f"       beta={r['beta']:<5.2f} max a={r['maxStrain']:>7.3f}  "
              f"band min |grad psi|={r['min'] if r['min'] is None else round(r['min'], 4)}"
              f"   -> {mark}")
    print()
    # State the verdict, including the falsifying outcome.
    if np.isfinite(slope):
        if abs(slope - 1.0) < 0.2:
            print(f"  VERDICT: slope {slope:.3f} ~ 1 -- consistent with "
                  f"g* = beta - a. `beta` chases a flow-set target, so refining "
                  f"the mesh cannot fix it.")
        elif abs(slope) < 0.2:
            print(f"  VERDICT: slope {slope:.3f} ~ 0 -- beta does NOT set the "
                  f"equilibrium level. The g* = beta - a mechanism is FALSIFIED; "
                  f"the offset comes from somewhere else.")
        else:
            print(f"  VERDICT: slope {slope:.3f} is neither ~1 nor ~0 -- beta "
                  f"moves the level but not one-for-one; the simple relaxation "
                  f"picture is incomplete.")

    # -------------------------------- CSV -----------------------------------
    csv_path = os.path.join(tabs, "sdpls_beta_target.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["beta", "h", "mean", "min", "max",
                                           "spread", "strain", "predicted",
                                           "residual", "maxStrain",
                                           "targetUnreachable"])
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})
    print(f"[beta_target] wrote {csv_path}")

    # ------------------------------- LaTeX ----------------------------------
    tex = os.path.join(tabs, "sdpls_beta_target.tex")
    with open(tex, "w") as fh:
        fh.write("% generated by workflow/scripts/make_beta_target_fig.py\n")
        fh.write("\\begin{tabular}{rrrrrrrc}\n\\toprule\n")
        fh.write(r"$\beta$ & band mean $|\nabla\psi|$ & spread & "
                 r"$\bar a$ & $\beta-\bar a$ & residual & $\max a$ & "
                 r"$g^*<0$ somewhere \\" "\n\\midrule\n")
        for r in rows:
            def c(x):
                return f"{x:.4f}" if x is not None else "--"
            fh.write(f"{r['beta']:.2f} & {c(r['mean'])} & {c(r['spread'])} & "
                     f"{c(r['strain'])} & {c(r['predicted'])} & "
                     f"{c(r['residual'])} & {c(r['maxStrain'])} & "
                     f"{'yes' if r['targetUnreachable'] else 'no'} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[beta_target] wrote {tex}")

    # ------------------------------- figure ---------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    # left: band mean vs h, one line per beta, with the R control
    for b in betas:
        v = data[b]
        ax[0].semilogx([t[0] for t in v], [t[1] for t in v], "o-",
                       label=fr"$\beta$={b:g}")
        ax[0].axhline(b, ls=":", lw=0.8, color="grey")
    if ctrl:
        ax[0].semilogx([h for h, _ in ctrl], [m for _, m in ctrl], "k^--",
                       label="R (control)")
    ax[0].axhline(1.0, color="k", lw=0.8)
    ax[0].set_xlabel("h"); ax[0].set_ylabel(r"band mean $|\nabla\psi|$ at $T/2$")
    ax[0].set_title("(i) h-independence\n(flat = a target, not a residual)",
                    fontsize=10)
    ax[0].legend(fontsize=7); ax[0].invert_xaxis()

    # middle: the decisive one-for-one test
    ax[1].plot(xs, ys, "o", ms=8, label="measured")
    if np.isfinite(slope):
        xx = np.linspace(min(xs), max(xs), 10)
        ax[1].plot(xx, slope*xx + icept, "-",
                   label=fr"fit, slope={slope:.3f}")
    ax[1].plot(xs, xs, "k:", label=r"slope 1 ($g^*=\beta-a$)")
    ax[1].set_xlabel(r"$\beta$")
    ax[1].set_ylabel(r"band mean $|\nabla\psi|$ at $T/2$")
    ax[1].set_title("(ii) one-for-one shift\n(decisive test)", fontsize=10)
    ax[1].legend(fontsize=8)

    # right: spread invariance + predicted-target residual
    sp = [r["spread"] for r in rows]
    rs = [r["residual"] for r in rows]
    if any(s is not None for s in sp):
        ax[2].plot(xs, [np.nan if s is None else s for s in sp], "s-",
                   label="band spread (max-min)")
    if any(r is not None for r in rs):
        ax[2].plot(xs, [np.nan if r is None else r for r in rs], "d-",
                   label=r"mean $-\,(\beta-a)$")
    ax[2].axhline(0, color="k", lw=0.8)
    ax[2].set_xlabel(r"$\beta$")
    ax[2].set_title("(iii) spread invariance and\n(iv) residual vs predicted target",
                    fontsize=10)
    ax[2].legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(figs, "sdpls_beta_target.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[beta_target] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
