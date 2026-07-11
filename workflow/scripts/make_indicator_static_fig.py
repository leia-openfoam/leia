#!/usr/bin/env python3
"""Static (t=0) phase-indicator volume convergence to the EXACT circle area.

Reads static_indicator_volume.csv (indicator,N,h,vol0 from leiaSetFields) and plots
the RELATIVE error of the initialised phase volume VOL_ALPHA_0 = sum(alpha_c V_c)
against the exact pseudo-2D cylinder volume V = pi r^2 * t_z. The relative error is
height-free (it equals the relative circle-area error pi r^2), so the pseudo-2D slab
height t_z divides out exactly. Both indicators use the same LLS plane, so their
initialised volumes coincide; the figure shows both and the measured convergence order.
"""
import csv, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 0.15          # circle radius (implicitSphere, cases/2Dvortex)
TZ = 0.1          # pseudo-2D slab thickness (z in [-0.05, 0.05])
A_EXACT = math.pi * R**2          # exact circle area
V_EXACT = A_EXACT * TZ            # exact cylinder volume

csv_path = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/tmaric/OpenFOAM/repos/leia/workflow/scripts/static_indicator_volume.csv"
import paths  # thematic docs layout
out = sys.argv[2] if len(sys.argv) > 2 else \
    os.path.join(paths.figs_dir("semi-lagrangian-level-set"),
                 "static_indicator_volume_convergence.png")

rows = list(csv.DictReader(open(csv_path)))
def series(ind):
    s = sorted([r for r in rows if r["indicator"] == ind], key=lambda r: float(r["h"]))
    hs = [float(r["h"]) for r in s]
    err = [abs(float(r["vol0"]) - V_EXACT) / V_EXACT for r in s]
    return hs, err

def order(hs, err):
    xs = [math.log(h) for h in hs]; ys = [math.log(e) for e in err]
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    return (n*sxy - sx*sy) / (n*sxx - sx*sx)

STYLE = {"geometric": ("#1e88e5", "o", "geometric  (polygon clip)"),
         "detrixheAslam": ("#e53935", "s", "detrixheAslam  (tet-fill)")}

fig, ax = plt.subplots(figsize=(7.6, 5.6))
ords = {}
for ind in ("geometric", "detrixheAslam"):
    hs, err = series(ind)
    if not hs:
        continue
    ords[ind] = order(hs, err)
    col, mk, lab = STYLE[ind]
    # small marker-size offset so the (coincident) series stay visible on top of each other
    ms = 10 if ind == "geometric" else 6
    ax.loglog(hs, err, marker=mk, color=col, lw=2, ms=ms, mfc="none" if ind=="geometric" else col,
              label=f"{lab}  (order {ords[ind]:.2f})")

# O(h^2) reference anchored to the coarsest geometric point
hs, err = series("geometric")
x0, y0 = max(hs), err[hs.index(max(hs))]
gx = sorted(hs); gy = [y0 * (x / x0)**2 for x in gx]
ax.loglog(gx, gy, "--", color="#777", lw=1.3, label=r"$O(h^2)$ guide")

ax.set_xlabel(r"mesh spacing  $h$", fontsize=12)
ax.set_ylabel(r"relative volume error  $|V_\alpha - \pi r^2 t_z|\,/\,\pi r^2 t_z$", fontsize=12)
ax.set_title("Phase-indicator initialisation converges to the exact area\n"
             r"circle $r=0.15$, pseudo-2D height $t_z=0.1$ divided out", fontsize=12.5, fontweight="bold")
ax.grid(True, which="both", ls=":", alpha=0.45)
ax.legend(fontsize=9.5, loc="lower right")
ax.annotate("geometric $\\equiv$ detrixheAslam\n(same LLS plane, identical to 8 digits)",
            (0.03, 0.06), xycoords="axes fraction", fontsize=9.5, color="#2e7d32", va="bottom")
fig.tight_layout()
fig.savefig(out, dpi=130)
print(f"exact area pi r^2 = {A_EXACT:.8f}, exact cyl vol = {V_EXACT:.8f}")
for ind, o in ords.items():
    print(f"  {ind:14s} measured order = {o:.3f}")
print("wrote", out)
