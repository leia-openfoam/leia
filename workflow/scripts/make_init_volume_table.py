#!/usr/bin/env python3
"""Initialisation accuracy: does the phase indicator's volume converge to the
ANALYTIC volume of the initialised shape?

WHY THIS GATE EXISTS. Every volume error the solvers report is referenced to the
DISCRETE initial volume -- `dropletInitialPhaseVolume` in the two-phase solver,
`alphaV0 = gSum(alpha0*V)` in the kinematic one. That is the right reference for
measuring conservation DRIFT, but it means the initialisation error cancels out
of the reported metric by construction and can never appear there. So the
accuracy of alpha(psi) at t = 0 has to be measured separately, against the
closed-form volume of the shape that was initialised. `leiaSetFields` already
writes it (`leiaSetFields.csv`: H, VOL_ALPHA_0); this script compares it with the
analytic value over a resolution ladder and fits the order.

Exact volumes, taken from the case's own `levelSet/implicitSurface` dictionary:

  * one-cell-thick mesh (nz = 1, the repo's 2D convention): the body is a
    cylinder of the given cross-section over the mesh thickness t
        implicitSphere                     ->  pi R^2 t
        signedDistanceEllipse / ellipsoid  ->  pi a b t
  * genuinely 3D mesh
        implicitSphere                     ->  4/3 pi R^3
        ellipsoid (any flavour)            ->  4/3 pi a b c

Writes `<study>_init_volume.csv` into the study dir and the theme's tables dir,
plus a log-log figure with the fitted order.

Usage:
  make_init_volume_table.py <study_dir> [--theme <theme>]
"""

import argparse
import csv
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import paths  # noqa: E402  (repo helper: theme -> data dirs)
except Exception:                                    # pragma: no cover
    paths = None


def _dict_block(text, name):
    """Body of the `name { ... }` sub-dictionary (first, non-nested match)."""
    i = text.find(name)
    while i != -1:
        j = text.find("{", i)
        if j != -1:
            depth, k = 1, j + 1
            while k < len(text) and depth:
                if text[k] == "{":
                    depth += 1
                elif text[k] == "}":
                    depth -= 1
                k += 1
            return text[j + 1:k - 1]
        i = text.find(name, i + 1)
    return ""


def surface_of(case):
    """(type, (a, b, c)) of the initialised implicit surface."""
    fv = os.path.join(case, "system", "fvSolution")
    if not os.path.isfile(fv):
        return None, None
    blk = _dict_block(open(fv).read(), "implicitSurface")
    if not blk:
        return None, None
    m = re.search(r"\btype\s+(\w+)\s*;", blk)
    typ = m.group(1) if m else None
    m = re.search(r"\baxes\s*\(\s*([^)]*)\)", blk)
    if m:
        v = [float(x) for x in m.group(1).split()]
        return typ, (v[0], v[1], v[2] if len(v) > 2 else v[1])
    m = re.search(r"\bradius\s+([0-9.eE+-]+)\s*;", blk)
    if m:
        r = float(m.group(1))
        return typ, (r, r, r)
    return typ, None


def mesh_extent_z(case):
    """(thickness, nz) from blockMeshDict; (None, None) if unreadable."""
    bm = os.path.join(case, "system", "blockMeshDict")
    if not os.path.isfile(bm):
        return None, None
    s = open(bm).read()
    # ONLY the vertices block: a bare triple regex over the whole file also
    # matches `simpleGrading (1 1 1)` and the hex cell counts, which put z = 1
    # into the list and inflated the thickness by ~3200x (caught by validating
    # against a hand-computed ladder).
    vblk = ""
    mv = re.search(r"\bvertices\s*\(", s)
    if mv:
        depth, k = 1, mv.end()
        while k < len(s) and depth:
            if s[k] == "(":
                depth += 1
            elif s[k] == ")":
                depth -= 1
            k += 1
        vblk = s[mv.end():k - 1]
    zs = [float(m.group(1)) for m in
          re.finditer(r"\(\s*[-0-9.eE+]+\s+[-0-9.eE+]+\s+([-0-9.eE+]+)\s*\)", vblk)]
    thick = (max(zs) - min(zs)) if zs else None
    m = re.search(r"hex\s*\([^)]*\)\s*\(\s*\S+\s+\S+\s+(\S+?)\s*\)", s)
    nz = None
    if m:
        tok = m.group(1)
        if tok.startswith("$"):
            mv = re.search(re.escape(tok[1:]) + r"\s+([0-9]+)\s*;", s)
            nz = int(mv.group(1)) if mv else None
        else:
            try:
                nz = int(tok)
            except ValueError:
                nz = None
    return thick, nz


def exact_volume(typ, axes, thick, nz):
    if not typ or not axes:
        return None
    a, b, c = axes
    if nz == 1:                      # one-cell-thick mesh: prism over the section
        if thick is None:
            return None
        return math.pi * a * b * thick
    return 4.0 / 3.0 * math.pi * a * b * c


def collect(study_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(study_dir, "*", "leiaSetFields.csv"))):
        case = os.path.dirname(f)
        recs = [r for r in csv.DictReader(open(f)) if r.get("VOL_ALPHA_0")]
        if not recs:
            continue
        h = float(recs[0]["H"])
        vol = float(recs[0]["VOL_ALPHA_0"])
        typ, axes = surface_of(case)
        thick, nz = mesh_extent_z(case)
        vex = exact_volume(typ, axes, thick, nz)
        if not vex or vol <= 0:
            continue
        rows.append(dict(case=os.path.basename(case), surface=typ, h=h, nz=nz or 0,
                         volume=vol, volumeExact=vex,
                         relError=abs(vol - vex) / vex))
    rows.sort(key=lambda r: (r["surface"] or "", -r["h"]))
    return rows


def with_orders(rows):
    """Per-level order within each surface flavour (coarse -> fine)."""
    for surf in {r["surface"] for r in rows}:
        sub = sorted([r for r in rows if r["surface"] == surf],
                     key=lambda r: -r["h"])
        prev = None
        for r in sub:
            r["order"] = ""
            if prev and r["relError"] > 0 and prev["relError"] > 0:
                r["order"] = round(
                    math.log(prev["relError"] / r["relError"])
                    / math.log(prev["h"] / r["h"]), 3)
            prev = r
    return rows


def figure(rows, out_png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for surf in sorted({r["surface"] for r in rows}):
        sub = sorted([r for r in rows if r["surface"] == surf], key=lambda r: r["h"])
        ax.loglog([r["h"] for r in sub], [r["relError"] for r in sub],
                  "o-", label=surf)
    hs = [r["h"] for r in rows]
    if hs:
        h0, e0 = max(hs), max(r["relError"] for r in rows)
        ax.loglog([h0, h0 / 4], [e0, e0 / 16], "k--", lw=1, label="second order")
    ax.set_xlabel("cell size $h$ [m]")
    ax.set_ylabel(r"$|V_\alpha - V_{\rm exact}| / V_{\rm exact}$")
    ax.set_title("Phase-indicator initialisation: volume vs analytic")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default=None)
    a = ap.parse_args()

    rows = with_orders(collect(a.study_dir))
    if not rows:
        print("[init-volume] no leiaSetFields.csv with a resolvable surface; skipped")
        return 0

    study = os.path.basename(os.path.normpath(a.study_dir))
    cols = ["case", "surface", "h", "nz", "volume", "volumeExact", "relError", "order"]
    local = os.path.join(a.study_dir, f"{study}_init_volume.csv")
    with open(local, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"[init-volume] wrote {local}")

    for r in rows:
        print(f"  {r['surface']:<24} h={r['h']:.5g}  rel={r['relError']:.4e}"
              f"  order={r['order']}")

    if a.theme and paths:
        try:
            tdir = paths.tables_dir(a.theme)
            fdir = paths.figs_dir(a.theme)
            os.makedirs(tdir, exist_ok=True)
            os.makedirs(fdir, exist_ok=True)
            import shutil
            shutil.copy(local, os.path.join(tdir, f"{study}_init_volume.csv"))
            figure(rows, os.path.join(fdir, f"{study}_init_volume.png"))
            print(f"[init-volume] curated into {tdir}")
        except Exception as e:                        # pragma: no cover
            print(f"[init-volume] theme export skipped: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
