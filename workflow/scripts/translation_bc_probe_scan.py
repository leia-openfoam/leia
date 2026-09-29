#!/usr/bin/env python3
"""Read the runs of translation_bc_probe.sh (STATUS 11.19, 2026-09-29).

For each probe directory: E_GEOM_ALPHA_REL and E_VOL_ALPHA_REL at T, and at every write time
the largest |psi - psi_exact| over the cells next to the right (outflow) patch and the count
of false zero-set cells (psi < 0 farther than 5 h from the exact circle). The mesh is the
uniform N x N blockMesh of cases/2Dtranslation (cell i = ix + N iy).

Usage: translation_bc_probe_scan.py <parent dir> <N> <probe dir name> [<probe dir name> ...]
"""
import csv, os, re, sys
import numpy as np


def internal(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
    n = int(m.group(1)); s = m.end(); e = txt.index(")", s)
    return np.array(txt[s:e].split(), dtype=float)


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 1
    parent, N = argv[0], int(argv[1]); h = 1.0 / N
    ix = np.tile(np.arange(N), N); iy = np.repeat(np.arange(N), N)
    xc, yc = (ix + 0.5) * h, (iy + 0.5) * h
    for name in argv[2:]:
        d = os.path.join(parent, name)
        r = list(csv.DictReader(open(os.path.join(d, "leiaSemiLagrangeLevelSetFoam.csv"))))[-1]
        print(f"{name}: E_GEOM_ALPHA_REL(T) {float(r['E_GEOM_ALPHA_REL']):.4e}"
              f"  E_VOL_ALPHA_REL(T) {float(r['E_VOL_ALPHA_REL']):.4e}")
        times = sorted((x for x in os.listdir(d) if re.fullmatch(r"[0-9.]+", x)), key=float)[1:]
        for tn in times:
            t = float(tn); psi = internal(os.path.join(d, tn, "psi"))
            ex = np.sqrt((xc - 0.25 - t)**2 + (yc - 0.5)**2) - 0.15
            err = np.abs(psi - ex)
            print(f"   t={t:.3f}  outflow-edge max|psi - exact| {err[ix == N - 1].max():.2e}"
                  f"  false zero-set cells {int(((psi < 0) & (ex > 5 * h)).sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
