#!/usr/bin/env python3
"""Gr1 of the row pointIso-davof: L1(alpha_f) against the exact plane (areaFractionError.dat), per case.
Reference of the contact line study: M. Fricke, T. Maric, D. Bothe, J. Comput. Phys. 407 (2020) 109221,
doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.   python3 gr1_an.py <leia dir>"""
import glob, json, os, sys
R = sys.argv[1] if len(sys.argv) > 1 else "."
for arm in ("normal", "oblique"):
    for c in sorted(glob.glob(os.path.join(R, "studies/davof/pi_gr1_plane_%s/planeTranslation3D_0*" % arm))):
        tk = json.load(open(c + "/case_params.json"))["tokens"]
        fs = sorted(glob.glob(c + "/postProcessing/areaFractionError/*/areaFractionError.dat"))
        rows = []
        if fs:
            for l in open(fs[-1]):
                if l.startswith("#") or not l.strip():
                    continue
                rows.append([float(x) for x in l.split()])
        log = open(c + "/log.davofAdvectorFoam", errors="replace").read() if os.path.isfile(c + "/log.davofAdvectorFoam") else ""
        done = "\nEnd" in log and "[stack trace]" not in log
        withif = [r for r in rows if r[5] > 0]
        l1 = [r[1] for r in withif]
        print("%s N %s target Co %s: steps with interface cells %d of %d rows, L1(alpha_f): max %.3e, median %.3e, last %.3e; status %s" % (
            arm, tk["N_CELLS"], tk["TARGET_ADVECTIVE_CO"], len(withif), len(rows), max(l1) if l1 else float("nan"),
            sorted(l1)[len(l1)//2] if l1 else float("nan"), l1[-1] if l1 else float("nan"), "completed" if done else "NOT completed"))
