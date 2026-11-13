#!/usr/bin/env python3
"""The configs of the bias test of the normal-order study (leia davof-normal-order,
2026-11-13; pre-registration in docs/davof/REPORT-swept.md "The bias test"). Writes
config/davof/nor_bias{Ob,Rot}_N{32,64,128}.yaml and nor_balance_N64.yaml.

    python3 workflow/scripts/nor_make_bias_configs.py

Row davof-planeSurf-ringTilt (the ring area normal rotated by theta_b = c h / R about
e_z x n, ringTiltNormal.H): TPF_TILT_MODE biased (the same sense in every cell) and
alternating (the sign by the parity of the cell's index sum), TPF_TILT_C in {0.125, 0.25,
0.5} (theta_b at N 32 of 0.026, 0.052, 0.104 against the ring normal's own L1 angle error
of 0.041 to 0.058), on the oblique translation (U = (0.3 0.5 0.8124), over 0.3) and the
rotation, N 32 / 64 / 128 on 8 ranks, OpenFOAM Co 0.2 and 1, with the position diagnostic
(positionDiag) and, on the translation, the per-step flux error. nor_balance_N64: the
oblique translation at N 64, OpenFOAM Co 0.2, the exact, ring, biased and alternating
rows (c = 0.25) for the step-resolved balance.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nor_make_configs as m

CS = ["0.125", "0.25", "0.5"]


def main():
    base = dict(m.SWEPT)
    base["TPF_POSITION_DIAG"] = ['"true"']
    for test, case, extra in (("Ob", "sphereTranslation3D", m.TRANS_OB), ("Rot", "sphereRotation3D", m.ROT)):
        for N, rt, mem in ((32, 480, 4000), (64, 1200, 4000), (128, 1800 if test == "Ob" else 2880, 8000)):
            axes = {"N_CELLS": [str(N)], "TPF_ROW": ["davof-planeSurf-ringTilt"],
                    "TPF_TILT_C": CS, "TPF_TILT_MODE": ["biased", "alternating"]}
            axes.update(extra)
            axes.update(base)
            if test == "Ob":
                axes["END_TIME"] = ['"0.3"']
                for k in ("TPF_WRITE_INTERVAL", "TPF_AFE_WRITE_INTERVAL", "TPF_VFE_WRITE_INTERVAL"):
                    axes[k] = ['"0.075"']
            m.write("nor_bias%s_N%d" % (test, N), case, 8, axes,
                    "The bias test, %s at N %d on 8 ranks: the ring normal rotated by a prescribed\n"
                    "# tilt theta_b = c h / R (c = 0.125, 0.25, 0.5; ringTiltNormal.H), the same sense in every\n"
                    "# cell (biased) or alternating by cell parity, OpenFOAM Co 0.2 and 1."
                    % ("the oblique translation" if test == "Ob" else "the rotation", N),
                    rt, mem)
    # the step-resolved balance: oblique translation N 64, OF Co 0.2
    axes = {"N_CELLS": ["64"], "TPF_ROW": ["davof-planeSurf-exact", "davof-planeSurf-ringArea", "davof-planeSurf-ringTilt"],
            "TPF_TILT_C": ["0.25"], "TPF_TILT_MODE": ["biased", "alternating"]}
    axes.update(m.TRANS_OB)
    axes["TARGET_ADVECTIVE_CO"] = ["0.124039"]
    axes.update(base)
    axes["END_TIME"] = ['"0.3"']
    for k in ("TPF_WRITE_INTERVAL", "TPF_AFE_WRITE_INTERVAL", "TPF_VFE_WRITE_INTERVAL"):
        axes[k] = ['"0.075"']
    m.write("nor_balance_N64", "sphereTranslation3D", 8, axes,
            "The step-resolved balance (oblique translation, N 64, OpenFOAM Co 0.2): per step the\n"
            "# flux error and the position error before and after CCS for the exact, ring, biased and\n"
            "# alternating (c = 0.25) normals.", 1200, 4000)


if __name__ == "__main__":
    main()
