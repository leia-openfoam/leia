#!/usr/bin/env python3
"""The configs of the signed flux error study (leia davof-normal-order, 2026-11-14;
pre-registration in docs/davof/REPORT-swept.md "Signed flux error of the biased tilt").
Writes config/davof/nor_sgn_N{32,64,128}_{base,bia,alt}.yaml:

    python3 workflow/scripts/nor_make_sgn_configs.py

The oblique translation (U = (0.3 0.5 0.8124), over 0.3), OpenFOAM Co 0.2, 8 ranks, with
the position diagnostic and the flux error diagnostics (sweptFluxErrorDiag). base: the
exact and the ring area rows; bia: ringTilt biased; alt: ringTilt alternating. N 64:
biased c in {-0.5, -0.25, -0.125, 0.125, 0.25, 0.5}, alternating c in {0.125, 0.25, 0.5};
N 32 and 128: biased c = -0.25, 0.25, alternating c = 0.25.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nor_make_configs as m

BIA = {64: ["-0.5", "-0.25", "-0.125", "0.125", "0.25", "0.5"], 32: ["-0.25", "0.25"], 128: ["-0.25", "0.25"]}
ALT = {64: ["0.125", "0.25", "0.5"], 32: ["0.25"], 128: ["0.25"]}


def main():
    base = dict(m.SWEPT)
    base["TPF_POSITION_DIAG"] = ['"true"']
    base["TPF_FLUX_ERROR_DIAG"] = ['"true"']
    for N, rt, mem in ((32, 480, 4000), (64, 1200, 4000), (128, 2400, 8000)):
        for grp in ("base", "bia", "alt"):
            axes = {"N_CELLS": [str(N)]}
            if grp == "base":
                axes["TPF_ROW"] = ["davof-planeSurf-exact", "davof-planeSurf-ringArea"]
            else:
                axes["TPF_ROW"] = ["davof-planeSurf-ringTilt"]
                axes["TPF_TILT_C"] = BIA[N] if grp == "bia" else ALT[N]
                axes["TPF_TILT_MODE"] = ["biased" if grp == "bia" else "alternating"]
            axes.update(m.TRANS_OB)
            axes["TARGET_ADVECTIVE_CO"] = ["0.124039"]
            axes.update(base)
            axes["END_TIME"] = ['"0.3"']
            for k in ("TPF_WRITE_INTERVAL", "TPF_AFE_WRITE_INTERVAL", "TPF_VFE_WRITE_INTERVAL"):
                axes[k] = ['"0.075"']
            m.write("nor_sgn_N%d_%s" % (N, grp), "sphereTranslation3D", 8, axes,
                    "The signed flux error study (oblique translation, N %d, OpenFOAM Co 0.2, 8 ranks;\n"
                    "# group %s): the flux error diagnostics (sweptFluxErrorDiag) and the position\n"
                    "# diagnostic. Pre-registration: docs/davof/REPORT-swept.md \"Signed flux error of\n"
                    "# the biased tilt\"." % (N, grp), rt, mem)


if __name__ == "__main__":
    main()
