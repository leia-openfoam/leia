#!/usr/bin/env python3
"""check_refined_band -- re-run the band check of leia_refine on an existing case.

    python3 check_refined_band.py [case] --mode hex|poly [--band-cells 6]
                                  [--alpha-name alpha.water] [--allow-pin-mismatch]

Needs 0/psi and 0/<alpha> as written by leiaSetFields on the FINAL mesh (and, for
hex, constant/polyMesh/cellLevel). Runs `postProcess -func writeCellVolumes` when
0/V is absent and removes it again. Prints the verdict, rewrites refinedBand.csv,
exits 2 on FAIL.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case", nargs="?", default=".")
    ap.add_argument("--mode", choices=("hex", "poly"), default=None,
                    help="default: hex when constant/polyMesh/cellLevel exists, else poly")
    ap.add_argument("--band-cells", type=int, default=None)
    ap.add_argument("--alpha-name", default="alpha.water")
    ap.add_argument("--allow-pin-mismatch", action="store_true")
    a = ap.parse_args(argv)
    toks = lr.tokens(a.case)
    mode = a.mode or ("hex" if os.path.isfile(os.path.join(a.case, "constant", "polyMesh", "cellLevel"))
                      else "poly")
    band_cells = a.band_cells if a.band_cells is not None else lr.tok_int(toks, "REFINE_BAND_CELLS", 6)
    had_v = os.path.isfile(os.path.join(a.case, "0", "V"))
    lr.ensure_cell_volumes(a.case)
    res = lr.band_check(a.case, mode, band_cells, a.alpha_name, toks,
                        allow_pin_mismatch=a.allow_pin_mismatch)
    res.update(lr.checkmesh_stats(os.path.join(a.case, "log.checkMesh")))
    lr.write_csv(os.path.join(a.case, "refinedBand.csv"), [res])
    lr.report(res)
    if not had_v:
        os.remove(os.path.join(a.case, "0", "V"))
    return 0 if res["status"] != "FAIL" else 2


if __name__ == "__main__":
    sys.exit(main())
