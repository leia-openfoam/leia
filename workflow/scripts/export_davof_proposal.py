#!/usr/bin/env python3
"""Copy the DAVOF proposal inputs out of the davof theme data.

    python3 workflow/scripts/export_davof_proposal.py --dest <proposal>/figures/davof

Copies tables/davof_normal_proposal_*.tex and figures/davof_plic_*.pdf (plus
the .png previews) from paths.tables_dir('davof') / paths.figs_dir('davof') to
--dest, creating it. The proposal \\input{}s and \\includegraphics{}es these
files by name, so a rerun of the studies followed by this export updates the
proposal's preliminary-results section without editing it.
"""
import argparse
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

PATTERNS = (("tables", "davof_normal_proposal_*.tex"),
            ("tables", "davof_regen_proposal_*.tex"),
            ("figures", "davof_plic_*.pdf"),
            ("figures", "davof_plic_*.png"))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", required=True, help="target directory (created)")
    ap.add_argument("--theme", default="davof")
    a = ap.parse_args(argv)
    srcs = {"tables": paths.tables_dir(a.theme), "figures": paths.figs_dir(a.theme)}
    os.makedirs(a.dest, exist_ok=True)
    n = 0
    for kind, pat in PATTERNS:
        for f in sorted(glob.glob(os.path.join(srcs[kind], pat))):
            shutil.copy2(f, os.path.join(a.dest, os.path.basename(f)))
            print(f"[export] {os.path.basename(f)} -> {a.dest}")
            n += 1
    if n == 0:
        print("[export] nothing found; run the davof studies first", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
