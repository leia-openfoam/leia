#!/usr/bin/env python3
"""Perturb the interface cells of an initial volume-fraction field by a fixed
amount (2026-10-30, the serial-versus-serial control of the decomposition
comparison of the parallel swept flux): every cell with 1e-8 < alpha < 1 - 1e-8
gets alpha + eps, the full and empty cells are left exact; the field is
rewritten with 17 significant digits. Runs in the case directory after the
field initialisation and before decomposePar (a study's setfields_args, e.g.
"; python3 ../../../../workflow/scripts/davof_perturb_alpha.py 0/alpha.water 1e-14").

    python3 davof_perturb_alpha.py <field file> <eps>
"""
import re
import sys


def main(argv):
    path, eps = argv[1], float(argv[2])
    text = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\(", text)
    if not m:
        sys.exit(f"davof_perturb_alpha: no nonuniform internalField in {path}")
    n = int(m.group(1))
    start = m.end()
    end = text.index(")", start)
    vals = text[start:end].split()
    if len(vals) != n:
        sys.exit(f"davof_perturb_alpha: {len(vals)} values, header says {n}")
    out, k = [], 0
    for v in vals:
        a = float(v)
        if 1e-8 < a < 1 - 1e-8:
            a += eps
            k += 1
        out.append(repr(a))
    text = text[:start] + "\n" + "\n".join(out) + "\n" + text[end:]
    open(path, "w").write(text)
    print(f"davof_perturb_alpha: {k} of {n} cells of {path} perturbed by {eps}")


if __name__ == "__main__":
    main(sys.argv)
