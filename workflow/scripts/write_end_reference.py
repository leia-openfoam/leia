#!/usr/bin/env python3
"""Write the exact end-time reference fields psiEnd and alphaEnd of a one-way translation.

Usage, in the case directory, after leiaSetFields and before decomposePar:
    python3 workflow/scripts/write_end_reference.py [--case DIR] [leiaSetFields arguments ...]

The kinematic solvers score the fields at T against psiEnd/alphaEnd when the case has them
(errorCalculation.H, READ_IF_PRESENT) and against the INITIAL fields otherwise, which is right
only for a flow that returns. This script reads `endReference` from system/fvSolution:
  none (or absent)  do nothing (exit 0);
  translate         shift the implicitSurface centre by velocityModel/velocity * endTime, run
                    leiaSetFields once more with the same arguments (the same mesh, profile and
                    phase indicator as the run) on the pristine inputs of 0.org, rename its
                    psi/alpha to psiEnd/alphaEnd, and restore 0/, system/fvSolution and
                    leiaSetFields.csv.
`translate` requires `velocityModel { type translation; oscillation off; }` and an
implicitSphere surface; anything else is refused (exit 2), so a reversed or unsupported flow
never gets a wrong reference. Written 2026-09-29 for cases/2Dtranslation (STATUS 11.19).
Needs the OpenFOAM environment (foamDictionary, leiaSetFields). Standard library only.
"""
import os
import re
import shutil
import subprocess
import sys


def fd(entry, path):
    r = subprocess.run(["foamDictionary", "-entry", entry, "-value", path],
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def vec(s):
    return [float(x) for x in s.strip().strip("()").split()]


def rename_field(src, dst, old_name, new_name):
    text = open(src, encoding="utf-8").read()
    text, n = re.subn(r"^(\s*object\s+)" + re.escape(old_name) + r"\s*;", r"\1" + new_name + ";",
                      text, count=1, flags=re.M)
    if n != 1:
        raise SystemExit(f"write_end_reference: no 'object {old_name};' header in {src}")
    open(dst, "w", encoding="utf-8").write(text)
    os.remove(src)


def main():
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "--case":
        os.chdir(args[1]); args = args[2:]
    fvsol, ctrl = "system/fvSolution", "system/controlDict"
    mode = fd("endReference", fvsol)
    if mode in (None, "", "none"):
        print("write_end_reference: endReference none, nothing written")
        return 0
    if mode != "translate":
        print(f"write_end_reference: unknown endReference '{mode}' (none | translate)")
        return 2
    vtype = fd("velocityModel/type", fvsol)
    osc = (fd("velocityModel/oscillation", fvsol) or "on").lower()
    stype = fd("levelSet/implicitSurface/type", fvsol)
    if vtype != "translation" or osc not in ("off", "false", "no", "0") or stype != "implicitSphere":
        print(f"write_end_reference: refused (velocityModel {vtype}, oscillation {osc}, "
              f"implicitSurface {stype}); translate needs translation, oscillation off, implicitSphere")
        return 2
    U = vec(fd("velocityModel/velocity", fvsol))
    T = float(fd("endTime", ctrl))
    c = vec(fd("levelSet/implicitSurface/center", fvsol))
    cT = [c[i] + U[i] * T for i in range(3)]
    alpha = "alpha"
    if "-alphaName" in args:
        alpha = args[args.index("-alphaName") + 1]
    # The second leiaSetFields call must see EXACTLY the inputs of the first one (0.org), only the
    # centre shifted: started from the first call's output, the indicator treats the new circle's
    # cells as bulk and writes 0/1 values there (MEASURED 2026-09-29: 1.6 % volume, 0.46 in a cell).
    # So 0/ is snapshot, rebuilt from 0.org for the call, and restored afterwards.
    snap = "0.endRefSave"
    if os.path.exists(snap):
        shutil.rmtree(snap)
    shutil.copytree("0", snap)
    saved = {p: p + ".endRefSave" for p in (fvsol, "leiaSetFields.csv") if os.path.isfile(p)}
    for p, s in saved.items():
        shutil.copy2(p, s)
    tmp_end = {}
    try:
        for f in ("psi", alpha):
            src = os.path.join("0.org", f)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join("0", f))
        r = subprocess.run(["foamDictionary", "-entry", "levelSet/implicitSurface/center", "-set",
                            "(%.17g %.17g %.17g)" % tuple(cT), fvsol], capture_output=True, text=True)
        if r.returncode != 0:
            print("write_end_reference: foamDictionary -set failed\n" + r.stderr)
            return 3
        with open("log.endReference", "w") as log:
            r = subprocess.run(["leiaSetFields"] + args, stdout=log, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print("write_end_reference: leiaSetFields failed, see log.endReference")
            return 3
        for f, end in (("psi", "psiEnd"), (alpha, "alphaEnd")):
            dst = end + ".endRefNew"
            rename_field(os.path.join("0", f), dst, f, end)
            tmp_end[end] = dst
    finally:
        shutil.rmtree("0")
        shutil.move(snap, "0")
        for p, s in saved.items():
            shutil.move(s, p)
    for end, dst in tmp_end.items():
        shutil.move(dst, os.path.join("0", end))
    print("write_end_reference: centre (%s) -> (%s) at T = %g, U = (%s); wrote 0/psiEnd, 0/alphaEnd"
          % (" ".join("%g" % x for x in c), " ".join("%g" % x for x in cT), T,
             " ".join("%g" % x for x in U)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
