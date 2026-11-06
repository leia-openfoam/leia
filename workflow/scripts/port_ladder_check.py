#!/usr/bin/env python3
"""The ladders of the port acceptance (docs/davof/PORT-DESIGN.md 4.3, 2026-11-06): the
rotation ladder (sphereRotation3D_hexRotFinal, N 32/64 serial and N 128 on 64 ranks) and
the sphere translation ladder (sphereTranslation3D_hexSweptLadderWispM9, N 32/64, and the
N 128 rungs), each run on the fork route (port_ladder_<name>_fork) and on the leia route
(port_ladder_<name>_leia) with the reference setting and the four round-4 keys:

  (a) the identity of the two routes per rung (port_identity_check.py: the written fields
      of the last time at writePrecision 12 of those cases, the function-object files, the
      converter CSV, the log numbers);
  (b) the ladders scored as the preprint scores them (davof_swept_ladders.py: the end
      values and run maxima of L1(alphaf), L1(n), L1(x), Eshape, the boundedness, the
      orders over N), both routes side by side.

Usage: port_ladder_check.py <studies_dir>
"""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sd = sys.argv[1]
LADDERS = {"rotation": ("port_ladder_rotation", "port_ladder_rotationN128"),
           "sphere": ("port_ladder_sphere", "port_ladder_sphereN128")}
allok = True
print("(a) the identity of the routes per study")
for name, studies in LADDERS.items():
    for st in studies:
        f, l = os.path.join(sd, st + "_fork"), os.path.join(sd, st + "_leia")
        if not (os.path.isdir(f) and os.path.isdir(l)):
            print(f"  {st}: missing"); allok = False; continue
        p = subprocess.run([sys.executable, os.path.join(HERE, "port_identity_check.py"), f, l], capture_output=True, text=True)
        lines = p.stdout.strip().split("\n")
        print(f"  {st}: {lines[-1]}")
        for x in lines[:-1]:
            if x.startswith(("IDENTICAL", "DIFFERENT")) or ("differs" in x) or ("first difference" in x) or x.startswith("    "):
                print("     " + x[:220])
        allok &= lines[-1] == "ALL IDENTICAL"
print("(b) the ladders scored (davof_swept_ladders.py), both routes")
for name, studies in LADDERS.items():
    for route in ("fork", "leia"):
        args = [sys.executable, os.path.join(HERE, "davof_swept_ladders.py")]
        for st in studies:
            d = os.path.join(sd, st + "_" + route)
            if os.path.isdir(d): args += ["--group", f"{name}_{route}", d]
        if len(args) > 2:
            p = subprocess.run(args, capture_output=True, text=True)
            print(f"--- {name}, {route} route")
            print(p.stdout.strip()[:6000])
            if p.returncode: print("  (scorer rc", p.returncode, p.stderr.strip()[-500:], ")")
print("LADDERS IDENTICAL" if allok else "LADDERS NOT IDENTICAL (see (a))")
sys.exit(0 if allok else 1)
