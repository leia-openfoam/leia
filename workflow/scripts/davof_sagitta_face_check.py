#!/usr/bin/env python3
"""Face-by-face check of the sagitta correction (study sphereTranslation3D_hexSagittaDiag,
2026-10-16): the plane row's per-face error e_plane = alphaf - alphaf_exact
(alphafError.<group> of the davof-exact case), the sagitta row's e_sag (the same field of the
davof-sagitta-exact case) and the correction dA/|S_f| the model applied (sagittaCorrection
.<group> of the sagitta case), all after ONE step from the identical state, so that
e_sag = e_plane + dA/|S_f| holds face by face if the stored value is what the model returns,
and dA/|S_f| = -e_plane up to O(h^2) if the correction is the right one.

    python3 davof_sagitta_face_check.py <plane case> <sagitta case> <time> [--group water]
"""
import argparse
import math
import os
import re
import sys


def read_scalar_list(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
    if not m:
        sys.exit(f"{path}: no nonuniform scalar internalField")
    n = int(m.group(1))
    body = txt[m.end():]
    end = body.index("\n)\n") if "\n)\n" in body else body.index(")\n;")
    vals = [float(v) for v in body[:end].split()]
    if len(vals) != n:
        sys.exit(f"{path}: {len(vals)} read, {n} declared")
    return vals


def stats(name, xs):
    n = len(xs)
    if n == 0:
        print(f"{name}: none")
        return
    mean = sum(xs)/n
    mabs = sum(abs(x) for x in xs)/n
    rms = math.sqrt(sum(x*x for x in xs)/n)
    print(f"{name}: n {n}, mean {mean:+.3e}, mean|.| {mabs:.3e}, rms {rms:.3e}, "
          f"min {min(xs):+.3e}, max {max(xs):+.3e}")


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("plane_case")
    ap.add_argument("sagitta_case")
    ap.add_argument("time")
    ap.add_argument("--group", default="water")
    a = ap.parse_args(argv)
    ep = read_scalar_list(os.path.join(a.plane_case, a.time, f"alphafError.{a.group}"))
    es = read_scalar_list(os.path.join(a.sagitta_case, a.time, f"alphafError.{a.group}"))
    dA = read_scalar_list(os.path.join(a.sagitta_case, a.time, f"sagittaCorrection.{a.group}"))
    if not (len(ep) == len(es) == len(dA)):
        sys.exit(f"sizes differ: {len(ep)} {len(es)} {len(dA)}")
    counted = [i for i in range(len(ep)) if ep[i] != 0 or es[i] != 0]
    corrected = [i for i in range(len(dA)) if dA[i] != 0]
    print(f"internal faces {len(ep)}; faces with an error entry {len(counted)}; faces corrected {len(corrected)}")
    stats("e_plane on counted faces", [ep[i] for i in counted])
    stats("e_sag   on counted faces", [es[i] for i in counted])
    stats("e_plane on corrected faces", [ep[i] for i in corrected])
    stats("e_sag   on corrected faces", [es[i] for i in corrected])
    stats("dA/|Sf| on corrected faces", [dA[i] for i in corrected])
    stats("e_sag - e_plane - dA (storage identity)", [es[i] - ep[i] - dA[i] for i in corrected])
    stats("e_plane + dA (what the correction leaves)", [ep[i] + dA[i] for i in corrected])
    # the correlation of dA with -e_plane
    xs = [dA[i] for i in corrected]
    ys = [-ep[i] for i in corrected]
    if xs:
        mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
        sxy = sum((x - mx)*(y - my) for x, y in zip(xs, ys))
        sxx = sum((x - mx)**2 for x in xs); syy = sum((y - my)**2 for y in ys)
        r = sxy/math.sqrt(sxx*syy) if sxx > 0 and syy > 0 else float("nan")
        slope = sxy/sxx if sxx > 0 else float("nan")
        print(f"correlation of dA with -e_plane: r = {r:.3f}, regression slope (-e_plane per dA) = {slope:.3f}")
    # corrected faces where the plane error is small but the correction is not, and vice versa
    big_wrong = sum(1 for i in corrected if abs(ep[i] + dA[i]) > 2*abs(ep[i]))
    print(f"corrected faces whose error GREW beyond twice the plane's: {big_wrong} of {len(corrected)}")
    uncorrected = [i for i in counted if dA[i] == 0]
    stats("e_plane on counted but uncorrected faces", [ep[i] for i in uncorrected])
    print("first 12 corrected faces: e_plane, dA, e_sag")
    for i in corrected[:12]:
        print(f"  {i:7d} {ep[i]:+.4e} {dA[i]:+.4e} {es[i]:+.4e}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
