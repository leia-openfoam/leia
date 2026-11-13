#!/usr/bin/env python3
"""Tally the NO PLANE dumps of the round-3 plane cases: per case the dumped faces by the state of
their vertex ring (all 0, all 1, mixed 0/1 pure, any partial alpha) and by the sweep."""
import glob, json, os, re, sys
S = "/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof"
DUMP = re.compile(r"NO PLANE in the vertex ring of face (\d+) at \(([^)]*)\) phi (\S+) sweeps (\d) V_swept (\S+) centroid \(([^)]*)\) centroid cell (-?\d+) alpha (\S+); ring:(.*)")
RING = re.compile(r"\[(\d+) a (\S+?)( P)?\]")
print("study,case,arm,N,co,steps,noPlaneSum,noPlaneMaxStep,dumps,ringAll0,ringAll1,ringMixedPure,ringPartial,ringWithPlane,dumpSweeps,dumpZeroFlux,limFacesSum,limVolumeSum,offFacesSum")
for d in sorted(glob.glob(os.path.join(S, "kiss_plane_*_np[18]"))):
    for c in sorted(glob.glob(os.path.join(d, "planeTranslation3D_*"))):
        lp = os.path.join(c, "log.davofAdvectorFoam")
        if not os.path.isfile(lp):
            continue
        t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        log = open(lp, errors="replace").read()
        nps = [int(x) for x in re.findall(r"none \(noPlaneFaces\) (\d+)", log)]
        lim = re.findall(r"limiter faces changed (\d+), cells (\d+), volume moved (\S+); faces of 0/1 cells with alphaf off 0/1 by more than 1e-12: (\d+)", log)
        k = dict(all0=0, all1=0, mixed=0, partial=0, plane=0, sweeps=0, zero=0)
        nd = 0
        for m in DUMP.finditer(log):
            nd += 1
            ring = [(float(a), bool(p)) for _, a, p in RING.findall(m.group(9))]
            al = [a for a, _ in ring]
            if any(p for _, p in ring): k["plane"] += 1
            if all(a == 0 for a in al): k["all0"] += 1
            elif all(a == 1 for a in al): k["all1"] += 1
            elif all(a in (0.0, 1.0) for a in al): k["mixed"] += 1
            else: k["partial"] += 1
            if m.group(4) == "1": k["sweeps"] += 1
            else: k["zero"] += 1
        print(",".join(str(x) for x in (os.path.basename(d), os.path.basename(c), t["TPF_VELOCITY_ARM"], t["N_CELLS"], t["TARGET_ADVECTIVE_CO"],
              len(re.findall(r"\nTime = ", log)), sum(nps), max(nps) if nps else "", nd, k["all0"], k["all1"], k["mixed"], k["partial"], k["plane"], k["sweeps"], k["zero"],
              sum(int(x[0]) for x in lim), sum(float(x[2]) for x in lim), sum(int(x[3]) for x in lim))))
