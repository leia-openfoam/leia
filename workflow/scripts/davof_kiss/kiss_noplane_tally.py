#!/usr/bin/env python3
"""Tally the NO PLANE dumps of the one-plane plane cases (python3 kiss_noplane_tally.py <glob>):
per case the dumped faces by ring state and the alpha of the planeless centroid cell
(the dumps are the first five cut faces per rank and step)."""
import glob, json, os, re, sys
S = "/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof"
PAT = sys.argv[1] if len(sys.argv) > 1 else "kiss_plane_*_np[18]"
DUMP = re.compile(r"NO PLANE in the vertex ring of face (\d+) at \(([^)]*)\) phi (\S+) sweeps (\d) V_swept (\S+) centroid \(([^)]*)\) centroid cell (-?\d+) alpha (\S+); ring:(.*)")
RING = re.compile(r"\[(\d+) a (\S+?)( P)?\]")
print("study,case,arm,N,co,steps,noPlaneSum,noPlaneCutSum,dumps,cenAlphaPartialN,cenAlphaMin,cenAlphaMax,ringPartialMax,limFacesSum,limVolumeSum,offFacesSum")
for d in sorted(glob.glob(os.path.join(S, PAT))):
    for c in sorted(glob.glob(os.path.join(d, "planeTranslation3D_*"))):
        lp = os.path.join(c, "log.davofAdvectorFoam")
        if not os.path.isfile(lp):
            continue
        t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        log = open(lp, errors="replace").read()
        nps = [int(x) for x in re.findall(r"none \(noPlaneFaces\) (\d+)", log)]
        cut = [int(x) for x in re.findall(r"of them cut (\d+)", log)]
        lim = re.findall(r"limiter faces changed (\d+), cells (\d+), volume moved (\S+); faces of 0/1 cells with alphaf off 0/1 by more than 1e-12: (\d+)", log)
        cen = []; rmax = 0.0; nd = 0
        for m in DUMP.finditer(log):
            nd += 1
            a = float(m.group(8))
            if 1e-12 < a < 1 - 1e-12:
                cen.append(min(a, 1 - a))
            for _, ra, _p in RING.findall(m.group(9)):
                ra = float(ra)
                if 1e-12 < ra < 1 - 1e-12:
                    rmax = max(rmax, min(ra, 1 - ra))
        print(",".join(str(x) for x in (os.path.basename(d), os.path.basename(c), t["TPF_VELOCITY_ARM"], t["N_CELLS"], t["TARGET_ADVECTIVE_CO"],
              len(re.findall(r"\nTime = ", log)), sum(nps), sum(cut), nd, len(cen), "%.3g" % min(cen) if cen else "", "%.3g" % max(cen) if cen else "", "%.3g" % rmax,
              sum(int(x[0]) for x in lim), sum(float(x[2]) for x in lim), sum(int(x[3]) for x in lim))))
