#!/usr/bin/env python3
"""The spec files of the decisions 1-3 test matrix (2026-11-07), written on the cluster
from the prepared studies (the swept agent's rendered cases are read only, copied by the
runner). Groups: G1 identity on the eight-case set (20 steps), G2 Enright N 32 Co 0.5
(serial, T = 3), G3 Enright N 64 Co 0.5 (8 ranks, T = 3), G4 wall clock rotation N 64 OF Co 1
(serial and 8 ranks, three repetitions), G5 parallel identity rotation N 32 np 1/2/4/8
(20 steps). Lines: name|source|np|steps|reps|keys.
    python3 enright_specs.py <out dir>
"""
import glob, json, os, sys

W = os.environ.get("ENRIGHT_W", "/work/scratch/tm83tomy/davof/agents/enright")
SW = "/work/scratch/tm83tomy/davof/agents/swept/leia/studies/davof"
out = sys.argv[1]
os.makedirs(out, exist_ok=True)


def cases(study):
    cs = []
    for c in sorted(glob.glob(study + "/*_0*")):
        if "." in os.path.basename(c) or not os.path.isfile(c + "/case_params.json"): continue
        cs.append((c, json.load(open(c + "/case_params.json"))["tokens"]))
    return cs


def pick(study, **cond):
    res = []
    for c, t in cases(study):
        ok = True
        for k, v in cond.items():
            if k == "co_max":
                continue
            if str(t.get(k)) != str(v): ok = False
        if ok: res.append((float(t.get("TARGET_ADVECTIVE_CO", 0)), c))
    res.sort()
    if not res: sys.exit("no case in %s for %s" % (study, cond))
    return res[-1][1] if cond.get("co_max", True) else res[0][1]


# the eight cases (name, source, np)
eight = [
    ("rotN32", pick(SW + "/sphereRotation3D_hexCmpRotSweptN32", TPF_ROW="davof-planeSurf-ringArea"), 1),
    ("rotN64", pick(SW + "/sphereRotation3D_hexCmpRotSweptN64", TPF_ROW="davof-planeSurf-ringArea"), 8),
    ("diagN32", pick(SW + "/sphereRotation3D_hexCmpLineDiagSweptN32", TPF_ROW="davof-planeSurf-ringArea"), 1),
    ("planeN40exact", pick(SW + "/planeTranslation3D_hexSweptWispM9Ccs2", N_CELLS="40", TPF_ROW="davof-planeSurf-exact", TPF_VELOCITY_ARM="normal"), 1),
    ("planeN20ringNormal", pick(SW + "/planeTranslation3D_hexSweptWispM9Ccs2", N_CELLS="20", TPF_ROW="davof-planeSurf-ringArea", TPF_VELOCITY_ARM="normal"), 1),
    ("planeN20ringOblique", pick(SW + "/planeTranslation3D_hexSweptWispM9Ccs2", N_CELLS="20", TPF_ROW="davof-planeSurf-ringArea", TPF_VELOCITY_ARM="oblique"), 1),
    ("planeN40ringNormal", pick(SW + "/planeTranslation3D_hexSweptWispM9Ccs2", N_CELLS="40", TPF_ROW="davof-planeSurf-ringArea", TPF_VELOCITY_ARM="normal"), 1),
    ("planeN40ringOblique", pick(SW + "/planeTranslation3D_hexSweptWispM9Ccs2", N_CELLS="40", TPF_ROW="davof-planeSurf-ringArea", TPF_VELOCITY_ARM="oblique"), 1),
]
CFG = [
    ("closure", "sweptPieceVolume=closure"),
    ("cTT", "sweptPieceVolume=closureThenTets;sweptPieceValidate=true"),
    ("sT", "sweptPieceVolume=signedTets;sweptPieceValidate=true"),
    ("cTT_scale", "sweptPieceVolume=closureThenTets;sweptBoundedFaceFraction=scale"),
    ("cTT_recompute", "sweptPieceVolume=closureThenTets;sweptBoundedFaceFraction=recompute"),
    ("cTT_scale_nAA", "sweptPieceVolume=closureThenTets;sweptBoundedFaceFraction=scale;interfaceCellRule=normalAndAlpha"),
    ("cTT_scale_nAAR", "sweptPieceVolume=closureThenTets;sweptBoundedFaceFraction=scale;interfaceCellRule=normalAndAlphaReset"),
]
with open(out + "/G1_identity.spec", "w") as fh:
    fh.write("# G1: the eight-case identity, 20 steps, writePrecision 17\n")
    for name, src, np_ in eight:
        for cname, keys in CFG:
            fh.write("%s__%s|%s|%d|20|1|%s\n" % (name, cname, src, np_, keys))
enr32 = pick(W + "/leia/studies/davof/sphereDeformation3D_hexSweptN32", TPF_ROW="davof-planeSurf-ringArea", co_max=False)
enr64 = pick(W + "/leia/studies/davof/sphereDeformation3D_hexSweptN64", TPF_ROW="davof-planeSurf-ringArea", co_max=False)
for tag, src, np_ in (("G2_enright32", enr32, 1), ("G3_enright64", enr64, 8)):
    with open(out + "/%s.spec" % tag, "w") as fh:
        fh.write("# %s: the deformation at Co 0.5, T = 3, the fragment log on\n" % tag)
        for cname, keys in CFG:
            fh.write("%s|%s|%d|end:3|1|%s;sweptFragmentLog=true\n" % (cname, src, np_, keys))
rot64 = pick(SW + "/sphereRotation3D_hexCmpRotSweptN64", TPF_ROW="davof-planeSurf-ringArea")
with open(out + "/G4_timing.spec", "w") as fh:
    fh.write("# G4: the rotation N 64 at OF Co 1, one turn, serial and 8 ranks, three repetitions, the order rotated\n")
    for r in range(3):
        order = [CFG[0], CFG[1], CFG[2]]
        order = order[r:] + order[:r]
        for cname, keys in order:
            for np_ in (1, 8):
                fh.write("%s_np%d_r%d|%s|%d|end:1|1|%s\n" % (cname, np_, r + 1, rot64, np_, keys.replace(";sweptPieceValidate=true", "")))
with open(out + "/G5_parallel.spec", "w") as fh:
    fh.write("# G5: the rotation N 32, np 1/2/4/8, 20 steps\n")
    for k in (1, 2, 4, 8):
        src = pick(SW + "/sphereRotation3D_hexSweptParN32np%d" % k, TPF_ROW="davof-planeSurf-ringArea")
        for cname, keys in (CFG[1], CFG[2], CFG[3], CFG[5], CFG[6]):
            fh.write("%s_np%d|%s|%d|20|1|%s\n" % (cname, k, src, k, keys.replace(";sweptPieceValidate=true", "")))
# the second set (2026-11-07, after G2 showed signedTets alone bounded to 1e-7 on the
# deformation at N 32): the bounded-fraction and rule variants on the signedTets base
CFG2 = [
    ("sT_scale", "sweptPieceVolume=signedTets;sweptBoundedFaceFraction=scale"),
    ("sT_recompute", "sweptPieceVolume=signedTets;sweptBoundedFaceFraction=recompute"),
    ("sT_scale_nAA", "sweptPieceVolume=signedTets;sweptBoundedFaceFraction=scale;interfaceCellRule=normalAndAlpha"),
    ("sT_scale_nAAR", "sweptPieceVolume=signedTets;sweptBoundedFaceFraction=scale;interfaceCellRule=normalAndAlphaReset"),
]
with open(out + "/G1b_identity.spec", "w") as fh:
    fh.write("# G1b: the eight-case identity of the signedTets-based settings, 20 steps\n")
    for name, src, np_ in eight:
        for cname, keys in CFG2:
            fh.write("%s__%s|%s|%d|20|1|%s\n" % (name, cname, src, np_, keys))
with open(out + "/G2b_enright32.spec", "w") as fh:
    fh.write("# G2b: the deformation N 32 Co 0.5, the signedTets-based settings\n")
    for cname, keys in CFG2:
        fh.write("%s|%s|1|end:3|1|%s;sweptFragmentLog=true\n" % (cname, enr32, keys))
with open(out + "/G3a_enright64.spec", "w") as fh:
    fh.write("# G3a: the deformation N 64 Co 0.5 on 8 ranks: closure, cTT, sT, cTT_scale, cTT_scale_nAA\n")
    for cname, keys in (CFG[0], CFG[1], CFG[2], CFG[3], CFG[5]):
        fh.write("%s|%s|8|end:3|1|%s;sweptFragmentLog=true\n" % (cname, enr64, keys))
with open(out + "/G3b_enright64.spec", "w") as fh:
    fh.write("# G3b: the deformation N 64 Co 0.5 on 8 ranks: the signedTets-based settings\n")
    for cname, keys in CFG2:
        fh.write("%s|%s|8|end:3|1|%s;sweptFragmentLog=true\n" % (cname, enr64, keys))
with open(out + "/G5b_parallel.spec", "w") as fh:
    fh.write("# G5b: the rotation N 32, np 1/2/4/8, 20 steps, the signedTets-based settings\n")
    for k in (1, 2, 4, 8):
        src = pick(SW + "/sphereRotation3D_hexSweptParN32np%d" % k, TPF_ROW="davof-planeSurf-ringArea")
        for cname, keys in (CFG2[0], CFG2[2], CFG2[3]):
            fh.write("%s_np%d|%s|%d|20|1|%s\n" % (cname, k, src, k, keys))
for f in sorted(glob.glob(out + "/*.spec")):
    print(f, sum(1 for l in open(f) if l.strip() and not l.startswith("#")), "runs")
