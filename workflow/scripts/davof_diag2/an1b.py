#!/usr/bin/env python3
"""Standing KISS results (read only): the pieces (recommended) cells at Co 0.2 and 0.5."""
import csv, glob, os, sys
sys.path.insert(0, "/work/scratch/tm83tomy/davof/agents/diag2")
import importlib.util
spec = importlib.util.spec_from_file_location("an1", "/work/scratch/tm83tomy/davof/agents/diag2/an1.py")
an1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(an1)

root = "/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof"
print("| family N | OF Co target | np | steps run | status | abort step (t) | min alpha before / after the conservative bounding (run) [step] | max alpha - 1 before / after the bounding (run) [step] | E_shape | wall (s) | keys (grazing / ccsEnd / ccsTol / tol) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for d in sorted(glob.glob(root + "/standing_*_variants")):
    f = glob.glob(d + "/*_database.csv")
    if not f:
        continue
    for r in csv.DictReader(open(f[0])):
        if r.get("TPF_ENRIGHT_SETTING") != "recommended" or r.get("TARGET_ADVECTIVE_CO") not in ("0.2", "0.5") or r.get("TPF_REPETITION") != "1":
            continue
        cd = "%s/%s_%s" % (d, os.path.basename(d).replace("standing_", "").split("_")[0] and r["case"], r["index"]) if False else None
        # the case directory name: <case>_<index>, case from the database row
        cdir = "%s/%s_%s" % (d, r["case"], r["index"]) if os.path.isdir("%s/%s_%s" % (d, r["case"], r["index"])) else None
        if cdir is None:
            cand = glob.glob(d + "/*_" + r["index"])
            cdir = cand[0] if cand else None
        if cdir is None or not os.path.exists(cdir + "/log.davofAdvectorFoam"):
            print("| %s | %s | %s | no log | | | | | | | |" % (os.path.basename(d), r["TARGET_ADVECTIVE_CO"], r["np"]))
            continue
        recs, abort = an1.parse(cdir + "/log.davofAdvectorFoam")
        mn = an1.ext(recs, "bmin", min); mx = an1.ext(recs, "bmax", max)
        mna = an1.ext(recs, "amin", min); mxa = an1.ext(recs, "amax", max)
        f2 = lambda v: ("%s [%s]" % (an1.fmt(v[0]), v[1])) if v else "-"
        print("| %s | %s | %s | %d | %s | %s | %s / %s | %s / %s | %s | %s | %s/%s/%s/%s |" % (
            os.path.basename(d).replace("standing_", "").replace("_variants", ""), r["TARGET_ADVECTIVE_CO"], r["np"], len(recs),
            "aborted" if abort else (r.get("advectorVoF.status") or r.get("logState")),
            ("%d (%.6g)" % (abort[0], abort[1])) if abort else "-",
            f2(mn), f2(mna), f2(mx), f2(mxa),
            (r.get("advectorVoF.Eshape") or "-")[:10], r.get("wallClockSolve"),
            r.get("TPF_GRAZING_CUT"), r.get("TPF_CCS_END_INTERVAL"), r.get("TPF_CCS_TOLERANCE"), r.get("TPF_POSITIONING_TOL")))
