#!/usr/bin/env python3
"""diag2 Part 1 analysis: one row per case log (study dir given on the command line as
studies/davof/<study>[/<case>]): status, abort step and time, run extremes of min(alpha) and
max(alpha) - 1 (after the whole step, the solver's line; and after the conservative bounding),
E_shape and wall clock from the study database."""
import csv, glob, json, os, re, sys

rx_time = re.compile(r"^Time = (\S+)")
rx_before = re.compile(r"davofSwept: before the bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)")
rx_after = re.compile(r"davofSwept: After conservative bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)")
rx_solver = re.compile(r"^t = (\S+),\s+sum\(alpha\*V\) = (\S+),\s+dev = (\S+)%,\s+1-max\(alpha1\) = (\S+) at cell (\d+),\s+min\(alpha1\) = (\S+) at cell (\d+)")
rx_abort = re.compile(r"ABORT at step (\d+), t = (\S+): (.*)")

def parse(path):
    step = 0
    recs = []
    cur = None
    abort = None
    for line in open(path, errors="replace"):
        m = rx_time.match(line)
        if m:
            step += 1
            cur = dict(step=step, t=float(m.group(1)), bmin=None, bmax=None, amin=None, amax=None, smin=None, smax=None)
            recs.append(cur)
            continue
        if cur is None:
            continue
        m = rx_before.search(line)
        if m:
            cur["bmin"], cur["bmax"] = float(m.group(1)), float(m.group(2)); continue
        m = rx_after.search(line)
        if m:
            cur["amin"], cur["amax"] = float(m.group(1)), float(m.group(2)); continue
        m = rx_solver.match(line)
        if m:
            cur["smax"] = -float(m.group(4)); cur["smin"] = float(m.group(6)); continue
        m = rx_abort.search(line)
        if m:
            abort = (int(m.group(1)), float(m.group(2)), m.group(3)[:90])
    return recs, abort

def ext(recs, key, fn):
    vals = [(r[key], r["step"]) for r in recs if r[key] is not None]
    if not vals:
        return None
    return fn(vals, key=lambda v: v[0])

def db_row(casedir):
    study = os.path.dirname(casedir)
    idx = os.path.basename(casedir).split("_")[-1]
    f = glob.glob(study + "/*_database.csv")
    if not f:
        return {}
    for r in csv.DictReader(open(f[0])):
        if r["index"] == idx:
            return r
    return {}

def fmt(v):
    return "-" if v is None else ("%.3e" % v)

if __name__ == "__main__":
    print("| case | keys (ccsEnd/ccsTol/grazing) | np | OF Co target | steps run | status | abort step (t) | min alpha before / after the conservative bounding (run) [step] | max alpha - 1 before / after the bounding (run) [step] | min alpha before / after at the last step | E_shape | wall (s) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for arg in sys.argv[1:]:
        dirs = sorted(glob.glob(arg + "/isoShear*_0*") + glob.glob(arg + "/isoDeformation*_0*") + glob.glob(arg + "/*_0*")) if os.path.isdir(arg) and not os.path.exists(arg + "/log.davofAdvectorFoam") else [arg]
        seen = set()
        for d in dirs:
            if d in seen or not os.path.exists(d + "/log.davofAdvectorFoam"):
                continue
            seen.add(d)
            recs, abort = parse(d + "/log.davofAdvectorFoam")
            r = db_row(d)
            keys = "%s/%s/%s" % (r.get("TPF_CCS_END_INTERVAL", "?"), r.get("TPF_CCS_TOLERANCE", "?"), r.get("TPF_GRAZING_CUT", "?"))
            mn = ext(recs, "bmin", min); mx = ext(recs, "bmax", max)
            mna = ext(recs, "amin", min); mxa = ext(recs, "amax", max)
            last = recs[-1] if recs else {}
            f2 = lambda v: ("%s [%s]" % (fmt(v[0]), v[1])) if v else "-"
            print("| %s | %s | %s | %s | %d | %s | %s | %s / %s | %s / %s | %s / %s | %s | %s |" % (
                os.path.relpath(d, "studies/davof") if os.path.isabs(d) is False else d.split("studies/davof/")[-1],
                keys, r.get("np", "?"), r.get("TARGET_ADVECTIVE_CO", "?"), len(recs),
                ("aborted" if abort else (r.get("advectorVoF.status") or r.get("logState") or "?")),
                ("%d (%.6g)" % (abort[0], abort[1])) if abort else "-",
                f2(mn), f2(mna), f2(mx), f2(mxa),
                fmt(last.get("bmin")), fmt(last.get("amin")),
                (r.get("advectorVoF.Eshape") or "-")[:10], r.get("wallClockSolve", "-")))
