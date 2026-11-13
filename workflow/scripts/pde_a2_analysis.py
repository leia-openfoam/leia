#!/usr/bin/env python3
"""pde_a2_analysis.py: the read-out of step A0.2 (the band-edge repair) of the PDE level-set plan.

    pde_a2_analysis.py --old <clone> --new <clone> [--celldist <file>] [--out <dir>]

1. Inertness, static: every legacy arm of aStatic2Dedge (new clone) against the same
   (VELOCITY_EXTENSION, FADE_MODE, N_CELLS) arm of aStatic2Dfade (old clone), tolerance 0, every
   common column of both CSVs of leiaTestVelocityExtension.
2. The fallback cells of closestPoint per arm (the last "fallback cells" count of the log; a
   diagnostic count, not a run state: the run state comes from foam_log_state.sh).
3. The seam check: psi at t = 0.02, a2Seam2Dvortexserial against a2Seam2Dvortex, per arm
   (compare_psi_fields.py), with the run state of both logs.
The advected identity (a0Identity2Dvortex) is read by pde_a0_analysis.py --old --new.
Writes <out>/a2_summary.txt (default <new>/studies/pdeLevelSet/a2). Stdlib only.
"""
import argparse, glob, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"


def tok(case):
    return {k: str(v) for k, v in json.load(open(os.path.join(case, "case_params.json")))["tokens"].items()}


def run(args):
    p = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    return p.returncode, p.stdout.strip().splitlines()


def state(log):
    p = subprocess.run([os.path.join(HERE, "foam_log_state.sh"), log], capture_output=True, text=True)
    return (p.stdout.strip().splitlines() or ["?"])[0]


def fallback(case):
    logs = glob.glob(os.path.join(case, "log.leia*"))
    n = None
    for lg in logs:
        for m in re.finditer(r"(\d+) foot-pointed band cells .*?, (\d+) fallback cells", open(lg, errors="ignore").read()):
            n = (int(m.group(1)), int(m.group(2)))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
    ap.add_argument("--celldist"); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "a2")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))

    say("1. Inertness (legacy arms of aStatic2Dedge against aStatic2Dfade, tolerance 0):")
    old = {}
    for c in glob.glob(os.path.join(a.old, "studies/pdeLevelSet/aStatic2Dfade", "2Dvortex_[0-9]*")):
        t = tok(c); old[(t["VELOCITY_EXTENSION"], t["FADE_MODE"], t["N_CELLS"])] = c
    n = bad = 0
    for c in sorted(glob.glob(os.path.join(a.new, "studies/pdeLevelSet/aStatic2Dedge", "2Dvortex_[0-9]*"))):
        t = tok(c)
        if t.get("VE_FADE_PLACEMENT") != "legacy":
            continue
        key = (t["VELOCITY_EXTENSION"], t["FADE_MODE"], t["N_CELLS"])
        o = old.get(key)
        for f in ("leiaTestVelocityExtension.csv", "leiaTestVelocityExtensionExact.csv"):
            if not o or not os.path.isfile(os.path.join(c, f)):
                v, d = "MISSING", []
            else:
                rc, d = run([os.path.join(HERE, "compare_metrics_csv.py"), os.path.join(o, f), os.path.join(c, f),
                             "--tol", "0", "--skip", SKIP])
                v = "PASS" if rc == 0 else "DIFF"
            n += 1; bad += v != "PASS"
            say(f"  {'/'.join(key):32s} {f:38s} {v}  {' | '.join(d[-2:])[:140]}")
    say(f"  -> {n - bad} of {n} identical: {'PASS' if n and not bad else 'FAIL'}")

    say()
    say("2. closestPoint: foot-pointed band cells and fallback cells (last call of each run):")
    for c in sorted(glob.glob(os.path.join(a.new, "studies/pdeLevelSet/aStatic2Dedge", "2Dvortex_[0-9]*"))):
        t = tok(c)
        if t["VELOCITY_EXTENSION"] != "closestPoint":
            continue
        fb = fallback(c)
        frac = f"{fb[1] / max(fb[0] + fb[1], 1):.2%}" if fb else "?"
        say(f"  N {t['N_CELLS']:>4s} {t['FADE_MODE']:9s} {t['VE_FADE_PLACEMENT']:18s} {fb} fallback share {frac}")

    say()
    say("3. Seam: psi at t = 0.02, serial against np 4, relative to max|psi|:")
    ser = {}
    for c in glob.glob(os.path.join(a.new, "studies/pdeLevelSet/a2Seam2Dvortexserial", "2Dvortex_[0-9]*")):
        t = tok(c); ser[(t["FADE_MODE"], t["VE_FADE_PLACEMENT"])] = c
    for p in sorted(glob.glob(os.path.join(a.new, "studies/pdeLevelSet/a2Seam2Dvortex", "2Dvortex_[0-9]*"))):
        t = tok(p); key = (t["FADE_MODE"], t["VE_FADE_PLACEMENT"]); s = ser.get(key)
        if not s:
            say(f"  {key}: serial twin missing"); continue
        say(f"  == {'/'.join(key)}: states {state(os.path.join(s, 'log.leiaLevelSetFoam'))} / "
            f"{state(os.path.join(p, 'log.leiaLevelSetFoam'))}; fallback serial {fallback(s)} np4 {fallback(p)}")
        args = [os.path.join(HERE, "compare_psi_fields.py"), s, p, "0.02"] + (["--celldist", a.celldist] if a.celldist else [])
        rc, d = run(args)
        for line in d[:4] + d[-1:]:
            say("     " + line)
    open(os.path.join(out, "a2_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
