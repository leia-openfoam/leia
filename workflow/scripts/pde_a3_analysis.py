#!/usr/bin/env python3
"""pde_a3_analysis.py: the read-out of step A0.3 (closestPoint repair 3) of the PDE level-set plan.

    pde_a3_analysis.py --old <clone a2> --new <clone a3> [--celldist <file>] [--out <dir>]

1. Inertness: aStatic2Dedge, a2Seam2Dvortex{,serial} and a0Identity2Dvortex of the new clone against
   the same case of the old clone, tolerance 0 (timing columns skipped).
2. Control: aStatic2Dfill (steadyUpwind, distance, cap 1.5) against the extendThroughFade/levelSet
   arm of aStatic2Dedge (new clone), tolerance 0.
3. The fill counters of the logs (closestPoint fallback cells; footPointTaylor filled/unreached).
4. Seam: a3Seam2Dvortexserial against a3Seam2Dvortex, psi at t = 0.02.
5. np 4 against serial (aStatic2DfillPar against aStatic2Dfill, cap 1.5, distance): the relative
   difference of E_EXT_L1_4H, E_EXT_L1_W and E_FLUX_L1_4H per N and fill.
Writes <out>/a3_summary.txt (default <new>/studies/pdeLevelSet/a3). Stdlib only.
"""
import argparse, csv, glob, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"
sys.path.insert(0, HERE)
from pde_a2_analysis import tok, run, state  # noqa: E402


def cmp(a, b):
    if not (os.path.isfile(a) and os.path.isfile(b)):
        return "MISSING", []
    rc, d = run([os.path.join(HERE, "compare_metrics_csv.py"), a, b, "--tol", "0", "--skip", SKIP])
    return ("PASS" if rc == 0 else "DIFF"), d


def counters(case):
    txt = "".join(open(f, errors="ignore").read() for f in glob.glob(os.path.join(case, "log.leia*")))
    fb = re.findall(r"(\d+) foot-pointed band cells \((\d+) via halo\), (\d+) fallback cells", txt)
    fl = re.findall(r"footPointTaylor: (\d+) cells filled in (\d+) sweeps, (\d+) unreached", txt)
    return (fb[-1] if fb else None), (fl[-1] if fl else None), len(fl)


def exact_row(case):
    f = os.path.join(case, "leiaTestVelocityExtensionExact.csv")
    return list(csv.DictReader(open(f)))[0] if os.path.isfile(f) and os.path.getsize(f) else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
    ap.add_argument("--celldist"); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "a3")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda root, st: os.path.join(root, "studies", "pdeLevelSet", st)

    say("1. Inertness (new clone against old clone, tolerance 0):")
    n = bad = 0
    for st, files in (("aStatic2Dedge", ("leiaTestVelocityExtension.csv", "leiaTestVelocityExtensionExact.csv")),
                      ("a2Seam2Dvortex", ("leiaLevelSetFoam.csv", "gradPsiError.csv")),
                      ("a2Seam2Dvortexserial", ("leiaLevelSetFoam.csv", "gradPsiError.csv")),
                      ("a0Identity2Dvortex", ("leiaLevelSetFoam.csv", "gradPsiError.csv"))):
        sn = sb = 0
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*"))):
            for f in files:
                v, d = cmp(os.path.join(S(a.old, st), os.path.basename(c), f), os.path.join(c, f))
                sn += 1; sb += v != "PASS"
                if v != "PASS":
                    say(f"  {st} {os.path.basename(c)} {f}: {v} {' | '.join(d[-2:])[:160]}")
        say(f"  {st}: {sn - sb} of {sn} identical")
        n += sn; bad += sb
    say(f"  -> {n - bad} of {n} identical: {'PASS' if n and not bad else 'FAIL'}")

    say()
    say("2. Control: aStatic2Dfill (steadyUpwind, distance, cap 1.5) against aStatic2Dedge extendThroughFade/levelSet:")
    edge = {}
    for c in glob.glob(os.path.join(S(a.new, "aStatic2Dedge"), "2Dvortex_[0-9]*")):
        t = tok(c)
        if t["VELOCITY_EXTENSION"] == "closestPoint" and t["FADE_MODE"] == "levelSet" and t["VE_FADE_PLACEMENT"] == "extendThroughFade":
            edge[t["N_CELLS"]] = c
    for c in sorted(glob.glob(os.path.join(S(a.new, "aStatic2Dfill"), "2Dvortex_[0-9]*"))):
        t = tok(c)
        if (t["CP_FILL"], t["CP_STEP"], t["CP_MAX_DIST"]) != ("steadyUpwind", "distance", "1.5"):
            continue
        for f in ("leiaTestVelocityExtension.csv", "leiaTestVelocityExtensionExact.csv"):
            v, d = cmp(os.path.join(edge.get(t["N_CELLS"], "/none"), f), os.path.join(c, f))
            say(f"  N {t['N_CELLS']:>4s} {f:38s} {v}")

    say()
    say("3. Counters (last call): (foot-pointed, via halo, fallback); footPointTaylor (filled, sweeps, unreached):")
    for st in ("aStatic2Dfill", "aStatic2DfillPar", "a3Seam2Dvortex", "a3Seam2Dvortexserial"):
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*"))):
            t = tok(c); fb, fl, nfl = counters(c)
            say(f"  {st:22s} N {t['N_CELLS']:>4s} {t.get('FADE_MODE', ''):8s} fill {t.get('CP_FILL'):15s} "
                f"step {t.get('CP_STEP'):8s} cap {t.get('CP_MAX_DIST'):4s} {fb} {fl} ({nfl} fill calls)")

    say()
    say("4. Seam: psi at t = 0.02, serial against np 4, relative to max|psi| (a2: levelSet 7.0e-8, layer 5.7e-8):")
    ser = {tok(c)["FADE_MODE"]: c for c in glob.glob(os.path.join(S(a.new, "a3Seam2Dvortexserial"), "2Dvortex_[0-9]*"))}
    for p in sorted(glob.glob(os.path.join(S(a.new, "a3Seam2Dvortex"), "2Dvortex_[0-9]*"))):
        fm = tok(p)["FADE_MODE"]; s = ser.get(fm)
        if not s:
            say(f"  {fm}: serial twin missing"); continue
        say(f"  == {fm}: states {state(os.path.join(s, 'log.leiaLevelSetFoam'))} / {state(os.path.join(p, 'log.leiaLevelSetFoam'))}")
        rc, d = run([os.path.join(HERE, "compare_psi_fields.py"), s, p, "0.02"] + (["--celldist", a.celldist] if a.celldist else []))
        for line in d[:4] + d[-1:]:
            say("     " + line)

    say()
    say("5. np 4 against serial (cap 1.5, distance): relative difference (np4 - serial)/serial:")
    ser = {}
    for c in glob.glob(os.path.join(S(a.new, "aStatic2Dfill"), "2Dvortex_[0-9]*")):
        t = tok(c)
        if t["CP_STEP"] == "distance" and t["CP_MAX_DIST"] == "1.5":
            ser[(t["CP_FILL"], t["N_CELLS"])] = exact_row(c)
    for c in sorted(glob.glob(os.path.join(S(a.new, "aStatic2DfillPar"), "2Dvortex_[0-9]*")), key=lambda c: (tok(c)["CP_FILL"], int(tok(c)["N_CELLS"]))):
        t = tok(c); r = exact_row(c); s = ser.get((t["CP_FILL"], t["N_CELLS"]), {})
        cells = []
        for k in ("E_EXT_L1_4H", "E_EXT_L1_W", "E_FLUX_L1_4H", "E_FADEFLUX_L1"):
            try:
                cells.append(f"{k} {float(r[k]):.3e} ({(float(r[k]) - float(s[k])) / float(s[k]):+.1%})")
            except (KeyError, ValueError, ZeroDivisionError):
                cells.append(f"{k} ?")
        say(f"  {t['CP_FILL']:15s} N {t['N_CELLS']:>4s}  " + "  ".join(cells))
    open(os.path.join(out, "a3_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
