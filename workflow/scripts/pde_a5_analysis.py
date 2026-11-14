#!/usr/bin/env python3
"""pde_a5_analysis.py: the read-out of step A proper, round 1 (q != 1 inputs, the ellipse).

    pde_a5_analysis.py --old <clone a3> --new <clone a5> [--out <dir>]

1. Inertness: aStatic2Dedge and a0Identity2Dvortex of the new clone against the old one (tolerance 0).
2. The input check: every non-distance case of aStatic2Dq, aStatic2DqPar, aStatic2Dellipse has the log
   line "testInput ... psi replaced".
3. Q4: closestPoint foot-pointed and fallback cells per (input, band) at each N (aStatic2Dq, newton).
4. np 4 against serial: aStatic2DqPar against the same arm of aStatic2Dq.
The order tables come from pde_a_static_analysis.py per study. Writes <out>/a5_summary.txt.
"""
import argparse, glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pde_a2_analysis import tok  # noqa: E402
from pde_a3_analysis import cmp, counters, exact_row  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "a5")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda root, st: os.path.join(root, "studies", "pdeLevelSet", st)

    say("1. Inertness (new clone against old clone, tolerance 0):")
    n = bad = 0
    for st, files in (("aStatic2Dedge", ("leiaTestVelocityExtension.csv", "leiaTestVelocityExtensionExact.csv")),
                      ("a0Identity2Dvortex", ("leiaLevelSetFoam.csv", "gradPsiError.csv"))):
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*"))):
            for f in files:
                v, d = cmp(os.path.join(S(a.old, st), os.path.basename(c), f), os.path.join(c, f))
                n += 1; bad += v != "PASS"
                if v != "PASS":
                    say(f"  {st} {os.path.basename(c)} {f}: {v} {' | '.join(d[-2:])[:160]}")
    say(f"  -> {n - bad} of {n} identical: {'PASS' if n and not bad else 'FAIL'}")

    say()
    say("2. The input was applied (log line 'testInput ... psi replaced'):")
    for st in ("aStatic2Dq", "aStatic2DqPar", "aStatic2Dellipse"):
        ok = miss = 0
        for c in glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*")):
            if tok(c).get("VE_TEST_INPUT", "distance") == "distance":
                continue
            txt = "".join(open(f, errors="ignore").read() for f in glob.glob(os.path.join(c, "log.leiaTestVelocityExtension")))
            if "psi replaced" in txt:
                ok += 1
            else:
                miss += 1; say(f"  MISSING in {st}/{os.path.basename(c)}")
        say(f"  {st}: {ok} applied, {miss} missing")

    say()
    say("3. closestPoint (newton) cells per input and band: (foot-pointed, via halo, fallback):")
    for st in ("aStatic2Dq", "aStatic2Dellipse"):
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*")),
                        key=lambda c: (tok(c)["VE_TEST_INPUT"], tok(c)["FADE_MODE"], int(tok(c)["N_CELLS"]))):
            t = tok(c)
            if t["VELOCITY_EXTENSION"] != "closestPoint" or t["CP_STEP"] != "newton":
                continue
            fb, fl, _ = counters(c)
            say(f"  {st:16s} {t['VE_TEST_INPUT']:15s} {t['FADE_MODE']:8s} N {t['N_CELLS']:>4s} {fb}")

    say()
    say("4. np 4 against serial (aStatic2DqPar against aStatic2Dq): (np4 - serial)/serial")
    ser = {}
    for c in glob.glob(os.path.join(S(a.new, "aStatic2Dq"), "2Dvortex_[0-9]*")):
        t = tok(c)
        if (t["VELOCITY_EXTENSION"], t["CP_STEP"], t["FADE_MODE"]) == ("closestPoint", "newton", "layer"):
            ser[(t["VE_TEST_INPUT"], t["N_CELLS"])] = exact_row(c)
    for c in sorted(glob.glob(os.path.join(S(a.new, "aStatic2DqPar"), "2Dvortex_[0-9]*")),
                    key=lambda c: (tok(c)["VE_TEST_INPUT"], int(tok(c)["N_CELLS"]))):
        t = tok(c); r = exact_row(c); s = ser.get((t["VE_TEST_INPUT"], t["N_CELLS"]), {})
        cells = []
        for k in ("E_EXT_L1_2H", "E_EXT_L1_4H", "E_FLUX_L1_4H"):
            try:
                cells.append(f"{k} {float(r[k]):.3e} ({(float(r[k]) - float(s[k])) / float(s[k]):+.1%})")
            except (KeyError, ValueError, ZeroDivisionError):
                cells.append(f"{k} ?")
        fb, _, _ = counters(c)
        say(f"  {t['VE_TEST_INPUT']:15s} N {t['N_CELLS']:>4s}  " + "  ".join(cells) + f"  cells {fb}")
    open(os.path.join(out, "a5_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
