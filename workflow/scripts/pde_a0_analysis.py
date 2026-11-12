#!/usr/bin/env python3
"""pde_a0_analysis.py: the read-out of the A0.1 gates of the PDE level-set plan.

    pde_a0_analysis.py --old <clone> --new <clone> [--out <dir>]

A0.1a (a0Identity2Dvortex): every case CSV of the old clone against the same case of the new
clone at tolerance 0. A0.1b (a0FluxForm2Dvortex, new clone): the inert projectionBlend in the
correction form against none (tolerance 0), in the interpolate form against none (the size of
the far-field flux confound at T/2 and T), and none with both token values (tolerance 0).
Writes <out>/a0_summary.txt (default <new>/studies/pdeLevelSet/a0). Stdlib only.
"""
import argparse, csv, glob, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"
CSVS = ("leiaLevelSetFoam.csv", "gradPsiError.csv")


def cmp(a, b, tol="0"):
    if not (os.path.isfile(a) and os.path.isfile(b)):
        return "MISSING", ""
    p = subprocess.run([sys.executable, os.path.join(HERE, "compare_metrics_csv.py"), a, b, "--tol", tol, "--skip", SKIP],
                       capture_output=True, text=True)
    return ("PASS" if p.returncode == 0 else "DIFF"), " | ".join(p.stdout.strip().splitlines()[-3:])


def tok(case):
    return {k: str(v) for k, v in json.load(open(os.path.join(case, "case_params.json")))["tokens"].items()}


def row_at(path, t):
    rows = list(csv.DictReader(open(path))) if os.path.isfile(path) and os.path.getsize(path) else []
    best = min(rows, key=lambda r: abs(float(r["TIME"]) - t), default={})
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "a0")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))

    say("A0.1a refactor identity (old clone against new clone, tolerance 0):")
    sub = os.path.join("studies", "pdeLevelSet", "a0Identity2Dvortex")
    n = bad = 0
    for case in sorted(glob.glob(os.path.join(a.old, sub, "2Dvortex_[0-9]*"))):
        name = os.path.basename(case)
        ve = tok(case).get("VELOCITY_EXTENSION")
        for f in CSVS:
            v, d = cmp(os.path.join(case, f), os.path.join(a.new, sub, name, f))
            n += 1; bad += v != "PASS"
            say(f"  {name} {ve:16s} {f:22s} {v}  {d[:160]}")
    say(f"  -> {n - bad} of {n} identical: {'PASS' if n and not bad else 'FAIL'}")

    say()
    say("A0.1b flux form (new clone): inert projectionBlend (PB_W1 0, PB_W2 1e-9) against none")
    sub = os.path.join(a.new, "studies", "pdeLevelSet", "a0FluxForm2Dvortex")
    cases = {}
    for case in sorted(glob.glob(os.path.join(sub, "2Dvortex_[0-9]*"))):
        t = tok(case)
        cases[(t["N_CELLS"], t["VELOCITY_EXTENSION"], t.get("VE_FLUX_FORM", "?"))] = case
    for N in sorted({k[0] for k in cases}, key=int):
        ref = cases.get((N, "none", "interpolate"))
        T = float(tok(ref)["END_TIME"]) if ref else 0.5
        for label, key, want in (("(3) none, correction token", (N, "none", "correction"), "identical"),
                                 ("(1) inert PB, correction", (N, "projectionBlend", "correction"), "identical"),
                                 ("(2) inert PB, interpolate", (N, "projectionBlend", "interpolate"), "the confound")):
            c = cases.get(key)
            if not (ref and c):
                say(f"  N {N} {label}: missing"); continue
            for f in CSVS:
                v, d = cmp(os.path.join(ref, f), os.path.join(c, f))
                say(f"  N {N} {label:28s} {f:22s} {v} (expected {want})  {d[:150]}")
            if key[2] == "interpolate" and key[1] != "none":
                for t in (T / 2, T):
                    r0, r1 = row_at(os.path.join(ref, "leiaLevelSetFoam.csv"), t), row_at(os.path.join(c, "leiaLevelSetFoam.csv"), t)
                    g0, g1 = row_at(os.path.join(ref, "gradPsiError.csv"), t), row_at(os.path.join(c, "gradPsiError.csv"), t)
                    say(f"      t = {t:g}: shape {r0.get('E_GEOM_ALPHA_REL')} -> {r1.get('E_GEOM_ALPHA_REL')}; "
                        f"volume {r0.get('E_VOL_ALPHA_REL')} -> {r1.get('E_VOL_ALPHA_REL')}; "
                        f"band grad L1 {g0.get('E_NARROW_L1_GRAD_PSI')} -> {g1.get('E_NARROW_L1_GRAD_PSI')}")
    open(os.path.join(out, "a0_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
