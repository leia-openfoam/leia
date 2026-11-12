#!/usr/bin/env python3
"""pde_step0_analysis.py: the read-out of step 0 of the PDE level-set plan (Eulerian transport alone).

    pde_step0_analysis.py [--studies studies/pdeLevelSet] [--out <dir>] [--celldist <file>]

Reads the studies of config/pdeLevelSet/s0*.yaml and writes, under <out> (default
<studies>/step0):
  s0_translation.csv   one row per case: N, CFL, nDefCorr, np, state, steps, psi solves and the
                       solves with a nonzero initial residual and no iteration, and at T/2 and T the
                       band gradient error (E_NARROW_L1/L2_GRAD_PSI, E_NARROW_MEAN_GRAD_PSI), and at T
                       the shape and volume errors against the exact translated state
  s0_orders.csv        the observed order between successive rungs, per (CFL, nDefCorr), np 4 only
  s0_seam.csv          serial against np 4, per case pair: the largest relative column difference
  s0_summary.txt       the pre-registered verdicts of s0Translation2D.yaml, and the seam probe
                       (compare_psi_fields.py at every written time of the vortex pair)
L1 and L2 only (CLAUDE.md: never an L_inf order). Stdlib only. The case state comes from
workflow/scripts/foam_log_state.sh, never from a grep of its own.
"""
import argparse, csv, glob, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SOLVER = "leiaLevelSetFoam"
GRAD_COLS = ("E_NARROW_L1_GRAD_PSI", "E_NARROW_L2_GRAD_PSI", "E_NARROW_MEAN_GRAD_PSI")
ITER_RE = re.compile(r"Solving for psi, Initial residual = ([-+0-9.eE]+), Final residual = [-+0-9.eE]+, No Iterations (\d+)")


def tokens(case):
    with open(os.path.join(case, "case_params.json")) as f:
        return {k: str(v) for k, v in json.load(f)["tokens"].items()}


def state(case):
    log = os.path.join(case, f"log.{SOLVER}")
    r = subprocess.run([os.path.join(HERE, "foam_log_state.sh"), log], capture_output=True, text=True)
    line = (r.stdout.strip().splitlines() or ["?"])[0]
    m = re.search(r"steps=(\d+)", line)
    return line.split()[0], int(m.group(1)) if m else 0


def iterations(case):
    """(psi solves, solves with a nonzero initial residual and zero iterations) from the log."""
    n = dead = 0
    log = os.path.join(case, f"log.{SOLVER}")
    if not os.path.isfile(log):
        return 0, 0
    for line in open(log, errors="ignore"):
        m = ITER_RE.search(line)
        if m:
            n += 1
            if int(m.group(2)) == 0 and float(m.group(1)) > 0.0:
                dead += 1
    return n, dead


def rows_of(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path) as f:
        return [r for r in csv.DictReader(f)]


def at_time(rows, t):
    """The row whose TIME is closest to t."""
    best = None
    for r in rows:
        try:
            d = abs(float(r["TIME"]) - t)
        except (KeyError, ValueError):
            continue
        if best is None or d < best[0]:
            best = (d, r)
    return best[1] if best else {}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def order(e1, e2, n1, n2):
    if not (e1 > 0 and e2 > 0):
        return float("nan")
    return math.log(e1 / e2) / math.log(n2 / n1)


def translation_rows(studies):
    out = []
    for study, npv in (("s0Translation2D", 4), ("s0Translation2Dserial", 1)):
        for case in sorted(glob.glob(os.path.join(studies, study, "2Dtranslation_[0-9]*"))):
            t = tokens(case)
            T = fnum(t.get("END_TIME", 0.5))
            st, steps = state(case)
            ns, dead = iterations(case)
            g = rows_of(os.path.join(case, "gradPsiError.csv"))
            a = rows_of(os.path.join(case, f"{SOLVER}.csv"))
            gh, gT, aT = at_time(g, T / 2), at_time(g, T), at_time(a, T)
            row = {"study": study, "case": os.path.basename(case), "np": npv, "N": int(t["N_CELLS"]),
                   "CFL": t["CFL"], "nDefCorr": t.get("N_DEF_CORR", "3"), "state": st, "steps": steps,
                   "psiSolves": ns, "psiSolvesNoIterNonzeroRes": dead, "T": T}
            for c in GRAD_COLS:
                row[c + "_Thalf"] = gh.get(c, "")
                row[c + "_T"] = gT.get(c, "")
            for c in ("E_GEOM_ALPHA_REL", "E_VOL_ALPHA_REL", "E_VOL_ALPHA_SIGNED_REL", "E_BOUND_ALPHA"):
                row[c + "_T"] = aT.get(c, "")
            out.append(row)
    return out


def write(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", default=os.path.join(REPO, "studies", "pdeLevelSet"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--celldist", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.studies, "step0")
    os.makedirs(out, exist_ok=True)
    lines = []
    say = lambda s="": (lines.append(s), print(s))

    rows = translation_rows(a.studies)
    write(os.path.join(out, "s0_translation.csv"), rows)
    say(f"STEP 0.1, one-way translation: {len(rows)} cases")
    bad = [r for r in rows if r["state"] != "COMPLETED"]
    dead = [r for r in rows if r["psiSolvesNoIterNonzeroRes"] > 0]
    say(f"  states: {len(rows) - len(bad)} COMPLETED, {len(bad)} other " + " ".join(f"{r['study']}/{r['case']}={r['state']}" for r in bad))
    say(f"  psi solves with a nonzero residual and no iteration: {sum(r['psiSolvesNoIterNonzeroRes'] for r in rows)} "
        f"in {len(dead)} cases (must be 0 before any metric is read)")

    # orders, np 4
    orows = []
    par = [r for r in rows if r["np"] == 4 and r["state"] == "COMPLETED"]
    for key in sorted({(r["CFL"], r["nDefCorr"]) for r in par}):
        arm = sorted([r for r in par if (r["CFL"], r["nDefCorr"]) == key], key=lambda r: r["N"])
        for m in ("E_NARROW_L1_GRAD_PSI_T", "E_NARROW_L2_GRAD_PSI_T", "E_NARROW_L1_GRAD_PSI_Thalf",
                  "E_NARROW_L2_GRAD_PSI_Thalf", "E_GEOM_ALPHA_REL_T", "E_VOL_ALPHA_REL_T"):
            for r1, r2 in zip(arm, arm[1:]):
                orows.append({"CFL": key[0], "nDefCorr": key[1], "metric": m, "N1": r1["N"], "N2": r2["N"],
                              "e1": r1[m], "e2": r2[m], "order": round(order(fnum(r1[m]), fnum(r2[m]), r1["N"], r2["N"]), 3)})
    write(os.path.join(out, "s0_orders.csv"), orows)
    say("  observed orders, np 4 (N 32-64 / 64-128 / 128-256):")
    for key in sorted({(o["CFL"], o["nDefCorr"]) for o in orows}):
        for m in ("E_NARROW_L1_GRAD_PSI_T", "E_NARROW_L2_GRAD_PSI_T", "E_GEOM_ALPHA_REL_T", "E_VOL_ALPHA_REL_T"):
            os_ = [o for o in orows if (o["CFL"], o["nDefCorr"]) == key and o["metric"] == m]
            fin = [o["e2"] for o in os_][-1:] or [""]
            say(f"    CFL {key[0]} nDefCorr {key[1]} {m:24s} " + " / ".join(f"{o['order']:.2f}" for o in os_) + f"   (finest {fin[0]})")

    # pre-registered verdict (s0Translation2D.yaml)
    say("  VERDICT (pre-registered: band gradient L1 and L2 at T order >= 0.8 over 64-128-256, shape order >= 1.0, seam <= 1e-12):")
    for key in sorted({(o["CFL"], o["nDefCorr"]) for o in orows}):
        def last2(m):
            return [o["order"] for o in orows if (o["CFL"], o["nDefCorr"]) == key and o["metric"] == m and o["N1"] >= 64]
        g1, g2, sh = last2("E_NARROW_L1_GRAD_PSI_T"), last2("E_NARROW_L2_GRAD_PSI_T"), last2("E_GEOM_ALPHA_REL_T")
        ok = lambda v, thr: bool(v) and all((not math.isnan(x)) and x >= thr for x in v)
        fals = (g1 and max(g1) <= 0.3) or (g2 and max(g2) <= 0.3)
        verdict = "FALSIFIED (transport defect)" if fals else ("PASS" if ok(g1, 0.8) and ok(g2, 0.8) and ok(sh, 1.0) else "FAIL (between the lines)")
        say(f"    CFL {key[0]} nDefCorr {key[1]}: grad L1 {g1} L2 {g2} shape {sh} -> {verdict}")

    # seam: serial against np 4, metric CSVs
    srows = []
    cmpr = os.path.join(HERE, "compare_metrics_csv.py")
    sers = {(r["N"], r["CFL"], r["nDefCorr"]): r for r in rows if r["np"] == 1}
    for r in rows:
        if r["np"] != 4:
            continue
        s = sers.get((r["N"], r["CFL"], r["nDefCorr"]))
        if not s:
            continue
        for f in (f"{SOLVER}.csv", "gradPsiError.csv"):
            A = os.path.join(a.studies, s["study"], s["case"], f)
            B = os.path.join(a.studies, r["study"], r["case"], f)
            if not (os.path.isfile(A) and os.path.isfile(B)):
                continue
            p = subprocess.run([sys.executable, cmpr, A, B, "--tol", "1e-12", "--skip", "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"],
                               capture_output=True, text=True)
            srows.append({"N": r["N"], "CFL": r["CFL"], "nDefCorr": r["nDefCorr"], "csv": f,
                          "verdict": "PASS" if p.returncode == 0 else "FAIL",
                          "detail": " | ".join(p.stdout.strip().splitlines()[-3:])})
    write(os.path.join(out, "s0_seam.csv"), srows)
    nf = sum(1 for x in srows if x["verdict"] != "PASS")
    say(f"  seam, serial against np 4 (tol 1e-12): {len(srows)} CSV pairs, {nf} FAIL")
    for x in [x for x in srows if x["verdict"] != "PASS"][:8]:
        say(f"    N {x['N']} CFL {x['CFL']} nDefCorr {x['nDefCorr']} {x['csv']}: {x['detail'][:200]}")

    # step 0.2: the vortex seam probe
    say()
    say("STEP 0.2, vortex seam probe (serial against np 4, psi cell by cell):")
    ser = sorted(glob.glob(os.path.join(a.studies, "s0Seam2Dvortexserial", "2Dvortex_[0-9]*")))
    pa = sorted(glob.glob(os.path.join(a.studies, "s0Seam2Dvortex", "2Dvortex_[0-9]*")))
    if ser and pa:
        def tdirs(c):
            out = set()
            for d in glob.glob(os.path.join(c, "[0-9]*")):
                try:
                    float(os.path.basename(d))
                except ValueError:
                    continue          # 0.org and the like are not time directories
                if os.path.isdir(d):
                    out.add(os.path.basename(d))
            return out
        times = sorted(tdirs(ser[0]) & tdirs(pa[0]), key=float)
        for t in times:
            if float(t) == 0.0:
                continue
            cmd = [sys.executable, os.path.join(HERE, "compare_psi_fields.py"), ser[0], pa[0], t]
            if a.celldist:
                cmd += ["--celldist", a.celldist]
            p = subprocess.run(cmd, capture_output=True, text=True)
            for l in (p.stdout + p.stderr).strip().splitlines():
                say("    " + l)
        for case in ser + pa:
            st, steps = state(case)
            say(f"    {os.path.relpath(case, a.studies)}: {st} steps={steps}")
        A = os.path.join(ser[0], f"{SOLVER}.csv"); B = os.path.join(pa[0], f"{SOLVER}.csv")
        p = subprocess.run([sys.executable, cmpr, A, B, "--tol", "1e-12", "--skip", "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"],
                           capture_output=True, text=True)
        say("    metric CSV: " + " | ".join(p.stdout.strip().splitlines()[-4:]))
    else:
        say("    the vortex pair is missing")
    open(os.path.join(out, "s0_summary.txt"), "w").write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
