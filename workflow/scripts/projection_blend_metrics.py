#!/usr/bin/env python3
"""Metrics and verdicts of the blended projection extension (config/projectionBlend/, STATUS 11.28).

    projection_blend_metrics.py --studies <study dir>... --out <dir>

Per case (leiaLevelSetFoam on a reversed kinematic flow): the band gradient error E_NARROW_L2_GRAD_PSI
and the band mean of |grad psi| at T/2 (gradPsiError.csv), the shape error E_GEOM_ALPHA_REL at T, the
volume error E_VOL_ALPHA_REL at T/2 and T (leiaLevelSetFoam.csv), the state (foam_log_state.sh).
Writes projectionBlend_cases.csv (one row per case), projectionBlend_orders.csv (least-squares orders per
arm and metric over the rungs) and projectionBlend_verdicts.csv (B1 to B4 of the config headers, applied
as written). L1 and L2 only."""
import argparse
import csv
import glob
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from static_source_metrics import read_csv, lsq_order, write_rows  # noqa: E402

METRICS = ("gradBand_T2", "qBand_T2", "geom_T", "vol_T2", "vol_T")


def at(rows, col, t):
    if not rows:
        return float("nan")
    r = min(rows, key=lambda r: abs(float(r["TIME"]) - t))
    try:
        return abs(float(r[col])) if col != "NARROW_MEAN_MAG_GRAD_PSI" else float(r[col])
    except (TypeError, ValueError):
        return float("nan")


def case_row(case):
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    T = float(tok["END_TIME"])
    adv = read_csv(os.path.join(case, "leiaLevelSetFoam.csv"))
    gp = read_csv(os.path.join(case, "gradPsiError.csv"))
    st = subprocess.run(["bash", os.path.join(HERE, "foam_log_state.sh"), os.path.join(case, "log.leiaLevelSetFoam")],
                        capture_output=True, text=True).stdout.strip()
    tLast = float(adv[-1]["TIME"]) if adv else float("nan")
    arm = f"{tok['VELOCITY_EXTENSION']}-{tok['PB_DISTANCE']}-{tok['SDPLS_SOURCE']}"
    return dict(study=os.path.basename(os.path.dirname(case)), case=os.path.basename(case),
                caseType=tok.get("case", os.path.basename(case).rsplit("_", 1)[0]), arm=arm,
                N=int(float(tok["N_CELLS"])), T=T, tLast=tLast, state=st.split(" ")[0] if st else "",
                gradBand_T2=at(gp, "E_NARROW_L2_GRAD_PSI", T / 2), qBand_T2=at(gp, "NARROW_MEAN_MAG_GRAD_PSI", T / 2),
                geom_T=at(adv, "E_GEOM_ALPHA_REL", T) if tLast >= T * (1 - 1e-9) else float("nan"),
                vol_T2=at(adv, "E_VOL_ALPHA_REL", T / 2),
                vol_T=at(adv, "E_VOL_ALPHA_REL", T) if tLast >= T * (1 - 1e-9) else float("nan"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cases = []
    for s in a.studies:
        for c in sorted(glob.glob(os.path.join(s, "*_0*"))):
            if os.path.isfile(os.path.join(c, "case_params.json")):
                cases.append((c, case_row(c)))
    rows = [r for _, r in cases]
    write_rows(os.path.join(a.out, "projectionBlend_cases.csv"), rows)

    orders, V = [], []
    by = {}
    for r in rows:
        by.setdefault((r["study"], r["arm"]), {})[r["N"]] = r
    for (study, arm), d in sorted(by.items()):
        Ns = sorted(d)
        for m in METRICS:
            vals = [d[N][m] for N in Ns]
            orders.append(dict(study=study, arm=arm, metric=m, N=" ".join(map(str, Ns)),
                               values=" ".join(f"{v:.4e}" for v in vals),
                               order=lsq_order([1 / N for N in Ns], vals) if m != "qBand_T2" else ""))
    write_rows(os.path.join(a.out, "projectionBlend_orders.csv"), orders)
    od = {(o["study"], o["arm"], o["metric"]): o for o in orders}

    def add(item, scope, value, threshold, verdict, note=""):
        V.append(dict(item=item, scope=scope, value=value, threshold=threshold, verdict=verdict, note=note))

    for study in sorted({r["study"] for r in rows}):
        base = by.get((study, "none-second-noSource"), {})
        for dist in ("first", "second"):
            for src in ("noSource", "gradientControl"):
                arm = f"projectionBlend-{dist}-{src}"
                ext = by.get((study, arm), {})
                Ns = sorted(set(base) & set(ext))
                if len(Ns) < 2:
                    continue
                fine = Ns[-2:]
                rg = [ext[N]["gradBand_T2"] / base[N]["gradBand_T2"] for N in fine]
                rs = [ext[N]["geom_T"] / base[N]["geom_T"] for N in fine]
                ob = od[(study, "none-second-noSource", "geom_T")]["order"]
                oe = od[(study, arm, "geom_T")]["order"]
                txt = (f"grad(T/2) ratio {' '.join(f'{x:.2f}' for x in rg)}; shape(T) ratio {' '.join(f'{x:.2f}' for x in rs)}; "
                       f"shape order {oe:.2f} against {ob:.2f}")
                if dist == "second" and src == "noSource":
                    fail = (ob - oe > 0.3) or (rs[-1] > 2)
                    ok = all(x <= 0.5 for x in rg) and all(x <= 1.2 for x in rs) and (oe >= ob - 0.3)
                    add("B1", f"{study} {arm}", txt, "grad <= 0.5x, shape <= 1.2x, order >= base - 0.3",
                        "FAIL" if fail else ("PASS" if ok else "OUT OF BAND"))
                elif src == "noSource":
                    add("B2", f"{study} {arm}", txt, "recorded", "RECORDED")
                else:
                    alone = by.get((study, f"projectionBlend-{dist}-noSource"), {})
                    r3 = [ext[N]["gradBand_T2"] / alone[N]["gradBand_T2"] for N in fine if N in alone]
                    add("B3", f"{study} {arm}", txt + f"; grad(T/2) against extension alone {' '.join(f'{x:.2f}' for x in r3)}",
                        "<= 1 against the extension alone", "PASS" if r3 and all(x <= 1 for x in r3) else "FAIL")
        # B4: the null pair and the states
        nulls = []
        for c1, r1 in cases:
            if r1["study"] != study or r1["arm"] != "none-first-noSource":
                continue
            c2 = next((c for c, r in cases if r["study"] == study and r["arm"] == "none-second-noSource" and r["N"] == r1["N"]), None)
            if c2:
                p = subprocess.run([sys.executable, os.path.join(HERE, "compare_metrics_csv.py"),
                                    os.path.join(c1, "leiaLevelSetFoam.csv"), os.path.join(c2, "leiaLevelSetFoam.csv"),
                                    "--tol", "0", "--skip", "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"], capture_output=True, text=True)
                nulls.append(p.returncode == 0)
        states = [r["state"] for r in rows if r["study"] == study]
        add("B4", study, f"null pairs {sum(nulls)}/{len(nulls)} bit-identical; states "
            + ", ".join(f"{s} {states.count(s)}" for s in sorted(set(states))), "all identical, all COMPLETED",
            "PASS" if nulls and all(nulls) and set(states) <= {"COMPLETED"} else "FAIL")
    write_rows(os.path.join(a.out, "projectionBlend_verdicts.csv"), V)
    for v in V:
        print(v["item"], v["verdict"], "|", v["scope"], "|", v["value"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
