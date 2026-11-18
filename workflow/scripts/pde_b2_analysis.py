#!/usr/bin/env python3
"""pde_b2_analysis.py: the read-out of the step-B discriminators (config bDiscr2Dstrain, D0 to D4).

    pde_b2_analysis.py --new <clone> [--ref <clone of round 1>] [--out <dir>]

Per case of bDiscr2Dstrain: E_Q_L1 at T/2 and T, the onset (the first row where E_Q_L1 > 2 x its t = 0
value; step index and time), and the mean growth factor per step from the onset to T. D0 compares the
base arm with the projectionBlend N 256 case of bMarker2Dstrain in --ref. Writes <out>/b2_summary.txt.
Stdlib only.
"""
import argparse, csv, glob, json, math, os, sys


def rows(case):
    f = os.path.join(case, "interfaceMarkers.csv")
    if not (os.path.isfile(f) and os.path.getsize(f)):
        return []
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(f))]


def at(rs, t):
    return min(rs, key=lambda r: abs(r["TIME"] - t))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--ref"); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "b2")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    res = {}
    for c in sorted(glob.glob(os.path.join(a.new, "studies", "pdeLevelSet", "bDiscr2Dstrain", "markerTest2D_*"))):
        t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        key = (float(t["CFL"]), str(t["PSI_DDT"]), str(t["DIV_GRAD"]), str(t["VE_FLUX_FORM"]))
        rs = rows(c)
        if not rs:
            say(f"{key}: no CSV"); continue
        T = float(t["END_TIME"]); e0 = rs[0]["E_Q_L1"]
        on = next((i for i, r in enumerate(rs) if r["E_Q_L1"] > 2 * e0), None)
        g = (rs[-1]["E_Q_L1"] / rs[on]["E_Q_L1"]) ** (1 / max(len(rs) - 1 - on, 1)) if on is not None else float("nan")
        res[key] = dict(case=c, rows=rs, eh=at(rs, T / 2)["E_Q_L1"], eT=at(rs, T)["E_Q_L1"],
                        ph=at(rs, T / 2)["E_POS_L1"], pT=at(rs, T)["E_POS_L1"],
                        on_step=on, on_t=rs[on]["TIME"] if on is not None else float("nan"),
                        g=g, steps=len(rs) - 1, tend=rs[-1]["TIME"])
    say("CFL    ddt               grad            flux         steps  E_Q_L1(T/2)  E_Q_L1(T)  E_POS_L1(T)  onset step / t   growth/step")
    for k in sorted(res):
        r = res[k]
        say(f"{k[0]:<6} {k[1]:17s} {k[2]:15s} {k[3]:12s} {r['steps']:5d}  {r['eh']:.3e}    {r['eT']:.3e}  {r['pT']:.3e}    "
            f"{str(r['on_step']):>5s} / {r['on_t']:.3f}   {r['g']:.4f}" + ("" if abs(r["tend"] - 0.5) < 1e-6 else f"  (ends at t {r['tend']:.3f})"))
    base = (0.5, "CrankNicolsonDdt", "CellLimitedLSQ", "correction")
    say()
    if a.ref and base in res:
        ref = None
        for c in glob.glob(os.path.join(a.ref, "studies", "pdeLevelSet", "bMarker2Dstrain", "markerTest2D_*")):
            t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
            if t["VELOCITY_EXTENSION"] == "projectionBlend" and int(t["N_CELLS"]) == 256:
                ref = rows(c)
        if ref:
            rb = res[base]["rows"]
            n = min(len(ref), len(rb))
            rel = max(abs(rb[i][k] - ref[i][k]) / abs(ref[i][k]) for i in range(n) for k in ("E_Q_L1", "E_POS_L1") if ref[i][k] != 0)
            say(f"D0 control against round 1: {len(rb)} and {len(ref)} rows, max rel. diff {rel:.2e} "
                f"[{'PASS' if rel <= 1e-12 and len(rb) == len(ref) else 'FAIL'}]")
    # D1: CFL sweep with the base schemes
    s = [res.get((c, "CrankNicolsonDdt", "CellLimitedLSQ", "correction")) for c in (0.5, 0.25, 0.125)]
    if all(s) and all(x["on_step"] is not None for x in s):
        st = [x["on_step"] for x in s]; tt = [x["on_t"] for x in s]
        ds = (max(st) - min(st)) / min(st); dt = (max(tt) - min(tt)) / min(tt)
        say(f"D1 onset step {st} (spread {ds:.0%}), onset time {[round(x, 4) for x in tt]} (spread {dt:.0%}): "
            + ("PER STEP" if ds < 0.3 else "") + (" PER TIME" if dt < 0.3 else "") + ("" if ds < 0.3 or dt < 0.3 else " neither"))
    else:
        say("D1: an arm of the CFL sweep has no onset (no growth beyond 2 x the t = 0 value): "
            + ", ".join(f"CFL {c}: onset {x['on_step'] if x else 'missing'}" for c, x in zip((0.5, 0.25, 0.125), s)))
    for tag, idx, alt, ref_v in (("D2 Euler", 1, "EulerDdt", "CrankNicolsonDdt"),
                                 ("D3 unlimited", 2, "UnlimitedLSQ", "CellLimitedLSQ"),
                                 ("D4 interpolate", 3, "interpolate", "correction")):
        rat = []
        for k, r in res.items():
            if k[idx] != ref_v:
                continue
            k2 = list(k); k2[idx] = alt; k2 = tuple(k2)
            if k2 in res:
                rat.append((k, res[k2]["eT"] / r["eT"]))
        say(f"{tag}: E_Q_L1(T) ratio (alternative / reference) over the other factors: "
            + ", ".join(f"{x:.3g}" for _, x in sorted(rat)))
    open(os.path.join(out, "b2_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
