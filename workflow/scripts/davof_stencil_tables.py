#!/usr/bin/env python3
"""The tables of the ring normal's stencil study (2026-11-08; the configs
config/davof/stencil_*.yaml, written by make_stencil_configs.py, carry the
predictions). One table per test with the variants (TPF_NORMAL_STENCIL) as rows,
printed as Markdown and written as CSV:

    python3 workflow/scripts/davof_stencil_tables.py --studies studies/davof --out <csv dir>
        [--port-plane <the port's port_enrightKeys_plane_leia study, for the identity G-0>]

Per case (the converter's advectorVoF.csv, one row per step; the log; .leia_launch):
the end values of L1(alphaf), L1(n) (LNormalDiff1), L1(x) (LCentreDiff1) and Eshape,
the run maxima of L1(alphaf) and L1(n), min alpha and the excess over 1 after the
bounding, the mass deviation at the end, the interface-cell count, the status, the
wall clock of the solve (.leia_launch end - start), the per-step median of the log's
timeStepTime, the stencil counters of the log (cells dropped by the angle filter, fit
fallbacks, remote rings; the first line per step). The growth factors of test 1b:
g_early = (L1 at step 6 / L1 at step 1)^(1/5), g_late = (L1 at the last step / L1 at
the middle step)^(1/(last - middle)), over the rows with TIME > 0.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import statistics

VARIANTS = ["production", "twoRings", "kernelGaussian", "kernelWendland", "angle30",
            "angle60", "planeFitMean", "planeFitOwn", "ringQuadric"]


def f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def order(ns, es):
    pts = [(1.0 / n, e) for n, e in zip(ns, es) if e is not None and e > 0]
    if len(pts) < 2:
        return None
    lh = [math.log(p[0]) for p in pts]
    le = [math.log(p[1]) for p in pts]
    mh, me = sum(lh) / len(lh), sum(le) / len(le)
    sxx = sum((a - mh) ** 2 for a in lh)
    return sum((a - mh) * (b - me) for a, b in zip(lh, le)) / sxx if sxx > 0 else None


def launch(c):
    d = {}
    p = os.path.join(c, ".leia_launch")
    if os.path.isfile(p):
        for line in open(p):
            if "=" in line:
                k, v = line.strip().split("=", 1)
                d[k] = v
    return d


_STEP = re.compile(r"timeStepTime = (\S+) s")
_EXEC = re.compile(r"^ExecutionTime = (\S+) s", re.M)
_STENCIL = re.compile(r"davof: ringArea stencil .*?: (\d+) dropped, (\d+) kept d, (\d+) fit fallbacks, (\d+) remote rings")
_MODEL = re.compile(r"normalModel (\w+) \(.*?\): (\d+) slivers, (\d+) replaced(?:, (\d+) kept m)?(?:, (\d+) fit fallbacks)?")


def measure(c):
    tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) else []
    rows = [r for r in rows if f(r.get("TIME")) is not None and f(r.get("TIME")) > 0]
    rows.sort(key=lambda r: f(r.get("TIME")))
    last = rows[-1] if rows else {}
    col = lambda k: [f(r.get(k)) for r in rows if f(r.get(k)) is not None]
    logf = os.path.join(c, "log.davofAdvectorFoam")
    log = open(logf).read() if os.path.isfile(logf) else ""
    steps = [f(x) for x in _STEP.findall(log)]
    steps = [x for x in steps if x is not None]
    if not steps:
        ex = [f(x) for x in _EXEC.findall(log)]
        steps = [b - a for a, b in zip(ex, ex[1:]) if a is not None and b is not None]
    dropped = [int(m[0]) for m in _STENCIL.findall(log)]
    fitfb = [int(m[2]) for m in _STENCIL.findall(log)]
    remote = [int(m[3]) for m in _STENCIL.findall(log)]
    ringfit_fb = [int(m[4] or 0) for m in _MODEL.findall(log) if m[0] == "ringFit"]
    la = launch(c)
    wall = None
    if la.get("start") and la.get("end"):
        wall = float(la["end"]) - float(la["start"])
    L1a, L1n = col("L1_alphaf"), col("LNormalDiff1")
    g_early = g_late = None
    if len(L1a) >= 6 and L1a[0] > 0 and L1a[5] > 0:
        g_early = (L1a[5] / L1a[0]) ** (1.0 / 5)
    if len(L1a) >= 4:
        mid = len(L1a) // 2
        if L1a[mid - 1] > 0 and L1a[-1] > 0 and len(L1a) - mid > 0:
            g_late = (L1a[-1] / L1a[mid - 1]) ** (1.0 / (len(L1a) - mid))
    nInt = col("nInterfaceAdv") or col("nInterface")
    return dict(
        study=os.path.basename(os.path.dirname(c)), case=os.path.basename(c),
        variant=tok.get("TPF_NORMAL_STENCIL", "production"), row=tok.get("TPF_ROW"),
        arm=tok.get("TPF_VELOCITY_ARM"), N=int(float(tok["N_CELLS"])), Co=tok.get("TARGET_ADVECTIVE_CO"),
        rep=tok.get("TPF_REPEAT", tok.get("TPF_REPETITION", "1")), np=la.get("np", "1"),
        steps=len(rows), status=last.get("status", "missing"), abortStep=last.get("abortStep", ""),
        L1a_end=f(last.get("L1_alphaf")), L1n_end=f(last.get("LNormalDiff1")),
        L1x_end=f(last.get("LCentreDiff1")), Eshape_end=f(last.get("Eshape")),
        L1a_max=max(L1a) if L1a else None, L1n_max=max(L1n) if L1n else None,
        L1a_first=L1a[0] if L1a else None, g_early=g_early, g_late=g_late,
        minAlpha=min(col("minAlphaBound")) if col("minAlphaBound") else None,
        maxExcess=max(col("maxAlphaBound")) if col("maxAlphaBound") else None,
        devAbsMax=max(abs(x) for x in col("devPercent")) if col("devPercent") else None,
        nInterface_end=nInt[-1] if nInt else None, nInterface_max=max(nInt) if nInt else None,
        wall_s=wall, step_med_s=statistics.median(steps) if steps else None,
        dropped_max=max(dropped) if dropped else 0, dropped_sum=sum(dropped) if dropped else 0,
        fitFallback_sum=(sum(fitfb) if fitfb else 0) + (sum(ringfit_fb) if ringfit_fb else 0),
        remote_max=max(remote) if remote else 0,
        L1a_series=L1a, host=la.get("host", ""),
    )


def cases(study_dir):
    out = []
    for c in sorted(glob.glob(os.path.join(study_dir, "*_0*"))):
        if "." in os.path.basename(c) or not os.path.isfile(os.path.join(c, "case_params.json")):
            continue
        out.append(measure(c))
    return out


def fmt(v, p="%.3e"):
    if v is None:
        return "-"
    if isinstance(v, float):
        return p % v
    return str(v)


def med(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def spread(vals):
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None
    m = statistics.median(vals)
    return (max(vals) - min(vals)) / m if m else None


def md(headers, lines):
    s = "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n"
    for l in lines:
        s += "| " + " | ".join(l) + " |\n"
    return s


def write_csv(path, recs, keys):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            w.writerow({k: (fmt(r.get(k), "%.6e") if isinstance(r.get(k), float) else r.get(k)) for k in keys})


KEYS = ["study", "case", "variant", "row", "arm", "N", "Co", "rep", "np", "steps", "status", "abortStep",
        "L1a_end", "L1n_end", "L1x_end", "Eshape_end", "L1a_max", "L1n_max", "L1a_first", "g_early", "g_late",
        "minAlpha", "maxExcess", "devAbsMax", "nInterface_end", "nInterface_max", "wall_s", "step_med_s",
        "dropped_max", "dropped_sum", "fitFallback_sum", "remote_max", "host"]


def identical(a, b):
    """Whether two cases' per-step L1(alphaf) and L1(n) series agree to the digit."""
    if not a["L1a_series"] or len(a["L1a_series"]) != len(b["L1a_series"]):
        return False
    return all(x == y for x, y in zip(a["L1a_series"], b["L1a_series"]))


def test_plane(st, out, gate_tol=1e-11, port_plane=None):
    text = ""
    gate = cases(os.path.join(st, "stencil_plane_gate"))
    if gate:
        write_csv(os.path.join(out, "stencil_plane_gate.csv"), gate, KEYS)
        lines = []
        prod = {(r["N"], r["Co"], r["arm"]): r for r in gate if r["variant"] == "production"}
        for v in VARIANTS:
            sel = [r for r in gate if r["variant"] == v]
            if not sel:
                continue
            npass = sum(1 for r in sel if r["L1a_max"] is not None and r["L1n_max"] is not None
                        and r["L1a_max"] < gate_tol and r["L1n_max"] < gate_tol)
            ncomp = sum(1 for r in sel if r["status"] == "completed")
            ident = sum(1 for r in sel if (r["N"], r["Co"], r["arm"]) in prod
                        and identical(r, prod[(r["N"], r["Co"], r["arm"])]))
            ratio = med([r["step_med_s"] / prod[(r["N"], r["Co"], r["arm"])]["step_med_s"]
                         for r in sel if (r["N"], r["Co"], r["arm"]) in prod and r["step_med_s"]
                         and prod[(r["N"], r["Co"], r["arm"])]["step_med_s"]])
            lines.append([v, "%d/%d" % (ncomp, len(sel)), "%d/%d" % (npass, len(sel)),
                          fmt(max(r["L1a_max"] for r in sel if r["L1a_max"] is not None) if any(r["L1a_max"] for r in sel) else None),
                          fmt(max(r["L1n_max"] for r in sel if r["L1n_max"] is not None) if any(r["L1n_max"] for r in sel) else None),
                          fmt(min(r["minAlpha"] for r in sel if r["minAlpha"] is not None) if any(r["minAlpha"] is not None for r in sel) else None),
                          fmt(max(r["maxExcess"] for r in sel if r["maxExcess"] is not None) if any(r["maxExcess"] is not None for r in sel) else None),
                          fmt(max(r["devAbsMax"] for r in sel if r["devAbsMax"] is not None) if any(r["devAbsMax"] is not None for r in sel) else None),
                          str(sum(r["dropped_sum"] for r in sel)), str(sum(r["fitFallback_sum"] for r in sel)),
                          "%d/%d" % (ident, len(sel)), fmt(ratio, "%.2f")])
        text += "\n### Test 1a: the plane gate (72 cases per variant: N 40/80, Co 0.2/0.5, both arms)\n\n"
        text += md(["variant", "completed", "gate < 1e-11", "max L1(alphaf)", "max L1(n)", "min alpha", "max excess",
                    "max |dev| %", "cells dropped", "fit fallbacks", "= production", "step / production"], lines)
        if port_plane:
            port = [r for r in cases(port_plane) if r["row"] == "davof-planeSurf-ringArea"]
            port = {(r["N"], r["Co"], r["arm"], json.load(open(os.path.join(port_plane, r["case"], "case_params.json")))["tokens"].get("TPF_ENRIGHT_SETTING")): r for r in port}
            n_id = n_tot = 0
            for k, r in prod.items():
                pr = port.get((k[0], k[1], k[2], "recommended"))
                if pr is None:
                    continue
                n_tot += 1
                n_id += identical(r, pr)
            text += "\nG-0, the production rows against the port's record (port_enrightKeys_plane_leia, recommended): %d of %d shared cells identical to the digit at every step.\n" % (n_id, n_tot)
    growth = cases(os.path.join(st, "stencil_plane_growth"))
    if growth:
        write_csv(os.path.join(out, "stencil_plane_growth.csv"), growth, KEYS)
        cells = sorted({(r["N"], r["Co"], r["arm"]) for r in growth})
        lines = []
        for v in VARIANTS:
            sel = {(r["N"], r["Co"], r["arm"]): r for r in growth if r["variant"] == v}
            if not sel:
                continue
            row = [v]
            for k in cells:
                r = sel.get(k)
                row.append("-" if r is None else "%s / %s / %s / %s" % (fmt(r["g_early"], "%.2f"), fmt(r["g_late"], "%.2f"), fmt(r["L1a_first"], "%.1e"), fmt(r["L1a_end"], "%.1e")))
            lines.append(row)
        text += "\n### Test 1b: the seeded perturbation (1e-8 on the partial faces): g_early / g_late / L1(alphaf) after step 1 / at the end\n\n"
        text += md(["variant"] + ["N %d Co %s %s" % k for k in cells], lines)
        first = med([r["L1a_first"] for r in growth])
        text += "\nL1(alphaf) at step 1 (median over the cases): %s; steps per cell: %s.\n" % (
            fmt(first), ", ".join("%s: %d" % (k, med([r["steps"] for r in growth if (r["N"], r["Co"], r["arm"]) == k])) for k in cells))
    return text


def test_sphere(st, out):
    recs = []
    for N in (32, 64, 128):
        recs += cases(os.path.join(st, "stencil_sphereN%d" % N))
    if not recs:
        return ""
    write_csv(os.path.join(out, "stencil_sphere.csv"), recs, KEYS)
    text = "\n### Test 2: the translated sphere (Co 0.2, one period; N 32/64/128 on 1/8/64 ranks; three repetitions)\n\n"
    lines = []
    prod = {}
    for v in VARIANTS:
        for N in (32, 64, 128):
            sel = [r for r in recs if r["variant"] == v and r["N"] == N and r["status"] == "completed"]
            if sel and v == "production":
                prod[N] = sel
    for v in VARIANTS:
        byN = {}
        for N in (32, 64, 128):
            sel = [r for r in recs if r["variant"] == v and r["N"] == N]
            ok = [r for r in sel if r["status"] == "completed"]
            byN[N] = (sel, ok)
        if not any(byN[N][0] for N in byN):
            continue
        ns = [N for N in (32, 64, 128) if byN[N][1]]
        L1n = {N: med([r["L1n_end"] for r in byN[N][1]]) for N in ns}
        p_n = order(ns, [L1n[N] for N in ns])
        L1a = {N: med([r["L1a_end"] for r in byN[N][1]]) for N in ns}
        p_a = order(ns, [L1a[N] for N in ns])
        L1x = {N: med([r["L1x_end"] for r in byN[N][1]]) for N in ns}
        p_x = order(ns, [L1x[N] for N in ns])
        Es = {N: med([r["Eshape_end"] for r in byN[N][1]]) for N in ns}
        p_e = order(ns, [Es[N] for N in ns])
        walls = {N: med([r["wall_s"] for r in byN[N][1]]) for N in ns}
        sp = {N: spread([r["wall_s"] for r in byN[N][1]]) for N in ns}
        ratio = {N: (med([r["step_med_s"] for r in byN[N][1]]) / med([r["step_med_s"] for r in prod[N]]))
                 if N in prod and med([r["step_med_s"] for r in prod[N]]) else None for N in ns}
        comp = "%d/%d" % (sum(len(byN[N][1]) for N in byN), sum(len(byN[N][0]) for N in byN))
        minA = min([r["minAlpha"] for N in ns for r in byN[N][1] if r["minAlpha"] is not None] or [None]) if ns else None
        exc = max([r["maxExcess"] for N in ns for r in byN[N][1] if r["maxExcess"] is not None] or [None]) if ns else None
        dev = max([r["devAbsMax"] for N in ns for r in byN[N][1] if r["devAbsMax"] is not None] or [None]) if ns else None
        lines.append([v, comp] + [fmt(L1n.get(N)) for N in (32, 64, 128)] + [fmt(p_n, "%.2f"),
                      fmt(L1a.get(128)), fmt(p_a, "%.2f"), fmt(L1x.get(128)), fmt(p_x, "%.2f"), fmt(Es.get(128)), fmt(p_e, "%.2f"),
                      fmt(minA), fmt(exc), fmt(dev)]
                     + ["%s (%s)" % (fmt(walls.get(N), "%.0f"), fmt(sp.get(N), "%.2f")) for N in (32, 64, 128)]
                     + [fmt(ratio.get(N), "%.2f") for N in (32, 64, 128)])
    text += md(["variant", "completed", "L1(n) N32", "L1(n) N64", "L1(n) N128", "p(n)", "L1(alphaf) N128", "p(alphaf)",
                "L1(x) N128", "p(x)", "E_shape N128", "p(E)", "min alpha", "max excess", "max |dev| %",
                "wall N32 s (spread)", "wall N64 s", "wall N128 s", "step/prod N32", "N64", "N128"], lines)
    # identity of the repetitions and of the angle filters
    idl = []
    for v in VARIANTS:
        for N in (32, 64, 128):
            sel = [r for r in recs if r["variant"] == v and r["N"] == N and r["status"] == "completed"]
            if len(sel) >= 2:
                same = all(identical(sel[0], r) for r in sel[1:])
                pr = [r for r in recs if r["variant"] == "production" and r["N"] == N and r["status"] == "completed"]
                eqp = identical(sel[0], pr[0]) if pr else None
                idl.append("%s N %d: reps %s, = production %s" % (v, N, "identical" if same else "DIFFER", eqp))
    text += "\nRepetitions and identity: " + "; ".join(idl) + ".\n"
    return text


def spearman(xs, ys):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for k, i in enumerate(s):
            r[i] = k
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    if n < 3:
        return None
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
    return 1 - 6 * d2 / (n * (n * n - 1))


def test_flow(st, out, study, title, Ns, rep_key="TPF_REPEAT"):
    recs = cases(os.path.join(st, study))
    if not recs:
        return ""
    write_csv(os.path.join(out, study + ".csv"), recs, KEYS)
    text = "\n### %s\n\n" % title
    lines = []
    prod = {N: [r for r in recs if r["variant"] == "production" and r["N"] == N and r["status"] == "completed"] for N in Ns}
    tab = {}
    for v in VARIANTS:
        for N in Ns:
            sel = [r for r in recs if r["variant"] == v and r["N"] == N]
            ok = [r for r in sel if r["status"] == "completed"]
            if not sel:
                continue
            st_ = "%d/%d" % (len(ok), len(sel))
            ab = ",".join(sorted({str(r["abortStep"]) for r in sel if r["status"] != "completed" and r["abortStep"]}))
            ratio = (med([r["step_med_s"] for r in ok]) / med([r["step_med_s"] for r in prod[N]])) if ok and prod[N] and med([r["step_med_s"] for r in prod[N]]) else None
            tab[(v, N)] = (med([r["L1n_end"] for r in ok]), med([r["Eshape_end"] for r in ok]), med([r["L1x_end"] for r in ok]))
            lines.append([v, str(N), st_ + ((" abort " + ab) if ab else ""),
                          fmt(min([r["minAlpha"] for r in sel if r["minAlpha"] is not None] or [None])),
                          fmt(max([r["maxExcess"] for r in sel if r["maxExcess"] is not None] or [None])),
                          fmt(max([r["devAbsMax"] for r in sel if r["devAbsMax"] is not None] or [None])),
                          fmt(med([r["Eshape_end"] for r in ok])), fmt(med([r["L1x_end"] for r in ok])),
                          fmt(med([r["L1n_end"] for r in ok])), fmt(med([r["L1a_end"] for r in ok])),
                          fmt(med([r["nInterface_end"] for r in ok]), "%.0f"), fmt(med([r["nInterface_max"] for r in sel]), "%.0f"),
                          str(max(r["dropped_max"] for r in sel)), str(sum(r["fitFallback_sum"] for r in sel)),
                          "%s (%s)" % (fmt(med([r["wall_s"] for r in ok]), "%.0f"), fmt(spread([r["wall_s"] for r in ok]), "%.2f")),
                          fmt(ratio, "%.2f")])
    text += md(["variant", "N", "completed", "min alpha", "max excess", "max |dev| %", "E_shape", "L1(x)", "L1(n)", "L1(alphaf)",
                "interface cells end", "max", "cells dropped (max/step)", "fit fallbacks", "wall s (spread)", "step/prod"], lines)
    for N in Ns:
        vs = [v for v in VARIANTS if (v, N) in tab and all(x is not None for x in tab[(v, N)])]
        if len(vs) >= 3:
            rho_e = spearman([tab[(v, N)][0] for v in vs], [tab[(v, N)][1] for v in vs])
            rho_x = spearman([tab[(v, N)][0] for v in vs], [tab[(v, N)][2] for v in vs])
            text += "\nSpearman rank correlation at N %d over %d variants: L1(n) vs E_shape %s, L1(n) vs L1(x) %s.\n" % (N, len(vs), fmt(rho_e, "%.2f"), fmt(rho_x, "%.2f"))
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", default="studies/davof")
    ap.add_argument("--out", default=".")
    ap.add_argument("--port-plane", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    text = test_plane(a.studies, a.out, port_plane=a.port_plane or None)
    text += test_sphere(a.studies, a.out)
    text += test_flow(a.studies, a.out, "stencil_rotation", "Test 3a: the rotation (one turn, Co 0.2, N 32/64 serial, three repetitions)", (32, 64))
    text += test_flow(a.studies, a.out, "stencil_deformation", "Test 3b: the 3D deformation (N 32, OpenFOAM Co 0.5, T = 3, three repetitions)", (32,), "TPF_REPETITION")
    print(text)
    with open(os.path.join(a.out, "stencil_tables.md"), "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    main()
