#!/usr/bin/env python3
"""Tables of the signed flux error study (leia davof-normal-order, 2026-11-14;
docs/davof/REPORT-swept.md "Signed flux error of the biased tilt").

    python3 workflow/scripts/sgn_tables.py studies/davof

Reads studies/davof/nor_sgn_N{32,64,128}_{base,bia,alt}/*/advectorVoF.csv and
case_params.json, prints the tables as markdown and writes sgn_cases.csv.
"""
import csv, glob, json, math, os, statistics, sys

root = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "sgn_cases.csv"


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


cases = []
for d in sorted(glob.glob(os.path.join(root, "nor_sgn_N*_*", "*"))):
    cp = os.path.join(d, "case_params.json")
    cs = os.path.join(d, "advectorVoF.csv")
    if not (os.path.exists(cp) and os.path.exists(cs)):
        continue
    t = json.load(open(cp))["tokens"]
    rows = list(csv.DictReader(open(cs)))

    def col(k):
        return [v for v in (fl(r.get(k)) for r in rows) if v is not None]

    def mean(k):
        x = col(k)
        return sum(x) / len(x) if x else None

    def med(k):
        x = col(k)
        return statistics.median(x) if x else None
    row = t["TPF_ROW"].replace("davof-planeSurf-", "")
    mode = t.get("TPF_TILT_MODE", "") if row == "ringTilt" else ""
    c = fl(t.get("TPF_TILT_C")) if row == "ringTilt" else 0.0
    key = {"exact": "exact", "ringArea": "ring"}.get(row, mode)
    last = rows[-1] if rows else {}
    cases.append(dict(
        N=int(t["N_CELLS"]), row=row, key=key, c=c, steps=len(col("fluxErrSigned")),
        status=last.get("status", ""),
        s_or=mean("fluxErrSigned"), s_flow=mean("fsgnFlow"), s_wet=mean("fsgnWet"),
        s_or2=mean("fsgnOrient"), fl1=mean("fluxErrL1"),
        ekL1=mean("ekL1"), ekL1med=med("ekL1"), ekL2=mean("ekL2"), ekInf=mean("ekInf"),
        ekSigned=mean("ekSigned"), ekSumMax=max([abs(v) for v in col("ekSum")] or [0]),
        rawAbs=mean("ekRawAbs"),
        sx=mean("sgnX"), sy=mean("sgnY"), sz=mean("sgnZ"),
        Eshape=fl(last.get("Eshape")), posB=med("posBefore"), posA=med("posAfter"),
        wall=fl(last.get("clockTime"))))

with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cases[0].keys()))
    w.writeheader()
    w.writerows(cases)


def g(x, f="%.3e"):
    return "-" if x is None else f % x


def find(N, key, c=None):
    for r in cases:
        if r["N"] == N and r["key"] == key and (c is None or abs(r["c"] - c) < 1e-12):
            return r


order = {"exact": 0, "ring": 1, "biased": 2, "alternating": 3}
print("## All cases (means over steps)\n")
print("| N | row | c | steps | status | s (orientation) | s (flow) | s (wetter) | L1 face | e_k L1 | e_k L2 | e_k Linf | e_k signed | max abs sum(e V) | sum abs raw | E_shape end | pos before | pos after | wall s |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(cases, key=lambda r: (r["N"], order[r["key"]], r["c"])):
    print("| %d | %s | %g | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        r["N"], r["key"], r["c"], r["steps"], r["status"], g(r["s_or"]), g(r["s_flow"]), g(r["s_wet"]),
        g(r["fl1"]), g(r["ekL1"]), g(r["ekL2"]), g(r["ekInf"]), g(r["ekSigned"]), g(r["ekSumMax"]),
        g(r["rawAbs"]), g(r["Eshape"]), g(r["posB"]), g(r["posA"]), g(r["wall"], "%.0f")))

# H1: s(c) at N 64, biased with c = 0 the ring row
print("\n## H1: s(c) at N 64, biased tilt (c = 0: the ring row), least-squares line\n")
pts = []
for r in cases:
    if r["N"] == 64 and r["key"] in ("biased", "ring") and r["s_or"] is not None:
        pts.append((r["c"], r["s_or"]))
pts.sort()
for c, s in pts:
    print("- c = %g: s = %s" % (c, g(s)))
if len(pts) >= 3:
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    s1 = sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx
    s0 = my - s1 * mx
    ss_res = sum((p[1] - (s0 + s1 * p[0])) ** 2 for p in pts)
    ss_tot = sum((p[1] - my) ** 2 for p in pts)
    print("\nfit s = s0 + s1 c: s0 = %s, s1 = %s, R^2 = %.4f, zero crossing c* = %s, "
          "max |residual| = %s; s(-0.25) measured %s, 2 s0 = %s" % (
              g(s0), g(s1), 1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
              g(-s0 / s1 if s1 else None, "%.3f"), g(max(abs(p[1] - (s0 + s1 * p[0])) for p in pts)),
              g(dict(pts).get(-0.25)), g(2 * s0)))
print("\nalternating, N 64: " + "; ".join("c = %g: s = %s" % (r["c"], g(r["s_or"]))
                                          for r in sorted(cases, key=lambda r: r["c"])
                                          if r["N"] == 64 and r["key"] == "alternating"))
print("exact N 64: s = %s" % g(find(64, "exact")["s_or"] if find(64, "exact") else None))

# H2: ratios
print("\n## H2: per-cell net error against E_shape (N 64 and all N)\n")
print("| N | row | c | e_k L1 | e_k L1 / ring | E_shape end / ring | sum abs raw | e_k L1 / (face L1) |")
print("|---|---|---|---|---|---|---|---|")
for N in (32, 64, 128):
    ring = find(N, "ring")
    for r in sorted([r for r in cases if r["N"] == N], key=lambda r: (order[r["key"]], r["c"])):
        if r["ekL1"] is None:
            continue
        print("| %d | %s | %g | %s | %s | %s | %s | %s |" % (
            N, r["key"], r["c"], g(r["ekL1"]),
            g(r["ekL1"] / ring["ekL1"] if ring and ring["ekL1"] else None, "%.2f"),
            g(r["Eshape"] / ring["Eshape"] if ring and ring["Eshape"] and r["Eshape"] else None, "%.2f"),
            g(r["rawAbs"]), g(r["ekL1"] / r["fl1"] if r["fl1"] else None, "%.3f")))

# split
print("\n## x/y/z split of the signed mean (orientation convention), N 64\n")
print("| row | c | x | y | z |")
print("|---|---|---|---|---|")
for r in sorted([r for r in cases if r["N"] == 64], key=lambda r: (order[r["key"]], r["c"])):
    print("| %s | %g | %s | %s | %s |" % (r["key"], r["c"], g(r["sx"]), g(r["sy"]), g(r["sz"])))

# orders
print("\n## Orders over N 32/64/128 (least squares of log y against log N, p = -slope)\n")
print("| row | c | p of s (|s|; sign) | p of e_k L1 | p of face L1 | p of E_shape |")
print("|---|---|---|---|---|---|")


def lsorder(vals):
    pts = [(math.log(N), math.log(abs(v))) for N, v in vals if v not in (None, 0)]
    if len(pts) < 3:
        return None
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    return -sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)


for key, c in (("exact", 0.0), ("ring", 0.0), ("biased", -0.25), ("biased", 0.25), ("alternating", 0.25)):
    rs = [(N, find(N, key, c if key in ("biased", "alternating") else None)) for N in (32, 64, 128)]
    if any(r is None for _, r in rs):
        continue
    signs = "".join("+" if r["s_or"] >= 0 else "-" for _, r in rs if r["s_or"] is not None)
    print("| %s | %g | %s (%s) | %s | %s | %s |" % (
        key, c, g(lsorder([(N, r["s_or"]) for N, r in rs]), "%.2f"), signs,
        g(lsorder([(N, r["ekL1"]) for N, r in rs]), "%.2f"),
        g(lsorder([(N, r["fl1"]) for N, r in rs]), "%.2f"),
        g(lsorder([(N, r["Eshape"]) for N, r in rs]), "%.2f")))
