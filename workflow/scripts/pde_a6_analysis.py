#!/usr/bin/env python3
"""pde_a6_analysis.py: the read-out of step A, round 2 (the model's own extension u(Y), per layer).

    pde_a6_analysis.py --old <clone of round 1 (leia-pde-a5)> --new <clone a6> [--out <dir>]

Reads leiaTestVelocityExtensionLayers.csv (per topological layer l = 0..8 and l = 0..3) of
aStatic2Dr2, aStatic2Dr2Ellipse and aStatic2Dr2Par, and checks the pre-registered read-outs R1 to R9
of the config headers. L1 and L2 only (never L_inf); orders between consecutive N. Writes
<out>/a6_summary.txt. Stdlib only.
"""
import argparse, csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pde_a2_analysis import tok  # noqa: E402
from pde_a3_analysis import cmp, counters  # noqa: E402

LAY = [f"L{l}" for l in range(9)] + ["L0TO3"]


def layers_row(case):
    f = os.path.join(case, "leiaTestVelocityExtensionLayers.csv")
    return {k: float(v) for k, v in list(csv.DictReader(open(f)))[0].items()} \
        if os.path.isfile(f) and os.path.getsize(f) else {}


def model(t):
    m = t["VELOCITY_EXTENSION"]
    return f"{m}/{t['PB_DISTANCE']}" if m == "projectionBlend" else f"{m}/{t['CP_STEP']}/{t['FADE_MODE']}"


def orders(v):
    return [math.log(e1 / e2) / math.log(n2 / n1) if e1 > 0 and e2 > 0 else float("nan")
            for (n1, e1), (n2, e2) in zip(v, v[1:])]


def load(study):
    """(model, input, N) -> (case, row); duplicates of one arm (an axis the model does not read) are
    checked for byte identity and kept once."""
    out, dup = {}, []
    for c in sorted(glob.glob(os.path.join(study, "2Dvortex_[0-9]*"))):
        t = tok(c)
        key = (model(t), t["VE_TEST_INPUT"], int(t["N_CELLS"]))
        if key in out:
            v, _ = cmp(os.path.join(out[key][0], "leiaTestVelocityExtensionLayers.csv"),
                       os.path.join(c, "leiaTestVelocityExtensionLayers.csv"))
            dup.append((key, v))
            continue
        out[key] = (c, layers_row(c))
    return out, dup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "a6")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda root, st: os.path.join(root, "studies", "pdeLevelSet", st)
    verdict = {}

    def judge(tag, ok):
        verdict.setdefault(tag, True)
        verdict[tag] = verdict[tag] and ok

    data = {}
    for st in ("aStatic2Dr2", "aStatic2Dr2Ellipse", "aStatic2Dr2Par"):
        data[st], dup = load(S(a.new, st))
        bad = [d for d in dup if d[1] != "PASS"]
        say(f"{st}: {len(data[st])} arms, {len(dup)} duplicate cases (an axis the model does not read), "
            f"{len(dup) - len(bad)} identical" + (f"; NOT identical: {bad}" if bad else ""))

    say()
    say("R6 / R6e INERTNESS: leiaTestVelocityExtensionExact.csv of every closestPoint case against round 1:")
    for st, old in (("aStatic2Dr2", "aStatic2Dq"), ("aStatic2Dr2Ellipse", "aStatic2Dellipse")):
        ref = {}
        for c in glob.glob(os.path.join(S(a.old, old), "2Dvortex_[0-9]*")):
            t = tok(c)
            if (t["VELOCITY_EXTENSION"], t["CP_STEP"], t["FADE_MODE"]) == ("closestPoint", "newton", "layer"):
                ref[(t["VE_TEST_INPUT"], int(t["N_CELLS"]))] = c
        n = ok = 0
        for (m, inp, N), (c, _) in sorted(data[st].items()):
            if not m.startswith("closestPoint"):
                continue
            f = "leiaTestVelocityExtensionExact.csv"
            v, d = cmp(os.path.join(ref.get((inp, N), "/nonexistent"), f), os.path.join(c, f))
            n += 1; ok += v == "PASS"
            if v != "PASS":
                say(f"  {st} {inp} N {N}: {v} {' | '.join(d[-2:])[:160]}")
        say(f"  {st} against {old}: {ok} of {n} identical")
        judge("R6", n > 0 and ok == n)

    say()
    say("Wrong-reason checks: the input applied; layers populated; delta = 0 in layers 0..3:")
    for st in data:
        miss = 0
        for (m, inp, N), (c, r) in data[st].items():
            if inp != "distance":
                txt = "".join(open(f, errors="ignore").read() for f in glob.glob(os.path.join(c, "log.leiaTestVelocityExtension")))
                miss += "psi replaced" not in txt
            if not r or any(r.get(f"N_L{l}", 0) <= 0 for l in range(4)):
                say(f"  {st} {m} {inp} N {N}: EMPTY layer or no CSV")
                judge("checks", False)
            if m.startswith("projectionBlend") and r:
                y, b = r["E_EXTY_L1_L0TO3"], r["E_EXTB_L1_L0TO3"]
                if abs(y - b) > 1e-12 * max(abs(y), 1e-300):
                    say(f"  {st} {m} {inp} N {N}: E_EXTB != E_EXTY on l <= 3 ({b:.6e} vs {y:.6e}): delta > 0 there")
                    judge("checks", False)
        say(f"  {st}: input missing in {miss} non-distance cases")
        judge("checks", miss == 0)

    def series(st, m, inp, col):
        return [(N, r[col]) for (mm, ii, N), (_, r) in sorted(data[st].items(), key=lambda kv: kv[0][2])
                if mm == m and ii == inp and col in r]

    def show(st, m, inp, col, tag=None, lo=None, hi=None):
        v = series(st, m, inp, col)
        o = orders(v)
        txt = "  ".join(f"N{n}:{e:.2e}" for n, e in v) + "   orders " + " / ".join(f"{x:.2f}" for x in o)
        res = ""
        if tag:
            # pre-registered over N 128-256-512: the orders of the pairs with N1 >= 128
            oo = [x for (n1, _), x in zip(v, o) if n1 >= 128]
            ok = bool(oo) and all((lo is None or x >= lo) and (hi is None or x <= hi) for x in oo)
            if v and all(e <= 1e-13 for _, e in v):
                ok = True  # identically zero (round-off) at every N: nothing to converge
            judge(tag, ok); res = f"   [{tag} {'PASS' if ok else 'FAIL'}]"
        say(f"    {col:20s} {txt}{res}")

    for st in ("aStatic2Dr2", "aStatic2Dr2Ellipse"):
        say()
        say(f"===== {st}: errors (L1, L2) and orders over N, layers 0..3 and per layer")
        sfx = "e" if "Ellipse" in st else ""
        for m, inp in sorted({(k[0], k[1]) for k in data[st]}):
            say(f"  == {m}  input {inp}")
            if m == "projectionBlend/first":
                show(st, m, inp, "E_EXTY_L1_L0TO3", "R1" + sfx, 1.8)
                show(st, m, inp, "E_FLUXY_L1_L0TO3", "R1" + sfx, 1.8)
                show(st, m, inp, "E_EXTY_L2_L0TO3")
                show(st, m, inp, "E_EXT_L1_L0TO3")
                for l in range(9):
                    show(st, m, inp, f"E_EXTB_L1_L{l}", "R2" + sfx, 1.8)
            elif m.startswith("closestPoint"):
                show(st, m, inp, "E_EXT_L1_L0TO3", "R3" + sfx, 1.8)
                show(st, m, inp, "E_EXTY_L1_L0TO3", "R3" + sfx, 1.8)
                show(st, m, inp, "E_EXTY_L2_L0TO3")
                show(st, m, inp, "E_FLUXY_L1_L0TO3")
            else:
                for col in ("E_EXTY_L1_L0TO3", "E_EXT_L1_L0TO3", "E_FLUXY_L1_L0TO3", "E_EXTY_L2_L0TO3"):
                    show(st, m, inp, col)
            if inp != "distance" and st == "aStatic2Dr2" and m.startswith("closestPoint"):
                show(st, m, inp, "D_YPI_L1_L0TO3", "R5", 0.8, 1.2)
            # per layer at the finest N
            v = {N: r for (mm, ii, N), (_, r) in data[st].items() if mm == m and ii == inp}
            if v:
                N = max(v); r = v[N]
                say(f"    per layer at N {N}:  l: E_EXTY_L1 / E_EXTB_L1 / E_FLUXY_L1 / E_RAWY_L1 / D_YPI / cells")
                for l in range(9):
                    t = f"L{l}"
                    say(f"      {l}: {r[f'E_EXTY_L1_{t}']:.2e} / {r[f'E_EXTB_L1_{t}']:.2e} / {r[f'E_FLUXY_L1_{t}']:.2e}"
                        f" / {r[f'E_RAWY_L1_{t}']:.2e} / {r[f'D_YPI_L1_{t}']:.2e} / {int(r[f'N_{t}'])}")
                if m == "projectionBlend/first" and inp == "distance" and st == "aStatic2Dr2":
                    ratio = r["E_EXTY_L1_L3"] / r["E_EXTY_L1_L0"] if r["E_EXTY_L1_L0"] > 0 else float("nan")
                    judge("R4", ratio >= 4)
                    say(f"    R4 E_EXTY_L1_L3 / E_EXTY_L1_L0 at N {N} = {ratio:.2f} [{'PASS' if ratio >= 4 else 'FAIL'}]")

    say()
    say("===== aStatic2Dr2Par against aStatic2Dr2 (np 4 against serial):")
    for (m, inp, N), (c, r) in sorted(data["aStatic2Dr2Par"].items()):
        s = data["aStatic2Dr2"].get((m, inp, N), (None, {}))[1]
        if not s or not r:
            say(f"  {m} {inp} N {N}: no serial or no parallel CSV"); judge("R7", False); continue
        same = all(r[f"N_{t}"] == s[f"N_{t}"] for t in LAY)
        judge("R7", same)
        cols = [k for k in r if k.endswith("_L0TO3") and k.startswith("E_")]
        rel = max(abs(r[k] - s[k]) / abs(s[k]) for k in cols if s[k] != 0)
        if m.startswith("projectionBlend"):
            judge("R8", rel <= 1e-8); tag = f"R8 {'PASS' if rel <= 1e-8 else 'FAIL'}"
        else:
            e = abs(r["E_EXT_L1_L0TO3"] - s["E_EXT_L1_L0TO3"]) / s["E_EXT_L1_L0TO3"]
            judge("R9", e <= 0.02); tag = f"R9 E_EXT_L1_L0TO3 {e:.2%} {'PASS' if e <= 0.02 else 'FAIL'}"
        fb, _, _ = counters(c)
        say(f"  {m:28s} {inp:15s} N {N:4d}: layer counts {'equal' if same else 'DIFFER'}, "
            f"max rel. diff of the E_*_L0TO3 columns {rel:.2e}  [{tag}]" + (f"  cells {fb}" if fb else ""))

    say()
    say("VERDICTS: " + ", ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in sorted(verdict.items())))
    open(os.path.join(out, "a6_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
