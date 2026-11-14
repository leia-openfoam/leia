#!/usr/bin/env python3
"""Analysis of the DIAG lines of the one-plane first-fault study (2026-11-15; docs/davof/REPORT-swept.md,
"One plane per face on curved interfaces: first fault"). Reads log.davofAdvectorFoam of a case and prints the
tables as markdown. Reading and printing only.

usage: diag1_an.py <log> [--label NAME]
"""
import re
import sys
from collections import OrderedDict

F = r"([-+0-9.eE]+|nan|inf|-inf)"
RE_STEP = re.compile(
    r"^DIAG step (\d+) t " + F + r" : unlimited update: cells outside \[-1e-10, 1\+1e-10\] (\d+), non-finite (\d+), "
    r"min alpha " + F + r", max alpha - 1 " + F + r", summed violation \(cell volumes\) " + F +
    r", cells with violation > 1e-3 (\d+) \| of the violating cells: with a no-plane cut face (\d+), "
    r"with a no-plane face \(cut or uncut\) (\d+), with an unprocessed face (\d+), all faces with planes (\d+), "
    r"of those with a face \|one - pieces\| > 0.1 V_k (\d+) \| no-plane faces uncut (\d+), cut (\d+) \| compare (\w+)")
RE_LIM = re.compile(
    r"^DIAG step (\d+) : limiter passes (\d+), faces changed (\d+), cells touched (\d+), volume moved \(cell volumes\) " + F)
RE_CMP = re.compile(
    r"^DIAG compare step (\d+) : .* over the (\d+) faces of interface cells: L1 \(mean\) " + F + r", Linf " + F +
    r", count > 0.01 (\d+), count > 0.1 (\d+) \(of those: donor plane (\d+), no plane cut (\d+), no plane uncut (\d+)\)")
RE_CELL = re.compile(
    r"^DIAG step (\d+) cell (\d+) at \(([^)]*)\) alpha\^n " + F + r" alpha after update " + F +
    r" hasPlane (\d) \|m\|/V\^\(2/3\) " + F + r"(?: outflow/V_k " + F + r" inflow/V_k " + F + r")?")
RE_FACE = re.compile(r"^    face (\d+)( \(boundary\))? flux out/V_k " + F + r" phi dt/V_k " + F + r" state (\S+) ringDepth " + F +
                     r" donor (-?\d+)(?: alpha\^n " + F + r" \|m\|/V\^\(2/3\) " + F + r")? donorIsCentroidCell (\d) prismWet " + F +
                     r" prismV/V_k " + F + r"(?: \|one-pieces\|/V_k " + F + r")?(?: other (-?\d+) alpha\^n " + F + r")?")
RE_TRAP = re.compile(r"^DIAG boundFlux .*")
RE_TIME = re.compile(r"^Time = ([0-9.eE+-]+)")


def parse(path):
    steps = OrderedDict()
    cmp_ = {}
    lim = {}
    cells = {}
    traps = []
    cur = None
    last_time = None
    abort = []
    with open(path, errors="replace") as fh:
        for line in fh:
            # parallel: a rank's partial line or its "[n] " prefix may precede the text
            idx = line.find("DIAG ")
            if idx > 0:
                line = line[idx:]
            m = RE_STEP.match(line)
            if m:
                g = m.groups()
                s = int(g[0])
                steps[s] = dict(t=float(g[1]), nViol=int(g[2]), nonFinite=int(g[3]), minA=float(g[4]), maxA1=float(g[5]),
                                sumViol=float(g[6]), nGt3=int(g[7]), cNoPlaneCut=int(g[8]), cNoPlane=int(g[9]),
                                cUnproc=int(g[10]), cAllPlanes=int(g[11]), cAllPlanesDiff=int(g[12]),
                                nF2=int(g[13]), nF3=int(g[14]), cmp=(g[15] == "yes"))
                continue
            m = RE_LIM.match(line)
            if m:
                lim[int(m.group(1))] = dict(passes=int(m.group(2)), faces=int(m.group(3)), cells=int(m.group(4)),
                                            vol=float(m.group(5)))
                continue
            m = RE_CMP.match(line)
            if m:
                g = m.groups()
                cmp_[int(g[0])] = dict(n=int(g[1]), L1=float(g[2]), Linf=float(g[3]), n01=int(g[4]), n1=int(g[5]),
                                       plane=int(g[6]), cut=int(g[7]), uncut=int(g[8]))
                continue
            m = RE_CELL.match(line)
            if m:
                g = m.groups()
                try:
                    cur = dict(step=int(g[0]), cell=int(g[1]), x=g[2], a0=float(g[3]), a1=float(g[4]), plane=int(g[5]),
                               m=float(g[6]), out=float(g[7] or "nan"), inn=float(g[8] or "nan"), faces=[])
                    cells.setdefault(cur["step"], []).append(cur)
                except ValueError:
                    cur = None   # a line glued to another rank's output (parallel)
                continue
            m = RE_FACE.match(line)
            if m and cur is not None:
                g = m.groups()
                try:
                    cur["faces"].append(dict(face=int(g[0]), flux=float(g[2]), phidt=float(g[3]), state=g[4], depth=float(g[5]),
                                             donor=int(g[6]), donorA=(float(g[7]) if g[7] else None),
                                             donorM=(float(g[8]) if g[8] else None), cen=int(g[9]), wet=float(g[10]),
                                             vtot=float(g[11]), dpiece=(float(g[12]) if g[12] else None),
                                             other=(int(g[13]) if g[13] else None),
                                             otherA=(float(g[14]) if g[14] else None)))
                except ValueError:
                    pass   # a line cut by another rank's output (parallel)
                continue
            if line.startswith("DIAG boundFlux"):
                traps.append(line.strip())
                continue
            m = RE_TIME.match(line)
            if m:
                last_time = m.group(1)
            if "ABORT" in line or "sigFpe" in line or "Floating point exception" in line:
                abort.append(line.strip()[:200])
    return steps, lim, cmp_, cells, traps, last_time, abort


def fmt(x):
    return f"{x:.3g}"


def main():
    path = sys.argv[1]
    label = sys.argv[sys.argv.index("--label") + 1] if "--label" in sys.argv else path
    steps, lim, cmp_, cells, traps, last_time, abort = parse(path)
    print(f"## {label}\n")
    if not steps:
        print("no DIAG lines\n")
        return
    last = max(steps)
    print(f"steps with a DIAG line: {len(steps)}; last step {last}, t = {steps[last]['t']:.6g}; "
          f"abort/trap lines: {abort[:3]}\n")

    def first(pred):
        for s, d in steps.items():
            if pred(s, d):
                return s
        return None

    sv = first(lambda s, d: d["nViol"] + d["nonFinite"] > 0)
    s_nc = first(lambda s, d: d["nF3"] > 0)
    s_ncv = first(lambda s, d: d["cNoPlaneCut"] > 0)
    s3 = first(lambda s, d: d["nGt3"] > 0)
    sl = first(lambda s, d: s in lim and lim[s]["vol"] > 1e-3)
    sl_any = first(lambda s, d: s in lim and lim[s]["vol"] > 0)
    print("### Timeline\n")
    print("| event | step | t |\n|---|---|---|")
    for name, s in [("first step with a pre-limiter violation outside [-1e-10, 1+1e-10]", sv),
                    ("first step with a no-plane cut face anywhere", s_nc),
                    ("first step with a violating cell that has a no-plane cut face", s_ncv),
                    ("first step with a violation > 1e-3", s3),
                    ("first step with any limiter volume moved", sl_any),
                    ("first step with limiter volume > 1e-3 cell volumes", sl),
                    ("last step printed (crash or abort or end)", last)]:
        print(f"| {name} | {s if s is not None else 'none'} | {steps[s]['t']:.6g} |" if s is not None else f"| {name} | none | |")
    print()
    if sv is not None:
        d = steps[sv]
        print(f"first violation step {sv}: {d['nViol']} cells, min alpha {fmt(d['minA'])}, max alpha - 1 {fmt(d['maxA1'])}, "
              f"summed violation {fmt(d['sumViol'])} cell volumes, no-plane faces uncut {d['nF2']} cut {d['nF3']}\n")

    # the first 50 steps with violations
    vs = [s for s, d in steps.items() if d["nViol"] + d["nonFinite"] > 0][:50]
    print("### The first 50 steps with violations (sums over those steps)\n")
    if vs:
        tot = lambda k: sum(steps[s][k] for s in vs)
        print(f"steps {vs[0]}..{vs[-1]} ({len(vs)} steps): violating cells {tot('nViol')}; "
              f"with a no-plane cut face {tot('cNoPlaneCut')}; with a no-plane face (cut or uncut) {tot('cNoPlane')}; "
              f"with an unprocessed face {tot('cUnproc')}; with all faces with planes {tot('cAllPlanes')}; "
              f"all faces with planes and a face |one - pieces| > 0.1 V_k {tot('cAllPlanesDiff')} "
              f"(compare steps among them: {sum(1 for s in vs if steps[s]['cmp'])}); "
              f"steps with max violation < 1e-3: {sum(1 for s in vs if steps[s]['nGt3'] == 0)}; "
              f"steps with a cut no-plane violating cell: {sum(1 for s in vs if steps[s]['cNoPlaneCut'] > 0)}\n")
        print("| step | cells out | min alpha | max alpha-1 | summed viol | cells > 1e-3 | cut no-plane | no-plane | unproc. | all planes | all planes, diff > 0.1 | limiter vol | compare Linf | compare n>0.1 |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for s in vs[:50]:
            d = steps[s]
            c = cmp_.get(s)
            print(f"| {s} | {d['nViol']} | {fmt(d['minA'])} | {fmt(d['maxA1'])} | {fmt(d['sumViol'])} | {d['nGt3']} | {d['cNoPlaneCut']} | "
                  f"{d['cNoPlane']} | {d['cUnproc']} | {d['cAllPlanes']} | {d['cAllPlanesDiff']} | "
                  f"{fmt(lim[s]['vol']) if s in lim else 'n/a'} | {fmt(c['Linf']) if c else 'n/a'} | {c['n1'] if c else 'n/a'} |")
        print()

    # the face classification of the first 20 violating cells at the first violation step and at the first > 1e-3 step
    for name, s in [("first violation step", sv), ("first step with a violation > 1e-3", s3), ("last step", last)]:
        if s is None or s not in cells:
            continue
        print(f"### The violating cells at the {name} (step {s}), first {len(cells[s])} by cell id\n")
        print("| cell | alpha^n | alpha new | has plane | outflow/V_k | inflow/V_k | faces: plane / noPlane uncut / noPlane CUT | donor = centroid cell | max prismWet | max flux out/V_k | max alpha^n of donors | face with max |one-pieces|/V_k |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for c in cells[s]:
            if not c["faces"]:
                continue
            st = {"plane": 0, "noPlane-uncut": 0, "noPlane-CUT": 0, "unprocessed": 0}
            for f in c["faces"]:
                st[f["state"]] += 1
            dp = [f["dpiece"] for f in c["faces"] if f["dpiece"] is not None]
            print(f"| {c['cell']} | {fmt(c['a0'])} | {fmt(c['a1'])} | {c['plane']} | {fmt(c['out'])} | {fmt(c['inn'])} | "
                  f"{st['plane']} / {st['noPlane-uncut']} / {st['noPlane-CUT']} (unproc. {st['unprocessed']}) | "
                  f"{sum(f['cen'] for f in c['faces'])} of {len(c['faces'])} | {fmt(max(f['wet'] for f in c['faces']))} | "
                  f"{fmt(max(f['flux'] for f in c['faces']))} | {fmt(max([f['donorA'] for f in c['faces'] if f['donorA'] is not None] or [0]))} | "
                  f"{fmt(max(dp)) if dp else 'n/a'} |")
        print()

    print("### Limiter per step (last 5 steps)\n")
    for s in sorted(lim)[-5:]:
        l = lim[s]
        print(f"- step {s}: passes {l['passes']}, faces changed {l['faces']}, cells touched {l['cells']}, volume moved {fmt(l['vol'])} cell volumes")
    print()
    print("### boundFlux division operands (denominator below 1e-12)\n")
    if traps:
        for t in traps[:6]:
            print("- " + t)
        print(f"- ({len(traps)} such lines)")
    else:
        print("none\n")
    print()
    print("### onePlane vs pieces flux difference |flux_onePlane - flux_pieces|/V_k over the faces of interface cells\n")
    print("| step | faces | L1 (mean) | Linf | count > 0.01 | count > 0.1 | of those: donor plane / no plane cut / no plane uncut |")
    print("|---|---|---|---|---|---|---|")
    want = set()
    for s in (sv, s3, last):
        if s is None:
            continue
        for k in (0, 10):
            if s - k >= 1:
                want.add(s - k)
    for s in sorted(cmp_):
        if s in want or s <= 2:
            c = cmp_[s]
            print(f"| {s} | {c['n']} | {fmt(c['L1'])} | {fmt(c['Linf'])} | {c['n01']} | {c['n1']} | {c['plane']} / {c['cut']} / {c['uncut']} |")
    if not cmp_:
        print("| (no comparison run) | | | | | | |")
    print()


if __name__ == "__main__":
    main()
