#!/usr/bin/env python3
"""Plane gate of a switch variant (2026-11-16): the Stage B ccsplane gate (kiss_an.mode_plane, ccsGate thresholds)
over studies swg_<variant>_plane_*_np[18], pass counts, the failing cases with the switch counters.
usage: gate_sw.py <variant e.g. S1S2low>"""
import argparse, csv, glob, io, os, re, sys, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kiss_an as ka

v = sys.argv[1]
_collect = ka.collect
def _wrap(*a):
    recs = _collect(*a)
    for r in recs:
        if str(r.get("setting", "")).startswith("onePlaneS"):
            r["setting"] = "onePlane"   # the gate applies its onePlane thresholds to the switch variant
    return recs
ka.collect = _wrap
_read = ka.read_case
def _safe(c):
    try:
        return _read(c)
    except ValueError as e:
        print("UNREADABLE LOG (garbled by a killed or duplicate process):", c, file=sys.stderr)
        return None
ka.read_case = _safe
out = os.path.join(ka.REPO, "studies", "davof", "sw_tables")
os.makedirs(out, exist_ok=True)
args = argparse.Namespace(studies_dir=os.path.join(ka.REPO, "studies", "davof"), out=out)
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ka.mode_plane(args, f"swg_{v}_plane_*_np[18]", f"swg_{v}", ka.IDENT_TOL_3B, True)
txt = buf.getvalue()
open(os.path.join(out, f"swg_{v}.txt"), "w").write(txt)
for l in txt.splitlines():
    if l.startswith("P1 (round 3)") or l.startswith("identity ("):
        print(l)
recs = ka.collect(args.studies_dir, f"swg_{v}_plane_*_np[18]")
rows = list(csv.DictReader(open(os.path.join(out, f"kiss_swg_{v}_table.csv"))))
print(f"cases {len(rows)}; statuses: " + ", ".join(f"{s} {sum(1 for r in rows if r['status'] == s)}" for s in sorted(set(r['status'] for r in rows))))
ok = [r for r in rows if r["verdict"].startswith("P1 ok")]
op = [r for r in rows if r["verdict"].startswith("OPEN")]
bad = [r for r in rows if r["verdict"].startswith("P1 FAIL") or (r["verdict"].startswith("OPEN") and "FAIL" in r["verdict"])]
print(f"P1 ok {len(ok)}, P1 FAIL {sum(1 for r in rows if r['verdict'].startswith('P1 FAIL'))}, OPEN cells (Co 1 normal arm) {len(op)} of which failing {sum(1 for r in op if 'FAIL' in r['verdict'])}")
cdir = {(r["study"], r["N"], round(r["co"], 5)): os.path.join(args.studies_dir, r["study"], r["case"]) for r in recs}
RE_SW = re.compile(r"switches acted on faces: S1 (\d+), S2 (\d+) \(upwind alpha < 0.5: (\d+), >= 0.5: (\d+)\)")
print("\n| study | N | Co | steps | status | L1(alpha_f) | Linf(alpha_f) | E_shape | min alpha | limiter faces | verdict | first step a switch acted | S1 faces total (max/step) | S2 faces total (max/step) | S2 alpha_up<0.5 total | S2 alpha_up>=0.5 total |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (r["study"], int(r["N"]), float(r["Co"]))):
    if r["verdict"].startswith("P1 ok"):
        continue
    d = cdir.get((r["study"], int(r["N"]), round(float(r["Co"]), 5)))
    sw = []
    if d and os.path.isfile(d + "/log.davofAdvectorFoam"):
        sw = [tuple(int(x) for x in m.groups()) for m in RE_SW.finditer(open(d + "/log.davofAdvectorFoam", errors="replace").read())]
    first = next((i + 1 for i, t in enumerate(sw) if t[0] + t[1] > 0), "none")
    print(f"| {r['study']} | {r['N']} | {r['Co']} | {r['steps']} | {r['status']} | {r['L1af']} | {r['Linfaf']} | {r['Eshape']} | {r['minAlpha']} | {r['limFaces']} | {r['verdict']} | {first} | "
          f"{sum(t[0] for t in sw)} ({max([t[0] for t in sw] or [0])}) | {sum(t[1] for t in sw)} ({max([t[1] for t in sw] or [0])}) | {sum(t[2] for t in sw)} | {sum(t[3] for t in sw)} |")
