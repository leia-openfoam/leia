#!/usr/bin/env python3
"""Curate a mass-momentum consistency table from a translating-droplet study.

Reads each arm's RESOLVED dictionaries (never the study config: a config states
intent, a rendered dictionary states what ran) plus its metrics CSV, and reports the
vector the consistency argument turns on:

  free-stream error   maxMagUPrime = max|U - U_ref|.  THE score for a translating
                      case. max|U| is dominated by U0 and cannot see a parasitic
                      current below it -- a 3D arm was once read as healthy while
                      dying because of exactly that.
  mass residual       ddt(rho) + div(rhoPhi), absolute and relative. For a droplet
                      translating at uniform U0 the momentum equation inherits U0
                      times this; it is not a gradient so the pressure cannot absorb
                      it, and it vanishes identically at U0 = 0.
  clip                what massFlux.boundRho had to remove to keep rho physical.
  volume / shape      the standing rule: never a headline metric alone.

Usage:  make_consistency_table.py <studyDir> [--csv out.csv]
"""
from __future__ import annotations
import argparse, csv, glob, os, re, sys

DT_TO_N = {"5.99991e-05": 32, "2.12129e-05": 64, "7.49989e-06": 128,
           "2.65161e-06": 256}


def strip_comments(text: str) -> str:
    return re.sub(r"//[^\n]*", "", re.sub(r"/\*.*?\*/", "", text, flags=re.S))


def dict_entry(body: str, key: str, default: str = "?") -> str:
    m = re.search(rf"\b{re.escape(key)}\s+([^;]+);", body)
    return m.group(1).strip() if m else default


def sub_dict(text: str, name: str) -> str:
    """Body of a named sub-dictionary.

    The closing brace sits at whatever indent the block was opened at -- top-level
    blocks like ddtSchemes close at column 0, nested ones like massFlux at four
    spaces -- so match the first line that is nothing but a closing brace rather
    than assuming a depth.
    """
    m = re.search(rf"(?m)^(\s*){re.escape(name)}\s*$\s*^\s*\{{", text)
    if m is None:
        m = re.search(rf"{re.escape(name)}\s*\{{", text)
        if m is None:
            return ""
        start = m.end()
    else:
        start = m.end()
    close = re.search(r"(?m)^\s*\}", text[start:])
    return text[start:start + close.start()] if close else text[start:]


def arm_settings(arm: str) -> dict:
    fv = strip_comments(open(os.path.join(arm, "system/fvSolution")).read())
    cd = open(os.path.join(arm, "system/controlDict")).read()
    sch = strip_comments(open(os.path.join(arm, "system/fvSchemes")).read())
    dt = dict_entry(cd, "deltaT")
    mf = sub_dict(fv, "massFlux")
    return {
        "N": DT_TO_N.get(dt, "?"),
        "deltaT": dt,
        "massFlux": dict_entry(mf, "type"),
        "boundRho": dict_entry(mf, "boundRho", "-"),
        "momentumDdt": dict_entry(sub_dict(sch, "ddtSchemes"), "default"),
        "rhoDdt": dict_entry(sub_dict(sch, "ddtSchemes"), "ddt(rho)", "(default)"),
        "curvatureExtension": dict_entry(sub_dict(fv, "curvatureExtension"), "type"),
        "trace": dict_entry(sub_dict(fv, "semiLagrangian"), "traceVelocity"),
    }


def arm_metrics(arm: str, solver: str) -> dict:
    path = os.path.join(arm, f"{solver}.csv")
    if not os.path.exists(path):
        return {}
    rows = [r for r in csv.DictReader(open(path)) if r.get("TIME")]
    if not rows:
        return {}

    def col(key):
        out = []
        for r in rows:
            v = r.get(key)
            if v in (None, ""):
                continue
            try:
                out.append(float(v))
            except ValueError:
                pass
        return out

    def last(key, default=float("nan")):
        v = col(key)
        return v[-1] if v else default

    def mx(key, default=float("nan")):
        v = col(key)
        return max(v) if v else default

    return {
        "rows": len(rows),
        "t_last": last("TIME"),
        "freeStreamLast": last("maxMagUPrime"),
        "freeStreamMax": mx("maxMagUPrime"),
        "massResL1": last("rhoMassResidualL1"),
        "massResRelL1": last("rhoMassResidualRelL1"),
        "clipL1": mx("rhoClipL1", 0.0),
        "clipFracMax": mx("rhoClipFraction", 0.0),
        "volErr": last("phaseVolumeRelError"),
        "shapeL2": last("zeroSetRadialL2"),
    }


def log_state(arm: str, solver: str) -> tuple[int, str]:
    log = os.path.join(arm, f"log.{solver}")
    if not os.path.exists(log):
        return 0, "NO LOG"
    text = open(log, errors="replace").read()
    steps = text.count("\nTime = ")
    end = re.search(r"^End$", text, re.M) is not None
    return steps, ("COMPLETED" if end else "did NOT finish")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--solver", default="leiaSemiLagrangianLevelSetTwoPhaseFoam")
    ap.add_argument("--csv")
    args = ap.parse_args()

    arms = sorted(d for d in glob.glob(os.path.join(args.study, "*/"))
                  if os.path.isfile(os.path.join(d, "system/fvSolution")))
    if not arms:
        print(f"no arms with a rendered fvSolution under {args.study}", file=sys.stderr)
        return 1

    recs = []
    for arm in arms:
        rec = {"arm": os.path.basename(arm.rstrip("/"))}
        rec.update(arm_settings(arm))
        steps, state = log_state(arm, args.solver)
        rec["steps"], rec["state"] = steps, state
        rec.update(arm_metrics(arm, args.solver))
        recs.append(rec)
    recs.sort(key=lambda r: (r["N"] if isinstance(r["N"], int) else 0,
                             r["massFlux"], r["momentumDdt"]))

    hdr = (f"{'N':<5}{'massFlux':<22}{'momDdt':<10}{'steps':<8}{'state':<16}"
           f"{'freeStream':<12}{'massResRelL1':<14}{'clipL1':<11}{'volErr':<11}{'shapeL2':<11}")
    print(hdr)
    print("-" * len(hdr))
    for r in recs:
        def f(key, w=12, p=4):
            v = r.get(key)
            return f"{v:<{w}.{p}g}" if isinstance(v, float) else f"{'-':<{w}}"
        print(f"{r['N']:<5}{r['massFlux']:<22}{r['momentumDdt']:<10}"
              f"{r['steps']:<8}{r['state']:<16}"
              f"{f('freeStreamLast')}{f('massResRelL1',14)}{f('clipL1',11)}"
              f"{f('volErr',11)}{f('shapeL2',11)}")

    if args.csv:
        fields = sorted({k for r in recs for k in r})
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(recs)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
