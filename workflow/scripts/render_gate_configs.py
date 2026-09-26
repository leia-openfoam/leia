#!/usr/bin/env python3
"""Render a method gate (config/gates/<gate>.yaml) and its candidates into one complete
study config per arm, for workflow/Snakefile (CLAUDE.md "Method gates").

A method enters a gate only as a CANDIDATE: a set of case tokens, from
config/candidates/<name>.yaml or from SET="line=...,TOKEN=value,..." on the command line.
This script merges gate + candidate into complete study configs with EXPLICIT collision
errors, so the shallow-merge trap of a base config (CLAUDE.md, "Never add a
config/best.yaml") cannot occur. It refuses:

  * a candidate token that is not a method token of the gate;
  * a token of the other solver line (SL_* in an Eulerian candidate, and the reverse);
  * a candidate token that collides with an arm token or a line token;
  * a candidate token that an arm's case does not render (no @!TOKEN!@ in its templates).
    Without this check foam_param drops the token silently and the arm runs the baseline
    under the candidate's name.

The METHOD settings are not in the gate file: the baseline is the .parameter layering,
the executable form of METHOD.md.

Outputs, under <studies>/<gate>[Smoke]_summary/<candidate>/:
  configs/<arm>.yaml    one study config per arm (and per seam decomposition)
  candidate.json        the resolved candidate (name, line, tokens, preRegistered)
  manifest.json         the arms: study name, config path, np, first, kind, case, ladder

Usage (the Makefile target `gate` and workflow/Snakefile.gate call it):
  render_gate_configs.py --gate config/gates/methodGate2D.yaml --candidates baseline+HL1z
      [--set "line=semiLagrangian,VELOCITY_EXTENSION=haloLimited"] [--smoke]
      [--studies-dir studies] [--preserve] [--list]
"""
import argparse
import copy
import datetime
import hashlib
import json
import os
import re
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import foam_param  # noqa: E402


class GateError(Exception):
    """A refused gate/candidate combination; main() prints it and exits with code 2."""


def load_yaml(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def rendered_tokens(case):
    """Every token the case's templates render, with the inputs of derived tokens
    (foam_param._DERIVED_INPUTS: DROPLET_OFFSET_X enters through DROPLET_CENTRE_X)."""
    ref = set(foam_param.referenced_tokens(os.path.join(REPO, "cases", case)))
    for tok in list(ref):
        ref.update(foam_param._DERIVED_INPUTS.get(tok, ()))
    return ref


def parse_set(set_string):
    """SET="line=semiLagrangian,TOKEN=value,GC_M_MU=1" -> (line, tokens)."""
    line, tokens = None, {}
    for item in [s for s in set_string.split(",") if s.strip()]:
        if "=" not in item:
            raise GateError(f"SET item {item!r} is not TOKEN=value")
        k, v = (x.strip() for x in item.split("=", 1))
        if k == "line":
            line = v
        else:
            tokens[k] = v
    return line or "semiLagrangian", tokens


def load_candidate(name, set_string=""):
    """A candidate from config/candidates/<name>.yaml, or the ad-hoc candidate of SET."""
    if name.startswith("adhoc"):
        line, tokens = parse_set(set_string)
        digest = hashlib.sha1(json.dumps([line, sorted(tokens.items())]).encode()).hexdigest()[:8]
        return {"candidate": f"adhoc-{digest}", "line": line, "tokens": tokens,
                "description": f"ad-hoc candidate from SET={set_string!r}",
                "preRegistered": False, "source": "SET"}
    path = os.path.join(REPO, "config", "candidates", f"{name}.yaml")
    if not os.path.isfile(path):
        raise GateError(f"no candidate file {os.path.relpath(path, REPO)}")
    cand = load_yaml(path) or {}
    if cand.get("candidate") != name:
        raise GateError(f"{os.path.relpath(path, REPO)} declares candidate "
                        f"{cand.get('candidate')!r}, not {name!r}")
    if not re.fullmatch(r"[A-Za-z0-9-]+", name):
        raise GateError(f"candidate name {name!r}: use letters, digits and '-' only "
                        f"(study names are <gate>_<candidate>_<arm>)")
    cand.setdefault("tokens", {})
    cand.setdefault("rates", {})
    cand["preRegistered"] = True
    cand["source"] = os.path.relpath(path, REPO)
    return cand


def split_rates(gate, cand):
    """Separate dimensionless rates (e.g. GC_M_MU) from plain tokens."""
    rates = dict(cand.get("rates") or {})
    tokens = {}
    for k, v in (cand.get("tokens") or {}).items():
        if k in gate.get("rates", {}):
            rates[k] = v
        else:
            tokens[k] = v
    return tokens, rates


def check_candidate(gate, cand, tokens):
    line = cand.get("line")
    if line not in gate["solvers"]:
        raise GateError(f"candidate {cand['candidate']}: line {line!r} is not one of "
                        f"{sorted(gate['solvers'])}")
    method = set(gate["methodTokens"])
    other = {t for ln, toks in gate.get("lineOf", {}).items() if ln != line for t in toks}
    for tok in tokens:
        if tok not in method:
            raise GateError(f"candidate {cand['candidate']}: {tok} is not a method token of "
                            f"{gate['gate']} (methodTokens); arm setup and METHOD.md settings "
                            f"are not a candidate's to change")
        if tok in other:
            raise GateError(f"candidate {cand['candidate']} (line {line}): {tok} belongs to "
                            f"the other solver line")


def render_arm(gate, arm_name, arm, cand, tokens, rates, smoke, np_override=None):
    line = cand["line"]
    case = arm["case"]
    if not os.path.isdir(os.path.join(REPO, "cases", case)):
        raise GateError(f"arm {arm_name}: case {case!r} does not exist under cases/")
    ref = rendered_tokens(case)
    solver = gate["solvers"][line][arm["kind"]]
    np_ = int(np_override if np_override is not None else arm.get("np", gate["np"]))
    if smoke and np_override is None and (arm.get("smoke") or {}).get("np"):
        np_ = int(arm["smoke"]["np"])       # a laptop smoke cannot start the cluster's ranks

    ladder = list((arm.get("smoke") or {}).get("N_CELLS", arm["N_CELLS"])) if smoke \
        else list(arm["N_CELLS"])
    arm_tokens = dict(arm.get("tokens") or {})
    if smoke:
        arm_tokens["END_TIME"] = arm["smoke"]["END_TIME"]
    line_tokens = dict((gate.get("lineTokens") or {}).get(line) or {})

    for tok in tokens:
        if tok in arm_tokens or tok in line_tokens:
            raise GateError(f"candidate {cand['candidate']}: {tok} collides with an arm or "
                            f"line token of arm {arm_name}")
        if tok not in ref:
            raise GateError(f"candidate {cand['candidate']} sets {tok}, which case {case} "
                            f"(arm {arm_name}) does not render. Add @!{tok}!@ to its "
                            f"templates with the inert default first.")
    for tok in arm_tokens:
        if tok not in ref:
            raise GateError(f"arm {arm_name}: token {tok} is not rendered by case {case}")
    axes = {"N_CELLS": ladder}
    for k, v in arm_tokens.items():
        axes[k] = [v]
    for k, v in line_tokens.items():
        if k in ref:            # e.g. ADVECTION exists only in the kinematic cases
            axes[k] = [v]
    for k, v in tokens.items():
        axes[k] = [v]
    for rate, target in (gate.get("rates") or {}).items():
        if rate in rates:
            if target not in ref:
                raise GateError(f"candidate {cand['candidate']} sets {rate}; case {case} "
                                f"does not render {target}")
            axes[target] = [f"{float(rates[rate])/float(arm['T_REF']):.10g}"]
    axes = {k: [str(x) for x in v] for k, v in axes.items()}

    cfg = {
        "study_name": None,           # set by the caller
        "case": case,
        "mesh": "hex",
        "mode": "parallel" if np_ > 1 else "serial",
        "np": np_,
        "solver": solver,
        "theme": gate["theme"],
        "solve_runtime": int(arm["solve_runtime"]),
        "axes_override": axes,
        "collapse_other_axes": True,
        "export_slides": False,
    }
    if arm.get("setfields_args"):
        cfg["setfields_args"] = arm["setfields_args"]
    return cfg


def write_if_changed(path, text):
    """Write only a changed file: an unchanged mtime keeps snakemake from re-running."""
    if os.path.isfile(path) and open(path).read() == text:
        return False
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return True


def render(gate_path, candidates, set_string="", smoke=False, studies_dir=None):
    gate = load_yaml(gate_path)
    gate_name = gate["gate"] + ("Smoke" if smoke else "")
    studies_dir = os.path.abspath(studies_dir or os.path.join(REPO, "studies"))
    summary_root = os.path.join(studies_dir, f"{gate_name}_summary")
    manifests = []
    for name in candidates:
        cand = load_candidate(name, set_string)
        tokens, rates = split_rates(gate, cand)
        check_candidate(gate, cand, tokens)
        cname = cand["candidate"]
        cdir = os.path.join(summary_root, cname)
        arms = []
        for arm_name, arm in gate["arms"].items():
            cfg = render_arm(gate, arm_name, arm, cand, tokens, rates, smoke)
            study = f"{gate_name}_{cname}_{arm_name}"
            cfg["study_name"] = study
            if smoke:
                cfg["studies_dir"] = studies_dir
            path = os.path.join(cdir, "configs", f"{arm_name}.yaml")
            header = (f"# RENDERED by workflow/scripts/render_gate_configs.py from "
                      f"{os.path.relpath(gate_path, REPO)} and candidate {cname} "
                      f"({cand['source']}). Do not edit; edit the gate or the candidate.\n")
            write_if_changed(path, header + yaml.safe_dump(cfg, sort_keys=False))
            arms.append({"arm": arm_name, "study": study, "config": path,
                         "np": cfg["np"], "first": bool(arm.get("first")),
                         "kind": arm["kind"], "case": arm["case"],
                         "ladder": cfg["axes_override"]["N_CELLS"]})
            seam = arm.get("seam") or {}
            for np_seam in ((seam.get("smokeNp") if smoke and seam.get("smokeNp") else seam.get("np")) or []):
                scfg = render_arm(gate, arm_name, arm, cand, tokens, rates, smoke,
                                  np_override=np_seam)
                scfg["axes_override"]["N_CELLS"] = [cfg["axes_override"]["N_CELLS"][0]]
                sstudy = f"{gate_name}_{cname}_seamNp{np_seam}"
                scfg["study_name"] = sstudy
                if smoke:
                    scfg["studies_dir"] = studies_dir
                spath = os.path.join(cdir, "configs", f"seamNp{np_seam}.yaml")
                write_if_changed(spath, header + yaml.safe_dump(scfg, sort_keys=False))
                arms.append({"arm": f"seamNp{np_seam}", "seamOf": arm_name, "study": sstudy,
                             "config": spath, "np": scfg["np"], "first": False,
                             "kind": arm["kind"], "case": arm["case"],
                             "ladder": scfg["axes_override"]["N_CELLS"]})
        resolved = {"candidate": cname, "line": cand["line"], "tokens": tokens,
                    "rates": rates, "preRegistered": cand["preRegistered"],
                    "source": cand["source"], "description": cand.get("description", ""),
                    "target": cand.get("target"), "gate": gate_name,
                    "gateFile": os.path.relpath(gate_path, REPO), "smoke": smoke}
        write_if_changed(os.path.join(cdir, "candidate.json"),
                         json.dumps(resolved, indent=2, sort_keys=True) + "\n")
        manifest = {"gate": gate_name, "candidate": cname, "line": cand["line"],
                    "summaryDir": cdir, "arms": arms}
        write_if_changed(os.path.join(cdir, "manifest.json"),
                         json.dumps(manifest, indent=2) + "\n")
        manifests.append(manifest)
    return gate, manifests


def preserve(manifests, studies_dir):
    """Rename every existing study of the rendered candidates with a dated suffix
    (CLAUDE.md "Provenance and preservation"). A rename, never a delete."""
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    for m in manifests:
        for a in m["arms"]:
            d = os.path.join(studies_dir, a["study"])
            if os.path.isdir(d):
                os.rename(d, f"{d}_pre-{stamp}")
                print(f"[gate] preserved {a['study']} -> {os.path.basename(d)}_pre-{stamp}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gate", required=True)
    ap.add_argument("--candidates", required=True,
                    help="'+'-separated candidate names; 'adhoc' uses --set")
    ap.add_argument("--set", default="", help="ad-hoc candidate tokens")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--studies-dir", default=None)
    ap.add_argument("--preserve", action="store_true")
    ap.add_argument("--list", action="store_true", help="print the study names")
    a = ap.parse_args(argv)
    names = [c for c in a.candidates.replace(",", "+").split("+") if c]
    if a.set and "adhoc" not in names:
        names.append("adhoc")
    try:
        _gate, manifests = render(a.gate, names, a.set, a.smoke, a.studies_dir)
    except GateError as e:
        print(f"[gate] REFUSED: {e}", file=sys.stderr)
        return 2
    studies_dir = os.path.abspath(a.studies_dir or os.path.join(REPO, "studies"))
    if a.preserve:
        preserve(manifests, studies_dir)
    for m in manifests:
        print(f"[gate] {m['gate']} candidate {m['candidate']} ({m['line']}): "
              f"{len(m['arms'])} studies -> {os.path.relpath(m['summaryDir'], REPO)}")
        if a.list:
            for arm in m["arms"]:
                print(f"    {arm['study']:48s} np {arm['np']:3d}  N {arm['ladder']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
