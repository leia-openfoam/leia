#!/usr/bin/env python3
"""The rendering identity of the two routes of the DAVOF advection studies
(docs/davof/PORT-DESIGN.md phase 3, criterion 2, 2026-11-06): every variation of a
study config rendered with the fork's solver (advectorVoF) and with leia's
(davofAdvectorFoam) must differ in exactly three lines of the case, the
`application` line of system/controlDict, the `advectionScheme` line of
system/fvSolution and the function objects' `libs` entries of system/controlDict,
and in nothing else (the setfields app is the Snakefile's, not the case's).

Usage:
    port_render_diff.py config/davof/<study>.yaml [more configs ...] [--keep DIR]

Mirrors the Snakefile's case construction (foam_param.build_token_grid, expand,
materialize.materialize) with DAVOF_SOLVER set to each route in turn. Exit 0 when
every variation of every config differs in the expected lines only.
"""
import argparse, os, re, sys, shutil, tempfile, difflib
import yaml
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import foam_param, materialize   # noqa: E402
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
EXPECTED = (re.compile(r"^\s*application\s"), re.compile(r"^\s*advectionScheme\s"), re.compile(r"^\s*libs\s"))

def render(cfg, solver, out_root):
    case = cfg["case"]; cases_dir = os.path.join(REPO, "cases")
    case_name = os.path.basename(case)
    mesh = cfg.get("mesh", "hex"); mode = cfg.get("mode", "serial"); np_ = int(cfg.get("np", 1))
    dims = 1 if case_name.startswith("1D") else (2 if case_name.startswith("2D") else 3)
    param = os.path.join(cases_dir, f"{case}_poly.parameter" if mesh in ("poly", "polyRefined") else f"{case}.parameter")
    base = os.path.join(cases_dir, case)
    axes, consts, _ = foam_param.build_token_grid(param, os.path.join(cases_dir, "default.parameter"), base,
                                                  axes_override=cfg.get("axes_override"),
                                                  collapse_other_axes=cfg.get("collapse_other_axes", False))
    vars_ = foam_param.expand(axes, consts)
    dirs = []
    for i, v in enumerate(vars_):
        d = os.path.join(out_root, f"{case_name}_{i:05d}")
        materialize.materialize(base, {**v, "DAVOF_SOLVER": solver}, d, np_, mesh, mode, dims, case_name, f"{i:05d}")
        dirs.append(d)
    return dirs

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("configs", nargs="+"); ap.add_argument("--keep", default=None)
    a = ap.parse_args()
    root = a.keep or tempfile.mkdtemp(prefix="port_render_")
    allok = True
    for c in a.configs:
        cfg = yaml.safe_load(open(c))
        fork = render(cfg, "advectorVoF", os.path.join(root, os.path.basename(c)[:-5], "fork"))
        leia = render(cfg, "davofAdvectorFoam", os.path.join(root, os.path.basename(c)[:-5], "leia"))
        for df, dl in zip(fork, leia):
            unexpected, expected_hits = [], set()
            for dp, dn, fn in os.walk(df):
                for f in fn:
                    pf = os.path.join(dp, f); rel = os.path.relpath(pf, df); pl = os.path.join(dl, rel)
                    if f == "case_params.json": continue
                    if not os.path.isfile(pl): unexpected.append(f"missing on leia: {rel}"); continue
                    A = open(pf, errors="replace").read().split("\n"); B = open(pl, errors="replace").read().split("\n")
                    if A == B: continue
                    for line in difflib.unified_diff(A, B, lineterm="", n=0):
                        if line.startswith(("---", "+++", "@@")): continue
                        body = line[1:]
                        if any(p.match(body) for p in EXPECTED): expected_hits.add(rel + ":" + body.strip().split()[0]); continue
                        unexpected.append(f"{rel}: {line[:120]}")
            ok = not unexpected and len(expected_hits) >= 3
            allok &= ok
            print(("OK  " if ok else "BAD ") + os.path.relpath(df, root) + "  expected lines: " + ", ".join(sorted(expected_hits)))
            for u in unexpected[:10]: print("     " + u)
    if not a.keep: shutil.rmtree(root, ignore_errors=True)
    print("RENDERING IDENTICAL (three lines apart)" if allok else "RENDERING DIFFERS BEYOND THE THREE LINES")
    sys.exit(0 if allok else 1)

if __name__ == "__main__":
    main()
