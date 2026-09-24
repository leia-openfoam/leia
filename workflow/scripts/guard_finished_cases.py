#!/usr/bin/env python3
"""Refuse a study launch that would re-run a FINISHED case.

    guard_finished_cases.py config/<study>.yaml --profile profiles/<name> [--studies-dir DIR]
                            [-- <extra snakemake arguments>]

The Makefile calls this before every `snakemake` it starts. It runs the same dry run the
launch would run (`--workflow-profile`, `--nolock`, `--rerun-triggers mtime`, `-n`), reads
the `solve` jobs it lists, and looks at each target case's solver log with
`workflow/scripts/foam_log_state.sh`. A case whose log is COMPLETED or DIVERGED holds a
result; the `solve` rule deletes the case CSV before it runs, so re-running it destroys
that result. If any listed `solve` job targets such a case, the guard prints them and
exits 2, and the launch does not happen.

Why. MEASURED 2026-09-22: a touch sweep over /work/scratch gave every file a new mtime,
after which the Makefile's `--rerun-triggers mtime` would have re-run 6 of the 42 finished
cases of sdplsConv2Dvortex and 2 of 14 of sdplsExpSource2Dvortex (STATUS.md section 9,
docs/plan-library-split-and-build-policy.md WP4). Snakemake's own view of "up to date"
is file times; the guard adds the one thing that matters: is there a finished result here.

Override: LEIA_ALLOW_RERUN=1 turns the refusal into a warning and exits 0. Preserve the
study directory FIRST (rename it with a dated suffix; CLAUDE.md, "Provenance and
preservation"), then set the variable for that one launch.

Exit codes: 0 no finished case would be re-run (or override set); 2 refusal;
the dry run's own code when the dry run fails (its last lines are printed)."""
import argparse, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FINISHED = {0: "COMPLETED", 3: "DIVERGED"}   # foam_log_state.sh exit codes


def solver_csv(solver):
    """The file the `solve` rule waits for -- the same rule as workflow/Snakefile."""
    return "leiaLevelSetFoam.csv" if solver == "leiaLevelSetTwoPhaseFoam" else f"{solver}.csv"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("config")
    ap.add_argument("--profile", required=True, help="workflow profile directory, as the Makefile passes it")
    ap.add_argument("--studies-dir", default=None, help="override studies_dir (default: the config's, else <repo>/studies)")
    ap.add_argument("extra", nargs="*", help="extra snakemake arguments after --")
    a = ap.parse_args()

    # the four top-level scalars this guard needs, read without PyYAML (stdlib only, so the
    # Makefile's plain `python3` works on every machine)
    cfg = {}
    for line in open(a.config):
        m = re.match(r"^(study_name|case|solver|studies_dir):\s*([^#\s]+)", line)
        if m:
            cfg[m.group(1)] = m.group(2).strip("\"'")
    study, case = cfg["study_name"], cfg["case"]
    solver = cfg.get("solver", "leiaLevelSetFoam")
    studies_dir = a.studies_dir or cfg.get("studies_dir") or os.path.join(REPO, "studies")

    cmd = ["snakemake", "--workflow-profile", a.profile, "--configfile", a.config,
           "--nolock", "--rerun-triggers", "mtime", "-n"] + a.extra
    env = dict(os.environ, PATH=os.path.expanduser("~/.local/bin") + ":" + os.environ.get("PATH", ""))
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, env=env)
    if r.returncode != 0:
        sys.stderr.write("guard: the dry run failed; its last lines:\n" +
                         "\n".join((r.stdout + r.stderr).strip().splitlines()[-15:]) + "\n")
        sys.exit(r.returncode)

    # job blocks of the dry run: "rule solve:" ... "wildcards: idx=00003"
    solve_idx = []
    rule = None
    for line in r.stdout.splitlines():
        m = re.match(r"^(?:local)?rule (\w+):", line)
        if m:
            rule = m.group(1); continue
        m = re.match(r"^\s+wildcards:.*\bidx=(\d+)", line)
        if m and rule == "solve":
            solve_idx.append(m.group(1))

    finished = []
    for idx in sorted(set(solve_idx)):
        d = os.path.join(studies_dir, study, f"{case}_{idx}")
        log = os.path.join(d, f"log.{solver}")
        csv = os.path.join(d, solver_csv(solver))
        if not os.path.isfile(log):
            continue
        s = subprocess.run([os.path.join(HERE, "foam_log_state.sh"), log], capture_output=True, text=True)
        if s.returncode in FINISHED:
            finished.append((idx, FINISHED[s.returncode], s.stdout.strip(), os.path.isfile(csv)))

    n = len(set(solve_idx))
    if not finished:
        print(f"guard: {study}: {n} solve job(s) listed, none over a finished case")
        return 0
    print(f"guard: {study}: {len(finished)} of {n} listed solve job(s) would RE-RUN a finished case:")
    for idx, state, line, has_csv in finished:
        print(f"    {case}_{idx}  {line}  csv={'present' if has_csv else 'absent'}")
    if os.environ.get("LEIA_ALLOW_RERUN") == "1":
        print("guard: LEIA_ALLOW_RERUN=1 -- launching anyway (the study directory must have been preserved first)")
        return 0
    print("guard: REFUSED. Preserve the study directory (rename it with a dated suffix), then\n"
          "       either mark it up to date (snakemake ... --touch) or set LEIA_ALLOW_RERUN=1 for this launch.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
