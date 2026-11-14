#!/usr/bin/env python3
"""diag2 config generator (2026-11-14): writes the Part 1 dump studies, the check study and the
Part 2 profile and scaling studies from the Stage B configs (config/davof/ccsB_*)."""
import re, sys, os

CFG = "config/davof"
HDR = "# diag2 (2026-11-14; docs/davof/REPORT-swept.md 'Pieces path: aborts and cost'): %s\n"

def base(fam, N, row):
    suffix = {"davof": "davof", "plic": "plic", "iso": "iso"}[row]
    return open(f"{CFG}/ccsB_{fam}_N{N}_{suffix}.yaml").read()

def write(name, text):
    open(f"{CFG}/{name}.yaml", "w").write(text)
    print("wrote", name)

def make(name, fam, N, row, targets, np=None, settings=None, end_time=None, profiling=True,
         old=False, runtime=None, descr=""):
    s = base(fam, N, row)
    # header comments: replace the leading comment block by one diag2 line and the sbatch line
    lines = [l for l in s.splitlines() if not l.startswith("#")]
    s = "\n".join(lines) + "\n"
    s = HDR % descr + f"# sbatch -J {name} --export=ALL,STUDY=davof/{name} run-studies.sbatch\n" + s
    s = re.sub(r"^study_name: .*$", f"study_name: {name}", s, flags=re.M)
    if np is not None:
        s = re.sub(r"^np: .*$", f"np: {np}", s, flags=re.M)
        s = re.sub(r"^mode: .*$", "mode: serial" if np == 1 else "mode: parallel", s, flags=re.M)
    s = re.sub(r"TARGET_ADVECTIVE_CO: \[.*\]", "TARGET_ADVECTIVE_CO: [" + ", ".join(str(t) for t in targets) + "]", s)
    extra = ""
    if end_time is not None:
        extra += f"  END_TIME: [{end_time}]\n"
    if profiling:
        extra += '  TPF_PROFILING: ["true"]\n'
    s = s.replace("  TPF_REPETITION: [1]\n", "  TPF_REPETITION: [1]\n" + extra)
    if row == "davof":
        s = re.sub(r"TPF_ENRIGHT_SETTING: \[.*\]", "TPF_ENRIGHT_SETTING: [" + ", ".join(settings) + "]", s)
        if old:
            s = s.replace("TPF_GRAZING_CUT: [openfoam]", "TPF_GRAZING_CUT: [fork]")
            s = s.replace("TPF_CCS_END_INTERVAL: [cubic]", "TPF_CCS_END_INTERVAL: [none]")
            s = s.replace("TPF_CCS_TOLERANCE: [relative]", "TPF_CCS_TOLERANCE: [absolute]")
            s = s.replace('TPF_POSITIONING_TOL: ["1e-12"]', 'TPF_POSITIONING_TOL: ["1e-14"]')
    if runtime is not None:
        s = re.sub(r"^solve_runtime: .*$", f"solve_runtime: {runtime}", s, flags=re.M)
    write(name, s)

T95 = {("deformation3D", 64): 0.95955, ("deformation3D", 128): 0.95834}

which = sys.argv[1]
if which == "part1":
    # the 1.84e+02 case with both key sets and the print-only dump; shear2D N 256, 8 ranks, Co 0.2
    for tag, old in (("new", False), ("old", True)):
        make(f"diag2_dump_shear2D_N256_{tag}", "shear2D", 256, "davof", [0.2], settings=["recommendedDump"],
             old=old, profiling=False, runtime=600,
             descr=f"shear2D N 256 (8 ranks), pieces + alphaDump, {'OLD' if old else 'NEW'} keys")
    # the old-key rerun of the third plain case is the dump study (print-only key): see the pre-registration
elif which == "check":
    make("diag2_check_shear2D_N64", "shear2D", 64, "davof", [0.5],
         settings=["recommended", "recommendedProfile", "recommendedDump"], profiling=False, runtime=120,
         descr="shear2D N 64, key off / phaseProfile / alphaDump: identity of the dictionaries and of the errors")
elif which == "part2":
    for fam, N, np, targets, et in (
        ("deformation3D", 64, 8, [0.2, T95[("deformation3D", 64)]], None),
        ("deformation3D", 128, 64, [0.2, T95[("deformation3D", 128)]], 0.6),
        ("shear3D", 64, 8, [0.2], None),
    ):
        rt = {64: 600, 128: 1400}[N]
        make(f"diag2_prof_{fam}_N{N}_davof", fam, N, "davof", targets, np=np, settings=["recommendedProfile"],
             end_time=et, runtime=rt, descr=f"{fam} N {N} ({np} ranks), pieces + phaseProfile, NEW keys")
        make(f"diag2_prof_{fam}_N{N}_plic", fam, N, "plic", targets, np=np, end_time=et, runtime=rt,
             descr=f"{fam} N {N} ({np} ranks), plicRDF + OpenFOAM profiling")
        make(f"diag2_prof_{fam}_N{N}_iso", fam, N, "iso", targets, np=np, end_time=et, runtime=rt,
             descr=f"{fam} N {N} ({np} ranks), isoAlpha + OpenFOAM profiling")
elif which == "scal":
    for np in (1, 8, 64):
        make(f"diag2_scal_deformation3D_N64_np{np}", "deformation3D", 64, "davof", [0.2], np=np,
             settings=["recommendedProfile"], end_time=0.25, runtime=600 if np == 1 else 120,
             descr=f"deformation3D N 64 strong scaling, pieces + phaseProfile, {np} ranks, first 200 steps (END_TIME 0.25)")
