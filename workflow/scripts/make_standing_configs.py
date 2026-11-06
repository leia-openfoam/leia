#!/usr/bin/env python3
"""The standing study matrix of the DAVOF advection in leia (Tomislav, 2026-11-07;
docs/davof/PORT-DESIGN.md section 5b): writes config/davof/standing/<family>_N<N>.yaml,
one study per case family and resolution, from the table below, so that the nine configs
agree in everything but the family and the rung. Run it after a change of the matrix
(the rows, the Courant numbers, the ranks, the limits) and commit the configs with it.

    python3 workflow/scripts/make_standing_configs.py            # writes the nine files
    python3 workflow/scripts/make_standing_configs.py --check    # exit 1 when a file differs

The families: shear2D (cases/davof/isoShear2D, the sheared disc of Rider and Kothe 1998 on
a one-cell-thick mesh, T = 8), shear3D (cases/davof/isoShear3D, the sphere in the shear
flow of Liovic et al. 2006, T = 3), deformation3D (cases/davof/sphereDeformation3D, the
LeVeque / Enright deformation, T = 3); the first two are the fork agent's templates
(origin/davof-swept, tokenised on 2026-11-07), the third the Enright agent's. Rows: the
DAVOF swept-face advection in its reference setting (davof-planeSurf-ringArea: the ring
normal, the plane surface, sweptFaceCells, the round-4 keys) and OpenFOAM v2606's
isoAdvector + plicRDF (plicRDF). OpenFOAM Co 0.2, 0.5, 1 (TARGET_ADVECTIVE_CO with the
case's REFERENCE_SPEED, which makes the target OpenFOAM's Courant number at t = 0 to within
a few per cent; the table reports the measured one). Ranks per rung fixed, the same
decomposition for every row (scotch, np of the config). Three repetitions per cell
(TPF_REPETITION 1, 2, 3: identical runs, each its own case and stamp; the table takes the
median wall clock and checks the three results are identical). The DAVOF rows run with
grazingCut fork (Tomislav's policy: the numbers are reproduced against the fork track's
record, which runs the same templates); the key is in every stamp and caption.
"""
import os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(REPO, "config", "davof", "standing")

FAMILIES = {
    # family: (case, setfields app on the fork route (the Snakefile maps the names on the leia route), rungs: N -> (np, runtime min, mem MB per cpu))
    "shear2D": ("davof/isoShear2D", "initDavofFields && initVortex2D",
                {64: (1, 240, 3000), 128: (1, 600, 3000), 256: (8, 900, 3000)}),
    "shear3D": ("davof/isoShear3D", "initDavofFields && initShearFlow3D",
                {32: (1, 240, 4000), 64: (8, 600, 4000), 128: (64, 900, 4000)}),
    "deformation3D": ("davof/sphereDeformation3D", "initDavofFields",
                      {32: (1, 240, 4000), 64: (8, 600, 4000), 128: (64, 900, 4000)}),
}
ROWS = ["davof-planeSurf-ringArea", "plicRDF"]
CO = [0.2, 0.5, 1]
REPS = [1, 2, 3]
GRAZING = "fork"

PREREG = {
    "shear2D": """# PRE-REGISTERED (2026-11-07, before any run; the swept-face advection has not run a 2D
# case before: the one-cell-thick mesh with two empty faces per cell is new ground for the
# swept flux, the explicit plane position was closed over the empty faces on 2026-11-07).
#   S2-1 Completion. plicRDF completes every cell bounded to 1e-9 (the fork agent's S-1).
#        The DAVOF row completes N 64 and 128 at OF Co 0.2 and 0.5 bounded to 1e-9 and
#        with |mass deviation| < 1e-10 (the swept flux is conservative by construction);
#        at OF Co 1 and at N 256 it may abort by the interface-cell guard (the fragmentation
#        of the thin filament at T/2: the 3D deformation aborted at N >= 64 on 2026-11-06).
#        Falsified if the DAVOF row aborts at N 64 OF Co 0.2 or loses the bound beyond 1e-6
#        at a cell it completes.
#   S2-2 Errors where both complete: E_shape(DAVOF) / E_shape(plicRDF) in [0.5, 1.5] at
#        every completed cell; L1(n)(DAVOF) / L1(n)(plicRDF) in [0.5, 2]. Falsified by a
#        ratio outside [0.3, 3].
#   S2-3 Orders over N (64 -> 128 -> 256) of E_shape between 1.5 and 2.5 for plicRDF at
#        every Co; for the DAVOF row over the completed rungs >= 1.5. Falsified by an
#        order below 1 on completed rungs.
#   S2-4 Wall clock per step at equal ranks: DAVOF / plicRDF in [1, 4] (the swept pieces
#        cost 1.5 to 3.6 x on the 3D cells of 2026-11-06; a one-cell-thick mesh has fewer
#        pieces per cell). Falsified by a ratio above 6.
#   S2-5 The three repetitions of a cell give identical results (every number of the
#        converter's last row), the wall clock within 5 % of the median. Falsified by a
#        differing result (the run would not be deterministic).
#   Also reported, not predicted: L1(alphaf) of the DAVOF row, the interface-cell counts
#   at the writes, the measured OF Courant number, the stamps.
""",
    "shear3D": """# PRE-REGISTERED (2026-11-07, before any run; the swept-face advection on the 3D shear
# flow is new: the fork agent's plane-clipped loop aborted on it at N >= 64 for OF Co >= 0.5
# on 2026-11-06, the swept-face advection has run the rotation, the translations and the
# Enright deformation).
#   S3-1 Completion. plicRDF completes every cell bounded to 1e-6 at most (it loses the
#        1e-9 bound where the Enright case did, N >= 64). The DAVOF row completes N 32 at
#        every Co and N 64 at OF Co 0.2 bounded to 1e-9 with |mass deviation| < 1e-10; it
#        aborts by the interface-cell guard at N 128 for OF Co >= 0.5 (the sheet of the
#        deformation case aborted at N >= 64 on 2026-11-06; the shear strains less).
#        Falsified if the DAVOF row aborts at N 32 OF Co 0.2 or completes N 128 at OF Co 1
#        bounded to 1e-9.
#   S3-2 Errors where both complete: E_shape(DAVOF) / E_shape(plicRDF) in [0.5, 1.5];
#        L1(n) ratio in [0.5, 2]. Falsified by a ratio outside [0.3, 3].
#   S3-3 Orders over N of E_shape for plicRDF between 1 and 2.5 at every Co; for the DAVOF
#        row over the completed rungs >= 1. Falsified by an order below 0.5.
#   S3-4 Wall clock per step at equal ranks: DAVOF / plicRDF in [1.5, 4] at N <= 64 and in
#        [1, 3] at N 128 on 64 ranks (the 2026-11-06 cells: 1.46 to 3.62 x). Falsified by
#        a ratio above 6.
#   S3-5 The three repetitions of a cell give identical results, the wall clock within 5 %
#        of the median (on 8 and 64 ranks the identity holds as it did for the Enright
#        cells with grazingCut fork). Falsified by a differing result.
#   Also reported, not predicted: L1(alphaf), the interface-cell and wisp counts at the
#   writes, the measured OF Courant number, the stamps.
""",
    "deformation3D": """# PRE-REGISTERED (2026-11-07, before any run; the Enright record of 2026-11-06, REPORT-swept
# "The 3D Enright deformation" and "Phase 3, the Enright identity", is the evidence).
#   D-1  Completion. plicRDF completes every cell (bounded to 1e-6 at N 32, unbounded to
#        1e-3 at N 128 OF Co 1 as on 2026-11-06). The DAVOF row with grazingCut fork
#        reproduces the fork's record bit for bit where the fork ran the same cell (N 32:
#        ends at T = 3 at OF Co 0.5 and 1; N 64: aborts at step 628 / 399 for OF Co 0.5 / 1;
#        N 128: aborts at 125 / 106) and completes N 32 at OF Co 0.2 bounded to 1e-9.
#        Falsified if an N 32 cell aborts or an N >= 64 cell completes.
#   D-2  Errors at N 32 (the completed rungs): E_shape(DAVOF) / E_shape(plicRDF) in
#        [1.2, 1.8] (2026-11-06: 1.137e-2 / 7.765e-3 = 1.46 at OF Co 0.5, 1.56 at OF Co 1),
#        L1(n) ratio in [1.5, 2.5]. Falsified by a ratio outside [1, 3].
#   D-3  Wall clock per step at equal ranks, N 32: DAVOF / plicRDF in [20, 40] (2026-11-06:
#        31 x at OF Co 0.5, 27 x at OF Co 1: the sheet's many interface cells and the swept
#        pieces of every one). Falsified by a ratio below 10 or above 60.
#   D-4  The three repetitions identical, the wall clock within 5 % of the median.
#   Also reported, not predicted: L1(alphaf), the interface-cell and wisp counts at the
#   writes (the fragmentation), the measured OF Courant number, the abort steps, the stamps.
""",
}


def config_text(family, N):
    case, setfields, rungs = FAMILIES[family]
    np_, runtime, mem = rungs[N]
    mode = "parallel" if np_ > 1 else "serial"
    rows = ", ".join(ROWS)
    return f"""# The standing study matrix (Tomislav, 2026-11-07; docs/davof/PORT-DESIGN.md 5b;
# workflow/scripts/make_standing_configs.py writes this file, edit the script): the
# family {family} ({case}) at N {N} on {np_} rank(s); the rows {rows}; OpenFOAM Co
# {', '.join(str(c) for c in CO)}; three identical repetitions per cell (TPF_REPETITION);
# the DAVOF rows with grazingCut {GRAZING} (Tomislav's policy of 2026-11-07: the numbers
# are reproduced against the fork track's record, which runs the same templates).
# Recompute only on change: workflow/scripts/case_stamp.py (README "Recompute only on
# change"); the tables: workflow/scripts/make_standing_tables.py {family}.
#
#   sbatch --parsable -J leia-standing-{family}-N{N} --export=ALL,STUDY=davof/standing/{family}_N{N} run-studies.sbatch
{PREREG[family]}study_name: standing_{family}_N{N}
studies_dir: studies/davof
case: {case}
mesh: hex
mode: {mode}
np: {np_}

solver: davofAdvectorFoam
theme:  davof

setfields_app: "{setfields}"
setfields_args: ""
post_solve: "python3 {{REPO}}/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv"

# The YAML booleans are QUOTED: they are rendered verbatim into the OpenFOAM
# dictionaries (false, not False).
axes_override:
  N_CELLS: [{N}]
  TPF_ROW: [{rows}]
  TARGET_ADVECTIVE_CO: [{', '.join(str(c) for c in CO)}]
  TPF_REPETITION: [{', '.join(str(r) for r in REPS)}]
  TPF_GRAZING_CUT: [{GRAZING}]
collapse_other_axes: true

solve_mem_mb: {mem}
solve_runtime: {runtime}

export_slides: false

env_preamble: "source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
"""


def main(argv):
    check = "--check" in argv
    os.makedirs(OUT, exist_ok=True)
    bad = 0
    for family, (case, setfields, rungs) in FAMILIES.items():
        for N in rungs:
            path = os.path.join(OUT, f"{family}_N{N}.yaml")
            text = config_text(family, N)
            old = open(path).read() if os.path.isfile(path) else None
            if check:
                if old != text:
                    print("differs:", os.path.relpath(path, REPO)); bad += 1
                continue
            if old != text:
                with open(path, "w") as fh:
                    fh.write(text)
                print("written:", os.path.relpath(path, REPO))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
