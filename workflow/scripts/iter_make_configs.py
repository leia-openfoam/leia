#!/usr/bin/env python3
"""The study configs of the iterative face-fraction corrector (leia davof-iter,
2026-11-09; src/davof/reconstruction/davof.H "The iterative face-fraction corrector",
the pre-registration in docs/davof/REPORT-swept.md "The iterative face-fraction
corrector"). Writes config/davof/iter_*.yaml; edit this script, not the files.

    python3 workflow/scripts/iter_make_configs.py

The rows: davof-planeSurf-ringArea (the DAVOF reference setting of the port: the swept
flux with the round-4 keys, the ring area normal, CCS positioning) with the corrector
none | mean1 | mean2 | mean4 (TPF_FACE_FRACTION_CORRECTOR: mean<k> = mean with k
iterations, materialize.py). The sphere rotation (sphereRotation3D: one turn, OF Co 0.2
= interface Co 0.0463, OF Co 1 = 0.2315) at N 32 / 64 / 128 on 1 / 8 / 64 ranks (the
standing protocol's ranks), also under OpenFOAM's isoAdvection with the DAVOF
reconstruction (TPF_ADVECTION_SCHEME isoAdvector; isoAdvectorScheme.H) at N 32 / 64;
the plane translation (planeTranslation3D, both arms, Co 0.2 / 0.5 / 1, N 20 / 40 / 80)
as the gate; the 3D deformation (isoDeformation3D, the standing matrix's template and
tokens, grazingCut fork so that the none row reproduces the standing record) at N 32 /
64, OF Co 0.5; the np 1 / np 8 identity pairs (20 steps, writePrecision 17) on the plane
and the rotation.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "config", "davof")

# the swept keys of the port's reference setting (port_cmpRotSweptN32_leia.yaml,
# port_enrightKeys_plane_leia.yaml): the round-4 keys on, the fork's reference defaults
SWEPT = {
    "TPF_VOLUME_RULE": ["kromerBothe"],
    "TPF_ANCHOR_METHOD": ["ccs"],
    "TPF_FLUX_INTEGRATION": ["sweptFaceCells"],
    "TPF_SWEPT_SCALE": ['"true"'],
    "TPF_SWEPT_FACE_UPDATE": ["pieces"],
    "TPF_SWEPT_COPLANAR": ["donorSide"],
    "TPF_SWEPT_CAP_SET": ["final"],
    "TPF_PLANE_SOURCE": ["anchor"],
    "TPF_KB_FAN": ['"auto"'],
    "TPF_SWEPT_WISP_SURFACE": ["matched"],
    "TPF_SWEPT_WISP_TOL": ["0"],
    "TPF_SWEPT_CONSISTENT": ['"true"'],
    "TPF_SWEPT_WISP_MATCH_TOL": ["1e-9"],
    "TPF_SWEPT_VOLUME_CORRECTION": ["scale"],
    "TPF_SWEPT_AXIS_PLANES": ['"true"'],
    "TPF_SWEPT_TET_CACHE": ['"true"'],
    "TPF_SWEPT_POINT_VELOCITY_LAZY": ['"true"'],
    "TPF_SWEPT_PRE_TEST_COPLANAR_SIDE": ['"true"'],
}
ROT = {
    "TPF_ROW": ["davof-planeSurf-ringArea"],
    "TPF_VELOCITY_ARM": ["rotation"],
    "TPF_ALPHA_BC": ["fixedValue"],
    "TPF_DAVOF_METHOD": ["exactSphere"],
    "TPF_AFE_CONSISTENCY": ['"false"'],
    "TPF_RECON_ERROR_EVERY_STEP": ['"false"'],
    "TPF_PROFILING": ['"false"'],
    "TPF_WRITE_PRECISION": ['"17"'],
    "END_TIME": ["1"],
    "TPF_WRITE_INTERVAL": ["0.25"],
}
CO_OF02 = "0.0462996198086"    # the interface Co of OpenFOAM Co 0.2 on the rotation
CO_OF1 = "0.231498099043"      # OpenFOAM Co 1

HEAD = """# {purpose}
# Written by workflow/scripts/iter_make_configs.py (leia davof-iter, 2026-11-09); edit the
# script. Pre-registration: docs/davof/REPORT-swept.md "The iterative face-fraction corrector".
#
#   sbatch --parsable -J leia-{name} --export=ALL,STUDY=davof/{name} run-studies.sbatch
study_name: {name}
studies_dir: studies/davof
case: davof/{case}
mesh: hex
mode: {mode}
np: {np}

solver: davofAdvectorFoam
theme:  davof

setfields_app: initDavofFields
setfields_args: ""
post_solve: "python3 {{REPO}}/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv"

# The YAML booleans and the numbers rendered into the dictionaries are QUOTED.
axes_override:
{axes}collapse_other_axes: true

solve_mem_mb: {mem}
solve_runtime: {runtime}

export_slides: false

env_preamble: "source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
"""


def write(name, case, np_, axes, purpose, runtime, mem=4000):
    lines = "".join("  %s: [%s]\n" % (k, ", ".join(v)) for k, v in axes.items())
    text = HEAD.format(purpose=purpose, name=name, case=case, mode=("serial" if np_ == 1 else "parallel"),
                       np=np_, axes=lines, mem=mem, runtime=runtime)
    path = os.path.join(OUT, name + ".yaml")
    with open(path, "w") as f:
        f.write(text)
    print("wrote", os.path.relpath(path, ROOT))


def main():
    corr = ["none", "mean1", "mean2", "mean4"]
    # the rotation: P2 (contraction), P3 (accuracy), P4 (cost), P5 (the first correction's order)
    for N, np_, cos, rt in ((32, 1, [CO_OF02, CO_OF1], 360), (64, 8, [CO_OF02, CO_OF1], 480),
                            (128, 64, [CO_OF02], 900)):
        axes = {"N_CELLS": [str(N)], "TARGET_ADVECTIVE_CO": cos, "TPF_FACE_FRACTION_CORRECTOR": corr}
        axes.update(ROT)
        axes.update(SWEPT)
        write("iter_rot_N%d" % N, "sphereRotation3D", np_, axes,
              "The corrector on the sphere rotation at N %d on %d rank(s): one turn, the DAVOF\n"
              "# reference setting (davof-planeSurf-ringArea, swept flux, round-4 keys), the corrector\n"
              "# none | mean with 1, 2, 4 iterations (P2-P5)." % (N, np_), rt)
    # the DAVOF reconstruction under OpenFOAM's isoAdvection (the predictor from the
    # translated planes, isoAdvectorScheme.H): the first-correction diagnostic of the
    # other flux
    for N, np_, rt in ((32, 1, 360), (64, 8, 480)):
        axes = {"N_CELLS": [str(N)], "TARGET_ADVECTIVE_CO": [CO_OF02],
                "TPF_FACE_FRACTION_CORRECTOR": ["none", "mean1", "mean4"],
                "TPF_ADVECTION_SCHEME": ["isoAdvector"]}
        axes.update(ROT)
        axes.update(SWEPT)
        write("iter_rotIso_N%d" % N, "sphereRotation3D", np_, axes,
              "The corrector under OpenFOAM's isoAdvection with the DAVOF reconstruction (the\n"
              "# predictor alphaf from the translated planes, isoAdvectorScheme.H) on the rotation at\n"
              "# N %d, OF Co 0.2: the first correction of the other flux (P5), the convergence (P2)." % N, rt)
    # the plane gate P1
    plane = {"N_CELLS": ["20", "40", "80"], "TARGET_ADVECTIVE_CO": ["0.2", "0.5", "1"],
             "TPF_VELOCITY_ARM": ["normal", "oblique"], "TPF_ROW": ["davof-planeSurf-ringArea"],
             "TPF_FACE_FRACTION_CORRECTOR": ["none", "mean4"],
             "TPF_WRITE_PRECISION": ['"17"'], "TPF_PROFILING": ['"false"']}
    plane.update(SWEPT)
    write("iter_plane", "planeTranslation3D", 1, plane,
          "The plane gate P1 of the corrector: both arms (U = n, U = (1 0 0)), Co 0.2 / 0.5 / 1,\n"
          "# N 20 / 40 / 80, the corrector none | mean4; the corrector must change nothing (0 or 1\n"
          "# iteration with max |dalphaf| <= 1e-12 every step; L1(alphaf), L1(n) <= 1e-12 as without).", 240)
    # the np 1 / np 8 identity pairs: 20 steps, every face fraction and alpha compared
    for np_ in (1, 8):
        axes = {"N_CELLS": ["40"], "TARGET_ADVECTIVE_CO": ["0.5"], "TPF_VELOCITY_ARM": ["oblique"],
                "TPF_ROW": ["davof-planeSurf-ringArea"], "TPF_FACE_FRACTION_CORRECTOR": ["mean4"],
                "END_TIME": ['"1e-3"'], "TPF_WRITE_INTERVAL": ['"5e-5"'],
                "TPF_WRITE_PRECISION": ['"17"'], "TPF_PROFILING": ['"false"']}
        axes.update(SWEPT)
        write("iter_par_plane_np%d" % np_, "planeTranslation3D", np_, axes,
              "The parallel identity of the corrector on the plane (N 40, oblique, Co 0.5, 20 steps =\n"
              "# 1e-3, dt 5e-5): np %d of the pair np 1 / np 8, alpha and alphaf at the end compared." % np_, 30)
        axes = {"N_CELLS": ["32"], "TARGET_ADVECTIVE_CO": [CO_OF02], "TPF_FACE_FRACTION_CORRECTOR": ["mean4"]}
        axes.update(ROT)
        axes.update(SWEPT)
        # 20 steps of the rendered dt 0.0005787452476 (= 0.0462996198086 (1/32)/2.5), every
        # step written so that the end state is on disk whatever the end-time rounding
        axes["END_TIME"] = ['"0.011574904952"']
        axes["TPF_WRITE_INTERVAL"] = ['"0.0005787452476"']
        write("iter_par_rot_np%d" % np_, "sphereRotation3D", np_, axes,
              "The parallel identity of the corrector on the rotation (N 32, OF Co 0.2, 20 steps):\n"
              "# np %d of the pair np 1 / np 8, alpha and alphaf at the end compared (the corrector\n"
              "# acts here, unlike on the plane)." % np_, 30)
    # the 3D deformation of the standing matrix (its template and tokens; grazingCut fork so
    # that the none row reproduces the standing record)
    for N, np_, rt in ((32, 1, 360), (64, 8, 720)):
        axes = {"N_CELLS": [str(N)], "TPF_ROW": ["davof-planeSurf-ringArea"], "TARGET_ADVECTIVE_CO": ["0.5"],
                "TPF_REPETITION": ["1"], "TPF_GRAZING_CUT": ["fork"], "TPF_FACE_FRACTION_CORRECTOR": corr}
        write("iter_def_N%d" % N, "isoDeformation3D", np_, axes,
              "The corrector on the 3D deformation (the standing matrix's isoDeformation3D at N %d on\n"
              "# %d rank(s), OF Co 0.5, grazingCut fork): none | mean 1, 2, 4 iterations (P3)." % (N, np_), rt)


if __name__ == "__main__":
    main()
