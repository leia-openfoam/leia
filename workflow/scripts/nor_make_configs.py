#!/usr/bin/env python3
"""The study configs of the normal-order study (leia davof-normal-order, 2026-11-12;
the pre-registration in docs/davof/REPORT-swept.md "Normal order versus normal bias").
Writes config/davof/nor_*.yaml; edit this script, not the files.

    python3 workflow/scripts/nor_make_configs.py

Three normals on identical initial states (initDavofFields exactSphere) and the same
swept-face scheme, the pieces path in the standing Co 1 setting (the round-4 keys, the
Enright setting recommended = sweptPieceVolume signedTets + sweptBoundedFaceFraction
scale + interfaceCellRule normalAndAlpha): davof-planeSurf-exact (the control),
davof-planeSurf-ringArea (the ring area normal d_k) and davof-planeSurf-youngs (the
Youngs normal from the Gauss linear gradient of alpha_k). Two tests: the sphere
translation (sphereTranslation3D, U = (1 0 0), one unit of travel; OpenFOAM Co =
TARGET_ADVECTIVE_CO for the uniform velocity) with the per-step flux error against the
exact flux (sweptExactFluxError), and the solid-body rotation (sphereRotation3D, one
turn; OpenFOAM Co 0.2 / 1 = interface Co 0.0463 / 0.2315). N 32 / 64 / 128 on 8 ranks,
OpenFOAM Co 0.2 and 1; N 32 also serial (the serial/parallel digit check).
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "config", "davof")

ROWS = ["davof-planeSurf-exact", "davof-planeSurf-ringArea", "davof-planeSurf-youngs"]

# the swept keys of the port's reference setting with the round-4 keys
# (port_cmpRotSweptN32_leia.yaml), the same on both tests
SWEPT = {
    "TPF_VOLUME_RULE": ["kromerBothe"],
    "TPF_ANCHOR_METHOD": ["ccs"],
    "TPF_ANCHOR_SKIP_CONVERGED": ['"true"'],
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
    "TPF_SWEPT_PRE_TEST": ["touching"],
    "TPF_SWEPT_STATE_SHORTCUT": ['"false"'],
    "TPF_SWEPT_TRACE_CACHE": ['"true"'],
    "TPF_SWEPT_AXIS_PLANES": ['"true"'],
    "TPF_SWEPT_TET_CACHE": ['"true"'],
    "TPF_SWEPT_POINT_VELOCITY_LAZY": ['"true"'],
    "TPF_SWEPT_PRE_TEST_COPLANAR_SIDE": ['"true"'],
    "TPF_POSITIONING_TOL": ["1e-14"],
    "TPF_GRAZING_CUT": ["openfoam"],
    "TPF_ENRIGHT_SETTING": ["recommended"],
    "TPF_YOUNGS_GRAD_SCHEME": ["Gauss_linear"],
    "TPF_DAVOF_METHOD": ["exactSphere"],
    "TPF_REF_METHOD": ["exactSphere"],
    "TPF_WRITE_PRECISION": ['"17"'],
    "TPF_PROFILING": ['"false"'],
    "TPF_RECON_ERROR_EVERY_STEP": ['"false"'],
    "TPF_AFE_CONSISTENCY": ['"false"'],
    "TPF_AFE_WRITE_CONTROL": ["adjustableRunTime"],
    "TPF_AFE_WRITE_INTERVAL": ["0.25"],
    "TPF_VFE_WRITE_CONTROL": ["adjustableRunTime"],
    "TPF_VFE_WRITE_INTERVAL": ["0.25"],
    "TPF_WRITE_INTERVAL": ["0.25"],
    "END_TIME": ["1"],
}
TRANS = {"TPF_VELOCITY_ARM": ["axial"], "TARGET_ADVECTIVE_CO": ["0.2", "1"],
         "TPF_EXACT_FLUX_ERROR": ['"true"']}
ROT = {"TPF_VELOCITY_ARM": ["rotation"], "TPF_ALPHA_BC": ["fixedValue"],
       "TARGET_ADVECTIVE_CO": ["0.0462996198086", "0.231498099043"],
       "TPF_EXACT_FLUX_ERROR": ['"false"']}

HEAD = """# {purpose}
# Written by workflow/scripts/nor_make_configs.py (leia davof-normal-order, 2026-11-12);
# edit the script. Pre-registration: docs/davof/REPORT-swept.md "Normal order versus
# normal bias".
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


def write(name, case, np_, axes, purpose, runtime, mem):
    lines = "".join("  %s: [%s]\n" % (k, ", ".join(v)) for k, v in axes.items())
    text = HEAD.format(purpose=purpose, name=name, case=case,
                       mode=("serial" if np_ == 1 else "parallel"), np=np_, axes=lines,
                       mem=mem, runtime=runtime)
    with open(os.path.join(OUT, name + ".yaml"), "w") as f:
        f.write(text)
    print("wrote config/davof/%s.yaml" % name)


def main():
    for test, case, extra in (("trans", "sphereTranslation3D", TRANS), ("rot", "sphereRotation3D", ROT)):
        for N, np_, rt, mem in ((32, 8, 240, 4000), (64, 8, 720, 4000), (128, 8, 1440, 8000), (32, 1, 600, 6000)):
            axes = {"N_CELLS": [str(N)], "TPF_ROW": ROWS}
            axes.update(extra)
            axes.update(SWEPT)
            name = "nor_%s_N%d%s" % (test, N, "_np1" if np_ == 1 else "")
            write(name, case, np_, axes,
                  "The normal-order study, %s at N %d on %d rank(s): the exact, ring area and\n"
                  "# Youngs normals on the pieces path (Enright setting recommended), OpenFOAM Co 0.2\n"
                  "# and 1." % ({"trans": "the sphere translation", "rot": "the sphere rotation"}[test], N, np_),
                  rt, mem)


if __name__ == "__main__":
    main()
