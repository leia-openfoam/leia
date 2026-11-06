#!/usr/bin/env python3
"""Write the study configs of the ring normal's stencil variants (2026-11-08, Tomislav's
task: the stability and the accuracy of the ring normal's stencil in the swept-face
advection; src/davof/normalModels/ringAreaNormal.H, tpf_rows.NORMAL_STENCILS):
config/davof/stencil_*.yaml. The predictions are registered in the headers before any
run; the tables: workflow/scripts/davof_stencil_tables.py.

    python3 workflow/scripts/make_stencil_configs.py
"""
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(REPO, "config", "davof")

VARIANTS = ["production", "twoRings", "kernelGaussian", "kernelWendland", "angle30",
            "angle60", "planeFitMean", "planeFitOwn", "ringQuadric"]

COMMON = """# THE RING NORMAL'S STENCIL (2026-11-08; workflow/scripts/make_stencil_configs.py writes this
# file, edit the script). The normal of cell k is the ring area normal d_k = sum over the ring
# R(k) of m_j, m_j = -sum_f alpha_f S_f (|m_j| = |Sigma_j| for a plane), n_k = d_k/|d_k|: the
# area-weighted average of the ring's unit normals, stable because the division is by the
# ring's total patch area (the own normal m_k/|m_k| grew 1.16 per step on the translated
# plane, tilt diagnostic of 2026-10-26). The variants, one word TPF_NORMAL_STENCIL
# (tpf_rows.NORMAL_STENCILS, the opt-in keys of normalModel ringArea; production = the block
# empty, bit for bit the model of every study before): twoRings (the union of the rings of
# the ring's cells, area weights), kernelGaussian / kernelWendland (one ring, w_j = |Sigma_j|
# K(|x_j - x_k|/h_k) with the cell centres; Gaussian exp(-q^2) cut at 2 h, Gingold and
# Monaghan 1977 doi:10.1093/mnras/181.3.375; Wendland C2 (1 - q/2)^4 (1 + 2 q), support 2 h,
# Wendland 1995 doi:10.1007/BF02123482), angle30 / angle60 (a ring cell whose own normal
# deviates from d_k by more than theta_max is dropped), planeFitMean / planeFitOwn (the
# weighted least-squares plane through the ring's patch centroids, the eigenvector of the
# smallest eigenvalue of the weighted covariance, Pearson 1901
# doi:10.1080/14786440109462720, Hoppe et al. 1992 doi:10.1145/133994.134011; own: centred
# on the cell's centroid with the moment-corrected weights sum_j w'_j (x_j - x_k) = 0, Liu,
# Jun, Zhang 1995 doi:10.1002/fld.1650200824), ringQuadric (the ring quadric of the
# preprint, normalModel ringFit). Every variant on the swept-face advection in the
# recommended setting (TPF_ENRIGHT_SETTING recommended: signedTets + scale + normalAndAlpha)
# with grazingCut fork (the identity policy of 2026-11-07), the four bit-for-bit keys of
# round 4, write precision 17.
"""

SWEPT_KEYS = """  TPF_ENRIGHT_SETTING: [recommended]
  TPF_GRAZING_CUT: [fork]
  TPF_VOLUME_RULE: [kromerBothe]
  TPF_ANCHOR_METHOD: [ccs]
  TPF_FLUX_INTEGRATION: [sweptFaceCells]
  TPF_SWEPT_SCALE: ["true"]
  TPF_KB_FAN: ["auto"]
  TPF_SWEPT_AXIS_PLANES: ["true"]
  TPF_SWEPT_TET_CACHE: ["true"]
  TPF_SWEPT_POINT_VELOCITY_LAZY: ["true"]
  TPF_SWEPT_PRE_TEST_COPLANAR_SIDE: ["true"]
  TPF_WRITE_PRECISION: ["17"]
"""

TAIL = """collapse_other_axes: true

solve_mem_mb: {mem}
solve_runtime: {runtime}

export_slides: false

env_preamble: "source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
"""


def body(name, case, mode, np_, setfields_args, axes, mem, runtime, extra_post=""):
    s = f"""study_name: {name}
studies_dir: studies/davof
case: davof/{case}
mesh: hex
mode: {mode}
np: {np_}

solver: davofAdvectorFoam
theme:  davof

setfields_app: initDavofFields
setfields_args: "{setfields_args}"
post_solve: "python3 {{REPO}}/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv"

# The YAML booleans are QUOTED: they are rendered verbatim into the OpenFOAM
# dictionaries (false, not False).
axes_override:
{axes}"""
    s += TAIL.format(mem=mem, runtime=runtime)
    return s


def write(name, header, text):
    path = os.path.join(OUT, name + ".yaml")
    with open(path, "w") as fh:
        fh.write(COMMON + header + "#\n#   sbatch --parsable -J leia-davof-" + name
                 + " --export=ALL,STUDY=davof/" + name
                 + ",PROFILE=profiles/slurm-small run-studies.sbatch\n" + text)
    print("[stencil]", path)


VAR = "  TPF_NORMAL_STENCIL: [" + ", ".join(VARIANTS) + "]\n"

# ---------------------------------------------------------------------------- test 1
PLANE_AXES = ("""  N_CELLS: [40, 80]
  TARGET_ADVECTIVE_CO: [0.2, 0.5]
  TPF_VELOCITY_ARM: [normal, oblique]
  TPF_ROW: [davof-planeSurf-ringArea]
""" + VAR + SWEPT_KEYS + """  TPF_RECON_ERROR_EVERY_STEP: ["true"]
  TPF_PROFILING: ["false"]
""")

GATE_HEADER = """# TEST 1a, THE PLANE GATE: the translated plane (cases/davof/planeTranslation3D, the exact
# Detrixhe-Aslam state, 0.2 L of travel = 0.2 N / Co steps) at N 40 / 80, OpenFOAM Co 0.2 /
# 0.5, both arms (normal: U = n, the sliver-amplifying arm; oblique: U = (1 0 0)), every
# variant: 2 x 2 x 2 x 9 = 72 serial cases. Read-out: the run maxima of L1(alphaf) and L1(n)
# (areaFractionError and reconstructionError every step), min alpha, the excess over 1, the
# mass deviation, the wall clock.
#
# PRE-REGISTERED (2026-11-08, before any run):
# G-0 (identity): the production rows reproduce the port's record bit for bit where the same
#     cell exists (port_enrightKeys_plane_leia, setting recommended, N 40, Co 0.5, both arms:
#     L1(alphaf) and L1(n) at every step equal to the digit). Falsified by any difference.
# G-1 (the gate): every variant keeps the plane at round-off: the run maxima of L1(alphaf)
#     and L1(n) below 1e-11 in all 72 cases. A plane is reproduced exactly by every stencil
#     (the sum of parallel area normals, a kernel-weighted sum, a filtered sum, the plane
#     through coplanar centroids, the quadric through coplanar planes). Falsified by any
#     case above 1e-11; the expected weak spots are planeFitOwn (the moment-corrected
#     weights can be negative) and ringQuadric (the Monge frame of m_k of a sliver cell).
# G-2 (identity of the angle filters on the plane): angle30 and angle60 drop no cell on the
#     plane (every ring normal is parallel to d_k) and reproduce production bit for bit.
"""

GROWTH_HEADER = """# TEST 1b, THE SEEDED PERTURBATION (the tilt-diagnostic approach of 2026-10-26 with a seed
# above round-off): the same 72 cells with the initial face fractions perturbed by a random
# field of amplitude 1e-8 on the partial faces (workflow/scripts/davof_perturb_alphaf.py,
# seed 20261108, the cell fractions exact), END_TIME 1e-3 (50 steps at N 40 Co 0.2, 100 at
# N 80 Co 0.2, 20 at N 40 Co 0.5, 40 at N 80 Co 0.5; 0.25 L of travel, the plane inside
# the corner, n.x <= 1.58 L < 1.6124 L). Read-out per case from the per-step L1(alphaf -
# alphaf_exact): the early factor g_early = (L1 at step 6 / L1 at step 1)^(1/5), the late
# factor g_late = (L1 at the last step / L1 at the middle step)^(1/(steps/2)), the end
# value; the growth is the tilt feedback alone (the perturbation enters the loop only
# through m_k = -sum_f alpha_f S_f and the planes positioned along the normal).
#
# PRE-REGISTERED (2026-11-08, before any run):
# S-1 (contraction): production, twoRings, kernelGaussian, angle30, angle60, planeFitMean
#     and ringQuadric contract: g_early < 1 in every cell and the end value at the gate's
#     floor (below 1e-11) at Co 0.2; at Co 0.5 g_early < 1 and the end value below 1e-9
#     (fewer steps). The ring sum's weights are bounded by 1/9 per cell (2026-10-26), the
#     wider and the smoother stencils damp more. Falsified by g_early >= 1 or an end value
#     above the seed (1e-8) for any of these.
# S-2 (the Wendland kernel is the own normal in disguise): with support 2 h the weights of
#     the face, edge and corner neighbours are 0.19, 0.03 and 0.001 of the cell's own, so
#     kernelWendland sits between the own normal (growth 1.16 per step at Co 0.2, 1.5-8.8 at
#     Co 0.5) and production: it GROWS on the normal arm at Co 0.5 (g_early > 1 at N 40 and
#     80) and contracts or stalls at Co 0.2. Falsified if it contracts in every cell.
# S-3 (the centred fit): planeFitOwn contracts at Co 0.2 (the eigenvector of a plane's
#     covariance is insensitive to a 1e-8 tilt of one point); its Co 0.5 behaviour on the
#     normal arm is the open question (the corrected weights are not bounded by 1), reported
#     as measured, not predicted.
# S-4 (the ranking is the stencil width): among the contracting one-ring variants the early
#     factor orders kernelGaussian > production > twoRings (the wider the stencil, the
#     smaller the factor). Falsified if twoRings' g_early exceeds production's in more than
#     half of the cells.
"""

write("stencil_plane_gate", GATE_HEADER,
      body("stencil_plane_gate", "planeTranslation3D", "serial", 1, "", PLANE_AXES, 4000, 60))
write("stencil_plane_growth", GROWTH_HEADER,
      body("stencil_plane_growth", "planeTranslation3D", "serial", 1,
           "; python3 {REPO}/workflow/scripts/davof_perturb_alphaf.py 0/alphaf.water 1e-8 20261108",
           "  END_TIME: [1e-3]\n" + PLANE_AXES, 4000, 90))

# ---------------------------------------------------------------------------- test 2
SPHERE_HEADER = """# TEST 2, THE ACCURACY ON THE TRANSLATED SPHERE (cases/davof/sphereTranslation3D along x,
# exactSphere state, R 0.15, one period L = 1 at Co 0.2 = N / 0.2 steps) at N {N} on {np}
# rank(s) ({mode}), every variant, three repetitions (TPF_REPEAT, the timing): 27 cases.
# Read-out at the end: L1(n) against the exact normal (reconstructionError at the end), its
# least-squares order over N 32 / 64 / 128 (the three studies stencil_sphereN32/64/128
# together, docs/davof tables by davof_stencil_tables.py), L1(alphaf), L1(x), E_shape
# (the fixed one against the exact sphere), min alpha, the excess over 1, the mass
# deviation, the interface-cell count, and the wall clock at equal ranks (the solve's
# total, three repetitions, median and spread; the per-step median from the log; nodes
# shared, not exclusive: the spread is reported with it).
#
# PRE-REGISTERED (2026-11-08, before any run; the record: the ring normal's L1(n) order 1.00
# in the paraboloid loop, 1.35 per period for the fractions; the ring quadric 1.52 at
# 1.6e-4, fifty times below plicRDF's; in the swept loop the ring row's L1(n) order is not
# on record):
# A-1 (the quadric is the accurate stencil): ringQuadric has the lowest L1(n) at every N
#     and an order >= 1.5 over 32 -> 128; production's order lies in [0.8, 1.2].
#     Falsified by ringQuadric above production at any N or an order below 1.2.
# A-2 (the area-weighted stencils are first order, whatever their width): twoRings,
#     kernelGaussian and kernelWendland have L1(n) orders in [0.8, 1.2] and end values
#     within a factor 2 of production's (the error is the asymmetry of the patch about the
#     cell, O(h/R), which a wider or a sharper stencil moves but does not remove).
#     Falsified by an order outside [0.6, 1.4] or a factor above 3.
# A-3 (the angle filters are inert on a resolved sphere): the normals of neighbouring
#     cells differ by about h/R (12 degrees at N 32, 3 at N 128), below both thresholds,
#     so angle30 and angle60 drop no cell and reproduce production bit for bit at every
#     N. Falsified by any dropped cell or any digit.
# A-4 (the centred plane fit): planeFitOwn's L1(n) order is >= 1.5 (the centring condition
#     cancels the first-order asymmetry of the patch about x_k; the eigenvector then sees
#     the second-order curvature term only), planeFitMean's stays in [0.8, 1.2]. This is
#     the hypothesis behind the centring condition; falsified by planeFitOwn's order
#     below 1.3.
# A-5 (cost): the per-step median of every one-ring variant is within 1.1 x production's,
#     twoRings within 1.3 x, planeFitMean / planeFitOwn / ringQuadric within 1.5 x (the
#     two positioning passes). Falsified by a ratio above 2 for any variant.
# A-6 (robustness): every variant completes every case bounded to 1e-9 (min alpha
#     >= -1e-9, excess <= 1e-9) and conserves the volume to 1e-10 %.
"""

SPHERE_AXES = """  N_CELLS: [{N}]
  TARGET_ADVECTIVE_CO: [0.2]
  TPF_ROW: [davof-planeSurf-ringArea]
  TPF_REPEAT: [1, 2, 3]
""" + VAR + SWEPT_KEYS + """  TPF_AFE_CONSISTENCY: ["false"]
"""

for N, np_, mode, mem, rt in ((32, 1, "serial", 4000, 60), (64, 8, "parallel", 3000, 120),
                              (128, 64, "parallel", 3000, 240)):
    name = f"stencil_sphereN{N}"
    write(name, SPHERE_HEADER.format(N=N, np=np_, mode=mode),
          body(name, "sphereTranslation3D", mode, np_, "", SPHERE_AXES.format(N=N), mem, rt))

# ---------------------------------------------------------------------------- test 3
ROT_HEADER = """# TEST 3a, THE ROTATION (cases/davof/sphereRotation3D: R 0.15 turned once about the axis
# through (0.75 . 0.5), omega 2 pi, END_TIME 1, Co 0.2 at the far point = 400 / 800 steps,
# exactSphere state, alpha fixedValue) at N 32 / 64 serial, every variant, three
# repetitions: 54 cases. Read-out at the end of the turn: bounded (min alpha, excess),
# conservative (mass deviation), E_shape (fixed), L1(x), L1(n), L1(alphaf), the
# interface-cell count, the wall clock (median of three).
#
# PRE-REGISTERED (2026-11-08, before any run; the record of the production row on the
# rotation: bounded, mass 1e-11 %, E_shape order 1.86, L1(x) 1.90):
# R-1 (robustness): every variant completes both rungs bounded to 1e-9 and conservative to
#     1e-10 %. Falsified by an abort or an excursion.
# R-2 (the ranking follows the normal): the variants' E_shape and L1(x) at N 64 order as
#     their L1(n) (Spearman rank correlation >= 0.7 over the nine variants); ringQuadric's
#     E_shape is below production's. Falsified by a correlation below 0.4.
# R-3 (identity): angle30 / angle60 reproduce production bit for bit (no cell dropped on a
#     sphere turned about an axis, as A-3).
"""

ROT_AXES = """  N_CELLS: [32, 64]
  TPF_ROW: [davof-planeSurf-ringArea]
  TARGET_ADVECTIVE_CO: [0.2]
  TPF_VELOCITY_ARM: [rotation]
  TPF_ALPHA_BC: [fixedValue]
  TPF_DAVOF_METHOD: [exactSphere]
  TPF_REPEAT: [1, 2, 3]
""" + VAR + SWEPT_KEYS + """  TPF_AFE_CONSISTENCY: ["false"]
  TPF_SWEPT_FACE_UPDATE: [pieces]
  TPF_SWEPT_COPLANAR: [donorSide]
  TPF_SWEPT_CAP_SET: [final]
  TPF_PLANE_SOURCE: [anchor]
"""

write("stencil_rotation", ROT_HEADER,
      body("stencil_rotation", "sphereRotation3D", "serial", 1, "", ROT_AXES, 8000, 240))

DEF_HEADER = """# TEST 3b, THE 3D DEFORMATION (cases/davof/isoDeformation3D, the standing template: LeVeque
# / Enright, sphere R 0.15 at (0.35 0.35 0.35), T = 3, OpenFOAM Co 0.5 = 480 steps) at N 32
# serial, every variant, three repetitions (TPF_REPETITION): 27 cases. Read-out at T against
# the initial sphere: bounded, conservative, E_shape, L1(n), L1(x), the interface-cell and
# wisp counts at T/2 and T, the wall clock; the variant that gains accuracy on the sphere
# but loses robustness on the sheet is seen here.
#
# PRE-REGISTERED (2026-11-08, before any run; the record of the recommended setting at
# N 32: completes bounded to 1e-7 after the bounding, E_shape 1.2-1.3 x plicRDF's):
# D-1 (robustness): production, kernelGaussian, kernelWendland, angle30, angle60 and
#     planeFitMean complete bounded to 1e-6. twoRings and ringQuadric are the risks: the
#     two-ring sum mixes the two sides of the sheet where it is thinner than 4 h (the
#     normals cancel, |d_k| collapses) and the quadric through a folding sheet is not a
#     Monge patch; prediction: both complete, twoRings with an E_shape above production's
#     (>= 1.2 x), ringQuadric with more fit fallbacks than the sphere. Falsified for the
#     first group by an abort, for twoRings by an E_shape below production's.
# D-2 (the angle filters act on the sheet): angle30 drops cells at T/2 (the folded sheet
#     holds neighbouring normals beyond 30 degrees) and differs from production in E_shape
#     by more than 1 %; angle60 drops fewer than a tenth of angle30's cells.
# D-3 (the centred fit on the sheet): planeFitOwn's E_shape is within [0.8, 1.2] of
#     production's (a plane fit through a thin sheet's centroids sees both sides alike).
"""

DEF_AXES = """  N_CELLS: [32]
  TPF_ROW: [davof-planeSurf-ringArea]
  TARGET_ADVECTIVE_CO: [0.5]
  TPF_REPETITION: [1, 2, 3]
""" + VAR + SWEPT_KEYS

write("stencil_deformation", DEF_HEADER,
      body("stencil_deformation", "isoDeformation3D", "serial", 1, "", DEF_AXES, 4000, 360))
