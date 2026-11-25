# Small transport test of the face-fraction corrector with the ring linear fit: tables (branch davof-iter-linfit)

Executor: Claude (a sub-session executing the designing session's plan). The pre-registration is the section "Pre-registration of the small transport test" of `REPORT-iterlin.md` (commit 17f76e01, before any transport run whose result was read). The PP gate and one probe preceded it (allowed by the order of work). Eight transport jobs (57253132-57253139) were submitted by mistake one command before the pre-registration commit; they were cancelled before any output was read and the matrix was resubmitted after the commit (57253296-57253303).

## Setup

- Case: leia `cases/davof/isoShear3D` rendered by leia's own `materialize.py` through `il_render.py` (route davofAdvectorFoam, `TPF_ADVECTION_SCHEME isoAdvector`, hexahedral 1 x 1 x 2 box with N x N x 2N cells, the sphere of radius 0.15 at (0.5 0.75 0.25), period T = 3 = END_TIME, fixed time step from the target OpenFOAM Co, `TPF_ALPHA_BC zeroGradient`), the initial alpha and alpha_f from davofInitFields, the flow from davofInitShearFlow3D. N 32 on 1 rank, N 64 on 8 ranks (scotch, one node). OpenFOAM Co 0.2 and 0.5.
- PLIC rows: `TPF_ROW plicRDF` (plicRDF, vofCutCell, normalModel cellLocal key inert). davof rows: `TPF_ROW davof-planeSurf-ringArea`, the row of the 2026-11-11 isoAdvector rows of the corrector study (davof, positioning ccs, wispTol 1e-3, applyTo all, interfaceCellRule normalOnly while the corrector is off and alphaOnly while it is on (forced by the scheme), all other keys the template's). AREA0 has `normalModel ringArea`. LIN0, LIN1, LIN2 have `normalModel ringLinear` with `ringLinearPasses 1; ringLinearCorrection true; ringLinearDiag true;` added by an edit of the rendered fvSolution because the template has no token for these keys (the only departure from the rendering; the keys exist in the port). LIN1 and LIN2 add `faceFractionCorrector mean; correctorIterations 1 | 2; correctorTolerance 0; correctorLog true;`. The full alpha dictionaries are at the end of this file.
- Metrics: E_shape, E_mass, E_bound at t = 3 from the last line of volumeFractionError.dat (written every 0.75, last line = end); L1 abs(n-n_e) = LNormalDiff1, L1(1-n.n_e) = LNormalDiff2^2/2, Linf abs(n-n_e) = LNormalDiffInf of reconstructionError at t = 3 (over the cells that hold a plane at the end; n_e the normal of the initial sphere, which is the exact solution at t = 3); min alpha and max alpha - 1 over the steps from the lines `After conservative bounding`; negative-weight and y = 0 counts from the `ringLinear pass 1` lines (reduced over the ranks; maximum over all reconstructions of the run); corrector changes from the `davof corrector (mean)` lines (per step, per iteration; the mean and the maximum over the steps). Interface cells are given twice: by alpha (1e-8 < alpha < 1 - 1e-8) in the written fields at t = 0, 1.5, 3, and as davof reports them in the log (first step, half the steps, last step; none for PLIC).
- Binaries (all transport jobs): libleiaDavof.so 84041d9ba17555a9, davofAdvectorFoam 547e31ad011dd062, libleiaDavofFunctionObjects.so 0017a107ff02820a, davofInitFields e7c5a00240b946bf, davofInitShearFlow3D 14a30a7b06a7fbb1; leia davof-iter-linfit; built in the workspace's own platforms/ directory.
- PP (parallel gate, job 57252508): davofTestCorrector, N 64 standard, 8 ranks against the serial run: ringLinear L1 abs(n-n_e) and L1(1-n.n_e) at passes 0, 1, 2 agree to a relative 1.7e-13 and 1.3e-11 at most (ringArea 4.6e-14 at most). PP held; N 64 runs on 8 ranks.

## How leia's isoAdvector scheme sets alpha_f in the davof rows (read from `src/davof/advection/isoAdvectorScheme.C`, `predictFaceFractions` and `afterAdvect`)

After isoAdvection's own step (its flux and bounding unchanged), every face that touches an interface cell of the step's reconstruction (the planes made at the START of the step, time level n, from alpha^n and alpha_f^n) takes the clip of the face polygon by the plane of its UPWIND cell translated by the cell velocity times dt (the plane at the END of the step, level n+1), by the other cell's translated plane when the upwind cell has none, and by the upwind cell's new alpha (clamped to [0, 1]) when neither cell has a plane; a face that touches no interface cell is not touched and keeps its alpha_f of the previous step, and a coupled face takes the half-sum of the two sides. The corrector acts directly afterwards in the same step (`afterAdvect` calls `predictFaceFractions()` and then `davof::correctFaceFractions()`): it reconstructs from the new alpha and this predicted alpha_f, replaces alpha_f by the mean of the two clips of the planes just reconstructed, and reconstructs again (so the next step starts from planes that match the corrected alpha_f), once or twice according to correctorIterations.

## Remarks on the numbers (not verdicts)

- The interface-cell counts by alpha of the rows without the corrector (AREA0, LIN0) grow from about 1 800 to 70 000 cells at N 64 by t = 3 (and to 14 000-17 000 at N 32), against about 1 800 for PLIC, LIN1 and LIN2: cells with alpha between 1e-8 and 1 - 1e-8 away from the interface (presumably wisps; not inspected) in the rows without the corrector, whose bounds are also the worst (min alpha down to -8e-3 and max alpha - 1 up to 7e-3; see the bounds columns). The counts davof reports in the log (its own interface-cell rule, wispTol 1e-3 in the rows without the corrector) stay close to PLIC's in every row (next column).
- The corrector's first iteration changes single faces by up to 1.0 (a face flips between 0 and 1); the mean change is 0.28-0.6 in the max norm and 0.005-0.014 in the L1.
- At N 32, LIN1 at Co 0.2 ends with L1 abs(n-n_e) 0.48 and Linf 1.99 although its E_shape is 1.09 times PLIC's: planes of near-empty cells at the end carry arbitrary normals; I did not investigate.
- Wall clock depends on the load of the node: eight jobs shared a node in the first matrix. PLIC against PLIC2 differs by up to 15 % (N 64, Co 0.2). A second set of jobs (the timing replicate) runs all five rows in one job per N 64 cell.

# Small transport test of the face-fraction corrector with the ring linear fit (isoAdvector, shear3D), tables

Generated by il_ttab.py from the logs and postProcessing of /work/scratch/tm83tomy/davof/agents/iterlin/transport/. Rows PLIC, AREA0, LIN0 ran in job A of a cell; PLIC2 (plicRDF again, the timing reference), LIN1, LIN2 in job B. Ratios are to the PLIC row of the SAME job.

## Run table (end of the period t = 3)

| N | Co | row | job | status | steps | E_shape | E_mass | E_bound | min alpha | max alpha - 1 | L1 abs(n-n_e) | L1 (1-n.n_e) | Linf abs(n-n_e) | interface cells t=0 / 1.5 / 3 (alpha in (1e-8, 1-1e-8)) | davof interface cells, first step / half / last step (log) | negW max | y=0 max | ExecutionTime s | ClockTime s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.2 | PLIC | A | completed | 720 | 0.004707 | 2.422e-12 | 2.311e-09 | -1.761e-06 | 4.441e-16 | 0.2925 | 0.05313 | 0.7363 | 416 / 936 / 468 | - |- | - | 20.48 | 21 |
| 32 | 0.2 | AREA0 | A | completed | 720 | 0.006291 | 1.452e-08 | 7.949e-08 | -0.0003753 | 7.305e-14 | 0.3694 | 0.08188 | 0.8009 | 416 / 8480 / 14305 | 416 / 948 / 503 |- | - | 13.05 | 14 |
| 32 | 0.2 | LIN0 | A | completed | 720 | 0.00496 | 3.25e-08 | 3.99e-06 | -0.00023 | 1.776e-14 | 0.2932 | 0.05286 | 0.6925 | 416 / 9440 / 15907 | 416 / 970 / 517 |173 | 5 | 14.75 | 16 |
| 32 | 0.2 | PLIC2 | B | completed | 720 | 0.004707 | 2.422e-12 | 2.311e-09 | -1.761e-06 | 4.441e-16 | 0.2925 | 0.05313 | 0.7363 | 416 / 936 / 468 | - |- | - | 21.84 | 22 |
| 32 | 0.2 | LIN1 | B | completed | 720 | 0.005119 | 1.619e-08 | 2.502e-08 | -0.0002901 | 4.441e-16 | 0.4793 | 0.2399 | 1.989 | 416 / 1011 / 565 | 416 / 1011 / 565 |198 | 63 | 21.77 | 23 |
| 32 | 0.2 | LIN2 | B | completed | 720 | 0.005086 | 1.562e-10 | 3.307e-14 | -3.483e-06 | 2.22e-16 | 0.3128 | 0.06578 | 1.943 | 416 / 1016 / 495 | 416 / 1016 / 495 |218 | 25 | 28.2 | 28 |
| 32 | 0.5 | PLIC | A | completed | 288 | 0.004603 | 5.02e-12 | 8.614e-07 | -0.0001409 | 4.441e-16 | 0.291 | 0.05157 | 0.7007 | 416 / 928 / 487 | - |- | - | 8.8 | 9 |
| 32 | 0.5 | AREA0 | A | completed | 288 | 0.00625 | 4.772e-08 | 7.27e-06 | -0.002151 | 1.51e-14 | 0.3657 | 0.07967 | 0.786 | 416 / 10843 / 17182 | 416 / 946 / 530 |- | - | 6.56 | 7 |
| 32 | 0.5 | LIN0 | A | completed | 288 | 0.004854 | 8.162e-08 | 5.03e-05 | -0.004721 | 1.71e-14 | 0.3497 | 0.1047 | 1.96 | 416 / 10396 / 16398 | 416 / 970 / 549 |190 | 12 | 8.56 | 8 |
| 32 | 0.5 | PLIC2 | B | completed | 288 | 0.004603 | 5.02e-12 | 8.614e-07 | -0.0001409 | 4.441e-16 | 0.291 | 0.05157 | 0.7007 | 416 / 928 / 487 | - |- | - | 9.82 | 10 |
| 32 | 0.5 | LIN1 | B | completed | 288 | 0.004856 | 2.915e-09 | 7.431e-06 | -0.0003156 | 2.353e-05 | 0.4126 | 0.1811 | 1.99 | 416 / 989 / 540 | 416 / 989 / 540 |188 | 15 | 10.15 | 11 |
| 32 | 0.5 | LIN2 | B | completed | 288 | 0.004859 | 1.432e-09 | 7.561e-06 | -0.0003114 | 3.306e-05 | 0.3219 | 0.08298 | 1.98 | 416 / 988 / 508 | 416 / 988 / 508 |192 | 21 | 12.08 | 12 |
| 64 | 0.2 | PLIC | A | completed | 1440 | 0.001342 | 4.094e-14 | 1.242e-09 | -3.022e-05 | 4.441e-16 | 0.1415 | 0.01318 | 0.361 | 1760 / 3954 / 1796 | - |- | - | 96.87 | 98 |
| 64 | 0.2 | AREA0 | A | completed | 1440 | 0.001881 | 4.133e-09 | 3.1e-06 | -0.0004714 | 0.0007962 | 0.1689 | 0.01798 | 0.4304 | 1760 / 40114 / 69491 | 1760 / 3924 / 1844 |- | - | 48.63 | 50 |
| 64 | 0.2 | LIN0 | A | completed | 1440 | 0.00164 | 6.247e-09 | 5.376e-06 | -0.000977 | 3.426e-10 | 0.15 | 0.01447 | 0.6163 | 1760 / 41212 / 73035 | 1760 / 3882 / 1789 |338 | 13 | 56.94 | 58 |
| 64 | 0.2 | PLIC2 | B | completed | 1440 | 0.001342 | 4.094e-14 | 1.242e-09 | -3.022e-05 | 4.441e-16 | 0.1415 | 0.01318 | 0.361 | 1760 / 3954 / 1796 | - |- | - | 81.89 | 83 |
| 64 | 0.2 | LIN1 | B | completed | 1440 | 0.001646 | 1.312e-13 | 9.785e-10 | -7.352e-05 | 4.039e-10 | 0.15 | 0.01444 | 0.383 | 1760 / 3909 / 1800 | 1760 / 3909 / 1800 |391 | 0 | 72.38 | 73 |
| 64 | 0.2 | LIN2 | B | completed | 1440 | 0.001645 | 7.233e-14 | 6.391e-10 | -1.048e-05 | 4.441e-16 | 0.1499 | 0.01443 | 0.3826 | 1760 / 3913 / 1802 | 1760 / 3913 / 1802 |394 | 0 | 94.88 | 96 |
| 64 | 0.5 | PLIC | A | completed | 576 | 0.001312 | 1.715e-11 | 6.682e-06 | -0.000343 | 0.0002226 | 0.1378 | 0.01252 | 0.3877 | 1760 / 3934 / 1877 | - |- | - | 44.88 | 45 |
| 64 | 0.5 | AREA0 | A | completed | 576 | 0.001856 | 1.755e-08 | 5.248e-05 | -0.008315 | 0.003891 | 0.1747 | 0.02287 | 1.813 | 1760 / 46657 / 79253 | 1760 / 3932 / 1977 |- | - | 23.73 | 24 |
| 64 | 0.5 | LIN0 | A | completed | 576 | 0.001608 | 2.076e-08 | 4.205e-05 | -0.006835 | 0.006903 | 0.159 | 0.02064 | 1.837 | 1760 / 46567 / 79716 | 1760 / 3885 / 1934 |378 | 9 | 27.2 | 28 |
| 64 | 0.5 | PLIC2 | B | completed | 576 | 0.001312 | 1.715e-11 | 6.682e-06 | -0.000343 | 0.0002226 | 0.1378 | 0.01252 | 0.3877 | 1760 / 3934 / 1877 | - |- | - | 44.79 | 45 |
| 64 | 0.5 | LIN1 | B | completed | 576 | 0.001607 | 9.872e-11 | 4.139e-07 | -0.0004048 | 4.547e-05 | 0.1472 | 0.01394 | 0.3771 | 1760 / 3910 / 1852 | 1760 / 3910 / 1852 |419 | 0 | 39.85 | 40 |
| 64 | 0.5 | LIN2 | B | completed | 576 | 0.001607 | 7.669e-11 | 4.277e-07 | -0.0004209 | 4.654e-05 | 0.1469 | 0.01391 | 0.3766 | 1760 / 3917 / 1849 | 1760 / 3917 / 1849 |419 | 0 | 49.51 | 50 |

## Ratios to PLIC of the same job

| N | Co | row | E_shape | E_shape / PLIC | ClockTime / PLIC | ExecutionTime / PLIC | L1(1-n.n_e) / PLIC |
|---|---|---|---|---|---|---|---|
| 32 | 0.2 | PLIC | 0.004707 | 1 | 1 | 1 | 1 |
| 32 | 0.2 | AREA0 | 0.006291 | 1.337 | 0.6667 | 0.6372 | 1.541 |
| 32 | 0.2 | LIN0 | 0.00496 | 1.054 | 0.7619 | 0.7202 | 0.995 |
| 32 | 0.2 | PLIC2 | 0.004707 | 1 | 1 | 1 | 1 |
| 32 | 0.2 | LIN1 | 0.005119 | 1.088 | 1.045 | 0.9968 | 4.515 |
| 32 | 0.2 | LIN2 | 0.005086 | 1.08 | 1.273 | 1.291 | 1.238 |
| 32 | 0.5 | PLIC | 0.004603 | 1 | 1 | 1 | 1 |
| 32 | 0.5 | AREA0 | 0.00625 | 1.358 | 0.7778 | 0.7455 | 1.545 |
| 32 | 0.5 | LIN0 | 0.004854 | 1.055 | 0.8889 | 0.9727 | 2.031 |
| 32 | 0.5 | PLIC2 | 0.004603 | 1 | 1 | 1 | 1 |
| 32 | 0.5 | LIN1 | 0.004856 | 1.055 | 1.1 | 1.034 | 3.512 |
| 32 | 0.5 | LIN2 | 0.004859 | 1.056 | 1.2 | 1.23 | 1.609 |
| 64 | 0.2 | PLIC | 0.001342 | 1 | 1 | 1 | 1 |
| 64 | 0.2 | AREA0 | 0.001881 | 1.401 | 0.5102 | 0.502 | 1.364 |
| 64 | 0.2 | LIN0 | 0.00164 | 1.222 | 0.5918 | 0.5878 | 1.098 |
| 64 | 0.2 | PLIC2 | 0.001342 | 1 | 1 | 1 | 1 |
| 64 | 0.2 | LIN1 | 0.001646 | 1.226 | 0.8795 | 0.8839 | 1.096 |
| 64 | 0.2 | LIN2 | 0.001645 | 1.226 | 1.157 | 1.159 | 1.095 |
| 64 | 0.5 | PLIC | 0.001312 | 1 | 1 | 1 | 1 |
| 64 | 0.5 | AREA0 | 0.001856 | 1.415 | 0.5333 | 0.5287 | 1.827 |
| 64 | 0.5 | LIN0 | 0.001608 | 1.226 | 0.6222 | 0.6061 | 1.649 |
| 64 | 0.5 | PLIC2 | 0.001312 | 1 | 1 | 1 | 1 |
| 64 | 0.5 | LIN1 | 0.001607 | 1.225 | 0.8889 | 0.8897 | 1.114 |
| 64 | 0.5 | LIN2 | 0.001607 | 1.225 | 1.111 | 1.105 | 1.112 |

## PLIC against PLIC2 (the same case in the two jobs of a cell: wall-clock scatter)

| N | Co | PLIC ClockTime / Execution | PLIC2 ClockTime / Execution | relative difference of ExecutionTime | E_shape PLIC / PLIC2 |
|---|---|---|---|---|---|
| 32 | 0.2 | 21 / 20.48 | 22 / 21.84 | 0.06641 | 0.004707 / 0.004707 |
| 32 | 0.5 | 9 / 8.8 | 10 / 9.82 | 0.1159 | 0.004603 / 0.004603 |
| 64 | 0.2 | 98 / 96.87 | 83 / 81.89 | 0.1546 | 0.001342 / 0.001342 |
| 64 | 0.5 | 45 / 44.88 | 45 / 44.79 | 0.002005 | 0.001312 / 0.001312 |

## E_shape orders N 32 - 64 (log2 of the ratio)

| Co | row | E_shape N 32 | E_shape N 64 | order |
|---|---|---|---|---|
| 0.2 | PLIC | 0.004707 | 0.001342 | 1.81 |
| 0.2 | AREA0 | 0.006291 | 0.001881 | 1.742 |
| 0.2 | LIN0 | 0.00496 | 0.00164 | 1.596 |
| 0.2 | LIN1 | 0.005119 | 0.001646 | 1.637 |
| 0.2 | LIN2 | 0.005086 | 0.001645 | 1.628 |
| 0.5 | PLIC | 0.004603 | 0.001312 | 1.811 |
| 0.5 | AREA0 | 0.00625 | 0.001856 | 1.752 |
| 0.5 | LIN0 | 0.004854 | 0.001608 | 1.594 |
| 0.5 | LIN1 | 0.004856 | 0.001607 | 1.595 |
| 0.5 | LIN2 | 0.004859 | 0.001607 | 1.596 |

## Corrector changes of alpha_f (per step: mean and maximum over the steps) and the identity defect

| N | Co | row | steps with the corrector | iteration | max abs(d alpha_f) mean | max abs(d alpha_f) max | L1 d alpha_f mean | L1 d alpha_f max | identity defect max |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.2 | LIN1 | 720 | 1 | 0.5835 | 1 | 0.01333 | 0.02051 | 5.581e-16 |
| 32 | 0.2 | LIN2 | 720 | 1 | 0.6364 | 1 | 0.01269 | 0.01874 | 4.908e-16 |
| 32 | 0.2 | LIN2 | 720 | 2 | 0.1908 | 0.9997 | 0.001599 | 0.003705 | 4.908e-16 |
| 32 | 0.5 | LIN1 | 288 | 1 | 0.5942 | 1 | 0.01371 | 0.01763 | 4.212e-16 |
| 32 | 0.5 | LIN2 | 288 | 1 | 0.6023 | 1 | 0.01387 | 0.01796 | 3.608e-16 |
| 32 | 0.5 | LIN2 | 288 | 2 | 0.257 | 1 | 0.001703 | 0.003247 | 3.608e-16 |
| 64 | 0.2 | LIN1 | 1440 | 1 | 0.2801 | 1 | 0.005262 | 0.00784 | 9.647e-16 |
| 64 | 0.2 | LIN2 | 1440 | 1 | 0.2798 | 1 | 0.005247 | 0.007806 | 8.741e-16 |
| 64 | 0.2 | LIN2 | 1440 | 2 | 0.07535 | 0.4194 | 0.0004856 | 0.001074 | 8.741e-16 |
| 64 | 0.5 | LIN1 | 576 | 1 | 0.3104 | 1 | 0.005899 | 0.008182 | 1.251e-15 |
| 64 | 0.5 | LIN2 | 576 | 1 | 0.3102 | 1 | 0.005887 | 0.007869 | 1.129e-15 |
| 64 | 0.5 | LIN2 | 576 | 2 | 0.08566 | 0.4735 | 0.0005257 | 0.001058 | 1.129e-15 |

## Verdicts (from the numbers only; N 64, both Co)

- PS1 (E_shape(LIN0)/E_shape(PLIC) > 1.05): Co 0.2 1.222, Co 0.5 1.226: **held**
- PS2 (E_shape(LIN1) < E_shape(LIN0)): Co 0.2: 0.001646 vs 0.00164; Co 0.5: 0.001607 vs 0.001608: **falsified**
- PS3 (E_shape(LIN1)/E_shape(PLIC) > 1.05): Co 0.2 1.226, Co 0.5 1.225: **held**
- PS4 (E_shape(LIN2) within 2 % of E_shape(LIN1)): Co 0.2 0.0439 %, Co 0.5 0.0132 %: **held**
- PS5 (wall(LIN0)/wall(PLIC) < 1 and wall(LIN1)/wall(PLIC) > 1, ClockTime): LIN0 0.5918 (Co 0.2), 0.6222 (Co 0.5); LIN1 0.8795, 0.8889: **falsified**
- Decision rule: accuracy clause (E_shape(LIN1)/E_shape(PLIC) <= 1.05 at both Co) fails; cost clause (wall(LIN1)/wall(PLIC) <= 0.8 at both Co) fails: **the accuracy clause fails: closed for this round**

## Timing replicate (N 64, all five rows in one job each, PLIC first; ratios inside the job)

| Co | row | status | E_shape | ClockTime s | ExecutionTime s | ClockTime / PLIC | ExecutionTime / PLIC |
|---|---|---|---|---|---|---|---|
| 0.2 | PLIC | completed | 0.001342 | 63 | 62.04 | 1 | 1 |
| 0.2 | AREA0 | completed | 0.001881 | 29 | 27.97 | 0.4603 | 0.4508 |
| 0.2 | LIN0 | completed | 0.00164 | 34 | 34.41 | 0.5397 | 0.5546 |
| 0.2 | LIN1 | completed | 0.001646 | 53 | 52.88 | 0.8413 | 0.8524 |
| 0.2 | LIN2 | completed | 0.001645 | 70 | 69.07 | 1.111 | 1.113 |
| 0.5 | PLIC | completed | 0.001312 | 38 | 37.41 | 1 | 1 |
| 0.5 | AREA0 | completed | 0.001856 | 13 | 12.63 | 0.3421 | 0.3376 |
| 0.5 | LIN0 | completed | 0.001608 | 16 | 15.86 | 0.4211 | 0.424 |
| 0.5 | LIN1 | completed | 0.001607 | 23 | 22.91 | 0.6053 | 0.6124 |
| 0.5 | LIN2 | completed | 0.001607 | 31 | 29.8 | 0.8158 | 0.7966 |

## Key sets of the alpha solver dictionary of every row (N 64, Co 0.5; N and Co change no key)


### PLIC

```
"alpha.water.*"
    {
```

### AREA0

```
"alpha.water.*"
    {
```

### LIN0

```
"alpha.water.*"
    {
```

### LIN1

```
"alpha.water.*"
    {
```

### LIN2

```
"alpha.water.*"
    {
```

## Jobs (all serial or 8 ranks on one node, deflt_short; none cancelled except the eight named below)

| job | purpose | final state |
|---|---|---|
| 57252508 | PP, N 64 standard, 8 ranks, ringLinear and ringArea, passes 0-2 | COMPLETED |
| 57252878 | probe, PLIC N 64 Co 0.5, 8 ranks (ClockTime 39 s) | COMPLETED |
| 57253132-57253139 | transport jobs submitted before the pre-registration commit, cancelled, never read | CANCELLED (by id) |
| 57253296 / 57253297 | N 64, Co 0.5, job A (PLIC, AREA0, LIN0) / job B (PLIC2, LIN1, LIN2) | COMPLETED |
| 57253298 / 57253299 | N 64, Co 0.2, A / B | COMPLETED |
| 57253300 / 57253301 | N 32, Co 0.5, A / B (1 rank) | COMPLETED |
| 57253302 / 57253303 | N 32, Co 0.2, A / B (1 rank) | COMPLETED |
| 57254046 / 57254047 | timing replicate, N 64, Co 0.5 / 0.2, all five rows in one job | COMPLETED |
