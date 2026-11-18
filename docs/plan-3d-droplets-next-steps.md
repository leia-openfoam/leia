# Plan: 3D droplets, next steps, and the data policy on scratch

Date: 2026-11-18. Branch `work/droplets3D`. Sources: STATUS 11.42 (the first 3D baseline),
STATUS 11.41 on the branch `fix/sl-remote-patch-faces`, STATUS 11.38.
Status of every decision below: **OPEN, the supervisor decides.** Claude's recommendation is
written next to each one. Nothing in this plan is started before the data clean-up of
section 3 (supervisor, 2026-11-18).

## 1. Where we are (MEASURED, STATUS 11.42)

| droplet (3D, N = 60, 78, 102, 32 ranks) | result |
|---|---|
| stationary | four of five proposed criteria hold; the spurious current (1.9e-6 m/s) does not fall with N |
| oscillating | period -1.4, -1.8, -2.7 % against Lamb (8.232 ms, corrected); volume error 2.4e-3 |
| translating | ends inside a blow-up; two separate defects (the curvature growth, the step-2 disturbance) |

## 2. Decisions

### D1. Merge `fix/sl-remote-patch-faces` into `development`

- Content: patch faces of another rank in the semi-Lagrangian stencil (base class), and the
  metric `MEAN_DELTA_X`.
- Evidence: P1 and P3 hold; P2 and P4 fail as written, at sizes that come from round-off
  (growth 1.25 per step with the non-default option `exclude`; summation order 4.7e-12).
- Recommendation: merge. The default paths are identical in all standard studies.
- If merged: every multi-rank result with `valueBound lipschitzCone` from an older tree
  carries the defect (1.6 to 1.9 % in the volume error of the polyhedral arms).

### D2. Merge `work/droplets3D` into `development`

- Content: the corrected 3D Lamb period (script, gate config, case comment, knowledge base),
  the tables of the first 3D baseline, three discriminator configs, this plan. No solver
  code.
- Recommendation: merge.

### D3. The accuracy criteria for the 3D droplets

- Proposed in the header of `config/gates/methodGate3D.yaml` (A-stationary, A-translating,
  A-oscillating). The reference period there is corrected to 8.232 ms.
- Open point: a reference for the damping rate. Lamb's rate for a viscous drop without an
  outer fluid is (n - 1)(2n + 1) nu / R^2 = 5.0 1/s (DERIVED; the air is neglected).
- Recommendation: accept the criteria as written, and add the damping reference as a
  reported value, not a criterion.

### D4. The translating droplet, defect 2: the disturbance at step 2

- MEASURED: 1.6e-2 m/s on N = 60 and 9.0e-2 m/s on N = 102 at the first three-level step;
  3.8e-4 m/s with Euler in both equations; a relative mass residual of 4.1e-3 at step 2 only.
- HYPOTHESIS: at step 2 the density of the level n-1 is the initial density, which is not
  consistent with the end-of-step reset that every later level has.
- Proposed step (cheapest discriminator first):
  1. Read-only: print the difference between the initial density and the reset density of
     step 0 on the N = 60 case (one short run).
  2. If the difference is not zero: a start-up option of the mass-flux family, selectable at
     run time and inert by default, that applies the same reset to the initial density.
  3. Gate: bit-identity with the option off; with the option on, the step-2 disturbance on
     N = 60 and N = 102 against the limit 1e-3 m/s; then the translating arm of the 3D gate.
- Recommendation: do this first. It is small, and it removes a disturbance of up to 1.8 U0.

### D5. The translating droplet, defect 1: the growth of the curvature error

- MEASURED: the band curvature deviation doubles about every 2 ms from t = 0.008 s; the same
  rate on N = 60 and 102, with BDF2 and Euler, and in a 6 mm and an 8 mm box until 0.016 s.
  The band gradient error grows linearly from the first step.
- Falsified: the BDF2 disturbance as the cause; a fixed outlet distance as the onset of the
  final blow-up (negative log, 2026-11-18).
- HYPOTHESIS: the loss of the signed-distance property under translation drives it.
- Options:
  1. A discriminator for the hypothesis: a redistancing INSTRUMENT on the 3D translating
     case. The 3D templates write `type noRedistancing`; the token must be added first
     (inert default, its own commit).
  2. The gradient-control candidates of `config/candidates/` on the 3D translating arm. All
     of them failed the 2D gate.
  3. Wait for the PDE level-set line with the velocity extension (another session, steps A
     and B in progress).
- Recommendation: option 1 as a measurement, then decide between 2 and 3 on its result.
- Also open: the cause of the final blow-up. A study that moves the outlet alone needs a
  droplet-offset token in the 3D template.

### D6. The oscillating droplet

- MEASURED: the period error grows under refinement. No discriminator ran.
- Proposed: one run of the N = 78 rung with the exact initial shape and a prescribed-free
  transport check is not available; the cheapest available discriminator is the horizon:
  the same ladder to 0.05 s (6 periods), to separate a start-up effect from a drift.
- Recommendation: after D4, because the oscillating droplet has no mean velocity and D4
  does not act on it; run it only if the supervisor wants the oscillating case before the
  translating one.

### D7. Not decided, listed for completeness

- The viscosity model of the 3D templates (`geo_lin`; the 2D best is `alg_lin`; no token).
- The deletion of the branches `integration/sl-pde`, `feature/sl-ghost-cells`,
  `fix/ghost-decomposition` (all contained in `development`).

## 3. Data on `/work/scratch/tm83tomy` (before any further run)

The account is at 24.9 of 25.0 million files (hard limit). The soft limit is 20 million.

### 3.1 Inventory of the directories of this session (MEASURED, job 56621795)

| directory | files | size | content |
|---|---|---|---|
| `leia-drop3d` | 138,917 | 54 GB | the 3D gate, the smoke, three discriminators, void runs |
| `leia-integration` | 226,157 | 289 GB | identity checks, the patch-face runs, the ghost gate re-run, 157 GB of tar files |
| `leia-integration-base` | 129,305 | 57 GB | the base tree of the identity checks, rounds 1 and 2 |
| `leia-ghost` | 28,930 | 94 GB | the ghost-cell gate and the HPCToolkit profile runs |
| `leia-amp-schemes`, `leia-amp-schemes-c`, `.leia-int-logs` | 144 | 1 MB | the scheme-matrix CSVs and logs |
| sum | 523,453 | 494 GB | |
| `leia-curvature` (older work of this project line, last change before 2026-11-10) | 583,452 | 281 GB | 569,886 files under `studies/` |

About 95 % of the files are in `processor*` directories. This session holds 2.1 % of the
files of the account (4.4 % with `leia-curvature`). The clean-up below cannot bring the
account below the soft limit; the other 23.8 million files belong to other work.

### 3.2 What can be deleted (proposal; nothing is deleted before the supervisor's word)

Every result of these directories is recorded in STATUS (11.35 to 11.42) and its curated
table is committed. A case can be simulated again from its committed config.

| item | files | size | reason |
|---|---|---|---|
| the eight tar files in `leia-integration/archives_20261118/` | 8 | 157 GB | raw data of finished identity rounds and locate runs |
| `leia-integration-base` (the whole clone) | 129,305 | 57 GB | reference tree of finished checks |
| `leia-integration` (the whole clone), after a copy of its CSVs and logs | 226,000 | 132 GB | identity re-runs and locate runs; both fix branches are on the hub |
| `leia-ghost` (the whole clone), after a copy of its CSVs, logs and the profile database | 28,930 | 94 GB | the branch is in `development`; the gate tables are committed |
| in `leia-drop3d`: the `_VOID_quota` and the smoke directories | 6,400 | 1.2 GB | void or not a measurement |
| in `leia-drop3d`: the `processor*` directories and all but the last written time of the finished studies | 125,000 | 50 GB | the CSVs, the logs and the last reconstructed time stay |
| `leia-amp-schemes`, `leia-amp-schemes-c`, `.leia-int-logs` | 144 | 1 MB | the CSVs are committed |

Result: about 517,000 files and 490 GB free; about 6,000 files stay (CSVs, logs, case
dictionaries, one reconstructed time per case of the 3D gate).
`leia-curvature`: the same rule (keep CSVs and logs, delete fields) frees about 560,000
files. The supervisor confirms first that no other session reads it.

### 3.3 How to reduce the output of future runs (proposal)

One 3D gate case on 32 ranks writes about 6,900 files now: 11 times, 33 directories
(32 processors and the reconstructed case), about 19 objects per time, ASCII.

| step | change | effect on one such case | gate |
|---|---|---|---|
| 1 | `purgeWrite 2` in the case templates, by a token with the present value as default | 2 written times in place of 11; the state before a failure is always there | none for the numbers; the report figures that need an intermediate time name their study |
| 2 | a workflow key `keep_fields: last` (default `all`): `reconstructPar -latestTime`, then delete `processor*`; for a failed case keep the `processor*` times and do not reconstruct | 1 directory in place of 33 | CSV identity on one study |
| 3 | `fileHandler collated` for parallel runs | one `processors32` directory while the run is live: 32 times fewer files during the run | 4-rank smoke and CSV identity |
| 4 | `writeFormat binary` | about half the size, the same file count | CSV identity |
| 5 | the solver writes only psi, U, p_rgh, alpha and kappa unless a `debugFields` switch is on (now 19 objects, with restart and diagnostic fields) | 5 objects in place of 19 | bit-identity of the CSVs; own commit |
| 6 | diagnostic jobs that write every step (`ghost_locate.sbatch`) pack their output into one tar at the end | 12,500 files become 1 | none |
| 7 | the gate removes the smoke directories after a pass | 2,900 files per gate | none |

With steps 1, 2 and 5 a finished case keeps about 10 field files and its CSVs and logs;
with step 3 a running case holds about 40 files in place of 6,900.
Recommendation: steps 1, 2, 6 and 7 first (workflow and templates only), then 3 and 4 with
their smoke, then 5 as a solver change with its gate.
