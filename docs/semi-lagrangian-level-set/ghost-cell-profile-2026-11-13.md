# HPCToolkit profile of the ghost-cell reconstruction (2026-11-13)

Solver `leiaSemiLagrangeLevelSetFoam`, case `3Dtranslation`, mesh `hexRefined`, N = 128
(1 269 760 cells, 296 steps, CFL 0.5, END_TIME 1), the N = 128 rung of
`config/gates/ghostCellRefinedBox3D.yaml`. Classes: `ghostCellQuadraticWeightedLeastSquares`
(ghost class) and `uncachedQuadraticWeightedLeastSquares` (baseline class).
Every number below is MEASURED and comes from a file, unless the line says DERIVED or
HYPOTHESIS. The build is `Opt` without `-g`: the attribution is to functions, not to lines.

## 1. Tool, commands, jobs

- HPCToolkit 2025.0.1 (spack, `~mpi`, `+papi`), `$HOME/hpctoolkit` ->
  `/work/home/tm83tomy/spack/opt/spack/linux-x86_64/hpctoolkit-2025.0.1-ky3i4lswv3ynbyvk6ygnunxtn6tsfa67`.
  The tools are called by path (`$HOME/hpctoolkit/bin/{hpcrun,hpcstruct,hpcprof}`); no module.
  `$HOME/Software/hpctoolkit-install` holds version 2017.06 (libraries missing); not used.
- Clone `/work/scratch/tm83tomy/leia-ghost`, branch `feature/sl-ghost-cells`. Commit
  a982f147 for the first four runs, 32f1cfe7 for the others (no source change between them).
- Preparation: `sbatch workflow/hpctoolkit/prep.sbatch` (job 56296196). The workflow runs
  `--until preprocess` for `config/gates/ghostCellProfileNp1.yaml` and `--until decompose`
  for `ghostCellProfileNp4.yaml`. The prepared cases are never solved.
- One run: `workflow/hpctoolkit/run.sbatch` copies a prepared case to
  `studies/ghostCellProfileRuns/<run>` and runs the solver command of the workflow there:
  - serial: `hpcrun -e REALTIME@5000 -o <run>/hpctoolkit-measurements leiaSemiLagrangeLevelSetFoam`
  - 4 ranks: `srun --ntasks=4 --overlap --cpu-bind=none --kill-on-bad-exit=1 hpcrun -e REALTIME@5000 -o ... leiaSemiLagrangeLevelSetFoam -parallel`
  - then `hpcstruct -j <cpus> hpctoolkit-measurements` and
    `hpcprof -j <cpus> -o hpctoolkit-database hpctoolkit-measurements`.
- Text extraction: `workflow/hpctoolkit/hpcdb_report.py <database>` reads `meta.db` and
  `profile.db` (format 4.0, `share/doc/hpctoolkit/FORMATS.md`). It ran on the laptop on
  copies of the two files (3 MB to 7 MB). Its outputs are in
  `docs/semi-lagrangian-level-set/ghost-cell-profile-2026-11-13/`.
- All jobs: account `special00005`, `--mem-per-cpu=12000`, partition `mem`.

| run | job | node (feature) | wall clock of the solver [s] | max. resident memory |
|---|---|---|---|---|
| `ghost_np1_hpc` | 56296324 | mpqc0001 (m01) | 1644 | 7.02 GB |
| `ghost_np1_plain` | 56296325 | mpqc0002 (m01) | 1774 | 7.00 GB |
| `base_np1_hpc` | 56296326 | mpqd0001 (m02) | 917 | 2.92 GB |
| `base_np1_plain` | 56296327 | mpqd0001 (m02) | 969 | 2.86 GB |
| `base_np1_hpc_m01` | 56296358 | mpqc0002 (m01) | 1642 | 2.90 GB |
| `base_np1_plain_m01` | 56296359 | mpqc0002 (m01) | 1665 | 2.89 GB |
| `ghost_np4_hpc` | 56297453 | mpqc0001 (m01) | 602 | 1.85 GB per rank |
| `ghost_np4_plain` | 56297454 | mpqc0001 (m01) | 583 | 1.80 GB per rank |

Wall clock: `ClockTime` of the last step in the solver log. Memory: `/usr/bin/time -v` for
the serial runs, `sacct` MaxRSS of the solver step for the 4-rank runs. The m01 nodes have
the Xeon Platinum 9242. The m02 node is 1.7 to 1.8 times faster for this solver (917 s
against 1642 s, same class, byte-identical metrics CSV). Compare runs only on one node type.

Measurements and databases (absolute):
`/work/scratch/tm83tomy/leia-ghost/studies/ghostCellProfileRuns/<run>/hpctoolkit-measurements`
and `.../<run>/hpctoolkit-database`, for `<run>` = `ghost_np1_hpc`, `base_np1_hpc`,
`base_np1_hpc_m01`, `ghost_np4_hpc`.

Profiler overhead: the difference between the run with hpcrun and the plain run is
-7 % (ghost serial, two nodes), -5 % (baseline, m02), -1 % (baseline, m01) and +3 %
(ghost, 4 ranks, same node, same time). The overhead is smaller than the run-to-run
difference on these shared nodes.

## 2. Serial, ghost class (`ghost_np1_hpc`, 1643.5 s sampled)

Inclusive time, call path `main -> slAdvection::advect -> pointValueScheme::advance ->
directCorrector::correct -> ghostCell...::update`:

| function | inclusive [s] | share | exclusive [s] |
|---|---|---|---|
| `slAdvection::advect` | 1217.0 | 74.1 % | |
| `directCorrector::correct` | 1178.2 | 71.7 % | 1.6 |
| `ghostCell...::update` | 1132.3 | 68.9 % | 2.8 |
| -- `fitCell` | 985.3 | 60.0 % | 184.3 |
| ---- `solveFit` (contains `choleskySolve`, inlined) | 566.7 | 34.5 % | 565.4 |
| ---- `basis` | 103.3 | 6.3 % | 103.3 |
| ---- heap: `List<double>::List` 40.6, `free` 28.1, `List<vector>::resize_nocopy` 25.6 | 94.3 | 5.7 % | |
| -- `collectData<tensor>` (coefficients of all stencil entries) | 109.4 | 6.7 % | 94.8 |
| -- `collectStencil` (`collectData<scalar>`, psi) | 16.3 | 1.0 % | |
| -- `collectData<scalar>` (mean of the children) | 15.7 | 1.0 % | |
| -- `build` (once) | 2.3 | 0.1 % | |
| `slCorrector::robustEvaluate` (departure-point value) | 41.0 | 2.5 % | 3.4 |
| -- `stencilRange` | 24.0 | 1.5 % | 24.0 |
| -- `evaluateRaw` | 10.4 | 0.6 % | |
| `fvc::grad<vector>` in `directCorrector::correct` | 32.4 | 2.0 % | |
| `slReconstruction::computeLimiters` | 0.35 | 0.02 % | |
| `slCorrector::footRadiusGuard` | 2.1 | 0.1 % | |
| function objects (`functionObjectList::execute`) | 220.3 | 13.4 % | |
| -- `writeIsoSurfaceTopo` (every step) | 103.6 | 6.3 % | |
| -- `gradPsiErrorCSV::write` | 51.5 | 3.1 % | |
| -- `fieldExpression::execute` (grad, gradPsiError) | 39.0 | 2.4 % | |
| -- field output of the function objects | 26.1 | 1.6 % | |
| `fvc::reconstruct<scalar>` in `main` | 83.0 | 5.1 % | |
| `reportErrors` in `main` | 37.7 | 2.3 % | |
| `Time::writeObject` (field output) | 30.6 | 1.9 % | |
| `slAdvection::New` (stencil construction, once) | 13.8 | 0.8 % | |

Top functions by exclusive time: `solveFit` 565.4 s (34.4 %), `fitCell` 184.3 (11.2 %),
`basis` 107.7 (6.6 %), `collectData<tensor>` 94.8 (5.8 %), `powf32x` 42.9 (2.6 %),
`_int_free` 31.6, `collectData<scalar>` 29.5, `_int_malloc` 26.0, `stencilRange` 24.0,
`isoSurfaceBase::getCellCutType` 22.7, `leastSquaresGrad::calcGrad` 22.7,
`outer<vector,vector>` 22.3, `malloc` 21.5, `volPointInterpolation::makeInternalWeights`
17.0, `fvc::surfaceSum<tensor>` 16.5. By library: `libleiaSemiLagrangian` 55.9 %,
`libfiniteVolume` 15.1 %, `libc` 11.6 %, solver executable 8.8 %.

## 3. Serial, baseline class, and the attribution of the difference

`base_np1_hpc_m01` (same node type as the ghost run, 1648.1 s sampled):
`update` 882.6 s (53.6 %): exclusive 596.1, `choleskySolve` 154.9, `basis` 102.1,
`collectStencil` 18.9. `robustEvaluate` 41.2. Function objects 422.8 s (25.7 %), of which
`writeIsoSurfaceTopo` 295.4 s. `fvc::reconstruct` 91.1 s.
`base_np1_hpc` (m02, 920.3 s) has the same shares: `update` 466.7 s (50.7 %), function
objects 250.4 s (27.2 %).

| part (m01, serial) | baseline [s] | ghost [s] | difference [s] |
|---|---|---|---|
| reconstruction `update` | 882.6 | 1132.3 | +249.7 |
| -- fit arithmetic (assembly, Cholesky, basis) | 853.1 | 857.2 | +4.1 |
| -- `collectData<tensor>` | 0 | 109.4 | +109.4 |
| -- heap allocation in `fitCell` | not visible | 94.3 | +94.3 |
| -- `collectData<scalar>` of the parent mean | 0 | 15.7 | +15.7 |
| -- psi stencil, `build`, rest | 29.5 | 55.7 | +26.2 |
| function objects | 422.8 | 220.3 | -202.5 |
| -- `writeIsoSurfaceTopo` | 295.4 | 103.6 | -191.8 |
| total | 1648.1 | 1643.5 | -4.6 |

1. The `update` of the ghost class costs 1.28 times the `update` of the baseline class.
   The fit arithmetic costs the same. The difference is the exchange of the coefficient
   tensor for every stencil entry, the heap allocation in `fitCell`, and the parent mean.
2. The serial totals are equal by coincidence. The baseline run has a large error on this
   mesh (`L_INF_E_PSI` = 2.36e6 at T = 1 in its CSV, gate result G1), and its iso-surface
   function object costs 3 times as much (`isoSurfaceTopo::removeInsidePoints` 147.0 s
   against 9.6 s). The wall-clock ratio of the two classes on this case measures the
   function objects too. Compare `update`, not the total.
3. Memory: 7.0 GB against 2.9 GB, ratio 2.4. DERIVED: the ghost class stores, per stencil
   entry, `stLevel_` (8 B), `stParentC_` (24 B), `stParentMean_` (8 B) and `stCoeff_`
   (72 B). With 27 entries per cell and 1.27e6 cells this is 3.8 GB, of which `stCoeff_` is
   2.5 GB. The measured difference is 4.1 GB. The profile has no memory event.

## 4. Four ranks, ghost class (`ghost_np4_hpc`)

- Wall clock 602 s with hpcrun, 583 s plain. Speed-up against the serial run: 2.73 (with
  hpcrun, 1644 s) and 3.04 (plain, 1774 s). Efficiency 0.68 and 0.76. The serial plain run
  was on another node, so the speed-up is between 2.7 and 3.0.
- Per rank (main thread): total 601.1, 601.2, 601.3, 601.3 s. Sum 2404.9 s against 1643.5 s
  serial: 761 s more.
- Time below MPI: 685.8 s = 28.5 % (per rank 171.9, 166.8, 176.1, 171.0 s). Two call sites
  hold 622.5 s of it:

| MPI call and call site | rank 0 | rank 1 | rank 2 | rank 3 | sum [s] |
|---|---|---|---|---|---|
| `MPI_Waitall` <- `UPstream::waitRequests` <- `GeometricBoundaryField<tensor>::evaluate` <- `ghostCell...::update` (`coeffFld_.correctBoundaryConditions()` between level 0 and level 1) | 2.9 | 148.6 | 146.8 | 0.5 | 298.7 |
| `MPI_Allreduce` <- `slCorrector::footRadiusGuard` <- `directCorrector::correct` (the first collective after `update`) | 155.8 | 0.6 | 9.5 | 157.9 | 323.8 |
| `MPI_Waitall` <- processor patches of `leastSquaresGrad` (function object gradPsiError) | 4.7 | 7.7 | 9.0 | 4.2 | 25.6 |
| all other MPI calls | | | | | 37.7 |

- This time is waiting, not communication. The decomposition (scotch) balances the cells
  (315 774 to 319 244 per rank) but not the levels: ranks 0 and 3 hold 196 116 and 196 523
  level-0 cells and 123 128 and 120 917 level-1 cells; ranks 1 and 2 hold 11 261 and 13 892
  level-0 cells and 306 041 and 301 882 level-1 cells (job 56297985). `update` fits level 0,
  then exchanges the coefficients, then fits level 1. Ranks 1 and 2 finish level 0 early and
  wait in the exchange. Ranks 0 and 3 finish level 1 early and wait in the next collective.
  `fitCell` itself is balanced: 260.4, 270.7, 263.2, 261.1 s.
- `update` inclusive per rank: 301.4, 456.8, 447.9, 299.4 s.
- DERIVED: without the 622.5 s the sum is 1782 s, 446 s per rank, speed-up 3.7.
- Work that grows with the rank count (sum over ranks against serial): function objects
  272.4 s against 220.3 s (`writeIsoSurfaceTopo` 127.1 against 103.6), heap allocation in
  `fitCell` 156.8 s against 94.3 s. `solveFit` is equal (563.6 against 566.7).
  No part is serial per step: no function has the serial time on every rank (examples per
  rank against serial: function-object output 20.6 s against 77.6 s, `Time::writeObject`
  7.0 s against 30.6 s, `reportErrors` 9.6 s against 37.7 s).
- Memory: 1.80 GB to 1.85 GB per rank, 7.2 GB to 7.4 GB in sum, against 7.0 GB serial.
- Three threads per rank exist (two OpenMPI helper threads). The report uses thread 0.

## 5. Serial against 4 ranks: the metrics CSV

`workflow/scripts/compare_metrics_csv.py ghost_np1.csv ghost_np4.csv --tol 0 --skip
ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`: FAIL, also with `--tol 1e-10`. 297 rows, equal step
counts. Largest relative difference per column:

| column | largest relative difference | where | first difference |
|---|---|---|---|
| `E_VOL_ALPHA` and its three variants | 6.06e-05 | row 213, t = 0.678 | row 8, 1.4e-12 |
| `E_GEOM_ALPHA`, `E_GEOM_ALPHA_REL` | 9.98e-07 | row 295, t = 0.9985 | row 195, 3.6e-12 |
| `L_INF_E_PSI` | 3.17e-07 | row 295 | row 237, 8.1e-12 |
| `TIME`, `DELTA_X`, `E_BOUND_ALPHA`, `MAX_CFL`, `MAX_ABS_SUM_VOL_FLUXES` | 0 | | |

The volume error at row 213 is 2.13e-07, so 6e-05 relative is 1.3e-11 absolute. The runs
with and without hpcrun are byte-identical in these columns, in serial and on 4 ranks, so
each decomposition is deterministic. The answer depends on the decomposition: the
difference starts at round-off size and grows to 3e-07 in the largest psi error in 60 steps.
`gradPsiError.csv` differs too: `MEAN_DELTA_X` by 8.5e-04 at t = 0 (the metric itself
depends on the decomposition), `E_L2_GRAD_PSI` by 1.0e-04 and `E_NARROW_L1_GRAD_PSI` by
9.2e-06 at T = 1.
HYPOTHESIS: the seed is the summation order of `syncTools::syncPointList` with `plusEqOp`
in the restriction (the mean of the children across a processor boundary), and the scheme
amplifies it. Not tested. The discriminator is a cell-by-cell difference of psi between the
serial case and the reconstructed 4-rank case at the first write time, located against the
processor boundaries and the two transitions.

## 6. Not measured, and what to change first

Not measured: line-level attribution (no `-g`; a rebuild was not necessary for these
conclusions). The time of `choleskySolve` inside `solveFit` and of the restriction loop
inside `update` (inlined). The time per refinement level inside `fitCell` (derived from the
cell counts only). Memory by allocation site. The 8-rank run of the gate. The baseline class
on 4 ranks, so the 4-rank ratio of the two classes is not known. The CPU model of the m02
node.

Changes, in this order:

1. Remove the level-by-level synchronization (MEASURED: 622.5 of 2404.9 rank-seconds,
   26 %). Only a level-1 cell with a coarser cell in its stencil needs the level-0
   coefficients. Fit every cell that has no coarser stencil entry first, on all ranks, then
   exchange the coefficients once, then fit the cells at the transition. HYPOTHESIS: the
   wait falls to the size of the transition set; the bound is speed-up 3.7 on 4 ranks.
2. Exchange and store the coefficient tensor only for the stencil entries that are coarser
   than their cell (MEASURED: `collectData<tensor>` 109.4 s serial, 6.7 %; DERIVED: 2.5 GB).
3. Remove the heap allocation from `fitCell`: four `DynamicList` and three `scalarField`
   per cell and step (MEASURED: 94.3 s serial, 5.7 %; 156.8 rank-seconds on 4 ranks).
4. The fit arithmetic is 853 s to 857 s serial, 52 %, in both classes. On a static mesh the data
   positions, the weights and the normal matrix of a cell do not change. HYPOTHESIS: a
   stored factor per cell removes the assembly and the factorization from every step.
5. Not the class, but in every run: `writeIsoSurfaceTopo` runs every step and builds
   `volPointInterpolation` every time (49.5 s of its 103.6 s); `pow` in `gradPsiErrorCSV`
   and `reportErrors` costs 45 s. Timing comparisons of reconstruction classes must not
   include the function objects.
6. Resolve the decomposition dependence of section 5 before any 4-rank or 8-rank result of
   the ghost class is used as a reference.

`run-studies.sbatch`, `run-sink-decay.sbatch` and the `cases/*/Allrun_*.sbatch` files still
name the account `special00004`. Only `profiles/slurm/config.yaml` changed.
