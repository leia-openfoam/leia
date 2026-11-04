# davof_prof: the HPCToolkit profiles of the swept-face advection (2026-11-04)

The scripts of the section "The cost of the swept-face advection: HPCToolkit
profiles" of docs/davof/REPORT-swept.md. They ran on Lichtenberg under
/work/scratch/tm83tomy/davof/agents/swept/prof with HPCToolkit 2025.0.1
(spack, ~/spack/opt/spack/linux-x86_64/hpctoolkit-2025.0.1-*).

- hp-build.sbatch: the profiling build, a detached worktree of the fork at
  davof-swept HEAD (720cfc6), WM_COMPILE_OPTION=Opt (-O3) plus
  FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer", its own platforms/.
- hp-run.sbatch: the runs (MODE=serial | par8 | par64) prepared from the
  comparison's study cases as davof-cmp-time.sbatch does; hpcrun with
  CPUTIME@1000 and ix86arch::LLC_MISSES@100000 on every rank; hpcstruct
  and hpcprof in the job.
- hpdb.py: a reader of HPCToolkit's database format v4 (meta.db,
  profile.db; FORMATS.md of the installation), used because the
  installation has no hpcviewer and no Python bindings.
- hpan.py: the analysis of one database: the flat profile (exclusive =
  HPCToolkit's "function" scope: the lexical subtree, not crossing a call or
  an inlined call), the per-line profile, the call tree, the per-rank
  table, the categories of rules_swept.py (a sample belongs to the nearest
  ancestor context matching a rule), the MPI time by call site, the callers.
- rules_swept.py: the categories (functions, source-line ranges of fork
  720cfc6, load modules).
- hp_all.sh, hp_merge.py: every database, and the CSVs
  (csv/hp_categories.csv, hp_ranks.csv, hp_hotspots.csv, hp_timing.csv on
  the cluster).
- logtimers.py, rankcells.py: the in-code timers of a comparison log and the
  per-rank partial-cell counts of a decomposed case.
