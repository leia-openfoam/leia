# Running leia on Lichtenberg (TU Darmstadt)

Canonical, agent-readable description of the laptop <-> cluster workflow. All
values below are **verified** against the live cluster (2026-07-27).

The KISS rule that dissolves every "two machines, one branch" headache:

> **Git carries code. rsync carries raw output. Neither carries the other.**
>
> The laptop (WSL) and Lichtenberg are two *equal clones* of the GitHub repo
> (`git@github.com:leia-openfoam/leia.git`). GitHub is the only hub — you never
> `scp` code and never `git push` into the cluster. Heavy simulation output
> (`studies/`, `runs/`) is git-ignored and lives on the cluster; it comes to the
> laptop by `rsync`, on demand, only to look at it.

## Why there is no "two commits for one milestone" problem

Laptop and cluster commit **disjoint paths**, so git merges them and
`git pull --rebase` keeps history linear:

| Machine    | Role                                    | Commits these paths                                           |
|------------|-----------------------------------------|---------------------------------------------------------------|
| Laptop/WSL | method dev (with Claude), viewing decks | solver/library `src/`, `*.template.html`, prose, plot scripts |
| Lichtenberg| heavy runs + BO orchestration           | regenerated results `docs/**/data/{figures,tables,mechanism}` |

`docs/**/data/*.{csv,png,svg}` are the small tracked results the decks + papers
consume (re-included past the global `*.csv`/`*.png` ignore). Built deck HTML is
git-ignored — rebuild it locally in seconds with `make docs`. So slides are
*authored* on the laptop but *fed* by data the cluster committed. Discipline:
**`git pull --rebase` before you start and before you push.**

## Verified cluster environment (Lichtenberg / HHLR)

- **User / login:** `tm83tomy@lclusterN.hrz.tu-darmstadt.de`. Helper `~/bin/licht N`
  opens `ssh -vvX tm83tomy@lclusterN...`. Login nodes: **1–7** (Phase 2 / i02),
  **13–19** (Phase 1 / i01); any is fine for submission. Passwordless (SSH key).
- **SLURM account:** `special00004` (current default for runs). Also available:
  `special00005`, `project00186/00450/00524/00727/01204/01456`.
- **SLURM gotcha:** the `job_submit` plugin **rejects any job without
  `--mem-per-cpu`**; partition is auto-routed from the time limit (default
  `deflt`, 24 h; `deflt_short` 30 min; `long` 7 d). The snakemake slurm profile
  supplies `--mem-per-cpu` via `mem_mb_per_cpu`.
- **Filesystems:** `$HOME` (=/home/tm83tomy) and `/work/home/tm83tomy` share a
  268 G quota (lots free); `/work/scratch/tm83tomy` has ~40 T for run output;
  the group tree is `/work/groups/da_mma_b`.
- **Toolchain modules:** `gcc/11.5.0-z7mc` (also 13.4.0, 14.3.0),
  `openmpi/4.1.8-6xzv`. `cmake 4.2`, `flex 2.6.4`, `bison 3.7.4`, `git 2.52`.
- **Login `.bashrc` noise:** it tries to `module load gcc/11.4.1 python/3.11.9`
  which don't exist on every node — harmless Lmod errors. Batch jobs ignore it
  (`module purge` + explicit loads in the job).
- **OpenFOAM:** built from source in **`$HOME/OpenFOAM/OpenFOAM-v2512`**
  (matches the WSL laptop version). Defaults `WM_COMPILER=Gcc`,
  `WM_COMPILER_TYPE=system`, `WM_MPLIB=SYSTEMOPENMPI` — uses the module gcc +
  module OpenMPI, no ThirdParty compiler/MPI build. A system spack module
  `openfoam/2512` also exists but does NOT populate the classic `wmake`
  environment, so we use the source build.

### How OpenFOAM-v2512 was built (reproduce with `$HOME/OpenFOAM/build-of2512-rebuild.sbatch`)

```bash
# on a login node (internet only exists here):
cd $HOME/OpenFOAM
curl -L -o OpenFOAM-v2512.tgz     https://dl.openfoam.com/source/v2512/OpenFOAM-v2512.tgz
curl -L -o ThirdParty-v2512.tar.gz https://dl.openfoam.com/source/v2512/ThirdParty-v2512.tar.gz   # NOTE: .tar.gz
tar xzf OpenFOAM-v2512.tgz && tar xzf ThirdParty-v2512.tar.gz
mv OpenFOAM-v2512/modules/OpenQBMM OpenFOAM-v2512/OpenQBMM.disabled   # see gotcha below
sbatch build-of2512-rebuild.sbatch   # -A special00004, -c 48, --mem-per-cpu=3600
```

The job does: `module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv`,
`source etc/bashrc`, then `./Allwmake -j 48 -s -l`. Result (verified 2026-07-27):
**270 apps, 128 libs, `foamInstallationTest` = "Critical systems ok".**

> **BUILD GOTCHA (cost one failed build).** Do **not** use `Allwmake`'s `-q`
> (queue / `wmakeCollect`) mode here. `-q` batches src+applications+modules into
> one parallel make, so when the bundled **OpenQBMM** module (population balance,
> unused by leia) fails a compile, `make` aborts scheduling and the standard
> solver *links* never run — you get all 128 libraries but `icoFoam`/`interFoam`
> missing. Fix = disable OpenQBMM **and** drop `-q` (build src -> applications ->
> modules in order, so a module failure can't poison the solvers). The other
> bundled modules (adios, visualization, external-solver) are optional and their
> ADIOS2/CGAL warnings are harmless.

### Runtime MPI (proven on the interconnect)

OpenMPI on Lichtenberg wants these (from the group's working job scripts); the
slurm profile exports them:

```bash
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
```

## Where leia lives on the cluster

**`/work/scratch/tm83tomy/leia`** — the checkout is on the **parallel file
system**, because that is where simulations must run (fast parallel I/O).
`/work/scratch` is periodically **purged**; that is fine because the repo is only
ever a clone of the GitHub hub and its heavy output is regenerable — if scratch
is wiped, re-clone. The OpenFOAM install and leia's compiled binaries
(`$FOAM_USER_APPBIN`) live in the persistent `$HOME`, so only the source + run
output are at risk.

## GitHub SSH from the cluster (git protocol, so push+pull work)

The cluster authenticates to GitHub with its **own** key
`~/.ssh/id_github_ed25519` (pinned for github.com in `~/.ssh/config` with
`IdentitiesOnly yes`, so the ~15 other keys in `~/.ssh` don't cause
"too many authentication failures"). The matching public key is registered on
GitHub. Regenerate + re-register only if it is ever lost:
```bash
ssh-keygen -t ed25519 -C "tm83tomy@lichtenberg-github" -f ~/.ssh/id_github_ed25519 -N ""
cat ~/.ssh/id_github_ed25519.pub   # add to GitHub -> Settings -> SSH keys
```

## One-time laptop setup

`~/.ssh/config` (so the Makefile's `CLUSTER=lichtenberg` and rsync helpers work):

```
Host lichtenberg
    HostName lcluster5.hrz.tu-darmstadt.de
    User tm83tomy
    IdentityFile ~/.ssh/id_ed25519
```

Clone on the cluster once (into scratch, over SSH):

```bash
ssh lichtenberg
cd /work/scratch/tm83tomy
git clone git@github.com:leia-openfoam/leia.git && cd leia
git checkout feature/velocity-extension
git submodule update --init --recursive        # pyFoamStudy (legacy)
python3 -m pip install --user --break-system-packages "snakemake>=8" snakemake-executor-plugin-slurm
source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake     # build leia against v2512
```

## Parallel-run status (verified 2026-07-28)

A single reversed 2D-vortex case with quadratic semi-Lagrangian advection
(`quadraticWeightedLeastSquares`) **runs to completion on 4 MPI ranks** on
Lichtenberg (`mpirun -np 4 leiaSemiLagrangeLevelSetFoam -parallel`, OpenFOAM-v2512
+ leia, account `special00004`). The proven pattern is a **standalone `#SBATCH -n 4`
job** (like the group's `isoAdv.sbatch`): `blockMesh -> decomposePar -force ->
mpirun -np {np} <solver> -parallel -> reconstructPar`.

### RESOLVED 2026-08-08 — parallel studies now run under `profiles/slurm`

The former gap ("`make studies PROFILE=profiles/slurm` does NOT drive parallel
solves") is **fixed**. The launcher in `profiles/slurm/config.yaml` is now:

```yaml
  - "mpi_launcher=srun --ntasks={np} --overlap"
```

Diagnosis, because the symptom points the wrong way. snakemake's
`slurm-jobstep` executor wraps every rule in an `srun` step with **one task**.
The group sbatch itself is allocated correctly — measured `NCPUS=4, ReqCPUS=4` —
but `mpirun -np 4`, launched *inside* that one-task step, sees a single slot and
aborts with *"There are not enough slots available in the system to satisfy the
4 slots that were requested"*. Because a plain `sbatch -n 4` runs the identical
case fine (verified: `SLURM_NTASKS=4`, whole case in < 60 s), this reads as an
allocation problem and is not one.

`--overlap` lets the inner step share the allocation rather than demand
exclusive resources — which is also what produced the original 2026-07-28
"CPU binding outside of job step allocation" failure with bare `srun`.

**The launcher alone was not enough.** With every rule still bundled into one
`group: "case"` sbatch, the group's allocation did not reliably carry the
solve's `tasks={np}`, and ~4% of jobs died with *"srun: error: Unable to create
step ...: More processors requested than permitted"* — scattered across all
arms, including `noSource` at N=32, so plainly not physics. The structural fix
is the rule layout below.

### Serial steps and the parallel step are separate jobs

Only the CFD solve is an MPI job:

| rule | parallelism | group | what it does |
|---|---|---|---|
| `mesh` | serial | `case_pre` | `blockMesh` / `pMesh` |
| `preprocess` | serial | `case_pre` | `leiaSetFields` on the **undecomposed** case |
| `decompose` | serial | `case_pre` | `decomposePar -force` |
| **`solve`** | **MPI, `tasks={np}`** | **ungrouped** | the solver — its own sbatch |
| `reconstruct` | serial | `case_post` | `reconstructPar -withZero` |
| `paraview_stub` | serial | `case_post` | `.foam` marker |

This is the standard OpenFOAM order — `blockMesh -> setFields -> decomposePar ->
solver -> reconstructPar` — and `leiaSetFields` running serially before
decomposition is also free of processor-boundary artefacts in the signed
distance and the narrow band. Because `solve` is ungrouped it gets its own
allocation with exactly `--ntasks={np}`, so the launcher always has the ranks it
asks for; the cheap serial steps still share one sbatch per case at each end.

#### The solve must also be pinned to ONE node (`resources: nodes=1`)

Ungrouping fixed np=4 but **not** np=8. SLURM satisfied `--ntasks=8` across
**two** nodes, so the `.batch` step held only 7 CPUs and `srun --ntasks=8`
inside it failed:

```
srun: error: Unable to create step for job 54042738:
      More processors requested than permitted
```

Measured allocations for that one rule before the pin: 8, 7, 6, 5 and 4 CPUs.
After `resources: nodes=1`: `ReqCPUS=8, NNodes=1` on every solve, and the study
went 12/18 -> 18/18.

> **Do NOT write this as `slurm_extra='--nodes=1'`.** The executor plugin
> rejects it — *"The --number-of-nodes option is not allowed in the
> 'slurm_extra' parameter. The number of nodes is set by the snakemake executor
> plugin"* — and the rejection is a **submission-time exception**, so no solve
> job is created at all and the study simply stalls with zero allocations. Use
> the plugin's own `nodes` resource.

#### The inner srun must also be given the CPUs (`--cpu-bind=none`)

**Added 2026-08-15, and it is the single largest performance fact on this page.**
`--overlap` lets the inner step share the outer jobstep's allocation -- including
its **CPU mask**. snakemake's slurm-jobstep wrapper creates that outer step with
one task, so every rank of the inner `srun --ntasks={np}` inherited a **1-CPU**
mask. OpenMPI's busy-wait progress engine then collapses: each rank spins its
whole timeslice waiting for ranks that are not scheduled.

There is **no error message**. The run is simply ~4 orders of magnitude slower,
and the CPU time looks plausible because spinning *is* CPU time.

Measured on `stationaryDroplet2D` N=64, np 8, everything else identical:

| launcher | steps/s | wall/CPU |
|---|---|---|
| `srun --ntasks=8 --overlap` | 0.0077 | 9.05 |
| `srun --ntasks=8 --overlap --cpu-bind=none` | **86** | **1.08** |

> **The diagnostic is the RATIO.** OpenFOAM prints `ExecutionTime` (CPU) and
> `ClockTime` (wall) every step. `wall/CPU ~= np` means the ranks are sharing one
> core; `~1` means they are not. Check this on any new cluster study before
> trusting its pace, and size time limits only from a ratio-1 measurement.

This is why a plain `sbatch -n 8` was always fast (60 steps in 1 s) while the
same case under snakemake needed 4 h for 111 steps -- which reads as a cluster
load or filesystem problem and is neither. Studies run under `profiles/slurm`
before 2026-08-15 were all affected; their results are unchanged (it is a pace
bug, not a numerics one) but their measured runtimes mean nothing.

#### A diverged solve is a result; a launch failure is not

`aggregate` needs every case's CSV, so one divergent arm used to leave a study
with no curated table. The solve rule therefore records a non-zero solver exit
(`.leia_solver_failed`, a `solverFailed` CSV column, blank metrics) instead of
stalling — **but only when the solver actually ran**. If the per-case
`log.<solver>` carries a launcher error the rule fails loudly with no marker and
no touched output. Without that discrimination the np=8 launch failures above
were logged as 11 "divergences", including the `euler` baseline at N=32, which
cannot fail for physical reasons. Infrastructure faults must never enter the
record as physics.

Measured on `sdplsStability` (18 parallel np=4 cases, one command):

| launcher | solver CSVs produced |
|---|---|
| `mpirun -np {np}` | **0 / 18** |
| `srun --ntasks={np} --overlap` | **18 / 18** |

`profiles/local` keeps `mpirun -np {np}` and needs no change — there is no
`srun` nesting there. The same study therefore runs both ways, switching only
`PROFILE`:

```bash
make studies-sdpls PROFILE=profiles/local     # laptop / WSL
make studies-sdpls PROFILE=profiles/slurm     # Lichtenberg, one sbatch per case
```

**On the cluster, submit the orchestrator as a job — do not run it on a login
node.** A login-node `make studies-euler` was SIGTERMed after ~5.5 h mid-sweep
(login nodes reap long-lived processes), losing the driver while its submitted
case jobs kept running:

```bash
sbatch run-studies.sbatch                                   # TARGET=studies-euler
sbatch --export=ALL,TARGET=studies-sdpls run-studies.sbatch
```

The orchestrator takes one core and only submits work. Resuming is free —
snakemake skips cases whose outputs exist — so re-submitting after any
interruption continues where it stopped.

### Study groups, and when they must be re-run together

`make studies-euler` runs **every** study whose level-set transport is a
finite-volume `div(phi,psi)` — the SDPLS line plus the GRL, velocity-extension
and remaining comparison lines. They share one discretization
(`div(phi,psi) Gauss linearUpwind grad(psi)`, `grad(psi) cellLimited leastSquares 1`,
`nDefCorr >= 3`), so **whenever that changes they must be re-run together**:
mixing discretizations inside one comparison table is exactly how the 2D/3D
SDPLS confound happened (2D ran `linearUpwind grad(psi)`, the newer 3D configs
ran `linear`, and nothing in the curated CSV recorded the difference).

The ~82 semi-Lagrangian studies have no `div(phi,psi)` term at all — that update
is a pointwise assignment — so they are unaffected and are not in the group.

What each run actually used is now recorded per case in the curated CSV
(`divPsi`, `divPsiGrad`, `divPsiGradScheme`, `gradPsiSdpls`, `gradUSdpls`,
`nDefCorr`), parsed from the **rendered** `system/fvSchemes` rather than from
study tokens — see `workflow/scripts/fvschemes.py` for why the tokens are not
trustworthy (2Dvortex hardcodes its div scheme and has no `DIV` token at all).

## The daily loop

```bash
# --- laptop: develop with Claude, then hand off ---
git pull --rebase && <edit / build / smoke-test in WSL> && git commit -am "..." && git push

# --- cluster: pull + run (this is all "run on Lichtenberg" means) ---
ssh lichtenberg
cd /work/scratch/tm83tomy/leia && git pull --rebase && ./Allwmake
make studies PROFILE=profiles/slurm                       # one sbatch per case
git commit -am "results: <study>" docs/**/data && git push  # only the tracked data

# --- laptop: see the new results ---
git pull --rebase && make docs        # rebuild decks + papers from docs/**/data
make pull-runs                        # OR rsync raw fields back to visualise
```

`make studies` swaps backend purely via `PROFILE` — identical study, local
(`mpirun`, WSL smoke test) vs Lichtenberg (`sbatch`/case).

## The three workflows

**(3) Method development (interactive, with Claude).** Edit + commit on the
laptop; smoke-test in WSL with `profiles/local`. Push. On the cluster
`git pull && ./Allwmake`, run the real resolution with `profiles/slurm`. `rsync`
raw fields back only to visualise (`make pull-runs`).

**(2) Automated BO runs.** Snakemake + SLURM fans out one iteration's candidate
simulations as independent `sbatch` jobs. Launch inside `tmux` on a login node so
it survives disconnect:
```bash
ssh lichtenberg && tmux new -s bo
cd /work/scratch/tm83tomy/leia && python3 workflow/scripts/bo_driver.py
```

**(1) Sequential BO (each sim triggered by the previous argmax).** BO is
inherently sequential, so the **outer loop is a thin Python driver**, not the
snakemake DAG. Per iteration: render candidate -> `snakemake --workflow-profile
profiles/slurm` (submit + wait) -> read objective from the study database CSV ->
fit GP -> argmax -> repeat. Snakemake still runs the inner simulation (reusing
caching/resume/launcher); it just doesn't decide the next point. DAG for fan-out,
Python for the feedback loop.

## rsync helpers (Makefile)

Raw output is git-ignored and stays on the cluster; pull it only to visualise:

```bash
make pull-runs                    # all of studies/ from the cluster
make pull-study STUDY=bulkVortexSL
```

`CLUSTER` (default `lichtenberg`) is the ssh alias; `REMOTE` the repo path.
Everything that feeds a deck or paper already travels via git as `docs/**/data`.

> **Check the size before pulling a 3D study.** These targets rsync the whole
> tree, raw fields included. `sdplsConv3Dshear` is **195 GB** — 17 GB for the
> single 203³ case — while everything the tables and figures need is **9 MB** of
> CSV. `make pull-study` there is a mistake, not a slow success.
>
> Aggregate in place and pull only the curated result:
>
> ```bash
> ssh $CLUSTER "cd $REMOTE && python3 -c \"
> import sys, glob; sys.path.insert(0,'workflow/scripts'); import aggregate
> d = sorted(glob.glob('studies/<study>/<case>_[0-9]*'))
> aggregate.build_database(d, 'studies/<study>/<study>_database.csv')\""
> rsync $CLUSTER:$REMOTE/studies/<study>/'<study>_*.csv' studies/<study>/
> ```
>
> Pull raw fields only for the one case you actually intend to visualise. The 3D
> iso-surface montages are generated locally regardless (no vtk on the cluster).
