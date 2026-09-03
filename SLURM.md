# Running leia on a SLURM cluster

Site-independent recipe: preparing and running OpenFOAM cases under SLURM, and driving
them with the Snakemake workflow in `workflow/`. Nothing here names a machine, an account
or a module version — every such value is something you look up once on your cluster and
put in one place (§1.3).

- **[CLUSTER.md](CLUSTER.md)** is the worked instance of this file for Lichtenberg
  (TU Darmstadt): verified hostnames, accounts, module versions, quotas and build recipe.
  If you run there, read that too; if you run elsewhere, this file plus your site's
  documentation is what you need.
- **[CLAUDE.md](CLAUDE.md)** carries the research discipline the studies are held to.

Two recipes: **(A)** running a case by hand, **(B)** the Snakemake workflow. Read A first
even if you only use B — B automates A, and when B breaks it breaks in A's terms.

---

## 1. Environment

### 1.1 The one ordering that matters

```bash
module purge
module load <compiler-module> <mpi-module>
source $WM_PROJECT_DIR/etc/bashrc            # e.g. $HOME/OpenFOAM/OpenFOAM-v2512
. $HOME/.leia_env                            # optional session overlay, AFTER the line above
```

A session overlay that pins your own binaries typically reads

```bash
export WM_PROJECT_USER_DIR=$HOME/OpenFOAM/<session-name>
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
```

> **ORDER TRAP — cost a whole failed launch.** `$WM_OPTIONS` is defined by `etc/bashrc`.
> Source the overlay first and it expands to nothing: `FOAM_USER_APPBIN` becomes
> `.../platforms//bin`, and that broken path is **exported to every job the shell
> launches**. Jobs then fail in `blockMesh` with no informative error. Verify every time:
>
> ```bash
> echo "$WM_OPTIONS"          # must be non-empty, e.g. linux64GccDPInt32Opt
> which blockMesh <solver>    # must both resolve
> ```

Why an overlay at all: when several sessions share one account, they also share
`$FOAM_USER_LIBBIN`. Pin anything you are measuring into a session-local
`WM_PROJECT_USER_DIR` so another session's rebuild cannot silently change your binary.

### 1.2 Filesystems

Put run output on the large scratch filesystem, never in `$HOME` — `$HOME` is small,
usually shared with the build tree, and often quota'd against it. Keep code in `$HOME`
(or wherever your clone lives) and point studies at scratch.

### 1.3 Site values to look up once

| value | how to find it |
|---|---|
| login host(s) | site documentation |
| account / project | `sacctmgr show assoc user=$USER format=account,partition` |
| partitions and limits | `sinfo -o "%P %l %c %m"` |
| compiler / MPI modules | `module avail` |
| scratch path and quota | site documentation; `df -h`, `quota` |
| whether `--mem-per-cpu` is mandatory | try a submission without it (§3.2) |

Record them in `profiles/slurm/config.yaml` and your `.leia_env`, and nowhere else.

---

## 2. Building

```bash
source $WM_PROJECT_DIR/etc/bashrc && ./Allwmake
```

> **STALE-LINK TRAP.** Rebuilding a library without relinking the solvers that link it
> yields a binary that segfaults at startup with an **empty log** — indistinguishable from
> a divergence. Rebuild the library *and* relink every solver against it.

> **BUILD-SYSTEM TRAP.** If you build OpenFOAM itself, avoid `Allwmake -q`
> (`wmakeCollect`). It batches src + applications + modules into one parallel make, so a
> failing optional module stops scheduling and the standard solver *links* never run —
> you get all the libraries and none of the solvers. Disable the offending module and drop
> `-q`, so a module failure cannot poison the solvers.

**Verify the artefact, never the exit code.** `ssh host 'bash -s' < build.sh` can exit 0
having printed nothing, because ssh drops long output; a build that never ran then looks
exactly like one that succeeded silently. Redirect on the remote side and read the file
back:

```bash
ssh host 'bash -s' < build.sh          # the script redirects into /tmp/build.log itself
ssh host 'cat /tmp/build.log'
ssh host 'ls -la $FOAM_USER_APPBIN/<solver>'                      # timestamp current?
ssh host 'strings $FOAM_USER_APPBIN/<solver> | grep -c <new-symbol>'
```

---

## 3. Recipe A — a case by hand

### 3.1 Prepare, then check the mesh

```bash
cd <scratch>/myRun
blockMesh
foamDictionary -entry entry0 -keywords constant/polyMesh/boundary   # what ACTUALLY exists
rm -rf 0 && cp -r 0.org 0
setFields
decomposePar -force
```

> **PATCH TRAP — voided an entire study.** OpenFOAM errors when a **mesh** patch is
> missing from a field, but **silently ignores a field entry matching no mesh patch**. A
> `blockMeshDict` putting all four sides into one `walls` patch, while every `0.org` field
> declared `inlet` and `outlet`, produced a case whose inflow condition was parsed and
> discarded on every run it ever did — a closed box that everyone believed had a through
> flow. `constant/polyMesh/boundary` is the authority; `boundaryField` entries are a wish
> list.

### 3.2 A minimal, correct job script

```bash
#!/bin/bash
#SBATCH -A <account>
#SBATCH -J myrun-hex-N128          # NAME IT: this is how you cancel it later
#SBATCH -n 16
#SBATCH --mem-per-cpu=4000         # many sites' submit filters reject jobs without this
#SBATCH -t 1440                    # minutes
#SBATCH -o log.slurm-%j

module purge
module load <compiler-module> <mpi-module>
source $WM_PROJECT_DIR/etc/bashrc
. $HOME/.leia_env                                  # AFTER etc/bashrc

# Site-specific MPI settings, if your interconnect needs them. Get these from your
# site's working job scripts -- wrong values hang or crash on the fabric.
# export OMPI_MCA_pml=... OMPI_MCA_btl=... OMPI_MCA_mtl=...

cd "$SLURM_SUBMIT_DIR"
mpirun -np $SLURM_NTASKS <solver> -parallel > log.<solver> 2>&1
reconstructPar -latestTime                        >> log.<solver> 2>&1
```

Submit and **write the id down**:

```bash
sbatch run.sbatch | tee -a .my_jobs
```

### 3.3 The 4-rank gate, before anything reaches the cluster

Run on **at least four MPI ranks locally** and *look at the result*.

```bash
decomposePar -force              # numberOfSubdomains >= 4; a serial case proves nothing
mpirun -np 4 <solver> -parallel
```

Serial success says nothing about the code paths that exist only in parallel: processor
patch values, halo exchange, and above all **collective calls**. A collective inside a
`Pstream::master()` guard deadlocks — rank 0 blocks in the reduction, the others never
enter it — and the job then sits **alive and silent** with a truncated log and no error.
Check the step count and the metrics CSV, never the exit code. Anything touching an
`fvMatrix` (an `fvOption`, a new term, a new solver) must clear this gate as part of its
inertness check, not as an extra.

---

## 4. Recipe B — the Snakemake workflow

One study = one `(case, mesh, mode)` swept over parameter axes. The backend is chosen
purely by the profile — `profiles/local*` runs locally, `profiles/slurm` submits one
`sbatch` per case. Nothing else in the study changes.

### 4.1 A study config

`config/<study>.yaml`:

```yaml
study_name: myStudy2D
case: translatingDroplet2D          # a directory under cases/
mesh: hex
mode: parallel
np: 16                              # ranks per arm
solver: leiaSemiLagrangianLevelSetTwoPhaseFoam
theme: method-comparison            # where curated output lands under docs/
setfields_args: "-alphaName alpha.water"
solve_runtime: 1380                 # minutes -> the SLURM time limit
axes_override:
  N_CELLS: [128, 256]               # every axis is a list; their product is the arm set
  END_TIME: [0.02]
collapse_other_axes: true
mpi_launcher: "mpirun -np {np}"
env_preamble: "source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc; [ -f $HOME/.leia_env ] && . $HOME/.leia_env || true"
```

Case directories hold `*.template` files with `@!TOKEN!@` placeholders. The materializer
renders every template under the case, substituting axis values and the defaults from
`cases/default.parameter`. `workflow/README.md` documents the committed studies.

**Write the header before you launch**: the prediction, the pass/fail threshold, and what
result would *falsify* the idea. Then the data decides rather than the reader.

**One variable per comparison, plus controls.** Audit every arm for a second difference
before launching. Arms are matched on mesh, time step *and its law* (`dt ~ h^p` hides an
h-dependence inside a dt sweep), and horizon. Put a null arm in the same matrix so it
validates itself.

### 4.2 Dry run, always

```bash
snakemake --workflow-profile profiles/slurm --configfile config/myStudy2D.yaml -n
```

Check the `solve` count equals the arm count you intended. A silent factor of two is
usually an axis you forgot to collapse.

### 4.3 Launch

```bash
cd <clone>
module purge; module load <compiler-module> <mpi-module>
source $WM_PROJECT_DIR/etc/bashrc
. $HOME/.leia_env
git pull --ff-only origin development
ls config/myStudy2D.yaml || exit 1          # fail loudly if the pull did not land
nohup snakemake --workflow-profile profiles/slurm --configfile config/myStudy2D.yaml \
      > studies_myStudy2D.log 2>&1 &
```

> **PATH TRAP — three failed submissions.** `mpirun` resolves the executable from the
> **orchestrating driver's PATH**, not from anything in the job script. Launch the driver
> from a shell that has not sourced your overlay and every rank silently runs whatever
> older binary sits on the default path. The one-line check:
>
> ```bash
> grep -m1 '^Exec' studies/myStudy2D/*/log.<solver>
> ```

> **PULL TRAP.** `git pull --ff-only` **aborts** when untracked or modified files would be
> overwritten. That is correct behaviour, and it means the config you are about to run may
> not be present. Always `ls` it after pulling and abort if missing, rather than launching
> against a stale clone. If files block the pull, move them aside — do not delete — and
> check whether they differ from the remote first.

### 4.4 Record what you started

```bash
squeue -h -o "%i|%Z" | awk -F'|' -v c="$PWD" '$2 == c {print $1}' | sort -u >> .my_jobs
```

---

## 5. Cancelling: only your own jobs, ever

**`scancel -u $USER` is forbidden** whenever an account is shared — an account-wide cancel
destroys other people's running work, not just yours. Measured: a 2x2x2 matrix was killed
19 minutes in by exactly this, together with five jobs belonging to another session;
`sacct` recorded the cancelling uid, which reads deceptively like an administrator action.

```bash
scancel <id> [<id> ...]          # explicit ids from your ledger: always safe
scancel -n myrun-hex-N128        # by JOB NAME, if the name is yours
```

Identify your work by **job name or id**, never by user (`squeue -u $USER` lists every
session's jobs on a shared account). A job you did not submit is not yours to kill even
when it looks stale. To replace a driver, cancel it by id and resubmit rather than
clearing the queue.

The same discipline applies to local cleanup: never `pkill -f <pattern>` (it matches its
own shell — exit 144) and never sweep by user. A background task from an earlier turn has
a different parent PID than the current shell, so "its parent is not my shell" does **not**
mean "not mine". Query your ledger by id.

---

## 6. Reading results

### Never write a new grep against a solver log

Every OpenFOAM solver prints, in its startup banner:

```
trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
```

so a grep for a bare `Floating point exception` — or bare `core dumped`, which SLURM
attaches to unrelated kill text — classifies **every healthy run as diverged**. This false
positive has bitten three times in this repository and once deleted a live 16-hour run.

```bash
workflow/scripts/foam_log_state.sh <log> [--stall S] [--wait [--poll S]]
```

returns `COMPLETED / RUNNING / STALLED / DIVERGED / LAUNCH_FAILURE / MISSING` with distinct
exit codes plus `steps=` and `age=`, using exactly the Snakefile's patterns. `--wait`
replaces hand-rolled completion loops. If the patterns ever change, change them in **both**
the helper and `workflow/Snakefile` in the same commit.

- **Completeness** comes from the terminating `^End$` or the solver's per-step CSV — never
  from a post-processing output that only exists after a run finishes.
- **Liveness** comes from log mtime and step count, never from `squeue` alone: an empty
  `squeue` during a controller outage reads exactly like "all jobs finished".
- A `DIVERGED` run is a **result** (a blow-up is data); a `LAUNCH_FAILURE` is not.

### Every waiter carries a timeout, and its exit condition must be reachable

An `until <cond>; do sleep N; done` with no timeout does not fail when its condition never
arrives — it polls until the session ends. Three such waiters were once still spinning
hours after the work they watched had finished. Two more rules, each from an orphan this
repository produced: **never grep for a pattern your own command line contains** (a waiter
conditioned on `! pgrep -f "configfile config/X"` can never exit, because the `pgrep`
matches the shell running it), and **ask whether the pattern can ever appear** — ssh may
have swallowed the output you are waiting on.

### Bringing output back

Code moves by git; raw simulation output moves by rsync.

```bash
make pull-runs                       # everything
make pull-study STUDY=myStudy2D      # one study
```

For a field inspection or an animation, pull only the fields you need:

```bash
rsync -a --include='0.0*/' --include='0.0*/U' --include='0.0*/alpha.water' \
      --include='0.0*/psi' --exclude='*' \
      <host>:<scratch>/studies/X/arm/ ./local/
```

### Preserve before you overwrite

Raw output that a re-run will overwrite is preserved **first** — rename the study
directory with a dated suffix; solvers truncate their metrics CSV on restart. Destructive
cleanup never runs while any driver may be alive, and liveness comes from
`foam_log_state.sh`, never from a declared output and never from `squeue` alone.

**A wrong setup voids its data.** When a case is found to be set up wrong, every number it
produced is void — including the ones that look unaffected. Rename with a
`_VOID_<reason>_<date>` suffix so it can never be curated by accident, fix, and re-run from
scratch. Do not reason about which metrics "should still be valid": that reasoning is
exactly as reliable as the setup was, and a partially-trusted table is worse than none.

---

## 7. Checklist before you submit

1. `echo $WM_OPTIONS` non-empty; `which <solver>` resolves; binary timestamp current.
2. `foamDictionary -entry entry0 -keywords constant/polyMesh/boundary` — the mesh has the
   patches the fields name.
3. Ran on **4 ranks** locally and the *result* was looked at, not the exit code.
4. `--mem-per-cpu` set (if your site requires it); job **named**.
5. `snakemake -n` shows the arm count you expect.
6. The config header states the prediction and what would falsify it.
7. Job ids appended to `.my_jobs`.
8. After launch: `grep -m1 '^Exec' <log>` confirms the binary every rank is running.
