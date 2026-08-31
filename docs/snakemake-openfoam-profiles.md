# Snakemake profiles for OpenFOAM: local + SLURM from one workflow

A portable recipe for driving OpenFOAM parameter studies with Snakemake so that
**the identical study runs on a laptop and on a SLURM cluster, switching backend
with `--workflow-profile` alone.** No per-study `.sbatch` scripts, no
`if cluster:` branches in the workflow.

Written for an agent porting this into another project. Everything marked
MEASURED is an observed failure in the `leia` repository (OpenFOAM-v2512,
Lichtenberg / TU Darmstadt, snakemake with `snakemake-executor-plugin-slurm`),
with the number that was observed. The failure modes are the valuable part of
this document — the YAML is short and obvious, and every non-obvious line in it
exists because one of these bit somebody.

---

## 1. The design in one sentence

The Snakefile knows **nothing** about SLURM. It declares, per rule, *how many
MPI ranks the rule needs* and *how to launch an MPI job*; a profile supplies the
executor and the launcher string. Backend choice is therefore a command-line
flag:

```bash
snakemake --workflow-profile profiles/local --configfile config/<study>.yaml
snakemake --workflow-profile profiles/slurm --configfile config/<study>.yaml
```

The payoff is not convenience, it is **comparability**: a study smoke-tested on
4 local ranks and then run at production resolution on the cluster is the same
DAG with the same commands, so a discrepancy between them is a real result and
not a difference in the driver.

---

## 2. The four things the Snakefile must do

Everything else in this document follows from these. Get them wrong and the
profiles cannot fix it.

### 2.1 Every rule declares `resources: tasks=N`

`tasks` is the MPI rank count. Serial rules declare `tasks=1`; the solve rule
declares `tasks=<np>`.

```python
NP       = int(config["np"])
PARALLEL = config["mode"] == "parallel"

rule solve:
    resources:
        tasks=(NP if PARALLEL else 1),
```

This single declaration does double duty:

* **On SLURM** the executor turns it into `--ntasks=N` on the generated `sbatch`.
* **Locally** it is a share of a global budget, which is what stops the machine
  being oversubscribed (§3).

### 2.2 The MPI launcher is a config string, not a hardcoded command

```python
MPI        = config.get("mpi_launcher", "mpirun -np {np}")
MPI_PREFIX = (MPI.format(np=NP) + " ") if PARALLEL else ""
PAR        = "-parallel" if PARALLEL else ""

# ... in the rule's shell:
#   f"{MPI_PREFIX}{SOLVER} {PAR}"
```

Note the `{np}` placeholder. `mpirun` needs an explicit rank count; `srun` does
**not**, because the SLURM allocation already carries `--ntasks`. Keeping it a
format string lets one mechanism serve both. This is the seam the whole design
turns on — see §4.3, where getting the launcher wrong cost a factor of 11 000.

### 2.3 An environment preamble runs at the top of every shell

OpenFOAM is a shell environment, not a set of binaries on `PATH`, and it must be
sourced inside each job because a cluster job does not inherit the submitting
shell.

```python
PREAMBLE = (config.get("env_preamble") or "").strip()

def sh(*lines):
    """Shell body: source OpenFOAM (lenient), then run the commands strictly."""
    pre = (PREAMBLE + "\n") if PREAMBLE else ""
    return "set +eu\n" + pre + "set -e\n" + "\n".join(lines)
```

`set +eu` around the preamble and `set -e` after it is **not stylistic**.
OpenFOAM's `etc/bashrc` reads unset variables (`WM_PROJECT_DIR` among them), so
sourcing it under `set -u` aborts every job before it starts. The commands you
actually care about still run strictly.

### 2.4 Serial and parallel steps are SEPARATE jobs

Group the cheap serial steps; leave the MPI solve **ungrouped**.

```python
rule mesh:        # blockMesh / pMesh
    resources: tasks=1
    group: "case_pre"

rule preprocess:  # field initialisation (setFields and friends)
    resources: tasks=1
    group: "case_pre"

rule decompose:   # decomposePar -force  -- SERIAL, on the undecomposed case
    resources: tasks=1
    group: "case_pre"

rule solve:       # THE ONLY MPI RULE -- deliberately ungrouped
    resources: tasks=NP, nodes=1

rule reconstruct: # reconstructPar -withZero
    resources: tasks=1
    group: "case_post"
```

Why this matters is §4.1. It is the first thing that goes wrong.

Two OpenFOAM-specific notes on the ordering above. **Field initialisation runs
before `decomposePar`, serially, on the whole case** — that is the standard
OpenFOAM order and it is also free of parallel artefacts, since a signed
distance field or a narrow band gets built with no processor-boundary exchange
at all. And `decomposePar` itself is serial: it is the step that *creates* the
parallel case, so it cannot be one.

---

## 3. `profiles/local/config.yaml`

```yaml
# Local profile: run every case on this machine, MPI via mpirun.
#   snakemake --workflow-profile profiles/local --configfile config/<study>.yaml
executor: local
cores: 24           # CPU cores Snakemake may use at once
jobs: 6             # max concurrent Snakemake jobs
resources:
  - tasks=24
printshellcmds: true
rerun-incomplete: true
```

The launcher and preamble come from the base `config/config.yaml`:

```yaml
mpi_launcher: "mpirun -np {np}"
env_preamble: "source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
```

**The `tasks=24` line is the whole trick.** Because `solve` declares
`resources: tasks=<np>`, a global budget of 24 bounds concurrent MPI solves to
`floor(24/np)` — automatically, for any `np`. At np=4 six solves run; at np=8,
three. The machine is never oversubscribed *regardless of what the study sets*,
which is why plain Snakemake suffices with no launcher script.

It also gives a one-flag override for a memory-bound tier, no config edit
needed:

```bash
snakemake --workflow-profile profiles/local --configfile config/<study>.yaml --resources tasks=8
```

### Provide more than one local profile

Fixed core budgets beat editing one file, especially when several sessions share
a workstation. The real profiles are `local` (24), `local20` (20 cores / 5
concurrent arms, leaving 4 free for interactive use) and `local8` (8 cores,
`jobs: 1` — one arm at a time, for when another session holds the rest of the
box).

Two hard-won constraints belong in those comments:

* **Memory, not cores, is usually the binding limit.** On a 15 GB laptop
  (~13 GB usable) the 3D ladder is impossible — 95³ is 857k cells, 120³ is
  1.73M — so laptop studies are 2D or small 3D. Sizing a local profile by core
  count alone gets you an OOM kill at the top rung.
* **On a heterogeneous CPU (Intel P+E cores, no SMT), keep `np` small and EQUAL
  across the arms of one study.** An MPI job runs at the pace of its slowest
  rank, so unequal `np` makes arms incomparable — and with `jobs > 1`, two
  concurrent solves contend for memory bandwidth, at which point the slower
  arm's wall time is not a property of the discretisation you are measuring.
  If wall time is part of the result, use `jobs: 1`.

---

## 4. `profiles/slurm/config.yaml`

```yaml
executor: slurm
jobs: 200                  # max sbatch jobs in flight
printshellcmds: true
rerun-incomplete: true

default-resources:
  slurm_account: "special00004"
  mem_mb_per_cpu: 4000
  runtime: 1440            # minutes (24 h)

set-resources:
  mesh:          {runtime: 120, mem_mb_per_cpu: 16000}
  preprocess:    {runtime: 120, mem_mb_per_cpu: 16000}
  decompose:     {runtime: 240, mem_mb_per_cpu: 24000}
  reconstruct:   {runtime: 240, mem_mb_per_cpu: 24000}

config:
  - "mpi_launcher=srun --ntasks={np} --overlap --cpu-bind=none"
  - "env_preamble=module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv; export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2; source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
```

Install the plugin first: `pip install snakemake-executor-plugin-slurm`.

The `config:` block **overrides the base config's values**, which is how one
Snakefile gets a different launcher and a different environment per backend
without a single conditional in the workflow.

### 4.1 Never group the MPI solve with serial steps

MEASURED. With `mesh`, `decompose` and `solve` all in one group, Snakemake
submits **one** `sbatch` for the group, and that allocation did not reliably
carry the solve's `tasks={np}`. The inner jobstep then failed with

```
srun: error: Unable to create step for job ...: More processors requested than permitted
```

or, with `mpirun`, *"There are not enough slots available in the system to
satisfy the {np} slots that were requested"*. Ungrouped, the solve gets its own
`sbatch` with exactly `--ntasks={np}` and the launcher always has the ranks it
asked for. The cheap serial steps still share one job at each end, so you do not
pay a scheduling round-trip per `blockMesh`.

### 4.2 Pin the solve to ONE node

MEASURED. SLURM satisfied `--ntasks=8` across **two** nodes, leaving the batch
step with only 7 CPUs, and `srun --ntasks=8` inside it failed to create a step.
Observed allocations for the same np=8 rule: 8, 7, 6, 5 and 4 CPUs. The ranks
must share a node for the jobstep to be creatable:

```python
rule solve:
    resources:
        tasks=NP,
        **({"nodes": 1} if PARALLEL else {}),
```

This **must** be the plugin's own `nodes` resource, not `slurm_extra`. The
executor rejects `--nodes` in `slurm_extra` outright, and the rejection is a
submission-time exception, so no job is created at all.

### 4.3 The launcher: `srun --ntasks={np} --overlap --cpu-bind=none`

Three failures in opposite directions produced this line. Read it before
choosing a launcher on a new cluster.

* **Bare `srun`** failed with *"CPU binding outside of job step allocation /
  Unable to satisfy cpu bind request"*. So the profile switched to `mpirun`.
* **`mpirun -np {np}`** then failed for *every* np > 1 study with *"not enough
  slots"* — even though the `sbatch` was allocated correctly (measured
  `NCPUS=4, ReqCPUS=4`). Cause: Snakemake's `slurm-jobstep` executor wraps every
  rule in an `srun` step with **one task**, and `mpirun` launched inside that
  step sees a single slot. A plain `sbatch -n 4` ran the identical case fine,
  which is exactly what made this look like an allocation problem rather than a
  nesting one. MEASURED on 18 np=4 cases: 0 solver CSVs with `mpirun`, 18 with
  `srun --overlap`.
* **`--cpu-bind=none`, and this is the expensive one.** `--overlap` lets the
  inner step share the outer step's resources — but it shares its **CPU mask**
  too, and that outer step has a single task. So all `{np}` ranks inherited a
  1-CPU mask, and OpenMPI's busy-wait progress engine collapsed: every rank
  spins its full timeslice waiting for ranks that are not scheduled.

  **There is no error message.** The run is simply slower, and the CPU time
  looks plausible because spinning *is* CPU. MEASURED, np=8, everything else
  identical:

  | launcher | steps/s | wall/CPU |
  |---|---|---|
  | `srun --ntasks=8 --overlap` | 0.0077 | 9.05 |
  | `srun --ntasks=8 --overlap --cpu-bind=none` | **86** | **1.08** |

  A factor of 11 000, silent, for months.

**The diagnostic to build into your workflow:** OpenFOAM prints `ExecutionTime`
(CPU) and `ClockTime` (wall) every step. **`wall/CPU ≈ np` means the ranks are
sharing one core; `≈ 1` means they are not.** Check this on the first parallel
job on any new cluster, before trusting a single timing.

Whether you need `--overlap` depends on whether your Snakemake plugin version
nests `srun`. Verify rather than copy: submit one np=4 job, read the generated
script, and check the ratio above.

### 4.4 Size the serial rules by MEMORY, not by their cost

Least obvious item here, and it fails only at the top of a ladder — after the
cheap rungs have all passed.

The serial rules are `tasks=1`, so `mem_mb_per_cpu` **is** the whole job's
memory, and being serial they do **not** get cheaper as the solve's `np` rises.
`blockMesh`, field initialisation and a scotch `decomposePar` each hold the
**entire mesh in one process**. `decompose` and `reconstruct` need the most:
scotch holds the cell connectivity graph on top of the mesh, and
`reconstructPar` merges `np` processor meshes plus every written time.

MEASURED: 3000 MB / 30 min is ample at 128³ (2.10M cells) and fails outright at
203³ (8.37M cells). Hence 16 GB for meshing/init and 24 GB for
decompose/reconstruct. On a fat node with no `MaxMemPerCPU` this is cheap — a
large single-core request pins one core, not a node.

### 4.5 Cluster-specific submission requirements

Lichtenberg's `job_submit` plugin **requires** `--mem-per-cpu` on every job and
auto-routes the partition from the runtime. Supply it via `mem_mb_per_cpu` in
`default-resources`. Check your target cluster's equivalent — a required flag
missing from `default-resources` fails at submission, which at least is loud.

Do not pad `runtime`. An oversized limit can route a job to a slower partition
and holds the allocation long after a diverged solve has stopped producing
steps. Allow a per-study override for the solve only:

```python
SOLVE_RUNTIME = config.get("solve_runtime")   # minutes
# in rule solve's resources:
**({"runtime": int(SOLVE_RUNTIME)} if SOLVE_RUNTIME else {}),
```

Absent, every study submits exactly as before; present, a job is sized to its
measured pace.

---

## 5. `np` must have exactly one source of truth

An OpenFOAM-specific trap with no Snakemake analogue: the rank count lives in
**two** places — the launcher, and `system/decomposeParDict`'s
`numberOfSubdomains`. If they disagree, `decomposePar` writes one number of
`processor*` directories and the solver is launched with another, and the job
dies with a confusing mismatch.

Fix it by **generating `decomposeParDict` from `np` when materialising the
case**, so the committed dictionary is never read:

```python
def _write_decompose_par_dict(case_dir, np_):
    """Scotch: only numberOfSubdomains is needed -- it partitions the cell
    connectivity graph, so there are no per-direction coeffs to keep consistent
    and it works for any np and any (2D/3D) mesh."""
    text = ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class dictionary;\n    object decomposeParDict;\n}\n\n"
            f"numberOfSubdomains {np_};\n\nmethod          scotch;\n")
    ...
```

Use **scotch**, not `simple` or `hierarchical`: those need per-direction
coefficients whose product must equal `np`, which is one more thing to keep in
sync and which constrains what `np` may be. Scotch takes any `np` for any mesh.

---

## 6. A diverged run is a RESULT; a launch failure is not

Physics studies deliberately include arms expected to blow up. If a non-zero
solver exit aborts the DAG, one divergent arm leaves the whole study with **no**
aggregated table.

So: run the solver with output captured per case, and on non-zero exit record
the code and ensure the output file exists.

```python
_solve = ["cd {params.d}", f"rm -f {SOLVER}.csv .leia_solver_failed"]
# ... run solver, capture rc, write rc to .leia_solver_failed, touch the CSV
```

Nothing is swallowed: the exit code lands in a marker file, the tail goes to the
job log, the full log stays in the case directory, and the aggregator emits a
`solverFailed` column. An empty CSV is treated as absent. Run with
`--keep-going`.

The distinction that matters: **a diverged solve is data; a launch failure is
not.** Never let the second be recorded as the first — that is how a broken
launcher becomes a published convergence order.

### Classify logs with one shared helper, never ad-hoc greps

Every OpenFOAM solver prints in its startup banner:

```
trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
```

so any grep for a bare `Floating point exception` classifies **every healthy
run** as diverged. Likewise `core dumped` appears in unrelated SLURM kill text.
This false positive has bitten three times in this repository and once deleted a
live 16-hour run.

Write **one** classifier, use it everywhere, and match on tightened patterns
(`sigFpe::sigHandler`, genuine `... (core dumped)` signal lines, `FOAM FATAL`),
returning distinct states: COMPLETED / RUNNING / STALLED / DIVERGED /
LAUNCH_FAILURE / MISSING. Two corollaries:

* **Completeness** comes from OpenFOAM's terminating `^End$` (or the solver's
  per-step CSV), never from a post-processing output that only exists after a
  run finishes.
* **Liveness** comes from log mtime and step count, never from `squeue` alone —
  an empty `squeue` during a controller outage reads exactly like *"all jobs
  finished"*.

---

## 7. Orchestrating a long sweep on the cluster

Do **not** run `snakemake` directly on a login node for a multi-hour sweep.
MEASURED: a login-node sweep was SIGTERMed after ~5.5 h, mid-run — login nodes
reap long-lived user processes.

Submit a **one-core orchestrator job** that runs Snakemake, which then submits
the real work:

```bash
#!/bin/bash
#SBATCH -J proj-studies
#SBATCH -A <account>
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=4000
#SBATCH -t 7-00:00:00
set -o pipefail   # NOTE: no `set -u` -- OpenFOAM etc/bashrc reads unset vars

cd /path/to/repo || exit 1
module purge && module load <toolchain>
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
snakemake --workflow-profile profiles/slurm --configfile config/"$STUDY".yaml \
          --nolock --keep-going --rerun-triggers mtime
```

```bash
sbatch --export=ALL,STUDY=<study> run-studies.sbatch
```

Resuming is free: Snakemake skips cases whose outputs exist, so re-submitting
after any interruption continues where it stopped. `tmux` on a login node also
works but does not survive a node reboot; the batch job does.

**Record every job id you submit** (`sbatch --parsable ... >> .my_jobs`) and
cancel only from that list. On a shared account, `scancel -u $USER` destroys
other people's work — MEASURED here as two hours of cluster time lost, where
`sacct` showed `CANCELLED by 64+`, the truncation of the shared account's own
uid rather than an administrator. Identify your jobs by **name**
(`squeue -u $USER -n <jobname>`), never by user.

---

## 8. Porting checklist

1. `pip install snakemake-executor-plugin-slurm`.
2. Give every rule `resources: tasks=N`; only the solve gets `N > 1`.
3. Add `nodes: 1` to the solve rule (plugin resource, **not** `slurm_extra`).
4. Group serial pre/post rules; leave the solve **ungrouped**.
5. Route the launcher through `config["mpi_launcher"]` with a `{np}` placeholder.
6. Route the environment through `config["env_preamble"]`, applied by an
   `sh()` helper that wraps it in `set +eu` / `set -e`.
7. Generate `decomposeParDict` from `np` at materialisation; use `scotch`.
8. Write `profiles/local/config.yaml` with `resources: [tasks=<cores>]`.
9. Write `profiles/slurm/config.yaml` with your account, memory, runtime, and
   fat `set-resources` for the serial mesh/decompose/reconstruct rules.
10. Make a non-zero solver exit a recorded result, not a DAG abort; run with
    `--keep-going`.
11. Write one log classifier; forbid ad-hoc greps.

## 9. Verify before trusting any number

Run the **same small study both ways** and compare, in this order:

1. `--workflow-profile profiles/local --configfile config/<smoke>.yaml -n`
   (dry run: does the DAG build?)
2. The same locally for real at np=4.
3. The same on SLURM at np=4. Then, on the generated job:
   * `sacct -j <id> -o JobID,NCPUS,ReqCPUS,NNodes` — did the solve get `np` CPUs
     on **one** node?
   * **`ClockTime / ExecutionTime` from the solver log — is it ≈ 1 and not ≈ np?**
     This is the check that would have caught the 11 000× regression.
4. Compare the two results field by field. On identical `np` and identical
   decomposition they should agree to solver tolerance; if they do not, you have
   a real parallel defect (see below), not a profile problem.
5. Only then scale up.

**One caveat that outranks the profile work.** A serial pass proves nothing
about the code paths that exist only in parallel — collective calls, halo
exchange, coupled-patch values. Before any new or changed solver goes to the
cluster, run it locally on **at least four ranks** and read the result, not the
exit code: a deadlock leaves the job alive with a truncated log and no non-zero
return anywhere. MEASURED in this repository: four cluster arms sat alive and
silent for 76 minutes inside the first momentum assembly, because collective
reductions were called from inside a rank-0-only guard. Every serial test passed.

The profiles above make the cluster easy to reach. That is exactly why the
4-rank local gate has to come first.
