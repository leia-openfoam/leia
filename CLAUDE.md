# leia — agent guide

> **CLAUDE.md and AGENTS.md are byte-identical duplicates of this guide.**
> Claude Code auto-loads `CLAUDE.md`; other agent tools auto-load `AGENTS.md`;
> neither one loads the other, and a markdown link is documentation, not a load
> instruction. Both files therefore carry the complete guide.
> **Any edit must be applied to both in the same commit** — `diff CLAUDE.md
> AGENTS.md` must print nothing. `CLUSTER.md` (cluster workflow) and
> `STATUS.md` (current state of the work) stay separate and are read on demand.

OpenFOAM level-set library + solvers, a Snakemake verification suite, and
thematic docs (reveal decks + Elsevier articles) fed by a single per-theme
`docs/<theme>/<slug>-article/data/`.
## Execution environment

- **All OpenFOAM / Python / Snakemake runs happen in WSL (Ubuntu), never Windows.**
  From a Windows shell, prefix commands with `wsl bash -lc '...'`. When passing a
  script to a remote host, pipe a file (`ssh host 'bash -s' < script.sh`) rather
  than a heredoc — nested Windows→WSL→ssh quoting expands `$VAR`/`$(...)` in the
  wrong shell.
- **Build:** `source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake`
  (OpenFOAM-v2512 is the current standard version, in `$HOME/OpenFOAM` on both
  the WSL laptop and Lichtenberg).

## Repo layout & git discipline

- Versioned: `src/`, `applications/`, `workflow/`, `cases/`, `config/`, and the
  small results `docs/**/data/{figures,tables,mechanism}` (decks/papers build
  from these without re-running studies).
- Git-ignored (never `git add`): `studies/`, `runs/`, `.snakemake/`, built deck
  `*.html` (keep `*.template.html`), article PDFs.
- Hub is GitHub `leia-openfoam/leia`. Code moves by git; raw simulation output
  moves by rsync. See **[CLUSTER.md](CLUSTER.md)**.

## Running studies

```bash
make studies PROFILE=profiles/local     # WSL, mpirun (smoke tests, small N)
make studies PROFILE=profiles/slurm     # Lichtenberg, one sbatch per case
snakemake --workflow-profile profiles/local --configfile config/<study>.yaml -n   # preview
```

One study = one `(case, mesh, mode)`; backend switches purely via `PROFILE`.
`workflow/README.md` documents every `config/*.yaml` study.

## Cluster (Lichtenberg, TU Darmstadt)

Full, verified workflow in **[CLUSTER.md](CLUSTER.md)**; the site-independent
recipe (any SLURM cluster, same traps) is **[SLURM.md](SLURM.md)**. Essentials:
`ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de` (passwordless; `~/bin/licht N` helper);
SLURM account `special00004`; jobs **must** set `--mem-per-cpu` (the slurm
profile does); OpenFOAM-v2512 source-built in `$HOME/OpenFOAM`; pull raw output
back with `make pull-runs` / `make pull-study STUDY=<name>`.

### Never cancel jobs you did not start

**`scancel -u $USER` is forbidden.** Several sessions share the `tm83tomy`
account, so an account-wide cancel destroys other people's running work, not
just your own. MEASURED 2026-08-28: the 2x2x2 coupled matrix was killed 19
minutes in by exactly this -- `sacct` recorded `CANCELLED by 64+`, which is the
truncation of uid 643395244 = `tm83tomy`, NOT an administrator; the same signal
took out five `ded-*` jobs of another session at 18:08. Two hours of cluster
time, lost to a cleanup that was aimed at somebody else's jobs.

**Keep a ledger and cancel from it.** Record every job id you submit -- the
studies append theirs to `.my_jobs` in the clone -- and cancel only those:

```bash
scancel <jobid> [<jobid> ...]        # explicit ids: always safe
scancel -n leia-curv                 # by JOB NAME, this session's own name
```

Corollaries: identify your work by **job name**, never by user (`squeue -u
$USER` lists every session's jobs on the shared account); a job you did not
submit is not yours to kill even when it looks stale; and if a driver must be
replaced, cancel it by id and resubmit rather than clearing the queue.

## Run it on 4 ranks before it leaves the laptop

**No new or changed algorithm is submitted to the cluster -- and above all not
as a parameter study -- until it has run on at least FOUR MPI RANKS locally and
the result has been looked at.** Serial success proves nothing about the code
paths that only exist in parallel: processor-patch values, halo exchange, and
above all COLLECTIVE calls.

MEASURED 2026-08-28. `semiImplicitCapillaryForce` was gated on bit-identity,
equilibrium behaviour and cost -- all in SERIAL, because the gate config is
`mode: serial, np: 1` -- and shipped. Its diagnostics block called `gAverage`,
`gMax` and `gSum`, every one of them COLLECTIVE, from inside a
`Pstream::master()` guard: rank 0 blocked in the reduction and the other ranks
never entered it. Four cluster arms sat ALIVE BUT SILENT for 76 minutes inside
the first momentum assembly before anyone noticed, and the failure is invisible
in every serial test that could have been run.

That is one instance of a class this repository has hit repeatedly, always the
same way -- code that is correct on one rank and wrong on several:
`setVelocity` writing FACE values into coupled patches, biasing `fvc::grad(U)`
by O(1) at every processor boundary; `interfaceExtension::updateFlux` assigning
the raw flux over the whole boundary field; the psi filter's `fvc::average`
inheriting `calculated` patch types so `L(psi)` was uncoupled across seams; a
narrow-band dilation looping internal faces only, making the filtered CELL SET
decomposition-dependent.

**The gate.** Before `sbatch`, on the laptop:

```bash
decomposePar -force            # numberOfSubdomains >= 4 -- a serial case
mpirun -np 4 <solver> -parallel  # decomposes to ONE domain and proves nothing
```

and then CHECK the result rather than the exit code: a deadlock leaves the job
alive with a truncated log and rc = 0 nowhere in sight, so read the step count
and the metrics CSV. Where the study is itself parallel, add the
serial-vs-np=4 equivalence check (`config/seamConsistency3D{serial,par4}.yaml`
is the committed pattern); a decomposition-dependent answer is a defect even
when both runs complete.

Specifically for anything that touches an `fvMatrix` -- an `fvOption`, a new
term, a new solver -- the parallel smoke test is part of the inertness gate,
not an extra.

## Time integration: BDF2 (`backward`) everywhere

**Momentum uses OpenFOAM's `backward` (BDF2) scheme in every leia two-phase
study — never Euler.** There is no "historical Euler default" to fall back to:
`MOMENTUM_DDT_SCHEME backward;` is the repository default and any new case
template must reference the token rather than hardcoding a scheme. The SL
foot-point trace is second-order in time and must not be fed a first-order
velocity; BDF2 evaluates fluxes and sources at t^{n+1}, matching the
interface-pipeline placement (force built from psi^{n+1}).

Measured basis: BDF2 vs Euler on matched windows moved the stationary-droplet
gain by +11.1/+2.9/-3.0 percent (sign-flipping, i.e. noise) with volume and
shape within 1.2 percent — BDF2 costs nothing and is formally right, so it is
mandatory, not optional.

## Reading OpenFOAM logs: use the classifier, never ad-hoc greps

Every OpenFOAM solver prints, in its STARTUP banner:

    trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).

so any grep for a bare `Floating point exception` (or bare `core dumped`, which
SLURM attaches to unrelated kill text) classifies EVERY healthy run as diverged.
This false positive has bitten three times in this repository -- the original
Snakefile divergence classifier, a dt-bisection probe, and a completion waiter
that declared three live studies finished; the same class of mistake once
deleted a live 16-hour interFoam arm.

**The rule: never write a new grep against a solver log.** Call

    workflow/scripts/foam_log_state.sh <log> [--stall S] [--wait [--poll S]]

which returns COMPLETED / RUNNING / STALLED / DIVERGED / LAUNCH_FAILURE /
MISSING (distinct exit codes, plus `steps=` and `age=`), using exactly the
Snakefile's tightened patterns (`sigFpe::sigHandler`, genuine
`... (core dumped)` signal lines, `FOAM FATAL`; launch failures separated so
they are never recorded as results). `--wait` replaces hand-rolled completion
loops. If the patterns ever need changing, change them in BOTH the helper and
workflow/Snakefile in the same commit.

Corollaries, learned the hard way:
- **Completeness** comes from OpenFOAM's terminating `^End$` (or the solver's
  per-step CSV), never from a post-processing output that only exists after a
  run finishes.
- **Liveness** comes from log mtime and step count (the helper's `age=`/
  `steps=`), never from `squeue` alone -- an empty `squeue` during a SLURM
  controller outage reads exactly like "all jobs finished".
- A DIVERGED run is a RESULT (a blow-up is data); a LAUNCH_FAILURE is not.

### Every waiter carries a timeout, and its exit condition must be reachable

Prefer `foam_log_state.sh --wait` over a hand-rolled loop. When a poll loop is
genuinely needed, three rules -- each written from an orphan this repository
actually produced:

- **`timeout` is mandatory.** An `until <cond>; do sleep N; done` with no
  timeout does not fail when its condition never arrives; it polls until the
  session ends. MEASURED 2026-08-31: THREE such waiters were still spinning
  hours after the work they watched had finished, and they were found only
  because someone asked whether anything was obsolete. An unbounded waiter is
  not a safety net, it is a leak.
- **Never grep for a pattern your own command line contains.** A waiter whose
  condition was `! pgrep -f "configfile config/<study>"` could never exit: the
  `pgrep` matched the shell running it, so the negation was permanently false.
  Same family as `pkill -f <pattern>` killing its own shell (exit 144), which
  this repository has hit repeatedly. Match on a file's CONTENT, or on a
  narrowly anchored process name, never on a string that appears in the
  invocation itself.
- **Ask whether the pattern can ever appear.** Two of the three orphans waited
  for text that the remote command never emitted, because ssh had swallowed it
  (below). A condition that is unreachable by construction is the same defect
  as an infinite loop, and it hides as "still running".

Corollary: identify the ORPHAN the way jobs are identified on the cluster --
from a ledger of what you started, by id. Never sweep by pattern or by user:
this laptop is shared with other sessions exactly as the `tm83tomy` account is,
and a `pkill`-style cleanup here takes out another session's container build or
job waiter. A background task from an earlier turn has a different parent PID
than the current shell, so "its parent is not my shell" does NOT mean "not
mine" -- that inference produced a wrong all-clear on 2026-08-31; query the
task list by id instead.

### ssh loses stdout: redirect on the remote side and read the file back

Long output from `ssh host 'bash -s' < script.sh` can vanish -- the command
exits 0 with nothing to show, so a build that never ran looks exactly like a
build that succeeded silently. MEASURED 2026-08-31: a cluster rebuild reported
nothing and was believed to have run; the binary was still eight days old, and
the study would have been submitted against stale code. It also cost two
orphaned waiters (above) that were waiting for text ssh had eaten.

Write the remote output to a file and read that file in a second call:

```bash
ssh host 'bash -s' < build.sh          # script redirects into /tmp/x.log
ssh host 'cat /tmp/x.log'
```

And verify the ARTEFACT, never the exit code: `ls -la` the binary for a current
timestamp and `strings <binary> | grep -c <a symbol you just added>` for the
change you expect. This is the same rule as reading a solver log rather than
`rc` -- exit 0 is not evidence that work happened.

## Method constraint: unstructured FVM only

**Never propose or pursue a method that depends on a structured or Cartesian
discretisation, on a named mode, or on a specific interface shape.** Everything
must work for the **unstructured Finite Volume Method as OpenFOAM implements it**
— or at minimum for unstructured methods that admit OpenFOAM meshes, use
**compact stencils**, and parallelise efficiently under **MPI domain
decomposition**.

Ruled out, non-exhaustively:

- filtering or projecting out named azimuthal/Fourier modes (e.g. "remove the
  m = 4, 8, 16 mesh-locked curvature content") — those modes exist only because
  the mesh is Cartesian and the interface is a circle or sphere;
- height-function constructions needing mesh-aligned columns;
- anything requiring a structured index space (i, j, k), a global transform, a
  wide or unbounded stencil, or a non-local solve that will not decompose cleanly
  across ranks;
- corrections tuned to one benchmark geometry.

**Why.** A mode- or mesh-specific correction can look excellent on the
stationary-droplet test and be worthless in general — and worse, it hides the
underlying defect behind a tuning parameter nobody can justify later. The target
is a method for arbitrary unstructured meshes and arbitrary interface topology; a
fix that does not survive that transition is not a fix.

**No partial solutions.** Do not bank a partial gain while the underlying defect
remains — e.g. do not raise the capillary time-step coefficient for a 2x speedup
when the growth it papers over is still there. That bakes a tuned number into
every case file for a problem that is about to be removed.

### No filtering in the production method

Filtering is a **research instrument only**. The psi filter (`biharmonicBand`),
curvature relaxation, and any smoothing of the level set, the curvature or the
force may be used to understand the defect; none of them may be what makes a
production run stable. The goal is a **Basilisk-like discretisation that is stable
on its own** — the stationary droplet relaxing to its numerical equilibrium with
the velocity falling to round-off because the operators are right.

A filter that stabilises hides the defect behind a coefficient nobody can justify,
and that coefficient then needs retuning per resolution. Measured: at matched
initial kick the theta = 0.2 band filter is 5.86x BETTER at R/h = 15.8 and 1.61x
WORSE at R/h = 10.0, turning damping into growth where there is no corrugation to
remove. That behaviour is a tuning knob, not a model.

Consequence: **score every candidate with all filters OFF**, and treat a filter's
benefit as a measurement of the defect's size rather than as a fix.

**How to apply.** When proposing a fix, state explicitly how it behaves on a
polyhedral unstructured mesh under MPI decomposition *before* proposing it.
Prefer formulations stated per cell or per face and summed — those are
mesh-agnostic by construction. Worked example that passes the test: the
variational capillary force, `f_c = -sigma * dA_h/dpsi_c` derived from a discrete
surface-energy functional `A_h = sum_c V_c |grad_h alpha_c|` — a sum over cells,
compact stencil, no structure assumed, no modes named, and zero net work around
any closed cycle by construction.

## Developing a numerical method for multiphase flow: the research loop

Six steps per experiment. They are cheap; the re-runs they prevent are not.
Every rule below is written from an incident in this repository's own history.

**1. Name the number, after decomposing the error.** State which measured
quantity the change must move, and its current value. Multiphase error metrics
almost always mix mechanisms, so decompose first — into *consistency* (the error
one step commits) and *stability* (how the scheme amplifies it), or into source
and amplifier. `max|U|(T) = u_0(h)*exp(G(h))` split a stalled question into a
converging kick and a growing amplifier, and showed that implicitness attacks
`G` while variational pairing attacks `u_0`. A proposal that cannot name its
target factor is not ready to run.

**2. Pre-register the read-out.** Prediction, pass/fail threshold, and what
result would *falsify* the idea — written in the config header BEFORE launching.
Then the data decides. Campaign-turning results came from gates written this way
(the 2x2 alignment probe, the varying-curvature ellipse gate that retracted an
adopted delivery); every expensive retraction came from a run launched without
one.

**3. Cheapest discriminator first.** The multiphase ladder is fixed:

    exact-solution unit gate    seconds  curvature of a sphere, alpha of a plane,
                                        closest point on an ellipsoid; the app
                                        exits nonzero on failure
    static field gate           minutes  no flow: estimator accuracy and its
                                        noise gain on a frozen interface
    kinematic transport         min-hrs  prescribed velocity, reversed flow:
                                        transport order with NO force in the loop
    well-balanced coupled gate  hours    coupled, but with EXACT or prescribed
                                        curvature: isolates the force balance
                                        from the estimator
    coupled 2D                  hours
    coupled 3D, then polyhedral core-hrs

Escalate only on a pass, and never spend a coupled run on a question a static
gate answers. **Anything new crosses the 4-rank gate before the cluster rung**
(see "Run it on 4 ranks before it leaves the laptop") -- a serial pass says
nothing about collectives, halo exchange or coupled patches. Two rungs are multiphase-specific and are the ones most often
skipped: **kinematics before dynamics** (a coupled failure is otherwise
ambiguous between transport and force), and the **well-balanced gate**
(prescribed exact curvature gave identically zero velocity in all six arms,
cleanly separating force-balance bugs from estimator error). Prefer a case with
a closed-form answer over one without.

**4. One variable, plus controls.** Audit every arm for a second difference
before launching — confounded verdicts have cost two full studies (a 3D verdict
confounded with a divergence scheme; a volume-conservation comparison confounded
with the psi time scheme). Multiphase arms are matched on mesh, time step AND
ITS LAW (`dt ~ h^p` hides an h-dependence inside a dt sweep), phase indicator,
and horizon. Put controls in the same matrix — a null arm, and where relevant a
re-run of the prior falsification — so the matrix validates itself.

**5. Check the measurement before believing it.** Is the perturbation mode
resolvable at this resolution (m <= 4 at N = 64)? Is `R/h >= 10`? Is the fit
window at least one capillary period? Are there enough ladder points — a fourth
point falsified a trend that three had made look like clean second order? Is
each metric read at the right instant (reversed flows: gradient at T/2, shape at
T, volume at both, because the reversal cancels errors at the endpoint)? Are the
compared runs at equal step counts (endpoint estimators are not comparable
across unequal horizons)? A number failing these is not a result.

**Is the interface still inside the domain?** Before any conclusion is drawn from
`t_blow`, compute where the interface IS at that step and how far it is from the
nearest boundary -- and check what that boundary actually IS, in
`constant/polyMesh/boundary` rather than in the field files. A blow-up whose step count
reproduces to better than a percent (9337/9331/9367/9324 on four independent occasions,
against the 5-38% scatter this campaign documents for genuine instabilities) is a
GEOMETRIC event, not an instability; the reproducibility is the tell. It can also
invert a refinement reading -- "finer blows first" was really "the finer mesh transports
the droplet more accurately, so it reaches the boundary sooner".

RETRACTED 2026-09-02, and kept here because the retraction is the lesson: this rule was
first written as "the leading edge reaches the OUTLET at t = 0.08". There was no outlet.
`translatingDroplet2D` had all four sides in one `walls` patch, so every measurement
behind that text came from a closed slip box (see "A wrong setup voids its data"). The
general point survives -- proximity to a boundary explains reproducible blow-ups -- but
the specific mechanism was wrong because nobody checked which patches the mesh actually
had. Check the mesh first, then the distance.

**6. Report the whole vector.** Never a headline metric alone. For an interface
method that vector is **shape/geometric error, volume conservation error, the
interface-profile diagnostic (`|grad psi|`, or boundedness of alpha), the
spurious-current level, and the pressure-jump error** — reported together. A
candidate improved volume order while failing the gradient; a single-metric view
would have called it a win. State where it ran, on which commit, and with which
binaries.

## Constraints that gate what may even be proposed

- **Unstructured FVM only** (see the section above), compact stencils,
  MPI-decomposable. A fix that does not survive a polyhedral mesh under domain
  decomposition is not a fix.
- **2D and 3D from one implementation.** A model is written once and must work
  in both; dimension-specific code paths are how a method acquires a defect that
  only one geometry ever sees.
- **Every new model is runtime-selectable and inert by default** (dictionary
  entry / case token, never a hardcoded switch), so it can be swept as a study
  axis and so every existing study is unaffected.
- **Regularisation is an instrument, not a fix** (see "No filtering in the
  production method"). Score every candidate with filters, smoothing and
  reinitialisation OFF.
- **No partial solutions.** Do not bank a tuned coefficient that papers over a
  defect that is still there.

## A wrong setup voids its data. Re-run, do not reuse

**When a case is found to be set up wrong, every number it produced is void --
including the ones that look unaffected, and including arms that were merely
"adjacent" to the defect.** Rename the study directory with a `_VOID_<reason>_<date>`
suffix so it can never be curated by accident, fix the setup, and re-run from
scratch. Aim at EXACTNESS AT THE COST OF CPU HOURS: cluster time is cheap next to a
conclusion built on a case that was not what anyone thought it was.

Do NOT reason about which metrics "should still be valid" under the defect. That
reasoning is exactly as unreliable as the setup was, it cannot be audited later, and
a partially-trusted table is worse than no table because nobody downstream knows
which rows to distrust.

MEASURED 2026-09-02. `cases/translatingDroplet2D` had all four sides of its
blockMeshDict in a single `walls` patch, so the mesh had NO `inlet` and NO `outlet`
-- while every field in `0.org` had always declared both. OpenFOAM errors when a MESH
patch is missing from a field, but SILENTLY IGNORES a field entry that matches no
mesh patch, so `inlet { type fixedValue; value uniform (0.05 0 0); }` was parsed and
discarded on every run this case ever did. The case was a CLOSED BOX with slip walls;
slip imposes `U.n = 0`, so on the x-normal faces the uniform stream the case is
initialised with is incompatible with its own boundaries and the pressure projection
annihilated it on step 1 (first-step local continuity error 1.0e-05 against 1.6e-12
by step 3). `maxMagUPrime = max|U - (U0,0,0)|` then measured the ANNIHILATED FREE
STREAM at ~2*U0 from the first step onward -- not a spurious current. An entire
`div(rhoPhi,U)` scheme comparison, a droplet-leaves-the-domain mechanism, and a
"scheme-independent kick" were read off that metric before anyone checked
`constant/polyMesh/boundary`.

The equal-density control then running on the same case was STOPPED rather than
allowed to finish, even though its isolation argument (`ddt(rho) = 0` and
`div(rhoPhi) = rho div(phi) = 0` at ratio 1) does not formally depend on the boundary
conditions. On the author's decision: a setup that was wrong in one way cannot be
assumed wrong in only that way.

**Corollary -- verify the mesh has the patches the fields talk about.** After
`blockMesh`, `constant/polyMesh/boundary` is the authority on what exists; the
boundaryField entries in `0.org` are a wish list. A one-line check belongs in any new
case's gate:

```bash
foamDictionary -entry entry0 -keywords constant/polyMesh/boundary   # what EXISTS
```

and any field patch name not in that list is silently doing nothing.

## Retraction is a first-class action

When data contradicts a previous claim, retract it explicitly and propagate the
retraction the same day: `STATUS.md`, the plan document, the slides, and any
curated table. The record has twice carried a retracted mechanism while the real
one sat only in a commit message — the reader then gets the wrong story with no
signal. Write "I was wrong about X, the measurement is Y" in plain terms: a
retraction that is easy to find is worth more than one that is well hedged.

Corollary: **never report a number you have not seen land.** If a run is
incomplete, say so and give the partial state; do not extrapolate a verdict from
arms that have not finished.

## Nothing is inert until measured

A new template token ships with its inert default **in the same commit** — a
valueless token once broke materialization of every study sharing that case. And
"this default changes nothing" is a claim to be TESTED, not asserted: a
`residualControl` block with `tolerance 0` is NOT equivalent to no block, because
its presence changes which corrector PIMPLE treats as final and therefore which
solver settings (`U`/`p_rgh` vs `*Final`) the last pass uses. The gate is a
bit-identity run against the pre-change state on a committed case, run before any
physics.

## Solver convergence in a cancellation-dominated problem

A balanced-force formulation is a **catastrophic-cancellation problem**: the net
force is a small difference of two large fluxes, and the parasitic current is
the fraction the pressure projection fails to absorb — measured here at
0.02-0.15%, i.e. 2.5e-3 m/s out of a 1.7 m/s capillary predictor at N = 128.
OpenFOAM normalises the linear-solver residual by
`sum(|A psi - A psibar| + |b - A psibar|)`, which at convergence is `2 sum|b|`:
THE LARGE QUANTITY. A relative tolerance is therefore safe only when

    tolerance  <<  (fraction of the source that carries the signal)

and that fraction must be MEASURED, not assumed. On uniform orthogonal hex the
margin is six orders and the exoneration is direct: a ~30x tighter pressure
solve moved max|U|, the velocity residual and the non-absorbable fraction by at
most 2.4e-4 RELATIVE over a whole run. **That result does not transfer to
non-orthogonal meshes** — the same measurement records that there the solver
residual CAN exceed the structural residual, and strict PCG once gained 18-29x.
Re-establish the margin before quoting solver-converged results on polyhedral or
strongly non-orthogonal meshes.

Two corollaries. Only the final corrector of the final outer iteration uses
`relTol 0`; with the interface re-advected inside the outer loop, an
under-converged intermediate velocity moves psi and kappa before any tight solve
happens. And a **diverging smoother destroys diagnosability**: a momentum
`smoothSolver` reaching 1e+98 at 1000 iterations makes a physical blow-up
indistinguishable from a linear-algebra one in the log (there it was a symptom —
max|U| had already grown 3.5x per step for six steps — but establishing that
took a timing analysis that a robust Krylov solver would have made unnecessary).

## Provenance and preservation

- Raw output that a re-run will overwrite is preserved FIRST (rename the study
  directory with a dated suffix); solvers truncate their metrics CSV on restart.
- Every curated number under `docs/**/data/` is regenerated by a committed
  script from committed configs. A number nobody can regenerate is a number
  nobody can defend.
- Destructive cleanup never runs while any driver may be alive, and liveness
  comes from log mtime and step count via `workflow/scripts/foam_log_state.sh` —
  never from a declared output, never from `squeue` alone.
- Rebuild the library **and relink every solver that links it**: a stale solver
  against a rebuilt library segfaults at startup with an EMPTY log, which reads
  exactly like a divergence. When sessions share `$FOAM_USER_LIBBIN`, pin
  binaries into a session-local directory for anything being measured.
- Environment workarounds are scoped to the tool that needs them: a
  study-global `LD_PRELOAD` for the mesher segfaulted the MPI solver.
