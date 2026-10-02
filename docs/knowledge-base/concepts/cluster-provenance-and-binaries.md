---
title: "Cluster provenance: binaries, ledger, ssh output"
description: "Every clone builds and runs its own binaries (etc/leia-env.sh, 2026-09-22) and every library prints its version stamp (2026-09-23); jobs are cancelled by id from the .my_jobs ledger, never by user (2026-08-28); a remote build is verified by its artefact, never by the exit code of ssh (2026-08-31)"
aliases: []
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-29
date_settled: 2026-09-24
decided_by: [author decision 2026-08-28, author decision 2026-08-31, docs/plan-library-split-and-build-policy.md]
code: [etc/leia-env.sh, etc/leia-stamp.sh, Allwmake, run-studies.sbatch, workflow/Snakefile, src/leiaLevelSet/leiaVersionRegistry.H, workflow/scripts/guard_finished_cases.py]
sources: [CLAUDE execution environment section, CLAUDE cluster and cancel sections, CLAUDE ssh section, CLAUDE provenance section, CLUSTER.md cancel and binaries sections, SLURM.md build and cancel sections, STATUS 9.3, STATUS 9.5, STATUS 10.1, STATUS 10.2, STATUS 10.4, STATUS 10.5, STATUS 11.12]
---
# Cluster provenance: binaries, ledger, ssh output

> Verdict (2026-09-29). A cluster result is about our code only after four checks. The checks are: which binary ran, which libraries it loaded, which jobs are ours, and whether the build happened. Every clone installs into its own `<clone>/platforms/`, because `etc/leia-env.sh` sets `WM_PROJECT_USER_DIR` to the clone root ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L38-L46)). Before that rule, a build from one clone landed in the shared account default. A second clone then found a library about 200 commits ahead of its own source on its path (2026-09-09, found 2026-09-22, [STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2773-L2780)). No job of that clone ran in that period ([STATUS 9.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2813-L2815)). Every leia library now carries a version stamp; every solver prints it in its banner and writes it to `<case>/leia.version` ([STATUS 10.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2980-L2994)). `scancel -u $USER` is forbidden on the shared account. One such cancel killed a 2x2x2 coupled matrix 19 minutes into its run, and five jobs of another session (2026-08-28, [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L213-L221)). Every job id goes into the ledger `.my_jobs`, and every cancel names ids from it. The ssh command can drop the output of a remote build and still exit 0. A rebuild once reported nothing, and the binary was still eight days old (2026-08-31, [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L374-L393)). The check is the artefact: the timestamp of the binary, and `strings` for a symbol that the change added.

## What it is

Before and after a launch:

1. In every shell that builds or runs a leia binary, source OpenFOAM's `etc/bashrc`, then `. ./etc/leia-env.sh` ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L38-L41)).
2. Build with `./Allwmake`. It prints the install directory, builds the nine leia libraries in dependency order, each with its stamp, and then relinks every application ([`Allwmake`](https://github.com/leia-openfoam/leia/blob/d1e3414/Allwmake#L3-L59)).
3. On a remote host, let the build script write its output into a file on that host. Read the file with a second ssh call ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L383-L388)).
4. Verify the artefact, not the exit code: `ls -la` the binary for a current timestamp, and `strings <binary or library> | grep -c <a symbol you added>` ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L390-L393), [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L501-L507)).
5. Submit from the clone root; `run-studies.sbatch` runs the clone that it was submitted from ([`run-studies.sbatch`](https://github.com/leia-openfoam/leia/blob/d1e3414/run-studies.sbatch#L26-L36)).
6. Record every job id in `.my_jobs`: `sbatch --parsable ... 2>/dev/null | tail -1 >> .my_jobs`. The site's submit plugin writes three lines before the id ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L115-L124)).
7. After the launch, read the `Exec` line of one solver log. It is the absolute path of the binary that ran ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L498-L503)).
8. Read the stamp lines of the banner. A line that differs from the others marks a library that was rebuilt from another tree ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L505-L511)).

To cancel:

1. Cancel by id from the ledger, or by the job name of this session (`scancel -n leia-curv`). Never cancel by user ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L223-L234)).
2. For a snakemake driver, read its `.err` for the lines `has been submitted with SLURM jobid N`. Cancel those child ids first, then the driver, and record every id ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L236-L241)).
3. Wait until `sacct -j <id>` shows CANCELLED or COMPLETED. Then rename or delete `studies/<study>`, and only then resubmit ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L134-L143)).

The mechanisms:

- `etc/leia-env.sh` sets `WM_PROJECT_USER_DIR`, `FOAM_USER_APPBIN` and `FOAM_USER_LIBBIN` to the clone. It removes every other OpenFOAM user directory from `PATH` and `LD_LIBRARY_PATH`, and it keeps OpenFOAM's own installation and ThirdParty ([`leia-env.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/etc/leia-env.sh#L20-L41)). The shared cfMesh install comes after the clone's directories, so it cannot hide a leia binary ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/etc/leia-env.sh#L43-L55)).
- The workflow's shell helper `sh()` sources the file in every job, after the profile's `env_preamble`. That preamble sources `etc/bashrc` again, which resets the user directory to the account default ([`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L127-L142)).
- `etc/leia-stamp.sh` writes `git describe --always` into a generated source file of each library. It adds `-dirty` when `src`, `applications`, `workflow`, `cases` or `config` has uncommitted changes, the same rule as the `gitCommit` column ([`leia-stamp.sh`](https://github.com/leia-openfoam/leia/blob/d1e3414/etc/leia-stamp.sh#L4-L19)). `aggregate.py` copies `leia.version` into the column `libStamps`, next to `gitCommit` ([STATUS 10.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2992-L2994)).
- A solver prints one stamp line per loaded leia library: at the split, two for `leiaSemiLagrangeLevelSetFoam` and seven for `leiaLevelSetFoam` ([STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3093-L3094)).

## Why it matters

- On 2026-09-01 fourteen arms of two 2D ladders ran an August 19 binary from the account default. They stopped at 0 to 1 steps with `Unknown fvOption type semiImplicitCapillaryForce` and were recorded as fourteen divergences. The PATH of the driver shell decided which binary `mpirun` launched ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L527-L532)).
- On 2026-09-09 the curvature clone's build wrote `libleiaLevelSet.so` into the account default at 11:00:40, 20 seconds after that clone pulled 0046961. Every SDPLS run from the SDPLS clone, at b3aa65e with object files of 2026-08-19, would then have loaded that library ([STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2773-L2780)). No SDPLS job ran between 2026-09-09 and 2026-09-22, so no result came from the mismatched pair ([STATUS 9.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2813-L2815)).
- A driver-side overlay did not reach the ranks. Job 54823248 resolved `leiaLevelSetFoam` to the account default, because the job preamble sourced `etc/bashrc` again (2026-09-22, [STATUS 9.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2795-L2800)).
- A stale solver against a rebuilt library segfaults at startup with an empty log, and that log reads as a divergence ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L815-L818), [SLURM.md](https://github.com/leia-openfoam/leia/blob/d1e3414/SLURM.md#L90-L92)).
- A study-global `LD_PRELOAD` for the mesher made the MPI solver segfault at startup (empty log, rc 139). The post-solve step recorded a "divergence"; the preload is now scoped to `pMesh` ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L1417-L1423), [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L819-L820)).
- A `--config` on the command line replaces the whole `config:` list of the profile. A baseline pass lost `mpi_launcher` and `env_preamble`, ran `mpirun` in one-task steps, and was voided (2026-09-23, [SLURM.md](https://github.com/leia-openfoam/leia/blob/d1e3414/SLURM.md#L270-L276), [STATUS 10.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3108-L3110)).
- A script with `set -e` stops without a message at the line that sources OpenFOAM's `etc/bashrc`. That file returns 1 on its last internal call ([STATUS](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2611-L2614)). The helper `sh()` therefore sources it under `set +eu` ([`Snakefile`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L142)).
- Between 10:59 and 11:12 on 2026-09-22 every file under `/work/scratch/tm83tomy` received a new mtime ([STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2781-L2784)). With `--rerun-triggers mtime`, 6 of 42 finished cases of `sdplsConv2Dvortex` and 2 of 14 of `sdplsExpSource2Dvortex` became due. The solve rule deletes a case CSV before it runs ([CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L392-L405)). Two measures followed: a `--touch` pass marked every study up to date on 2026-09-24, and the launch guard now runs before every launch ([[concepts/log-classifier-and-waiters]]).

## Where in the code

- `etc/leia-env.sh` ([L1-L55](https://github.com/leia-openfoam/leia/blob/d1e3414/etc/leia-env.sh#L1-L55)), created in [019ed08](https://github.com/leia-openfoam/leia/commit/019ed08) (2026-09-22); the cfMesh line in [b7e444c](https://github.com/leia-openfoam/leia/commit/b7e444c) (2026-09-23).
- `Allwmake`: removes a stale monolith `libleiaLevelSet.so`, runs `etc/leia-check-deps.py`, stamps and builds the nine libraries, then `wmake all applications` ([L9-L59](https://github.com/leia-openfoam/leia/blob/d1e3414/Allwmake#L9-L59)).
- `etc/leia-stamp.sh` ([L1-L28](https://github.com/leia-openfoam/leia/blob/d1e3414/etc/leia-stamp.sh#L1-L28)) and `src/leiaLevelSet/leiaVersionRegistry.H` ([L1-L13](https://github.com/leia-openfoam/leia/blob/d1e3414/src/leiaLevelSet/leiaVersionRegistry.H#L1-L13)).
- `workflow/Snakefile`, the helper `sh()` ([L127-L142](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile#L127-L142)).
- `run-studies.sbatch`, the orchestrator job of a study group ([L1-L45](https://github.com/leia-openfoam/leia/blob/d1e3414/run-studies.sbatch#L1-L45)).
- `workflow/scripts/guard_finished_cases.py`, the launch guard ([L1-L31](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/guard_finished_cases.py#L1-L31)).
- The site-independent recipe: [SLURM.md](https://github.com/leia-openfoam/leia/blob/d1e3414/SLURM.md#L84-L110) (build), [same](https://github.com/leia-openfoam/leia/blob/d1e3414/SLURM.md#L278-L306) (ledger and cancel).

## Evidence

| claim | number | where |
|---|---|---|
| the account-wide cancel | the 2x2x2 coupled matrix killed 19 minutes in; five `ded-*` jobs of another session; about two hours of cluster time; `CANCELLED by 64+` is uid 643395244, the account itself (2026-08-28) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L215-L221), [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L105-L113), MEASURED |
| a cancelled driver leaves its children | a 96-core solver child ran 53 minutes as an orphan under a UUID job name (2026-09-05) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L236-L241), MEASURED |
| a cancelled driver is not gone at once | it ran two more hours and re-queued dead arms as new 32-core jobs (2026-09-05) | [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L134-L143), MEASURED |
| ssh dropped the output of a rebuild | the binary was still eight days old (2026-08-31) | [CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L376-L381), MEASURED |
| a shared user directory | a library about 200 commits ahead of the clone's source; the `strings` output of it has 2 lines with `compositeFlux` (2026-09-09) | [STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2773-L2780), [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L533-L537), MEASURED |
| a stale binary recorded as divergences | 14 arms, 0 to 1 steps (2026-09-01) | [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L527-L532), MEASURED |
| the touch sweep | 10:59 to 11:12 on 2026-09-22; 6 of 42 and 2 of 14 finished cases due for a re-run | [STATUS 9.3](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2781-L2784), [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L394-L397), MEASURED |
| the touch pass | 43 studies in the SDPLS clone and 98 in the curvature clone, 0 dry-run failures (2026-09-24) | [STATUS 10.4](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3168-L3170), MEASURED |
| per-clone binaries change no number | WP1: 12 cases, both CSVs at tolerance 0; WP6, after the retirement: 42 CSV pairs per clone | [STATUS 10.1](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2935-L2941), [STATUS 10.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3223-L3226), MEASURED |
| the stamp shows an edited clone that was not rebuilt | an untracked file made `gitCommit` read `cc79df4-dirty` while the banner read `...-gcc79df4` | [STATUS 10.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3021-L3025), MEASURED |
| the pre-launch check in practice | the clone's solver on `PATH`, stamps `g4015404` without `-dirty`, the new symbols in the libraries, a dry run of 70 jobs (2026-09-26) | [STATUS 11.12](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3524-L3528), MEASURED |
| a cancel from the ledger | the `leia-gcls` gate cancelled by 209 child ids from its ledger (2026-09-26) | [STATUS 11.12](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3545-L3547), MEASURED |
| OpenFOAM-v2606 on Lichtenberg needs two build passes | pass 1: 151 applications, 113 libraries, `icoFoam` missing, job exit code 0 (FFTW built after `etc/bashrc` was sourced); pass 2: 270 and 129, "Critical systems ok" | MEASURED, [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4818-L4848), [CLUSTER.md](https://github.com/leia-openfoam/leia/blob/a177025d/CLUSTER.md#L64-L86) |
| The development line runs on Lichtenberg with v2606 | clone `leia-dev`: the CI script passes; the repair gate through `profiles/slurm` matches the laptop's v2512 run in every column except the four volume sums (2.7e-12 column-scaled) | MEASURED, [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4818-L4848) |

## Decisions

- 2026-08-28: never cancel by user; cancel by id from the ledger ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L213-L234)).
- 2026-08-31: verify the artefact, never the exit code ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L390-L393)).
- 2026-09-22: every clone builds and runs its own binaries, WP1 of the build-policy plan ([STATUS 10.1](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2908-L2929)).
- 2026-09-23: a version stamp in every library, WP2 ([STATUS 10.2](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L2980-L2994)).
- 2026-09-24: the old shared install folders retired with a dated suffix, WP6 ([STATUS 10.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3198-L3218)).

## Open questions

1. Delete the retired folders on or after 2026-10-01, if both sessions have run a study from the new folders by then ([STATUS 10.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3231-L3234)).
2. `CLUSTER.md` says that a banner lists eight stamp lines ([L508-L510](https://github.com/leia-openfoam/leia/blob/d1e3414/CLUSTER.md#L508-L510)). The count is the number of loaded leia libraries, and since 2026-09-26 there are nine libraries ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L47-L51)).

## Related

- Hub: [[hubs/verification]].
- Siblings: [[concepts/log-classifier-and-waiters]] (liveness, the launch guard, waiters that waited for text that ssh dropped), [[concepts/bit-identity-and-inertness-gates]] (the WP1 to WP6 gates), [[concepts/wrong-setup-voids]] (stale arms in a resubmitted sweep), [[concepts/data-archive-per-version]], [[concepts/seam-checks-and-decomposition-invariance]].
- Models: [[models/semi-implicit-capillary-force]] (the fvOption that the stale binary of 2026-09-01 did not know).

## Log

### 2026-09-29
Created.

### 2026-09-30
OpenFOAM-v2606 on Lichtenberg (two passes, FFTW) and the new clone `leia-dev` on the development line (STATUS 11.21).
