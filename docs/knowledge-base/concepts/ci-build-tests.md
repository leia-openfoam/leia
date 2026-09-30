---
title: "Build Tests: the CI of leia on OpenFOAM-v2606"
description: "Since 2026-09-30 every push and pull request to main and development builds leia on OpenFOAM-v2606 in the image ghcr.io/leia-openfoam/openfoam-v2606_ubuntu-noble, checks every Make/files target and runs the one-way 2Dtranslation case with checks on its numbers; the same script runs on a workstation"
aliases: [CI, Build Tests, ci-build-and-smoke, CI image]
kind: concept
status: settled
part: verification
tags: [concept, part/verification]
date: 2026-09-30
date_settled: 2026-09-30
decided_by: [author request 2026-09-30 ("make the new image for Docker; we will need CI to work")]
code: [.github/workflows/build.yml, .github/workflows/ci-image.yml, .github/docker/openfoam-v2606/Dockerfile, .github/scripts/ci-build-and-smoke.sh, etc/leia-check-build.py, cases/2Dtranslation/Allrun.sh]
sources: [STATUS 11.21]
---
# Build Tests: the CI of leia on OpenFOAM-v2606

> Verdict (2026-09-30). `Build Tests` builds leia on OpenFOAM-v2606 on every push and pull request to `main` and `development`, in the image `ghcr.io/leia-openfoam/openfoam-v2606_ubuntu-noble` ([`build.yml`](https://github.com/leia-openfoam/leia/blob/f010ed1a/.github/workflows/build.yml), [`Dockerfile`](https://github.com/leia-openfoam/leia/blob/f010ed1a/.github/docker/openfoam-v2606/Dockerfile)). One script does the work and runs the same way on a workstation ([`ci-build-and-smoke.sh`](https://github.com/leia-openfoam/leia/blob/f010ed1a/.github/scripts/ci-build-and-smoke.sh)). The first run on `development` passed in 8 minutes ([STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4795-L4811)). It replaces the image `tmaric/openfoam-v2206_ubuntu-focal`, with which the job failed on `main` on 2026-09-29.

## What it is

The script does four things and stops at the first failure:

1. It sources the OpenFOAM `etc/bashrc` given as its argument (default: `$LEIA_FOAM_BASHRC` of the image) and runs `./Allwmake`.
2. It runs `etc/leia-check-build.py`, which checks that every `EXE` and `LIB` target of every `src/**/Make/files` and `applications/**/Make/files` exists in the clone's `platforms/` ([`leia-check-build.py`](https://github.com/leia-openfoam/leia/blob/f010ed1a/etc/leia-check-build.py)). `wmake all` does not always exit non-zero when one application fails, so the exit code of `./Allwmake` is not evidence.
3. It runs `cases/2Dtranslation/Allrun.sh` (serial, N = 128, CFL 0.5) and classifies the log with `foam_log_state.sh` ([[concepts/log-classifier-and-waiters]]).
4. It checks the numbers of that run: `E_GEOM_ALPHA_REL` = 2 at t = 0 (the start and end circles are disjoint, so the exact end reference is in use), below 3e-3 at T; `E_VOL_ALPHA_REL` below 1e-12 at t = 0 and 1e-3 at T; `E_BOUND_ALPHA` = 0 ([[retractions/reversed-2dtranslation]]).

On a workstation:

```bash
.github/scripts/ci-build-and-smoke.sh $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc
```

The image is Ubuntu 24.04 with the openfoam.com package `openfoam2606-dev`, `openmpi-bin`, git, python3, bc, procps and make. The workflow `CI image` builds and pushes it when the Dockerfile or the workflow changes. The package is private to the organisation; the build job pulls it with its own token.

## Why it matters

Several sessions push to `development`, each from its own worktree or clone. A shared change can break another method's build, and until 2026-09-30 nothing ran the build outside a session. The CI runs the build on the standard OpenFOAM version at every push and gives every session the same pass or fail.

## Where in the code

- [`build.yml`](https://github.com/leia-openfoam/leia/blob/f010ed1a/.github/workflows/build.yml): the triggers, the container and the one step.
- [`ci-image.yml`](https://github.com/leia-openfoam/leia/blob/f010ed1a/.github/workflows/ci-image.yml): the image build, push and check.
- [`ci-build-and-smoke.sh`](https://github.com/leia-openfoam/leia/blob/a177025d/.github/scripts/ci-build-and-smoke.sh): the four steps; `fail()` prints a GitHub error annotation with the tail of the failing log, and a pass prints a notice with the counts and the numbers.

## Evidence

| claim | number | where |
|---|---|---|
| The script passes on the laptop (OpenFOAM-v2606 source build) | 26 executables and 14 libraries present, 0 missing; 168 steps; `E_GEOM_ALPHA_REL` 2.000000 -> 1.6090e-03; `E_VOL_ALPHA_REL` 2.27e-14 -> 5.7984e-04 | MEASURED, [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4795-L4811) |
| The script passes on Lichtenberg (clone `leia-dev`, job 55179979) | the same counts and the same numbers to every printed digit | MEASURED, [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4818-L4848) |
| The same end values as the study | `kinematicTranslation2D`, N = 128, CFL 0.5: 1.6090e-03 and 5.7984e-04 | MEASURED, STATUS 11.19 item 5 |
| The job passes on GitHub | run 36706765044, `development` f010ed1a, 8 minutes, the build step 7.5 minutes | MEASURED, [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4795-L4811) |
| No official development image of v2606 | OpenCFD's newest `openfoam-dev` tag is 2512; 2606 exists only as `openfoam-run`, without compilers | MEASURED (Docker Hub, 2026-09-30), [STATUS 11.21](https://github.com/leia-openfoam/leia/blob/a177025d/STATUS.md#L4786-L4794) |

## Why it failed, or why we think so

Three failures on the way, all fixed in the script or the image:

1. Sourcing OpenFOAM's `etc/bashrc` passes the script's own arguments to it, and it sources every readable file among them. A path to the bashrc itself in `$1` recursed until bash segfaulted (exit 139). The script runs `set --` first.
2. `cases/2Dtranslation/Allrun.sh` has mode 644 in git (exit 126). The script runs `bash ./Allrun.sh`.
3. The first image check called `foamVersion`, which the openfoam2606 packages do not install (exit 127). The check now prints `WM_PROJECT_VERSION`.

## Decisions

- The image is our own, built by GitHub, because the laptop has no Docker and OpenCFD publishes no v2606 development image.
- The CI reads its results through annotations: the public check-runs API returns them without login, the job logs need admin rights.

## Open questions

1. A parallel smoke run (np 2 on the 4-core runner) and the seam check against the serial run: not in the CI yet.
2. The actions `checkout@v4`, `upload-artifact@v4` and the Docker actions target Node.js 20; GitHub runs them on Node.js 24 with a warning.

## Related

- [[hubs/verification]]
- [[concepts/cluster-provenance-and-binaries]]
- [[concepts/log-classifier-and-waiters]]
- [[retractions/reversed-2dtranslation]]

## Log

### 2026-09-30
Created with the CI image, the build job and the script; the Lichtenberg run of the same script added (STATUS 11.21).
