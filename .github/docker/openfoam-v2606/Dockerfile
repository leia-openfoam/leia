# The CI image of leia: Ubuntu 24.04 (noble) with the openfoam.com package openfoam2606-dev
# (OpenFOAM-v2606, the leia standard since 2026-09-28, STATUS 10.6) and the tools the build and
# the smoke run need. Built and pushed by .github/workflows/ci-image.yml to
# ghcr.io/leia-openfoam/openfoam-v2606_ubuntu-noble; used by .github/workflows/build.yml.
# It replaces tmaric/openfoam-v2206_ubuntu-focal (OpenFOAM-v2206, retired 2026-09-30).
FROM ubuntu:24.04

LABEL org.opencontainers.image.source="https://github.com/leia-openfoam/leia"
LABEL org.opencontainers.image.description="OpenFOAM-v2606 (openfoam.com packages) on Ubuntu 24.04 for the leia CI"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

ARG DEBIAN_FRONTEND=noninteractive

# 1. The openfoam.com Debian repository (its own script adds the key and the source list).
# 2. openfoam2606-dev: the binaries, the headers and wmake; its dependencies bring the compiler,
#    flex and the OpenMPI development files.
# 3. openmpi-bin (mpirun for parallel runs), git (the library version stamps read
#    `git describe`), python3 (etc/leia-check-deps.py, etc/leia-check-build.py, the workflow
#    scripts), bc and procps for the case scripts.
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl gnupg \
 && curl -fsSL https://dl.openfoam.com/add-debian-repo.sh | bash \
 && apt-get update \
 && apt-get install -y --no-install-recommends openfoam2606-dev openmpi-bin git python3 bc procps make \
 && rm -rf /var/lib/apt/lists/*

# The OpenFOAM environment of this image: source it in every step, then etc/leia-env.sh.
ENV LEIA_FOAM_BASHRC=/usr/lib/openfoam/openfoam2606/etc/bashrc

# A check at build time: the environment loads and wmake is on PATH.
RUN bash -c '. "$LEIA_FOAM_BASHRC" && echo "WM_PROJECT_VERSION=$WM_PROJECT_VERSION WM_OPTIONS=$WM_OPTIONS" && command -v wmake && command -v blockMesh && command -v mpirun'
