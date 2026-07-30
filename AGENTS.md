# Agent guide

This repo's agent instructions live in **[CLAUDE.md](CLAUDE.md)** (execution
environment, git discipline, running studies) and the cluster workflow in
**[CLUSTER.md](CLUSTER.md)** (Lichtenberg / TU Darmstadt: SSH, SLURM account
`special00004`, OpenFOAM-v2512, rsync helpers).

Read both before making changes or running simulations. Key rules:

- Run OpenFOAM/Python/Snakemake in **WSL**, not Windows (`wsl bash -lc '...'`).
- Code moves by **git** (hub: GitHub `leia-openfoam/leia`); raw simulation output
  (`studies/`, `runs/`) is git-ignored and moves by **rsync** — never `git add` it.
- Build against `$HOME/OpenFOAM/OpenFOAM-v2512`; run studies with
  `make studies PROFILE=profiles/{local,slurm}`.
