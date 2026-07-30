# leia — agent guide

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

Full, verified workflow in **[CLUSTER.md](CLUSTER.md)**. Essentials:
`ssh tm83tomy@lcluster5.hrz.tu-darmstadt.de` (passwordless; `~/bin/licht N` helper);
SLURM account `special00004`; jobs **must** set `--mem-per-cpu` (the slurm
profile does); OpenFOAM-v2512 source-built in `$HOME/OpenFOAM`; pull raw output
back with `make pull-runs` / `make pull-study STUDY=<name>`.
