# Secondary data of the pre-print, software version `shared-method-config-2026-09-01-192-g1150e68`

This folder holds every number of the pre-print `gclsLevelSet.tex` in this version. The folder
name is the version stamp of the leia libraries (`libleiaCore` in each case's `leia.version`)
that produced the headline results, the fixed 2D method gate. Nothing here is typed by hand:
`figures/make_archive.py` writes it from the raw output, and
`figures/make_result_figures.py` and `workflow/scripts/make_gate_tables.py` make the figures
and the tables of the paper from it. `MANIFEST.csv` lists every file with its row count, size
and SHA-256.

## Content

| folder | what | code version of the runs (libleiaCore stamp) | raw output |
|---|---|---|---|
| `gate/` | the 2D method gate of the paper: `summaries/<candidate>/{summary,orders,vsBaseline,seam}.csv`, `verdict.txt`, `candidate.json`, `exact1D.csv`; `cases.csv` (every case: arm, N, ranks, state, steps, stamp); `histories_<arm>.csv` (reduced time histories, about 300 rows per case) | `...-192-g1150e68` (both corrections of section 5.7) | Lichtenberg, `/work/scratch/tm83tomy/leia/studies/methodGate2D_*` (orchestrator 55048916) |
| `prefix/` | the same gate before the corrections of section 5.7; `translating_endTime0p1.csv`: the 0.1 s translating runs | `...-173-g79b5a92` | Lichtenberg, `studies/methodGate2D_*_pre-20260927-*`, `methodGate2D_summary_pre-20260927-020856`, `methodGate2D_*_translating_endTime0p1_20260927` |
| `laptop/` | `translating_runs.csv`, `translating_histories.csv` (section 7.8), `seam_np4_vs_serial.csv`, `seam_face_density.csv` (section 5.7), `eulerian_mass_flux.csv` (Table 4) | per run in `translating_runs.csv`; the corrected runs used libraries stamped `...-183-g935fd4e-dirty` with solver binaries built from b1798c3, whose solver sources equal those of 1150e68 | laptop, `runs/gcls-laptop-20260927` (git-ignored), copied to Lichtenberg `/work/scratch/tm83tomy/leia/runs/gcls-laptop-20260927` |

## Regenerate

```bash
# on Lichtenberg, in the clone (gate and pre-fix record)
python3 docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_archive.py gate \
    --studies-dir studies --summary-dir studies/methodGate2D_summary --out <archive>/gate
python3 docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_archive.py prefix \
    --studies-dir studies --summary-dir studies/methodGate2D_summary_pre-20260927-020856 --out <archive>/prefix
# on the laptop (the discriminator runs)
python3 docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_archive.py laptop \
    --runs runs/gcls-laptop-20260927 --out <archive>/laptop
python3 docs/gradient-controlled-level-set/gcls-level-set-article/figures/make_archive.py manifest --out <archive>
# the figures and tables of the paper
cd docs/gradient-controlled-level-set/gcls-level-set-article
python3 figures/make_result_figures.py --archive <archive> --figures data/figures --tables data/tables
```

A new version of the pre-print gets a new folder named after the stamp of its data; an older
folder is never overwritten.
