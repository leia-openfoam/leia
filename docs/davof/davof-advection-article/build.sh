#!/usr/bin/env bash
# Build of the advection preprint (LaTeX only; nothing else computes on the laptop).
# The table bodies table_cmp_*.tex are written by make_table_cmp.py from data/swept/cmp_table.csv
# and the panels fig_sweptflux_{a,b}.tex by make_fig_sweptflux.py, both run on Lichtenberg;
# they are committed and never edited by hand.
set -e
cd "$(dirname "$0")"
mkdir -p _build
BIBINPUTS=.. latexmk -pdf -outdir=_build -interaction=nonstopmode ms.tex
pdfinfo _build/ms.pdf | grep Pages
