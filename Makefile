# =============================================================================
# leia -- one-command reproduction of all results, figures and documents.
#
# Prerequisite (once per shell): source OpenFOAM and have the Python deps
# (snakemake, foamlib, vtk, numpy, matplotlib) + pdflatex on PATH:
#     source $HOME/OpenFOAM/OpenFOAM-v2506/etc/bashrc
#
# Targets:
#     make build      compile the leia library + solvers (Allwmake)
#     make studies    run the convergence studies (snakemake); each study's results
#                     agglomerate into ITS theme's docs/<theme>/<slug>-article/data/
#     make docs       (re)build the reveal decks + compile the articles from that data
#     make all        build + studies + docs
#     make clean      remove the (regenerable) built decks + article PDFs
#
# One data source per theme -- docs/<theme>/<slug>-article/data/{figures,tables} --
# feeds BOTH the standalone reveal deck and the Elsevier article, so slides and
# paper always reflect the newest results and each is independently shareable.
# =============================================================================
SHELL   := /bin/bash
PROFILE ?= profiles/local
TASKS   ?= 12
SNAKE    = PATH=$$HOME/.local/bin:$$PATH snakemake --workflow-profile $(PROFILE) \
             --nolock --keep-going --resources tasks=$(TASKS) --rerun-triggers mtime

# Publication studies (override on the command line, e.g. make studies SL_STUDIES="uncachedConv2Dvortex").
SL_STUDIES ?= uncachedConv2Dvortex uncachedConv3Dshear uncachedConv3Ddeformation \
              uncachedConv3DshearPoly uncachedConv3DdeformationPoly
VE_STUDIES ?= bulkVortex steadyVortex2D staticExtension

ART_SL := docs/semi-lagrangian-level-set/sl-level-set-article

.PHONY: all build studies studies-sl studies-ve docs decks articles article-sl clean help
.DEFAULT_GOAL := help

help:
	@echo "leia reproduction (source OpenFOAM first):"
	@echo "  make build     - compile the leia library + solvers"
	@echo "  make studies   - run the convergence studies -> docs/<theme>/.../data"
	@echo "  make docs      - build the reveal decks + compile the articles from data"
	@echo "  make all       - build + studies + docs"
	@echo "  make clean     - remove regenerable built decks + article PDFs"

build:
	./Allwmake

studies-sl:
	@for cfg in $(SL_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies-ve:
	@for cfg in $(VE_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies: studies-sl studies-ve

decks:
	bash docs/build-decks.sh

article-sl:
	@command -v latexmk >/dev/null && (cd $(ART_SL) && latexmk -pdf -interaction=nonstopmode -halt-on-error *.tex) \
	  || echo "[skip] latexmk not found; install a LaTeX toolchain to build the article PDF"
articles: article-sl

docs: decks articles

all: build studies docs

clean:
	rm -f docs/*/*-presentation/*.html docs/*/*-presentation/*-linear.html \
	      docs/*/*-article/*.pdf docs/*/*-article/*.aux docs/*/*-article/*.log \
	      docs/*/*-article/*.bbl docs/*/*-article/*.blg docs/*/*-article/*.out
