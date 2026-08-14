# =============================================================================
# leia -- one-command reproduction of all results, figures and documents.
#
# Prerequisite (once per shell): source OpenFOAM and have the Python deps
# (snakemake, foamlib, vtk, numpy, matplotlib) + pdflatex on PATH:
#     source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc
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
# `--resources tasks=N` is a LOCAL cap on concurrent MPI ranks -- it stops a
# laptop oversubscribing its cores. On SLURM the scheduler already does that,
# and the cap would instead serialise the solves (an np=8 study would run one
# case at a time). Lift it for the cluster profile.
ifeq ($(PROFILE),profiles/slurm)
TASKS   ?= 10000
else
TASKS   ?= 12
endif

# Cluster sync (see CLUSTER.md). CLUSTER is an ssh alias; REMOTE_DIR is the leia
# checkout on the cluster's PARALLEL file system (/work/scratch, where sims run).
# Raw output (studies/, runs/) is git-ignored and travels by rsync; docs/**/data
# travels by git.
CLUSTER    ?= lichtenberg
REMOTE_DIR ?= /work/scratch/tm83tomy/leia
REMOTE     ?= $(CLUSTER):$(REMOTE_DIR)
SNAKE    = PATH=$$HOME/.local/bin:$$PATH snakemake --workflow-profile $(PROFILE) \
             --nolock --keep-going --resources tasks=$(TASKS) --rerun-triggers mtime

# Publication studies (override on the command line, e.g. make studies SL_STUDIES="uncachedConv2Dvortex").
SL_STUDIES ?= uncachedConv2Dvortex uncachedConv3Dshear uncachedConv3Ddeformation \
              uncachedConv3DshearPoly uncachedConv3DdeformationPoly
# LINEAR semi-Lagrangian method line (linearTaylor reconstruction): an
# INDEPENDENT method with its own studies, docs theme and paper
# (docs/linear-semi-lagrangian-level-set/).
SL_LINEAR_STUDIES ?= linearConv2Dvortex linearConv3Dshear linearConv3Ddeformation \
                     linearConv3DshearPoly linearConv3DdeformationPoly
# Two-phase (Navier-Stokes) studies: stationary-droplet parasitic currents, the
# standalone mean-curvature accuracy study (no flow solve), and the field-montage
# visualisation job (velocity-glyph-over-alpha + curvature across resolution).
DROPLET_STUDIES ?= stationaryDroplet2D curvatureDroplet2D dropletFields2D
VE_STUDIES ?= bulkVortex steadyVortex2D staticExtension
# Geometrically-redistanced level set (leiaRedistancedLevelSetFoam): static
# circle convergence + trigger ablation + advection studies (2D + 3D).
BENCH_STUDIES ?= benchVortexEulerT2 benchVortexEulerT8 benchVortexSLT2 benchVortexSLT8 benchVortexVET2 benchVortexVET8 \
                 benchVortexSLimproved benchVortexSLimprovedPerturbed benchVortexSLflux benchVortexGRLfrozen
GRL_STUDIES ?= redistanceStatic2D redistanceCircle2D vortexTriggerGRL vortexBoundsGRL bulkVortexGRL 3DshearGRL 3DdeformationGRL
# SDPLS source line: the 2D reversed-vortex arm matrix (noSource/R/beta x both
# admissible linearizations), the beta-target sweep that separates a residual
# from a wrong target, and the 3D shear/deformation companions.
SDPLS_STUDIES ?= sdplsStability benchVortexEulerT2 benchVortexEulerT8 \
                 sdplsBetaSweep sdplsConv2Dvortex sdplsOrderAblation \
                 sdplsConv3Dshear sdplsConv3Ddeformation \
                 sdplsConvMoll2Dvortex sdplsConvMoll3Dshear \
                 sdplsConvMoll3Ddeformation \
                 sdplsBand2Dvortex sdplsBand3Dshear sdplsBand3Ddeformation
# EVERY study whose psi transport is an FV div(phi,psi) -- i.e. solver is
# leiaRedistancedLevelSetFoam / leiaLevelSetTwoPhaseFoam, or leiaLevelSetFoam
# without ADVECTION=semiLagrangian. These share ONE discretization
# (div(phi,psi) Gauss linearUpwind grad(psi), grad(psi) cellLimited leastSquares 1,
# nDefCorr >= 3) and must be re-run together whenever it changes, or the
# comparison tables mix discretizations -- exactly the 2D/3D SDPLS confound.
# The semi-Lagrangian studies have no div(phi,psi) term and are unaffected.
EULER_STUDIES ?= $(SDPLS_STUDIES) \
                 benchVortexGRLfrozen benchVortexVET2 benchVortexVET8 \
                 bulkVortex bulkVortexHighRes bulkVortexGRL \
                 steadyVortex2D steadyVortex2DHighRes shearVortex \
                 phaseIndicatorConvergence \
                 vortexTriggerGRL vortexBoundsGRL \
                 vortexThresholdGRL vortexUnclippedGRL vortexSDPLSGRL \
                 3DshearGRL 3DdeformationGRL

ART_SL := docs/semi-lagrangian-level-set/sl-level-set-article
ART_LSL := docs/linear-semi-lagrangian-level-set/lsl-level-set-article
ART_GRL := docs/geometrically-redistanced-levelset/grl-level-set-article

.PHONY: all build studies studies-sl studies-sl-linear studies-droplet studies-ve studies-grl studies-sdpls studies-euler studies-one print-sdpls-studies print-euler-studies check-discretization docs decks articles article-sl article-lsl article-sdpls article-grl comparison sl-quadratic sl-linear curvature curvature-mode-gate pressure-workflow pressure-compatibility-gate pressure-nonorthogonal-sweep pressure-operator-pair-gate pressure-rauf-gate pressure-tolerance-gate pressure-solver-gate clean help pull-runs pull-study
.DEFAULT_GOAL := help

help:
	@echo "leia reproduction (source OpenFOAM first):"
	@echo "  make build     - compile the leia library + solvers"
	@echo "  make studies   - run the convergence studies -> docs/<theme>/.../data"
	@echo "  make docs      - build the reveal decks + compile the articles from data"
	@echo "  make comparison- run the 6 benchVortex studies + build the comparison deck + paper"
	@echo "                   (= snakemake -s workflow/Snakefile.comparison --workflow-profile profiles/local)"
	@echo "  make sl-quadratic - QUADRATIC SL method line end-to-end: studies + figures + deck + paper"
	@echo "                   (= snakemake -s workflow/Snakefile.sl-quadratic --workflow-profile profiles/local)"
	@echo "  make sl-linear - LINEAR SL method line end-to-end: studies + figures + deck + paper"
	@echo "                   (= snakemake -s workflow/Snakefile.sl-linear --workflow-profile profiles/local)"
	@echo "  make curvature - curvature-only convergence study (static circle, no flow) -> table + figure"
	@echo "                   (= snakemake -s workflow/Snakefile.curvature --workflow-profile profiles/local)"
	@echo "  make curvature-mode-gate - perturbed-circle m=2,3,4 transfer gate on uniform/perturbed N=32,64,128 meshes"
	@echo "  make pressure-workflow - Snakemake DAG: nonorthogonal + compatibility + operator + rAUf + pressure algebra"
	@echo "  make pressure-compatibility-gate - alias for the corresponding pressure Snakemake rule"
	@echo "  make pressure-nonorthogonal-sweep - alias for the corresponding pressure Snakemake rule"
	@echo "  make pressure-operator-pair-gate - alias for the corrected/uncorrected pressure Snakemake rule"
	@echo "  make pressure-rauf-gate - alias for the paired pressure-operator/constant-rAUf Snakemake rule"
	@echo "  make pressure-tolerance-gate - alias for the pressure-algebra Snakemake rule"
	@echo "  make pressure-solver-gate - alias for the GAMG/PCG pressure Snakemake rule"
	@echo "  make studies-sdpls - SDPLS source line (2D arm matrix + 3D shear/deformation)"
	@echo "  make studies-euler - EVERY FV div(phi,psi) study; re-run together when the"
	@echo "                   discretization changes.  PROFILE=profiles/slurm on Lichtenberg."
	@echo "  make check-discretization - assert all FV studies ran ONE discretization"
	@echo "  make all       - build + studies + docs"
	@echo "  make clean     - remove regenerable built decks + article PDFs"

build:
	./Allwmake

studies-sl:
	@for cfg in $(SL_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies-sl-linear:
	@for cfg in $(SL_LINEAR_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies-droplet:
	@for cfg in $(DROPLET_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies-ve:
	@for cfg in $(VE_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies-grl:
	@for cfg in $(GRL_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
# SDPLS source line: the 2D reversed-vortex arm matrix (both linearizations x
# {noSource,R,beta}) plus the 3D shear/deformation companions.
# Enumerate a study group (scripts and humans: `make -s print-euler-studies`).
# Gate: every FV div(phi,psi) case ran ONE discretization. This is the check
# that would have caught the 2D-vs-3D confound. Exits non-zero on a mismatch.
check-discretization:
	@python3 workflow/scripts/check_discretization.py $(STUDY)

print-sdpls-studies:
	@echo $(SDPLS_STUDIES)
print-euler-studies:
	@echo $(EULER_STUDIES)

# Run ONE study by name: make studies-one STUDY=sdplsConv3Dshear
studies-one:
	@test -n "$(STUDY)" || { echo "usage: make studies-one STUDY=<name>"; exit 1; }
	$(SNAKE) --configfile config/$(STUDY).yaml

studies-sdpls:
	@for cfg in $(SDPLS_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
# EVERY study whose psi transport is a finite-volume div(phi,psi). These share
# one discretization (linearUpwind grad(psi), cell-limited grad(psi), nDefCorr
# >= 3), so they must be re-run TOGETHER whenever it changes -- otherwise rows
# from different discretizations end up in the same comparison table, which is
# exactly how the 2D/3D SDPLS confound happened. The 82 semi-Lagrangian studies
# have no div(phi,psi) term and are unaffected.
studies-euler:
	@for cfg in $(EULER_STUDIES); do echo ">>> $$cfg"; $(SNAKE) --configfile config/$$cfg.yaml; done
studies: studies-sl studies-sl-linear studies-droplet studies-ve studies-grl

decks:
	bash docs/build-decks.sh

article-sl:
	@command -v latexmk >/dev/null && (cd $(ART_SL) && latexmk -pdf -interaction=nonstopmode -halt-on-error *.tex) \
	  || echo "[skip] latexmk not found; install a LaTeX toolchain to build the article PDF"
article-sdpls:
	cd docs/sdpls-level-set/sdpls-article && latexmk -pdf sdplsLevelSet.tex

# One-command method comparison: runs the six benchVortex studies, harvests the
# figures + decision table, rebuilds the comparison deck AND compiles the paper.
# (Convenience alias for the documented snakemake entry point.)
comparison:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.comparison \
	  --workflow-profile $(PROFILE) --nolock --keep-going --resources tasks=$(TASKS)

# One-command METHOD LINES: each runs its studies, harvests all figures/tables,
# dual-copies into article + presentation data/, rebuilds deck + paper.
sl-quadratic:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.sl-quadratic \
	  --workflow-profile $(PROFILE) --nolock --keep-going --resources tasks=$(TASKS)
sl-linear:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.sl-linear \
	  --workflow-profile $(PROFILE) --nolock --keep-going --resources tasks=$(TASKS)

# One-command CURVATURE convergence study (static circle, no flow): runs
# config/curvatureDroplet2D.yaml across N=32..512 and harvests the model-comparison
# table + log-log figure. See workflow/Snakefile.curvature. The study is a few seconds
# per resolution, so it is recomputed from scratch every time (removes the gitignored
# studies/curvatureDroplet2D) -- guaranteeing the result always reflects the current
# `make build` binary, without any dependence on stale snakemake mtimes.
curvature:
	rm -rf studies/curvatureDroplet2D
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.curvature \
	  --workflow-profile $(PROFILE) --nolock --keep-going --forceall --resources tasks=$(TASKS)

curvature-mode-gate:
	python3 workflow/scripts/run_curvature_mode_transfer_gate.py --fresh

pressure-workflow:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-compatibility-gate:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility pressure_compatibility_gate \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-nonorthogonal-sweep:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility nonorthogonal_sweep \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-operator-pair-gate:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility pressure_operator_pair_gate \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-rauf-gate:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility pressure_rauf_gate \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-tolerance-gate:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility pressure_tolerance_gate \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

pressure-solver-gate:
	PATH=$$HOME/.local/bin:$$PATH snakemake -s workflow/Snakefile.pressure-compatibility pressure_solver_gate \
	  --workflow-profile $(PROFILE) --nolock --resources tasks=$(TASKS)

article-lsl:
	@command -v latexmk >/dev/null && (cd $(ART_LSL) && latexmk -pdf -interaction=nonstopmode -halt-on-error *.tex) \
	  || echo "[skip] latexmk not found; install a LaTeX toolchain to build the article PDF"

article-grl:
	@command -v latexmk >/dev/null && (cd $(ART_GRL) && latexmk -pdf -interaction=nonstopmode -halt-on-error *.tex) \
	  || echo "[skip] latexmk not found; install a LaTeX toolchain to build the article PDF"
articles: article-sl article-lsl article-grl

docs: decks articles

all: build studies docs

clean:
	find docs -path '*-presentation/*.html' ! -name '*.template.html' -delete
	rm -f docs/*/*-article/*.pdf docs/*/*-article/*.aux docs/*/*-article/*.log \
	      docs/*/*-article/*.bbl docs/*/*-article/*.blg docs/*/*-article/*.out

# --- cluster sync (git-ignored raw output only; see CLUSTER.md) --------------
# Pull ALL study output from the cluster to inspect/visualise locally.
pull-runs:
	rsync -avz --progress $(REMOTE)/studies/ studies/
# Pull ONE study: make pull-study STUDY=bulkVortexSL
pull-study:
	@test -n "$(STUDY)" || { echo "usage: make pull-study STUDY=<name>"; exit 1; }
	rsync -avz --progress $(REMOTE)/studies/$(STUDY)/ studies/$(STUDY)/
