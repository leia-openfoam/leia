#!/bin/bash
# The arXiv package of the DAVOF preprint (2026-10-24, build flow of 2026-10-27): the sources
# the build needs, flat, with the bbl of the build in _build/, into _build/arxiv/ and
# _build/davof-preprint-arxiv.tar.gz; the package is test-compiled from scratch in a temporary
# folder with pdflatex + bibtex (what arXiv runs); the checklist is printed and written to
# _build/arxiv/CHECKLIST.md. Nothing is submitted here.
set -e
cd "$(dirname "$0")"
mkdir -p _build
python3 prepare_tables.py > _build/prepare_tables.log 2>&1
BIBINPUTS=.. latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=_build ms.tex > _build/latexmk.log 2>&1 || { echo "build failed, see _build/latexmk.log"; exit 1; }
OUT=_build/arxiv
rm -rf "$OUT"; mkdir -p "$OUT"
cp ms.tex _build/ms.bbl "$OUT"/
# every file ms.tex inputs or includes as a graphic (flat names)
for f in $(grep -o '\\input{[^}]*}' ms.tex | sed 's/\\input{//; s/}//'); do
  if [ -f "$f" ]; then cp "$f" "$OUT"/; elif [ -f "$f.tex" ]; then cp "$f.tex" "$OUT"/; else echo "MISSING input: $f"; fi
done
for f in $(grep -o '\\includegraphics\[[^]]*\]{[^}]*}' ms.tex | sed 's/.*\]{//; s/}//' | tr -d '{}'); do
  found=0
  for ext in "" .pdf .png .jpg; do [ -f "$f$ext" ] && cp "$f$ext" "$OUT"/ && found=1 && break; done
  [ $found -eq 1 ] || echo "MISSING graphic: $f"
done
# files the input files themselves input (the TikZ figures use tikz_styles.tex; table bodies are flat)
for f in tikz_styles.tex; do [ -f "$f" ] && cp "$f" "$OUT"/; done
# test compile from scratch, as arXiv does (pdflatex x3 + bibtex), in a temporary folder
TMP=$(mktemp -d)
cp "$OUT"/* "$TMP"/
( cd "$TMP" && pdflatex -interaction=nonstopmode ms.tex > c1.log 2>&1 && pdflatex -interaction=nonstopmode ms.tex > c2.log 2>&1 && pdflatex -interaction=nonstopmode ms.tex > c3.log 2>&1 ) \
  && TEST="compiles from the package alone: $(pdfinfo "$TMP/ms.pdf" | awk '/Pages/{print $2}') pages, $(grep -c 'undefined' "$TMP/ms.log" || true) undefined-reference warnings" \
  || TEST="DOES NOT COMPILE from the package alone, see $TMP"
ABS=$(python3 - <<'PY'
import re
t = open("ms.tex", encoding="utf-8").read()
a = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", t, re.S).group(1).strip()
a = a.replace("\\,\\%", " %").replace("\\%", "%").replace("\\", "")
print(len(a))
PY
)
PAGES=$(pdfinfo _build/ms.pdf | awk '/Pages/{print $2}')
PLACEHOLDERS=$(grep -c "XXXXXXXX\|YYYYYYYY" literature.bib _build/ms.bbl | awk -F: '{s+=$2} END{print s}')
cat > "$OUT/CHECKLIST.md" <<EOF
# arXiv submission checklist (DAVOF preprint), generated $(date -I)

- Sources: ms.tex, ms.bbl (bibliography pre-built; arXiv does not run bibtex on .bib), the input table bodies and TikZ files, the figures, flat in this folder. $TEST.
- Build: pdflatex + bibtex (elsarticle, preprint 12pt, in TeX Live on arXiv); $PAGES pages in the last build.
- Abstract in ms.tex: $ABS characters. arXiv's METADATA abstract is limited to 1920 characters: paste a shortened version there (the PDF keeps the full one).
- Title/author: sole author Tomislav Marić, Mathematical Modeling and Analysis, Department of Mathematics, TU Darmstadt; footer "Tomislav Marić, preprint, 27 October 2026".
- Primary category: physics.flu-dyn; cross-list: math.NA, physics.comp-ph.
- License: CC BY 4.0 (or arXiv's non-exclusive licence if a journal version is planned with a restrictive publisher).
- Placeholders still in the bibliography (must be real DOIs or the entries removed before submission): $PLACEHOLDERS lines matching XXXXXXXX/YYYYYYYY (the figshare data and code deposits).
- Code: TwoPhaseFlow fork https://github.com/tmaric/TwoPhaseFlow branch davof at the commit of the tables (7c06c57); leia https://github.com/leia-openfoam/leia branch development. Tag both (e.g. davof-preprint-v1) before submission so the citation is stable.
- NOT SUBMITTED: Tomislav submits. Before submission: his read of the PDF, the leia commit of docs/davof/davof-article (ms.tex, figures, scripts, reviews/), the push, the tags.
EOF
( cd _build && tar czf davof-preprint-arxiv.tar.gz arxiv )
echo "package: $OUT ($(ls "$OUT" | wc -l) files), _build/davof-preprint-arxiv.tar.gz; abstract $ABS characters; $PAGES pages"
cat "$OUT/CHECKLIST.md"
