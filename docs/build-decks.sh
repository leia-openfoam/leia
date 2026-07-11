#!/usr/bin/env bash
#
# Build the self-contained reveal.js decks (one per theme) from their editable
# sources. Each theme owns a *-presentation/ folder with ONE *.template.html whose
# <img src> point at that theme's single data source
# (../<theme>-article/data/figures/...). This inlines those figures as base64 +
# the CDN assets (reveal.js, MathJax) so the built *.html opens offline and is
# shareable standalone. Built *.html are NOT tracked (regenerated in seconds).
#
#   docs/semi-lagrangian-level-set/sl-level-set-presentation/sl.template.html    -> sl.html
#   docs/velocity-extension/velocity-extension-presentation/index.template.html  -> index.html
#
# Regenerating the FIGURES/tables is a separate, heavier step (the Snakemake
# studies) -- see the repo Makefile / README.
#
#   bash docs/build-decks.sh            # build every theme deck
#   bash docs/build-decks.sh --linear   # + the flat *-linear.html variants
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # repository root
command -v python3 >/dev/null || { echo "error: python3 not found on PATH" >&2; exit 1; }

shopt -s nullglob
built=0
for tpl in "$repo"/docs/*/*-presentation/*.template.html; do
    dir="$(dirname "$tpl")"; name="$(basename "$tpl")"
    python3 "$repo/workflow/scripts/export_html.py" "$dir" "$name" "$@"
    out="$dir/${name/.template.html/.html}"
    python3 - "$out" <<'PY'
import sys
t = open(sys.argv[1], encoding="utf-8").read()
o, c = t.count("<section"), t.count("</section>")
assert o == c, f"unbalanced <section> in {sys.argv[1]}: {o} vs {c}"
print(f"[ok] {sys.argv[1]}: {o} balanced sections, {len(t)//1024} KiB")
PY
    built=$((built + 1))
done
[ "$built" -gt 0 ] || { echo "[warn] no docs/*/*-presentation/*.template.html found" >&2; exit 1; }
echo "[ok] built $built theme deck(s)."
