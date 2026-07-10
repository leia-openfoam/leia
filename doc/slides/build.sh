#!/usr/bin/env bash
#
# Build the self-contained reveal.js decks from their editable sources.
#
#   doc/slides/index.template.html   (edit THIS: the VELOCITY-EXTENSION deck)
#   doc/slides/sl.template.html      (edit THIS: the SEMI-LAGRANGIAN deck)
#        │  build.sh  →  workflow/scripts/export_html.py
#        ▼
#   doc/slides/index.html            (velocity-extension deck; tracked)
#   doc/slides/sl.html               (semi-Lagrangian deck; tracked)
#   doc/slides/{index,sl}-linear.html (flat reading order; ONLY with --linear, not tracked)
#
# Each output is ONE self-contained file: local figures/*.png are inlined as base64
# and the CDN assets (reveal.js, theme, MathJax) are fetched and inlined so the deck
# opens offline. Requires only python3 (standard library) + a network connection ON
# THE FIRST BUILD to inline the CDN assets; without network the deck still builds but
# keeps CDN <link>/<script> tags (so it then needs the internet to render).
# Regenerating the FIGURES is a separate, heavier step -- see README.md.
#
# Run from anywhere:
#   bash doc/slides/build.sh            # index.html + sl.html
#   bash doc/slides/build.sh --linear   # + the flat *-linear.html variants
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../doc/slides
repo="$(cd "$here/../.." && pwd)"                       # repository root

command -v python3 >/dev/null || { echo "error: python3 not found on PATH" >&2; exit 1; }

for tpl in index.template.html sl.template.html; do
    [ -f "$here/$tpl" ] || { echo "[skip] no $tpl"; continue; }
    python3 "$repo/workflow/scripts/export_html.py" "$here" "$tpl" "$@"
    out="$here/${tpl/.template.html/.html}"
    # sanity: the reveal structure must be balanced (unbalanced <section> = broken deck)
    python3 - "$out" <<'PY'
import sys
t = open(sys.argv[1], encoding="utf-8").read()
o, c = t.count("<section"), t.count("</section>")
assert o == c, f"unbalanced sections in {sys.argv[1]}: {o} <section> vs {c} </section>"
print(f"[ok] {sys.argv[1]} built: {o} balanced sections, {len(t)//1024} KiB")
PY
done
echo "[ok] open doc/slides/index.html (velocity extension) or sl.html (semi-Lagrangian)."
