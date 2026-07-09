#!/usr/bin/env bash
#
# Build the self-contained reveal.js deck from its editable source.
#
#   doc/slides/index.template.html   (edit THIS: external figure + CDN refs)
#        │  build.sh  →  workflow/scripts/export_html.py
#        ▼
#   doc/slides/index.html            (vertical section stacks; the tracked deck)
#   doc/slides/index-linear.html     (flat, front-to-back reading; ONLY with --linear,
#                                     not tracked in git -- generated on demand)
#
# Each output is ONE self-contained file: local figures/*.png are inlined as base64
# and the CDN assets (reveal.js, theme, MathJax) are fetched and inlined so the deck
# opens offline. Requires only python3 (standard library) + a network connection ON
# THE FIRST BUILD to inline the CDN assets; without network the deck still builds but
# keeps CDN <link>/<script> tags (so it then needs the internet to render).
# Regenerating the FIGURES is a separate, heavier step -- see README.md.
#
# Run from anywhere:
#   bash doc/slides/build.sh            # index.html
#   bash doc/slides/build.sh --linear   # index.html + index-linear.html
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../doc/slides
repo="$(cd "$here/../.." && pwd)"                       # repository root

command -v python3 >/dev/null || { echo "error: python3 not found on PATH" >&2; exit 1; }

python3 "$repo/workflow/scripts/export_html.py" "$here" "$@"

# sanity: the reveal structure must be balanced (unbalanced <section> = broken deck)
python3 - "$here/index.html" <<'PY'
import sys
t = open(sys.argv[1], encoding="utf-8").read()
o, c = t.count("<section"), t.count("</section>")
assert o == c, f"unbalanced sections in {sys.argv[1]}: {o} <section> vs {c} </section>"
n = t.count('mjx-container')  # 0 in source; MathJax renders at load, not here
print(f"[ok] {sys.argv[1]} built: {o} balanced sections, {len(t)//1024} KiB")
PY
echo "[ok] open doc/slides/index.html in a browser (offline-capable)."
