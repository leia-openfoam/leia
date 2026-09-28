#!/usr/bin/env bash
# Build the knowledge base with Quartz v5 into build/kb/quartz/public: the same steps as
# .github/workflows/knowledge-base.yml.
#
#   bash docs/knowledge-base/.quartz/build.sh            # build
#   bash docs/knowledge-base/.quartz/build.sh --serve    # build, then serve on http://localhost:8080
#   KB_DECKS=0 ...                                       # skip the reveal decks
#   KB_PREPRINTS=1 ...                                   # also compile the pre-prints (latexmk)
#
# Quartz v5 needs Node >= 22. Order of preference: node on PATH, nvm's Node 22, a Node 22
# tarball downloaded once into build/node (no change to the system), Docker's node:22-slim.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
kb="$repo/docs/knowledge-base"; out="$repo/build/kb"; q="$out/quartz"; tag="${QUARTZ_TAG:-v5.0.0}"
serve=""; [ "${1:-}" = "--serve" ] && serve="--serve"

python3 "$kb/.quartz/check_kb.py" "$kb"
python3 "$kb/.quartz/build_graph.py" "$kb"

mkdir -p "$out"
if [ ! -d "$q/.git" ]; then
    git clone --depth 1 --branch "$tag" https://github.com/jackyzha0/quartz.git "$q"
fi
rm -rf "$q/content"; mkdir -p "$q/content"
rsync -a --exclude .obsidian --exclude .quartz --exclude templates "$kb"/ "$q/content/"
cp "$kb/.quartz/quartz.config.yaml" "$q/quartz.config.yaml"

mkdir -p "$q/content/record"
cp "$repo/STATUS.md" "$repo/METHOD.md" "$repo/CLAUDE.md" "$q/content/record/"
cp "$repo"/docs/plan-*.md "$repo/docs/capillary-level-set-research-roadmap.md" \
   "$repo/docs/IMPROVEMENTS.md" "$q/content/record/" 2>/dev/null || true

if [ "${KB_DECKS:-1}" = 1 ]; then
    bash "$repo/docs/build-decks.sh"
    mkdir -p "$q/content/decks"
    find "$repo/docs" -path '*-presentation/*.html' ! -name '*.template.html' -exec cp {} "$q/content/decks/" \;
fi
if [ "${KB_PREPRINTS:-0}" = 1 ]; then
    mkdir -p "$q/content/preprints"
    for tex in "$repo"/docs/*/*-article/*.tex "$repo"/docs/*/*-report/*.tex; do
        [ -f "$tex" ] || continue
        dir=$(dirname "$tex"); name=$(basename "$tex" .tex)
        if (cd "$dir" && latexmk -pdf -interaction=nonstopmode -halt-on-error "$name.tex" > "latexmk-$name.log" 2>&1); then
            cp "$dir/$name.pdf" "$q/content/preprints/"
        else
            echo "[kb] skip $tex (see $dir/latexmk-$name.log)"
        fi
    done
fi

node_ok() { command -v node >/dev/null 2>&1 && [ "$(node -p 'process.versions.node.split(".")[0]')" -ge 22 ]; }
if ! node_ok && [ -s "$HOME/.nvm/nvm.sh" ]; then
    . "$HOME/.nvm/nvm.sh"; nvm use 22 >/dev/null 2>&1 || nvm install 22
fi
if ! node_ok; then
    nd="$repo/build/node"
    if [ ! -x "$nd/bin/node" ]; then
        f=$(curl -fsSL https://nodejs.org/dist/latest-v22.x/SHASUMS256.txt | grep -o 'node-v22[0-9.]*-linux-x64\.tar\.xz' | head -1)
        [ -n "$f" ] || { echo "[kb] cannot find a Node 22 tarball name at nodejs.org" >&2; exit 1; }
        echo "[kb] downloading $f into $nd (no change to the system)"
        mkdir -p "$nd"
        curl -fsSL "https://nodejs.org/dist/latest-v22.x/$f" | tar -xJ -C "$nd" --strip-components=1
    fi
    export PATH="$nd/bin:$PATH"
fi
if node_ok; then
    (cd "$q" && npm ci && npx quartz plugin install && npx quartz plugin resolve \
        && npx quartz build -d content -o public $serve)
elif command -v docker >/dev/null 2>&1; then
    echo "[kb] no Node >= 22 on PATH: building in node:22-slim"
    docker run --rm -it -v "$q:/work" -w /work -p 8080:8080 -p 3001:3001 node:22-slim \
        bash -c "npm ci && npx quartz plugin install && npx quartz plugin resolve && npx quartz build -d content -o public $serve"
else
    echo "[kb] neither Node >= 22 nor Docker is available" >&2; exit 1
fi
echo "[kb] site: $q/public"
