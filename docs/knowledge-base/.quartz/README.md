# Building the knowledge base site

`docs/knowledge-base/` is an Obsidian vault. Quartz v5 (pinned to `QUARTZ_TAG`, never vendored)
builds it into a static site; GitHub Actions (`.github/workflows/knowledge-base.yml`) deploys it
to <https://leia-openfoam.github.io/leia/> on every push that touches `docs/**`, `STATUS.md`,
`METHOD.md` or `CLAUDE.md`, together with the reveal decks (`decks/`) and the pre-print PDFs
(`preprints/`).

- `check_kb.py <vault>`: the gate (frontmatter, wikilinks, log coverage, names). Run it before
  every commit: `make kb-graph` runs it and regenerates the 3D graph data.
- `build_graph.py <vault>`: writes `graph3d/graph.json` (git-ignored) for the three.js page.
- `quartz.config.yaml`: the site configuration, copied over Quartz's own at build time.
- `build.sh [--serve]`: the local build (`make kb`, `make kb-serve`). Needs Node >= 22: it uses
  the node on PATH, nvm's Node 22, a Node 22 tarball downloaded into `build/node`, or Docker,
  in that order.

Preview the 3D graph without Node: `make kb-graph`, then
`python3 -m http.server -d docs/knowledge-base 8000` and open
<http://localhost:8000/graph3d/graph.htm?local=1>. The page is `graph.htm`, not `index.html`: Quartz
strips the `.html` extension of a non-Markdown file, and GitHub Pages would not serve the result as a
directory index.
