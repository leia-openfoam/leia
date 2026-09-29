# Brief for writing notes of the leia knowledge base (2026-09-28)

You write Markdown notes into the Obsidian vault `docs/knowledge-base/` of the repository
`/home/tmaric/OpenFOAM/repos/leia-gcls` (branch `feature/gradient-controlled-level-set`, HEAD
`8867581`). The vault is published with Quartz to https://leia-openfoam.github.io/leia/ and is
read by humans and by future Claude sessions to find what is relevant, what was decided, what was
retracted, and why something failed or why we think so.

## Hard rules
0. YAML: write `title` and `description` ALWAYS in double quotes (a colon or a hash inside an
   unquoted value breaks the site build); escape an inner double quote as \\". Quote a list
   item that contains a colon. `check_kb.py` flags these.
1. Read `docs/knowledge-base/conventions.md` and the template of your kind in
   `docs/knowledge-base/templates/` first. Follow the frontmatter schema exactly (flow lists
   `[a, b]`; no nested maps). `## Log` is the last section of every note.
2. Write ONLY the files assigned to you in the manifest (`.quartz/manifest.md`). Do not edit any other
   file. Do not run git commands that change the repository (no add, commit, stash, checkout).
3. Write in ASD-STE100 Simplified Technical English: short sentences (max 25 words), active
   voice, one meaning per word, no idioms, no analogies, numbered procedures. Use only known
   CFD/FVM/multiphase terms; give a formula for a new quantity.
4. Every number carries a link to where it lives (a STATUS/METHOD line range pinned to commit
   8867581, an article section, a curated table, a config header). Never invent a number. If you
   cannot find a number in the sources, write "not recorded" rather than a guess. Never report an
   L_inf error as a result (L2 and L1 only); where the record only has L_inf, say so.
5. Link format: another note is `[[folder/slug]]` (always with the folder; the alias form
   `[[folder/slug|text]]` is allowed; inside a table cell escape the pipe: `[[folder/slug\|text]]`).
   Use only slugs that exist in the manifest (yours or another owner's). Root notes: `[[index]]`,
   `[[decision-log]]`, `[[retraction-log]]`, `[[conventions]]`.
   GitHub links: `https://github.com/leia-openfoam/leia/blob/<commit>/<path>#L<a>-L<b>` (pinned,
   drift-proof), with the path in backticks next to the link. `<commit>` is the commit whose tree
   you READ (`git rev-parse HEAD` in the clone, and the file must be unmodified there: check
   `git status --short <path>`); line numbers of a modified working-tree file do not match any
   commit. The notes of 2026-09-28 pin to d1e3414 (the writers) and 8867581 (the session lead). Article sections: the same blob
   link at the `\section` line plus the `\label` in backticks. Deck slides: the site URL
   `https://leia-openfoam.github.io/leia/decks/<name>.html#/<h>` (h = the horizontal slide index,
   counted from 0; add `/<v>` for a vertical slide). Pre-prints: `https://leia-openfoam.github.io/leia/preprints/<texname>.pdf`.
   Literature: the DOI.
6. Each note: half a page to two pages. Verdict first (one dated paragraph, every clause linked).
   Then the template sections: What it is / Why it matters / Where in the code / Evidence (a table
   claim | number | where) / Why it failed, or why we think so (only if something failed) /
   Decisions / Open questions / Related / Log. Empty sections may be dropped except Related and
   Log. `Related` lists the hub(s) of the part and the sibling notes.
7. `status`: settled (a measured, current result), open (unmeasured or contested), retracted,
   voided (the setup was wrong), candidate (not yet gated). `part`: advection | viscosity |
   surface-tension | mass-flux | gradient-control | verification | all. `tags`: `[<kind>, part/<part>]`.
   `sources`: short locators such as `STATUS 4.1`, `METHOD 8.1 row SL_CLIP`, `SL article sec:visc`.
   `code`: repository paths of the code the note is about.
8. Mark the reliability of each claim in the Evidence table's "where" column: a line reference is
   MEASURED; a derivation is DERIVED; a reading without a run is HYPOTHESIS.
9. Sources to read (all under the repository root): `STATUS.md` (the lab notebook, ~4000 lines;
   section 0 READ THIS FIRST, section 4 state of the measurements, section 11 gradient control),
   `METHOD.md` (the best configuration; section 8.1 is the decision table), `CLAUDE.md` (the
   rules), `docs/plan-*.md`, `docs/capillary-level-set-research-roadmap.md`, the articles under
   `docs/<theme>/<slug>-article/*.tex`, the decks `docs/<theme>/<slug>-presentation/*.template.html`,
   `cases/default.parameter` (the token defaults with their decision comments),
   `config/gates/methodGate2D.yaml`, `workflow/README.md`. The file `kb-raw-material.md` next to
   this brief is an extracted index of decisions, findings and line references (tags: S = STATUS,
   M = METHOD, C = CLAUDE, PCS = plan-curvature-stabilization, PSH = plan-shannon-parasitic-currents,
   PCT = plan-combined-source-terms, PHL = plan-halo-limited-gradient-control, RM = roadmap,
   DP = cases/default.parameter, G2 = methodGate2D.yaml). Use it to find the lines, then READ the
   lines: quote the number as the source states it.
10. When you finish, run `python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base`
    from the repository root. Fix every message that names one of YOUR files, except
    "link [[...]] does not resolve" for a slug that the manifest assigns to another owner (that
    note is being written in parallel). Report: the list of files you wrote, the check output for
    your files, and any number you could not source.
