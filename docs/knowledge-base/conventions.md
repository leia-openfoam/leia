---
title: "Conventions"
description: "How to write in the leia knowledge base, on one page"
kind: index
status: settled
part: all
tags: [index]
date: 2026-09-28
---
# Conventions

[[index]] <- back

This vault is the concise, cross-linked record of the leia methods. It is a curated layer over
`STATUS.md` (the lab notebook), `METHOD.md` (the best configuration), the plan documents and the
pre-prints. A number is written first there, or in a curated CSV, and quoted here with its link.

## Kinds and folders

| kind | folder | one note per |
|---|---|---|
| hub | `hubs/` | moving part (advection, viscosity, surface tension, mass flux, gradient control, verification), and one for the method lines |
| concept | `concepts/` | idea or mechanism |
| model | `models/` | runtime-selectable family; the members are the rows of its table |
| decision | `decisions/` | settled setting or question (mirrors `METHOD.md` section 8.1) |
| retraction | `retractions/` | retracted, voided or corrected claim |
| case | `cases/` | benchmark family |
| study | `studies/` | pre-print or deck theme, or a campaign |
| session | `sessions/` | handover: `current.md` is rewritten each sitting; the dated notes are frozen |

## Frontmatter

Every note has `title`, `description` (one sentence; search, popovers and the graph show it),
`kind`, `status` (`settled`, `open`, `retracted`, `voided`, `candidate`), `part` (`advection`,
`viscosity`, `surface-tension`, `mass-flux`, `gradient-control`, `verification`, `all`), `tags`
(the kind, and `part/<part>`) and `date`. Optional: `aliases`, `date_settled`, `decided_by`
(config paths, or `author decision YYYY-MM-DD`), `code` (repository paths), `sources` (short
locators such as `STATUS 11.15`, `METHOD 8.1 row SL_CLIP`, `SL article sec:surften`,
`SL deck #/3/2`). Lists are flow lists `[a, b]` or block lists; no nested maps. Write `title`
and `description` in double quotes: a colon, a hash or a leading special character breaks the
YAML parser of the site (found by the first build, 2026-09-28).

## Body

The verdict comes first: one paragraph, dated, every clause linked. Then the sections of the
template of the kind (`templates/`). `## Log` is the last section: append-only, one
`### YYYY-MM-DD` entry per change, newest at the bottom.

## Links

- Another note: `[[folder/slug]]`, always with the folder, for example `[[hubs/advection]]` or
  `[[decisions/sl-clip-and-value-bound-off|SL_CLIP]]`. Root notes: `[[index]]`, `[[decision-log]]`.
- An article section: a GitHub blob URL pinned to a commit, at the `\section` line, with the file
  path and the `\label` in backticks.
- A pre-print PDF: `https://leia-openfoam.github.io/leia/preprints/<name>.pdf`.
- A deck slide: `https://leia-openfoam.github.io/leia/decks/<name>.html#/h/v`.
- `STATUS.md`, `METHOD.md`, `CLAUDE.md`, a plan: a GitHub blob URL pinned to a commit with the
  line range (`#L120-L135`); the link text is the locator, for example `STATUS 11.15`. A heading
  anchor on the branch is allowed when the section is short.
- Code: a GitHub blob URL pinned to a commit, with `#Lstart-Lend`, the path in backticks.
- Data: the curated table or figure on the branch; an archive folder by its stamp.
- Literature: the DOI.

## Rules

1. A claim without a `status` does not exist here.
2. Every number carries a link to the place where it lives.
3. Write in ASD-STE100 Simplified Technical English.
4. Never rename or delete a note. Repurpose it and add an alias.
5. A decision, a retraction, a closed or opened question, a finished gate: write the note or the
   log line in the same commit as the result.
6. Keep a note to one screen where possible. Link to the pre-print, the deck and `STATUS.md` for
   the details.
7. No file named `log.*` and no folder named with one digit: the repository `.gitignore` ignores
   them silently.
8. `python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base` must pass before a commit.

## Log

### 2026-09-28
Created with the vault.
