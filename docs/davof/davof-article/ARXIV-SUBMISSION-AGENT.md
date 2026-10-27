# Creating the arXiv submission of the DAVOF preprint (instructions for an agent)

Written 2026-10-27 for an agent that operates a browser (or guides Tomislav step by step)
to submit version 1.0 of the preprint to arXiv. The agent does not change the manuscript;
if anything below does not match what it finds, it stops and reports instead of improvising.

## 0. Guardrails

- The submitter of record is Tomislav Marić (tomislav.maric@tu-darmstadt.de); the agent acts
  on his account only with his session open in the browser. Never create an account, never
  change account settings, never add co-authors.
- Upload only the archive named in section 1. Do not edit, rename or add files inside it. If the
  archive does not process, report the AutoTeX log; do not patch sources on the arXiv site.
- The final "Submit" click is Tomislav's unless he has said in writing, in the same session,
  that the agent may press it. Everything up to the preview is the agent's job.
- No other side effects: no emails, no posts, no changes to the repositories.
- Record every identifier arXiv shows (submission id `submit/NNNNNNN`, later the arXiv id
  `YYMM.NNNNN`) in the report, verbatim.

## 1. What goes in (verify before starting)

| item | value |
|---|---|
| source archive | `C:\Users\TomislavMaric\Proton Drive\tomi.maric\My files\Documents\research\projects\2026-DFG-AVOF\_build\davof-preprint-arxiv-v1.0.tar.gz` (also `~/OpenFOAM/repos/davof/leia/docs/davof/davof-article/_build/davof-preprint-arxiv.tar.gz` in WSL) |
| archive contents | 20 files at the top level: `ms.tex`, `ms.bbl`, five `fig_*.tex`, eight `table_*_pp.tex`, `tikz_styles.tex`, `plic_ellipsoid.pdf`, `onering_convergence.png`, `fig_ladders.png`, `davof_advection_history_sphereTranslation3D_hex_N64_Co0.2.png`; no subfolder, no `.bib`, no PDF of the paper |
| expected build | pdflatex, document class `elsarticle` (in arXiv's TeX Live), 45 pages, no undefined references |
| metadata source | `~/OpenFOAM/repos/davof/leia/docs/davof/davof-article/arxiv-metadata.md` (the texts below are copied from it; if the two differ, that file wins and the agent reports the difference) |
| code and data | tag `davof-v1.0` in https://github.com/leia-openfoam/leia and https://github.com/tmaric/TwoPhaseFlow; both tag URLs must answer HTTP 200 before submission |

Pre-flight checks the agent runs (shell, read-only):

```bash
tar tzf "<archive>" | wc -l                 # 21 lines: "./" plus 20 files
tar tzf "<archive>" | grep -c "/.*/"        # 0: nothing in a subfolder
curl -s -o /dev/null -w "%{http_code}\n" https://github.com/leia-openfoam/leia/tree/davof-v1.0   # 200
curl -s -o /dev/null -w "%{http_code}\n" https://github.com/tmaric/TwoPhaseFlow/tree/davof-v1.0  # 200
```

If any check fails: stop, report, do not submit.

## 2. Metadata to enter (verbatim)

**Title**

Dual Area-Volume-of-Fluid (DAVOF): co-advecting volume and face-area fractions for a closed-form, cell-local interface normal and position on unstructured meshes

**Authors**

Tomislav Marić (Mathematical Modeling and Analysis, Department of Mathematics, Technical University of Darmstadt). arXiv's author field takes "Tomislav Marić"; the affiliation goes into the affiliation field if the form offers one, otherwise it is omitted (it is on the PDF).

**Abstract** (1913 characters, within arXiv's 1920 limit; paste as plain text, no LaTeX macros):

Geometric Volume-of-Fluid (VOF) on unstructured polyhedral meshes reconstructs the interface orientation from the volume fractions of a stencil of neighbours; on irregular cells and in the presence of wisps it degrades, and the curvature built on it drives the parasitic currents of capillary flows. We transport, consistently to the volume fraction of every cell, the wetted fraction of every face area. From the area and volume fractions, closed-form equations recover the interface normal and position of a cell, for the first time in VOF: the normal from one sum over the cell's own faces, the divergence theorem on the liquid sub-volume, and the plane's position from a second, without neighbours or iteration. Both identities are exact for a plane on hexahedral, perturbed and polyhedral cells. On a sphere and an ellipsoid the normal has the order of the face fractions, the position is second order for every state and mesh, and a quadric through one ring of recovered normals gives a converging curvature where the reference models show none. Transported with the plane clipping, the ring normal, the same identity on the union of the cell's ring, stabilises the state and carries a sphere over 3.3 diameters bounded and conservative, first order in the face fractions as the clipping law predicts, at plicRDF's position and shape errors. A preliminary transport implementation with one paraboloid per cell, its frame from the normal, its curvature tensor from the ring and its position from the volume fraction, carries the sphere bounded and conservative over the resolution ladder with face fractions of order 2.8 per period with the exact normal and 2.1 with the ring quadric's normal over N = 32 to 128. The loop has the order of the normal it is fed; a transported normal of second order is the research goal. This document is the basis on which the further development of the method is documented.

**Comments**

45 pages, 9 figures, 18 tables. Code and data at tag davof-v1.0: https://github.com/leia-openfoam/leia and https://github.com/tmaric/TwoPhaseFlow

**Categories**

Primary: physics.flu-dyn (Fluid Dynamics). Cross-lists: math.NA (Numerical Analysis), physics.comp-ph (Computational Physics).

**Other fields**

- Report number: none. Journal reference: none. DOI: none. MSC class: optional, `65M08, 76T06` if the form asks. ACM class: none.
- License: CC BY 4.0 (`http://creativecommons.org/licenses/by/4.0/`). If Tomislav says in the session that a journal with a restrictive policy is planned, choose arXiv's non-exclusive licence to distribute instead.

## 3. The submission, step by step (arxiv.org, logged in as Tomislav)

1. Open https://arxiv.org/submit and choose "Start New Submission".
2. **Verify contact information**: name Tomislav Marić, e-mail tomislav.maric@tu-darmstadt.de, affiliation TU Darmstadt. Confirm; change nothing.
3. **Licence**: select CC BY 4.0 (section 2). Continue.
4. **Policy / originality** statements: read; tick the affirmations only if they are true as stated (sole author, original work, not under embargo). Continue.
5. **Upload files**: upload the archive of section 1 as one file. arXiv unpacks it. The file list it shows must be the 20 files of section 1 and nothing else. If arXiv flags `ms.bbl` as unused or asks for `.bib`, that is wrong: `ms.tex` has `\bibliography{literature}` and arXiv uses the supplied `ms.bbl`; continue.
6. **Process**: let AutoTeX run. Open the log. Required outcome: "pdflatex" build succeeded, no errors; warnings about overfull boxes are acceptable. Download or open the produced PDF and check:
   - 45 pages (44 or 46 is acceptable if the TeX Live version differs; more than that means a missing input, report it);
   - page 1 shows the title of section 2, "Tomislav Marić", the affiliation, and the footer "Tomislav Marić, preprint, 27 October 2026";
   - no "??" anywhere (`pdftotext file.pdf - | grep -c "??"` must print 0);
   - the four graphics are present (page with the PLIC ellipsoid renderings, the one-ring convergence plot, the history plot, the ladders figure);
   - the appendix "Implementation record of the transport studies" and the bibliography are present at the end.
   If the build fails because `elsarticle.cls` or `elsarticle-num-names.bst` is not found, report it; the fix (adding the class file to the archive) is made by the maintainer, not on the site.
7. **Metadata**: enter title, authors, abstract, comments and categories from section 2 exactly. Leave journal reference, DOI and report number empty.
8. **Preview**: compare the preview page with section 2 field by field (title, author, abstract first and last sentence, comments, categories, licence). Take a screenshot for the report.
9. **Submit**: Tomislav presses it, or the agent only if his written permission for this click is in the session. After submission, record the submission id and the announced date shown on the confirmation page.
10. **After the announcement** (arXiv sends the e-mail; the id appears on the "My Articles" page): record the arXiv id `YYMM.NNNNN` and version `v1` in the report.

## 4. What to report back

- Pre-flight: the four checks of section 1 with their outputs.
- The AutoTeX result: page count, warnings worth knowing, the "??" count.
- The preview screenshot and the field-by-field comparison with section 2 (any difference is listed, even a punctuation change arXiv made).
- The submission id, and later the arXiv id with its version.
- Anything the agent did not do, and why.

## 5. After the arXiv id exists (not the submission agent's job, listed so that nothing is forgotten)

- `arxiv-metadata.md` and `FINDINGS.md` get the id.
- The animation's title card and YouTube description get the id; the preprint's next version cites the animation.
- The research proposal replaces "arXiv preprint, in preparation, 2026" by the id (the proposal session does this).
- A version 2 of the preprint, when it comes, cites `arXiv:YYMM.NNNNN` as its own v1 and bumps the repository tag.
