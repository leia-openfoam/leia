<!-- Completed 2026-09-29. The record of what was done is STATUS.md section 11.18 and
docs/knowledge-base/sessions/current.md. Kept verbatim as approved on 2026-09-28. -->

# Plan: a linked knowledge base for the leia methods, the SL/source-term separation, and a critical technical report on gradient control

(The previous plan of this file, the halo-limited gradient-control implementation, is committed as
`docs/plan-halo-limited-gradient-control.md` and is complete through its first 2D gate campaign.)

## Context

**Why.** leia now has six method lines (quadratic and linear semi-Lagrangian, SDPLS, normal-projected SL,
geometric redistancing, gradient control) with one pre-print and one deck each, plus STATUS.md (4000
lines), METHOD.md, five plan documents and a research roadmap. The reasoning behind decisions, the
retractions, and the "why it failed" of every negative result are spread over these files. A new session
cannot find what is relevant without reading all of them. The user wants (1) a concise, human-readable,
cross-linked knowledge base in `docs/` that tracks the reasoning and decisions of method development,
built from the slides and pre-prints, that Obsidian opens directly and that GitHub Pages publishes
(Quartz, one source, no duplicate); (2) the semi-Lagrangian topic separated from the source-term topic:
the SL findings of the gradient-control session move out of the gcls pre-print into the SL article and a
handover note for the separate SL session; (3) a critical technical report (LaTeX, PDF) on the source-term
and velocity-extension approach: why each candidate failed and which methods to proceed with, summary
first; (4) CLAUDE.md/AGENTS.md rules: the knowledge base is the point of reference, and the working mode
"user = supervisor, Claude = expert numerical-methods developer; propose, assess critically, decide
together" for SL and every other leia method line.

**Decisions taken by the user (2026-09-28).** Quartz (single Markdown source in `docs/knowledge-base/`,
published by a GitHub Actions workflow to the leia repository's GitHub Pages); SL findings move to the SL
article; the technical report is LaTeX with a PDF.

**Facts that shape the design (measured this session).**
- Quartz is at v5.0.0 (2026-03-14; v4.5.2 is the last v4). v5 configures with `quartz.config.yaml`,
  installs its components as plugins (`npx quartz plugin install --from-config`), needs Node >= 22 (the
  laptop has Node 18.19.1, no nvm), builds `content/` to `public/`, and its documented GitHub Pages
  workflow uses Node 24, `fetch-depth: 0`, `npm ci`, `npx quartz plugin install && npx quartz build`,
  `upload-pages-artifact` of `public`, `deploy-pages`. It supports wikilinks with aliases and headings,
  transclusion, callouts, Obsidian frontmatter (title, aliases, tags, draft, date, description), block
  references, Mermaid, KaTeX math, full-text search, backlinks, and a global and a local force-directed
  graph view (node size = link count; `showTags`, `removeTags`, `depth`, `enableRadial`). `baseUrl` takes
  a subpath for repository sites (`user.github.io/repo`).
- The repository already has `.github/workflows/docs.yml`, a Pages deployment of Doxygen from
  `doc/Doxygen`, a path that no longer exists (the API docs moved to `docs/api`); it runs on `main` only.
  One repository has ONE Pages site, so the knowledge-base workflow replaces this workflow and can host the
  API docs under `/api/` later.
- The user's own Obsidian vault (`repos/orga`) uses YAML frontmatter, `[[wikilinks]]`, `status:` fields,
  hub notes with `tags: [area]`, and an append-only `## Log` with dated `### YYYY-MM-DD` entries. The
  knowledge base follows these conventions so the two vaults feel the same.
- The velocity-extension article is a stub (four TODO sections); the SDPLS article holds the earlier
  coupled-source result (the strain-cancelling source amplifies the mesh-locked mode-4 current 260x; with
  exact curvature the loop has no seed); the SL article has no source-term content and no mass-flux
  section; the gcls pre-print (written 2026-09-27) holds SL-baseline content that must move.

## 1. The knowledge base (`docs/knowledge-base/`, an Obsidian vault, published with Quartz)

### 1a. Layout, kinds, templates, links, logs (design settled 2026-09-28; verbatim artefacts in Appendix D)

**Layout: kind folders, `part` as frontmatter and nested tags** (not one folder per moving part: many
notes belong to two parts, tags allow several, and both graphs colour or filter by tag; the user's own
vault sorts by kind and gives topics a hub). The vault root is the Quartz content root:

```
docs/knowledge-base/
  index.md            front door: the six hubs, the two logs, sessions/current, the 3D graph
  conventions.md      how to write here (the rules of section 4, condensed)
  decision-log.md     chronological, one line per settled decision, month headings ## YYYY-MM
  retraction-log.md   chronological, one line per retraction / void / correction
  hubs/               6 notes, one per moving part (+ verification)
  concepts/           the ideas and mechanisms (~13)
  models/             one note per runtime-selectable FAMILY, members as rows (~14)
  decisions/          one note per settled setting or question, mirrors METHOD §8.1 row groups (~16)
  retractions/        one note per retracted, voided or corrected claim (~12)
  cases/              one note per benchmark family (~8)
  studies/            one note per pre-print/deck theme and per campaign (~11)
  sessions/           current.md (the living handover, rewritten each sitting) + dated handovers
  graph3d/            index.html (three.js page) + graph.json (generated, git-ignored)
  assets/             a few decisive PNG/SVG (<= 200 kB); default is to link figures
  templates/          Obsidian templates (excluded from the site)
  .obsidian/          committed app/appearance/core-plugins/graph/templates .json; workspace ignored
  .quartz/            quartz.config.yaml, build.sh, check_kb.py, build_graph.py, README.md
```
Every folder has an `index.md` (title + two lines; Quartz uses it as the folder page). Slugs are
kebab-case; decision notes are named by the setting (`decisions/sl-clip-and-value-bound-off.md`), not by
the date; a note is never renamed or deleted (repurpose, add `aliases`); no file named `log.*.md` and no
single-digit folder (the repo's `.gitignore` traps).

**Frontmatter (all notes):** `title`, `aliases`, `tags` (the kind and `part/<part>`), `kind`
(hub|concept|model|decision|retraction|case|study|session|index), `status`
(settled|open|retracted|voided|candidate), `part` (advection|viscosity|surface-tension|mass-flux|
gradient-control|verification|all), `date`, `date_settled`, `decided_by` (config paths or "author
decision YYYY-MM-DD"), `code` (repo paths), `sources` (short locators: "STATUS 11.15", "METHOD 8.1 row
SL_CLIP", "gcls article sec:gate-verdict", "SL deck #/3/2"). Body sections of a concept/model/case/study
note: one-paragraph verdict with date (shown in popovers and search) / What it is / Why it matters /
Where in the code / Evidence (table: claim, number, where) / Why it failed, or why we think so / Decisions
/ Open questions / Related / `## Log` (append-only, dated, newest at the bottom). Decision note: the
question, the measurement that decided it (table with pre-registered read-out link), what it does not
cover, related, log. Retraction note: the claim as stated, where it lived, why it was wrong, what
survives, the propagation checklist (STATUS, METHOD, plan, deck/table, every KB note, the log line),
related, log. Hub: the question this part answers, current verdict (dated, every clause linked), map
table (kind, note, status, one line), open in order, log. Session note: branch, sha, stamp, last
handover, what is being worked on, open in order, traps, where the numbers live. STE prose; one screen
where possible; every number carries a link; a claim without a `status` does not exist in the KB.

**Links, one convention per target:** KB note → `[[folder/slug]]` (Obsidian `newLinkFormat: absolute`,
Quartz `markdownLinkResolution: absolute`; `[[x\|alias]]` inside tables). Article section → GitHub blob
URL pinned to a commit with the `\section` line anchor, text "SL article sec:surften", plus the repo path
and `\label` in backticks (greppable, drift-proof). Pre-print PDF → site URL `preprints/<name>.pdf`
(document level). Deck slide → site URL `decks/<name>.html#/h/v` (reveal `hash: true`; new decks add
section ids). STATUS/METHOD/CLAUDE/plans/roadmap section → GitHub blob on the branch with the heading
anchor, text "STATUS 11.15" (live, because these files receive CORRECTED notes in place). Code → GitHub
blob pinned to a commit with `#Lstart-Lend`, path in backticks. Gate/candidate config → pinned blob with
the pre-registered header lines. Curated table/figure → blob on the branch, file name as text; a decisive
figure may be copied to `assets/`. Data archive → tree URL on the branch, the stamp as text. Literature →
DOI.

**The two logs** (both, plus per-note `## Log`): global lines are never edited (a correction is a new
line linking the old one); month headings, newest month at the bottom; verbs SETTLED / CORRECTED /
REOPENED (decision log) and RETRACTED / VOIDS / FALSIFIED / CORRECTED (retraction log); line format
`date · [[note]] · VERB claim: number · gate/config · see [[...]]`. Every note in `decisions/` has a line
in `decision-log.md` and every note in `retractions/` in `retraction-log.md` (checked by `check_kb.py`).

**Quartz setup (single source, nothing vendored):** `.github/workflows/knowledge-base.yml` replaces the
dead `docs.yml` (Doxygen from a path that no longer exists): on push to `main`/`development` touching
`docs/**`, STATUS/METHOD/CLAUDE, and `workflow_dispatch`; Node 22, Python 3.12; `check_kb.py` gate;
`build_graph.py`; clone Quartz at `QUARTZ_TAG: v5.0.0` (never committed); `npm ci && npx quartz plugin
install && npx quartz plugin resolve`; rsync the vault into `content/` (excluding `.obsidian`, `.quartz`,
`templates`), copy the committed `.quartz/quartz.config.yaml`; stage `record/` copies of STATUS, METHOD,
CLAUDE, plans, roadmap (full-text search only; links point at GitHub); build the decks into
`content/decks/` (`docs/build-decks.sh`, python3 only); compile the pre-prints into `content/preprints/`
(texlive install, `continue-on-error`); `npx quartz build -d content -o public`; `upload-pages-artifact@v3`
+ `deploy-pages@v4`. `baseUrl: leia-openfoam.github.io/leia`. Enabled plugins: note-properties,
created-modified-date (frontmatter), syntax-highlighting, obsidian-flavored-markdown, GFM, TOC,
crawl-links (absolute, external links in new tab), description, latex (KaTeX), remove-draft,
alias-redirects, content-index (sitemap), favicon, content/folder/tag pages, explorer (left), graph
(right; local depth 1, global radial, showTags), search, backlinks, article-title, content-meta,
tag-list, page-title, darkmode, reader-mode, breadcrumbs, footer (GitHub, STATUS.md, METHOD.md),
stacked-pages. **Local build:** `make kb` / `make kb-serve` run `.quartz/build.sh` (same steps; uses
nvm's Node 22 if present, else the `node:22-slim` Docker image; Docker Desktop exists on this machine;
the system Node 18 is not upgraded). **`.obsidian/app.json`:** `newLinkFormat: absolute`, wikilinks,
`attachmentFolderPath: assets`; `graph.json` colour groups by kind tag; `.gitignore` ignores
`workspace*.json`, the plugin data and `graph3d/graph.json`, and un-ignores `assets/*.png|svg`.
**README.md:** the badge and the documentation line point at the site (`https://leia-openfoam.github.io/
leia/`); Doxygen stays local (`docs/api/Allwmake`).

**`check_kb.py` (CI gate, stdlib):** frontmatter fields and vocabularies; every wikilink resolves
(target, heading, alias, asset); every decision/retraction note has its log line; a settled decision has
`date_settled` and `decided_by`; no ignored file names; `## Log` is the last section. **Relation to the
record files:** a curated layer, not a replacement: STATUS.md stays the lab notebook, METHOD.md the
executable best configuration, the plans the programme documents, CLAUDE.md the rulebook; a number is
written first there or in a curated CSV, then quoted in the KB with its link. Later (phase 2): STATUS §0,
§1, §2, §3, §7 shrink to a pointer paragraph to `sessions/current.md` and the hubs.

**Public site:** the repository is public and Pages already serves it, so the knowledge base is public.
The notes describe the methods (which the pre-prints already publish) and cite the two dossiers by
title, author and date; no dossier text is copied into the vault.

### 1b. Coverage: the notes (about 75), what each records, and its sources

The coverage below maps onto the kind folders of 1a: hubs → `hubs/`; the topic notes → `concepts/`
(ideas and mechanisms) or `models/` (one per RTS family, members as rows); the "settled" items → also a
`decisions/` note (16, mirroring METHOD §8.1 row groups: sl-reconstruction-uncached-qwls,
sl-fit-normal-equations, sl-trace-velocity-projected-flux, sl-clip-and-value-bound-off, psi-filter-none,
mass-flux-rholent, mass-flux-alphaf-donor-plane, mass-flux-bound-rho (open), phase-indicator-
detrixhe-aslam, surface-tension-reconstructed-curvature, curvature-extension-cell-centre-inverse
(case-dependent), curvature-inverse-gaussian, viscosity-face-model-alg-lin (3D open),
momentum-schemes-bdf2-upwind, mesh-family-hexahedral, process-gates-2d-first-and-no-best-yaml); the
retracted items → `retractions/` (12: closed-box-translating-droplet, distance-cone-bound-as-transport-
bound, clip-damage-is-the-narrow-band, polyhedral-popinet-3d-mesh-defect, t-blow-baseline,
psi-filter-seam-bug, gradu-coupled-patch-contamination, cell-mean-delivery-adoption, force-at-n-not-n-
plus-1, advection-orders-3-2-factor, late-translating-instability-is-the-outlet, mass-momentum-
consistency-dominant-term); cases → `cases/` (stationary-droplet, translating-droplet, oscillating-
droplet, popinet-translating-droplet, ellipse-and-foliation-gates, kinematic-advection-cases,
exact-1d-stretch, curvature-static-gates); studies → `studies/` (one per theme: sl-quadratic-pre-print,
sl-linear-pre-print, grl-pre-print, sdpls-pre-print, npsl-design, velocity-extension-pre-print,
method-comparison, gcls-pre-print; per campaign: curvature-stabilization-campaign, shannon-parasitic-
currents-campaign, poly3d-roadmap); sessions → `sessions/current.md`, `sessions/2026-09-27-gcls-first-
campaign.md` (the SL handover of section 2 is a section of it and a note of its own,
`sessions/sl-session-handover.md`). First commit: index, conventions, both logs, 6 hubs, current, all 16
decisions, all 12 retractions, the 8 model notes with a settled decision, 6 concepts, 4 cases, 3 studies;
the rest follow in the same campaign.

Every note is concise (half a page to two pages), states what the thing is, where it is in the code,
the evidence with numbers, why it failed or why we think so, the decisions, the open questions, and links
to the article section, the deck slide group and the STATUS section that hold the details. The sources
below are from the docs inventory (paths under `docs/`; `SL` = `semi-lagrangian-level-set/sl-level-set-
article/semiLagrangianLevelSet.tex`, `SDPLS` = `sdpls-level-set/sdpls-article/sdplsLevelSet.tex`, `MC` =
`method-comparison/.../methodComparison.tex`, `GCLS` = the gcls pre-print, `VE deck` = the velocity-
extension deck, `SL-neg` / `SDPLS-neg` / `GRL-neg` = the negative-results decks).

**Hubs (8).** `Home` (map of content: the four moving parts, the method lines, how to read the graph,
the rules for writing here); `advection`, `viscosity`, `surface-tension`, `mass-flux`,
`gradient-control`, `verification` (one hub per part: the questions, the settled answers, the open ones,
the notes); `method-lines` (the seven lines with status: quadratic SL (production), linear SL, normal-
projected SL (closed), Eulerian SDPLS (kinematic positive, coupled negative), geometric redistancing
(closed), velocity extension (dominated; stub article), gradient control (first campaign failed)).

**Interface advection and phase indicator (16).** `semi-lagrangian-quadratic-transport` (SL §2, deck
"Semi-Lagrangian advection"; orders 2.97/2.59, 2.95/3.28, 1.36/1.46); `departure-foot-ab2-centring` (SL
§2.1; taylor vs rk2 equal on the coupled droplet); `quadratic-wlsq-reconstruction` (SL §2.2 and ¶429
admissibility; value fit, no constant term, pivot tolerance 0.3, CPC vs CFC stencils; QR bit-identical
to Cholesky); `polyhedral-fit-amplification` (Lambda 1.2608 on pMesh vs 1.05 on every hex family; the
pivot census; `sl_fit_amplification.csv`); `linear-semi-lagrangian` (LSL article; nestedLSQ ~1.1,
linearTaylor O(1e9) gradient defect; two-phase workhorse); `normal-projected-sl` (nPSL article §"Measured
verification record": trace clean, write-back diverges x1.7 per 10 steps, corrugation re-enters through
the fitted normals; what survives); `idec-defect-correction-failure` (SL supplementary: rho > 1);
`eulerian-fv-transport` (MC §"Faults", limiter drops order 3.0 -> 0.9; flux-form loses 17 % volume);
`trace-velocity-projected-flux` (METHOD 8.1, STATUS 2026-08-31: the reconstruct operator carries 70 %;
-52 vs +118 1/s); `psi-outer-correctors` (default yes since 5cbfaaa; 3 frozen = 12 re-advected);
`value-bounds-and-clips` (SL_CLIP, stencilBounds, lipschitzCone falsified; the mesh-noise floor
diagnostic); `redistancing-geometric-grl` (GRL article and GRL-neg: planeFootWave second order, PDE
reinit "stability bomb", frozen-band injurious, line closed by MC §304); `phase-indicator-detrixhe-
aslam` (SL §2.7, ¶1531: order 2, geometric = DA to 8 digits); `narrow-band` (signChange band, halo
dilation incident); `static-local-refinement` (SL §"Local static refinement": 3.4-6.9x fewer core-hours,
hanging nodes measured); `advection-regression-set` (CLAUDE.md: hex 2D, hex 3D, poly 3D; mesh-noise
floor 0.0 / 0.0 / 14.2 %).

**Viscosity (2).** `face-viscosity-models` (SL §"The viscous term and its face viscosity" ¶794-875: the
six-name grid, alg_lin the only positive-order model, "viscosity spreads the spurious current rather than
damping it", the frozen-muf fix); `viscosity-open-items` (3D templates without the token; geometric
sharp alpha_f consistent with the plane).

**Surface tension (14).** `balanced-force-csf-flux` (SL deck "Two-phase flow": momentum matrix, forces
in flux space, what survives the projection; SL §"What the pressure projection can absorb": 99.85-99.98 %
absorbed, residual structural); `curvature-from-the-fit` (METHOD §4, SL §"Curvature accuracy": symbolic
kappa, the parallel-curve offset correction, 11.5 % -> 0.46 %, second order on a clean circle);
`cell-centre-inverse-curvature` (cellCentreInverse + K-aware inverse; order +2.01 non-gradient content
on constant curvature; 4.60/3.81/1.55/1.49x lower residual; constant-curvature caveat; 3D h^1.95 vs
h^1.02); `face-curvature-deliveries` (plan-curvature-stabilization §8-14: stabilizedFootPointFace,
cutCell (blow-up 3.3x sooner), cellMean, footPointEvaluated, symmetricFaceMean, connectedInterface; gain
G h^2 vs accuracy dissociation; the ellipse gate collapses one-value-per-cell deliveries to first order);
`parasitic-current-mechanism` (SL §"stationary droplet" and §"translating droplet", plan-shannon §0d:
source = curvature estimator (step-1 kick independent of U0), amplifiers = translation and density
ratio, the two-factor law max|U| = u0(h) exp(G(h)), t_blow ~ N^-3/2, exact-curvature arms bounded);
`curvature-corrugation-and-the-fit` (m > 4 modes, grid-scale aliasing, the psi filter as instrument only:
5.86x better / 1.61x worse); `semi-implicit-capillary-force` (Hysing/Raessi/SAAMPLE fvOption; the
collective-in-master-guard deadlock; -45.3 vs -52.0 1/s: not needed with projectedFlux); `integral-
surface-tension-cst` (SL-neg "Force delivery": sign trap, decaying equilibrium 5e-7 at N=64, diverges at
N=128; the meta-law "better static balance, higher dynamic gain"); `kang-gfm-and-sharp-heaviside`
(rejected: 4.4e-3 rising vs 1.3e-5); `force-time-centring` (MC deck "Time centring": the amplification
matrix, endStep vs midpoint; the 2026-09-27 midpoint test 32 % earlier); `capillary-time-step`
(0.2323 Brackbill; the dt sweep: ~90 % of the growth dt-proportional at N=128); `pressure-projection-
and-linear-solvers` (roadmap 2026-07-28 gates: operator pair, rAUf, tolerance, solver; non-orthogonal
caveat); `variational-capillary-force` (SL §Outlook: the proposal, no code); `well-balanced-exact-
curvature-gate` (Delta p = 145.470 vs 145.48; spurious velocity <= 3.8e-10).

**Mass flux and density-ratio consistency (9).** `rholent-mass-flux` (METHOD §6, SL ¶2015, roadmap
"Gate 2": the auxiliary density equation, reset after the loop, residual 1e-13; stationary
+1.0/-22/+0.1 %); `mass-flux-comparison-modes` (interpolatedDensity, geometricFaceDensity 4.8e-2
residual); `bound-rho` (active with rhoLENT; the only evidence VOID; clip reported, rhoClipFraction
counts round-off); `alphaf-source-donor-plane` (donorPlane by the author's instruction; averagedPlanes,
donorPlaneAdvected, trapezoid time level); `mass-flux-projection` (curl-free correction; compatibility
needs zero mean); `ddt-scheme-pairing-bdf2` (the BDF2 rule; matching argument; the pairing tables VOID);
`eulerian-solver-mass-flux-port` (frozen rho until 2026-09-27; 0.42 vs 1.0001 travelled fraction; the
shared headers); `coupled-face-density-defect` (2026-09-27: rho_f 90 % apart across seams; the fix;
1e-4 -> 1e-8); `closed-box-void-2026-09-02` (the wrong setup, what it voided, the free-stream gate).

**Gradient control: source terms and velocity extension (13).** `gradient-control-overview` (the
direction, the two dossiers (not in git), the plan, the pre-print, the technical report);
`source-law-family` (Table of laws, strain weights; F(1) = 0; the fixed point 1 - a/mu);
`sdpls-source-eulerian` (SDPLS article: R converges kinematically (+0.74 vs -0.26, 31x), beta's target is
wrong (fixed point beta - a), Rdiv withdrawn (diverges kinematically), exponentialImplicit repairs
order; coupled: mode-4 current amplified 260x, seed = curvature error, exact kappa gives 0.000 in 12/12);
`sl-source-step` (the exponential update, band 3h, clamp 30; the zero-set shift delta = h ab/(a+b)^2 dt
(F_A - F_B); the S1/HL1 gate results); `combined-source-note` (F = a + lambda(1 - q): invariant interval;
subsumed); `halo-limited-extension` (design: S_R, w, Y, flux correction, stencilFit sampler; unit gates
pass 2e-16; HL0 order 1.39, coupled arms destroyed; R = h not mesh-resolved); `closest-point-extension`
(FP0; decomposition dependent 45 %; MC: VE "dominated"; the DISPUTED 2185x claim);
`legacy-velocity-extensions` (VE deck verdicts for none, anisotropicDiffusion, pseudoTime, steadyUpwind,
steadyUpwindLinear (divergent cascade), meshWave; static (n.grad)Uext = 0 test; one-sided errors);
`material-form-transport` (the Sp(div Phi) correction; Test 0a; what SDPLS Rdiv taught);
`method-gate-2d-campaign-2026-09` (the six candidates, verdicts table, per-arm numbers, seam checks);
`why-the-candidates-failed` (the mechanisms with their evidence and status: hypothesis / measured; links
to the technical report); `gradient-control-next-experiments` (the ranked discriminators from the
report, each with its pre-registered prediction); `gradient-control-open-decisions`. Added after the
red-team memo (2026-09-28): `source-discretisation-defects` (the centred-q eigenmode at rate μ, the hard
band edge under the accumulated O(1) factor, the sublattice-incompatible fixed point, the soft wall's
40σ slope, the clamp that is never reached; each with its mark and its discriminator E1.1-E1.3, E1.7,
E2.1, E2.2); `extension-strain-relocation` (the K(t) profile of the capped/fractionReached extension, the
shell at 1-2h with the 27 % overshoot, the 1D `qBandMean` 0.46-0.59 against the 0.49 prediction, the
radius cap that blocks the R ladder, the tangential shear of every normal-constant extension); and the
gate's two blind spots in `method-gates` (the vacuous exact-1D closed form for candidates; the centred
band metric) with the repairs listed in the report.

**Verification methodology (11).** `method-gates` (2D/3D definitions, arms, criteria 1-4, INVALID);
`error-vector-and-read-out-instants` (L2/L1 only, never L_inf; gradient at T/2, shape at T, volume at
both; the oscillating arm by period and damping); `richardson-ladders-and-orders` (three rungs minimum,
h ratios, Celik/GCI, the fourth-rung falsification); `seam-checks-and-decomposition-invariance`
(kinematic and coupled checks; the class of one-rank-correct defects: setVelocity, updateFlux, the psi
filter patches, the band dilation, the face density, the metrics; the gradU contamination of
2026-08-26); `wrong-setup-voids` (closed box, tilted wall faces, algebraic psi, frozen rho; rename and
re-run); `retraction-log` (chronological, one line per retraction with the note it belongs to);
`decision-log` (chronological, one line per decision: date, decision, decided by, note);
`bit-identity-and-inertness-gates` (compare_metrics_csv, DICT_MODE, the pre-extpoints set, "nothing is
inert until measured"); `cluster-provenance-and-binaries` (per-clone binaries, stamps, ledger, the
touch sweep, the ssh-output rule); `log-classifier-and-waiters` (foam_log_state.sh; the trapFpe false
positive); `data-archive-per-version` (the pre-print archive convention).

**Cases and studies (2).** `benchmark-cases` (one table: every case under `cases/` with geometry,
closed form or reference, which parts it tests, its gate figure); `studies-index` (config/*.yaml grouped
by topic, with the docs table each feeds).

**Sessions and handover (3).** `session-gradient-control-2026-09` (what the session built and changed,
by commit); `sl-session-handover` (what the SL session must know: the two parallel fixes, the metric
fixes, the Eulerian port, the translating-droplet box study, the oscillating drift at N = 200, the gate
infrastructure and the seam check, the rhoClipFraction and oscillating scoring corrections, the record
corrections of METHOD.md, and the open decisions listed for the author); `how-to-write-here` (the rules
of the knowledge base for humans and sessions).

### 1c. The interactive 3D knowledge graph (three.js, on the same GitHub Pages site)

Decided by the user (2026-09-28): the interactive three.js graph is part of the first deployment, not a
later addition. It is one static page inside the vault, so Quartz publishes it with the notes (Quartz
copies every non-Markdown file under `content/` to `public/`), and no second tool chain exists.

- **`docs/knowledge-base/.quartz/build_graph.py`** (Python standard library only; committed; shares the
  frontmatter and wikilink parser of `check_kb.py`): reads every
  note, parses the frontmatter (`title`, `part`, `status`, `tags`, `description`) and the `[[wikilinks]]`
  (aliases `[[note|text]]` and heading links `[[note#h]]` resolve to the note), and writes
  `docs/knowledge-base/graph3d/graph.json`: `nodes: [{id, title, part, status, tags, description, degree,
  url}]` and `links: [{source, target}]`. Unresolved links are reported (they are errors of the vault: the
  script exits nonzero, so the CI build fails on a broken link). `url` is the Quartz URL of the note
  (extensionless path relative to the site root; GitHub Pages serves `<path>.html` for it) and `file` the
  local `.md` path, so the page works both on the site and in a local preview.
- **`docs/knowledge-base/graph3d/index.html`** (committed, plain HTML + JS): loads `3d-force-graph` from
  jsDelivr (`https://cdn.jsdelivr.net/npm/3d-force-graph`; it bundles three.js), fetches `graph.json`, and
  renders: node colour by `part` (advection, viscosity, surface-tension, mass-flux, gradient-control,
  verification, hub: one fixed palette, shown in a legend); node size proportional to the square root of
  the degree; `status` by opacity and ring (settled solid, open outlined, retracted and voided grey);
  hover: title, part, status and description; click: opens the note (site URL, or the `.md` in local mode
  `?local=1`); a legend with checkboxes that filter parts and statuses; a search box that focuses the
  camera on a node; a "hub only" toggle that hides leaf notes. Directed links are drawn with a small
  arrow. The page has a `<title>`, works on a phone width, and needs no build step of its own.
- **In the vault:** `Home.md` links to `graph3d/` ("the graph, interactive"), and Quartz's own 2D graph
  stays in the right sidebar of every note; the 3D page links back to `Home`.
- **In the CI workflow:** `python3 docs/knowledge-base/.quartz/build_graph.py docs/knowledge-base` runs before the Quartz build
  (the JSON is generated, not committed; `.gitignore` gets `docs/knowledge-base/graph3d/graph.json`).
- **Locally (no Node needed):** `make kb-graph` regenerates `graph.json` and prints the command
  `python3 -m http.server -d docs/knowledge-base 8000` so the 3D page opens at
  `http://localhost:8000/graph3d/?local=1` from the raw vault.
- **Verification:** the script's link check passes on the vault; the page renders the graph with all notes
  as nodes (count equals the number of notes), the filters hide and show parts, a click opens the note on
  the deployed site; checked in a browser before the first deployment (rendered pages inspected).

## 2. Separating the semi-Lagrangian topic from the source-term topic

The SL article (`docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex`) holds
no source-term content; the mixing is in the gcls pre-print, which carries four SL-baseline findings.
They move, with their data, so that each paper owns one topic:

| gcls pre-print section (moves out) | SL article destination | data |
|---|---|---|
| §5.6 "Consistency under domain decomposition" (the face-density defect, the diagnostics defect, Tables `gcls_face_density`, `gcls_seam`, Figure `gcls_seam.pdf`) | new subsection "Consistency under domain decomposition" at the end of §"Software design and reproducibility" (L1101), with the earlier one-rank-correct incidents (setVelocity, updateFlux, gradU contamination) as the class it belongs to | the tables and figure move to the SL article `data/` (the archive folder stays with the gcls theme and both READMEs cite it) |
| §5.6 ¶"The Eulerian two-phase solver" (Table `gcls_eulerian`) | new subsection "Mass–momentum-consistent flux (rhoLENT)" in §"Two-phase flow coupling" (the SL article has no mass-flux subsection; METHOD.md §6 is the source of the formulation), with the Eulerian port as one paragraph | `gcls_eulerian.tex` moves |
| §7.8 "The translating droplet at long times" (Table `gcls_translating_runs`, Figure `gcls_translating_boxes.pdf`) | extends §"Two-phase verification: the translating droplet" (L1901) with a paragraph "The late instability: the outlet and the interior growth" (the box-length study, the start-position test, the 40 mm ladder at two rungs) and corrects §"Limitations" accordingly | table and figure move |
| §7.1 ¶"Oscillating droplet" (baseline drift at N = 200; Figure `gcls_oscillating_drift.pdf`) | a paragraph in §"Limitations and outlook" (the gradient drift without redistancing turns unstable at N = 200 within ten periods) | figure moves |

What stays in the gcls pre-print: the baseline as the reference of the gate (its error vector table
stays, with one sentence per arm and a citation of the SL article for the long-time and parallel
findings), the source terms, the extension, the gate, the candidate results, the discussion. The
pre-print's abstract, contributions (items 4 and 5 removed) and conclusions are rewritten accordingly;
`figures/make_result_figures.py` keeps generating the moved tables and figures but writes them into the SL
article's `data/` (one script, two output folders, both from the same archive), and the SL article's
`REPRODUCE`-style provenance names the archive.

**The pre-print's discussion is corrected to the report's assessment** (section 3): the "coupled loop
through the velocity at rate μ" reading for HL1q/HL1z is withdrawn (the linear laws read no velocity:
`gradientControlLaw.H:182-185` F = w·a + G, w = 0 for every candidate, no σ in linearQ/Z; the measured
growth rates 245-276/s match μ = 270/s and the flow loop is the second stage); the crossing shift is
demoted to an O(h) contribution (the shape errors are O(1) in h, orders 0.02-0.06); the `reconstruct`
argument to one O(h|∇u|) contribution (FP0 has no shell and loses the same order); the
sampler-across-the-jump reading is kept for the oscillating arm only, and the translating arm's reading
becomes the unprojected injection of the spurious current. A CORRECTED line goes into the retraction log
(`retractions/gcls-coupled-loop-reading`), and the KB note `why-the-candidates-failed` carries both
readings with their marks and the experiments that decide them.

The handover for the separate SL session is the knowledge-base note `sl-session-handover` (section 1b)
plus a short pointer at the top of STATUS.md §11.13-11.15 ("SL findings of this session: see the
handover note"). Nothing is deleted from STATUS.md.

## 3. The technical report (LaTeX, PDF)

**Path:** `docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex` (class
`article`, 11 pt, the packages of the pre-print, `\input` of the same `data/tables/*.tex` and
`\includegraphics` of the same `data/figures/*.pdf` by relative path, so the report and the pre-print read
one archive), built by `make report-gcls` (latexmk, like `article-gcls`; PDF git-ignored). Title:
"Gradient control by source terms and velocity extension in the leia semi-Lagrangian level set: why the
first campaign failed, and what to do next. Technical report, 2026-09-28." Written as the expert
developer's assessment addressed to the supervisor; ASD-STE100 style; every number from the archive.

**Structure.** Sections 3-6 were rewritten on 2026-09-28 after the red-team memo. Its code-path claims
were verified against the source of `leia-gcls` at 8867581: `gradientControlLaw.H:182-185` F = w(q²)·a +
G(q², σ), w = 0 for every candidate, no σ in linearQ/linearZ, so the linear laws READ NO VELOCITY;
`slGradientControlSource.C:146` q from `fvc::grad(psi, "gradPsiSource")`, which the case templates set to
`leastSquares` (a centred gradient); `:170` the hard band switch `|ψ_c| > 3 h_c q_c → skip`; `:185-194`
clamp and exponential update; `make_gate_summary.py:190-198` `q_exact = None` for every candidate;
`slAlphaEqn.H:71` extension evaluated before the ψⁿ restore at `:149`; `haloLimited.C:108-114` the
radius cap. Every claim in the report carries a mark: MEASURED (a number in the archive), DERIVED (a
closed-form result, re-derived by hand in the appendix when the report is written) or HYPOTHESIS (a
discriminating experiment named).

1. **Summary** (two pages, first). The verdict in one paragraph: the first campaign failed on (i) the
   source DISCRETISATION, not the laws; (ii) an extension that relocates the normal strain into a shell
   instead of removing it; (iii) a gate with two blind spots. One line per candidate: what it is, how it
   failed, the mechanism held responsible, its mark. The three baseline findings that bound any
   candidate. The recommendation: abandon / repair / pursue, and the next experiments in order with
   pre-registered predictions and cost.
2. **What was tested and how.** The gate in half a page. What the dossiers prescribed and was skipped: the
   Tier-I rigid-motion and homogeneous-strain tests with q0 ≠ 1, and the zero-flow static test. Two gate
   defects found on the data: the exact-1D arm has no closed form for any candidate (`qError = None`, so
   criterion "target" is vacuous for candidates); the band-gradient metric (`gradPsiMetric leastSquares`,
   centred) is blind to the cell-scale mode that destroyed the stationary droplet. The kinematic arms
   carry no curvature column. The soft wall was inactive on the stationary droplet (σ ≈ 1e-3/s and
   |q − 1| < δ_s), so "S1 = baseline to four digits" says nothing about coupled safety.
3. **The source-term route.**
   3.1 The source's own map, with no flow in it (DERIVED; discriminator E1.1). Linearise ψ_c ←
   ψ_c exp(−μ dt (q_c − 1)) about ψ = d: δ_c ← δ_c − μ dt d_c (n·∇_h δ)_c. With a centred q the mode
   δ_c = (−1)^c E d_c has (n·∇_h δ)_c = (−1)^{c+1} E and grows by (1 + μ dt) per step: rate μ,
   independent of dt, h and the flow (2D: (−1)^{i+j} E d). A one-sided (Godunov/Rouy–Tourin) q damps
   the same mode at 2μ d/h. At u ≈ 0 the SL step is the identity (`evaluateRaw` at the cell centre
   returns ψ_c), so nothing damps the mode; at CFL 0.5 the quadratic fit halves a checkerboard per step,
   which is why the shear arm hides it. SDPLS never saw it: R reads a, not q; beta is linearQ at μ = 1/s.
   3.2 The band edge under an accumulated O(1) factor (DERIVED; discriminators E1.2, E1.3). A law that
   reaches its target multiplies |ψ| inside the band by e^{∫F dt} = q_far/q_target ≈ 2-3 (shear arm) and
   by 1 outside: |ψ| is non-monotone at 3h, independent of h, and spurious crossings follow (MEASURED:
   HL1z band error 10.5 and volume error 1.6 at T; HL1q band error 0.44 at t = 1.9 s → 6.58 at T). With
   a centred gradient the fixed point "q = 1 in the band, outside data c·d" splits into two sublattices;
   for c = 1.5 the odd sublattice's fixed point has the wrong sign at the interface cells, the
   sign-preserving update drives |ψ| → 0 there, and the zero set moves O(0.5h) with u = 0. The z-law's
   −μq²/2 is unbounded for the large q this produces (runaway at t ≈ 0.7-1.0 s in shear).
   3.3 HL1q/HL1z on the stationary droplet (MEASURED rates; DERIVED mechanism). The linear laws read no
   velocity, so the pre-print's loop ψ → κ → f → u → ψ → q → F cannot be primary. The curvature error
   e-folds at 245/s and the spurious current at 276/s against μ = 1/T_ref = 270/s; the mesh-scale
   capillary rate (~1e5/s) is not seen; HL1q and HL1z agree to three digits until 0.02 s; the curvature
   error is 10x baseline at 4.3 ms while the velocity is 1.4x; the centred band metric stays within 2x
   baseline until ~10 ms while the curvature error is 100x; the smooth part of the q error decays at μ
   as the continuum says (3.6e-4 → 2.2e-4 in 2.3 ms, e^{−0.62} = 0.54); min/max q leave [0.99, 1.01]
   symmetrically (0.83/1.17 at 20 ms, 0.33/2.29 at 30 ms). Seed: the cell-to-cell part of the initial q
   error, O(h²/R²) ≈ 6e-4 at N = 200, reaches O(1) at ln(1/6e-4)/270 = 27 ms; observed 25-35 ms. The
   flow loop is the second stage (t > 10 ms). The clamp |dt F| = 30 is never reached (HL1q at q = 190:
   dt F = −0.2); for linearQ q → 0 gives F → μ, bounded.
   3.4 S1 (MEASURED seam; DERIVED slope). The soft wall's slope |dG/dq| reaches 1.25σ·2.57/δ_s ≈ 40σ at
   |q − 1| ≈ 0.069 (`gradientControlLaws.C:210-220`), so dt·dG/dq ≈ 0.5 per step in the shear arm
   (σ ≈ 3/s, dt 4.1e-3 s): a nearly discontinuous per-cell map. The 76 % seam failure is its signature;
   HL1q/HL1z pass the same coupled seam check at 4.7e-6, so "round-off in q_c" is not their mechanism.
   In the coupled arms the soft wall reads σ of the spurious current: that, and only that, is the SDPLS
   loop of the pre-print.
   3.5 What the pre-print's mechanisms explain and do not. The crossing shift δ = h·ab/(a+b)²·dt·
   (F_A − F_B) is first order in h: summed over the shear run |δ| ≤ (h/4)·μT·max|Δq| ≈ 0.2h, while the
   measured HL1q shape error 0.107 is a ~1.1h displacement of the whole perimeter with order 0.02 (flat
   in N); S1's 0.265/0.327/0.254 is flat too. The band-edge kink at 3h lies outside the interface cells'
   CPC stencils (reach ≈ 2.4h). The nSL result (per-cell write-back from the FITTED gradient) does not
   transfer: the source reads `fvc::grad`, not the fit. Kept: the shift as an O(h) contribution; the seam
   signature for the soft wall.
4. **The velocity-extension route.**
   4.1 R = h relocates the strain, it does not remove it (DERIVED; MEASURED on the 1D arm). For `capped`
   travel and `fractionReached` β = 1: S = c d, w = c, c = (1 + t⁴)^{−1/4}, t = d/R; the transmitted
   normal strain K(t) = 1 − (1 − t⁴)(1 + t⁴)^{−3/2} = 0, 0.14, 1.00, 1.27, 1.21, 1.11 at t = 0, 0.5, 1,
   1.5, 2, 3, with ∫(1 − K) dd = R²/D → 0: a shell at d ≈ 1-2h with a 27 % overshoot, inside the
   geometry-evaluation neighbourhood (band 3h, stencil 1.4h), which the dossier itself forbids ("the
   transition should be broad, located outside the geometry-evaluation neighborhood"). Measured:
   exact-1D `qBandMean` HL0 0.46-0.59 against e^{−1} = 0.368 (baseline), 0.96-1.09 (FP0) and the
   K-profile prediction 0.49; shear after 5 steps HL0 band error 0.020 vs 0.035 (baseline) vs 0.001
   (FP0). The pre-registered HL0 read-out conflated the interface with the band. The q profile's O(1)
   variation over one cell (ψ''' ~ 0.6/h²) makes the quadratic reconstruction error at the feet ~0.1h
   instead of ~4e-4 h: first-order transport whatever the trace velocity (HL0/baseline shape ratio
   4.3 → 8.4 → 14x over the rungs).
   4.2 Within-cell variation of the flux correction (the `reconstruct` argument of the pre-print) is
   neither necessary (FP0 has no shell, `phiExt = interpolate(Uext)·Sf`, and loses the same order 1.38)
   nor sufficient (HL0's excess volume error 2.3e-4 vs 5.5e-5 appears by t = 0.02 s, before any q-kink
   exists). Kept as one O(h|∇u|) contribution.
   4.3 The translating droplet has no physical gradient to fit across (rigid translation, D = 0), and HL0
   fails 64x there (order −3.1): the sampler reads the spurious current u' (cell-scale, ~6e-4 m/s) and
   injects w[u_h(Y) − u_h(x)] ~ u' into the trace flux UNPROJECTED, while the baseline traces
   `reconstruct(phi)`, the smoothest velocity the solver has. HL0 tracks the baseline to 0.02 s
   (κ-error 74 vs 62 /m), then κ-error and velocity grow together at 183/s and 178/s (MEASURED). The
   oscillating droplet is where the sampler-across-the-jump reading is plausible (HYPOTHESIS): the
   air-side boundary layer δ ≈ √(ν T) ≈ 8h at N = 200 gives ∂_n u_t ≈ 250/s, the fit extrapolates the
   water side by up to 25x, and HL0's band gradient error grows at 450/s, five times the physical
   strain; Y_f can lie 1.5h from a cell centre, outside the 3x3 hull (`evaluateRaw` has no trust
   region). The ratio-1 test of the pre-print cannot separate the two readings (it shrinks u' and the
   kink at once); sampling `reconstruct(phi)` can (E1.6).
   4.4 FP0: no shell, the same order loss, seam-dependent (45 %), collapses at the shear tail (0.1-0.4 s:
   the closest-point search fails, the fallback steady solve, the distance function's kinks within 3h).
   Any normal-constant extension shears the band tangentially: level sets at distance d rotate at
   Ω(d) = ω(R + d(1 − c²))/(R + d) (HL) or ωR/(R + d) (FP); a circle is invariant, a slotted disc is
   not (the Tier-I rigid-rotation prediction; HYPOTHESIS until run).
   4.5 Time-level mismatch on outer passes 2-3 (coupled solver only): `slVelExt->correct()` runs before
   the ψⁿ restore, so the extension geometry (d, e, w, S, Y) comes from ψ^{n+1,(k−1)} while ψⁿ is
   transported; with R = h a CFL·h shift in d changes w by up to 0.3 (w = 0.84 at h, 0.50 at 2h).
   Discriminator E1.5.
   4.6 The R ladder is blocked by construction: `radiusCells > 1` is a FatalIOError because the sampler
   evaluates the face's own cells' models; R = 2-3h needs a sampler that locates the cell containing Y
   with a 2-3 layer halo, the halo the design set out to avoid; enlarging R makes HL FP-like in cost.
   4.7 Checked and consistent: AB2/UextOld, the departure-centred foot, the physical patch values, the
   source's time level.
5. **The baseline's own limits and the gate's blind spots.** Curvature of the moving interface does not
   converge; oscillating drift at N = 200; the translating interior growth; the outlet. The gate defects
   of section 2 restated as repairs.
6. **Assessment and recommendation.**
   Abandon as formulated: HL at R = h as a gradient-control device; HL2; the soft wall as parameterised
   (p = 5, δ_s = 0.08 is bang-bang; safe on the stationary droplet only because it was off); FP0 as a
   production candidate (reference only); the Taylor/weighted-SDPLS strain route (weighted SDPLS at the
   PDE level, dossier lines 1596-1607; a weight w ≤ 1 reduces the coupled gain by at most w and w = 1 in
   the geometry cells by design; the SDPLS coupled failure is seeded on the force side, which no source
   fixes).
   Pursue, in order: (1) repair the source discretisation before testing any law: a one-sided
   (Rouy–Tourin) q on unstructured meshes with coupled-patch neighbours (E2.1), a smooth band taper or no
   band (E2.2), a source-CFL guard μ·3h·dt/h < 1/2; verified by E1.1-E1.3 and the non-vacuous 1D closed
   forms (E0.1); then the stationary droplet at M_μ ∈ {0.1, 1, 10} with the prediction "no growth at
   rate μ" (a slower residual growth is the force-side seed and no source fixes it); (2) the
   velocity-free linear law (linearQ, strainWeight none) as the only family with a clean linearisation
   (offset q* = 1 − a/μ, so μ ≈ 10/T_ref in strained arms: the discretisation must be monotone first);
   (3) F evaluated at the algebraic foot x_c − d n and copied along the normal (E2.4) if the O(h) shift
   remains measurable after (1); (4) S2 (scalar Q memory) only after (1)-(3), in the variant that filters
   q_h with the upwind q, and with its a_h/σ_h term off first; (5) ALG (the uncapped algebraic sample,
   w = 1 in the band, the dossier's K taper with ℓ₀ = 3h, ℓ₁ = 6h, K_p ≈ 1.8) only if an extension is
   still wanted for droplets (it does not help the vortex); (6) V1 last, only if (1)-(4) fail for lack
   of orientation memory.
   Gate repairs regardless: the candidate closed forms in the 1D arm; a checkerboard-sensitive band
   diagnostic (L2 of the second difference of ψ, or min/max of a one-sided q); a curvature column in the
   kinematic arms; the zero-flow static arm as a gate arm; the extension's time level on outer passes
   ≥ 2. The rule stands: no candidate re-enters the coupled gate before the kinematic ladder (E1.1-E1.3,
   the Tier-I cases) passes at the baseline's order.
   The experiment ladder, cheapest first, each with its prediction and the outcome that falsifies it:

   | id (cost) | experiment | prediction | falsifies |
   |---|---|---|---|
   | E0.1 (hours) | candidate closed forms in the 1D arm: dq/dt = q(F − αK(d/R)), dd/dt = αd(1 − c²) per band cell | HL0 band mean ≈ 0.49 (measured 0.46-0.59); S1 → 0.925 (measured 0.998 → 0.697); HL1q with μ = α: band mean → ~0.5 (measured 0.80 → 0.53) | nothing yet; it makes the 1D gate real |
   | E1.1 (dictionary) | zero-flow static test: kinematic solver, u = 0, circle SDF, N = 100, linearQ, bandCells 3, μ ∈ {27, 270, 2700}/s, μ dt ≤ 1e-2 | max\|ψ_c − d_c\| grows as e^{μt} with pattern (−1)^{i+j} d; fit curvature error grows at μ; centred metric stays O(seed) until amplitude O(h); rates 1:10:100 | A: if nothing grows at u = 0 the source map is stable and the flow-loop reading returns |
   | E1.2 | E1.1 with ψ₀ = 1.5d (dossier B1, Tier-I test 4) | \|ψ\| → 0 in interface cells on one sublattice within ~5/μ; zero-set displacement O(0.5h); c = 1.08: 0.1-0.2h | B, if the zero set stays put to O(h²) |
   | E1.3 | E1.1/E1.2 with bandCells 1000 | B disappears, A remains | separates A from B |
   | E1.4 | HL1q stationary at GC_M_MU 0.1 and 10 | survives 0.1 s / dies by ~3.5 ms; rate ∝ μ | A as the primary mechanism |
   | E1.5 | HL0 translating with psiOuterCorrectors no | the departure from the baseline moves later or vanishes if D matters | D |
   | E1.6 | HL0 translating and oscillating with the extension sampling reconstruct(phi) instead of U (`velocityExtension.C:49-50`) | translating improves, oscillating does not; then ratio 1 for the oscillating arm | 4.3 against the jump reading |
   | E1.7 | S1 with δ_s = 0.3, p = 1 (slope ≈ 6σ) | seam ≤ 1e-5, shear shape error ≥ 10x lower, band-error gain shrinks | the bang-bang reading, if the seam still fails |
   | E2.1 (days) | upwind q in `slGradientControlSource::apply` + the source-CFL guard | E1.1's mode damped at 2μd/h; HL1q stationary at baseline level for 0.1 s | A |
   | E2.2 | smooth band taper F·τ(\|ψ\|/(3hq)), τ a C² cutoff, or no band | HL1q/S1 shear shape orders rise from ~0 to ≥ 1 | B |
   | E2.3 | kinematic `cellCentred` trace with Uext (`leiaSemiLagrangeLevelSetFoam.C:218`) | HL0 shear constant drops, order stays < 2 (the shell remains); order 3 would mean the trace alone was the cause | 4.1 |
   | E2.4 | F at the algebraic foot, copied along the normal | the O(h) shift vanishes at first order; A and B remain unless E2.1/E2.2 | 3.5 |
   | E3 (weeks) | extended sampler with `findCell` and a 3-layer halo for R = 3h; S2 | only if HL survives E2.3 | — |
7. **Appendix:** the data (archive folder, MANIFEST), the commit ids, the exact configs, and the two
   closed-form derivations (K(t) of the capped/fractionReached extension; the checkerboard eigenmode of
   the centred exponential source), one page each.

The report's numbers are the archive's; the mechanisms marked DERIVED and HYPOTHESIS are stated with the
evidence for and against and with the discriminating experiment that decides them, not as facts.

## 4. CLAUDE.md and AGENTS.md (byte-identical)

Two additions, after "Repo layout & git discipline":

1. **"The knowledge base is the point of reference."** `docs/knowledge-base/` (an Obsidian vault,
   published by Quartz) holds the concise, cross-linked record of every method line and moving part:
   what was decided, on which measurement, what was retracted, why something failed or why we think so,
   and the open questions. Start a session there (`Home.md`, then the hub of the part you work on, then
   the notes it links), and follow its links to the pre-prints, decks, STATUS.md sections and code for
   details. When a session takes a decision, retracts a claim, closes or opens a question, or finishes a
   gate, it writes the note (or the log line) in the same commit as the result; the `## Log` of a note is
   append-only. The graph (`graph3d/`, and Quartz's graph view) shows what relates to what; a note
   without links is not finished.
2. **"How we work: supervisor and expert developer."** The user is the supervisor; Claude works as an
   expert developer of numerical methods for multiphase flow. For every method line of leia (the
   semi-Lagrangian transport, the Eulerian lines, the surface-tension and mass-flux models, the
   gradient-control direction, and any new one): propose with the mechanism named and the prediction
   pre-registered, run the cheapest discriminator first, assess the results critically against the
   knowledge base (what was tried, what failed and why), state which routes to abandon and which to
   pursue with reasons, and write that assessment as a technical report in the theme's docs folder
   (summary first, details after) whenever a campaign closes; decisions are taken together, and the
   knowledge base records them.

## 5. Verification

- Vault: every note has frontmatter with the required fields; `python3 docs/knowledge-base/.quartz/
  check_kb.py docs/knowledge-base` and `.quartz/build_graph.py` exit 0 (no unresolved wikilink, every
  decision and retraction has its log line), the node count equals the note count; Obsidian opens
  the folder as a vault and shows the graph (checked by the user); a spot check of ten notes against
  their cited sources (numbers and section references).
- Site: `make kb` builds Quartz locally (with a Node >= 22 from nvm, installed on request) or the CI
  workflow builds it on a push to the branch (workflow_dispatch first); the deployed site shows the notes,
  backlinks, search, the 2D graph and the 3D graph page; links to decks and PDFs resolve.
- Papers: `make article-sl article-gcls report-gcls` build without undefined references; the moved
  tables and figures render in the SL article; the gcls pre-print no longer contains the moved sections;
  `git grep` finds no dangling `\ref` or `\input`.
- Report: every MEASURED number is checked against the archive CSVs named in the MANIFEST (the memo's
  numbers are its reading of the archive and are re-read, not copied); every DERIVED claim is re-derived
  by hand in the appendix (K(t); the eigenmode growth (1 + μ dt); the soft-wall slope 40σ; the summed
  crossing shift ≤ 0.2h); every code citation is checked against the pinned commit, line by line.
- CLAUDE.md and AGENTS.md: `diff CLAUDE.md AGENTS.md` prints nothing.
- Commit and push the branch; pull the Lichtenberg clone (no job of this session runs there now).

## 6. Record corrections found by the audit (small, done in the same campaign)

The extraction agent found 13 inconsistencies between the record files (Appendix A, section H). The
knowledge base states the current truth for each; the following are also corrected in place:
- `config/gates/methodGate2D.yaml:69` still cites the VOID `rhoDdtGate2D` as the basis for boundRho; and
  its 4.8e-2 residual comment comes from a pre-fix run (say so).
- `METHOD.md:35` names `libleiaLevelSet` (nine libraries now); `METHOD.md` §8 table and §9 items 1, 5, 6
  describe 2026-07-31 (mark CORRECTED with pointers, as done for §4.1/§10); `METHOD.md:375` names the
  wrong gate for projectedFlux (the evidence is STATUS 2026-08-31 and `default.parameter`).
- `STATUS.md:6` "Last updated 2026-09-10" (the file runs to 2026-09-28); STATUS 470-501, 520-535, 775-844
  get inline "superseded by" markers; STATUS:3305 "Open" is closed by b1798c3.
- `cases/default.parameter:57-76` still calls the clip "REQUIRED on polyhedra" (refuted); its lines
  817-869 (alphaFTimeLevel, donorPlaneAdvected, projectMassFlux, +5-8 % inflation) predate the
  closed-box fix and are VOID but unmarked; add them to the VOID list in STATUS §0.
- `printMethodBanner.H:139` prints "alg_lin (default)" while the code default is `geo_lin`; the
  `baselineEulerian` candidate header still says the Eulerian solver freezes rho; `workflow/README.md`
  says seven method libraries (eight). One-line fixes.
- `docs/IMPROVEMENTS.md` links three briefs that do not exist; the nPSL note's status line is stale;
  `plan-library-split` says "nothing executed" (it was). Pointers, not rewrites.
- Gate blind spots found by the red-team memo: `workflow/scripts/make_gate_summary.py:190-198` sets the
  exact-1D closed form to `None` for every candidate (criterion "target" vacuous for candidates), and
  the band-gradient metric is centred. Recorded in the report §2, the KB `method-gates` note and a
  STATUS §11 line as OPEN gate repairs; the code fixes (E0.1 closed forms, the checkerboard-sensitive
  diagnostic, the curvature column, the zero-flow arm) are the first items of the next campaign, not
  part of this plan. `config/candidates/HL0.yaml` gets a CORRECTED comment: its pre-registered read-out
  ("zero normal strain on the interface, so the band error falls") conflated the interface with the
  band.
Everything else stays where it is; the knowledge base links to it.

---

## Appendix A. Raw material: the decision and finding log extracted from the record (agent, 2026-09-28)

(Verbatim from the extraction agent. Tags: S = STATUS.md, M = METHOD.md, C = CLAUDE.md, I = docs/
IMPROVEMENTS.md, PCS = plan-curvature-stabilization.md, PSH = plan-shannon-parasitic-currents.md, PCT =
plan-combined-source-terms.md, PHL = plan-halo-limited-gradient-control.md, PLS = plan-library-split,
RM = capillary-level-set-research-roadmap.md, DP = cases/default.parameter, G2 = methodGate2D.yaml,
T2 = stationaryDroplet2D fvSolution.template; citations are tag:line-range.)

### STATUS.md section list
- 14–73 §0 READ THIS FIRST — the translating droplet case was a CLOSED BOX (2026-09-02); 74–129 The source is the curvature estimator; translation and density ratio are amplifiers; 130–189 ANCHORED: Popinet's benchmark reproduced (2026-09-04); 190–228 DECIDED: the instability does not survive a perfect capillary force; 229–247 What SURVIVES the retraction; 248–259 The curvature chain is exonerated (RETRACTED 2026-09-27); 260–269 Leading open defect (RETRACTED 2026-09-27); 270–289 New, runtime-selectable, default-off.
- 290–350 §1 What is being worked on; 351–364 §2 Where the writing lives; 365–388 §3 Where the numbers live; 389–551 §4 State of the measurements, with subsections: 552 Mesh alignment EXONERATED (08-26); 580 gradU contamination CLOSED for the 2D vortex (08-27); 619 INVALIDATION: parallel kinematic SL baselines predate the gradU fix (08-26); 718 Domain size is an axis, 6R certified (08-18); 775 The wide ladders (08-19); 846 K is EXONERATED (08-19); 892 INVALIDATION: filtered results predate a psi-filter seam bug (08-19); 924 The 3D "instability under refinement" may not be an h effect; 949 Two further defects fixed; 981 The two-factor law (08-20); 1034 Two process failures; 1050 Full-horizon stability gate and t_blow RETRACTION (08-31); 1103 ANSWERED: projectedFlux win is the RECONSTRUCT OPERATOR (08-31); 1166 RUNNING: refinement; 1197 DECIDED: hanging nodes (09-04); 1260 DECIDED: static refinement equivalent (09-04); 1309 RUNNING: static local refinement; 1514 RUNNING: Popinet 3D poly (09-05); 1738 RESULT and VOID: tilted wall faces (09-05); 1846 RESULT: poly ladder diverges late, inlet/outlet transport defect (09-05); 2008 DECIDED: quasi-monotone clip (09-08); 2041 DECIDED: SI units (09-08); 2121 MEASURED: the amplification bound of the fit (09-08); 2184 CONFIRMED: not a time-step effect (09-09); 2233 DECIDED: polyhedral cells, not cfMesh (09-09); 2312 G4 FALSIFIES the candidate (09-09); 2389 RETRACTED AND REPLACED: the clip's damage (09-09); 2497 capillary time step on polyhedra audited (09-08).
- 2525 §5 Lichtenberg running; 2581 §6 how to run; 2643 §7 Next; 2690 §8 sync; 2703 §9 SDPLS thread sync (09-22, 9.1–9.7); 2885 §10 build policy (10.1–10.5); 3218 §11 gradient control (11.1–11.16, 3228–4041).

### A. Interface advection and phase indicator
- A1 quadratic SL (production): the update is an assignment (M:52–67); departure-centred AB2 foot (M:69–96; why: arrival form leaves +dt^2 d_t u, 2–4 % early, 35–47 % by t=0.02 in the oscillating droplet); constant-free quadratic WLS fit, w=1/|d|, Cholesky (M:98–121; value fit and degree >= 2 are load-bearing: Taylor variants diverge, a linear value fit drives ||grad psi|-1|| to 1e21); CPC on hex, CFC on poly (M:123–125; RM:1390); orders 2.97/2.59, 2.95/3.28, 1.36/1.46, 20x vs Eulerian at half the wall clock (M:343–347), contaminated 08-26, re-established 08-27 for the 2D vortex 2.84/3.30 (PCS:52–61); SL_RECONSTRUCTION uncached per case (M:373); SL_FIT normalEquations: QR blows up identically, pivot 0.757 at Lambda 1.2608 (M:377; S:2269–2278); quadraticPivotTol 0.3 (09-05→08): cfMesh boundary-layer cells cond 3e7–6e12, tolerance history 0/1e-3/1e-2 diverge at steps 8/12/16, bit-inert on hex (S:1550–1722); SL_STENCIL_BOUNDARY_FACES include, inflowOnly selectable (S:1874–1990); 1/3/6 non-orth correctors inert (S:1941–1957); SI units (S:2041–2081); translation none saturates at N=256 (M:648–653, OPEN); the transport operator amplifies on every mesh, rho(B)=1.00441 hex / 1.01028 pMesh, growth per physical time ~1/h (M:402–484), acceptance rho<=1 met nowhere (M:767–776).
- A2 linear SL: coupled transport-order axis flat (S:503–518, filtered → subject to the filter seam bug); linearTaylor grows under refinement kinematically (S:514–518); quadraticTaylor needs the clip (RM:22–28).
- A3 nSL: strain mode anti-convergent (order −1.24; PCT:78–95); geometric write-back diverges x1.7 per 10 steps, dt-independent, delta d ~ psi delta|g|/|g|^2 (PCT:167–187); normalProjection trajectory latest runaway but physically invalid by t=0.03, "do not promote" (RM:136–155, 302–327, 525–534).
- A4 Eulerian SDPLS transport: benchVortex N=512: band gradient noSource 1.46e-4, R 7.87e-1, beta 5.16e-1 (PCT:47–61); beta+explicit diverges at every CFL (PCT:66–71); R over-flattens near the resolution limit (PCT:72–76); R amplifies solver-tolerance noise 7 orders over 130 steps, seam gate needs psi tol 1e-14 (S:2836–2856); Rdiv N=128 diverged (DP:186–189).
- A5 redistancing retired (PCT:599–604; PCS:378–381; M:777–790; RM:31–32, 536–550): plane-based band rewrite O(h^2 kappa) one-signed compounding; central-Hamiltonian PDE redistancer divergent (DP:495–497); frozen-band gated reset never run (PCS:302–344, 531–538, 1585–1587).
- A6 clips/value bounds: sigma=0 control fails on the poly Popinet mesh at step 508 (S:1822–1844); global clip removes the far-field failure (S:2008–2039); band-aware clip RETRACTED (S:2389–2495): clip fires in six stencil-extremum cells, extremum exemption cuts damage 30x, apex tie; G4 falsifies (S:2312–2387: 59.2 % of cells to bound are extrema); SL_CLIP false, limiters collapse the order (BJ 3.0→0.1, Venk 3.0→0.9) (M:378–380, 777–790); slValueBound family, sentinel fromClipSwitch, bit-identity 8 arms x 1563 steps (S:280–286; DP:94–158); lipschitzCone: exact gate passes, N=64 gain RETRACTED the same day by the advection ladder (2.3–3.6x worse translation, 19x vortex, 9x 3D shear), converged ladders (8.6x→189.7x), resolution ladder (floor 7.1e-3, centroid +119 %) (M:488–729); the 2D orders of METHOD 8.3.7 were 3/2 too high (S:3237–3242), CSVs to regenerate OPEN.
- A7 projectedFlux: fullHorizonStability2D: off+projectedFlux best, −52 vs +118 1/s (S:1050–1088); the win is the reconstruct operator 70 %, extension 0 %, solenoidality 30 % (S:1103–1164); projFlux ladders recorded only in DP:1029–1035; default since c935883 (M:375); under a frozen uniform stream identical to cellCentred (M:437–455); traceFlux physical|extension added 09-26 (S:3300–3303); SAAMPLE warning: reconstruct diverges at a gradient jump (PSH:791–805); traceTranslating3D disturbance columns VOID (S:3716–3723).
- A8 foot integrator taylor, no gate (M:376); kernel order 3.00 verified (PCS:1364–1372); rk2 inside the scatter on the translating droplet (S:3736).
- A9 psi outer correctors: frozen-force lag exonerated (PCS:1352–1363); gain3D +0.4 %/−0.1 % (PSH:207–236); default yes 08-28 on consistency grounds, m=2 rate inert (DP:403–449; M:314–317).
- A10 Lambda and rho(B): Lambda_c=|1−Σg|+Σ|g|, hex 1.0527, pMesh 1.2608, a smaller step cannot help (S:2121–2182); dt sweep confirms within 6.5 % (S:2184–2231); polyhedral cells, not cfMesh or layers (S:2233–2310); Lambda is not a proxy for rho (M:423–435).
- A11 gradU contamination (08-26): setVelocity wrote face values into processor patches; 31 parallel kinematic studies contaminated (S:619–640); closed for the 2D vortex (S:580–617); regression bit-identical (M:680–690).
- A12 phase indicator detrixheAslam, NOT INERT via c935883 (DP:5–9; M:129–140); the SDF assumption survives in the indicator's first-order offset, unmeasured (PCS:863–878, OPEN); determinant guard failure at N=512 fixed (RM:1339–1352); DA volume drift ~6 % at N=32 under translation (RM:170–177); algebraic psi in oscillatingDroplet2D → DROPLET_SURFACE token (S:3271–3279, OPEN void); band metrics use the unlimited gradPsiMetric (S:3261–3269).

### B. Viscosity
- VISCOSITY_FACE_MODEL alg_lin decided 09-03 on the 36-arm mufGrid2D ladder (DP:701–763; M:397): only model with positive order in all three metrics at both jumps; L1 1.400e-3, p 1.10 at ratio 1000. Grid {alg,geo}x{lin,harm,blend}, blend weight |n.Sf|/|Sf| parameter-free; harmonic and blend worse everywhere (DP:761–762). RETRACTED switch to geo_lin: muf was frozen for the algebraic arms (bug, 39e59b3). Popinet's face properties are alg_lin (S:139–143). Viscosity damps in the Laplace sweep (S:164–169). Eulerian solver froze muf, fixed c094bd8. OPEN: 3D templates have no token. Flag: mufGrid2D ran at np 8 before the 09-27 parallel fixes.

### C. Surface tension
- C1 balanced-force flux, constant kappa absorbed (~3e-11) (M:243–272); exact kappa cuts the step-1 kick 6 orders (S:97–102); amplifierGate: exact-kappa arms bounded, "the curvature estimator is the whole story" (S:190–228); projection converged, absorbs 99.85–99.98 %, 30x tighter solve <=2.4e-4 (PCS:546–602); hanging nodes do not break the balance (S:1197–1258); skewed meshes: exact kappa leaves 1.39e-3, levers do not compose, suspect reconstruct (RM:1028–1295); interface-mean replay quiet → spatial variation drives (RM:724–777); pressure jump over pure phases, the alpha=1/2 partition biased 3 % (M:321–337); Young–Laplace t=0 solve fixed (S:960–966).
- C2 estimators/deliveries: arithmetic h^1.13, ellipse 0.97, G h^2 0.64; psiOverGradPsi offset correction 0.477 % second order, superseded by cellCentreInverse; stabilizedFootPointFace h^2.04, ellipse 1.98, the only delivery meeting the criterion, production until c935883; cutCell retired (ellipse 1.02, blows sooner); cellMean retired (lowest gain, prefactor win); symmetric face mean 1.10; footPointEvaluatedFace coupled verdict falsified (blows earlier, volume 5–15x worse; score on r(A2h)); closestPointNewton normal extension of kappa dead end (40–85x amplification); harmonicLaplace dead end 2; height function 2D only; connectedInterface Helmholtz halves m=2, fails coupled; FVM div(grad psi/|grad psi|) h^1.16; interFoam operator swap runs away (PSH:658–678); cellCentreInverse +2.01 non-gradient content, default since c935883, stationary residual 4.60/3.81/1.55/1.49x lower, translating undecided, 93.5 % unfilled-cell inversion bug fixed, second order on constant curvature only; Kang 58x better statics but earlier blow-up; sharpHeaviside "58x/0.07 s" (T2:215–217 only; attribution flag); acceptance criterion G h^2 <= 0.65 and ellipse order >= 1.9 (S:2685–2688).
- C3 K-aware inverse: sphere h^1.95 vs h^1.02 (RM:1362–1396); K off → 962 vs 3.2, two arms diverge (S:846–890); needs a parallel foliation not an SDF; torus test open.
- C4 semi-implicit force: dropped 08-19, revived, deadlocked on the cluster 08-28 (collectives in a master guard), +10.5 vs +118 but the winner uses none (S:1080–1096); untested on the translating droplet after the fix (OPEN). endStep centring, midpoint spectrally identical, midpoint 32 % earlier divergence (S:3737).
- C5 integral surface tension: 2D static equilibrium at N=64, blows at N=128; "integral reformulations smooth corrupted geometry, they do not create missing signal" (PCS:371–373); translating runaway ~1.4e-3 s (RM:157–256); conormal prototype 49.45 % residual (RM:411–442).
- C6 psi filter biharmonicBand: delay device (PCS:90–97); "combination works" WITHDRAWN 08-19 (S:892–922 seam defects f83a1ab); damping flips sign with resolution (S:986–1014); no filtering in production (C:382–398).
- C7 variational force: the worked example (C:403–407); interFoam bounded by operator–transport–state pairing, not accuracy (PSH:585–656); kappaClamp and freezeKappa: the pump is the live curvature refill (PSH:807–889); no implementation.
- C8 mechanism: max|U|(T)=u0(h) exp G(h), u0 = curvature error independent of U0, G grows with U0 and density ratio (S:74–125); two-factor law u0~h^3.5, G~h^−3.27 (3D) (S:981–1005; PSH:244–324); r·dt gain changes sign with resolution (S:924–947); G h^2 orders deliveries but is incomplete (PCS:731–767); loop model gamma ~ 0.647 sigma dt^2/(rho h^3), r=r0+c dt (PCS:1300–1372); the unstable mode is m=2, the instability is the explicit capillary coupling (PCS:1485–1592); Popinet reframing: the m=8 capillary wave anti-damped +100..+250 vs −128 physical, net work per cycle of a mesh-locked curvature error (PSH:398–583); corrugation baseline 0.209h at 16h displacement (PCT:144–149); psi side supplies ~21 %, band renormalisation 3.4x worse (PCS:1155–1298); mesh alignment exonerated (S:552–578); 6R box equivalent to 10R, 4R fails (S:718–773); 3D poly ladder currents grow after decay while geometry converges (S:1474–1512); old t_blow RETRACTED (S:1050–1068); Shannon target g <~ 4e-5 vs 4.4e-4 (PSH:1003–1021).

### D. Mass flux
- D1 closed-box VOID (S:14–65; C:661–705; 440107f): VOID lists S:55–62, 3566–3572.
- D2 rhoLENT (Liu 2023 eq. 40) vs geometricFaceDensity vs interpolatedDensity: history c935883 → 82ca995; stationary evidence stands (+1.0/−22/+0.1 %), translating VOID; dec002f FALSIFIED mass–momentum consistency as the dominant term (9 orders in the residual, <2x in the excess; the density ratio is the mechanism); geometricFaceDensity residual 0.04–0.56 harmless on the stationary droplet; interpolatedDensity unmeasured.
- D3 boundRho true since 28a1383, ACTIVE with rhoLENT; basis VOID → OPEN; rhoClipFraction counts round-off.
- D4 donorPlane by instruction, unmeasured; averagedPlanes VOID/diverged; donorPlaneAdvected ~1 %, trapezoid marginal, projectMassFlux 47x on the reducible part but stops the droplet — all pre-fix (VOID by the rule, unmarked: flag 8).
- D5 ddt pairing backward/backward by the matching argument; tables VOID; flip-flop history f7307b5/28a1383/b60e3df; MOMENTUM_DIV upwind inert on the stationary droplet.
- D6 density ratio is the mechanism (kickOriginGate2D; ratio 1 completes; keep rho1+rho2 fixed; Popinet at ratio 1 4.86 % of U).
- D7 late translating instability after the fix: horizon cut to 0.05 s (stated as post-hoc), not a decomposition artefact, outlet triggers, 20 mm box diverges by 0.21 s, start-position test refuted the box-centre reading, 40 mm box interior degradation shrinks with h (two rungs) — OPEN.
- D8 coupled-face density defect 28d13f0 (rho_f 90 % apart; 1e-5..5e-4 → 1e-8..1e-10); which earlier parallel studies to re-run OPEN; RM:510–511 had flagged it in July.
- D9 Eulerian frozen rho and the port c094bd8 (travel 1.00012 vs 0.420); void the five Eulerian studies OPEN; curvature dispatch and droplet CSV remain.

### E. Gradient control
- E1 SDPLS continuum diagnosis (PCT:97–115), noise gain <= 1 + lambda dt pi W, design window lambda dt <~ 1/(pi W) (PCT:117–165); combined = gradientControl linearQ + full, SUBSUMED (PCT:3–7); WP0 delta_h diagnostic never run; dead ends PCT:594–634.
- E2 laws: eight + strain weights; boundedGradient not a class; 109 unit checks; SW_GAMMA = artanh 0.9; GC_EPS placeholder until E0 (S:3373–3427).
- E3 slSource: exponential update, clamp 30; q and n from fvc::grad instead of the geometry fit (deviation from plan); 49 checks, halo mutant fails on 4 ranks (S:3415–3435).
- E4 sdplsGradientControl equals R bit for bit (189 checks).
- E5 extensions: none identity; roadmap N=32 translating survival none 0.0103, meshWave 0.0131, closestPoint 0.0258, steadyUpwind 0.0388 s (RM:121–134); UEXT_DIV upwind (deferred correction diverges on the steady equation, DP:16–18); extensions worst accuracy at 12–27x cost for pure advection (PCT:623–625); closestPoint decomposition-dependent seam FAIL 0.45/1.11; haloLimited D4 (S:3437–3470); HL0 FAIL, first suspect stencilFit (S:3972–3990); HL0 at R=2h not run (S:3487–3491); steadyUpwindLinear, pseudoTime, anisotropicDiffusion: no record.
- E6 extension points Phase C: 53 cases bit-identical, 250/250 renders (S:3318–3401).
- E7 gate campaign: pre-fix all FAIL on the shear arm alone; fixed gate kinematic arms byte-identical (177 pairs); scoring corrections; baseline oscillating unstable at N=200; final verdicts (FP0 0.939/3/8/5/5, HL0 0.828/16/13/0/0, HL1q 0.204/22/18/0/0, HL1z 6.490/23/19/0/0, HL2 0.754/7/5/0/3, S1 0.707/12/5/1/1); coupled seam S1 0.76 FAIL, FP0 1.11 FAIL, HL2 n.c.; HL1q 5x better gradient, 173x worse shape; 3D gate does not run.

### F. Verification methodology
- Research loop C:409–508; method gates C:603–642, PHL:448–595; mesh-convergence rule C:510–537; regression set C:539–601; L_inf rule S:196–198, PHL:554–560; seams and 4 ranks C:199–246 with the defect list (gradU 08-26, updateFlux, psi filter f83a1ab, band dilation, one-sided kappa_f at seams, semiImplicit collectives, 28d13f0, b1798c3); wrong-setup voids (closed box, tilted wall faces, stale arms, vacuous snappy arm, algebraic psi, 3D translating reference velocity; rule C:661–705); retraction protocol C:707–718; bit-identity C:720–729, S:3010–3113; measurement traps (log classifier C:263–296; empty squeue S:1034–1048; absolute thresholds S:2226–2231; a finished 0/ C:778–788; per-step counters S:2475–2481; shallow config merge M:25–29; mis-set thresholds S:1220–1246; t_blow as proxy PCS:1311–1313); cancellation-dominated solver convergence C:731–759; provenance C:761–788.

### G. OPEN questions and author decisions
- PHL:695–704: (1) curvature extension for the translating and oscillating arms; (2) void the algebraic-psi studies; (3) void the frozen-rho studies. PCT:670–703 seven open questions. PCS:121, 422–434. RM:118–119 gates 4–5. S:3175 WP6 deletion 2026-10-01; S:3241 regenerate CSVs; S:3584 boundRho; S:3587 3D viscosity token; S:3671 which parallel studies to re-run; S:3782 longer box token; S:3968 oscillating horizon; S:1101 metric; S:2306–2310; S:2383–2387; S:2519–2521; S:2661–2668 indicator offset; S:3487–3491 HL0 at R=2h; S:3620–3622/3862 untried levers; S:3946 N=200 rung; S:3988–3990 HL0 discriminator; S:3703–3710 Eulerian steps; S:3257 D-c; S:3713–3714 diagnostics; M:648–653; M:677–678; M:752–799; C:748–750; RM:1271–1287.

### H. Audit flags (13)
1. G2:69 cites VOID rhoDdtGate2D for boundRho; G2:66/DP:668 quote a pre-fix 4.8e-2 residual. 2. M:375 wrong gate for projectedFlux. 3. METHOD §8 table and §9 items 1/5/6 describe 2026-07-31; item 6 contradicts M:234–239. 4. M:35 names libleiaLevelSet. 5. S:6 "Last updated 2026-09-10". 6. S:470–501, 520–535, 775–844 lack withdrawal markers. 7. DP:57–76 calls the clip REQUIRED on polyhedra. 8. DP:817–869 studies of 2026-09-02 morning are VOID by the rule but unmarked. 9. mufGrid2D at np 8 predates the parallel fixes. 10. PHL:698–700 proposes none for translating; gate uses cellCentreInverse. 11. "58x/0.07 s" attributed to both Kang and sharpHeaviside. 12. RM:101–102, 525–527 superseded July rules. 13. S:3305 "Open" closed by b1798c3 but not marked.

## Appendix B. Raw material: the docs inventory (agent, 2026-09-28; paths under docs/)

Gaps: IMPROVEMENTS.md links three briefs that do not exist (`semi-lagrangian-level-set/improvement-interface-offset.md`, `velocity-extension/improvement-metric-footpoint.md`, `normal-projected-semi-lagrangian/improvement-drift-gate.md`); the velocity-extension article (88 lines), the GRL article (advected sections TODO) and the GCLS deck (5 slides, no images) are skeletons; built decks and PDFs are git-ignored.

1. **semi-lagrangian-level-set** — article `sl-level-set-article/semiLagrangianLevelSet.tex` (2555 lines), title "A high-order, matrix-free semi-Lagrangian level-set method for incompressible two-phase flow on arbitrary polyhedral meshes". Sections: 122 Intro; 213 Method (216 SL update, 231 second-order foot, 304 WLS value reconstruction with ¶418 why a value fit and ¶429 admissibility, 488 mesh-adaptive stencil, 542 boundedness, 550 cache-free, 621 phase indicator); 715 Two-phase coupling (723 NS, 759 viscous term ¶794/¶830/¶875, 919 surface tension from reconstructed curvature, 985 static refinement ¶1000/1038/1075/1090); 1101 Software design; 1161 Verification (1164 setup, 1200 2D convergence, 1250 3D hex/poly, 1331 static refinement measured ¶1379 WB gate ¶1457 production ladder ¶1493 poly, 1531 phase-indicator accuracy, 1556 alpha evolution, 1577 interface evolution, 1603 curvature accuracy, 1658 stationary droplet ¶1754 stabilisation study ¶1780 delivery study ¶1852 cost, 1901 translating droplet ¶1931 free stream ¶1944 curvature error source ¶1991 amplifier ¶2015 mass–momentum not the mechanism ¶2030 net propulsion ¶2058 exact force, 2107 Popinet comparison, 2200 flux-space residual); 2275 Outlook variational force (2282, 2304, 2358, 2395); 2448 Limitations; 2508 Conclusions; 2543 code and data. Key numbers: ¶830 alg_lin only model with positive orders (L1 3.49e-3 at N=512, p 0.77/1.10); 1200 shape 2.97/2.59, gradient order ~0.5; 1250 shear 2.95/3.28, deformation 1.36/1.46; 1331 refined 51 640 vs 216 000 cells, Delta p 145.470 vs 145.48, orders agree to <=0.09, 3.4–6.9x fewer core-hours (13.4x two levels; poly 2.7x); 1531 indicator order ~2.0, geometric = DA to 8 digits; 1603 curvature O(h^1.2) 35 %→1 %; 1658 max|u| drops 46x N=32→64, t_blow 0.47/0.44/0.11/0.03 s, CST floor 5e-7 at N=64 diverges at N=128 t~0.047, "better static balance, higher dynamic gain"; 1901 one-step kick independent of U0 to 0.6 %, exact kappa removes it 2e4–7e4x, kErr 7 % at R/h 12.8, drift 13 % of U0 by 60 ms, exact-kappa plateau 3e-5 vs 3.0e-2; 2107 Popinet Linf 4.86 % vs ~5 %, orders 0.49/0.91/1.71; 2200 projection absorbs 99.85–99.98 %. Supplementary `supplementary/iDEC_failure_report.tex` (iDEC diverges: 39→1.4e4→4.7e9; rho>1). Decks: `quadratic-semi-lagrangian-level-set.template.html` (117 sections: roadmap; level set & phase indicator; SL advection incl. QWLS/stencil/UQWLSR/linear line/2D+3D convergence; two-phase flow incl. segregated step, Rhie–Chow, balanced-force theorem, curvature estimation/delivery, parallel-surface inverse, cut-cell blow-up 3.3x sooner, stationary droplet, static refinement, parasitic-current instability, lecture series slides, fields; conclusions; roadmap capillary coupling incl. solver residual, projectedFlux; software design) and `...-negative-results.template.html` (28 sections: time discretization; curvature post-processing; alternative estimators; force delivery Kang/CST/meta-law; transport & reinit (mass–momentum flux, plane-anchored redistancing preliminary); process traps (band filter, PL pipeline partial positive, methodology traps, index of failed levers)). Data: 63 tables (droplet_parasitic, curvature_error, capillary_flux_residual, refined_*, popinetTranslating, stationaryDroplet3D*, sl_fit_pivot_census, ...), 20 figures, 5 mechanism CSVs.
2. **sdpls-level-set** — `sdpls-article/sdplsLevelSet.tex` (3389 lines) "The SDPLS level-set method for two-phase flows: source-term control of the signed-distance property without interface displacement". Sections: 124 Method (125 source, 152 linearization, 167 defect correction); 176 Two defects (sign, use-after-free, wrong operator metric, one discretization); 232 Verification (¶241 which instant, ¶299 measurement, ¶350 3D shear, ¶370 deformation, ¶394 result, ¶406 where the order is not, 441 error estimate, 684 where the error lives); 842 Two-way coupling (865 setup, 920 production curvature destroys the droplet, 1041 time discretization not the driver, 1127 drain seeded by curvature error, 1263 Rdiv and a retraction, 1791 difference of two upwind operators, 2076 two invariants); 2162 Two models in response (2169 exponentialImplicit, 2507 global volume correction); 2678 Two coupled-patch defects (2690 face-centre value, 2776 raw flux, 2805 consequences); 2862 Negative results; 3172 Conclusions. Key numbers: R band gradient order +0.74 vs −0.26 (31x) on the non-reversing vortex; 3D shear +0.668 vs −0.094 but R 3.7x worse shape, 18x worse volume; order ceiling ~1.2 under second-order transport; coupled: R drains N=32 (−1.000), diverges N=64/128, mode-4 current amplified 260x, exact kappa 0.000 in 12/12; Rdiv diverges kinematically (−3.8/−1.3/−1.5), withdrawn; exponentialImplicit shape +1.603 vs +1.360; volume correction is a crossover; np=8 seam still 2.8e-3; "closestPoint 2185x worse" DISPUTED (serial re-measure 0.172/0.267); beta's target structurally wrong (band mean 1.516→1.481 flat). Decks: `sdpls-level-set.template.html` (65 sections) and `...-negative-results.template.html` (36 sections: why beta fails, retracted R near resolution limit, kinematic gates are blind, Rdiv, temporal corrections, volume correction crossover, coupled-patch defects, dead `limit` entry). Note `improvement-interface-defects.md` (interfaceDefects function object brief, pending). Data: 28 tables, 6 figures.
3. **method-comparison** — `method-comparison-article/methodComparison.tex` (354 lines) "Which level-set advection method, given sufficient resolution?": 93 convergence (T=8 filament under-resolved at N<=128), 115 cost vs accuracy (SL 20x more accurate at half the wall clock at 512², T=8), 124 decision table (SL 1.49e-5 at 595 s; Eulerian 3.04e-4 at 1170 s; VE closestPoint 1.44e-3 at 1.61e4 s), 147 flux-form volume loss (−17 % by t=6.5; band |grad psi| 0.72 vs 0.96), 258 measured SL improvements (limiter 3.0→0.9, QR bit-identical, flux-form conserves ∫psi to 2e-14 but 3–13x volume error), 304 frozen-band redistancing injurious (E_vol 0.017→1.77), 194 verdict (SL quadratic; Eulerian robust second; SDPLS only for coarse-mesh volume; VE dominated). Deck `level-set-method-comparison.template.html` (79 sections incl. a long "Capillary balance" track: one enforceable flux contract, only exact curvature balanced, five models fail before 4 ms, spatial curvature variation is the defect, connected face curvature, mode gates, oracles, pressure gates; and "Time centring": amplification matrix). Data: 133 tables (+ VOID_closedBox_20260902/README), 84 figures, 1 animation; `REPRODUCE.md`.
4. **gradient-controlled-level-set** — `gcls-level-set-article/gclsLevelSet.tex` (988 lines; sections as written 2026-09-27: 141 transport and q; 170 source terms; 258 extension; 294 discretisation incl. 397 consistency under decomposition; 477 gate; 614 results incl. 789 translating at long times; 846 discussion; 958 appendix data); deck skeleton (5 slides); `figures/` TikZ cases + make_case_figures.py, make_archive.py, make_result_figures.py; data/tables (41), data/figures (5 PDFs), data/archive/shared-method-config-2026-09-01-192-g1150e68 (README, MANIFEST, gate/, prefix/, laptop/); README.md (theme, dossiers not in git).
5. **linear-semi-lagrangian-level-set** — `lsl-level-set-article/linearSemiLagrangianLevelSet.tex` (679 lines): nestedLSQ shape ~1.1, volume ~1.5, band gradient ~2.1 at CFL ½; CFL 1 destabilises at N>=128; linearTaylor gradient defect O(1e9); 3D deformation 2.28/1.52, hex shear stable only to N=50; UQWLSR ~3 vs linear ~1.1. Deck 34 sections. Data 10 tables, 11 figures.
6. **geometrically-redistanced-levelset** — `grl-level-set-article/geometricallyRedistancedLevelSet.tex` (484 lines, draft): planeFootWave band Linf 4.48e-5 at h=1/256 (order 2.0), anchoredEikonal 3.53e-4, PDE reinit increases error; one step preserves a plane to 1e-14, circle volume change O(h^2) 1.5 %→0.02 %; advected sections TODO; negative results: foot-cloud scalloping, anchored-Eikonal first order, PDE reinit "stability bomb". Decks: main (23 sections) and negative results (29 sections: scalloping, Eikonal fill, PDE reinit two failures, measurement pitfalls, scoreboard). Data 10 tables, 6 figures.
7. **velocity-extension** — article stub (TODO sections, placeholder figure). Deck `velocity-extension.template.html` (162 sections): models (none, anisotropicDiffusion, pseudoTime, steadyUpwind, steadyUpwindLinear divergent cascade, closestPoint, meshWave) each with motivation/model/discretization/code/verdict; static verification (n·grad Uext = 0, static convergence, why legacy models cannot converge, anchoring, one-sided errors); advected verification (reversed vortex convergence per model, volume loss at large T, SD quality, normal-constancy channel, reversibility bias, non-reversing steady vortex "none never wins", winding & annihilation); conclusions "which extension, by problem"; appendix atlases (alpha and |grad psi|−1 for 6 models x 4 T). Data 60 figures, no tables.
8. **combined-source-terms** — `levelset_combined_source_note.tex` (1355 lines, 2026-08-06; analysis only: F = a + lambda(1−|grad phi|) gives the logistic interfacial law and an invariant interval); `improvement-sdpls-combined.md` SUBSUMED 2026-09-26.
9. **normal-projected-semi-lagrangian** — `npsl-article/normalProjectedSemiLagrangian.tex` (505 lines): trace clean (0.017h one step), write-back diverges (1.107h after 100 steps, x1.7/10 steps), corrugation hypothesis falsified (0.209h vs 0.223h), deforming flow orders 2.9 vs −1.24/−0.17, status: value path most accurate; deck 12 flat slides; notes `normal-projected-semi-lagrangian.md` (521 lines, stale status line) and `stable-foot-point-3d.md` (185 lines).
10. **docs/api** Doxygen setup (Doxyfile INPUT=../../src). **Top-level:** IMPROVEMENTS.md (briefs index, 2026-08-25); capillary-level-set-research-roadmap.md (1414 lines, gates 0–5, dated records 07-27..08-07); gradU-coupled-patch-contamination.md (234 lines, handover 08-26); plan-combined-source-terms.md (SUBSUMED); plan-curvature-stabilization.md (1592 lines, v0.2, §6–18 measured, contamination notice); plan-halo-limited-gradient-control.md (approved 09-26); plan-library-split-and-build-policy.md (stale "nothing executed"); plan-shannon-parasitic-currents.md (1100 lines, §0–0g evidence, Pólya §1–6, execution phases); snakemake-openfoam-profiles.md (499 lines). **Build:** `docs/build-decks.sh` (propagate_data.py, export_html.py inlines figures + CDN assets for every `docs/*/*-presentation/*.template.html`, 11 templates; `--linear`); Makefile: ART_SL/LSL/GRL/GCLS, `decks`, `article-sl` (latexmk *.tex), `article-sdpls` (not in `articles`), `article-lsl`, `article-grl`, `article-gcls`, `articles`, `docs = decks articles`, `comparison`, `sl-quadratic`, `sl-linear`, `gate`; no target for velocityExtension.tex, normalProjectedSemiLagrangian.tex, the combined note, the iDEC supplementary; methodComparison.tex only through `make comparison`.

## Appendix C. Raw material: the code map (agent, 2026-09-28)

Discrepancies found: `printMethodBanner.H:139` prints "alg_lin (default)" while `leiaLevelSetTwoPhaseFoam/createFields.H:92` defaults to `geo_lin`; base classes of phaseIndicator, profile, surfaceTensionForce, velocityModel declare TypeName("none") but never register it (velocityModel::New defaults to "none" → fatal without a type); `config/candidates/baselineEulerian.yaml` header still says rho/rhoPhi frozen (fixed by c094bd8; "writes no droplet metrics" still true); `workflow/README.md` says seven method libraries (eight); code default phaseIndicator `geometric` vs token `detrixheAslam`; slReconstruction code default `quadraticWeightedLeastSquares` vs gate `uncached...`.

**Libraries** (src/leiaLevelSet/*/Make/files): libleiaCore (phaseIndicator, narrowBand, profile, velocityModel, leiaVersionRegistry; `schemes/levelSetBlended` not compiled), libleiaAdvection, libleiaGradientControl, libleiaRedistancer, libleiaSdplsSource, libleiaSemiLagrangian, libleiaSurfaceTension (incl. the semi-implicit fvOption), libleiaVelocityExtension, libleiaVolumeCorrection; plus liblevelSetImplicitSurfaces, libleiafiniteVolume, libleiaFunctionObjects. LINKS (etc/leia-check-deps.py): GradientControl {}, SdplsSource {GradientControl}, SemiLagrangian {GradientControl}, Redistancer/VolumeCorrection/SurfaceTension {}, VelocityExtension {SemiLagrangian}, Advection {SdplsSource, SemiLagrangian, VelocityExtension}.

**Families → models (dictionary key):** phaseIndicator (`levelSet.phaseIndicator.type`, default geometric): heaviside, sharpJump, geometric, detrixheAslam (geometrySource levelSetField|analyticImplicitSurface). narrowBand (`levelSet.narrowBand.type`, default none): none, empty, signChange, neighbours, distance, phaseIndicator. profile (`levelSet.profile.type`): signedDistance, tanh. velocityModel (top-level `velocityModel.type`): shear2D, deformation3D, shear3D, translation, rotation, vortex2D, periodic2D, uniaxialStrain (+ fluxCorrection helper). levelSetAdvection (`levelSet.advection.type`, default eulerian): eulerian (nDefCorr, builds velocityExtension + sdplsSource), semiLagrangian. sdplsSource (`levelSet.sdplsSource.type`, default noSource): noSource, R, beta, Rdiv, RdivStrictSp, gradientControl; strategies discretization {none, explicit, simpleLinearImplicit, strictNegativeSpLinearImplicit, exponential, exponentialImplicit}, gradPsi {fvc, narrowLS}, mollifier {none, m1, band}. gradientControlLaw (`...law.type`, required): none, linearQ, linearZ, cubicQ, cubicZ, twoThirdsZReg, saturatedLinearZ, softWall; strainWeight {none, full, omega}. slReconstruction (`levelSet.semiLagrangian.reconstruction`, default quadraticWeightedLeastSquares; keys stencil point|face, stencilBoundaryFaces include|exclude|inflowOnly, slopeLimiter, clipToStencilBounds, clipRegion, clipKeepExtrema; `levelSet.geometryFit.type` second fit): linearTaylor, linearWeightedLeastSquares, signedDistanceLinearWeightedLeastSquares, quadraticTaylor, quadraticWeightedLeastSquares (alias quadraticWLSQ), uncachedQuadraticWeightedLeastSquares, signedDistanceQuadraticWeightedLeastSquares, bandQuadraticWeightedLeastSquares, defectCorrectedIDW. slValueBound (`...valueBound`, sentinel fromClipSwitch): none, stencilBounds ("falsified as a fix"), lipschitzCone (lipschitzMode unity|stencil, onInadmissible cellOnly|none). slCorrector (`...correction`, default direct): direct, deferredCorrection. slScheme (`...scheme`, default pointValue): pointValue (footIntegrator taylor|rk2, trajectoryVelocity input|normalProjection|normalClosestPoint), fluxForm, normalProjected (renormalization geometric|strain|none, offsetEngine rayRoot|footPoint). slSource (`...source.type`, default none): none, gradientControl (bandCells, law). velocityExtension (`levelSet.velocityExtension.type`, default none): none, anisotropicDiffusion, pseudoTime, steadyUpwind, steadyUpwindLinear, closestPoint, meshWave, haloLimited (radiusCells <= 1; strategies travel capped, weight fractionReached, direction levelSet, sampler stencilFit); interfaceExtension intermediate base (nLayers, nDescent, nAnchorLayers, projectFlux, fadeMode, maxScale). redistancer (`levelSet.redistancer.type`, default noRedistancing; trigger interval|gradPsiThreshold|signedDistanceBounds): noRedistancing, PDE, anchoredEikonal, planeFootWave. volumeCorrection (default noVolumeCorrection): noVolumeCorrection, newtonShift (alias globalShift). surfaceTensionForce (`levelSet.surfaceTensionForce.type`, required; faceCurvatureSource model|registered): constantCurvatureSurfaceTension, constantCurvaturePressurePotential, curvaturePressurePotential, correctionKang, divGradAlphaSnGradAlpha, divGradPsiSnGradAlpha, traceGradGeoNormalSnGradAlpha, traceGradGradPsiSnGradAlpha, integralSurfaceTension, integralConormalSurfaceTension, isoCurvature, reconstructedCurvature. fv::option semiImplicitCapillaryForce (form off|value|increment, coeff, laplaceBeltrami). implicitSurface (`levelSet.implicitSurface.type`): implicitPlane, hesseNormalPlane, implicitSlab, implicitSphere, slottedSphere, implicitEllipsoid, signedDistanceEllipsoid, signedDistanceEllipse, implicitSinc. Others: gradScheme noBvGrad; functionObjects gradPsiError, gradPsiErrorCSV, psiConservationCSV, writeIsoSurfaceTopo. Not RTS: the psi filter (`levelSet.psiFilter.type none|biharmonicBand`, theta).

**Solvers:** leiaLevelSetFoam (kinematic unified; redistancer, narrowBand, phaseIndicator, velocityModel, levelSetAdvection, volumeCorrection); leiaRedistancedLevelSetFoam (Eulerian advective + criterion-gated redistancing); leiaSemiLagrangeLevelSetFoam (kinematic SL; traceVelocity cellCentred|projectedFlux, traceFlux physical|extension; velocityExtension only with extension); leiaLevelSetTwoPhaseFoam (Eulerian two-phase, interIsoFoam-based; shares createMassFluxFields.H/updateFaceDensity.H/updateMassFlux.H since c094bd8; viscosityFaceModel default geo_lin; psiOuterCorrectors; writes leiaLevelSetFoam.csv); leiaSemiLagrangianLevelSetTwoPhaseFoam (SL two-phase; reuses the sibling's createFields/UEqn/pEqn/YoungLaplaceEqn; createSLFields builds slAdvection; createTransportFields builds slVelExt + mass flux). Non-RTS switches: massFlux {type interpolatedDensity|geometricFaceDensity|rhoLENT, alphaFSource donorPlane|averagedPlanes|donorPlaneAdvected (legacy upwind/central), boundRho, massResidualDiagnostic, alphaFTimeLevel new|trapezoid, projectMassFlux}; viscosityFaceModel {alg_lin, alg_harm, alg_blend, geo_lin, geo_harm, geo_blend} (legacy interpolated/geometric/geometricHarmonic); curvatureExtension.type {none → slAdvection::meanCurvatureNoExtension, closestPointNewton → meanCurvatureClosestPoint, fvm, harmonicLaplace → harmonicMeanCurvatureExtensionEqn.H, footPointHeightFunction → footPointCurvature.H, connectedInterface → connectedInterfaceCurvature.H (fitHalfWidth, tangentialRegularization, estimator connectedFit|analyticImplicitSurface|rdfQuadratic), interfaceMean → interfaceMeanCurvature.H, cellCentreInverse → cellCentreInverseCurvature.H (gaussianCurvature), stabilizedFootPointFace | cutCellFootPointFace | cellMeanFootPointFace | symmetricFaceMeanFootPointFace (faceSmoothing) | footPointEvaluatedFace | cellFootPointEvaluatedFace → stabilizedFootPointFaceCurvature.H, all six require offsetCorrection none; relax}; capillaryForceCentring endStep|midpoint; PIMPLE.psiOuterCorrectors; traceVelocity cellCentred|projectedFlux|reconstructedU|reconstructedMomentum; footIntegrator taylor|rk2 (pointValueScheme.C:63); psiRenormalization none|footPointGradient; psiFilter; curvatureClamp; dropletReferenceVelocity; env SL_FREEZE_KAPPA, SL_FREEZE_RHOPHI, SL_SMOOTH_RHO; fvOptions tokens SEMI_IMPLICIT_CAPILLARY(+_COEFF, _LB).

**Tests (applications/test):** leiaTestConnectedCurvatureModes, CurvatureNoiseGain, DeparturePoint, FoamGeometry, FoliationResidual, GradScheme, GradientControlLaw, HaloLimited (cases/haloLimitedUnit), LevelSet, MeanCurvature, ParallelSurfaceInverse, Redistance, SLReconstruction, SdplsSource (cases/sdplsSourceUnit), SignedDistanceEllipsoid, SlSource (cases/slSourceUnit), TransportSpectrum, VelocityExtension.

**Workflow:** study config keys (study_name, case, mesh hex|perturbed|poly|hexRefined|polyRefined, mode, np, solver, theme, axes_override, collapse_other_axes, setfields_args/app, post_solve, solve_runtime, env_preamble, mpi_launcher, export_slides); token precedence axes_override > cases/<case>.parameter values{} (…_poly.parameter for poly) > cases/default.parameter; derived tokens in materialize._with_derived_tokens; gate files (config/gates/methodGate2D|3D.yaml: solvers, lineTokens, methodTokens, lineOf, rates, twoPhaseCoupling, verdict; arms exact1D, shear(+seam), stationary, translating(+translatingSeamNp1), oscillating; candidates baseline, baselineEulerian, FP0, HL0, HL1q, HL1z, HL2, S1); Snakefile.gate rules all/run_arm/exact1d_check/summarize/compare; scripts render_gate_configs.py, make_gate_summary.py, make_gate_tables.py, foam_log_state.sh (COMPLETED 0/RUNNING 1/STALLED 2/DIVERGED 3/LAUNCH_FAILURE 4/MISSING 5), compare_metrics_csv.py, guard_finished_cases.py, aggregate.py, materialize.py, richardson.py; other Snakefiles: comparison, curvature, sl-linear, sl-quadratic, pressure-compatibility.

**Cases:** kinematic 1Dstretch, 1DredistanceTest, 2DgradTest, 2DredistanceCircle, 2DredistanceStatic, 2Dtranslation, 2Dvortex, 2Dcontactline-{periodic,vortex,vortex-oscillate}, 3Dcontactline, 3Ddeformation, 3Drotation, 3Dshear, 3Dtranslation; static curvature ellipseDroplet2D, ellipsoidDroplet3D; coupled stationaryDroplet2D/3D, translatingDroplet2D/3D, oscillatingDroplet2D/3D, popinetTranslating2D/3D, isoStaticDroplet2D, linWLSdroplet2D, qadvFvmDroplet2D, oscISO/ISTDroplet2D, transISO/ISTDroplet2D; references interFoamDroplet2D/3D, interFlowDroplet2D/3D; units haloLimitedUnit, sdplsSourceUnit, slSourceUnit. **Configs (352)** grouped: transport ladders (advConv*, bench*, bulkVortex*, *Conv*, nsl/npsl, uncached*, trace*, ...), stationary 2D/3D droplet studies, curvature gates, mass-flux gates (massFluxComparison2D, matchedBDF2Translating2D, rhoBoundGate2D, rhoDdtGate2D, rhoLENTStationary2D, projFlux*, momentumDivScheme*, volumeCorrection*, wellBalancedTranslating2D), viscosity (viscGate2D_*, mufDecide2D_*, mufGrid2D_*, viscousHorizon2D), translating (translating*, resLadder2D_*, kickOriginGate2D, amplifierGate*, bestConfig*, transIST*), oscillating, popinet*, seam checks, gates, sdpls*.

## Appendix D. Verbatim artefacts of the vault design (design agent, 2026-09-28; measured against Quartz v5.0.0 sources)

Facts: `github.com/leia-openfoam/leia` is public (HTTP 200); Pages serves `https://leia-openfoam.github.io/leia/` (stale Doxygen; README line 9 badge and line 97 point at it); `.github/workflows/docs.yml` is dead (`doc/Doxygen/` does not exist). `.gitignore` ignores `docs/**/*.html` (templates excepted), `*.pdf`, `*.csv`, `*.png`, `log`, `log.*`, `[0-9]/`, `build/`. Decks: reveal `hash: true`, sections without `id`, vertical stacks → anchors `#/h` or `#/h/v`. Every article section has a `\label`. Quartz v5.0.0: Node >= 22, npm >= 10.9.2; one `quartz.config.yaml` with `configuration:`, `plugins:`, `layout:`; plugins in `.quartz/plugins/` installed by `npx quartz plugin install` (from `quartz.lock.json`) and `npx quartz plugin resolve` (config entries the lockfile lacks); no `--from-config` flag in v5.0.0; non-Markdown files under `content/` are emitted at their path; `folder/index.md` titles a folder; frontmatter keys `title, description, tags, aliases, permalink, date/created, modified, published, publish, draft, enableToc`; CrawlLinks `markdownLinkResolution` must equal Obsidian's link format. Laptop: Node 18.19.1, no nvm, Docker Desktop on WSL2, python 3.13, latexmk.

### D1. `.github/workflows/knowledge-base.yml` (replaces docs.yml, deleted in the same commit)

```yaml
name: Knowledge base
on:
  push:
    branches: [main, development]
    paths: ["docs/**", "STATUS.md", "METHOD.md", "CLAUDE.md", ".github/workflows/knowledge-base.yml"]
  workflow_dispatch:
permissions: {contents: read, pages: write, id-token: write}
concurrency: {group: pages, cancel-in-progress: false}
env: {QUARTZ_TAG: v5.0.0, KB_DIR: docs/knowledge-base}
jobs:
  build:
    runs-on: ubuntu-24.04
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with: {node-version: 22}
      - uses: actions/setup-python@v5
        with: {python-version: "3.12"}
      - name: Check the knowledge base (frontmatter, wikilinks, log coverage)
        run: python3 "$KB_DIR/.quartz/check_kb.py" "$KB_DIR"
      - name: Build the 3D graph data
        run: python3 "$KB_DIR/.quartz/build_graph.py" "$KB_DIR"
      - name: Clone Quartz at the pinned tag
        run: git clone --depth 1 --branch "$QUARTZ_TAG" https://github.com/jackyzha0/quartz.git quartz
      - name: Install Quartz and its plugins
        working-directory: quartz
        run: |
          npm ci
          npx quartz plugin install
          npx quartz plugin resolve
      - name: Stage the content (single source: docs/knowledge-base)
        run: |
          rm -rf quartz/content && mkdir -p quartz/content
          rsync -a --exclude ".obsidian" --exclude ".quartz" --exclude "templates" "$KB_DIR"/ quartz/content/
          cp "$KB_DIR/.quartz/quartz.config.yaml" quartz/quartz.config.yaml
      - name: Stage the record files for full-text search (links still point at GitHub)
        run: |
          mkdir -p quartz/content/record
          cp STATUS.md METHOD.md CLAUDE.md quartz/content/record/
          cp docs/plan-*.md docs/capillary-level-set-research-roadmap.md docs/IMPROVEMENTS.md quartz/content/record/
      - name: Build the reveal decks (python3 only) into content/decks
        run: |
          bash docs/build-decks.sh
          mkdir -p quartz/content/decks
          for f in docs/*/*-presentation/*.html; do case "$f" in *.template.html) continue ;; esac; cp "$f" quartz/content/decks/; done
      - name: Compile the pre-prints into content/preprints (a LaTeX failure does not block the site)
        continue-on-error: true
        run: |
          sudo apt-get update -q
          sudo apt-get install -y --no-install-recommends latexmk texlive-latex-base texlive-latex-recommended \
            texlive-latex-extra texlive-publishers texlive-science texlive-pictures texlive-fonts-recommended cm-super
          mkdir -p quartz/content/preprints
          for tex in docs/*/*-article/*.tex; do
            dir=$(dirname "$tex"); name=$(basename "$tex" .tex)
            if (cd "$dir" && latexmk -pdf -interaction=nonstopmode -halt-on-error "$name.tex"); then cp "$dir/$name.pdf" quartz/content/preprints/; else echo "[skip] $tex"; fi
          done
      - name: Build the site
        working-directory: quartz
        run: npx quartz build -d content -o public
      - uses: actions/upload-pages-artifact@v3
        with: {path: quartz/public}
  deploy:
    needs: build
    runs-on: ubuntu-latest
    environment: {name: github-pages, url: "${{ steps.deployment.outputs.page_url }}"}
    steps:
      - id: deployment
        uses: actions/deploy-pages@v4
```
(The technical report's PDF is compiled by the same loop when its folder is named `*-report/*.tex`: add `docs/*/*-report/*.tex` to the glob, output into `content/preprints/`.)

### D2. `docs/knowledge-base/.quartz/quartz.config.yaml` (derived from quartz.config.default.yaml at v5.0.0)

```yaml
configuration:
  pageTitle: leia knowledge base
  pageTitleSuffix: " - leia"
  enableSPA: true
  enablePopovers: true
  analytics: null
  locale: en-GB
  baseUrl: leia-openfoam.github.io/leia
  ignorePatterns: [private, templates, .obsidian, .quartz]
  defaultDateType: created
  theme:
    fontOrigin: googleFonts
    cdnCaching: true
    typography: {header: Schibsted Grotesk, body: Source Sans Pro, code: IBM Plex Mono}
    colors:
      lightMode: {light: "#faf8f8", lightgray: "#e5e5e5", gray: "#b8b8b8", darkgray: "#4e4e4e", dark: "#2b2b2b", secondary: "#284b63", tertiary: "#84a59d", highlight: "rgba(143, 159, 169, 0.15)", textHighlight: "#fff23688"}
      darkMode: {light: "#161618", lightgray: "#393639", gray: "#646464", darkgray: "#d4d4d4", dark: "#ebebec", secondary: "#7b97aa", tertiary: "#84a59d", highlight: "rgba(143, 159, 169, 0.15)", textHighlight: "#b3aa0288"}
plugins:
  - {source: "github:quartz-community/note-properties", enabled: true, order: 5, options: {includeAll: false, includedProperties: [status, part, date_settled, decided_by, tags, aliases], excludedProperties: [], hidePropertiesView: false, delimiters: "---", language: yaml}, layout: {position: beforeBody, priority: 15, display: all}}
  - {source: "github:quartz-community/created-modified-date", enabled: true, order: 10, options: {priority: [frontmatter]}}
  - {source: "github:quartz-community/syntax-highlighting", enabled: true, order: 20, options: {theme: {light: github-light, dark: github-dark}, keepBackground: false}}
  - {source: "github:quartz-community/spacer", enabled: true, order: 25, options: {}, layout: {position: left, priority: 25, display: mobile-only}}
  - {source: "github:quartz-community/obsidian-flavored-markdown", enabled: true, order: 30, options: {enableInHtmlEmbed: false, enableCheckbox: true}}
  - {source: "github:quartz-community/github-flavored-markdown", enabled: true, order: 40}
  - {source: "github:quartz-community/table-of-contents", enabled: true, order: 50}
  - {source: "github:quartz-community/crawl-links", enabled: true, order: 60, options: {markdownLinkResolution: absolute, prettyLinks: true, openLinksInNewTab: true, externalLinkIcon: true}}
  - {source: "github:quartz-community/description", enabled: true, order: 70}
  - {source: "github:quartz-community/latex", enabled: true, order: 80, options: {renderEngine: katex}}
  - {source: "github:quartz-community/remove-draft", enabled: true}
  - {source: "github:quartz-community/alias-redirects", enabled: true}
  - {source: "github:quartz-community/content-index", enabled: true, options: {enableSiteMap: true, enableRSS: false}}
  - {source: "github:quartz-community/favicon", enabled: true}
  - {source: "github:quartz-community/content-page", enabled: true}
  - {source: "github:quartz-community/folder-page", enabled: true}
  - {source: "github:quartz-community/tag-page", enabled: true}
  - {source: "github:quartz-community/explorer", enabled: true, options: {title: Knowledge base, folderClickBehavior: collapse, folderDefaultState: open, useSavedState: true}, layout: {position: left, priority: 50}}
  - {source: "github:quartz-community/graph", enabled: true, options: {localGraph: {drag: true, zoom: true, depth: 1, scale: 1.1, repelForce: 0.5, centerForce: 0.3, linkDistance: 30, fontSize: 0.6, opacityScale: 1, removeTags: [], showTags: true, enableRadial: false}, globalGraph: {drag: true, zoom: true, depth: -1, scale: 0.9, repelForce: 0.5, centerForce: 0.3, linkDistance: 30, fontSize: 0.6, opacityScale: 1, removeTags: [], showTags: true, focusOnHover: true, enableRadial: true}}, layout: {position: right, priority: 10}}
  - {source: "github:quartz-community/search", enabled: true, layout: {position: left, priority: 20, group: toolbar, groupOptions: {grow: true}}}
  - {source: "github:quartz-community/backlinks", enabled: true, layout: {position: right, priority: 30}}
  - {source: "github:quartz-community/article-title", enabled: true, layout: {position: beforeBody, priority: 10}}
  - {source: "github:quartz-community/content-meta", enabled: true, layout: {position: beforeBody, priority: 20}}
  - {source: "github:quartz-community/tag-list", enabled: true, layout: {position: beforeBody, priority: 30}}
  - {source: "github:quartz-community/page-title", enabled: true, layout: {position: left, priority: 10}}
  - {source: "github:quartz-community/darkmode", enabled: true, layout: {position: left, priority: 30, group: toolbar}}
  - {source: "github:quartz-community/reader-mode", enabled: true, layout: {position: left, priority: 35, group: toolbar}}
  - {source: "github:quartz-community/breadcrumbs", enabled: true, layout: {position: beforeBody, priority: 5, condition: not-index}}
  - {source: "github:quartz-community/footer", enabled: true, options: {links: {GitHub: "https://github.com/leia-openfoam/leia", STATUS.md: "https://github.com/leia-openfoam/leia/blob/development/STATUS.md", METHOD.md: "https://github.com/leia-openfoam/leia/blob/development/METHOD.md"}}}
  - {source: "github:quartz-community/stacked-pages", enabled: true, layout: {position: afterBody, priority: 50, display: all}}
layout:
  groups: {toolbar: {priority: 35, direction: row, gap: 0.5rem}}
  byPageType:
    "404": {positions: {beforeBody: [], left: [], right: []}}
    content: {}
    folder: {exclude: [reader-mode], positions: {right: []}}
    tag: {exclude: [reader-mode], positions: {right: []}}
    canvas: {}
    bases: {}
```
(Written in block style in the repo; the flow style above is for compactness. Omitted plugins: citations, hard-line-breaks, ox-hugo, roam, explicit-publish, og-image, cname, canvas-page, bases-page, comments, recent-notes, encrypted-pages. Re-read `quartz.config.default.yaml` at the pinned tag before the first build and adjust option names if the tag differs.)

### D3. `.quartz/build.sh` and the Makefile targets

```makefile
kb:        ## knowledge base site (docs/knowledge-base -> build/kb/public); Quartz v5 needs Node >= 22
	bash docs/knowledge-base/.quartz/build.sh
kb-serve:
	bash docs/knowledge-base/.quartz/build.sh --serve
kb-graph:  ## regenerate graph3d/graph.json and print the local preview command
	python3 docs/knowledge-base/.quartz/check_kb.py docs/knowledge-base && python3 docs/knowledge-base/.quartz/build_graph.py docs/knowledge-base
	@echo "python3 -m http.server -d docs/knowledge-base 8000   # then open http://localhost:8000/graph3d/?local=1"
```
```bash
#!/usr/bin/env bash
# Build the knowledge base with Quartz v5 into build/kb/public. Same steps as the workflow.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
kb="$repo/docs/knowledge-base"; out="$repo/build/kb"; q="$out/quartz"; tag="${QUARTZ_TAG:-v5.0.0}"
python3 "$kb/.quartz/check_kb.py" "$kb"; python3 "$kb/.quartz/build_graph.py" "$kb"
mkdir -p "$out"; [ -d "$q/.git" ] || git clone --depth 1 --branch "$tag" https://github.com/jackyzha0/quartz.git "$q"
rm -rf "$q/content"; mkdir -p "$q/content"
rsync -a --exclude .obsidian --exclude .quartz --exclude templates "$kb"/ "$q/content/"
cp "$kb/.quartz/quartz.config.yaml" "$q/quartz.config.yaml"
bash "$repo/docs/build-decks.sh" && mkdir -p "$q/content/decks" \
  && find "$repo/docs" -path '*-presentation/*.html' ! -name '*.template.html' -exec cp {} "$q/content/decks/" \;
serve=""; [ "${1:-}" = "--serve" ] && serve="--serve"
node_ok() { command -v node >/dev/null && [ "$(node -p 'process.versions.node.split(".")[0]')" -ge 22 ]; }
if ! node_ok && [ -s "$HOME/.nvm/nvm.sh" ]; then . "$HOME/.nvm/nvm.sh"; nvm use 22 >/dev/null || nvm install 22; fi
if node_ok; then (cd "$q" && npm ci && npx quartz plugin install && npx quartz plugin resolve && npx quartz build -d content -o public $serve)
else echo "[kb] node >= 22 not found: building in node:22-slim (Docker Desktop)"
  docker run --rm -it -v "$q:/work" -w /work -p 8080:8080 -p 3001:3001 node:22-slim \
    bash -c "npm ci && npx quartz plugin install && npx quartz plugin resolve && npx quartz build -d content -o public $serve"; fi
echo "[kb] site: $q/public"
```
nvm install (on request, once): `curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash && nvm install 22`. The system Node 18 is not touched.

### D4. `.quartz/check_kb.py` (stdlib; exit 1 on any failure, prints `file: message`)
1. every `*.md` outside `templates/` has frontmatter with `title`, `kind` in {hub, concept, model, decision, retraction, case, study, session, index}, `status` in {settled, open, retracted, voided, candidate}, `part` in the six parts or `all`, `tags` containing the kind; 2. every `[[target]]`, `[[target#h]]`, `[[target|alias]]`, `![[target]]` resolves to `<target>.md` (or an alias, or a file under `assets/`) relative to the vault root; 3. every note in `decisions/` appears in `decision-log.md` and every note in `retractions/` in `retraction-log.md`; 4. every settled decision has `date_settled` and non-empty `decided_by`; 5. no file name matches `log.*` and no single-digit folder; 6. `## Log` is the last `##` section of every note. `build_graph.py` (section 1c) reuses the same parser and fails on an unresolved link.

### D5. `.obsidian/` (committed) and `.gitignore`
`app.json`: `{"newLinkFormat": "absolute", "useMarkdownLinks": false, "alwaysUpdateLinks": true, "attachmentFolderPath": "assets", "showUnsupportedFiles": true, "strictLineBreaks": true}`; `templates.json`: `{"folder": "templates", "dateFormat": "YYYY-MM-DD"}`; `core-plugins.json`: the orga vault's file with canvas, daily-notes, bases, sync false; `appearance.json`: `{}`; `graph.json`: the orga file plus `colorGroups` for `tag:#hub` (2899536), `#decision` (3066993), `#retraction` (15158332), `#model` (15844367), `#case` (10181046), `#study` (3447003). `.gitignore` additions: `docs/knowledge-base/.obsidian/workspace.json`, `.../workspace-mobile.json`, `.../cache`, `.../plugins/*/data.json`, `docs/knowledge-base/graph3d/graph.json`, and `!docs/knowledge-base/assets/*.png`, `!docs/knowledge-base/assets/*.svg`.

### D6. Templates (the exact Markdown lives in `templates/`; body outlines)
- `templates/note.md` (concept | model | case | study): frontmatter (title, aliases, tags [kind, part/<part>], kind, status, part, date, date_settled, decided_by, code, sources); `# title`; `> one-paragraph verdict with date`; `## What it is`; `## Why it matters`; `## Where in the code` (pinned blob links); `## Evidence` (table claim | number | where); `## Why it failed, or why we think so` (only when status is retracted/voided or open with a measured failure); `## Decisions`; `## Open questions`; `## Related`; `## Log` (`### YYYY-MM-DD` entries).
- `templates/decision.md`: title "TOKEN: <setting in one line>", aliases [TOKEN]; `> TOKEN <value> in the <layer> since YYYY-MM-DD. Decided by config/<gate>.yaml: <number>`; `## The question`; `## The measurement that decided it` (table arm | metric | value | where; the pre-registered read-out link); `## What it does not cover`; `## Related`; `## Log` with "SETTLED ... Entered in [[decision-log#YYYY-MM]]".
- `templates/retraction.md`: title = the claim as stated; `> RETRACTED YYYY-MM-DD. The claim was "...". The measurement is .... Scope of the void: ...`; `## The claim, and where it lived`; `## Why it was wrong, or why we think so`; `## What survives`; `## Propagation (checklist, same commit)` (STATUS, METHOD, plan, deck/table, every KB note, the log line); `## Related`; `## Log`.
- `templates/hub.md`: `[[index]] <- back`; `## The question this part answers`; `## Current verdict (date)`; `## Map` (table kind | note | status | one line); `## Open, in order`; `## Log`.
- `templates/session.md`: `> Branch, sha, stamp; last handover; read [[decision-log#YYYY-MM]] and [[retraction-log#YYYY-MM]] from that date`; `## What is being worked on`; `## What is open, in order`; `## Traps that cost time`; `## Where the numbers live`; `## Log`.

### D7. Log line formats
decision-log: `- 2026-09-09 · [[decisions/sl-clip-and-value-bound-off]] · SETTLED \`SL_CLIP false\`: G4 fails at step 506 against 527 unclipped; +30.4 % volume error on 2D hex translating · \`config/popinet3D_poly_sigma0_clipGate.yaml\``. retraction-log: `- 2026-09-02 · [[retractions/closed-box-translating-droplet]] · VOIDS "rhoLENT + matched momentum scheme makes the translating droplet run" and every translatingDroplet2D result before 440107f: the mesh had no inlet and no outlet · propagated to STATUS 0, METHOD 8.1 rows ...`.

### D8. Implementation order of the vault (each step one commit)
1. Skeleton: folders and index files, conventions.md, empty logs with month headings, five templates, `.obsidian/`, `.quartz/{quartz.config.yaml, build.sh, check_kb.py, build_graph.py, README.md}`, `graph3d/index.html`, the workflow, the Makefile targets, `.gitignore` lines, README badge and sentence, delete `docs.yml`. Verify: Obsidian opens the vault; `make kb` builds (Docker or nvm); a `workflow_dispatch` deploys. 2. Six hubs and `sessions/current.md` from STATUS §0, §1, §7 and METHOD §9. 3. Sixteen decision notes from METHOD §8.1 and the decision log. 4. Twelve retraction notes and the retraction log. 5. Models, concepts, cases, studies (Appendix A–C are the raw material). 6. The CLAUDE.md/AGENTS.md sections; pointer lines at the top of STATUS.md and METHOD.md.
