# The PDE level-set direction (D. Bothe and K. Soga) on the Eulerian FVM line

Started 2026-11-12 by the author's decision. The model:

    d(phi)/dt + u_e . grad(phi) = beta phi g(|grad phi|),   u_e(x) = u(Pi(x))

a velocity extension (Bothe and Soga, J. Differential Equations 468, 2026) plus a gradient source
(Bothe, Fricke and Soga, Math. Ann. 390, 2024), in `leiaLevelSetFoam` and
`leiaLevelSetTwoPhaseFoam`. The semi-Lagrangian line is out of scope here.

| what | where |
|---|---|
| the plan (sub-problems 0, A0, A to E, each with an exact reference) | `docs/plan-velocity-extension-verification.md` |
| the progress deck for the supervisor (editable source) | `pde-level-set-presentation/pde-level-set.template.html` |
| build the deck (self-contained HTML, git-ignored) | `python3 workflow/scripts/export_html.py docs/pde-level-set/pde-level-set-presentation pde-level-set.template.html` |
| the record of results | `STATUS.md` from section 11.29 |
| the knowledge-base note | `docs/knowledge-base/concepts/pde-level-set-verification.md` |

The dossiers and D. Bothe's handwritten notes of 2026-10-29 stay outside git (`agent-input/`); the
deck carries a transcription of the notes and their reading.
