# VOID: every curated table of cases/2Dtranslation before 2026-09-29

`cases/2Dtranslation` was meant to translate a circle one way, left to right, at U = (1 0 0).
It was made from `cases/2Dvortex` (commit cd97e6da, 2026-09-01) without the vortex's line
`oscillation @!OSCILLATION!@;`, so the velocity model's default (`oscillation` on, `tau` =
`endTime`, `src/leiaLevelSet/velocityModel/velocityModel.C` lines 49-50) multiplied U by
cos(pi t / tau). The circle moved at most tau/pi = 0.16 and returned to its start at T, and every
error metric compared T with the initial fields. Confirmed from the recorded runs: U =
(-0.99991 0 0) at t = 0.4979, the centroid back at x = 0.2502, and step counts that follow the
cosine law exactly.

By the rule "A wrong setup voids its data" (CLAUDE.md), every number these tables carry is void,
including the ones that look unaffected:

| file | study | used in |
|---|---|---|
| `kinematicTranslation2D_errors.csv` | `config/kinematicTranslation2D.yaml` (18 arms) | the config's RESULT block; the projectedFlux null control |
| `value_bound_advection_2Dtranslation.csv` | `config/coneBoundMesh2Dtranslation.yaml` + `workflow/scripts/advect_bound_arm.sh` | METHOD.md 8.3.4, the 2Dtranslation row and finding 1 |
| `advConv2Dtranslation_convergence.csv` | `config/advConv2Dtranslation.yaml` | METHOD.md 8.3.7 (the translation part), 8.3.8, the hex 2D regression rung |

The study trees are renamed `*_VOID_reversedTranslation_20260929` (the laptop's
`~/OpenFOAM/repos/leia/studies/` and Lichtenberg's `leia-curvature/studies/`). The fixed case
(oscillation off, an exact end reference, psi zeroGradient on all four patches) re-runs every
study; the new tables replace these at their old paths. (CORRECTED 2026-09-29: this README first
said "exact boundary data"; exact psi values on the outflow and on the inflow patch made the SL
update unstable, and those two re-runs are void too, STATUS.md 11.19 item 3.) STATUS.md 11.19, the knowledge-base note
`retractions/reversed-2dtranslation`.
