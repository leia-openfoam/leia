## Pieces path: serial performance fixes F1 to F5 (diag2, 2026-11-18/19)

Branch `davof-ccs-pieces-perf` off `davof-ccs-pieces` 240f3974, all changes in `src/davof/geometry/davofSweptPieces.H` (sweptPieceExactVolume) and the
call sites in `davofSweptAdvectionFlux.C` (one member `davof::pieceWorkspace pieceWs_` of `davofSweptAdvection`, declared in `davofSweptAdvection.H`).
Nothing merged. Builds: one worktree per commit under `agents/diag2/perf/<tag>` (-O3, own platforms/).

### Set-up facts

- Step times are the mean timeStepTime of the first 100 steps of deformation3D N 64, OF Co 0.5 (the case of the 2026-11-16 profile), 1 rank,
  an exclusive node of the type of that profile (feature `i01`, Intel Xeon Platinum 9242, 96 cores; the profile baseline 3.39 s ran on `mpsc0345`),
  three repetitions, median of the three repetition means. A first G4 series on `i02` nodes (Xeon Platinum 8470Q, 1.6 s for the baseline) was
  discarded as the wrong node type (kept in the logs `logs/d2-g4-*.out`).
- The default constructors of FixedList and Vector do not initialise (`= default`); the 7.5 % of the profile at the clipAll entry was therefore not a zero
  initialisation but the entry of the per-polygon loop with the polygon temporary and the copy (`q = out`); after F1 that category is 1.1 to 2.0 %.
- G3 compares every printed per-step quantity of the solver log (alpha extremes before and after the bounding, mass changes, sum(alpha V), bounded-face
  counters, Courant numbers, the final areaFractionError) and E_shape, E_mass of `advectorVoF.csv` with the build of 240f3974: shear2D N 64 (1331
  steps), shear3D N 32 (288), deformation3D N 32 (480), 1 rank, OF Co 0.5, the cases `diag2_c_*` index 00001. E_shape of the base build is the stored
  value of the 112c284a rows (1.4764013347392773e-04, 6.0140199357620163e-03, 9.3147214126802838e-03).
- G1: `davofTestPieceClip` (the two dumped tets) compared line by line with the base build (40 lines).
- G2: the four `diag2_gate_*` studies (80 pieces cases), thresholds of `kiss_an.py ccsplane`, run for one build after the other.
- A first attempt at G2 for F1 and F2 at once, and the first G3 runs of deformation3D N 32, hit "Disk quota exceeded" (the F1 studies alone were 13 GB at that
  moment); those runs were discarded, the field and mesh directories of the invalid studies and of the 2026-11-15 gate set were removed (files created by
  this task and by the earlier diag2 gate), the gates were rerun one build at a time and cleaned after the analysis.

### The table

| fix | commit | step time (s, median of 3) | vs 3.339 s | G1 | G3 (printed digits; E_shape relative difference) | G2 (80 cases) |
|---|---|---|---|---|---|---|
| base | 240f3974 | 3.339 (3.339, 3.370, 3.337) | 0 | reference | reference | 80 of 80 (2026-11-15, 112c284a) |
| F1 reused workspace | 25601230 | 2.549 (2.537, 2.560, 2.549) | -23.7 % | 40 of 40 lines identical | all printed lines identical, E_shape 0 on the three cases | 80 of 80 P1 ok; maxima identical to 112c284a (L1 1.076e-12, Linf 4.03e-11, E_shape 1.55e-22, min alpha -3.69e-13); np1/np8 pair (oblique N 160 Co 1) differs 7.04e-3 as in 112c284a |
| F2 classification before clipping | 8e7fc2dc | 1.885 (1.882, 1.885, 1.888) | -43.5 % | identical | identical, E_shape 0 | 80 of 80 P1 ok; maxima L1 1.078e-12, Linf 4.24e-11, E_shape 1.55e-22, min alpha -5.61e-13; the np1/np8 pair differs 4.89e-3 (7.04e-3 before) |
| F3 vertex distances once | 3cea1cab | 1.792 (1.794, 1.792, 1.791) | -46.3 % | identical | identical, E_shape 0 | 80 of 80 P1 ok; same maxima and pair as F2 |
| F4 pseudo-angle order (reverted 4cfee6c9) | 27742435 | 1.552 (1.552, 1.555, 1.551) | -53.5 % | identical | FAILS identity: first differing step 2 (the CCS residual of the davof line, 2.6228e-15 against 2.6229e-15; the step-1 alpha extremes are identical); E_shape relative 2.6e-2 shear2D N 64, 1.3e-12 shear3D N 32, 3.2e-9 deformation3D N 32; shear3D step 4 sets 8965 faces against 9019 | not run |
| F5 no cap for the PLIC clip, on F4 (reverted with F4) | dee2bc8f | 1.565 (1.565, 1.568, 1.564) | -53.1 % | the dumped PLIC volumes are exactly 0 (2 of 40 lines differ; 3.5e-28 before) | first differing step 1 (alpha extremes at 1e-15); E_shape relative 1.9e-13, 1.4e-13, 3.0e-7 | not run |
| F5 alone, on F3 (build 4cfee6c9; reverted 04ea39a1) | dee2bc8f after 3cea1cab | 1.747 (1.747, 1.752, 1.747) | -47.7 % | PLIC volumes exactly 0 (2 of 40 lines) | FAILS the 1e-8 prediction: first differing step 1; E_shape relative 7.4e-3 shear2D N 64, 1.7e-12 shear3D N 32, 3.0e-8 deformation3D N 32 | not run |
| final (HEAD) | 04ea39a1 | the sources equal F3 (`git diff 3cea1cab HEAD` empty): 1.79 s | -46.3 % | identical | identical | as F3 |

The reverts are commits (4cfee6c9, 04ea39a1); F4 and F5 stay in the history and can be cherry-picked. The perturbation experiment of F4, F5 and F4+F5 shows that E_shape
of shear2D N 64 changes by 2.6e-2, 7.4e-3 and 1.9e-13 for three different last-bit changes of the cap arithmetic: the case is sensitive to round-off through the
bounding and limiter decisions (face counts differ at step 4 of shear3D N 32), so a gate of identity to the digit is the only discriminating one for these three cases.

### Profile shares (HPCToolkit CPUTIME, serial, deformation3D N 64, OF Co 0.5, 100 steps; `perf/prof/<tag>/an.*`)

| item | baseline 240f3974 (3.42 s/step, CPUTIME) | after F2 (8e7fc2dc; profiled on an i02 node) | final = F3 tree (1.80 s/step, CPUTIME) |
|---|---|---|---|
| sweptPieceExactVolume, inclusive | 82.3 % | not read | 66.0 % |
| clipOne (distances, crossing, copy into out), category | 40.0 % | 25.5 % | 25.5 % |
| polygon copies (std::copy of FixedList<point, 24>, exclusive) | 15.3 % | about 0.5 % | about 1.3 % (memcpy only) |
| clipAll per-polygon loop (temporary, `q = out`; the "zero-init" line) | 4.5 % as a category (7.5 % on the entry line) | 1.1 % | 2.0 % |
| cap (frame, angle, sort, copy and orientation) | 13.3 % | 22.3 % | 21.5 % |
| of which atan2 | 7.1 % | 13.0 % | 12.4 % |
| classification (clipPlane, F2) | n/a | 5.4 % | 5.1 % |
| areaOf (all sites) | 3.7 % | 4.6 % | 6.4 % |
| area clips (sweptPieceArea) | 1.7 % | 3.0 % | 3.5 % |
| function objects | 4.6 % | 9.8 % | 9.0 % |

The shares are of each run's own total, so the cap and the other untouched parts grow as the step shrinks (the cap's own time is 45.5 s of 342 s in the
baseline and 38.7 s of 180 s at the end). The largest item of the final profile is still sweptPieceExactVolume (66.0 %); inside it the clip (clipOne, 25.5 %) is first, the cap
second (21.5 %, atan2 12.4 %).

### 64 ranks, deformation3D N 128, OF Co 0.5, first 200 steps, once (exclusive `i01` node, production builds, the decomposition of the 2026-11-16 profile)

| build | mean step (s) | median step (s) |
|---|---|---|
| base 240f3974 | 5.081 | 5.110 |
| final (sources of 3cea1cab) | 2.572 | 2.565 |

The factor is 1.98 (serial N 64: 1.86). No decomposition work was done; the step still follows the most loaded rank.

### Pre-registered predictions (facts; the verdicts are Tomislav's)

- PF1 (step falls 15 to 22 % after F1): 23.7 %.
- PF2 (a further 25 % or more of the baseline after F2): 19.9 % of the baseline (2.549 to 1.885 s).
- PF3 (at most 2.0 s after F1 to F5, falsifier above 2.4 s): 1.79 s with F1 to F3 kept (F4 and F5 reverted); 1.55 s with F4 kept.
- PF4 (sweptPieceExactVolume still the largest item): yes, 66.0 %.

### Jobs

Builds 56591392 (base), 56591393 (F1), 56591394 (F2), 56591395 (profiled F2), 56690826 (F3), 56690827 (F4), 56690828 (F5), 56690829 (profiled F5), 56692120 (F3 + F5),
56692121 (profiled F3 + F5), 56693341 (profiled F3); G1/G3 56593059, 56593060, 56593251, 56691582, 56691584, 56691586, 56692321 and the deformation3D reruns 56690830 to 56690832;
G4 on i01 56690854 to 56690856, 56691583, 56691585, 56691587, 56692322 (the i02 series 56593248, 56593250, 56593340 discarded); profiles 56593412 (F2, i02), 56691588 (F4 + F5),
56692323 (F3 + F5), 56693351 (F3, final); 64 ranks 56693349, 56693350; G2 orchestrators of F1 56691605 to 56691608 cancelled (duplicate chain) and rerun as 56692489 to 56692492, F2 and F3 by
the chain jobs 56694336 and 56694440 (results in `perf/g2/`).
