# Manifest of the leia knowledge base notes (slugs are final; owners: A, B, C, D, E, F, ME)

Format: `path` | title | scope (what the note records) | sources to read | owner

## hubs/ (owner ME)
hubs/advection | Interface advection and the phase indicator | hub | | ME
hubs/viscosity | Viscosity | hub | | ME
hubs/surface-tension | Surface tension | hub | | ME
hubs/mass-flux | Mass flux and density-ratio consistency | hub | | ME
hubs/gradient-control | Gradient control: source terms and velocity extension | hub | | ME
hubs/verification | Verification methodology | hub | | ME
hubs/method-lines | The method lines and their status | hub | | ME

## models/ (one per runtime-selectable family; the members are rows of a table with: member, dictionary word, status, one-line verdict, evidence link)
models/sl-reconstruction | slReconstruction: the semi-Lagrangian value fit | members linearTaylor, linearWeightedLeastSquares, signedDistanceLinearWLS, quadraticTaylor, quadraticWeightedLeastSquares, uncachedQuadraticWLS (production), signedDistanceQuadraticWLS, bandQuadraticWLS, defectCorrectedIDW; keys stencil point/face, stencilBoundaryFaces, quadraticPivotTol 0.3, SL_FIT normalEquations vs QR; why a value fit, degree >= 2 load-bearing | M:98-125, M:373-377, S:1550-1722, S:1874-1990, S:2269-2278, SL article sec 304-548, kb-raw A1 | A
models/sl-scheme | slScheme: pointValue, fluxForm, normalProjected | pointValue (production; footIntegrator taylor/rk2; trajectoryVelocity input/normalProjection/normalClosestPoint), fluxForm (conserves int psi to 2e-14 but 3-13x volume error), normalProjected (closed) | M:376, MC article 258-304, nPSL article, kb-raw A3, A8 | A
models/sl-value-bound | slValueBound and the clips | none (production), stencilBounds (falsified as a fix), lipschitzCone (gate passes, ladders falsify), SL_CLIP history, limiters collapse the order | M:378-380, M:488-729, M:777-790, S:1822-1844, S:2008-2039, S:2312-2495, DP:57-158, kb-raw A6 | A
models/level-set-advection | levelSetAdvection: eulerian vs semiLagrangian | the two transport lines in the kinematic solver; where each is production; orders | MC article, SL article sec 1200-1330, kb-raw A1, A4 | A
models/phase-indicator | phaseIndicator: heaviside, sharpJump, geometric, detrixheAslam | detrixheAslam production (order ~2, geometric = DA to 8 digits), the SDF assumption in the offset (open), determinant guard, DA volume drift 6 % at N=32 | M:129-140, DP:5-9, SL article sec 621 and 1531, PCS:863-878, RM:1339-1352, RM:170-177, kb-raw A12 | A
models/narrow-band | narrowBand: none, empty, signChange, neighbours, distance, phaseIndicator | which is used where; the halo-dilation incident (band decomposition-dependent); band metrics use gradPsiMetric | C (seam section), S:3261-3269 | A
models/redistancer | redistancer: noRedistancing, PDE, anchoredEikonal, planeFootWave | line closed: planeFootWave second order static, PDE reinit divergent ("stability bomb"), frozen-band injurious (E_vol 0.017 -> 1.77), one-signed compounding | GRL article, MC article 304, PCT:599-604, PCS:378-381, M:777-790, DP:495-497, kb-raw A5 | A
models/volume-correction | volumeCorrection: noVolumeCorrection, newtonShift | global volume correction is a crossover, not a fix | SDPLS article sec 2507, kb-raw B2 | A
models/surface-tension-force | surfaceTensionForce: the twelve models | constantCurvature*, curvaturePressurePotential, correctionKang, divGrad*, traceGrad*, integralSurfaceTension, integralConormal, isoCurvature, reconstructedCurvature (production); verdict per member | SL article sec 919 and 1658-1852, SL negative-results deck, kb-raw C2, C5 | B
models/curvature-extension | curvatureExtension: the delivery words | none, closestPointNewton, fvm, harmonicLaplace, footPointHeightFunction, connectedInterface, interfaceMean, cellCentreInverse (production, case-dependent), stabilizedFootPointFace, cutCell, cellMean, symmetricFaceMean, footPointEvaluated, cellFootPointEvaluated; each with order, gain, verdict | PCS sections 6-18, M:4 and 8.1, S:2685-2688, kb-raw C2, C3 | B
models/semi-implicit-capillary-force | fv::option semiImplicitCapillaryForce | Hysing/Raessi/SAAMPLE forms; the collective-in-master-guard deadlock (2026-08-28); -45.3 vs -52.0 1/s: not needed with projectedFlux; untested after the fix on the translating droplet | C (4-rank section), S:1080-1096, kb-raw C4 | B
models/viscosity-face-model | viscosityFaceModel: alg_lin, alg_harm, alg_blend, geo_lin, geo_harm, geo_blend | alg_lin decided 2026-09-03 (36-arm mufGrid2D ladder; only positive orders; L1 1.400e-3, p 1.10 at ratio 1000); retracted switch to geo_lin (frozen muf bug 39e59b3); Popinet's face properties; Eulerian solver froze muf until c094bd8; 3D templates lack the token; printMethodBanner prints the wrong default | DP:701-763, M:397, S:139-143, S:164-169, SL article sec 759-875, kb-raw B | B
models/mass-flux | massFlux: interpolatedDensity, geometricFaceDensity, rhoLENT | the family and its sub-switches (alphaFSource, boundRho, alphaFTimeLevel, projectMassFlux); rhoLENT production; what is measured and what is VOID | M:6, M:8.1 rows, S:0, S:11.13-11.16, DP:817-869, kb-raw D | C
models/gradient-control-law | gradientControlLaw: the eight laws and the strain weights | | | ME
models/sl-source | slSource: none, gradientControl | | | ME
models/sdpls-source | sdplsSource: noSource, R, beta, Rdiv, RdivStrictSp, gradientControl | | | ME
models/velocity-extension | velocityExtension: the eight models | | | ME

## concepts/ advection (owner A)
concepts/departure-foot-ab2-centring | The departure foot: AB2, departure-centred | why the arrival form fails (+dt^2 d_t u; 2-4 % early, 35-47 % by t=0.02 s oscillating), taylor vs rk2 equal, kernel order 3.00 | M:69-96, M:376, PCS:1364-1372, S:3736, SL article sec 231 | A
concepts/polyhedral-fit-amplification | The amplification bound of the fit on polyhedra | Lambda_c = |1-sum g| + sum|g|, hex 1.0527 vs pMesh 1.2608; the pivot census; not a dt effect; polyhedral cells not cfMesh; Lambda is not a proxy for rho(B) | S:2121-2310, M:402-484, M:423-435, sl_fit_amplification.csv, SL article | A
concepts/idec-defect-correction-failure | iDEC defect correction diverges | 39 -> 1.4e4 -> 4.7e9; rho > 1 | SL supplementary iDEC_failure_report.tex, kb-raw A1 | A
concepts/trace-velocity-projected-flux | The trace velocity: projectedFlux and the reconstruct operator | fullHorizonStability2D: -52 vs +118 1/s; 70 % reconstruct, 0 % extension, 30 % solenoidality; identical to cellCentred under a frozen uniform stream; traceFlux physical/extension; SAAMPLE warning | S:1050-1164, M:375, M:437-455, S:3300-3303, PSH:791-805, DP:1029-1035 | A
concepts/psi-outer-correctors | psiOuterCorrectors: re-advect psi in every outer iteration | frozen-force lag exonerated; gain3D +0.4/-0.1 %; default yes since 5cbfaaa on consistency grounds; 3 frozen = 12 re-advected | PCS:1352-1363, PSH:207-236, DP:403-449, M:314-317 | A
concepts/static-local-refinement | Static local refinement around the interface | hexRefined/polyRefined; 51 640 vs 216 000 cells; Delta p 145.470 vs 145.48; orders agree to 0.09; 3.4-6.9x fewer core-hours; hanging nodes do not break the balance | SL article sec 985-1090 and 1331-1530, S:1197-1310, workflow/README | A
concepts/linear-semi-lagrangian | The linear semi-Lagrangian line | nestedLSQ ~1.1 shape, 1.5 volume, 2.1 band gradient; CFL 1 unstable at N >= 128; linearTaylor gradient defect O(1e9); two-phase workhorse; 3D 2.28/1.52 | LSL article, S:503-518, RM:22-28, kb-raw A2 | A
concepts/normal-projected-sl | The normal-projected semi-Lagrangian line (closed) | trace clean 0.017h; write-back diverges x1.7/10 steps; corrugation hypothesis falsified 0.209h vs 0.223h; orders 2.9 vs -1.24/-0.17; do not promote | nPSL article, PCT:78-95, PCT:167-187, RM:136-155, RM:302-327, RM:525-534 | A
concepts/eulerian-fv-transport | Eulerian FV transport of the level set | limiter drops the order 3.0 -> 0.9; flux-form loses 17 % volume; SL 20x more accurate at half the wall clock at 512^2; decision table | MC article sec 93-194, 258 | A
concepts/redistancing-geometric-grl | Geometric redistancing (GRL), the closed line | planeFootWave 4.48e-5 at h=1/256 order 2.0; anchoredEikonal first order; PDE reinit increases the error; foot-cloud scalloping; frozen-band injurious | GRL article and decks, MC article 304, PCS:302-344 | A
concepts/value-bounds-and-clips | Clips and value bounds: what was falsified | the clip on polyhedra, G4 falsifies, lipschitzCone retracted the same day, mesh-noise floor diagnostic (0.0/0.0/14.2 %), limiters collapse the order | S:1822-1844, S:2008-2039, S:2312-2495, M:488-729, C (regression set) | A

## concepts/ surface tension and viscosity (owner B)
concepts/balanced-force-csf-flux | The balanced-force capillary flux | forces in flux space; constant kappa absorbed to 3e-11; the projection absorbs 99.85-99.98 %; the flux-space residual; 30x tighter solve moves at most 2.4e-4 | M:243-272, SL article sec 2200, PCS:546-602, SL deck two-phase flow | B
concepts/curvature-from-the-fit | Curvature from the quadratic fit | symbolic kappa, parallel-curve offset correction 11.5 % -> 0.46 %, second order on a clean circle, O(h^1.2) on the moving interface 35 % -> 1 % | M:4, SL article sec 1603 | B
concepts/cell-centre-inverse-curvature | cellCentreInverse with the K-aware inverse | +2.01 non-gradient content on constant curvature; 4.60/3.81/1.55/1.49x lower residual; 93.5 % unfilled-cell inversion bug fixed; second order on the signed-distance ellipse and ellipsoid (2026-09-29; 2.10 with K, 1.00 without); the sphere h^1.95 vs h^1.02 is the per-face inverse; needs a parallel foliation | M:4.1 (CORRECTED), S:846-890, RM:1362-1396, kb-raw C2, C3 | B
concepts/face-curvature-deliveries | Face curvature deliveries and their gain | stabilizedFootPointFace h^2.04 ellipse 1.98; cutCell blows 3.3x sooner; cellMean; footPointEvaluated falsified; symmetricFaceMean 1.10; connectedInterface halves m=2 fails coupled; acceptance criterion G h^2 <= 0.65 and ellipse order >= 1.9; the ellipse gate collapses one-value-per-cell deliveries to first order | PCS sections 6-18, S:2685-2688, kb-raw C2 | B
concepts/parasitic-current-mechanism | The parasitic-current mechanism: source and amplifiers | max|U|(T) = u0(h) exp G(h); u0 = curvature error independent of U0 (step-1 kick to 0.6 %); exact kappa removes it 2e4-7e4x; G grows with U0 and density ratio; u0 ~ h^3.5, G ~ h^-3.27 (3D); t_blow ~ N^-3/2; loop model gamma ~ 0.647 sigma dt^2/(rho h^3); the m=2 mode; Popinet reframing | S:74-125, S:190-228, S:981-1005, PSH:244-324, PSH:398-583, PCS:1300-1372, PCS:1485-1592, SL article sec 1658-2058 | B
concepts/curvature-corrugation-and-the-fit | Curvature corrugation, aliasing and the psi filter as instrument | m > 4 modes; 0.209h corrugation at 16h displacement; psi filter 5.86x better / 1.61x worse; psi side supplies ~21 %; band renormalisation 3.4x worse; filtered results predate the seam bug | PCT:144-149, S:892-922, S:986-1014, PCS:1155-1298, C (no filtering) | B
concepts/integral-surface-tension-cst | Integral surface tension (CST) and the meta-law | 2D static equilibrium 5e-7 at N=64; diverges at N=128 t~0.047; "better static balance, higher dynamic gain"; conormal prototype 49.45 % residual; translating runaway 1.4e-3 s; sign trap | PCS:371-373, RM:157-256, RM:411-442, SL article sec 1658, SL negative deck | B
concepts/kang-gfm-and-sharp-heaviside | Kang GFM correction and the sharp Heaviside | 58x better statics, earlier blow-up (4.4e-3 rising vs 1.3e-5); the "58x/0.07 s" attribution flag | T2:215-217, kb-raw C2, H11 | B
concepts/force-time-centring | Time centring of the capillary force | endStep vs midpoint; the amplification matrix; midpoint spectrally identical, 32 % earlier divergence (2026-09-27); force at n+1 not n | S:3737, MC deck time centring, kb-raw C4 | B
concepts/capillary-time-step | The capillary time step | 0.2323 Brackbill; dt = 0.010861 N^-1.5; ~90 % of the growth dt-proportional at N=128; dt sweep within 6.5 %; polyhedral audit | S:2184-2231, S:2497-2523, C (no partial solutions) | B
concepts/pressure-projection-and-linear-solvers | The pressure projection and the linear solvers | operator pair, rAUf, tolerance, solver gates (roadmap 2026-07-28); cancellation-dominated convergence; non-orthogonal caveat (strict PCG 18-29x); a diverging smoother destroys diagnosability | RM gates, C (solver convergence section), PCS:546-602 | B
concepts/variational-capillary-force | The variational capillary force (proposal) | f_c = -sigma dA_h/dpsi_c; interFoam bounded by operator-transport-state pairing, not accuracy; kappaClamp/freezeKappa: the pump is the live curvature refill; no implementation | SL article sec 2275-2447, PSH:585-656, PSH:807-889, C:403-407 | B
concepts/well-balanced-exact-curvature-gate | The well-balanced gate with exact curvature | Delta p = 145.470 vs 145.48; spurious velocity <= 3.8e-10; identically zero velocity in all six arms; separates force balance from estimator | SL article sec 1379, S:190-228, C (research loop step 3) | B
concepts/viscosity-open-items | Viscosity: open items | 3D templates without the token; geometric sharp alpha_f consistent with the plane; mufGrid2D ran at np 8 before the parallel fixes; the banner default | S:3587, kb-raw B, H9 | B

## concepts/ mass flux (owner C)
concepts/rholent-mass-flux | rhoLENT: the mass-momentum-consistent flux | Liu 2023 eq. 40; auxiliary density equation, reset after the loop, residual 1e-13; stationary +1.0/-22/+0.1 %; translating VOID; formulation from METHOD 6 | M:6, SL article 2015, RM gate 2, S:0, kb-raw D2 | C
concepts/bound-rho | boundRho: active, evidence void | true since 28a1383; basis rhoDdtGate2D VOID; rhoClipFraction counts round-off (not scored); OPEN | S:3584, G2:66-69, DP:668, kb-raw D3 | C
concepts/alphaf-source-donor-plane | alpha_f source: donorPlane and the alternatives | donorPlane by author instruction (unmeasured); averagedPlanes VOID/diverged; donorPlaneAdvected ~1 %; trapezoid marginal; projectMassFlux 47x on the reducible part but stops the droplet -- all pre-closed-box-fix, VOID | DP:817-869, kb-raw D4 | C
concepts/mass-flux-projection | Projecting the mass flux (curl-free correction) | compatibility needs zero mean; what it did pre-fix (VOID) | DP:817-869, kb-raw D4 | C
concepts/ddt-scheme-pairing-bdf2 | Pairing the ddt schemes: BDF2 everywhere | the matching argument; the pairing tables VOID; flip-flop history f7307b5/28a1383/b60e3df; MOMENTUM_DIV upwind inert on the stationary droplet; BDF2 vs Euler +11.1/+2.9/-3.0 % (noise) | C (BDF2 section), kb-raw D5, M:8.1 | C
concepts/eulerian-solver-mass-flux-port | The Eulerian two-phase solver: frozen rho and the rhoLENT port | frozen rho until 2026-09-27; travelled fraction 0.420 -> 1.00012; np4 = serial 1e-8; shared headers createMassFluxFields/updateFaceDensity/updateMassFlux; the five Eulerian studies to void OPEN | S:11.16 (S:3703-3714), c094bd8, kb-raw D9 | C
concepts/coupled-face-density-defect | The coupled-face density defect (fixed 28d13f0) | rho_f 90 % apart across seams; seam residual 1e-5..5e-4 -> 1e-8..1e-10; the metrics counted processor faces (b1798c3); which earlier parallel studies to re-run OPEN; RM:510-511 flagged it in July | S:11.13-11.14 (S:3620-3680), kb-raw D8 | C
concepts/density-ratio-amplifier | The density ratio is the amplifier, not the mass-flux consistency | dec002f FALSIFIED mass-momentum consistency as the dominant term (9 orders in the residual, <2x in the excess); kickOriginGate2D: ratio 1 completes; Popinet at ratio 1 4.86 % of U | kb-raw D2, D6, S:74-125, SL article sec 2015 | C

## concepts/ verification (owner C)
concepts/method-gates | The method gates (2D and 3D) | arms exact1D, shear (+seam), stationary, translating (+translatingSeamNp1), oscillating; criteria completion / regression 10 % / order -0.3 / target / seam; INVALID; the two blind spots found 2026-09-28: exact-1D closed form None for candidates (make_gate_summary.py:190-198), centred band metric; the repairs | C:603-642, PHL:448-595, config/gates/methodGate2D.yaml, workflow/README, S:11 | C
concepts/error-vector-and-read-out-instants | The error vector and the read-out instants | shape, volume, |grad psi|, spurious current, pressure jump together; L2/L1 never L_inf; reversed flows: gradient at T/2, shape at T, volume at both; equal step counts; the oscillating arm by period and damping | C (research loop steps 5-6), S:196-198, PHL:554-560, memory rules | C
concepts/richardson-ladders-and-orders | Richardson ladders and observed orders | three rungs minimum; the fourth rung falsified a trend; matched horizon and dt law; report the order next to the error; the 2D orders of METHOD 8.3.7 were 3/2 too high | C (mesh convergence section), S:3237-3242, workflow/scripts/richardson.py | C
concepts/seam-checks-and-decomposition-invariance | Seam checks: code that is right on one rank and wrong on several | the class: setVelocity (2026-08-26), updateFlux, psi filter patches f83a1ab, band dilation, one-sided kappa_f, semiImplicit collectives, face density 28d13f0, metrics b1798c3; the 4-rank gate; serial-vs-np4 configs | C:199-246, S:552-640, docs/gradU-coupled-patch-contamination.md | C
concepts/wrong-setup-voids | A wrong setup voids its data | closed box 2026-09-02; tilted wall faces; stale arms; vacuous snappy arm; algebraic psi; 3D translating reference velocity; frozen rho; rename _VOID_ and re-run; check constant/polyMesh/boundary | C:661-705, S:14-65, S:1738-1845, kb-raw F | C
concepts/bit-identity-and-inertness-gates | Bit-identity and inertness gates | compare_metrics_csv.py --skip timing columns; DICT_MODE; "nothing is inert until measured" (residualControl tolerance 0); a finished 0/ is not the initial state | C:720-729, C:778-788, S:3010-3113 | C
concepts/cluster-provenance-and-binaries | Cluster provenance: binaries, ledger, ssh output | per-clone platforms/; the shared-account scancel incident; the .my_jobs ledger; the touch sweep 2026-09-22; ssh loses stdout; verify the artefact | C (cluster and ssh sections), CLUSTER.md, kb-raw F | C
concepts/log-classifier-and-waiters | Reading solver logs and waiting on them | foam_log_state.sh states and exit codes; the trapFpe false positive (three strikes); every waiter has a timeout; pgrep traps | C:263-296, workflow/scripts/foam_log_state.sh | C
concepts/data-archive-per-version | The data archive per software version | docs/<theme>/<slug>-article/data/archive/<stamp>/ with README and MANIFEST (SHA-256); make_archive.py; every curated number regenerated by a committed script | docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/, C (provenance) | C
concepts/advection-regression-set | The standing advection regression set | hex 2D (2Dvortex or 2Dtranslation), hex 3D (3Dshear), poly 3D (3Dshear); the mesh-noise floor 0.0/0.0/14.2 %; the slValueBound lesson | C:539-601 | C

## cases/ (owner C except where noted)
cases/stationary-droplet | The stationary droplet (2D and 3D) | geometry, La, the metrics, what it tests (source of the parasitic current), the ladders, the 6R box, exact-curvature arms | SL article sec 1658, S:718-773, cases/stationaryDroplet2D | B
cases/oscillating-droplet | The oscillating droplet | period and damping read-out; algebraic psi void (DROPLET_SURFACE token); the baseline drift at N=200 within ten periods; the gate arm | S:3271-3279, S:3946-3968, gcls article results, cases/oscillatingDroplet2D | B
cases/curvature-static-gates | The static curvature gates: circle, ellipse, ellipsoid | the ellipse gate collapses one-value-per-cell deliveries; the K-aware inverse on the sphere; curvatureDroplet2D | cases/ellipseDroplet2D, cases/ellipsoidDroplet3D, PCS, workflow/Snakefile.curvature | B
cases/translating-droplet | The translating droplet (2D and 3D) | the closed-box void; inlet/outlet; the free-stream gate; the late instability (outlet trigger, interior growth ~30 1/s shrinking with h; 20 mm vs 40 mm boxes; start-position test); END_TIME 0.05 in the gate | S:0, S:11.15 (S:3720-3790), gcls article sec 789 (moves to the SL article), cases/translatingDroplet2D | C
cases/popinet-translating-droplet | Popinet's translating droplet benchmark | reproduced 2026-09-04: Linf 4.86 % vs ~5 %, orders 0.49/0.91/1.71; the poly 3D mesh defect (tilted wall faces) VOID; poly ladder diverges late | S:130-189, S:1514-1845, SL article sec 2107 | C
cases/kinematic-advection-cases | The kinematic advection cases | 2Dvortex, 2Dtranslation (the only O(1)-displacement gate), 3Dshear, 3Ddeformation, 3Drotation, 3Dtranslation, contact-line cases; the boundary-stencil trap (u = 0 on walls) | C (step 5 boundary paragraph), cases/*, workflow/README | A
cases/exact-1d-stretch | The exact 1D stretch (1Dstretch) | u = alpha x, q = exp(-alpha t) closed form; the gate arm; the vacuous closed form for candidates | cases/1Dstretch, methodGate2D.yaml, make_gate_summary.py:185-200 | C
cases/benchmark-cases | Every case under cases/: one table | geometry, closed form or reference, which parts it tests, gate figure | cases/, workflow/README, kb-raw C (Cases) | C

## studies/ (owner F)
studies/sl-quadratic-pre-print | The quadratic semi-Lagrangian pre-print and decks | sections, claims with numbers, decks (main and negative results), data folders, what moved in on 2026-09-28 (parallel consistency, rhoLENT, late instability, oscillating drift) | SL article, kb-raw B1 | F
studies/sl-linear-pre-print | The linear semi-Lagrangian pre-print | | LSL article, kb-raw B5 | F
studies/grl-pre-print | The geometrically redistanced level set pre-print (draft) | | GRL article and decks, kb-raw B6 | F
studies/sdpls-pre-print | The SDPLS pre-print and decks | | SDPLS article, kb-raw B2 | F
studies/npsl-design | The normal-projected SL design note and pre-print | | nPSL article and notes, kb-raw B9 | F
studies/velocity-extension-pre-print | The velocity extension deck and the stub article | the 162-section deck: models, static verification, advected verification, verdicts by problem | VE deck, kb-raw B7 | F
studies/method-comparison | The method comparison article and deck | decision table; capillary balance track; time centring | MC article and deck, kb-raw B3 | F
studies/gcls-pre-print | The gradient-controlled level set pre-print (2026-09-27) and its archive | sections; the data archive per version; what moved out on 2026-09-28; the technical report | gcls article, its README, kb-raw B4 | F
studies/curvature-stabilization-campaign | The curvature stabilization campaign (plan v0.2, 2026-08) | sections 6-18 measured; the contamination notice; what it decided | docs/plan-curvature-stabilization.md | F
studies/shannon-parasitic-currents-campaign | The Shannon parasitic-currents campaign | sections 0-0g evidence; the Polya sections; execution phases; what it decided | docs/plan-shannon-parasitic-currents.md | F
studies/poly3d-roadmap | The capillary level-set research roadmap (gates 0-5, 2026-07/08) | dated records 07-27..08-07; gates; what stands and what was superseded | docs/capillary-level-set-research-roadmap.md, kb-raw H12 | F
studies/method-gate-2d-campaign-2026-09 | The first 2D method-gate campaign (2026-09-26/27) | | | ME

## decisions/ (owner D; one note per row group of METHOD 8.1 plus the process decisions; each with a line in decision-log.md under its month)
decisions/sl-reconstruction-uncached-qwls | SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares | | M:8.1, M:98-125, M:373 | D
decisions/sl-fit-normal-equations | SL_FIT normalEquations (QR bit-identical, blows up identically) | | M:377, S:2269-2278 | D
decisions/sl-trace-velocity-projected-flux | SL_TRACE_VELOCITY projectedFlux | | S:1050-1164, M:375 (wrong gate cited: say so), DP:1029-1035 | D
decisions/sl-clip-and-value-bound-off | SL_CLIP false and SL_VALUE_BOUND none | | M:378-380, S:2312-2495, M:488-729 | D
decisions/psi-filter-none | PSI_FILTER none (no filtering in production) | | C:382-398, S:892-922, S:986-1014 | D
decisions/mass-flux-rholent | MASS_FLUX rhoLENT | | M:6, M:8.1, S:0 | D
decisions/mass-flux-alphaf-donor-plane | ALPHAF_SOURCE donorPlane (author decision) | | DP:817-869, M:8.1 | D
decisions/mass-flux-bound-rho | BOUND_RHO true (open: the basis is void) | status open | S:3584, G2:66-69 | D
decisions/phase-indicator-detrixhe-aslam | PHASE_INDICATOR detrixheAslam | | M:129-140, DP:5-9 | D
decisions/surface-tension-reconstructed-curvature | SURFACE_TENSION reconstructedCurvature | | M:8.1, SL article sec 919 | D
decisions/curvature-extension-cell-centre-inverse | CURVATURE_EXTENSION cellCentreInverse (case-dependent: none for the Popinet family) | | M:4.1, M:8.1, DP, PHL:695-704, H10 | D
decisions/curvature-inverse-gaussian | The K-aware (Gaussian-curvature) inverse | | RM:1362-1396, S:846-890 | D
decisions/viscosity-face-model-alg-lin | VISCOSITY_FACE_MODEL alg_lin (3D open) | | DP:701-763, M:397 | D
decisions/momentum-schemes-bdf2-upwind | MOMENTUM_DDT_SCHEME backward, MOMENTUM_DIV upwind | | C (BDF2 section), kb-raw D5 | D
decisions/mesh-family-hexahedral | The mesh family: hexahedral production, polyhedral as the amplification rung | | S:2233-2310, M:402-484 | D
decisions/process-gates-2d-first-and-no-best-yaml | Process: 2D gate before 3D; no config/best.yaml; METHOD.md updated in the same commit | | C (best configuration section), M:25-29 | D

## retractions/ (owner E; each with a line in retraction-log.md under its month)
retractions/closed-box-translating-droplet | "The translating droplet ran": the case was a closed box (2026-09-02) | | S:14-65, C:661-705, 440107f | E
retractions/distance-cone-bound-as-transport-bound | The distance-cone (lipschitzCone) bound as a transport gain (2026-09-10) | | M:488-729, C (mesh convergence section) | E
retractions/clip-damage-is-the-narrow-band | "The clip's damage is the narrow band" (2026-09-09) | | S:2389-2495, S:2312-2387 | E
retractions/polyhedral-popinet-3d-mesh-defect | The polyhedral Popinet 3D results: tilted wall faces (2026-09-05) | | S:1738-1845 | E
retractions/t-blow-baseline | The t_blow baseline as a proxy (2026-08-31) | | S:1050-1068, PCS:1311-1313 | E
retractions/psi-filter-seam-bug | Filtered results predate the psi-filter seam bug (2026-08-19) | | S:892-922 | E
retractions/gradu-coupled-patch-contamination | 31 parallel kinematic studies contaminated by gradU on coupled patches (2026-08-26) | | S:619-640, S:580-617, docs/gradU-coupled-patch-contamination.md | E
retractions/cell-mean-delivery-adoption | The adopted cellMean delivery, retracted by the varying-curvature ellipse gate | | PCS sections 8-14 | E
retractions/force-at-n-not-n-plus-1 | The force built from psi^n (fixed: psi^{n+1}) | | C (BDF2 section), kb-raw C4 | E
retractions/advection-orders-3-2-factor | The 2D advection orders of METHOD 8.3.7 were 3/2 too high | | S:3237-3242 | E
retractions/late-translating-instability-is-the-outlet | "The late translating instability is the outlet" (corrected 2026-09-27: interior growth remains) | | S:11.15 (S:3720-3790), gcls article sec 789 | E
retractions/mass-momentum-consistency-dominant-term | Mass-momentum consistency as the dominant term (falsified dec002f) | | kb-raw D2, S:74-125 | E
retractions/gcls-coupled-loop-reading | The "coupled loop through the velocity at rate mu" reading of HL1q/HL1z (corrected 2026-09-28) | the linear laws read no velocity (gradientControlLaw.H:182-185); measured rates 245-276/s = mu 270/s; the crossing shift is O(h); see the technical report | src/leiaLevelSet/gradientControl/gradientControlLaw.H, gcls article sec 846, the memo summary in kb-manifest appendix | E
retractions/reversed-2dtranslation | "2Dtranslation is a uniform one-way translation": the case was a reversed flow (2026-09-29) | the dropped oscillation line and the cos(pi t/T) default; the silent error-reference fallback; the two fixed-psi setups found during the repair; the one-way re-run (orders 2.90/2.16/1.86; cone 14.1x and 9.5x at N = 64; null control 5.0e-12) | S:4167-4319 (11.19), S:3090-3091, S:3386-3387, M:405, M:414, M:608-741, C:528-532, C:598 (all aaa0a7dd); old text M:373-653 at 8867581 | session lead

## sessions/ (owner ME)
sessions/current | Current handover | | | ME
sessions/2026-09-27-gcls-first-campaign | Session 2026-09-26/28: gradient control, first campaign | | | ME
sessions/sl-session-handover | Handover to the semi-Lagrangian session | | | ME

## concepts/ gradient control (owner ME)
concepts/gradient-control-overview, concepts/sdpls-source-eulerian, concepts/sl-source-step, concepts/combined-source-note, concepts/halo-limited-extension, concepts/closest-point-extension, concepts/material-form-transport, concepts/why-the-candidates-failed, concepts/source-discretisation-defects, concepts/extension-strain-relocation, concepts/gradient-control-next-experiments, concepts/gradient-control-open-decisions

## Appendix: the red-team memo summary for retractions/gcls-coupled-loop-reading (owner E)
The gcls pre-print (2026-09-27, sec 846 "Discussion") explained the HL1q/HL1z stationary-droplet
destruction by the coupled loop psi -> kappa -> f -> u -> psi -> q -> F at rate mu = 1/T_ref, the
same loop that amplified the SDPLS current 260x. On 2026-09-28 a re-analysis showed: the linear
laws read no velocity (src/leiaLevelSet/gradientControl/gradientControlLaw.H:182-185: F = w(q^2)*a
+ G(q^2, sigma), strainWeight none gives w = 0, and linearQ/linearZ have no sigma); the measured
growth rates (archive gate/histories_stationary.csv, N = 200, window 10 to 20 ms) are about 290/s (curvature error) and 320 to 350/s (spurious current), of the order of mu = 270/s and far below the mesh-scale capillary rate of about 1e5/s (NOT the 245/276 of an earlier reading); HL1q and
HL1z agree to three digits until t = 0.02 s; the curvature error leads the velocity by ~5 ms. The
primary mechanism is now read as a linear instability of the explicit exponential source with a
CENTRED q (fvc::grad leastSquares): the mode delta_c = (-1)^c E d_c grows by (1 + mu dt) per step
(DERIVED; discriminator: a zero-flow static test). The crossing-shift formula is first order in h
and cannot produce the O(1) shape errors (orders 0.02-0.06). The reconstruct argument for HL0 is
neither necessary (FP0 has no shell, same order loss) nor sufficient. The sampler-across-the-jump
reading survives for the oscillating arm only. Scope: the pre-print's discussion and the
knowledge-base note why-the-candidates-failed are corrected; no data is void. Details: the
technical report docs/gradient-controlled-level-set/gcls-technical-report/ (written 2026-09-28).
