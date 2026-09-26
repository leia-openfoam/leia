# plan-halo-limited-gradient-control.md

Status: approved 2026-09-26. Branch `feature/gradient-controlled-level-set`. Progress is recorded in
STATUS.md section 11. The source dossiers stay out of git (`agent-input/`, see
`docs/gradient-controlled-level-set/README.md`).


## Context

**Why this work exists.** D. Bothe's technical dossier (revision 6, 2026-09-26,
`agent-input/ModifiedLevelSet_TechnicalDossier_HaloLimitedDirectionalExtension_20260926_115845_CEST.tex`)
proposes a one-scalar level-set method for unstructured, MPI-decomposed finite-volume meshes. The target
is a band 0.9 <= |grad psi| <= 1.1 near the interface. It is not exact signed distance and it uses no
redistancing. The method has three parts:

1. Material-form transport with the ACTUAL transport flux: `ddt(psi) + div(Phi,psi) - Sp(div(Phi),psi) = psi F`.
2. A halo-limited directional velocity extension: a sample point `Y = x - S_R(d) e` toward the interface,
   capped at a radius R, blended with the local velocity, `u_H = (1-w) u(x) + w u(Y)`. The blend has zero
   normal strain on the interface. The flux is `Phi^H = Phi^NS + w [U(Y) - U(x)] . S_f`.
3. A weak scalar source F: none, `lambda (1-q)`, or a smooth soft wall `-kappa tanh[gamma((q-1)/delta_s)^p]`.

The companion dossier (Bothe, Maric, Soga, `agent-input/Bothe-Maric-Soga-OpenFOAM-FVM-footpoint-source-dossier.tex`)
adds source laws in q and in z = q^2 with closed-form reference solutions and a test matrix.

**What the user asked for.**
1. Implement the method in leia and keep the runtime-selection (RTS) modularity.
2. Add a rule to CLAUDE.md and AGENTS.md: the level-set method must be extensible without modification of
   existing code (RTS and similar mechanisms). Eulerian methods run in the Eulerian solvers,
   semi-Lagrangian methods in the semi-Lagrangian solvers.
3. One Snakemake gate that runs four 2D cases on Lichtenberg at np 4: 2D shear (equidistant), stationary,
   translating and oscillating droplet. It reports the error vector of the existing studies (shape,
   gradient, volume, phase-indicator bounds) plus the MPI rank count and the wall-clock time. Reason: a
   method once gave good advection and bad hydrodynamics.
4. An equivalent coarse 3D gate on equidistant meshes with three resolutions.
5. Richardson best practice: the cell count doubles per rung (h ratio sqrt(2) in 2D, 2^(1/3) in 3D; the
   request said "cubic root (3)", but doubling the cell count in 3D gives 2^(1/3) = 1.26). The first rung
   has R/h of about 10, not less.
6. All studies run on Lichtenberg. A method change runs the 2D gate and a report. The 3D gate runs only
   after a 2D pass.
7. A docs theme folder for the pre-print.

**Intended outcome.** Every new source law, extension, travel law, weight law, direction or sampler is one
new class and one dictionary word. The 2D and the 3D gates run with one command each. Every gate result
carries a matched baseline on the same commit, the whole error vector, and the observed orders with
Richardson/GCI estimates.

**Decisions taken by the user (2026-09-26).**
1. The first campaign runs the semi-Lagrangian (SL) line. The Eulerian coupled arms join after Phase F
   repairs the Eulerian two-phase solver (frozen density, missing droplet metrics).
2. 3D ladder: h ratio of at least 1.3 per rung (Celik et al. 2008), about 2.2 times the cells per rung.
   2D keeps sqrt(2) = 1.414 (cells x2), which is above 1.3.
3. Docs theme: `gradient-controlled-level-set`, slug `gcls-level-set`.
4. The two dossiers stay out of git. The theme README cites them by title, authors and date.
5. ONE 2D study and ONE 3D study test every method; there are no per-method study configs. A method
   enters the gates only as a set of case tokens (a candidate), so a method that cannot be configured by
   tokens alone is not modular and cannot be tested. The gates reuse the existing cases and the existing
   `workflow/Snakefile` unchanged. The smoke run, the exact 1D check, the static gradient floor and the
   decomposition check are parts of the same two gate studies, not separate studies.

---

## 1. What exists today

### 1.1 Code: what the method can reuse, and what blocks it

| Need | Exists | Gap |
|---|---|---|
| Material-form FV transport | `applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H:20-29` (flux `phi`); `src/leiaLevelSet/advection/eulerianAdvection.C:172-204` (flux `velExt_->phi()`) | none for the flux form |
| Velocity-extension family | RTS `velocityExtension` (`none`, `closestPoint`, `steadyUpwind`, ...); `correct()` produces `Uext` and `phiExt` (`velocityExtension.C:118-124`) | the Eulerian two-phase solver has no call site and does not link the library; the SL two-phase solver uses the extension only with `traceVelocity cellCentred`, and the default `projectedFlux` traces `reconstruct(phi)` (`createTransportFields.H:142-156`); the kinematic SL solver has no extension |
| Flux contract | `interfaceExtension::updateFlux` (`interfaceExtension.C:400-470`) replaces the whole internal flux by `interpolate(Uext) & Sf` | the dossier needs the correction form, which keeps the pressure-corrected flux where w = 0 |
| Eulerian source family | RTS `sdplsSource` with the strategies `discretization`, `gradPsi`, `mollifier`; `nonLinearPart()` is the extension point | the selector lists the literal names `"Rdiv"`, `"RdivStrictSp"` (`sdplsSource.C:83-84`); `sdplsRdiv` finds its flux by the registry name `"phi"` (`sdplsRdiv.C:95-96`), which is not the flux that advects psi when an extension is active |
| SL source | none | no source hook in any SL solver; `docs/plan-combined-source-terms.md` WP3 designs `slSource` composed in `slAdvection`, not implemented |
| Foot points | `closestPoint` (Newton descent, known-vicinity walk, zoneDistribute halo) keeps its foot points private; `slReconstruction::footPointDistance` (byte-frozen) returns a signed distance only | no public foot-point API; `closestPoint` falls back to a steady solve when a query fails, and its halo Taylor data makes it decomposition dependent at seams |
| Partition-independent sampling | `extendedCentredCellToCellStencil`, used by `slReconstruction` (`slReconstruction.C:92-99`), seam-consistent to about 1e-12 | no velocity (vector) stencil fit yet |
| Eulerian two-phase physics | `leiaLevelSetTwoPhaseFoam` | `rho` and `rhoPhi` are frozen at t = 0 except on mesh refinement (`leiaLevelSetTwoPhaseFoam.C:149`; `IMPROVEMENTS.md` item 4, open); the solver writes no droplet metrics (no max|U|, no pressure jump) |

### 1.2 Studies and workflow

| Need | Exists | Gap |
|---|---|---|
| One command for four cases | `workflow/Snakefile` runs one case and one solver per config; meta-Snakefiles (`workflow/Snakefile.comparison` and others) run several configs | the meta-Snakefiles hard-code `profiles/local` and `--resources tasks=12` and skip `guard_finished_cases.py`; a base config is forbidden (Snakemake merges config files shallowly) |
| 2D arms | `2Dvortex` (the shear2D case; no `2Dshear` exists), `stationaryDroplet2D`, `translatingDroplet2D`, `oscillatingDroplet2D`; the closest configs are `traceKinematic2Dvortex`, `stationaryLadder2Dshared`, `translatingLadder2D`, `oscillatingLadder2Dshared` (all np 4) | most ladders start below R/h = 10 (N = 32, 64); `oscillatingDroplet2D` initialises an ALGEBRAIC level set (`implicitEllipsoid`: value = sum (x_i - c_i)^2/a_i^2 - 1, so |grad psi| is about 1.8e3 at the interface; `fvSolution.template:382`); `signedDistanceEllipse` exists and `oscISTDroplet2D` uses it |
| 3D arms | `3Dshear`, `stationaryDroplet3D` and `translatingDroplet3D` on the 6R box; 2^(1/3) ladder 60/76/95 | no 3D oscillating droplet case |
| Error vector | advection CSV in every solver (`E_GEOM_ALPHA_REL`, `E_VOL_ALPHA_REL`, `E_BOUND_ALPHA`, `ELAPSED_CLOCK_TIME`); `gradPsiError.csv` (kinematic); droplet CSV of the SL two-phase solver (`maxMagU`, `l2MagUPrime`, `phaseVolumeRelError`, `zeroSetRadialL2`, `pLaplace`, `kErrL2Band`, `gradPsiL2ErrorBand`, `m2Amplitude`, ...) | no exact pressure-jump column; the real rank count is only in the log; no Richardson extrapolation or GCI anywhere (`make_convergence_table.py` fits a least-squares order) |
| Cluster launch | `run-studies.sbatch` (orchestrator job), `make studies-one`, the launch guard, `.my_jobs`, `make pull-study` | none |

### 1.3 Defects to fix before a gate table can be trusted

| ID | Where | Defect | Effect |
|---|---|---|---|
| D-a | `workflow/scripts/aggregate.py:383-391` | the token copy sits inside `if isfile(leia.version)` | a case without `leia.version` loses N_CELLS, END_TIME and h |
| D-b | `workflow/scripts/advection_convergence_table.py:107` | h_eff = nCells^(-1/3) also in 2D | every 2D order in `advConv2D*_convergence.csv` is 1.5 times too high; METHOD.md section 8.3.7 (lines 584-618) quotes them (5.70 is really 3.80, 4.15 is 2.76, 4.85 is 3.23) |
| D-c | `workflow/Snakefile:79` | dims from the case-name prefix | every `*Droplet2D` case records dims 3 |
| D-e | `cases/oscillatingDroplet2D/system/fvSolution.template:382` | algebraic level set | the gradient columns and every band criterion in psi units have no meaning |
| D-f | `applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H:195` | the band |grad psi| uses the limited `fvc::grad(psi)` | the template requires the unlimited `gradPsiMetric` (`fvSchemes.template:104-108`); every two-phase template defines it |
| D-i | `aggregate.py:68` | completeness = t_last >= 0.99 END_TIME | short fixed-dt runs lose their final metrics (measured in `seamConsistency3Dpar4`) |
| D-j | `materialize.py:331`, `aggregate.py:329` | `np` = the configured value | a serial run records np > 1 |
| D-m | `leiaLevelSetTwoPhaseFoam.C:149` | rho and rhoPhi frozen at t = 0 | Eulerian coupled translating and oscillating results are wrong physics |

---

## 2. Rules added to CLAUDE.md and AGENTS.md (same commit, byte-identical)

New section after "Repo layout & git discipline":

> ## Extension without modification
>
> The level-set method is open to extension and closed to modification. New behaviour is a new
> runtime-selectable (RTS) class. It is not an edit of an existing model class, of a solver equation or
> of a selector.
>
> 1. A new model goes into the RTS family that owns the behaviour: one class with
>    `addToRunTimeSelectionTable` and one line in that family's `Make/files`. A dictionary word selects
>    it. A case token with an inert default renders that word.
> 2. If no family owns the behaviour, add a new family. Its base class is the inert default model. It
>    gets its own `Make/`, a row in `etc/leia-check-deps.py`, an entry in `Allwmake`, and an `-l` entry in
>    the `EXE_LIBS` of every solver that can select it.
> 3. An alternative inside a model is a strategy, and a strategy is also a RTS family (the pattern is the
>    `discretization`, `gradPsi` and `mollifier` of `sdplsSource`). An `if` chain over dictionary words is
>    not an extension point.
> 4. A solver holds composition roots only. It constructs each family with `New()` and calls the family
>    interface. A solver or a selector never tests the type name of a model.
> 5. A solver without a composition root for a family gets one as a one-time infrastructure change. The
>    change ships with the inert default, in its own commit, gated by a bit-identity run against the
>    pre-change state.
> 6. A model that needs data the family interface does not give extends the base interface once, with an
>    inert default, gated as in item 5. A model never looks up a field by a name that another class
>    chooses.
> 7. The Eulerian solvers and the semi-Lagrangian solvers expose the same families where the method
>    allows it. An Eulerian model runs in an Eulerian solver. A semi-Lagrangian model runs in a
>    semi-Lagrangian solver.
>
> MEASURED 2026-09-26: the rule was broken in four places. The Eulerian two-phase solver cannot select a
> velocity extension. The SL two-phase solver ignores any extension with the default `projectedFlux`
> trace. `sdplsSource::New` lists `"Rdiv"` and `"RdivStrictSp"` by name. `sdplsRdiv` finds its flux by
> the registry name `"phi"`.

New section after "The standing advection regression set":

> ## Method gates: 2D first, 3D on a pass
>
> Any change to the level-set method (a new model, a changed default, a composition root that is not
> inert) runs the 2D method gate before any other coupled study:
> `make gate GATE=methodGate2D CANDIDATES=<name> PROFILE=profiles/slurm`. The 3D method gate runs only
> after a 2D pass. The gates run on Lichtenberg. The laptop runs only the unit tests and the 4-rank smoke
> that "Run it on 4 ranks before it leaves the laptop" requires.
>
> 1. Every gate run carries the `baseline` candidate on the same commit and binaries. A verdict against a
>    baseline from another commit is not a verdict.
> 2. The resolution ladder is a Richardson ladder. In 2D the cell count doubles per rung (h ratio
>    sqrt(2) = 1.414). In 3D the h ratio is at least 1.3 per rung (Celik et al., J. Fluids Eng. 130,
>    078001, 2008), about 2.2 times the cells, because doubling (2^(1/3) = 1.26) puts the rungs too close
>    for a stable order estimate. Three rungs at least. The first rung has R/h >= 10. All rungs have the
>    same parity of N.
> 3. The gate reports the whole vector per arm (shape, gradient band error, volume, phase-indicator bounds,
>    spurious currents, pressure-jump error, curvature, ranks, wall clock), the observed order of every L2
>    and L1 metric, and the Richardson/GCI estimate for a quantity without an exact value. It never reports
>    an L_inf order.

---

## 3. Target architecture

### 3.1 Families

| Family | New or existing | Library | Dictionary | Models in the first release |
|---|---|---|---|---|
| `gradientControlLaw` | NEW | `libleiaGradientControl` (NEW, `src/leiaLevelSet/gradientControl/`) | `law { type ...; }` inside a source dictionary | `none`, `linearQ`, `linearZ`, `cubicQ`, `cubicZ`, `twoThirdsZReg`, `saturatedLinearZ`, `boundedGradient` (S0), `softWall` (S1) |
| `strainWeight` (strategy of the law) | NEW | `libleiaGradientControl` | `law { strainWeight { type ...; } }` | `none`, `full` (the SDPLS strain a), `omega` (S1c) |
| `sdplsSource` | existing | `libleiaSdplsSource` | `levelSet.sdplsSource` | adds `gradientControl` (the Eulerian consumer) |
| `slSource` | NEW | `libleiaSemiLagrangian` (`semiLagrangian/source/`) | `levelSet.semiLagrangian.source` | `none`, `gradientControl` (the SL consumer) |
| `velocityExtension` | existing | `libleiaVelocityExtension` | `levelSet.velocityExtension` | adds `haloLimited` |
| `extensionTravel` (strategy of `haloLimited`) | NEW | `libleiaVelocityExtension` (`velocityExtension/haloLimited/`) | `travel { type ...; }` | `capped` (S_R, parameter m); later `full` |
| `extensionWeight` | NEW | same | `weight { type ...; }` | `fractionReached` ((S/d)^beta); later `compact`, `designedK` |
| `extensionDirection` | NEW | same | `direction { type ...; }` | `levelSet` (d = psi/q, e = grad psi/q); later `fitFootPoint` |
| `extensionSampler` | NEW | same | `sampler { type ...; }` | `stencilFit` (weighted least-squares fit of U over the CPC stencil); later `cellPoint` |

Interface sketches:

```cpp
// src/leiaLevelSet/gradientControl/gradientControlLaw.H
class gradientControlLaw
{
public:
    TypeName("none");
    declareRunTimeSelectionTable(autoPtr, gradientControlLaw, Dictionary,
        (const dictionary& lawDict), (lawDict));
    static autoPtr<gradientControlLaw> New(const dictionary& lawDict);
    explicit gradientControlLaw(const dictionary& lawDict);   // builds strainWeight_
    virtual ~gradientControlLaw() = default;

    //- True if rate() reads sigma = |symm(grad U)|_F. The consumer computes sigma only then.
    virtual bool needsRate() const { return false; }
    //- G [1/s] in D_t psi = psi*F. q2 = |grad psi|^2. Zero at q2 = 1 for every law.
    virtual scalar rate(const scalar q2, const scalar sigma) const { return 0; }
    //- dG/dq, for the stiffness number Co_B = dt |q dG/dq|.
    virtual scalar dRateDq(const scalar q2, const scalar sigma) const { return 0; }
    //- The complete coefficient F = strainWeight(q2)*a + rate(q2, sigma).
    scalar F(const scalar q2, const scalar sigma, const scalar a) const;

protected:
    autoPtr<strainWeight> strainWeight_;
};
```

```cpp
// src/leiaLevelSet/semiLagrangian/source/slSource.H
class slSource
{
public:
    TypeName("none");
    declareRunTimeSelectionTable(autoPtr, slSource, Mesh, (const fvMesh& mesh), (mesh));
    static autoPtr<slSource> New(const fvMesh& mesh);   // reads levelSet.semiLagrangian.source
    //- psi holds the SL arrival value on entry and psi^{n+1} on exit.
    virtual void apply(volScalarField& psi, const volVectorField& Utrajectory,
                       slReconstruction& geometryFit, const scalar dt) {}
};
```

`haloLimited` derives from `velocityExtension` directly. It needs no alpha seed band and no steady solve,
and it must not inherit the full-replacement flux of `interfaceExtension`.

### 3.2 Composition roots (one-time edits; after them a new model needs no solver edit)

| Solver | Extension root | Source root |
|---|---|---|
| `leiaLevelSetFoam` (kinematic, `ADVECTION eulerian` or `semiLagrangian`) | exists in `eulerianAdvection` | `eulerianAdvection` passes `velExt_->phi()` through `setTransportFlux`; the SL path gets `slSource` through `slAdvection` |
| `leiaSemiLagrangeLevelSetFoam` (kinematic SL) | NEW: `velocityExtension::New` and `traceFlux physical|extension` | through `slAdvection` |
| `leiaLevelSetTwoPhaseFoam` (Eulerian two-phase) | NEW header `velocityExtensionFieldsEuler.H`; `alphaEqn.H` uses `velExt->phi()` in `div` and `Sp` | `source->setTransportFlux(velExt->phi())` |
| `leiaSemiLagrangianLevelSetTwoPhaseFoam` (SL two-phase) | NEW `traceFlux physical|extension` inside `slTraceField` | through `slAdvection` |

### 3.3 Link graph

- `etc/leia-check-deps.py`: `LIB["gradientControl"] = "leiaGradientControl"`; `LINKS`:
  `leiaGradientControl: {}`, `leiaSdplsSource: {leiaGradientControl}`,
  `leiaSemiLagrangian: {leiaGradientControl}`; the other rows do not change.
- `Allwmake` (and `Allwclean`): `leiaGradientControl:src/leiaLevelSet/gradientControl` after `leiaCore`.
- `EXE_LIBS`: all five solvers and the new test apps add `-lleiaGradientControl`; the kinematic SL solver
  and the Eulerian two-phase solver add `-lleiaVelocityExtension`.
- `README.md`: library table and link matrix get the new rows.

---

## 4. Work packages, in order

Every step names its files, its inert default and the gate that closes it. Every commit is inert or gated.

### Phase A: rules, plan document, docs theme (no code)

- **A1.** Add the two CLAUDE.md sections of section 2 to CLAUDE.md and AGENTS.md. Gate: `diff CLAUDE.md AGENTS.md` prints nothing.
- **A2.** Commit this plan as `docs/plan-halo-limited-gradient-control.md`. Put a "SUBSUMED by
  docs/plan-halo-limited-gradient-control.md" line at the top of `docs/plan-combined-source-terms.md`
  (WP2, WP3) and `docs/combined-source-terms/improvement-sdpls-combined.md`: `sdplsCombined` is
  `gradientControl` with `linearQ` and `strainWeight full`; WP3's `slSource` is adopted here.
- **A3.** Docs theme skeleton (section 7). Gate: `latexmk` builds the empty article; `paths.theme_of("<theme>")` resolves.

### Phase B: workflow repairs and gate infrastructure (no solver change)

- **B1.** Fix D-a, D-b, D-c, D-i, D-j.
  - D-a: move the token copy out of the `leia.version` branch.
  - D-b: h = nCells^(-1/dims), with dims from `case_params.json`.
  - D-c: a per-case `DIMS` entry in `cases/<case>.parameter` wins over the name prefix. Add `DIMS (2);` to
    every 2D case whose name does not start with `2D` (the 2D droplet family, `popinetTranslating2D`, the
    ISO/IST 2D cases).
  - D-i: completeness comes from `workflow/scripts/foam_log_state.sh` (COMPLETED), not from 0.99 END_TIME.
  - D-j: the `solve` rule writes a sidecar `<case>/.leia_launch` (np, SLURM_JOB_ID, host, start and end
    epoch seconds, exit code). `foam_log_state.sh` gets an additive `nprocs=` output field from the log
    header (no pattern change). `aggregate.py` reports `nRanks` (log) and `wallClockSolve` (sidecar).
  - Gate: re-aggregate three preserved studies (one without `leia.version`, `seamConsistency3Dpar4`, one
    2D droplet study) and diff the `_errors.csv`. Only the named columns change. Solver CSVs are not touched.
- **B2.** Retraction, on the same day as the D-b fix. Regenerate
  `docs/method-comparison/method-comparison-article/data/tables/advConv2D*_convergence.csv`. Correct
  METHOD.md section 8.3.7 (lines 584-618) with an explicit sentence: "the 2D orders were computed with
  h = nCells^(-1/3); the correct orders are 2/3 of the published values". Add a STATUS.md entry. Search the
  decks and articles for the old numbers and correct them.
- **B3.** Fix D-f: `writeDropletMetrics.H:195` uses `fvc::grad(psi, "gradPsiMetric")`. Gate: the SL
  droplet CSV of `stationaryLadder2Dshared` at N = 64, np 4, is identical at tolerance 0 except the four
  `*GradPsi*Band` columns (`compare_metrics_csv.py --skip` those four and the two clock columns).
- **B4.** Fix D-e: a token `DROPLET_SURFACE` with the inert default `implicitEllipsoid` in
  `oscillatingDroplet2D`. The gates pin `signedDistanceEllipse` (block pattern from
  `cases/oscISTDroplet2D/system/fvSolution.template:302-310`). Gate: render `oscillatingLadder2Dshared`
  before and after; the dictionaries are identical.
- **B5.** Gate infrastructure (section 5): `workflow/Snakefile.gate`,
  `workflow/scripts/render_gate_configs.py`, `workflow/scripts/make_gate_summary.py`,
  `workflow/scripts/richardson.py`, `config/gates/{methodGate2D,methodGate3D}.yaml`,
  `config/candidates/baseline.yaml`, the Makefile targets `gate` and `studies-one-file`, and a
  `workflow/README.md` section. Gates:
  - `python3 workflow/scripts/richardson.py --self-test` recovers p = 1, 2, 3 from synthetic
    f = f0 + C h^p with non-constant r to 1e-10, and flags oscillatory convergence. It exits nonzero on failure.
  - The 4-rank laptop smoke `make gate GATE=methodGate2D SMOKE=1 CANDIDATES=baseline PROFILE=profiles/local`
    completes all four arms, and the summary has every column filled (check `steps=`, not the exit code).
- **B6.** New case `cases/oscillatingDroplet3D`: a clone of `stationaryDroplet3D` with
  `signedDistanceEllipsoid`, axes a = 1.1 R and b = c = R/sqrt(1.1) (volume of the sphere of radius R).
  Check whether `m2Amplitude` applies in 3D; if not, add axis-extent columns to the droplet CSV (additive).
  Gate: the t = 0 volume and curvature check (the `initVolumeEllipsoid3D` pattern) and the 4-rank smoke.
- **B7.** Run the 2D gate with `baseline` only on Lichtenberg (section 5.6) and report it. This is the
  matched reference of the current production method and the first test for "good advection, bad
  hydrodynamics". Then run the 3D gate with `baseline` once.

### Phase C: open the extension points (one-time, bit-identical)

Before C1, preserve the pre-change gate set of `docs/plan-library-split-and-build-policy.md` WP3 (unit,
`sdpls1Dstretch`, the advection regression set, `stationaryDroplet3DbitIdentity`, the N = 32 arm of
`sdplsPsiBudgetDroplet2D`, the seam pair, one np 4 case per solver) as `studies/<study>.pre-extpoints-<date>`.

- **C1.** `sdplsSource`: replace the literal-name exemption by a virtual trait `usesDiscretization()`
  (true in the base, false in `sdplsRdiv` and `sdplsRdivStrictSp`). The selector constructs the model
  first and then checks. Add `setTransportFlux(const surfaceScalarField&)` and a protected
  `transportFlux()` accessor that is fatal if no flux was set. Gate: `leiaTestSdplsSource` (89
  assertions), `cases/sdplsSourceUnit`, `sdpls1Dstretch` 12 of 12 at tolerance 0, the Rdiv arm of
  `sdplsPsiBudgetDroplet2D`.
- **C2.** `slSource` base (`none`) composed in `slAdvection::advect` after `scheme_->advance`, with the
  scheme's effective dt (the Stormer-Verlet halves). `postAdvect` is a no-op for the production
  reconstruction, so the order is safe. Gate: the advection regression set (default arms, three rungs,
  tolerance 0), `stationaryDroplet3DbitIdentity`, the seam pair.
- **C3.** Eulerian two-phase extension root: header `velocityExtensionFieldsEuler.H`, included only from
  `leiaLevelSetTwoPhaseFoam.C` (the shared `createFields.H` cannot hold it: the SL solver already builds
  one `velocityExtension`, and a mesh holds only one). `alphaEqn.H` calls `velExt->correct()` on every
  outer corrector and uses `velExt->phi()` in `div` and `Sp`. Link `-lleiaVelocityExtension`. With `none`,
  `phiExt` is a copy of `phi` named "phi", so the `div(phi,psi)` scheme applies. Gate: the Eulerian N = 32
  arm of `sdplsPsiBudgetDroplet2D` and one np 4 Eulerian droplet case at tolerance 0.
- **C4.** SL two-phase `traceFlux physical|extension` (token `SL_TRACE_FLUX`, default `physical`):
  `slTraceField` reconstructs `slVelExt->phi()` for `extension`. A new word, not a changed meaning of
  `projectedFlux`: existing configs that set an extension with `projectedFlux` keep their meaning. Gate:
  `stationaryDroplet3DbitIdentity` and `stationaryLadder2Dshared` N = 64 np 4 at tolerance 0.
- **C5.** Kinematic SL solver extension root (`velocityExtension::New`, `traceFlux`, link). For the
  reversed flow the root extends the base (unscaled) velocity and flux once per step and scales both trace
  levels by the oscillation factor. This is exact, because the flux correction is linear in U at fixed
  geometry. Gate: `advConv2Dvortex`, `advConv2Dtranslation`, `advConv3DshearHex`, `advConv3DshearPoly`
  default arms at tolerance 0, and a 4-rank `2Dvortex` arm.
- **C6.** Tokens and template blocks with inert defaults (table below) in `cases/default.parameter` and
  in the templates of the gate cases and of the kinematic gate cases (`2Dvortex`, `2Dtranslation`,
  `3Dshear`, `1Dstretch`, `stationaryDroplet2D/3D`, `translatingDroplet2D/3D`, `oscillatingDroplet2D/3D`).
  Also correct the dead value `meshWaveExt` in `default.parameter:13` and the two-phase template comment
  (the type name is `meshWave`). Gate: render every committed config that uses a touched template before
  and after; the diff shows only the new inert entries. Then repeat the C1-C5 bit-identity runs on the
  final state.

| Token | Inert default | Lands in |
|---|---|---|
| `VELOCITY_EXTENSION` | `none` (exists) | `levelSet.velocityExtension.type` |
| `HL_RADIUS_CELLS` | `1` | `haloLimited.radiusCells` (the model refuses a value above 1) |
| `HL_M`, `HL_BETA` | `2`, `1` | `travel.m`, `weight.beta` |
| `HL_DIRECTION`, `HL_SAMPLER` | `levelSet`, `stencilFit` | `direction.type`, `sampler.type` |
| `SL_TRACE_FLUX` | `physical` | `levelSet.semiLagrangian.traceFlux` |
| `SL_SOURCE`, `SL_SOURCE_BAND_CELLS` | `none`, `3` | `levelSet.semiLagrangian.source.{type,bandCells}` |
| `SDPLS_SOURCE` | `noSource` (exists) | `levelSet.sdplsSource.type` |
| `GC_LAW`, `GC_STRAIN_WEIGHT` | `none`, `none` | `law.type`, `law.strainWeight.type` (in both source blocks) |
| `GC_MU` | `1` | `law.mu` [1/s]; read only by a law that needs it |
| `GC_EPS`, `GC_SAT_C` | `0.02`, `1` | `law.eps` (placeholder until E0 measures the floor), `law.c` |
| `SW_C_KAPPA`, `SW_DELTA_S`, `SW_P`, `SW_GAMMA`, `SW_EPS_D` | `1.25`, `0.08`, `5`, `1.4722`, `0` | `softWall` entries (the dossier reference values) |
| `OMEGA_BETA`, `OMEGA_M` | `1`, `2` | `strainWeight omega` |

### Phase D: the new models (each with unit tests and a 4-rank run)

- **D1. `libleiaGradientControl`**: `gradientControlLaw`, `strainWeight` and the nine laws. Every law
  parameter is REQUIRED in the dictionary (no code default); the tokens carry the dossier values.
  `softWall` refuses an even p. New app `applications/test/leiaTestGradientControlLaw`, exit nonzero on failure:
  - F(q = 1) = 0 exactly; the sign of F is opposite to the sign of q - 1.
  - The q and z pairs have equal slopes at q = 1 for equal mu.
  - Soft wall: Psi(0.08) = 0.9; the dossier table 0.011/0.139/0.336/0.638/0.900/0.990/0.9997 at 3-10 %;
    Theta* = 0.9 for (delta_A, eps_q, delta_s, p, gamma) = (0.10, 0.02, 0.08, 5, artanh 0.9).
  - RK4 integration of q' = q F reproduces the logistic solution (`linearQ`), the z-logistic solution
    `z(t) = z0/(z0 + (1 - z0) exp(-mu t))` (`linearZ`) and the cubic decay.
- **D2. `sdplsGradientControl`** (`TypeName("gradientControl")`) in `libleiaSdplsSource`.
  `nonLinearPart(R, psi, U)` returns `w(q2) R + G(q2, sigma)`. q2 comes from the `gradPsi` strategy
  (unlimited `gradPsiSdpls`); sigma = |symm(fvc::grad(U, "gradUSdpls"))| only if `needsRate()`. The
  `discretization`, `gradPsi` and `mollifier` strategies work unchanged. The model warns if the velocity
  extension is not `none` and the strain weight is not `none` (double cancellation of the strain).
  - Coverage of older ideas: SDPLS `R` = law `none` + `full`; `sdplsCombined` = `linearQ` + `full`;
    S0 = `boundedGradient`.
  - Tests: extend `leiaTestSdplsSource` (equality with `sdplsR` and `sdplsBeta` to round-off on a static
    field). The 1D closed forms run in the `exact1D` arm of the 2D gate (no new study config): with
    u = alpha x and psi = c(t) x, c' = c (G(c) - alpha) for strain weight `none` and c' = c G(c) for
    `full`; the summary integrates the same ODE independently (RK4) and reports the error.
- **D3. `slGradientControlSource`** (`TypeName("gradientControl")`) in `libleiaSemiLagrangian`.
  - Flat band gate |psi|/q <= bandCells h; outside the band F = 0 exactly.
  - q2 and n from the geometry fit of the arrival field; sigma and a = n . symm(grad U) . n from
    `fvc::grad(Utrajectory, "gradUSource")`.
  - Update psi_c <- psi_c exp(dt F_c); |dt F| is clamped at 30; every clamp is counted, logged with
    `gSum` (outside any `Pstream::master()` guard) and written to a field.
  - New app `leiaTestSlSource`: one step reproduces psi exp(dt F) exactly for planar psi = g0 x with
    U = 0; cells outside the band stay bit-identical; no value changes sign; np 4 equals serial.
  - The SL line's 1D closed form runs in the same `exact1D` arm (`leiaLevelSetFoam` with
    `ADVECTION semiLagrangian`, the same `slAdvection` path as the SL solvers).
- **D4. `haloLimited`** velocity extension. Per `correct()`:
  1. g = fvc::grad(psi, "gradPsiExtension") (unlimited leastSquares); psi_f and g_f by linear interpolation.
  2. Direction `levelSet`: d = psi/sqrt(Q2), e = g/sqrt(Q2), with Q2 = q2 + zeta(q2) and a flat
     regulariser zeta that is zero for q2 >= 0.25.
  3. R = radiusCells h_f. h_f is the coupled-aware cell size (the `interfaceExtension::computeCellSize`
     definition), independent of the decomposition.
  4. S = travel.S(d, R), w = weight.w(d, R), Y_f = x_f - S e_f. By construction |Y_f - x_f| < R.
  5. Sampler `stencilFit`: a weighted least-squares fit of each U component over the owner's and the
     neighbour's CPC stencil (`extendedCentredCellToCellStencil`). U_h(Y_f) and U_h(x_f) are the mean of
     the two fits. On processor faces the two ranks swap and average with `syncTools`.
  6. Flux, correction form: phiExt_f = phi_f + w_f [U_h(Y_f) - U_h(x_f)] . S_f. Physical patches keep phi.
  7. Uext_c = (1 - w_c) U_c + w_c U_h(Y_c), for the SL `cellCentred` trace and for diagnostics.
  8. Fields (AUTO_WRITE): `hlWeight`, `hlDistance`, `hlReach` (|Y - x|/R, always < 1), `hlCorrection`.
     Per step, collective and outside any master guard: the maximum of `hlReach` and the band L2 of
     div(phiExt) - div(phi).
  - Tests in `leiaTestVelocityExtension`: S(0) = 0, S'(0) = 1, |S| < R up to |d| = 1e6 R, S odd; w(0) = 1
    and 1 - w = O((d/R)^(2m)); planar signed distance with affine U: n . grad_h(u_H) . n on the interface
    faces is zero to 1e-12 relative (the linear fit is exact for affine U); uniform U: phiExt equals phi bit
    for bit; the dossier's interfacial test psi = (1.08 + 0.1 y) x + 0.2 x^2 (DY . n = 0 on x = 0); the
    dossier's nonlinear velocity u = (alpha x + beta x^2, -(alpha + 2 beta x) y): velocity blend minus
    point blend = beta w (1 - w) x^2; serial equals np 4 to 1e-12.
  - New app `leiaTestMaterialTransport` (dossier Test 0a): psi = const with a manufactured flux of nonzero
    divergence, and with phiExt; the material operator leaves psi unchanged to 1e-14.
- **D5. `closestPoint`** stays unchanged as the full foot-point reference (FP0, FP1). It replaces the whole
  flux and falls back to a steady solve, so it is a reference. It is not a promotion candidate unless it
  passes the seam gate.

After each D step: a 4-rank laptop run of the unit apps and of one gate arm at N = 32; check the steps and
the CSV; then push.

### Phase E: the campaign (Lichtenberg)

- **E0.** Static gradient floor (dossier Test 0), read from the t = 0 rows of the baseline gate run: the
  band L2 of |q - 1| and of |zeta| on the exact signed-distance initial fields of every arm and rung.
  This sets `SW_DELTA_S = 0.10 - eps_q` and `GC_EPS = C_eps E_zeta` with C_eps in {0.5, 1, 2}. The gates
  use hex meshes only; other mesh families are out of their scope.
- **E1.** The `exact1D` arm of every gate run, 2D and 3D (closed forms of D2 and D3). The other arms
  start only after it passes.
- **E2.** The `seam` sub-arm of every gate run (2D: np 1, 4, 8 at the coarsest shear rung; 3D: np 16 and
  32). Pass: every CSV column equal to 1e-10 relative (for an Eulerian `R`-type strain weight the psi
  solver tolerance is 1e-14, per STATUS.md section 9).
- **E3.** The 2D gate campaign (section 6). **E4.** The 3D gate on a 2D pass. **E5.** Promotion to
  METHOD.md (section 6).

### Phase F: the Eulerian two-phase line (after the first SL campaign; decision 1)

C3 already opens the extension point in the Eulerian two-phase solver, inert. The Eulerian candidates
enter the coupled gate only after F1 and F2. The Eulerian line uses E1 and E2 (kinematic, 1D and seam)
before that, because those gates do not need the two-phase solver.

- **F1.** Update rho and rhoPhi on every outer corrector (pattern: `slAlphaEqn.H:408-433`). Not inert by
  design. Gates: the translating droplet at density ratio 1 and 840 (the droplet moves at U0), the mass
  residual, the stationary droplet. Before F1 the author decides whether the Eulerian coupled studies are
  void (`sdplsDropletNS2D`, `sdplsDropletMechanism2D`, `sdplsDropletBdf2Droplet2D`,
  `sdplsPsiBudgetDroplet2D`, `sdplsRdivDroplet2D`); if void, rename them `_VOID_frozenRho_<date>`.
- **F2.** Droplet metrics parity: move the droplet metric computation into a function object in
  `src/functionObjects` (field names as dictionary entries), so any solver writes the droplet CSV. Gate:
  the SL droplet CSV stays bit-identical; the Eulerian solver writes the same columns.
- **F3.** The Eulerian line joins the 2D gate (`line: eulerian` candidates).

---

## 5. The method gates

### 5.0 What the gates reuse, and what they replace

| Existing asset | Use in the gates |
|---|---|
| cases `1Dstretch`, `2Dvortex`, `stationaryDroplet2D`, `translatingDroplet2D`, `oscillatingDroplet2D`, `3Dshear`, `stationaryDroplet3D`, `translatingDroplet3D` | the gate arms, unchanged except for new inert tokens; only `oscillatingDroplet3D` is new, because no 3D oscillating case exists |
| `workflow/Snakefile` | runs every arm; the gate adds no rule to it |
| `stationaryLadder2Dshared`, `translatingLadder2D`, `oscillatingLadder2Dshared`, `traceKinematic2Dvortex`, `traceStationary3Dhex` | source of the pinned tokens; the configs stay as the historical record and are not re-run |
| the advection regression set and the WP3 bit-identity set (`docs/plan-library-split-and-build-policy.md`) | the inertness check of a refactor, compared against their preserved baselines; no new config |
| `guard_finished_cases.py`, `foam_log_state.sh`, `compare_metrics_csv.py`, `aggregate.py`, `run-studies.sbatch` | used as they are, by every arm |

What the gates replace (planned earlier as separate studies, now parts of the gates):
- the smoke gate: `SMOKE=1` on the same gate (coarse N, about 20 steps, `profiles/local`);
- the exact 1D checks: the `exact1D` arm, the FIRST arm of BOTH gates (`1Dstretch`, closed form, seconds);
  the other arms start only after it passes, so every gate run carries its own exact check;
- the static gradient metrology: the t = 0 rows of every arm (exact signed distance on the gate meshes);
- the seam checks: the `seam` sub-arm (the coarsest shear rung at np 1 and np 8 next to np 4 in 2D; np 16
  next to np 32 in 3D), reported as the maximum relative CSV difference.

### 5.1 Mechanism: one command per gate

- **Gate definition** `config/gates/<gate>.yaml`: the arms (case, kind, N ladder, T_REF, arm tokens,
  metric instants, solve_runtime), np, the solver per method line, the `fixed` tokens (a copy of
  METHOD.md section 8.1; a candidate cannot change them), the `methodTokens` (the only tokens a candidate
  may set), the verdict thresholds.
- **Candidate** `config/candidates/<name>.yaml`: `line` (`semiLagrangian` or `eulerian`), the method
  tokens (single values), dimensionless rates (`GC_M_MU`, with mu = M_mu/T_REF per arm, constant across
  rungs, never tied to dt), the pre-registered read-out in the header, and a target criterion.
- **On the fly**: `SET="TOKEN=value,..."` on the command line defines an ad-hoc candidate without a file.
  The renderer names it `adhoc-<hash of the tokens>` and marks it `preRegistered: false` in every summary,
  so an exploratory run can never pass for a pre-registered one.
- **Renderer** `workflow/scripts/render_gate_configs.py` writes one complete study config per candidate
  and arm to `studies/<gate>_summary/<candidate>/configs/<arm>.yaml`. It refuses a token that is not a
  method token, a token of the other line, and a collision with a fixed token. It omits `mpi_launcher` and
  `env_preamble`, so the profile supplies them. It prints the token difference between the candidate and
  `baseline`. The merge is done by a committed script with explicit collision errors, so the shallow-merge
  trap of a base config cannot occur.
- **Meta-workflow** `workflow/Snakefile.gate`: `render` -> `run_arm` (per candidate and arm:
  `make studies-one-file CFG=<rendered> PROFILE=<profile>`, which runs the guard and the same `snakemake`
  command line as `make studies-one`) -> `summarize` (per candidate) -> `compare` (candidate vs baseline).
  The outer rules run locally in the orchestrator job; the inner drivers submit to SLURM.
- **Names** [proposal]: `<gate>_<candidate>_<arm>`, for example `methodGate2D_HL1z_translating`; arms
  `shear`, `stationary`, `translating`, `oscillating`.
- **Makefile**: `make gate GATE=methodGate2D CANDIDATES=baseline+HL1z PROFILE=profiles/slurm [PRESERVE=1] [DRYRUN=1]`.
  `baseline` is always added. `PRESERVE=1` renames existing studies of the named candidates with a dated
  suffix before the launch (a rename, never a delete).

### 5.2 2D gate `methodGate2D` (np 4, hex, equidistant)

| Arm | Case | SL-line solver | Eulerian-line solver | N | Cells | h ratio | R/h | END_TIME | Steps |
|---|---|---|---|---|---|---|---|---|---|
| exact1D | `1Dstretch` (uniaxial strain, closed form) | `leiaLevelSetFoam` (`ADVECTION semiLagrangian`) | `leiaLevelSetFoam` (`ADVECTION eulerian`) | as the existing `sdpls1Dstretch` ladder | 1D | | | as `sdpls1Dstretch` | seconds |
| shear | `2Dvortex` (shear2D, reversed, T = 2) | `leiaSemiLagrangeLevelSetFoam` | `leiaLevelSetFoam` (`ADVECTION eulerian`) | 68/96/136 | 4 624 / 9 216 / 18 496 | 1.412, 1.417 | 10.2 / 14.4 / 20.4 | 2 | about 270 / 380 / 550 (CFL 0.5) |
| seam | `2Dvortex` at N = 68 | same | same | 68 at np 1 and np 8 | 4 624 | | 10.2 | 2 | about 270 |
| stationary | `stationaryDroplet2D` (L = 10 mm, R = 1 mm) | `leiaSemiLagrangianLevelSetTwoPhaseFoam` | `leiaLevelSetTwoPhaseFoam` | 100/142/200 | 10 000 / 20 164 / 40 000 | 1.420, 1.408 | 10.0 / 14.2 / 20.0 | 0.1 | 9 207 / 15 580 / 26 042 |
| translating | `translatingDroplet2D` (offset -2.5 mm, U = 0.05 m/s) | same | same | 100/142/200 | same | same | same | 0.1 | same |
| oscillating | `oscillatingDroplet2D` (mode 2, signed-distance ellipse) | same | same | 100/142/200 | same | same | same | 0.1 (about 10.5 periods of 9.51 ms) | same |

- The droplet arms start only after the `exact1D` arm passes (cheapest discriminator first). The summary
  reads the t = 0 row of every arm as the static gradient floor (E0).
- Cell counts double per rung (ratios 1.99-2.02). dt = 10.861 h^1.5 (0.2323 of the Brackbill limit) at
  every rung. All N are even: the droplet centre and the vortex centre sit on mesh vertices at every rung.
- solve_runtime: 240 min (shear), 600 min (droplets; `stationaryLadder2Dshared` ran N = 256 in that limit).
- T_REF: shear 1 s (L/U); stationary 3.7 ms (sqrt(rho R^3/sigma)); translating 20 ms (R/U); oscillating
  1.51 ms (1/omega).
- Fixed tokens, SL line: METHOD.md section 8.1 (`SL_RECONSTRUCTION uncachedQuadraticWeightedLeastSquares`,
  `SL_TRACE_VELOCITY projectedFlux`, `SL_FOOT_INTEGRATOR taylor`, `SL_FIT normalEquations`, `SL_CLIP false`,
  `PSI_FILTER none`, `VOLUME_CORRECTION noVolumeCorrection`, `MASS_FLUX rhoLENT`,
  `PHASE_INDICATOR detrixheAslam`, `SURFACE_TENSION_FORCE reconstructedCurvature`,
  `FACE_CURVATURE_SOURCE model`, `MOMENTUM_DDT_SCHEME backward`, `RHO_DDT_SCHEME backward`),
  `OSCILLATION` on for the shear arm, `CURVATURE_EXTENSION` per arm (see section 9).
- Fixed tokens, Eulerian line: the shared Eulerian discretization of `config/sdplsConv2Dvortex.yaml`
  (`make check-discretization` must pass).
- Cost: the finest droplet rung is about 1.04e9 cell-steps. The `stationaryLadder2Dshared` time limit gives
  an upper bound of about 4 h for it on 4 ranks. The baseline run measures the real number.

### 5.3 3D gate `methodGate3D` (np 32, hex, equidistant, h ratio >= 1.3)

| Arm | Case | N | Cells | h ratio | R/h | END_TIME | Steps |
|---|---|---|---|---|---|---|---|
| exact1D | `1Dstretch` (the same arm as in the 2D gate) | as `sdpls1Dstretch` | 1D | | | as `sdpls1Dstretch` | seconds |
| shear | `3Dshear` (unit cube, R = 0.15, reversed) | 68/90/118 | 314 432 / 729 000 / 1 643 032 | 1.324, 1.311 | 10.2 / 13.5 / 17.7 | 3 (CFL 0.3) | CFL-controlled |
| stationary | `stationaryDroplet3D`, 6R box (L = 6 mm) | 60/78/102 | 216 000 / 474 552 / 1 061 208 | 1.300, 1.308 | 10.0 / 13.0 / 17.0 | 0.025 | 2 302 / 3 412 / 5 102 |
| translating | `translatingDroplet3D`, 6R box | 60/78/102 | same | same | same | 0.02 | 1 842 / 2 730 / 4 082 |
| oscillating | `oscillatingDroplet3D` (NEW, B6) | 60/78/102 | same | same | same | 0.025 (about 2.1 periods of 11.65 ms) | 2 302 / 3 412 / 5 102 |

- Cells grow by 2.20-2.32 per rung. All N are even.
- The finest droplet rung is about 5.4e9 cell-steps per arm. solve_runtime 1440 min. The baseline 3D run
  measures the real core-hours.

### 5.4 Error vector and orders

| Quantity | shear | stationary | translating | oscillating | Instant |
|---|---|---|---|---|---|
| Shape error | `E_GEOM_ALPHA_REL` | `zeroSetRadialL2`/R | `zeroSetRadialL2`/R (about x0 + U t) | period and damping rate from `m2Amplitude`, with the Richardson reference; the Lamb period as a check | T |
| Gradient band error (L2 of |q - 1|, unlimited gradient) | `E_NARROW_L2_GRAD_PSI` (`gradPsiError.csv`) | `gradPsiL2ErrorBand` (after B3) | same | same | T/2 for shear, T for droplets |
| Volume error | `E_VOL_ALPHA_REL` | `phaseVolumeRelError` | same | same | T/2 and T |
| Phase-indicator bounds | `E_BOUND_ALPHA` | `E_BOUND_ALPHA`, `rhoClipFraction` | same | same | maximum over t |
| Spurious currents | none | `maxMagU` (maximum over t, and at T), `l2MagUPrime` | `maxMagUPrime`, `l2MagUPrime` | (the damping rate above) | |
| Pressure-jump error | none | abs(pLaplace - sigma/R)/(sigma/R) | same | none | T |
| Curvature error | none | `kErrL2Band` | same | none | T |
| Travelled fraction | none | none | centroid displacement/(U T) | none | T |
| Cost | `nRanks`, `wallClockSolve`, `solverClockTime` (`ELAPSED_CLOCK_TIME`), `coreSeconds`, `steps`, `secondsPerStep` | same | same | same | |
| State, provenance | COMPLETED or DIVERGED (classifier), `endTimeReached`, `gitCommit`, `libStamps`, `runDate` | same | same | same | |

- `E_BOUND_ALPHA` is 0 by construction for `detrixheAslam` (METHOD.md). The column stays, because the user
  asked for it and because another indicator can be selected.
- No L_inf metric gets an order.
- Orders (`workflow/scripts/richardson.py`): for a metric with an exact value, the two pairwise orders and
  the least-squares order over the three rungs. For a quantity without an exact value (oscillation period,
  damping rate), the procedure of Celik et al., J. Fluids Eng. 130, 078001 (2008): the apparent order p by
  fixed-point iteration for non-constant r, the extrapolated value, GCI_fine with F_s = 1.25, the
  convergence type from sign(eps32/eps21) (monotone, oscillatory, divergent), and the asymptotic-range
  ratio GCI_32/(r21^p GCI_21), which is about 1 in the asymptotic range.
- Outputs: `studies/<gate>_summary/<candidate>/{summary,orders,vsBaseline}.csv` and copies in
  `docs/<theme>/<slug>-article/data/tables/<gate>_<candidate>_{summary,orders,vsBaseline}.csv`.

### 5.5 Verdict, pre-registered in `config/gates/<gate>.yaml` (the thresholds are proposals)

1. Every baseline case is COMPLETED. A candidate case that DIVERGED where the baseline COMPLETED is a FAIL,
   and a result.
2. No regression: every vector metric at the finest rung <= 1.10 x baseline, and every least-squares order
   >= baseline order - 0.3.
3. The candidate's own target from its header (for example: the gradient band error at the finest rung
   <= 0.5 x baseline in every arm).
4. The cost ratio is reported, not scored.

The summary prints PASS or FAIL per criterion. The author reads the whole vector.

### 5.6 Launch on Lichtenberg

```bash
cd /work/scratch/tm83tomy/<clone> && git pull --rebase
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc; . ./etc/leia-env.sh; ./Allwmake > /tmp/leia-build.log 2>&1
make gate GATE=methodGate2D CANDIDATES=baseline+HL1z PROFILE=profiles/slurm DRYRUN=1
sbatch --parsable -J leia-gate2D --export=ALL,TARGET=gate,GATE=methodGate2D,CANDIDATES=baseline+HL1z run-studies.sbatch 2>/dev/null | tail -1 >> .my_jobs
```

Then: add the child job ids from the driver's `.err` to `.my_jobs`; check liveness with
`workflow/scripts/foam_log_state.sh`; verify the binary (`strings` for a new symbol, `grep -m1 '^Exec'` in
one log); pull the results with `make pull-study STUDY=methodGate2D_summary` (only CSVs for 3D). Job names
[proposal]: `leia-gate2D`, `leia-gate3D`.

---

## 6. Campaign: candidates and decision logic

Run order on the SL line. Every run carries `baseline` on the same commit.

| Candidate | Extension | Source | Question |
|---|---|---|---|
| `baseline` | none | none | the matched reference (METHOD.md) |
| `S1` | none | `softWall` | does scalar control alone keep the band? |
| `HL0` | `haloLimited`, `traceFlux extension` | none | how much does the extension alone protect q? (the dossier's first control) |
| `HL1q`, `HL1z` | same | `linearQ` / `linearZ`, M_mu = 1 | the paired q/z comparison (dossier 2: Bq vs Bz) |
| `HL2` | same | `softWall` | band protection when HL1 would need a large mu |
| `FP0` | `closestPoint`, `traceFlux extension` | none | full foot-point extension vs halo limiting |
| later | | `twoThirdsZReg`; `linearZ` at M_mu in {0.5, 2, 5} | after E0 sets eps |

Decision rules (dossier D1-D6, mapped to leia):
- HL0 passes and holds the band: no source is needed ("the source is justified only by measured drift").
- HL1q and HL1z agree within the ladder's GCI: keep `linearZ` (no square root) only if it also costs less;
  else keep `linearQ`.
- S1 is as good as HL1 at lower cost: the extension does not pay for itself; report that.
- FP0 is much better than HL0: improve the direction (`fitFootPoint`) before any source tuning (dossier D6).
- Promotion to METHOD.md needs all of: a 2D gate PASS, a 3D gate PASS, the seam gate E2, the advection
  regression set with the order reported (the change is not inert), the cost ratio. METHOD.md changes in
  the same commit, with the configs and the numbers. A later contradiction starts the retraction protocol
  of CLAUDE.md.

---

## 7. Docs theme and pre-print skeleton

Theme key `gradient-controlled-level-set`, slug `gcls-level-set` (decided). The dossiers stay out of git;
the README cites them:

```
docs/gradient-controlled-level-set/
  README.md                          purpose, source dossiers (title, date, location), status
  gcls-level-set-article/
    gclsLevelSet.tex                 elsarticle; preamble from docs/sdpls-level-set/sdpls-article/sdplsLevelSet.tex
    refs.bib                         seeded from the dossier bibliographies
    data/figures/.gitkeep
    data/tables/.gitkeep
  gcls-level-set-presentation/
    gcls-level-set.template.html     reveal.js, from an existing theme template
```

- Registration: one `_THEMES` row in `workflow/scripts/paths.py`; a Makefile target `article-gcls`.
- Case figures (user request, 2026-09-26): one TikZ/PGF figure per gate case with the geometry, the
  initial condition and the boundary conditions, generated by
  `gcls-level-set-article/figures/make_case_figures.py` from the case files (streamlines, the
  orthographic 3D view and the hidden edges are computed), iterated on rendered pages to publication
  quality; a label that would cross a line gets a small white rectangular background.
- The existing notes (`docs/combined-source-terms/`, `docs/velocity-extension/improvement-metric-footpoint.md`)
  stay where they are; the article cites them.
- Section outline (the source of each section's numbers in brackets):
  1. Introduction: D_t q = -a q; why redistancing is unwanted; the admissible band.
  2. Prior art (dossier 1 section 4): reinitialisation, conservative level set, variational methods,
     Sabelnikov et al., Hamamuki, SDPLS, the bounded-gradient equation, gradient-augmented level set,
     velocity extension, OpenFOAM and LEIA.
  3. Continuum model: D_t q = q (F - a); the halo-limited extension and the interfacial cancellation; the
     law family; the robust inward-pointing estimate.
  4. Discretisation: the material-form operator and the two divergences; the flux-correction form; the
     partition-independent stencil sampler; the exponential source update; the SL realisation.
  5. Implementation: RTS families, composition roots, extension without modification.
  6. Verification [D1-D4, E0, E1, E2 tables].
  7. Coupled results [2D and 3D gate tables with orders and GCI].
  8. Parallel cost and decomposition invariance [cost columns, E2].
  9. Discussion: what was falsified; limits (volume conservation, transition strain).
  10. Conclusions. Appendices: derivations, pseudo-code.
- refs.bib keys: OsherSethian1988, SussmanSmerekaOsher1994, SussmanFatemi1999, AdalsteinssonSethian1999,
  OlssonKreiss2005, OlssonKreissZahedi2007, LiEtAl2010, KeesEtAl2011, SabelnikovEtAl2014, MaricEtAl2015,
  HamamukiNtovoris2016, ToureSoulaimani2016, GibouEtAl2018, Hamamuki2019, QuezadaKuzminKees2019,
  ZhangYue2019, LyrasEtAl2020, KarakusEtAl2022, FrickeEtAl2022, ShaoEtAl2023, LiuEtAl2023,
  NaveRosalesSeibold2010, BockmannVartdal2014, MullerRuggeri1998, BotheFrickeSoga2024, Shakoor2025,
  FerroEtAl2025, BotheSoga2026, WellerEtAl1998, Reitzel2023, LEIA, and Celik2008 (GCI).
- STATUS.md: a new section "11. Halo-limited extension and gradient-control sources (from 2026-09-26)".
  Every gate result is a dated entry with the config, the commit, the binary stamps, the numbers, PASS or
  FAIL against the pre-registered threshold, and any retraction.

---

## 8. Verification, end to end

1. `./Allwmake` passes `etc/leia-check-deps.py`; `ldd` of each solver shows `libleiaGradientControl`, and
   `libleiaVelocityExtension` where it was added.
2. The unit apps exit 0: `leiaTestGradientControlLaw`, `leiaTestSdplsSource`, `leiaTestSlSource`,
   `leiaTestVelocityExtension`, `leiaTestMaterialTransport`; and `python3 workflow/scripts/richardson.py --self-test`.
3. Bit identity after each C step, and after C6 on the final state:
   `compare_metrics_csv.py A.csv B.csv --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME` for every pair
   of the preserved gate set. PASS = every pair identical.
4. The 4-rank laptop smoke: `make gate GATE=methodGate2D SMOKE=1 CANDIDATES=baseline PROFILE=profiles/local`.
   Check `steps=` and the summary columns, not the exit code.
5. Lichtenberg: the gate dry run, then the baseline 2D gate. The summary CSV has every column filled, the
   orders computed, nRanks = 4 from the log, and the wall clock present.
6. The seam gate E2 and the exact 1D gates E1 pass before any candidate enters the 2D gate.

---

## 9. Risks and open author decisions

Open decisions (they do not block Phases A-D; each is needed before the step named):
- **CURVATURE_EXTENSION for the translating and oscillating arms** (before B7). METHOD.md has no value for
  these cases. Proposal: `none` for translating (as `translatingLadder2D`), `cellCentreInverse` for
  oscillating. Record both in the per-case `.parameter` layer and in METHOD.md.
- **Algebraic psi in past oscillating studies (D-e)** (before B4 lands). Void them with `_VOID_algebraicPsi_<date>`,
  or keep them with a caveat in STATUS.md.
- **Frozen rho in the Eulerian coupled studies (D-m)** (before F1). Void them with `_VOID_frozenRho_<date>`,
  or keep them with a caveat.

Risks:
- **`stencilFit` accuracy.** A least-squares fit smooths U. The unit tests measure the order of the
  cancellation error on a circle. If it is below second order, add the `cellPoint` sampler with a
  seam-consistent halo.
- **The halo-limited weight never reaches zero**, so every face samples. The cost is one polynomial
  evaluation per face; the gate records `secondsPerStep`.
- **One `velocityExtension` object per mesh** (registered `Uext`): every composition root constructs
  exactly one.
- **Double cancellation.** `haloLimited` together with `strainWeight full` removes the strain twice. The
  Eulerian consumer warns; the gate candidates never combine the two.
- **Cost.** Sections 5.2 and 5.3 give cell-steps. The baseline gate runs measure the real core-hours before
  any candidate campaign is sized.
