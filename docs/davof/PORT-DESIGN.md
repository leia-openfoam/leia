# The port of DAVOF into leia: delta analysis and design (phase 0, 2026-11-06)

Tomislav Maric, 2026-11-06: "I do not want to make TwoPhaseFlow 'famous' by my work,
I want to have it in my own project." The DAVOF method (the reconstruction, the
swept-face advection, the paraboloid dual advection, the volume rules, the CCS
positioning, the function objects and the apps) moves from his TwoPhaseFlow fork into
leia as leia's own library, built against OpenFOAM v2606 only. The port must reproduce
the fork bit for bit on the verification set. This document is the phase-0 deliverable:
(a) the delta between the fork's VoF base classes and OpenFOAM v2606's geometricVoF,
(b) the library layout, (c) the options dropped, (d) the acceptance test, and the
phase 3 added on 2026-11-06 (the transfer of the running studies). Nothing was coded.

## 0. Sources, versions, machines (measured 2026-11-06)

| what | where | state |
|---|---|---|
| the fork (read only) | laptop `/home/tmaric/OpenFOAM/repos/davof/agents/swept/TwoPhaseFlow`, cluster `/work/scratch/tm83tomy/davof/agents/enright/TwoPhaseFlow` | branch `davof-swept`, HEAD `c5c00c5` (the reference setting `407f969` as defaults) |
| the joint branch | `tmaric/davof-swept-joint` `6a12d7f`, rebased on `407f969`, not on `c5c00c5` | `sweptHexClip` lives in `df125bb` together with the joint-clip variants (8 files, +927 lines; `davofSweptPieces.H` +142, `isoAdvectionSwept.C`) |
| leia | `/home/tmaric/OpenFOAM/repos/davof/leia` (main checkout on `development`, where the preprint folder is); the port works in the worktree `/home/tmaric/OpenFOAM/repos/davof/agents/port/leia`, branch `davof-port` from `davof-swept` `930d2328` | `davof-swept` is 90 commits ahead of `development`, `development` 22 ahead of `davof-swept` (the article folder) |
| OpenFOAM | laptop `/home/tmaric/OpenFOAM/OpenFOAM-v2606`; cluster `$HOME/OpenFOAM/OpenFOAM-v2606` | both the leia and the fork builds on the cluster are stamped `v2606` (`platforms/linux64GccDPInt32Opt/.openfoam-version`); the laptop leia build too |
| the fork's base library | the fork builds against v2606's `finiteVolume` already: it has NO `zoneDistribute` of its own (`find src -name 'zoneDistribute*'`: none), `interpolationCellPoint` and `volPointInterpolation` are OpenFOAM's | the delta is exactly libVoF (fork) against libgeometricVoF (v2606) |
| Eigen 3 (the overlap library) | cluster `$HOME/OpenFOAM/thirdparty/eigen/eigen-3.4.0`, laptop `/usr/include/eigen3` (`etc/davof-env.sh` lines 67-76) | leia needs the same `EIGEN_INC` logic in `etc/leia-env.sh` |
| cluster workspace | `/work/scratch/tm83tomy/davof/agents/port` does not exist yet; created in phase 1 with its own `.my_jobs` | the Enright agent's jobs (`davof-enright-N*`) and the swept agent's jobs are running: not mine, never cancelled |

The DAVOF code in the fork (`wc -l`, `c5c00c5`): `src/VoF/davof/` 11 013 lines (of which
`thirdParty/overlap/overlap.hpp` 1 990, MIT, Strobl, Formella, Poeschel 2016),
`reconstructionSchemes/plicSchemes/davof/` 1 736, the advection member functions inside
`isoAdvection`: `isoAdvectionSwept.C` 2 080, `isoAdvectionCurved.C` 2 742, and the hooks
edited into `isoAdvection.C/.H/Templates.C` (2 862 inserted lines over the fork's master),
the three function objects 1 884, the five apps 2 282 (`advectorVoF` 466 with its include
files, `initDavofFields` 229, `staticFaceFractions` 572, `planeVolumeConsistency` 1 111,
`testSweptPieces` 140). About 22 000 lines move; about 1 500 of them (the hooks in the
plane-clipped loop) do not.

## 1. (a) The delta: the fork's base classes against v2606's geometricVoF

Method: `diff -u` of every class the DAVOF code uses, v2606 against the fork at `c5c00c5`,
and separately v2606 against the fork's `master` (TwoPhaseFlow's own divergence from
OpenFOAM, which is not DAVOF's doing) and `master..davof-swept` (the DAVOF edits). The
diffs are kept in this session's scratchpad; the counts are changed lines.

### 1.1 What v2606 lacks that the DAVOF code uses

| the fork's symbol | used by | v2606 | consequence for the port |
|---|---|---|---|
| `Foam::advectionSchemes` (RTS base of the advection: `advect(Sp, Su)`, `surf()`, `alphaPhi()`, `getRhoPhi`, registered in the object registry as `"advectionSchemes"`) | `advectorVoF`, `areaFractionError`, `volumeFractionError` (`mesh_.foundObject<advectionSchemes>("advectionSchemes")`, `advPtr->surf()`) | absent: v2606 ships `Foam::isoAdvection` alone, no base, no RTS | leia defines its own small base (section 2.3) |
| `isoAdvection` members added for DAVOF: `alphafPtr_`/`alphafOwnPtr_`/`alphafSetFace_`, `davofPtr_`, `planeFluxUpCell_`, `faceClippingPtr_`, the curved-surface state (`cellSurfaces_`, `volumeRulePtr_`, `wispSurfaces_`, `sweptRemote_`, ...), 255 data members against v2606's 29; the functions `initAreaFractions`, `setDefaultAreaFraction`, `defaultAreaFractions`, `publishWetFaces`, `storeAreaFraction`, `sweepSecondRing`, `boundCellContent`, `planeSpeaks`, `planeFlux` (in `isoAdvection.C`), `sweptTriangles`, `sweptFaceFlux`, `sweptFaceCellsFlux`, `writeFaceDump`, `sweptConsistencyProbe`, `packSweptCell`, `sweptExchange`, `classifyPreImage` (`isoAdvectionSwept.C`), `buildCellSurfaces`, `anchorTensorCentroids`, `curvedFlux`, `reportCurved`, `ensureAllPointVelocities`, `ensurePointVelocity` (`isoAdvectionCurved.C`), `syncProcPatches` with the `alphaf` argument, the second-layer bounding in `limitFluxes`, the DAVOF branch of `advect` | the whole advection | absent (v2606's class is closed: private members, no virtuals) | the advection becomes leia's own class carrying its copy of the parts of `isoAdvection` it uses (section 2.3); nothing is derived from v2606's `isoAdvection` |
| `cutFaceAdvect::subFaceAreaEnd()`, `faceSwept()`, `keepChord()`, `chordEnd()`, `edgeCrossings()` (the moving plane's submerged area at `t = dt`, 64 + 105 lines) | `isoAdvection.C` only (13 uses: `planeFlux`, `sweepSecondRing`, `planeSpeaks`), i.e. the plane-clipped loop | absent | not needed: the plane-clipped loop stays in the fork (section 3) |
| `cutCellImpFunc`, `cutFaceImpFunc` | `volumeFractionError.C` line 27 includes `cutCellImpFunc.H` but never uses the class (its exact state is `exactSphere` through `davofInit` or `cutCellIso`, line 89, 159, 195, 360) | absent | the include is dropped; `cutCellIso` exists in v2606 |
| `markInterfaceRegion` | `advectorVoF.C` line 43 includes it; not used | absent | dropped |
| `reconstruction::plicRDF` of the fork (`RDF_` a reference to a shared `reconstructedDistanceFunction`, `centreAndNormalBC` for `alphaContactAngleTwoPhaseFvPatchScalarField`) | the normal model `davofNormalModels::plicRDF` (a helper reconstruction, `plicRDFNormal.C` line 62) and the comparison arm | present as `Foam::reconstruction::plicRDF`, constructed from `(alpha1, phi, U, dict)` the same way | the helper uses v2606's; the numbers of that normal model change where v2606's plicRDF differs (section 1.3); no acceptance case uses it |
| `inline namespace twoPhaseFlow` around every fork base class | all | absent | the DAVOF code refers to `Foam::reconstructionSchemes`, `Foam::cutCellPLIC`, ... unqualified, so nothing changes in the DAVOF sources |
| `reconstructionSchemes::modelDict()` etc. without `noexcept`; `normal_`/`centre_` with `AUTO_WRITE` | `davof.C` | `noexcept` accessors; `writeFields` key (default `false`: `NO_WRITE`) | the davof scheme writes its own fields (`alphaf`, `areaNormalDavof`, ...) itself; where a study reads `normal`/`centre` from a time directory the case sets `writeFields true` (checked per study in phase 1; the converter `tpf_advection_metrics.py` reads the function-object CSVs, not the fields) |
| `implicitFunction` (`libalphaFieldFunctions`) | `davofInit`, `exactNormal`, the three function objects, `initDavofFields`, `staticFaceFractions` | present: `applications/utilities/preProcessing/setAlphaField/alphaFieldFunctions` builds `$(FOAM_LIBBIN)/libalphaFieldFunctions` | the same library (the exact signed distances must be the same code for the identity) |
| `functionObjects::setFlow` | `advectorVoF` | present (`src/functionObjects/field/setFlow`) | the same |
| `immiscibleIncompressibleTwoPhaseMixture` | `advectorVoF/createFields.H` (reads `transportProperties`) | present | the same |
| `zoneDistribute` (`setUpCommforZone`, `updateStencil`, `getStencil`, `getDatafromOtherProc`, `getValue`, `globalNumbering`) | `davof.C`, `davofNormalModel.C`, `ringFitNormal.C`, `isoAdvectionSwept.C` (`sweptExchange`) | the fork USES v2606's (`src/finiteVolume/fvMesh/zoneDistribute`) | no delta |

### 1.2 The shared classes: API identical, three arithmetic differences

`cutCell`, `cutCellPLIC`, `cutCellIso`, `cutFace`, `cutFacePLIC`, `surfaceIteratorPLIC`,
`reconstructionSchemes`: every member the DAVOF code calls (`calcSubCell`, `VolumeOfFluid`,
`cutValue`, `faceArea`, `faceCentre`, `facePoints`, `subCellCentre`, `subCellVolume`,
`cellStatus`, `vofCutCell`, `surfaceArea`, `surfaceCentre`, `subFaceArea`, `subFacePoints`,
`normal()`, `centre()`, `interfaceCell()`, `interfaceLabels()`, `modelDict()`,
`reconstruct()`) exists in v2606 with the same signature (the fork writes the accessors
out of line, v2606 inline `noexcept`). The arithmetic differs in three places, all
TwoPhaseFlow's divergence from OpenFOAM, not DAVOF's:

1. `cutCellPLIC::calcSubCell`: a grazing cut is classified empty or full when
   `mag(faceArea_) < 10*SMALL` (fork, 1e-14) against `< ROOTVSMALL` (v2606, 1e-150).
   A cut polygon with area between 1e-150 and 1e-14 is a status-0 cell with a tiny
   `VOF_` in v2606 and a status +-1 cell with `VOF_` 0 or 1 in the fork. The DAVOF code
   calls `calcSubCell` in the CCS positioning (`davofCcsPositioning.H`), the reconstruction
   (`davof.C`, 7 uses), `planeVolumeConsistency` and `staticFaceFractions`. THE ONE RISK
   TO THE BIT-FOR-BIT IDENTITY in the shared classes. Plan: run the identity with v2606's
   class (the instruction); if a case differs and the first differing cell is such a
   grazing cut, the port carries a 210-line copy `davofCutCellPLIC` (GPL, Scheufler and
   Roenby 2019 cited) with the fork's threshold and the reason in its header; Tomislav
   decides between the copy and accepting the (then documented) difference.
2. `cutCell::calcIsoFacePointsFromEdges`: `mag(faceArea) == 0` (fork) against
   `< VSMALL` (1e-300): no case produces an area in (0, 1e-300).
3. `cutFace::calcSubFaceCentreAndArea`, the triangle branch: the fork writes
   `(1.0/3.0)*(p0 + p1 + p2)` and `0.5*((p1 - p0)^(p2 - p0))`; v2606 calls
   `triPointRef::centre` and `areaNormal`, which are (`triangleI.H` lines 162-190)
   `(1.0/3.0)*(a_ + b_ + c_)` and `0.5*((b_ - a_)^(c_ - a_))`: the same expressions in the
   same order, bit-identical. `sortedOrder` by value against by reference: identical.

### 1.3 The comparison arm: v2606's isoAdvection and plicRDF against the fork's

The Co 1 comparison and the Enright study run plicRDF + isoAdvector. In leia that is
OpenFOAM proper's `Foam::isoAdvection` and `Foam::reconstruction::plicRDF` of v2606,
whose numbers will not equal the fork's rows in `REPORT-swept.md` (the fork's
`isoAdvection` is TwoPhaseFlow's, DLR 2019-2020 on DHI 2016-2017). The differences,
from the v2606-against-fork-master diffs, to be stated with the numbers in phase 2:

- `limitFluxes`: the overshoot tolerance `aTol` is `100*SMALL` (v2606) against `1e-12`
  (fork); `boundFlux`: `100*SMALL` against `10*SMALL`. Which cells are bounded, and when
  the bounding loop stops, differ at round-off; the per-face correction arithmetic is
  the same (`alpha -= dVfcorr/V` in one form, `+= -dVfcorr/V` in the other: identical).
- `timeIntegratedFlux`: v2606 computes `Un0` for all interface cells up front (the same
  interpolation and product), handles upwind cyclic patches, and divides by porosity when
  enabled (off here). No arithmetic change for the cases at hand.
- `applyBruteForceBounding`: `clamp_range` against `min(1, max(0, .))`: identical.
- `plicRDF::interpolateNormal` and `calcResidual`: v2606 fetches neighbour centres with
  `getPosition(..., getCyclicPatches(...))` (a transformation only on cyclic patches; the
  cases have none) where the fork takes `getValue`; `clamp` against `max(min())`: identical;
  the fork adds `centreAndNormalBC` (the contact angle, inert) and holds `RDF_` by
  reference to a shared `reconstructedDistanceFunction`; v2606 owns one.
  `reconstructedDistanceFunction.C` differs by 448 lines, NOT read line by line in phase 0:
  the phase-2 report diffs it where the plicRDF numbers differ beyond round-off.
- The fork's `isoAdvection` holds `alphaPhi_` in the base; v2606's computes `getRhoPhi`
  from `dVf_/deltaT`: the same quantity.

Conclusion for (d): the acceptance cannot ask for the comparison arm's identity with the
fork; it asks for the identity of the DAVOF arm and a stated, explained difference of the
comparison arm (section 4.4).

### 1.4 The DAVOF edits that touch v2606 code paths the port keeps

The port's advection class carries its own copies (not derivations) of these fork
functions, so the identity does not depend on v2606's versions: `advect` (the template,
with the DAVOF branch `curvedFlux` and the mass report), `limitFluxes` (with
`boundSecondLayer`), `boundFlux`, `applyBruteForceBounding`, `syncProcPatches(dVf, phi,
alphaf)`, `setProcessorPatches`, `checkIfOnProcPatch`, `extendMarkedCells`, `netFlux`,
`faceValue`/`setFaceValue`, `setDownwindFaces`. These are the parts of `isoAdvection` the
swept and paraboloid loops use (the fork's `advect` takes `curvedFlux` whenever
`clippingSurface != none`; `timeIntegratedFlux`, the plane-clipped loop, is never
entered on those paths: `isoAdvectionCurved.C` and `isoAdvectionSwept.C` do not call
`planeFlux`, `sweepSecondRing`, `boundCellContent`, `planeSpeaks`, `faceClippingPtr_`
or `advectFace_`; measured by grep, 0 uses). Every copied function keeps the DHI/OpenCFD
/DLR copyright lines of its origin and cites Roenby, Bredmose, Jasak 2016 and Scheufler,
Roenby 2019.

## 2. (b) The layout in leia

### 2.1 One library: `libleiaDavof` under `src/davof` (decision 1 for Tomislav)

leia already has `libleiaDavof` (`src/leiaLevelSet/davof`: `davofState`, the DAVOF state
constructed from a level set or an implicit surface, and the ring-quadric curvature models;
4 671 lines; linked by `leiaTestDavofNormal` and listed in `etc/leia-check-deps.py`). The
port IS the DAVOF method, so it should carry the plain name. Proposal: the port builds
`libleiaDavof` from `src/davof/`, and the existing library is renamed `libleiaDavofState`
in the same phase-1 commit (four mechanical edits: `Allwmake`'s list, `leia-check-deps.py`'s
`LIB`/`LINKS` tables, `leiaTestDavofNormal/Make/options`, `leia-stamp.sh`'s name; the
sources do not move). Alternative if Tomislav prefers no rename: the port is
`libleiaDavofVof`. The design assumes the rename.

The three geometry headers `davofSphereGeometry.H`, `davofSimplexGeometry.H`,
`davofQuadraticFaceGeometry.H` exist in both places and are identical except for one
comment line in the fork ("copied from leia src/leiaLevelSet/davof, commit f13a529"): the
port keeps leia's copies in `src/leiaLevelSet/davof` as the single home and includes them
by path (`-I$(LEIA)/src/leiaLevelSet/davof`, header-only, no link, documented in
`Make/options`); nothing is duplicated, nothing deleted. `davofInit` (the fork, the state
the studies start from: `detrixheAslam`, `quadraticFaces`, `exactSphere`, `areaNormal`,
`wispThreshold`) overlaps `davofState` in purpose but not in code; it moves as is, bit for
bit, and the merge of the two state constructors is a later item, not the port's.

Namespaces and type names stay as in the fork, because the dictionary keys and the
selection names are what the configs render: `Foam::reconstruction::davof` (TypeName
`davof`), `Foam::davofNormalModels::*`, `Foam::davofVolumeRules::*` (keys `volumeRule
kromerBothe | exactSphere | splitQuadrature`), `Foam::davof::{cellSurface, velocityAt,
remoteCell, halfSpace, ...}`, `Foam::davofInit`, `Foam::davofCcsPositioning`.

```
src/davof/                                   libleiaDavof (one library, one RTS table set)
  Make/files, Make/options                   -lgeometricVoF -lfiniteVolume -lmeshTools -lsurfMesh
                                             -lsampling -lalphaFieldFunctions -ltwoPhaseMixture ...;
                                             -isystem $(EIGEN_INC); --no-undefined (leia's rule)
  reconstruction/davof.{C,H}                 from plicSchemes/davof (: Foam::reconstructionSchemes of v2606, RTS)
  reconstruction/davofCcsPositioning.H       the CCS positioning (Maric 2021, JCP:X 11 100093)
  normalModels/{davofNormalModel, cellLocalNormal, youngsNormal, plicRDFNormal, exactNormal,
                ringAreaNormal, ringFitNormal}.{C,H}
  volumeRules/{davofVolumeRule, kromerBotheVolumeRule, exactSphereVolumeRule,
               splitQuadratureVolumeRule}.{C,H}, davofKromerBothe.H
  thirdParty/overlap/{overlap.hpp, LICENSE, README.md, CITATION.cff}   (MIT, unchanged)
  geometry/{davofSweptPieces, davofSweptRemote, davofVelocityAt, davofOrbit,
            davofCurvedGeometry}.H
  init/davofInit.{C,H}
  advection/davofAdvectionScheme.{C,H}       the base (section 2.3)
  advection/davofSweptAdvection.{C,H}        ONE class: the dual advection with its own flux assembly
  advection/davofSweptAdvectionFlux.C        from isoAdvectionSwept.C (sweptFaceCells, the exchange)
  advection/davofSweptAdvectionSurfaces.C    from isoAdvectionCurved.C (cell surfaces, anchor, curvedFlux, prism)
  advection/davofSweptAdvectionTemplates.C   advect, limitFluxes, boundFlux (the fork's versions)
  advection/isoAdvectorScheme.{C,H}          the comparison arm: a wrapper of v2606's Foam::isoAdvection
src/functionObjects/davof/                   libleiaDavofFunctionObjects (-lleiaDavof)
  areaFractionError/, volumeFractionError/, reconstructionError/
applications/solvers/davofAdvectorFoam/      from apps/benchmark/advectorVoF (+ its four include files)
applications/utilities/preProcessing/davofInitFields/     from apps/initDavofFields
applications/test/davofStaticFaceFractions/              from apps/benchmark/staticFaceFractions
applications/test/davofPlaneVolumeConsistency/           from apps/benchmark/planeVolumeConsistency
applications/test/davofTestSweptPieces/                  from apps/test/testSweptPieces
```

The apps take names distinct from the fork's binaries on purpose: `etc/leia-env.sh`
appends the fork's `platforms/.../bin` (`LEIA_TPF_DIR`) after leia's own, so both routes
coexist on one `PATH` for the acceptance comparison; a same-named `initDavofFields` would
shadow the fork's silently. The function objects keep their type names
(`areaFractionError`, `volumeFractionError`, `reconstructionError`): they are selected per
case through the `libs` entry, and a leia case never loads `libpostProcess.so` (loading
it would pull the fork's `libVoF` with a second `reconstructionSchemes` RTS table into
the same process).

### 2.2 Headers, citations, licence

leia's `LICENSE.md` is the GPL v3. Every ported file gets leia's header form (as
`davofState.C`: `leia :: davof :: <name>`) with `Copyright (C) 2026 Tomislav Maric, TU
Darmstadt` and the GPL block; files copied from OpenFOAM/TwoPhaseFlow keep their original
copyright lines above his (DHI 2016-2017, OpenCFD, DLR 2019-2020, Johan Roenby) and cite
Roenby, Bredmose, Jasak, R. Soc. Open Sci. 3 (2016) 160405, doi:10.1098/rsos.160405 and
Scheufler, Roenby, J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009 where
geometricVoF's classes are used (the reconstruction, the advection base, the function
objects' use of `cutCellIso`). The algorithm citations already in the fork's headers
(Kromer-Bothe 2023 doi:10.1016/j.jcp.2022.111840; Strobl-Formella-Poeschel 2016
doi:10.1016/j.jcp.2016.02.003; Maric-Marschall-Bothe 2018 doi:10.1016/j.jcp.2018.03.048
(UFVFC); Maric 2022 CPC 273 108249 doi:10.1016/j.cpc.2021.108249 (triSurfaceImmersion);
Maric 2021 JCP:X 11 100093 doi:10.1016/j.jcpx.2021.100093 (CCS); Detrixhe-Aslam;
Dowell-Jarratt (Illinois); Rodrigues 1840) move with the code and are checked file by file
in phase 1 against the rule of `feedback-full-citation-in-every-implementation` (authors,
title, journal, volume, year, pages, DOI; arXiv id with the journal DOI in one citation).
The overlap library keeps its MIT `LICENSE` and `CITATION.cff` next to `overlap.hpp`.

### 2.3 The advection class

`Foam::davofAdvectionScheme` (`advection/davofAdvectionScheme.H`): an `IOdictionary`
(`regIOobject` name `davofAdvectionScheme`, so a function object finds it by
`foundObject<davofAdvectionScheme>("davofAdvectionScheme")`), constructed from `(alpha1,
phi, U)`, holding `surf_` (`reconstructionSchemes::New(alpha1, phi, U, dict)` of v2606,
so every row still selects `reconstructionScheme davof | plicRDF | gradAlpha | isoAlpha`),
with the RTS table on the key `advectionScheme` of the alpha solver dictionary and the
virtuals `advect(const volScalarField::Internal& Sp, const volScalarField::Internal& Su)`,
`surf()`, `alphaPhi()`, `getRhoPhi(rho1, rho2)`, `advectionTime()`. Two concrete schemes:

- `davofSwept` (`davofSweptAdvection`): the dual advection. Dictionary keys as in the fork
  (section 3 lists what is dropped), with `surface plane | paraboloid | exactSphere`
  replacing `clippingSurface` (`none` no longer exists: the plane-clipped loop is gone; the
  fork's `clippingSurface` is accepted as an alias of `surface` for one release so that
  today's templates render) and `flux sweptFaceCells | prism` replacing `fluxIntegration`
  (same alias rule). The class holds the state of the fork's `isoAdvection` that the two
  loops use (section 1.4), the alphaf storage (`alphaf.<group>` registered, `alphafSetFace_`,
  the wet-face list for the davof scheme), `volumeRulePtr_`, `cellSurfaces_`,
  `wispSurfaces_`, `sweptRemote_`, the caches of the round-4 keys, the diagnostics kept
  (section 3). Its `advect` is the fork's `advect` with the `clippingSurface != none`
  branch made unconditional. It requires `surf()` to be `reconstruction::davof` (the fork's
  `davofPtr_`, a `dynamic_cast` at construction with a `FatalError` otherwise).
- `isoAdvector` (`isoAdvectorScheme`): owns a `Foam::isoAdvection` of v2606 and forwards
  `advect`, `surf()`, `alphaPhi()` (from `getRhoPhi` with rho 1, 0, or `dVf/deltaT`),
  `advectionTime()`. The comparison arm (plicRDF + isoAdvector) and nothing else.

The solver `davofAdvectorFoam` is `advectorVoF` with `advectionSchemes::New` replaced by
`davofAdvectionScheme::New`, the unused `markInterfaceRegion` include removed, the abort
guard, `setFlow`, `period`/`reverseTime`/`spirallingFlow` and the `isoAdvector` dictionary
kept as they are (the templates render them).

### 2.4 The workflow

The study route is a token, not a code path: `solver: advectorVoF` (the fork, kept for
the acceptance comparison and nothing else afterwards) or `solver: davofAdvectorFoam`
(leia). `workflow/Snakefile` treats both as the "static" advection app (`_STATIC_APPS`,
the `post_solve` converter `tpf_advection_metrics.py`, the `rm -rf postProcessing`
before a solve). `tpf_rows.py` gains one column per row, the advection scheme of the leia
route (`davofSwept` for every davof row, `isoAdvector` for `plicRDF`/`isoAlpha`), rendered
through a new token `DAVOF_ADVECTION_SCHEME` whose fork-route value is `isoAdvection`;
`setfields_app: initDavofFields` becomes `davofInitFields` on the leia route (a token
too). The case templates (`cases/davof/{sphereRotation3D, sphereTranslation3D,
planeTranslation3D, sphereDeformation3D}/system/{fvSolution, controlDict}.template`)
change in exactly three places: the `advectionScheme` line, the `application` line and the
function objects' `libs` entries (`libpostProcess.so` against
`libleiaDavofFunctionObjects.so`), each a token. Everything else renders identically on
both routes (the phase-3 criterion). The scoring and comparison scripts under
`workflow/scripts/` (`tpf_advection_metrics.py`, `davof_cmp/`, `davof_prof/`,
`davof_swept_ladders.py`, `make_davof_advection_table.py`, the rescore scripts) read
CSVs and logs, not binaries; the runners that name a binary (`davof_cmp`'s timed runner,
`davof_prof`'s profile runners, the enright runner) take the app name from the same token.
`etc/leia-env.sh` gains the `EIGEN_INC` default of the fork's `davof-env.sh` and keeps
`LEIA_TPF_DIR` (the fork route).

## 3. (c) What is dropped, what is kept

Kept (the production path and what the studies rely on): `surface plane` with
`planeSource anchor | exact`, `surface paraboloid` (the ring-tensor paraboloid, `fitAnchor`
`faces | volume`, `fitDof six | one | none`, `tensorCentroid`, `tensorSource ring`,
`fitPriorWeight/fitTol/fitMaxIter`, the curvature tensor of the davof scheme that the
paraboloid prior needs), `surface exactSphere` (the control), `flux sweptFaceCells` with
every reference-setting key (`sweptWispSurface/-Tol/-MatchTol`, `sweptConsistent`,
`sweptVolumeCorrection`, `positioning ccs | vofCutCell`, `positioningTol/-MaxIter`,
`anchorSkipConverged`, `sweptPreTest`, `sweptTraceCache`, `sweptCapSet`, `sweptCoplanar`,
`sweptCellsFaceUpdate`, `sweptScaleToPhi`, `sweptPlaneTol`, `anchorMethod ccs | illinois |
newton | compare`, `anchorTolRel`), the four round-4 bit-for-bit keys (`sweptAxisPlanes`,
`sweptTetCache`, `sweptPointVelocityLazy`, `sweptPreTestCoplanarSide`), the two round-off
opt-ins (`sweptFaceUpdateFromVolume`, `sweptSliverClosedForm`), `sweptHexClip` from the
joint branch, `flux prism` (the paraboloid's flux and the Kromer-Bothe rule), the volume
rules `kromerBothe` (+ `kromerBotheFan`, `kromerBotheMinSplit`, `kromerBothePrismApex`,
`volumeRuleFallbackTolRel/-MaxDepth`), `exactSphere`, `splitQuadrature` (the capped
fallback of kromerBothe, `curvedTolRel`, `curvedMaxDepth`), the normal models (all seven),
`faceState scalar | vector` (the warped-face state the reconstruction reads), `boundSecondLayer`,
the bounding keys (`nAlphaBounds`, `snapTol`, `clip`, `isoFaceTol`, `surfCellTol`,
`wispTol`), and the diagnostics the studies and the profiles read: `sweptFaceLog`,
`sweptClassStats`, `sweptClassify`, `sweptPreTestCheck`, `sweptConsistencyLog`,
`sweptFaceDump`, `sweptDebugPoint`, `anchorDiagnostics`, `volumeRuleDebug`, `curvedDebug`,
`defectEvery/-Levels`, `writeCurvedFields`, `writeCurvatureTensor`, the mass report and the
parallel wall-time line (both printed by `advect`), the abort guard.

Dropped, with the reason:

| option | reason |
|---|---|
| the plane-clipped loop (`clippingSurface none`: `timeIntegratedFlux` with `planeFlux`, `sweepSecondRing`, `boundCellContent`, `planeSpeaks`, `contentBound`, `defaultAreaFractions`' second-ring part, the `cutFaceAdvect` additions) | the preprint v1 baseline; stays in the fork at `davof-v1.0` and `c5c00c5`; the port's `advect` is the curved branch alone |
| the face-clipping models `faceClipping plane | planeSagitta` (`davofFaceClipping`, `planeClipping`, `planeSagittaClipping`, the keys `sagittaMinCos/-MinChord/-MaxFraction/-Form`, `writeSagittaCorrection`) | they correct the plane-clipped loop's face update (`faceClippingPtr_` has 14 uses in `isoAdvection.C`, 0 in the swept and curved files); the curvature tensor they share with the paraboloid stays in the reconstruction (`updateCurvatureTensor`, `requestCurvatureTensor`), requested by the paraboloid alone |
| `fluxIntegration sweptFace` (one donor) | falsified on the rotation 2026-10-29 (every pre-image leaves the donor; unbounded at N 64); not a fallback of anything kept |
| `fluxIntegration trapezoidal | simpson | simpsonSplit | simpsonSplitS | gaussSplit` | the time-quadrature fluxes of branch `davof-fluxint` (2026-10-26): the first three fail (vertex passages); `gaussSplit` holds the plane gate at 1.5-1.8x the prism's cost and is blocked by the one-surface-per-face selection; none is the tested fallback of `prism` or `sweptFaceCells` |
| `sweptStateShortcut` | "not achieved" (README 2026-11-04); broke the plane gate when defaulted on in `3aaf49d` |
| `sweptNoFluxLayer` | an opt-in, default `false` in `c5c00c5` (`isoAdvection.C` line 225), so no study and no acceptance key depends on it; round-4 gain 1 (`fbbfaf9`) is a path the reference setting does not take |
| `sweptDonorFastPath` | "useless under rotation (vertices always leave laterally)", round 3 |
| `sweptClip joint | hybrid | jointGuarded | jointConfig`, `sweptDonorCheck`, `sweptJointDumpSteps` | negative result 2026-11-05 (the flux must be single-valued per region); stay in the fork's joint branch |
| `sweptHexClip` (the joint branch, `df125bb`) | NOT a second key: it is the same device as `c5c00c5`'s round-4 key `sweptAxisPlanes` (gain 2, `fbbfaf9`..`c5c00c5`: `halfSpace::axis`/`sgn`, `setAxis`, `dist` = the one coordinate difference, bit for bit the dot product, applied to the local and the remote cell planes), implemented independently on the `407f969` base by the joint agent (`halfSpace::axis`, `axisOf`, `sdist`, `remoteCellPlanes(.., hexClip)`); found in phase 1 (2026-11-06) when extracting it. The port carries `sweptAxisPlanes`; the acceptance item 4.1 "sweptHexClip on against off" is run as `sweptAxisPlanes` on against off |
| `faceClipping`-dependent rows in `tpf_rows.py` (`davof-sagitta-*`), the `clippingSurface none` rows (`davof-cellLocal`, `-ringArea`, `-ringFit`, `-exact`, `-plicRDFnormal*`) | they select the dropped loop; the rows stay in `tpf_rows.py` marked fork-only (the fork route renders them, the leia route refuses them with a message) |
| `apps/benchmark/reconstructInterface`, `apps/test/testCutFaceAdvectEndArea`, `src/postProcessing/{interface, interfaceRegion, superHeated, interfaceEnergyFluxes, wallHeatFluxSolid}` | not DAVOF; `testCutFaceAdvectEndArea` tests the dropped `cutFaceAdvect` additions |
| `markInterfaceRegion`, `cutCellImpFunc` includes | unused |

## 4. (d) The acceptance test (phase 2, on Lichtenberg, workspace `agents/port`)

The fork side of every comparison is the fork at `c5c00c5` built in the port workspace
(`agents/port/TwoPhaseFlow`, a clone, never committed to), the leia side the `davof-port`
branch built against the same OpenFOAM v2606; one `env.sh` sources both as the swept
workspace does. Every run is a leia study (`config/davof/port_*.yaml`, `studies/davof/`)
submitted through the workflow, so that the case rendering is the same on both routes
except the three tokenised lines. "Bit for bit" means: the written fields of the last
time (`alpha.water`, `alphaf.water`, `areaNormalDavof` where written) with
`writePrecision 17` are byte-identical (`cmp`), the function objects' CSVs are
byte-identical, and the per-step `Info` numbers of the log (the mass report, the
min/max lines, the fallback counts) are identical as text.

### 4.1 The identity set (eight cases, serial, 20 steps, the reference setting + four round-4 keys)

| case | arm | N | row | state | Co |
|---|---|---|---|---|---|
| sphereRotation3D | rotation | 32 | davof-planeSurf-ringArea | exactSphere | 0.2 |
| sphereRotation3D | rotation | 64 | davof-planeSurf-ringArea | exactSphere | 0.2 |
| sphereRotation3D | lineDiagXZ | 32 | davof-planeSurf-ringArea | exactSphere | 0.2 |
| planeTranslation3D | normal | 40 | davof-planeSurf-exact | detrixheAslam | 0.5 |
| planeTranslation3D | normal / oblique | 20, 40 | davof-planeSurf-ringArea | detrixheAslam | 0.5 |

Keys on both sides: the reference defaults of `407f969` plus `sweptAxisPlanes true`,
`sweptTetCache true`, `sweptPointVelocityLazy true`, `sweptPreTestCoplanarSide true`
(written explicitly, since the fork's defaults for the four are `false`). Criterion: all
eight bit for bit. Then `sweptAxisPlanes false` on the leia side only (the other three round-4 keys on): all
eight bit for bit against the run with it on (the fork's claim for gain 2; the joint
branch's `sweptHexClip` is the same device, section 3). A difference in any case is traced
to its first differing cell and step (`sweptFaceLog` and the per-cell print exist for
that) and explained in the report before anything else runs; section 1.2 item 1 is the
first suspect.

### 4.2 The parallel identity

`sphereRotation3D` rotation N 32, ring row, `END_TIME 0.05`, `writePrecision 17`, np 1, 2,
4, 8, 16 (the configs `sphereRotation3D_hexSweptConsParN32np*` with the port tokens).
Criterion: leia at np k bit for bit against the fork at np k with the same decomposition
(`decomposePar` of the same case gives the same `processor*`); leia np k against leia np 1
within 1e-12 of L1(alphaf), as the fork's parallel was validated (`REPORT-swept.md`, the
par merge of 2026-11-02).

### 4.3 The ladders

`sphereRotation3D_hexRotFinal` (rotation, N 32/64 serial; N 128 on 64 ranks, the config
`sphereRotation3D_hexRotFinalN128`) and `sphereTranslation3D_hexSweptLadderWispM9` (N 32/64
serial; N 128 on 64 ranks), rows planeSurf-exact and planeSurf-ringArea, with the port
tokens. Criterion: the scored numbers (`davof_swept_ladders.py`: L1(alphaf), L1(n), L1(x),
E_shape, boundedness, mass) identical in every printed digit to `csv/ref_setting.csv`
(2026-11-05) in serial, within 1e-10 relative in parallel (the fork's own parallel scatter).

### 4.4 The Co 1 comparison rerun

The 18 cells of 2026-11-03 (`sphereRotation3D_hexCmp{Rot, LineX, LineDiag}{Swept, RDF}N{32,
64, 128}`, OF Co 1 and 0.2, same state, mesh, dt and ranks, exclusive nodes, three
repetitions of the timed cells) with the leia apps: the DAVOF arm against the fork's
numbers (identical digits serial, 1e-10 parallel), the comparison arm `plicRDF +
isoAdvector` of v2606 against the fork's plicRDF + TwoPhaseFlow-isoAdvector rows: every
difference stated with the version (v2606 `geometricVoF` against TwoPhaseFlow master's
`libVoF`), the settings (`nAlphaBounds 5`, `snapTol 0`, `clip false`, `interpolateNormal`),
and the numbers, with section 1.3's list as the explanation; a difference beyond what the
tolerance change in `limitFluxes` explains is diffed down to the line
(`reconstructedDistanceFunction.C`). Wall clock at equal ranks, three repetitions, the
medians of the per-step mean (`davof_cmp`'s timed runner), reported next to the fork's
(`cmp-time`), with the expectation that the leia binary is within the run-to-run scatter
of the fork's (same code, same compiler flags).

### 4.5 Build and gate on both machines

`./Allwmake` of leia under v2606 only (the gate: `grep -r TwoPhaseFlow src/davof
applications/*/davof* src/functionObjects/davof` returns nothing; no `-I` or `-l` to the
fork; `ldd` of every leia binary shows no `libVoF`, no `libpostProcess`); the laptop does
compile checks only; every run on the cluster. The fork's testsuite is not rerun (the fork
is read-only).

## 5. Phase 3: the transfer of the running studies (added 2026-11-06)

(1) Every improvement on the fork's branches comes over: the round-4 keys (in the
acceptance already), the two round-off opt-ins `sweptFaceUpdateFromVolume` and
`sweptSliverClosedForm` as opt-in keys (verified by one rotation N 64 run each against the
fork with the key on: identical digits), `sweptHexClip` (section 4.1), the diagnostics
of section 3 (kept). The joint-clip variants stay in the fork.

(2) The leia study route switches to the leia apps: every config under `config/davof` and
every template under `cases/davof` renders for the leia app through the tokens of section
2.4 (`DAVOF_ADVECTION_SCHEME`, the app name, the `libs` entry); `workflow/scripts/`'s
runners point at the leia binaries; the fork route stays selectable by `solver:
advectorVoF` for the acceptance comparison only. Criterion: for every config, the case
rendered on the leia route differs from the case rendered on the fork route in exactly
the `advectionScheme` line of `fvSolution`, the `application` line and the `libs` entries
of `controlDict`, and in nothing else (`diff -r` of two renderings, scripted over all
configs under `config/davof`). The coordinator's criterion names the solver line alone;
the `libs` entries must differ as well (section 2.1, last paragraph), so the criterion
is stated with the three lines.

(3) The studies running in the fork are rerun in leia, in this order: the 3D Enright
deformation first (the Enright agent's `cases/davof/sphereDeformation3D` and the six
configs `sphereDeformation3D_hex{Swept, RDF}N{32, 64, 128}` on leia branch
`davof-swept-enright`, fork branch `davof-swept-enright`: taken over unchanged except the
tokenised lines; the DAVOF arm's numbers reproduced), then the reference-setting ladders
and the Co 1 comparison as the acceptance already does. Criterion for the Enright
transfer: the DAVOF arm's scored numbers (E_shape at T, L1(alphaf), boundedness, mass,
the per-step CSV) identical within round-off to the fork's (identical digits serial,
1e-10 relative on 8 and 64 ranks); the comparison arm's difference stated as in 4.4.
The Enright agent's work is not touched; its results in the fork are the reference.

(4) From the port on, new results are recorded against leia commits (the version stamp of
`libleiaDavof` through `etc/leia-stamp.sh`, printed by every app and written to
`<case>/leia.version`, replaces the fork's commit in the CSVs' version column); the
fork's results stay as the reference in `REPORT-swept.md`, under the heading "The port
into leia (2026-11-06)" appended in phase 2.

## 6. Effort estimate

| phase | work | estimate |
|---|---|---|
| 1, the library | move 22 000 lines, rename the base class and remove the plane-clipped hooks (about 1 500 lines in `isoAdvection.C/.H/Templates.C` to disentangle by hand, the one place where new code is written), the two scheme classes, the function-object and app Make files, the `EIGEN_INC` logic, the rename of the state library, the citation pass over every header, the workflow tokens, compile on the laptop, build on the cluster | 2 to 3 agent working days, the disentangling of `advect`/`limitFluxes`/the alphaf storage being the uncertain part (half a day to a day) |
| 1, `sweptHexClip` | extract from `df125bb` by hand | half a day, verified by 4.1 |
| 2, acceptance | 8 + 5 + 6 + 18 x 3 runs; the N 128 ladders about 13 min each on 64 ranks; the queue is the unknown | 1 day of runs, half a day of analysis; one more day if the identity fails and section 1.2 item 1 must be settled |
| 3, transfer | the tokens over all configs and templates, the renderer diff, the Enright rerun (N 128 on 64 ranks, T = 3: hours per cell), the ladders | 1 to 2 days plus the run time |
| total | | 5 to 7 agent days, of which 2 to 3 are waiting on the cluster |

## 7. Decisions for Tomislav (asked at the phase-0 hand-back)

1. The library name: `libleiaDavof` for the port with the existing `src/leiaLevelSet/davof`
   renamed `libleiaDavofState` (recommended), or `libleiaDavofVof` with no rename.
2. The grazing-cut threshold of `cutCellPLIC` (section 1.2 item 1): accept v2606's and
   document a difference if one appears, or carry the 210-line copy with the fork's
   threshold. The design runs the identity first and asks only if it fails there.
3. The app names (`davofAdvectorFoam`, `davofInitFields`, `davofStaticFaceFractions`,
   `davofPlaneVolumeConsistency`, `davofTestSweptPieces`) and the keys `surface` and `flux`
   with the fork's `clippingSurface`/`fluxIntegration` as aliases for one release.
4. The push of `davof-port` to `origin`: attempted at the end of phase 0; if the
   permission classifier refuses it (as it refused the tag pushes of 2026-10-27), Tomislav
   pushes.
