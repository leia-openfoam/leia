/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
    Modified code Copyright (C) 2019 Johan Roenby
    Modified code Copyright (C) 2019-2020 DLR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "davofSweptAdvection.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "syncTools.H"
#include "profiling.H"
#include "davof.H"
#include <cstring>

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(davofSweptAdvection, 0);
    addToRunTimeSelectionTable(davofAdvectionScheme,davofSweptAdvection, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofSweptAdvection::davofSweptAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    davofAdvectionScheme
    (
        typeName,
        alpha1,
        phi,
        U
    ),
    // The reconstruction scheme, as the fork's base class constructed it:
    // before every other member, from the alpha solver dictionary
    surf_(reconstructionSchemes::New(alpha1, phi, U, modelDict())),
    // General data
    mesh_(alpha1.mesh()),
    alpha1In_(alpha1.ref()),
    dVf_
    (
        IOobject
        (
            "dVf_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol, Zero)
    ),
    advectionTime_(0),
    timeIndex_(-1),

    // Tolerances and solution controls
    nAlphaBounds_(modelDict().lookupOrDefault<label>("nAlphaBounds", 3)),
    isoFaceTol_(modelDict().lookupOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().lookupOrDefault<scalar>("surfCellTol", 1e-8)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0),

    // DAVOF face fractions
    storeAreaFractions_
    (
        modelDict().getOrDefault<bool>
        (
            "storeAreaFractions",
            surf_->type() == "davof"
        )
    ),
    alphafPtr_(nullptr),
    alphafOwnPtr_(nullptr),
    alphafSetFace_(),
    alphafSetFaces_(0),
    alphafPrevSetFaces_(0),
    alphafFallbackFaces_(0),
    defaultsValid_(false),
    alphaDefault_(0),
    phiPos_(),
    phiEventNo_(-1),
    wetFace_(),
    wetFaces_(0),
    davofPtr_(nullptr),
    nAlphafClipped_(0),
    boundSecondLayer_(modelDict().getOrDefault<bool>("boundSecondLayer", true)),
    nBoundedBeyondRing_(0),
    clippingSurface_
    (
        // the key surface (plane | paraboloid | exactSphere); the fork's
        // clippingSurface accepted as an alias (deprecated, printed below)
        modelDict().found("surface")
      ? modelDict().get<word>("surface")
      : modelDict().getOrDefault<word>("clippingSurface", "plane")
    ),
    faceState_(modelDict().getOrDefault<word>("faceState", "scalar")),
    wetVecPtr_(nullptr),
    tensorCentroid_(modelDict().getOrDefault<word>("tensorCentroid", "cut")),
    nCentroidFail_(0),
    curvedTolRel_(modelDict().getOrDefault<scalar>("curvedTolRel", 1e-11)),
    curvedMaxDepth_(modelDict().getOrDefault<label>("curvedMaxDepth", 12)),
    fitPriorWeight_(modelDict().getOrDefault<scalar>("fitPriorWeight", 1e-2)),
    fitTol_(modelDict().getOrDefault<scalar>("fitTol", 1e-10)),
    fitMaxIter_(modelDict().getOrDefault<label>("fitMaxIter", 50)),
    fitAnchor_(modelDict().getOrDefault<word>("fitAnchor", "faces")),
    fitDof_(modelDict().getOrDefault<word>("fitDof", "six")),
    anchorTolRel_(modelDict().getOrDefault<scalar>("anchorTolRel", 1e-11)),
    defectEvery_(modelDict().getOrDefault<label>("defectEvery", 10)),
    curvedDebug_(modelDict().getOrDefault<bool>("curvedDebug", false)),
    anchorDiagnostics_(modelDict().getOrDefault<bool>("anchorDiagnostics", false)),
    anchorMethod_(modelDict().getOrDefault<word>("anchorMethod", "ccs")),
    anchorSkipConverged_(modelDict().getOrDefault<bool>("anchorSkipConverged", true)),
    fluxIntegration_
    (
        // the key flux (sweptFaceCells | prism); the fork's fluxIntegration
        // accepted as an alias (deprecated, printed below)
        modelDict().found("flux")
      ? modelDict().get<word>("flux")
      : modelDict().getOrDefault<word>("fluxIntegration", "prism")
    ),
    sweptPlaneTol_(modelDict().getOrDefault<scalar>("sweptPlaneTol", 1e-6)),
    nSweptFaces_(0),
    nSweptFlips_(0),
    nSweptOutOfDonor_(0),
    nSweptClipTris_(0),
    sweptClosureMax_(0),
    sweptClassify_(modelDict().getOrDefault<bool>("sweptClassify", false)),
    sweptClass_(),
    nSweptClass_(label(0)),
    sweptEventReported_(false),
    dVfBeforeBound_(),
    nSweptPieces_(0),
    nSweptBeyond_(0),
    sweptBeyondMax_(0),
    nSweptRankFallback_(0),
    nSweptPartialCells_(0),
    sweptPieceClosureMax_(0),
    nSweptDegenerate_(0),
    sweptCellsFaceUpdate_(modelDict().getOrDefault<word>("sweptCellsFaceUpdate", "pieces")),
    sweptCoplanar_(modelDict().getOrDefault<word>("sweptCoplanar", "donorSide")),
    nSweptCoplanar_(0),
    sweptPreTest_(modelDict().getOrDefault<word>("sweptPreTest", "touching")),
    nSweptPreSkipped_(0),
    nSweptPreClipped_(0),
    sweptStateShortcut_(modelDict().getOrDefault<bool>("sweptStateShortcut", false)),
    sweptPreTestCheck_(modelDict().getOrDefault<bool>("sweptPreTestCheck", false)),
    nSweptPreViol_(0),
    nSweptPreSlab_(0),
    nSweptShortViol_(0),
    nSweptShortNotInside_(0),
    cellBoxes_(),
    sweptTraceCache_(modelDict().getOrDefault<bool>("sweptTraceCache", true)),
    traceCell_(),
    traceU_(),
    traceDone_(),
    nSweptTraceHits_(0),
    nSweptTraceMisses_(0),
    nSweptShortFull_(0),
    nSweptShortEmpty_(0),
    nSweptShortPieces_(0),
    nSweptShortFullTol_(0),
    nSweptShortEmptyTol_(0),
    sweptClassStats_(modelDict().getOrDefault<bool>("sweptClassStats", false)),
    sweptDonorFastPath_(modelDict().getOrDefault<bool>("sweptDonorFastPath", false)),
    nSweptDonorFast_(0),
    sweptFaceUpdateFromVolume_(modelDict().getOrDefault<bool>("sweptFaceUpdateFromVolume", false)),
    sweptPreTestCoplanarSide_(modelDict().getOrDefault<bool>("sweptPreTestCoplanarSide", false)),
    sweptNoFluxLayer_(modelDict().getOrDefault<bool>("sweptNoFluxLayer", false)),
    sweptAxisPlanes_(modelDict().getOrDefault<bool>("sweptAxisPlanes", false)),
    sweptPointVelocityLazy_(modelDict().getOrDefault<bool>("sweptPointVelocityLazy", false)),
    sweptTetCache_(modelDict().getOrDefault<bool>("sweptTetCache", false)),
    sweptSliverClosedForm_(modelDict().getOrDefault<bool>("sweptSliverClosedForm", false)),
    nSweptSliverClosed_(0),
    sweptExactCaps_(modelDict().getOrDefault<bool>("sweptExactCaps", false)),
    sweptReinitFaceFractions_(modelDict().getOrDefault<bool>("sweptReinitFaceFractions", false)),
    sweptFragmentLog_(modelDict().getOrDefault<bool>("sweptFragmentLog", false)),
    fragFallbackFaces_(),
    fragFallbackFacesPrev_(),
    fragBound2Cells_(),
    fragBound2CellsPrev_(),
    fragWispCells_(),
    fragWispCellsPrev_(),
    fragLogged_(),
    nFragLogged_(0),
    nFragFallback_(0),
    nFragBound2_(0),
    nFragWisp_(0),
    nFragNone_(0),
    fragStream_(),
    nSweptExactUsed_(0),
    nSweptExactChecked_(0),
    sweptExactDiffMax_(0),
    sweptCoverageCheck_(modelDict().getOrDefault<bool>("sweptCoverageCheck", false)),
    nSweptCovChecked_(0),
    sweptCovClosureMax_(0),
    sweptCovExactMax_(0),
    sweptPieceDumpTol_(modelDict().getOrDefault<scalar>("sweptPieceDumpTol", -1)),
    nSweptPieceDumps_(0),
    fragBoundFaces_(),
    fragBoundFacesPrev_(),
    fragBoundMag_(),
    fragBoundMagPrev_(),
    nFragAlphafChanged_(0),
    nFragBoundFaces_(0),
    fragAlphafCopy_(),
    sweptPieceVolumeWord_
    (
        modelDict().getOrDefault<word>
        (
            "sweptPieceVolume",
            modelDict().getOrDefault<bool>("sweptExactCaps", false)
          ? word("closureThenTets") : word("closure")
        )
    ),
    sweptPieceMode_(0),
    sweptPieceValidate_(modelDict().getOrDefault<bool>("sweptPieceValidate", false)),
    nSweptClosurePieces_(0),
    sweptBoundedFaceFraction_(modelDict().getOrDefault<word>("sweptBoundedFaceFraction", "none")),
    nBoundedFaces_(0),
    nBoundedZero_(0),
    nBoundedPartial_(0),
    maxOneMinusS_(0),
    nRecomputed_(0),
    nRecomputeNoPlane_(0),
    nRecomputeNotBracketed_(0),
    recomputeEvals_(0),
    recomputeEvalsMax_(0),
    sweptFaceTris_(),
    sweptFacePlaneN_(),
    sweptFacePlaneX0_(),
    sweptFaceHasPlane_(),
    upLazy_(),
    upDone_(),
    anyPatchPoint_(),
    upWeights_(16),
    nUpLazy_(0),
    cellTets_(),
    nClsFaces_(label(0)),
    nClsTouched_(label(0)),
    nClsPieces_(label(0)),
    nClsVolClips_(label(0)),
    nClsAreaClips_(label(0)),
    nClsTris_(label(0)),
    tClsVol_(scalar(0)),
    tClsArea_(scalar(0)),
    tClsTrace_(scalar(0)),
    tClsPre_(scalar(0)),
    sweptCapSet_(modelDict().getOrDefault<word>("sweptCapSet", "final")),
    sweptConsistencyLog_(modelDict().getOrDefault<bool>("sweptConsistencyLog", false)),
    sweptConsistent_(modelDict().getOrDefault<bool>("sweptConsistent", true)),
    sweptConsistentTol_(modelDict().getOrDefault<scalar>("sweptConsistentTol", 1e-8)),
    nConsProbed_(0),
    nConsReach_(0),
    nConsReachInterface_(0),
    nConsReachWisp_(0),
    nConsReachOpposite_(0),
    nConsProbeFallback_(0),
    consAlgMisSum_(0),
    consAlgMisAbs_(0),
    consAlgMisMax_(0),
    nConsWispPieces_(0),
    consWispMisSum_(0),
    consWispMisAbs_(0),
    consWispMisMax_(0),
    sweptFaceConsDelta_(0),
    nConsCandAdded_(0),
    nConsNoSurface_(0),
    nConsFallback_(0),
    nConsBoundaryKept_(0),
    sweptNoDonorSurface_(false),
    sweptWispSurface_(modelDict().getOrDefault<word>("sweptWispSurface", "matched")),
    sweptWispTol_(modelDict().getOrDefault<scalar>("sweptWispTol", 0)),
    wispSurfaces_(),
    hasWispSurface_(),
    nSweptWispClipped_(0),
    nSweptWispNoNeighbour_(0),
    nSweptWispAnchorFail_(0),
    nSweptWispMatched_(0),
    sweptWispMatchTol_(modelDict().getOrDefault<scalar>("sweptWispMatchTol", 1e-9)),
    consResidueStep_(),
    consSummary_(),
    consResiduePtr_(),
    planeSource_(modelDict().getOrDefault<word>("planeSource", "anchor")),
    sweptFaceDump_(modelDict().getOrDefault<bool>("sweptFaceDump", false)),
    dumpFaces_(),
    dumpPieces_(),
    sweptScaleToPhi_(modelDict().getOrDefault<bool>("sweptScaleToPhi", false)),
    sweptVolumeCorrection_
    (
        // pyramid when no key (the reference setting, 2026-11-04); auto:
        // the mode of sweptScaleToPhi, the builds before the key
        [this]() -> word
        {
            const word w =
                modelDict().getOrDefault<word>("sweptVolumeCorrection", "pyramid");
            if (w != "auto") return w;
            return
                modelDict().getOrDefault<bool>("sweptScaleToPhi", false)
              ? word("scale") : word("none");
        }()
    ),
    sweptVolumeKeyGiven_
    (
        modelDict().getOrDefault<word>("sweptVolumeCorrection", "pyramid") != "auto"
    ),
    sweptPyrDeltaMax_(0),
    sweptVolDefectMax_(0),
    sweptDispAxisMax_(0),
    nSweptTwisted_(0),
    nSweptPyrSkipped_(0),
    sweptRemote_(),
    sweptRemoteRing_(),
    nSweptRemoteCells_(0),
    nSweptRemoteFaces_(0),
    nSweptEulerTrace_(0),
    sweptFaceRemote_(false),
    sweptFaceFallback_(false),
    sweptExchangeTime_(0),
    sweptCurvedTime_(0),
    sweptExchangeParts_(Zero),
    sweptStencilDone_(),
    sweptMap_(),
    sweptPrevNeed_(),
    sweptPrevCompact_(),
    sweptFaceLog_(modelDict().getOrDefault<bool>("sweptFaceLog", false)),
    sweptDebug_(modelDict().found("sweptDebugPoint")),
    sweptDebugPoint_(modelDict().getOrDefault<point>("sweptDebugPoint", point(GREAT, GREAT, GREAT))),
    sweptFaceLogFile_(),
    debugR0Max_(0),
    debugCellR0_(-1),
    debugSurface_(),
    debugAlpha_(0),
    nAnchorFail_(0),
    nAnchorLeaves_(0),
    nAnchorEvals_(0),
    nAnchorEdge_(0),
    nAnchorStartCells_(0),
    anchorStartResMax_(0),
    anchorStartResSum_(0),
    nAnchorStartBelowTol_(0),
    nAnchorSkipped_(0),
    anchorKbVsCutMax_(-1),
    nCompare_(0),
    nCompareDiff_(0),
    nCompareDiffEdge_(0),
    nCompareDiffResolved_(0),
    anchorShiftMax_(0),
    anchorShiftSum_(0),
    defectLevels_(modelDict().getOrDefault<label>("defectLevels", 3)),
    writeCurvedFields_(modelDict().getOrDefault<bool>("writeCurvedFields", false)),
    cellSurfaces_(0),
    hasSurface_(),
    volumeRulePtr_(nullptr),
    nRuleFallback_(0),
    nCurvedFaces_(0),
    nCurvedSecondRing_(0),
    nCurvedBoundary_(0),
    nCurvedSkipped_(0),
    nCurvedFallback_(0),
    nCurvedLeaves_(0),
    nFitCells_(0),
    nFitConverged_(0),
    nFitStalled_(0),
    nFitSingular_(0),
    nFitNoBracket_(0),
    nFitNoFrame_(0),
    fitResMax_(0),
    fitResSum_(0),
    defectMax_(0),
    defectSum_(0),
    nDefectCells_(0),
    fitResidualPtr_(nullptr),
    paraboloidDefectPtr_(nullptr)
{

    // decisions 1 and 2 (2026-11-07): the piece-volume algorithm and the
    // bounded face fractions
    if (sweptPieceVolumeWord_ == "closure") sweptPieceMode_ = 0;
    else if (sweptPieceVolumeWord_ == "oneDonorFallback")
    {
        sweptPieceMode_ = 0;
        Info<< "davofSwept: sweptPieceVolume oneDonorFallback is a DEPRECATED"
            << " alias of closure" << endl;
    }
    else if (sweptPieceVolumeWord_ == "closureThenTets") sweptPieceMode_ = 1;
    else if (sweptPieceVolumeWord_ == "signedTets") sweptPieceMode_ = 2;
    else
    {
        FatalErrorInFunction
            << "sweptPieceVolume " << sweptPieceVolumeWord_
            << ": closure | closureThenTets | signedTets" << exit(FatalError);
    }
    if (modelDict().found("sweptExactCaps"))
    {
        Info<< "davofSwept: sweptExactCaps is DEPRECATED: use sweptPieceVolume"
            << " closureThenTets" << endl;
    }
    if
    (
        sweptBoundedFaceFraction_ != "none" && sweptBoundedFaceFraction_ != "scale"
     && sweptBoundedFaceFraction_ != "recompute"
    )
    {
        FatalErrorInFunction
            << "sweptBoundedFaceFraction " << sweptBoundedFaceFraction_
            << ": none | scale | recompute" << exit(FatalError);
    }
    Info<< "davofSwept: sweptPieceVolume " << sweptPieceVolumeWord_
        << ", sweptPieceValidate " << sweptPieceValidate_
        << ", sweptBoundedFaceFraction " << sweptBoundedFaceFraction_ << endl;

    initAreaFractions();

    // The scheme is the dual advection of the DAVOF state: it needs the davof
    // reconstruction (the face fractions alphaf.<group> it transports)
    if (!(storeAreaFractions_ && alphafPtr_ && isA<reconstruction::davof>(*surf_)))
    {
        FatalIOErrorInFunction(modelDict())
            << "advectionScheme davofSwept needs reconstructionScheme davof"
            << " with stored face fractions (storeAreaFractions true)"
            << exit(FatalIOError);
    }
    davofPtr_ = &refCast<reconstruction::davof>(*surf_);

    // The fork's key names (TwoPhaseFlow branch davof-swept) are accepted as
    // aliases for one release: clippingSurface for surface, fluxIntegration
    // for flux (docs/davof/PORT-DESIGN.md, 2026-11-06)
    if (modelDict().found("clippingSurface") && !modelDict().found("surface"))
    {
        Info<< "davofSwept: the key clippingSurface is deprecated, use surface"
            << endl;
    }
    if (modelDict().found("fluxIntegration") && !modelDict().found("flux"))
    {
        Info<< "davofSwept: the key fluxIntegration is deprecated, use flux"
            << endl;
    }

    {
        // The dual advection with ONE surface per cell for the volume flux
        // and the face update (davofSweptAdvectionSurfaces.C)
        {
            // plane (2026-10-28): the paraboloid path with a zero curvature
            // tensor, i.e. the plane in the frame of the used normal
            // positioned by the volume anchor (fitDof none or one, fitAnchor
            // volume), the surface of the swept-face flux
            if
            (
                clippingSurface_ != "exactSphere" && clippingSurface_ != "paraboloid"
             && clippingSurface_ != "plane"
            )
            {
                FatalIOErrorInFunction(modelDict())
                    << "surface '" << clippingSurface_
                    << "' unknown. Valid: plane, paraboloid, exactSphere."
                    << exit(FatalIOError);
            }
            if (clippingSurface_ == "paraboloid")
            {
                davofPtr_->requestCurvatureTensor(true);
            }
            Info<< "davofSwept: surface " << clippingSurface_
                << " (the volume flux and the face update from one surface"
                << " per cell; curvedTolRel " << curvedTolRel_
                << ", curvedMaxDepth " << curvedMaxDepth_;
            if (planeSource_ != "anchor" && planeSource_ != "exact")
            {
                FatalIOErrorInFunction(modelDict())
                    << "planeSource '" << planeSource_
                    << "' unknown. Valid: anchor, exact." << exit(FatalIOError);
            }
            if (planeSource_ == "exact")
            {
                Info<< ", planeSource exact (the tangent plane of the exact"
                    << " sphere per interface cell, no anchor)";
            }
            if (clippingSurface_ == "plane" && fitDof_ == "six")
            {
                FatalIOErrorInFunction(modelDict())
                    << "surface plane needs fitDof none or one (the"
                    << " six-coefficient fit would set a curvature tensor)"
                    << exit(FatalIOError);
            }
            if (clippingSurface_ != "exactSphere")
            {
                if (fitAnchor_ != "faces" && fitAnchor_ != "volume")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitAnchor '" << fitAnchor_
                        << "' unknown. Valid: faces, volume." << exit(FatalIOError);
                }
                if (fitDof_ != "six" && fitDof_ != "one" && fitDof_ != "none")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitDof '" << fitDof_
                        << "' unknown. Valid: six, one, none." << exit(FatalIOError);
                }
                if (fitDof_ == "none" && fitAnchor_ != "volume")
                {
                    FatalIOErrorInFunction(modelDict())
                        << "fitDof none needs fitAnchor volume (the offset"
                        << " comes from nowhere else)" << exit(FatalIOError);
                }
                if
                (
                    anchorMethod_ != "ccs" && anchorMethod_ != "newton"
                 && anchorMethod_ != "illinois" && anchorMethod_ != "compare"
                )
                {
                    FatalIOErrorInFunction(modelDict())
                        << "anchorMethod '" << anchorMethod_
                        << "' unknown. Valid: ccs, newton, illinois, compare."
                        << exit(FatalIOError);
                }
                Info<< ", fitPriorWeight " << fitPriorWeight_ << ", fitTol "
                    << fitTol_ << ", fitMaxIter " << fitMaxIter_
                    << ", fitAnchor " << fitAnchor_ << ", fitDof " << fitDof_
                    << ", anchorTolRel " << anchorTolRel_ << ", anchorMethod "
                    << anchorMethod_ << ", defectEvery "
                    << defectEvery_ << ", anchorDiagnostics "
                    << (anchorDiagnostics_ ? "on" : "off")
                    << ", anchorSkipConverged "
                    << (anchorSkipConverged_ ? "true" : "false");
            }
            // The flux of a face: prism (the liquid volume of the prism between
            // the face's pre-image and the face under the upwind surface, the
            // paraboloid's flux) or sweptFaceCells (the piecewise swept-face
            // flux, planes only). The fork's one-donor sweptFace and its
            // time-quadrature fluxes (trapezoidal, simpson, simpsonSplit,
            // simpsonSplitS, gaussSplit; branch davof-fluxint) are not ported
            // (PORT-DESIGN.md section 3); the one-donor clip remains as the
            // rank-deficient-cap fallback of sweptFaceCells (sweptFaceFlux)
            if (fluxIntegration_ != "prism" && fluxIntegration_ != "sweptFaceCells")
            {
                FatalIOErrorInFunction(modelDict())
                    << "flux '" << fluxIntegration_
                    << "' unknown. Valid: prism, sweptFaceCells."
                    << exit(FatalIOError);
            }
            const bool swept = (fluxIntegration_ == "sweptFaceCells");
            if (swept && clippingSurface_ == "exactSphere")
            {
                FatalIOErrorInFunction(modelDict())
                    << "flux " << fluxIntegration_ << " needs a plane surface"
                    << " (surface plane, or paraboloid with a zero"
                    << " tensor); exactSphere is a sphere" << exit(FatalIOError);
            }
            Info<< ", flux " << fluxIntegration_;
            if (swept)
            {
                // the volume correction (2026-11-02): the mode replaces the
                // flag; sweptScaleToPhi stays the flag of the mode scale
                if
                (
                    sweptVolumeCorrection_ != "none" && sweptVolumeCorrection_ != "scale"
                 && sweptVolumeCorrection_ != "pyramid"
                )
                {
                    FatalIOErrorInFunction(modelDict())
                        << "sweptVolumeCorrection '" << sweptVolumeCorrection_
                        << "' unknown. Valid: auto (sweptScaleToPhi decides), none,"
                        << " scale, pyramid." << exit(FatalIOError);
                }
                sweptScaleToPhi_ = (sweptVolumeCorrection_ == "scale");
                Info<< " (sweptPlaneTol " << sweptPlaneTol_ << ", sweptScaleToPhi "
                    << (sweptScaleToPhi_ ? "true" : "false");
                if (sweptVolumeKeyGiven_)
                {
                    Info<< ", sweptVolumeCorrection " << sweptVolumeCorrection_
                        << " (V_f = phi_f dt: phi of the half step, setFlow at"
                        << " t + dt/2, paired with the midpoint trace)";
                }
                if (fluxIntegration_ == "sweptFaceCells")
                {
                    if (sweptCellsFaceUpdate_ != "pieces" && sweptCellsFaceUpdate_ != "donor")
                    {
                        FatalIOErrorInFunction(modelDict())
                            << "sweptCellsFaceUpdate '" << sweptCellsFaceUpdate_
                            << "' unknown. Valid: pieces, donor." << exit(FatalIOError);
                    }
                    if (sweptCoplanar_ != "donorSide" && sweptCoplanar_ != "both")
                    {
                        FatalIOErrorInFunction(modelDict())
                            << "sweptCoplanar '" << sweptCoplanar_
                            << "' unknown. Valid: donorSide, both." << exit(FatalIOError);
                    }
                    if (sweptCapSet_ != "final" && sweptCapSet_ != "crossing")
                    {
                        FatalIOErrorInFunction(modelDict())
                            << "sweptCapSet '" << sweptCapSet_
                            << "' unknown. Valid: final, crossing." << exit(FatalIOError);
                    }
                    if
                    (
                        sweptWispSurface_ != "none" && sweptWispSurface_ != "neighbour"
                     && sweptWispSurface_ != "anchored" && sweptWispSurface_ != "matched"
                    )
                    {
                        FatalIOErrorInFunction(modelDict())
                            << "sweptWispSurface '" << sweptWispSurface_
                            << "' unknown. Valid: none, neighbour, anchored, matched."
                            << exit(FatalIOError);
                    }
                    if (sweptPreTest_ != "none" && sweptPreTest_ != "strict" && sweptPreTest_ != "touching")
                    {
                        FatalIOErrorInFunction(modelDict())
                            << "sweptPreTest '" << sweptPreTest_
                            << "' unknown. Valid: none, strict, touching."
                            << exit(FatalIOError);
                    }
                    Info<< ", sweptCellsFaceUpdate " << sweptCellsFaceUpdate_
                        << ", sweptCoplanar " << sweptCoplanar_
                        << ", sweptCapSet " << sweptCapSet_
                        << ", sweptPreTest " << sweptPreTest_
                        << ", sweptStateShortcut " << sweptStateShortcut_
                        << ", sweptTraceCache " << sweptTraceCache_;
                    if (sweptPreTestCheck_) Info<< ", sweptPreTestCheck true";
                    if (sweptClassStats_) Info<< ", sweptClassStats true";
                    Info<< ", sweptDonorFastPath " << sweptDonorFastPath_
                        << ", sweptFaceUpdateFromVolume " << sweptFaceUpdateFromVolume_
                        << ", sweptPreTestCoplanarSide " << sweptPreTestCoplanarSide_
                        << ", sweptNoFluxLayer " << sweptNoFluxLayer_
                        << ", sweptAxisPlanes " << sweptAxisPlanes_
                        << ", sweptPointVelocityLazy " << sweptPointVelocityLazy_
                        << ", sweptTetCache " << sweptTetCache_
                        << ", sweptSliverClosedForm " << sweptSliverClosedForm_;
                    // the wisp rule (2026-11-02), printed only when not none
                    if (sweptWispSurface_ != "none")
                    {
                        Info<< ", sweptWispSurface " << sweptWispSurface_;
                    }
                    Info<< ", sweptWispTol " << sweptWispTol_
                        << ", sweptWispMatchTol " << sweptWispMatchTol_;
                    // the consistency keys (2026-11-02), printed only when on
                    if (sweptConsistent_)
                    {
                        Info<< ", sweptConsistent true (sweptConsistentTol "
                            << sweptConsistentTol_ << ")";
                    }
                    if (sweptConsistencyLog_)
                    {
                        if (Pstream::parRun())
                        {
                            FatalIOErrorInFunction(modelDict())
                                << "sweptConsistencyLog true is a serial"
                                << " diagnostic (it probes the algebraic faces"
                                << " in the local vertex rings)"
                                << exit(FatalIOError);
                        }
                        Info<< ", sweptConsistencyLog true";
                    }
                }
                Info<< ")";
            }
            Info<< ", defectLevels " << defectLevels_ << ")" << endl;

            // The rule for the volumes and face fractions under a surface
            volumeRulePtr_ = davofVolumeRule::New(mesh_, modelDict());
        }

        // The face state on warped faces (2026-10-25): the wetted area
        // vector next to the scalar fraction, produced by the curved face
        // update only
        if (faceState_ != "scalar" && faceState_ != "vector")
        {
            FatalIOErrorInFunction(modelDict())
                << "faceState '" << faceState_
                << "' unknown. Valid: scalar, vector." << exit(FatalIOError);
        }
        if (tensorCentroid_ != "cut" && tensorCentroid_ != "anchored")
        {
            FatalIOErrorInFunction(modelDict())
                << "tensorCentroid '" << tensorCentroid_
                << "' unknown. Valid: cut, anchored." << exit(FatalIOError);
        }
        if (faceState_ == "vector")
        {
            wetVecPtr_ = davofPtr_->wettedAreaPtr();
            if (!wetVecPtr_)
            {
                FatalIOErrorInFunction(modelDict())
                    << "faceState vector: the davof scheme holds no wetted"
                    << " area vector field" << exit(FatalIOError);
            }
            Info<< "davofSwept: faceState vector (the wetted area vector "
                << wetVecPtr_->name() << " transported with alphaf)" << endl;
        }
        if (tensorCentroid_ == "anchored" && clippingSurface_ != "paraboloid")
        {
            FatalIOErrorInFunction(modelDict())
                << "tensorCentroid anchored needs surface paraboloid"
                << " (the ring tensor is fitted there only)" << exit(FatalIOError);
        }
        if (tensorCentroid_ == "anchored")
        {
            Info<< "davofSwept: tensorCentroid anchored (the ring tensor's"
                << " centroids re-anchored by the volume rule)" << endl;
        }
    }

    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        setProcessorPatches();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::davofSweptAdvection::setProcessorPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    surfaceCellFacesOnProcPatches_.clear();
    surfaceCellFacesOnProcPatches_.resize(patches.size());

    // Append all processor patch labels to the list
    procPatchLabels_.clear();
    forAll(patches, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchi])
            && patches[patchi].size() > 0
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}

void Foam::davofSweptAdvection::extendMarkedCells
(
    bitSet& markedCell
) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace: the marked faces are visited
    // (the bitSet iterates its set bits, skipping the zero words)
    // instead of testing every face of the mesh; the resulting SET of
    // cells is the same, and its consumers (boundFlux, the second-layer
    // count) iterate the bitSet in cell order whatever the order the
    // bits were set in
    const label nInt = mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    for (const label facei : markedFace)
    {
        markedCell.set(own[facei]);
        if (facei < nInt)
        {
            markedCell.set(nei[facei]);
        }
    }
}

void Foam::davofSweptAdvection::initAreaFractions()
{
    const word name(IOobject::groupName("alphaf", alpha1_.group()));

    if (mesh_.foundObject<surfaceScalarField>(name))
    {
        alphafPtr_ = &mesh_.lookupObjectRef<surfaceScalarField>(name);
    }
    else if (storeAreaFractions_)
    {
        Info<< "davofSwept: creating the face fractions " << name
            << " (upwind) for the DAVOF face update" << endl;
        alphafOwnPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                upwind<scalar>(mesh_, phi_).interpolate(alpha1_)
            )
        );
        alphafPtr_ = alphafOwnPtr_.get();
    }
}


void Foam::davofSweptAdvection::setDefaultAreaFraction
(
    const label facei,
    const scalarField& alpha
)
{
    const labelList& own = mesh_.faceOwner();

    if (mesh_.isInternalFace(facei))
    {
        const scalar phif = phi_.primitiveField()[facei];
        const scalar af =
            (phif >= 0) ? alpha[own[facei]] : alpha[mesh_.faceNeighbour()[facei]];
        alphafPtr_->primitiveFieldRef()[facei] = af;
        wetFace_.set(facei, (af != 0) ? 1u : 0u);
        if (wetVecPtr_)
        {
            // the planar default of the vector state: alpha_up S_f (the
            // fan triangles' area vectors sum to S_f)
            wetVecPtr_->primitiveFieldRef()[facei] = af*mesh_.faceAreas()[facei];
        }
        return;
    }

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];
    fvsPatchScalarField& pf = alphafPtr_->boundaryFieldRef()[patchi];
    if (pf.size())
    {
        const scalar af = alpha[own[facei]];
        pf[facei - pbm[patchi].start()] = af;
        wetFace_.set(facei, (af != 0) ? 1u : 0u);
        if (wetVecPtr_)
        {
            wetVecPtr_->boundaryFieldRef()[patchi][facei - pbm[patchi].start()] =
                af*mesh_.faceAreas()[facei];
        }
    }
}


void Foam::davofSweptAdvection::defaultAreaFractions()
{
    addProfilingInFunction(geometricVoF);
    const scalarField& alphaIn = alpha1_.primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const label nInt = mesh_.nInternalFaces();
    const label nFaces = mesh_.nFaces();

    // The target of every face no plane set in this step is alpha of its
    // upwind cell (boundary faces: the owner), with alpha AFTER the
    // bounding. The invariant kept between the steps: such a face holds
    // alphaDefault_ of that cell, alphaDefault_ being a copy of alpha at
    // the last evaluation, as long as nothing has written the face since
    // (a plane, the sync or a fallback writes a face in the step that
    // marks or lists it). The stored value can then differ from the
    // target only on (a) a face whose upwind direction changed, (b) a
    // face set in the previous step, (c) a fallback of this step (it
    // holds the entry alpha of its upwind cell), (d) a face of a cell
    // whose alpha differs bitwise from alphaDefault_. Every other unset
    // face already holds its target; evaluating a face always writes the
    // target, so visiting a face twice is harmless.

    if (!defaultsValid_)
    {
        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
        phiPos_.resize(nInt);
        phiPos_.reset();
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (phiIn[facei] >= 0)
            {
                phiPos_.set(facei);
            }
        }
        phiEventNo_ = phi_.eventNo();
        alphaDefault_ = alphaIn;

        // The wet faces from the final values (the set faces included)
        wetFace_.resize(nFaces);
        wetFace_.reset();
        const scalarField& afIn = alphafPtr_->primitiveField();
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (afIn[facei] != 0) wetFace_.set(facei);
        }
        const surfaceScalarField::Boundary& bf = alphafPtr_->boundaryField();
        forAll(bf, patchi)
        {
            const label start = mesh_.boundaryMesh()[patchi].start();
            const fvsPatchScalarField& pf = bf[patchi];
            forAll(pf, i)
            {
                if (pf[i] != 0) wetFace_.set(start + i);
            }
        }
        defaultsValid_ = true;
        return;
    }

    // (a) phi written since the last evaluation (its event number moved:
    // every non-const access of a GeometricField records one): re-scan
    // the upwind direction, the faces whose direction changed
    if (phi_.eventNo() != phiEventNo_)
    {
        for (label facei = 0; facei < nInt; ++facei)
        {
            const bool pos = (phiIn[facei] >= 0);
            if (pos == phiPos_.test(facei)) continue;
            if (pos) phiPos_.set(facei); else phiPos_.unset(facei);
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
        phiEventNo_ = phi_.eventNo();
    }

    // (b) the faces set in the previous step and (c) the fallbacks
    for (const label facei : alphafPrevSetFaces_)
    {
        if (!alphafSetFace_.test(facei))
        {
            setDefaultAreaFraction(facei, alphaIn);
        }
    }
    for (const label facei : alphafFallbackFaces_)
    {
        if (!alphafSetFace_.test(facei))
        {
            setDefaultAreaFraction(facei, alphaIn);
        }
    }

    // (d) the faces of the cells whose alpha changed, compared bitwise
    // (-0 and +0 differ in the written field; NaN never occurs here)
    const cellList& cells = mesh_.cells();
    forAll(alphaIn, celli)
    {
        if (std::memcmp(&alphaIn[celli], &alphaDefault_[celli], sizeof(scalar)) == 0)
        {
            continue;
        }
        alphaDefault_[celli] = alphaIn[celli];
        for (const label facei : cells[celli])
        {
            if (!alphafSetFace_.test(facei))
            {
                setDefaultAreaFraction(facei, alphaIn);
            }
        }
    }
}


void Foam::davofSweptAdvection::publishWetFaces()
{
    wetFaces_ = wetFace_.toc();
    if (davofPtr_)
    {
        davofPtr_->setWetFaces(wetFaces_);
    }
}


void Foam::davofSweptAdvection::setDownwindFaces
(
    const label celli,
    DynamicLabelList& downwindFaces
) const
{
    DebugInFunction << endl;

    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[celli];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label facei = c[fi];
        const scalar phi = faceValue(phi_, facei);

        if (own[facei] == celli)
        {
            if (phi >= 0)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < 0)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


Foam::scalar Foam::davofSweptAdvection::netFlux
(
    const surfaceScalarField& dVf,
    const label celli
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[celli];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == celli)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


Foam::DynamicList<Foam::label>  Foam::davofSweptAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces,
    surfaceScalarField* alphaf,
    surfaceVectorField* wetVec
)
{
    DynamicLabelList syncedFaces(0);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;

            if (alphaf)
            {
                const UIndirectList<scalar> afPatch
                (
                    alphaf->boundaryField()[patchi],
                    surfCellFacesOnProcPatch
                );
                toNbr << afPatch;
            }
            if (wetVec)
            {
                const UIndirectList<vector> wvPatch
                (
                    wetVec->boundaryField()[patchi],
                    surfCellFacesOnProcPatch
                );
                toNbr << wvPatch;
            }
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;
            List<scalar> nbrAf;
            if (alphaf)
            {
                fromNeighb >> nbrAf;
            }
            List<vector> nbrWv;
            if (wetVec)
            {
                fromNeighb >> nbrWv;
            }
            if (returnSyncedFaces)
            {
                List<label> syncedFaceI(faceIDs);
                for (label& faceI : syncedFaceI)
                {
                    faceI += procPatch.start();
                }
                syncedFaces.append(syncedFaceI);
            }

            if (debug)
            {
                Pout<< "Received at time = " << mesh_.time().value()
                    << ": surfCellFacesOnProcPatch = " << faceIDs << nl
                    << "Received at time = " << mesh_.time().value()
                    << ": dVfPatch = " << nbrdVfs << endl;
            }

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
                if (alphaf)
                {
                    // The upwind side's face fraction, no sign change
                    alphaf->boundaryFieldRef()[patchi][facei] = nbrAf[i];
                    if (wetVec)
                    {
                        // the wetted area vector is oriented with each
                        // side's own outward patch face: the sign flips,
                        // as the flux's does
                        wetVec->boundaryFieldRef()[patchi][facei] = -nbrWv[i];
                    }
                    const label meshFacei = procPatch.start() + facei;
                    wetFace_.set(meshFacei, (nbrAf[i] != 0) ? 1u : 0u);
                    if (!alphafSetFace_.test(meshFacei))
                    {
                        alphafSetFaces_.append(meshFacei);
                    }
                    alphafSetFace_.set(meshFacei);
                }
                if (debug && mag(localFlux[facei] + nbrdVfs[i]) > ROOTVSMALL)
                {
                    Pout<< "localFlux[facei] = " << localFlux[facei]
                        << " and nbrdVfs[i] = " << nbrdVfs[i]
                        << " for facei = " << facei << endl;
                }
            }
        }

        if (debug)
        {
            // Write out results for checking
            forAll(procPatchLabels_, patchLabeli)
            {
                const label patchi = procPatchLabels_[patchLabeli];
                const scalarField& localFlux = dVf.boundaryField()[patchi];
                Pout<< "time = " << mesh_.time().value() << ": localFlux = "
                    << localFlux << endl;
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }

    return syncedFaces;
}


void Foam::davofSweptAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[facei - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}




void Foam::davofSweptAdvection::applyBruteForceBounding()
{
    addProfilingInFunction(geometricVoF);
    bool alpha1Changed = false;

    scalar snapAlphaTol = modelDict().lookupOrDefault<scalar>("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    if (modelDict().lookupOrDefault("clip", true))
    {
        alpha1_ = min(scalar(1), max(scalar(0), alpha1_));
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}

// ************************************************************************* //
