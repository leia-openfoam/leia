/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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

#include "davof.H"
#include "emptyPolyPatch.H"
#include "davofInit.H"
#include "davofCcsPositioning.H"
#include "zoneDistribute.H"
#include "syncTools.H"
#include "cutFacePLIC.H"
#include "coupledFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(davof, 0);
    addToRunTimeSelectionTable(reconstructionSchemes, davof, components);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceScalarField>
Foam::reconstruction::davof::readAlphafIfAbsent(const volScalarField& alpha1)
{
    const fvMesh& mesh = alpha1.mesh();
    const word name(IOobject::groupName("alphaf", alpha1.group()));

    if (mesh.foundObject<surfaceScalarField>(name))
    {
        return nullptr;
    }

    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!io.typeHeaderOk<surfaceScalarField>(true))
    {
        FatalErrorInFunction
            << "The face-fraction field " << name
            << " is neither registered nor in " << mesh.time().timeName()
            << "/ of the case." << nl
            << "The davof reconstruction needs the face fractions of the"
            << " DAVOF state: run initDavofFields (system/initAlphaFieldDict"
            << " with davofMethod) to write " << name << "."
            << exit(FatalError);
    }

    Info<< "davof: reading " << name << " from " << mesh.time().timeName()
        << endl;
    return autoPtr<surfaceScalarField>::New(io, mesh);
}


Foam::autoPtr<Foam::surfaceVectorField>
Foam::reconstruction::davof::readWetVecIfAbsent
(
    const volScalarField& alpha1,
    const word& faceState
)
{
    if (faceState != "vector")
    {
        return nullptr;
    }
    const fvMesh& mesh = alpha1.mesh();
    const word name(IOobject::groupName("wettedAreaDavof", alpha1.group()));

    if (mesh.foundObject<surfaceVectorField>(name))
    {
        return nullptr;
    }

    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!io.typeHeaderOk<surfaceVectorField>(true))
    {
        FatalErrorInFunction
            << "faceState vector: the wetted area vector field " << name
            << " is neither registered nor in " << mesh.time().timeName()
            << "/ of the case." << nl
            << "Run initDavofFields with faceState vector in"
            << " system/initAlphaFieldDict to write " << name << "."
            << exit(FatalError);
    }

    Info<< "davof: reading " << name << " from " << mesh.time().timeName()
        << endl;
    return autoPtr<surfaceVectorField>::New(io, mesh);
}


bool Foam::reconstruction::davof::positionPlane
(
    const label celli,
    const vector& n
)
{
    // The RETURNED status: vofCutCell returns -1 (alpha below the
    // tolerance) or 1 (above 1 - tolerance) before it cuts anything, and
    // cellStatus() then still holds the previous cell's cut (surfaceArea
    // and surfaceCentre with it). An interface cell is selected here by its
    // face fractions, so its alpha can be at 0 or 1; reading cellStatus()
    // gave such a cell the plane of the previously positioned cell (a
    // phantom surface cell whose plane sits elsewhere, which isoAdvection
    // then fluxed; found on the sphere 2026-10-13).
    if (positioning_ == "ccs")
    {
        // Consecutive Cubic Spline positioning on the same cut, from the
        // explicit position pPlane_[celli] along n (2026-11-04)
        const scalar nMag = mag(n);
        if (nMag == 0)
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
            return false;
        }
        const scalar ac = alpha1_[celli];
        Foam::davof::ccsRecord rec;
        const label status = Foam::davof::ccsPositionPlane
        (
            mesh_, cutCellPos_, celli, ac, n/nMag, pPlane_[celli],
            positioningTol_, positioningMaxIter_, rec
        );
        const bool atEnd = (ac < isoFaceTol_) || (1 - ac < isoFaceTol_);
        ++nPosCells_;
        nPosEvals_ += rec.nEvals;
        maxPosEvals_ = max(maxPosEvals_, rec.nEvals);
        ++posEvalHist_[min(rec.nEvals, posEvalHist_.size() - 1)];
        if (rec.exitCode == 2) ++nPosFloor_;
        if (rec.exitCode == 4) ++nPosMaxIter_;
        if (atEnd) ++nPosEnd_;
        const scalar res = mag(rec.rFinal);
        positioningResidual_[celli] = rec.rFinal;   // signed, of V
        maxPosRes_ = max(maxPosRes_, res);
        sumPosRes_ += res;

        if (status == 0)
        {
            normal_[celli] = cutCellPos_.faceArea();
            centre_[celli] = cutCellPos_.faceCentre();
            if (mag(normal_[celli]) != 0)
            {
                return true;
            }
        }
        if (atEnd) ++nPosEndNoPlane_;
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        return false;
    }

    // No normal (an alphaOnly interface cell whose ring sum vanishes,
    // 2026-11-09): no plane. Never taken on the default path, where every
    // selected normal has a positive magnitude.
    if (mag(n) == 0)
    {
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        return false;
    }

    const label status = sIterPLIC_.vofCutCell
    (
        celli,
        alpha1_[celli],
        isoFaceTol_,
        100,
        n
    );

    if (status == 0)
    {
        normal_[celli] = sIterPLIC_.surfaceArea();
        centre_[celli] = sIterPLIC_.surfaceCentre();
        if (mag(normal_[celli]) != 0)
        {
            return true;
        }
    }

    normal_[celli] = Zero;
    centre_[celli] = Zero;
    return false;
}


void Foam::reconstruction::davof::updateHasValue()
{
    if (hasValue_.size() == mesh_.nFaces() && !mesh_.topoChanging())
    {
        return;
    }

    // Faces that carry a value: the internal faces and those of patches
    // whose fvsPatchField has the patch size (empty patches have none).
    hasValue_.resize(mesh_.nFaces());
    hasValue_ = false;
    emptyFace_.resize(mesh_.nFaces());
    emptyFace_ = false;
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        hasValue_[facei] = true;
    }
    const surfaceScalarField::Boundary& bf = alphaf_.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (bf[patchi].size() == pp.size())
        {
            forAll(pp, i)
            {
                hasValue_[pp.start() + i] = true;
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                emptyFace_[pp.start() + i] = true;
            }
        }
    }
}


void Foam::reconstruction::davof::formAreaNormal()
{
    addProfilingInFunction(geometricVoF);
    const vectorField& Sf = mesh_.faceAreas();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const label nCells = mesh_.nCells();

    vectorField& m = areaNormal_.primitiveFieldRef();

    // The per-cell sums of a pass over every face in face order, the
    // dry faces (alpha_f = 0) skipped (davofInit::areaNormal: 0 S_f adds
    // +-0, which changes no sum, and the sums start at +0). The wet list
    // holds exactly the faces with alpha_f != 0, ascending (the internal
    // faces, then the patch faces in patch order, as the full pass
    // visits them), so every cell receives the same terms in the same
    // order. A cell no wet face touches keeps the +0 it holds: the cells
    // of the previous sum are zeroed, every other cell was never written
    // since the last full zeroing.
    const bool useList =
    (
        wetFacesPtr_
     && wetFacesNFaces_ == mesh_.nFaces()
     && !mesh_.topoChanging()
     && touchedCell_.size() == nCells
     && m.size() == nCells
    );

    if (useList)
    {
        for (const label celli : touchedCell_)
        {
            m[celli] = Zero;
        }
        touchedCell_.reset();

        const scalarField& afIn = alphaf_.primitiveField();
        const vectorField* wvIn = wetVecPtr_ ? &wetVecPtr_->primitiveField() : nullptr;
        for (const label facei : *wetFacesPtr_)
        {
            const scalar af = (facei < nInt) ? afIn[facei] : faceFraction(facei);
            if (af == 0) continue;
            // faceState vector: the wetted area vector in place of af S_f
            const vector contrib =
                wvIn
              ? ((facei < nInt) ? (*wvIn)[facei] : faceWettedArea(facei))
              : af*Sf[facei];
            m[own[facei]] += contrib;
            touchedCell_.set(own[facei]);
            if (facei < nInt)
            {
                m[nei[facei]] -= contrib;
                touchedCell_.set(nei[facei]);
            }
        }
        return;
    }

    m = Zero;
    touchedCell_.resize(nCells);
    touchedCell_.reset();

    const scalarField& afIn = alphaf_.primitiveField();
    const vectorField* wvIn = wetVecPtr_ ? &wetVecPtr_->primitiveField() : nullptr;
    for (label facei = 0; facei < nInt; ++facei)
    {
        const scalar af = afIn[facei];
        if (af == 0) continue;
        const vector contrib = wvIn ? (*wvIn)[facei] : af*Sf[facei];
        m[own[facei]] += contrib;
        m[nei[facei]] -= contrib;
        touchedCell_.set(own[facei]);
        touchedCell_.set(nei[facei]);
    }
    const surfaceScalarField::Boundary& bf = alphaf_.boundaryField();
    forAll(bf, patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const fvsPatchScalarField& pf = bf[patchi];
        forAll(pf, i)
        {
            const label facei = pp.start() + i;
            if (pf[i] == 0) continue;
            m[own[facei]] +=
                wetVecPtr_
              ? wetVecPtr_->boundaryField()[patchi][i]
              : pf[i]*Sf[facei];
            touchedCell_.set(own[facei]);
        }
    }
}


Foam::scalar Foam::reconstruction::davof::stateTime()
{
    const Time& runTime = mesh_.time();

    if (runTime.timeIndex() != lastTimeIndex_)
    {
        lastTimeIndex_ = runTime.timeIndex();
        nReconstructInIndex_ = 0;
    }
    ++nReconstructInIndex_;

    // The constructor (start time index) and every reconstruction after
    // the advection step's own (the first one of a time index) see the
    // state at the current time; the advection step's sees the state the
    // step starts from.
    if
    (
        runTime.timeIndex() == runTime.startTimeIndex()
     || nReconstructInIndex_ > 1
    )
    {
        return runTime.value();
    }

    return runTime.value() - runTime.deltaTValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::davof::davof
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),
    wispTol_(modelDict().getOrDefault<scalar>("wispTol", 1e-3)),
    surfCellTol_(modelDict().getOrDefault<scalar>("surfCellTol", 1e-8)),
    isoFaceTol_(modelDict().getOrDefault<scalar>("isoFaceTol", 1e-8)),
    interfaceCellRule_(modelDict().getOrDefault<word>("interfaceCellRule", "normalOnly")),
    nRuleOutside_(0),
    nRuleOutsidePlane_(0),
    nRuleResetCells_(0),
    nRuleResetFaces_(0),
    ruleResetMaxDaf_(0),
    faceFractionCorrector_
    (
        modelDict().getOrDefault<word>("faceFractionCorrector", "none")
    ),
    correctorIterations_(modelDict().getOrDefault<label>("correctorIterations", 1)),
    correctorTolerance_(modelDict().getOrDefault<scalar>("correctorTolerance", 1e-12)),
    correctorLog_(modelDict().getOrDefault<bool>("correctorLog", false)),
    correctedTimeIndex_(-1),
    corrIterationsUsed_(0),
    corrMaxChange_(),
    corrL1Change_(),
    corrFaces_(0),
    corrNoPlaneCells_(0),
    corrInterfaceCells_(0),
    corrIdentity_(0),
    positioning_(modelDict().getOrDefault<word>("positioning", "ccs")),
    positioningTol_(modelDict().getOrDefault<scalar>("positioningTol", 1e-14)),
    positioningMaxIter_(modelDict().getOrDefault<label>("positioningMaxIter", 30)),
    grazingCut_(modelDict().getOrDefault<word>("grazingCut", "openfoam")),
    cutCellPos_(alpha1.mesh(), grazingCut_ == "fork"),
    positioningResidual_
    (
        IOobject
        (
            IOobject::groupName("positioningResidual", alpha1.group()),
            alpha1.mesh().time().timeName(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    nPosCells_(0),
    nPosEvals_(0),
    maxPosEvals_(0),
    nPosFloor_(0),
    nPosMaxIter_(0),
    nPosEnd_(0),
    nPosEndNoPlane_(0),
    maxPosRes_(0),
    sumPosRes_(0),
    posEvalHist_(8, Zero),
    lastTimeIndex_(-1),
    nReconstructInIndex_(0),
    sIterPLIC_(mesh_, surfCellTol_),
    cutCell_(mesh_, grazingCut_ == "fork"),
    alphafPtr_(readAlphafIfAbsent(alpha1)),
    alphaf_
    (
        mesh_.lookupObjectRef<surfaceScalarField>
        (
            IOobject::groupName("alphaf", alpha1.group())
        )
    ),
    faceState_(modelDict().getOrDefault<word>("faceState", "scalar")),
    wetVecOwnPtr_(readWetVecIfAbsent(alpha1, faceState_)),
    wetVecPtr_(nullptr),
    hasValue_(),
    emptyFace_(),
    wetFacesPtr_(nullptr),
    wetFacesNFaces_(-1),
    touchedCell_(),
    areaNormal_
    (
        IOobject
        (
            IOobject::groupName("areaNormalDavof", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimArea, Zero)
    ),
    sliverFlag_
    (
        IOobject
        (
            IOobject::groupName("sliverFlag", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    pPlane_
    (
        IOobject
        (
            IOobject::groupName("pPlane", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero)
    ),
    planeVolDefect_
    (
        IOobject
        (
            IOobject::groupName("planeVolDefect", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    planeAreaDefect_
    (
        IOobject
        (
            IOobject::groupName("planeAreaDefect", alpha1.group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    maxPlaneVolDefect_(0),
    maxPlaneAreaDefect_(0),
    normalModel_(nullptr),
    usedNormal_(0),
    tensorRequested_(false),
    tensorAlignTol_(modelDict().getOrDefault<scalar>("tensorAlignTol", 0.2)),
    tensorMinCells_(modelDict().getOrDefault<label>("tensorMinCells", 6)),
    writeCurvatureTensor_
    (
        modelDict().getOrDefault<bool>("writeCurvatureTensor", false)
    ),
    polyMoment_(0),
    curvatureTensorPtr_(nullptr),
    curvatureOffsetPtr_(nullptr),
    mRingPtr_(nullptr),
    mRingFilled_(),
    tensorCells_(),
    tensorCentrePtr_(nullptr),
    tensorCentresValid_(false),
    nTensorCells_(0),
    nTensorFallback_(0)
{
    if (faceState_ == "vector")
    {
        wetVecPtr_ = &mesh_.lookupObjectRef<surfaceVectorField>
        (
            IOobject::groupName("wettedAreaDavof", alpha1.group())
        );
        Info<< "davof: faceState vector (m_c and the face moment from the"
            << " wetted area vector " << wetVecPtr_->name() << ")" << endl;
    }
    else if (faceState_ != "scalar")
    {
        FatalIOErrorInFunction(modelDict())
            << "faceState '" << faceState_
            << "' unknown. Valid: scalar, vector." << exit(FatalIOError);
    }

    if (positioning_ != "vofCutCell" && positioning_ != "ccs")
    {
        FatalIOErrorInFunction(modelDict())
            << "positioning '" << positioning_
            << "' unknown. Valid: vofCutCell, ccs." << exit(FatalIOError);
    }
    if (positioning_ == "ccs")
    {
        Info<< "davof: positioning ccs (Consecutive Cubic Spline on the"
            << " cutCellPLIC cut from the explicit position, positioningTol "
            << positioningTol_ << " of V, positioningMaxIter "
            << positioningMaxIter_ << ")" << endl;
    }
    if (grazingCut_ != "openfoam" && grazingCut_ != "fork")
    {
        FatalIOErrorInFunction(modelDict())
            << "grazingCut '" << grazingCut_
            << "' unknown. Valid: openfoam, fork." << exit(FatalIOError);
    }
    if (grazingCut_ == "fork")
    {
        Info<< "davof: grazingCut fork (the cell cuts classify a cut polygon"
            << " below 10*SMALL as no cut, the TwoPhaseFlow threshold; the"
            << " default openfoam is v2606's cutCellPLIC, ROOTVSMALL)" << endl;
    }

    // The iterative face-fraction corrector (2026-11-09)
    if (faceFractionCorrector_ != "none" && faceFractionCorrector_ != "mean")
    {
        FatalIOErrorInFunction(modelDict())
            << "faceFractionCorrector '" << faceFractionCorrector_
            << "' unknown. Valid: none, mean." << exit(FatalIOError);
    }
    if (faceFractionCorrector_ != "none")
    {
        if (correctorIterations_ < 1)
        {
            FatalIOErrorInFunction(modelDict())
                << "correctorIterations " << correctorIterations_
                << ": at least 1" << exit(FatalIOError);
        }
        if (faceState_ == "vector")
        {
            FatalIOErrorInFunction(modelDict())
                << "faceFractionCorrector " << faceFractionCorrector_
                << " is not implemented for faceState vector"
                << exit(FatalIOError);
        }
        Info<< "davof: faceFractionCorrector " << faceFractionCorrector_
            << " (correctorIterations " << correctorIterations_
            << ", correctorTolerance " << correctorTolerance_
            << ", correctorLog " << correctorLog_ << ")";
        if (interfaceCellRule_ != "alphaOnly")
        {
            Info<< "; interfaceCellRule " << interfaceCellRule_
                << " replaced by alphaOnly (the corrector's interface cells"
                << " are the cells with surfCellTol < alpha < 1 - surfCellTol)";
            interfaceCellRule_ = "alphaOnly";
        }
        Info<< endl;
    }

    normalModel_ = davofNormalModel::New(modelDict(), *this);

    Info<< "davof: normalModel " << normalModel_->type() << ", applyTo "
        << normalModel_->applyTo() << ", sliverTol "
        << normalModel_->sliverTol() << ", maxAngle "
        << normalModel_->maxAngle() << endl;

    reconstruct(true);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::davof::reconstruct(bool forceUpdate)
{
    addProfilingInFunction(geometricVoF);
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    // The face-fraction corrector ran at the end of the previous step and
    // left the planes of the state (alpha, alphaf) it ends with; nothing has
    // changed them since, so the advection step's reconstruction (the first
    // one of this time index) is the corrector's last (2026-11-09)
    if
    (
        faceFractionCorrector_ != "none"
     && correctedTimeIndex_ >= 0
     && mesh_.time().timeIndex() == correctedTimeIndex_ + 1
     && mesh_.time().timeIndex() != lastTimeIndex_
     && !mesh_.topoChanging()
    )
    {
        lastTimeIndex_ = mesh_.time().timeIndex();
        nReconstructInIndex_ = 1;
        if (correctorLog_)
        {
            Info<< "davof: reconstruction skipped, the face-fraction"
                << " corrector's planes of the previous step are current"
                << endl;
        }
        return;
    }

    if (mesh_.topoChanging())
    {
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }

    // The interface cells of the previous reconstruction: the planes and
    // the per-cell diagnostics are set in interface cells only and are
    // zero everywhere else (zero at construction), so they are reset
    // there instead of over the whole mesh.
    const labelList prevLabels(interfaceLabels_);
    interfaceCell_ = false;
    interfaceLabels_.clear();

    const scalar tState = stateTime();

    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const labelUList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();

    updateHasValue();

    // m_c = sum_f alpha_f S_f^out, formed in areaNormal_ itself over the
    // wet faces (the list isoAdvection keeps) or over every face
    formAreaNormal();
    vectorField& m = areaNormal_.primitiveFieldRef();

    for (const label celli : prevLabels)
    {
        normal_[celli] = Zero;
        centre_[celli] = Zero;
        sliverFlag_[celli] = 0;
        pPlane_[celli] = 0;
        planeVolDefect_[celli] = 0;
        planeAreaDefect_[celli] = 0;
        positioningResidual_[celli] = 0;
    }
    maxPlaneVolDefect_ = 0;
    maxPlaneAreaDefect_ = 0;
    nPosCells_ = nPosEvals_ = maxPosEvals_ = nPosFloor_ = nPosMaxIter_ = 0;
    nPosEnd_ = nPosEndNoPlane_ = 0;
    maxPosRes_ = sumPosRes_ = 0;
    posEvalHist_ = Zero;

    // A cell whose face fractions are all equal holds no interface (all
    // wet, all dry): its sum_f S_f vanishes only to round-off. Zero it.
    // A cell no wet face touched has m exactly +0 (the sums start at +0
    // and skip the dry faces): all dry, or an exact cancellation, in
    // either case not an interface cell and m already what the rule
    // below would leave (the former full sum of 0 S_f terms gave +0 as
    // well, and the all-dry rule assigned +0), so it is skipped without
    // reading its faces: the cells the sum touched are visited, in cell
    // order (the bitSet iterates ascending), the others hold +0.
    // decision 3 (2026-11-07): the alpha condition of the interface-cell rule
    if
    (
        interfaceCellRule_ != "normalOnly" && interfaceCellRule_ != "normalAndAlpha"
     && interfaceCellRule_ != "normalAndAlphaReset" && interfaceCellRule_ != "alphaOnly"
    )
    {
        FatalErrorInFunction
            << "interfaceCellRule " << interfaceCellRule_
            << ": normalOnly | normalAndAlpha | normalAndAlphaReset | alphaOnly"
            << exit(FatalError);
    }
    // alphaOnly (2026-11-09, the face-fraction corrector's rule): the
    // interface cells are the cells with tol < alpha < 1 - tol, whatever
    // their face fractions; a cell among them with m_c = 0 gets its normal
    // from the model (the ring sum) or no plane
    const bool alphaOnly = (interfaceCellRule_ == "alphaOnly");
    const scalar tolA = surfCellTol_;
    auto alphaInside = [&](const label celli) -> bool
    {
        const scalar a = alpha1_[celli];
        return (a > tolA && a < 1 - tolA);
    };
    nRuleOutside_ = nRuleOutsidePlane_ = nRuleResetCells_ = nRuleResetFaces_ = 0;
    ruleResetMaxDaf_ = 0;
    if (interfaceCellRule_ == "normalAndAlphaReset")
    {
        // the faces of a dry or full cell with |m_c| above the threshold are
        // set to its state; a face shared with a neighbour changes the
        // neighbour's m_c as well; no mass moves (alpha is untouched)
        labelList resetFlag(mesh_.nFaces(), 0);   // 1: set to 0, 2: set to 1
        bool any = false;
        for (const label celli : touchedCell_)
        {
            if (m[celli] == vector::zero) continue;
            if (alphaInside(celli)) continue;
            if (!(mag(m[celli]) > davofInit::wispThreshold(mesh_, celli, wispTol_))) continue;
            const scalar state = (alpha1_[celli] <= tolA) ? 0 : 1;
            const labelList& cf = cells[celli];
            bool changed = false;
            forAll(cf, i)
            {
                const label f = cf[i];
                if (!hasValue_[f]) continue;
                const scalar af = faceFraction(f);
                if (af == state) continue;
                ruleResetMaxDaf_ = max(ruleResetMaxDaf_, mag(af - state));
                setFaceFraction(f, state);
                resetFlag[f] = (state == 0) ? 1 : 2;
                ++nRuleResetFaces_;
                changed = true;
            }
            if (changed)
            {
                ++nRuleResetCells_;
                any = true;
            }
        }
        if (returnReduce(any, orOp<bool>()))
        {
            if (Pstream::parRun())
            {
                syncTools::syncFaceList(mesh_, resetFlag, maxEqOp<label>());
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
                forAll(pbm, patchi)
                {
                    if (!pbm[patchi].coupled()) continue;
                    forAll(pbm[patchi], i)
                    {
                        const label f = pbm[patchi].start() + i;
                        if (resetFlag[f] == 1) setFaceFraction(f, 0);
                        else if (resetFlag[f] == 2) setFaceFraction(f, 1);
                    }
                }
            }
            // the wet-face list no longer describes alphaf: the full pass
            wetFacesNFaces_ = -1;
            formAreaNormal();
        }
    }

    for (const label celli : touchedCell_)
    {
        if (m[celli] == vector::zero) continue;

        const labelList& cf = cells[celli];
        scalar lo = GREAT, hi = -GREAT;
        forAll(cf, i)
        {
            if (!hasValue_[cf[i]]) continue;
            const scalar af = faceFraction(cf[i]);
            lo = min(lo, af);
            hi = max(hi, af);
        }
        if (lo == hi)
        {
            m[celli] = Zero;
        }
        else if (mag(m[celli]) > davofInit::wispThreshold(mesh_, celli, wispTol_))
        {
            const bool inside = alphaInside(celli);
            if (!inside) ++nRuleOutside_;
            if (alphaOnly)
            {
                // the set is formed from alpha below; the count of the
                // cells with wet faces and alpha outside is kept
                continue;
            }
            if (interfaceCellRule_ == "normalOnly" || inside)
            {
                if (!inside) ++nRuleOutsidePlane_;
                interfaceCell_[celli] = true;
                interfaceLabels_.append(celli);
            }
        }
    }
    if (alphaOnly)
    {
        const label nCells = mesh_.nCells();
        for (label celli = 0; celli < nCells; ++celli)
        {
            if (alphaInside(celli))
            {
                interfaceCell_[celli] = true;
                interfaceLabels_.append(celli);
            }
        }
    }

    // The faces of one cell in ascending order: primitiveMesh::cells()
    // lists the owned faces first, then the faces the cell neighbours,
    // while the face moment used to be accumulated by one pass over all
    // faces in face order; the same order keeps its partial sums.
    DynamicList<label> sortedFaces(16);

    // The sliver flag, the explicit position and the realizability defects
    // of the explicit plane (always from m_c; cutCell_ cuts the explicit
    // plane, not the positioned one). The models that read the planes
    // positioned by alpha along m_c (exact, ringFit) get them here, in a
    // first pass; for the other models every plane is positioned once,
    // below, along the normal the model selects. The final plane is the
    // same either way: vofCutCell starts from the cell's vertex values
    // and carries nothing over from an earlier cut, so positioning along
    // the selected normal once gives the plane the former two passes
    // ended with (the second pass re-positioned along that normal, the
    // first pass's plane along m_c was only overwritten).
    const bool twoPass = normalModel_->needsPositionedPlanes();

    // The curvature tensor's per-cell records (requested by the face
    // update's model only; otherwise untouched)
    usedNormal_.resize(interfaceLabels_.size());
    if (tensorRequested_)
    {
        polyMoment_.resize(interfaceLabels_.size());
    }

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        const scalar mMag = mag(m[celli]);
        // mMag > 0 on every path but alphaOnly (the interface cells are
        // selected by |m_c| above a positive threshold otherwise)
        const vector n = (mMag > 0) ? m[celli]/mMag : vector(Zero);

        if (normalModel_->isSliver(celli, mMag))
        {
            sliverFlag_[celli] = 1;
        }

        // The face moment sum_f alpha_f (x_f - x_c).S_f^out of the explicit
        // position, over this cell's faces only (it is read for interface
        // cells only), in face order with the per-face terms of the former
        // pass over all faces. Computed before the first-pass positioning
        // (2026-11-04: positioning ccs starts from pPlane_); the two do not
        // depend on each other, so the order leaves vofCutCell's result as
        // it was.
        sortedFaces.clear();
        sortedFaces.append(cells[celli]);
        Foam::sort(sortedFaces);
        scalar faceMoment = 0;
        for (const label facei : sortedFaces)
        {
            // A face of an empty patch carries no value (the one-cell-thick
            // mesh of the 2D tests): its cross-section is the cell's, so
            // its fraction is alpha_c and its term closes the Gauss sum
            // (the fork's davof-isoAdvector 3ff7ac7, 2026-11-07; a closed
            // cell has no such face and keeps its sum bit for bit)
            const bool empty = !hasValue_[facei];
            if (empty && !emptyFace_[facei]) continue;
            const scalar af = empty ? alpha1_[celli] : faceFraction(facei);
            if (af == 0) continue;
            // faceState vector: (x_f - x_c).a_f, the exact first moment of
            // the fan face (every fan triangle contains x_f)
            const scalar term =
                (wetVecPtr_ && !empty)
              ? ((Cf[facei] - C[celli]) & faceWettedArea(facei))
              : af*((Cf[facei] - C[celli]) & Sf[facei]);
            if (own[facei] == celli)
            {
                faceMoment += term;
            }
            else
            {
                faceMoment -= term;
            }
        }

        const scalar p =
            (mMag > 0) ? (faceMoment - 3.0*alpha1_[celli]*V[celli])/mMag : 0;
        pPlane_[celli] = p;

        if (twoPass)
        {
            const bool ok = positionPlane(celli, m[celli]);
            if (tensorRequested_)
            {
                polyMoment_[i] = ok ? cutPolygonMoments() : symmTensor(Zero);
            }
        }

        if (mMag > 0)
        {
            const label status = cutCell_.calcSubCell(celli, p, n);
            const scalar alphaPlane =
                (status == 0) ? cutCell_.VolumeOfFluid() : ((status < 0) ? 1 : 0);
            const vector areaPlane = (status == 0) ? cutCell_.faceArea() : vector(Zero);

            planeVolDefect_[celli] = mag(alphaPlane - alpha1_[celli]);
            planeAreaDefect_[celli] =
                mag(areaPlane - m[celli])/Foam::pow(V[celli], 2.0/3.0);
            maxPlaneVolDefect_ = max(maxPlaneVolDefect_, planeVolDefect_[celli]);
            maxPlaneAreaDefect_ = max(maxPlaneAreaDefect_, planeAreaDefect_[celli]);
        }
        else
        {
            // no explicit plane without an area normal (alphaOnly)
            planeVolDefect_[celli] = 0;
            planeAreaDefect_[celli] = 0;
        }
    }

    if (twoPass)
    {
        // The boundary values of the planes for the ring models' exchange
        normal_.correctBoundaryConditions();
        centre_.correctBoundaryConditions();
    }

    // The model's normals; the plane of every interface cell positioned by
    // alpha along the selected normal (nUse[i] is m_c itself where the
    // model does not replace it), or, with two passes, re-positioned where
    // the normal changed
    vectorField nUse;
    boolList changed;
    normalModel_->select(tState, m, nUse, changed);

    forAll(interfaceLabels_, i)
    {
        if (!twoPass || changed[i])
        {
            const bool ok = positionPlane(interfaceLabels_[i], nUse[i]);
            if (tensorRequested_)
            {
                polyMoment_[i] = ok ? cutPolygonMoments() : symmTensor(Zero);
            }
        }
        const scalar nMag = mag(nUse[i]);
        usedNormal_[i] = (nMag > 0) ? nUse[i]/nMag : vector(Zero);
    }

    maxPlaneVolDefect_ = returnReduce(maxPlaneVolDefect_, maxOp<scalar>());
    maxPlaneAreaDefect_ = returnReduce(maxPlaneAreaDefect_, maxOp<scalar>());
    const label nInterface = returnReduce(interfaceLabels_.size(), sumOp<label>());

    normal_.correctBoundaryConditions();
    centre_.correctBoundaryConditions();
    areaNormal_.correctBoundaryConditions();
    sliverFlag_.correctBoundaryConditions();
    pPlane_.correctBoundaryConditions();
    planeVolDefect_.correctBoundaryConditions();
    planeAreaDefect_.correctBoundaryConditions();

    Info<< "davof: " << nInterface << " interface cells (wispTol " << wispTol_
        << "), ";
    normalModel_->report(Info);
    if (normalModel_->type() == "exact")
    {
        Info<< " at t = " << tState;
    }
    Info<< ", explicit plane defects max |alpha_plane - alpha| = "
        << maxPlaneVolDefect_ << ", max |A n - m|/V^(2/3) = "
        << maxPlaneAreaDefect_;
    if (positioning_ == "ccs")
    {
        label nC = nPosCells_, nE = nPosEvals_, mE = maxPosEvals_;
        label nF = nPosFloor_, nM = nPosMaxIter_, nEnd = nPosEnd_, nEndNo = nPosEndNoPlane_;
        scalar rMax = maxPosRes_, rSum = sumPosRes_;
        labelList hist(posEvalHist_);
        if (Pstream::parRun())
        {
            reduce(nC, sumOp<label>()); reduce(nE, sumOp<label>());
            reduce(mE, maxOp<label>()); reduce(nF, sumOp<label>());
            reduce(nM, sumOp<label>()); reduce(nEnd, sumOp<label>());
            reduce(nEndNo, sumOp<label>());
            reduce(rMax, maxOp<scalar>()); reduce(rSum, sumOp<scalar>());
            forAll(hist, k) reduce(hist[k], sumOp<label>());
        }
        Info<< "; positioning ccs: " << nC << " positionings, cuts per cell "
            << scalar(nE)/max(nC, label(1)) << " (max " << mE
            << ", histogram 1.." << hist.size() - 1 << "+ " << hist
            << "), residual max " << rMax << " mean " << rSum/max(nC, label(1))
            << ", floor exits " << nF << ", maxIter exits " << nM
            << ", end-of-interval cells " << nEnd << " (" << nEndNo
            << " without a plane)";
    }
    {
        label nO = nRuleOutside_, nOP = nRuleOutsidePlane_;
        label nRc = nRuleResetCells_, nRf = nRuleResetFaces_;
        scalar dMax = ruleResetMaxDaf_;
        if (Pstream::parRun())
        {
            reduce(nO, sumOp<label>()); reduce(nOP, sumOp<label>());
            reduce(nRc, sumOp<label>()); reduce(nRf, sumOp<label>());
            reduce(dMax, maxOp<scalar>());
        }
        Info<< "; interface-cell rule " << interfaceCellRule_
            << ": cells with |m_c| above the threshold and alpha outside (tol, 1 - tol) "
            << nO << ", with a plane " << nOP;
        if (interfaceCellRule_ == "normalAndAlphaReset")
        {
            Info<< ", reset cells " << nRc << ", faces " << nRf
                << ", max |dalphaf| " << dMax;
        }
    }
    Info<< endl;
}


Foam::symmTensor Foam::reconstruction::davof::cutPolygonMoments()
{
    // The second moments (1/A) int (x - c)(x - c)^T dA of the cut polygon
    // about its centroid c, by the fan of triangles (c, p_i, p_{i+1}): for
    // a triangle with the vertices 0, v1, v2 the integral of x x^T is
    // (A_t/12)(2 v1 v1^T + 2 v2 v2^T + v1 v2^T + v2 v1^T)
    const bool ccs = (positioning_ == "ccs");
    const DynamicList<point>& pts = ccs ? cutCellPos_.facePoints() : sIterPLIC_.facePoints();
    const point& c = ccs ? cutCellPos_.faceCentre() : sIterPLIC_.surfaceCentre();
    const label n = pts.size();

    symmTensor M(Zero);
    scalar A = 0;
    for (label i = 0; i < n; ++i)
    {
        const vector v1 = pts[i] - c;
        const vector v2 = pts[(i + 1) % n] - c;
        const scalar At = 0.5*mag(v1 ^ v2);
        M += (At/6.0)*(sqr(v1) + sqr(v2) + symm(v1*v2));
        A += At;
    }
    if (A > 0)
    {
        return M/A;
    }
    return symmTensor(Zero);
}


void Foam::reconstruction::davof::updateCurvatureTensor()
{
    addProfilingInFunction(geometricVoF);

    if (!tensorRequested_ || polyMoment_.size() != interfaceLabels_.size())
    {
        FatalErrorInFunction
            << "The curvature tensor was not requested before the"
            << " reconstruction (requestCurvatureTensor(true))."
            << exit(FatalError);
    }

    if (!curvatureTensorPtr_)
    {
        const IOobject::writeOption wo =
            writeCurvatureTensor_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE;
        curvatureTensorPtr_.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    IOobject::groupName("curvatureTensorDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    wo
                ),
                mesh_,
                dimensionedSymmTensor(inv(dimLength), Zero)
            )
        );
        curvatureOffsetPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("curvatureOffsetDavof", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    wo
                ),
                mesh_,
                dimensionedScalar(dimLength, Zero)
            )
        );
        mRingPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("davofRingTensor", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(dimArea, Zero)
            )
        );
    }

    volSymmTensorField& Bf = *curvatureTensorPtr_;
    volScalarField& z0f = *curvatureOffsetPtr_;
    symmTensorField& B = Bf.primitiveFieldRef();
    scalarField& z0 = z0f.primitiveFieldRef();

    // Zero at the cells of the previous update (zero elsewhere)
    for (const label celli : tensorCells_)
    {
        B[celli] = Zero;
        z0[celli] = 0;
    }
    tensorCells_ = interfaceLabels_;

    // The ring exchange (as the ring models): m_c, the plane centroids and
    // normals of the interface cells across the processors
    davofNormalModel::fillInterfaceNormals(*this, *mRingPtr_, mRingFilled_);
    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    const boolList& zone = interfaceCell_;
    exchange.setUpCommforZone(zone, true);

    // The centroids: the cut polygons' (centre_) or the ones isoAdvection
    // re-anchored for this update (tensorCentroid anchored)
    const bool useAnchored = (tensorCentresValid_ && tensorCentrePtr_);
    const volVectorField& xField = useAnchored ? *tensorCentrePtr_ : centre_;
    tensorCentresValid_ = false;

    const Map<vector> mapM(exchange.getDatafromOtherProc(zone, *mRingPtr_));
    const Map<vector> mapX(exchange.getDatafromOtherProc(zone, xField));
    const Map<vector> mapN(exchange.getDatafromOtherProc(zone, normal_));
    const labelListList& stencil = exchange.getStencil();

    const vectorField& m = areaNormal_.primitiveField();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();

    DynamicList<vector> xs(32);
    DynamicList<vector> ns(32);
    DynamicList<scalar> ws(32);

    nTensorCells_ = interfaceLabels_.size();
    nTensorFallback_ = 0;

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        // The ring's positioned planes, weighted by |m_k|
        xs.clear(); ns.clear(); ws.clear();
        scalar wSum = 0;
        for (const label gblIdx : stencil[celli])
        {
            const vector mk = exchange.getValue(*mRingPtr_, mapM, gblIdx);
            const scalar mkMag = mag(mk);
            if (!(mkMag > 0)) continue;       // not an interface cell
            const vector nk = exchange.getValue(normal_, mapN, gblIdx);
            const scalar nkMag = mag(nk);
            if (!(nkMag > 0)) continue;       // no positioned plane
            const vector xk = exchange.getValue(xField, mapX, gblIdx);
            xs.append(xk);
            ns.append(nk/nkMag);
            ws.append(mkMag);
            wSum += mkMag;
        }

        // The frame of the used normal (into the liquid) at the cell's
        // plane centroid
        vector e3 = usedNormal_[i];
        if (!(mag(e3) > 0))
        {
            const scalar mMag = mag(m[celli]);
            e3 = (mMag > 0) ? m[celli]/mMag : vector(Zero);
        }
        const point o = (mag(normal_[celli]) > 0) ? xField[celli] : C[celli];
        const scalar h = Foam::cbrt(V[celli]);

        vector e1, e2;
        scalar q[6];
        if
        (
            !(mag(e3) > 0)
         || xs.size() < tensorMinCells_
         || !(wSum > 0)
         || !davofNormalModel::mongeFit
            (
                e3, o, h, xs, ns, ws, wSum, tensorAlignTol_, tensorMinCells_,
                q, e1, e2
            )
        )
        {
            ++nTensorFallback_;
            continue;
        }

        // The Hessian of z = h (q0 + q1 x/h + ... + q3 x^2/h^2 + q4 x y/h^2
        // + q5 y^2/h^2) in the frame (e1, e2): [[2 q3, q4], [q4, 2 q5]]/h,
        // mapped to the global frame; the vertex offset from the volume
        // match over the cut polygon
        const scalar B11 = 2.0*q[3]/h;
        const scalar B12 = q[4]/h;
        const scalar B22 = 2.0*q[5]/h;
        const symmTensor Bc =
            B11*sqr(e1) + (2.0*B12)*symm(e1*e2) + B22*sqr(e2);
        B[celli] = Bc;
        z0[celli] = -0.5*(Bc && polyMoment_[i]);
    }

    nTensorCells_ = returnReduce(nTensorCells_, sumOp<label>());
    nTensorFallback_ = returnReduce(nTensorFallback_, sumOp<label>());

    Bf.correctBoundaryConditions();
    z0f.correctBoundaryConditions();

    Info<< "davof: curvature tensor on " << nTensorCells_
        << " cells, " << nTensorFallback_ << " tensor fallbacks"
        << (useAnchored ? " (anchored centroids)" : "") << endl;
}


Foam::vectorField& Foam::reconstruction::davof::tensorCentres()
{
    if (!tensorCentrePtr_)
    {
        tensorCentrePtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("davofTensorCentre", alpha1_.group()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(dimLength, Zero)
            )
        );
    }
    tensorCentresValid_ = true;
    return tensorCentrePtr_->primitiveFieldRef();
}


// ************************************************************************* //


void Foam::reconstruction::davof::setFaceFraction(const label facei, const scalar v)
{
    const label nInt = mesh_.nInternalFaces();
    if (facei < nInt)
    {
        alphaf_.primitiveFieldRef()[facei] = v;
        return;
    }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label patchi = pbm.whichPatch(facei);
    alphaf_.boundaryFieldRef()[patchi][facei - pbm[patchi].start()] = v;
}


// * * * * * * * * the iterative face-fraction corrector (2026-11-09) * * * * //

void Foam::reconstruction::davof::inducedFaceFractions
(
    scalarField& afNew,
    bitSet& set,
    scalar& dMax,
    scalar& dSum,
    label& nSet,
    label& nNoPlane
)
{
    addProfilingInFunction(geometricVoF);
    const label nFaces = mesh_.nFaces();
    const label nInt = mesh_.nInternalFaces();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& aIn = alpha1_.primitiveField();
    const vectorField& nIn = normal_.primitiveField();
    const vectorField& cIn = centre_.primitiveField();
    const scalar tolA = surfCellTol_;

    afNew.resize(nFaces);
    set.resize(nFaces);
    set.reset();
    for (label facei = 0; facei < nFaces; ++facei)
    {
        afNew[facei] = hasValue_[facei] ? faceFraction(facei) : 0;
    }

    // The clip of a face polygon by a plane: the wet area over the face
    // area (cutFacePLIC, H. Scheufler, J. Roenby, J. Comput. Phys. 383
    // (2019) 1-23, doi:10.1016/j.jcp.2019.01.009; the plane's normal points
    // into the liquid, as the cut assumes: (x - base).n > 0 is liquid).
    // The UNIT normal, as cutCellPLIC passes it: calcSubFace snaps a point's
    // signed value (x - base).n to zero below SMALL (1e-15, absolute) and
    // calls a face without a positive point a gas face; with the cut
    // polygon's area vector normal_ (|normal_| down to 1e-12 for a sliver
    // cut) every point of a fully wet face snapped and the clip returned 0,
    // which gave the 0.5 changes of the plane gate of 2026-11-09 (located
    // 2026-11-10 by the unit check on the exact planar state).
    cutFacePLIC cutF(mesh_);
    auto clip = [&](const label facei, const vector& nrm, const point& ctr) -> scalar
    {
        cutF.calcSubFace(facei, nrm/mag(nrm), ctr);
        const scalar af = mag(cutF.subFaceArea())/mag(Sf[facei]);
        return min(max(af, scalar(0)), scalar(1));
    };
    // The contribution of a cell to the mean on a face: its plane's clip
    // when it has one, its (clamped) alpha otherwise; an interface cell
    // without a plane (the ring sum vanished, or the positioning failed)
    // is counted
    label nNoPlaneLocal = 0;
    auto contrib = [&]
    (
        const label facei,
        const bool interfaceCell,
        const vector& nrm,
        const point& ctr,
        const scalar a
    ) -> scalar
    {
        if (mag(nrm) > 0)
        {
            return clip(facei, nrm, ctr);
        }
        if (interfaceCell)
        {
            ++nNoPlaneLocal;
        }
        return min(max(a, scalar(0)), scalar(1));
    };

    for (label facei = 0; facei < nInt; ++facei)
    {
        const label o = own[facei], n = nei[facei];
        const bool io = interfaceCell_[o], in = interfaceCell_[n];
        if (!(io || in)) continue;
        afNew[facei] =
            0.5*
            (
                contrib(facei, io, nIn[o], cIn[o], aIn[o])
              + contrib(facei, in, nIn[n], cIn[n], aIn[n])
            );
        set.set(facei);
    }

    // The coupled patches: the other side's plane and alpha from the
    // patch-neighbour fields (the remote cell's values), its interface
    // flag from its alpha (the alphaOnly rule, the same on both sides);
    // the non-coupled patches are left as the transport set them
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    bitSet coupledSet(nFaces);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled()) continue;
        if (alphaf_.boundaryField()[patchi].size() != pp.size()) continue;
        const coupledFvPatchField<vector>& pn =
            refCast<const coupledFvPatchField<vector>>(normal_.boundaryField()[patchi]);
        const coupledFvPatchField<vector>& pc =
            refCast<const coupledFvPatchField<vector>>(centre_.boundaryField()[patchi]);
        const coupledFvPatchField<scalar>& pa =
            refCast<const coupledFvPatchField<scalar>>(alpha1_.boundaryField()[patchi]);
        const tmp<vectorField> tnN(pn.patchNeighbourField());
        const tmp<vectorField> tcN(pc.patchNeighbourField());
        const tmp<scalarField> taN(pa.patchNeighbourField());
        const vectorField& nN = tnN();
        const vectorField& cN = tcN();
        const scalarField& aN = taN();
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label o = own[facei];
            const bool io = interfaceCell_[o];
            const bool in = (aN[i] > tolA && aN[i] < 1 - tolA);
            if (!(io || in)) continue;
            afNew[facei] =
                0.5*
                (
                    contrib(facei, io, nIn[o], cIn[o], aIn[o])
                  + contrib(facei, in, nN[i], cN[i], aN[i])
                );
            set.set(facei);
            coupledSet.set(facei);
        }
    }
    if (Pstream::parRun() || coupledSet.any())
    {
        // one value per coupled face: the half-sum of the two sides' means
        // (the sum is symmetric bit for bit; the two clips of one face by
        // one plane differ at round-off through the face's orientation)
        scalarList both(nFaces, Zero);
        for (const label facei : coupledSet)
        {
            both[facei] = afNew[facei];
        }
        syncTools::syncFaceList(mesh_, both, plusEqOp<scalar>());
        for (const label facei : coupledSet)
        {
            afNew[facei] = 0.5*both[facei];
        }
    }

    dMax = 0;
    dSum = 0;
    scalar nSetW = 0;
    label worst = -1;
    for (const label facei : set)
    {
        const scalar d = mag(afNew[facei] - faceFraction(facei));
        const scalar w = coupledSet.test(facei) ? 0.5 : 1.0;
        if (d > dMax)
        {
            dMax = d;
            worst = facei;
        }
        dSum += w*d;
        nSetW += w;
    }
    const scalar dMaxLocal = dMax;
    reduce(dMax, maxOp<scalar>());
    reduce(dSum, sumOp<scalar>());
    reduce(nSetW, sumOp<scalar>());
    nSet = label(nSetW + 0.5);
    nNoPlane = returnReduce(nNoPlaneLocal, sumOp<label>());

    // The face of the largest change above the tolerance (correctorLog;
    // printed by the rank that holds it): the located cell pair, not a
    // guess (2026-11-10)
    if (correctorLog_ && worst >= 0 && dMaxLocal == dMax && dMax > correctorTolerance_)
    {
        const label o = own[worst];
        const bool internal = (worst < nInt);
        Pout<< "davof corrector: worst face " << worst
            << (internal ? "" : " (boundary)") << " centre " << mesh_.faceCentres()[worst]
            << " alphaf " << faceFraction(worst) << " -> " << afNew[worst]
            << "; owner " << o << " alpha " << aIn[o] << " interface " << interfaceCell_[o]
            << " |n| " << mag(nIn[o]);
        if (internal)
        {
            const label n = nei[worst];
            Pout<< "; neighbour " << n << " alpha " << aIn[n] << " interface "
                << interfaceCell_[n] << " |n| " << mag(nIn[n]);
            if (mag(nIn[o]) > 0) Pout<< "; clip owner " << clip(worst, nIn[o], cIn[o]);
            if (mag(nIn[n]) > 0) Pout<< " clip neighbour " << clip(worst, nIn[n], cIn[n]);
        }
        Pout<< endl;
    }
}


bool Foam::reconstruction::davof::correctFaceFractions()
{
    addProfilingInFunction(geometricVoF);
    if (faceFractionCorrector_ == "none")
    {
        return false;
    }

    // the remote cells' alpha on the processor patches (the bounding may
    // have changed alpha after the last evaluation)
    if (Pstream::parRun())
    {
        alpha1_.correctBoundaryConditions();
    }

    // iteration 0: the planes of the predictor state (the transported
    // face fractions; the advection's wet-face list describes them)
    reconstruct(true);

    corrIterationsUsed_ = 0;
    corrMaxChange_.clear();
    corrL1Change_.clear();
    corrFaces_ = 0;
    corrNoPlaneCells_ = 0;
    bool changedAny = false;

    scalarField afNew;
    bitSet set;
    for (label it = 0; it < correctorIterations_; ++it)
    {
        scalar dMax = 0, dSum = 0;
        label nSet = 0, nNoPlane = 0;
        inducedFaceFractions(afNew, set, dMax, dSum, nSet, nNoPlane);
        corrMaxChange_.append(dMax);
        corrL1Change_.append(dSum/max(nSet, label(1)));
        corrFaces_ = nSet;
        corrNoPlaneCells_ = max(corrNoPlaneCells_, nNoPlane);
        if (dMax <= correctorTolerance_)
        {
            break;
        }
        for (const label facei : set)
        {
            setFaceFraction(facei, afNew[facei]);
        }
        ++corrIterationsUsed_;
        changedAny = true;
        // the wet-face list no longer describes alphaf: the full pass
        wetFacesNFaces_ = -1;
        reconstruct(true);
    }
    correctedTimeIndex_ = mesh_.time().timeIndex();
    corrInterfaceCells_ = returnReduce(interfaceLabels_.size(), sumOp<label>());

    // The identity sum_k m_k = sum_b alpha_f S_f over the non-coupled
    // boundary faces (the internal and coupled faces cancel when every face
    // carries ONE value; m_k here is the code's sum_f alpha_f S_f^out, into
    // the liquid, i.e. minus the area normal of the notation)
    {
        const vectorField& m = areaNormal_.primitiveField();
        const vectorField& Sf = mesh_.faceAreas();
        vector sumM(Zero), sumB(Zero);
        scalar sumMag = 0;
        forAll(m, celli)
        {
            sumM += m[celli];
            sumMag += mag(m[celli]);
        }
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        forAll(pbm, patchi)
        {
            if (pbm[patchi].coupled()) continue;
            const fvsPatchScalarField& pf = alphaf_.boundaryField()[patchi];
            forAll(pf, i)
            {
                sumB += pf[i]*Sf[pbm[patchi].start() + i];
            }
        }
        reduce(sumM, sumOp<vector>());
        reduce(sumB, sumOp<vector>());
        reduce(sumMag, sumOp<scalar>());
        corrIdentity_ = mag(sumM - sumB)/max(sumMag, VSMALL);
    }

    if (correctorLog_)
    {
        Info<< "davof corrector (" << faceFractionCorrector_ << "): iterations "
            << corrIterationsUsed_ << " of " << correctorIterations_
            << ", interface cells (alpha) " << corrInterfaceCells_
            << ", faces " << corrFaces_
            << ", cells without a plane " << corrNoPlaneCells_
            << "; max |dalphaf| per iteration (";
        forAll(corrMaxChange_, k)
        {
            Info<< (k ? " " : "") << corrMaxChange_[k];
        }
        Info<< "), L1 |dalphaf| (";
        forAll(corrL1Change_, k)
        {
            Info<< (k ? " " : "") << corrL1Change_[k];
        }
        Info<< "); identity |sum_k m_k - sum_b alphaf S_f| / sum_k |m_k| = "
            << corrIdentity_ << endl;
    }
    return changedAny;
}
