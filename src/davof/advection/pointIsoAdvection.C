/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2022 OpenCFD Ltd.
    Modified code Copyright (C) 2019-2020 DLR
    Modified code Copyright (C) 2018, 2021 Johan Roenby
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

#include "pointIsoAdvection.H"
#include "davof.H"
#include "coupledFvPatchField.H"
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
#include "mathematicalConstants.H"
#include "pointIsoPlaneFit.H"
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointIsoAdvection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointIsoAdvection::pointIsoAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    mesh_(alpha1.mesh()),
    dict_(mesh_.solverDict(alpha1.name())),
    alpha1_(alpha1),
    alpha1In_(alpha1.primitiveFieldRef()),
    phi_(phi),
    U_(U),
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
    alphaPhi_
    (
        IOobject
        (
            "alphaPhi_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol/dimTime, Zero)
    ),
    advectionTime_(0),
    timeIndex_(-1),

    // Tolerances and solution controls
    nAlphaBounds_(dict_.getOrDefault<label>("nAlphaBounds", 3)),
    isoFaceTol_(dict_.getOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(dict_.getOrDefault<scalar>("surfCellTol", 1e-8)),
    writeIsoFacesToFile_(dict_.getOrDefault("writeIsoFaces", false)),

    // The face fractions (before surf_)
    faceFractions_(dict_.getOrDefault<bool>("faceFractions", false)),
    wetFraction_(makeWetFraction(alpha1, faceFractions_)),

    // Cell cutting data
    surfCells_(label(0.2*mesh_.nCells())),
    surf_(reconstructionSchemes::New(alpha1_, phi_, U_, dict_)),
    advectFace_(alpha1.mesh(), alpha1),
    bsFaces_(label(0.2*mesh_.nBoundaryFaces())),
    bsx0_(bsFaces_.size()),
    bsn0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),

    // Porosity
    porosityEnabled_(dict_.getOrDefault<bool>("porosityEnabled", false)),
    porosityPtr_(nullptr),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0),

    // pointIsoAdvection
    pointVelocity_(dict_.getOrDefault<word>("pointVelocity", "edgeLinear")),
    compareIso_(dict_.getOrDefault<bool>("pointIsoCompare", false)),
    printPoints_(dict_.getOrDefault<bool>("pointIsoPrintPoints", false)),
    bndPoint_(mesh_.nPoints()),
    bndEdge_(mesh_.nEdges()),
    afSet_(),
    wetFaces_(),
    davofPtr_(nullptr),
    pointU_()
{
    davofPtr_ = dynamic_cast<reconstruction::davof*>(&surf_());
    if (faceFractions_ && !davofPtr_)
    {
        FatalErrorInFunction
            << "faceFractions needs the davof reconstruction, the scheme is "
            << surf_().type() << exit(FatalError);
    }
   cutFaceAdvect::debug = debug;

    if (pointVelocity_ != "edgeLinear")
    {
        FatalErrorInFunction
            << "pointVelocity " << pointVelocity_ << " is not implemented;"
            << " only edgeLinear" << exit(FatalError);
    }

    // The boundary points: the points and the edges of the boundary faces of
    // the patches that prescribe a velocity (walls, inlets, fixedValue and
    // the like; not empty, wedge, symmetry, cyclic or processor patches)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const volVectorField::Boundary& Ub = U_.boundaryField();
        forAll(pbm, patchi)
        {
            if (pbm[patchi].empty() || !Ub[patchi].fixesValue()) continue;

            const label start = pbm[patchi].start();
            forAll(pbm[patchi], pfi)
            {
                const label facei = start + pfi;
                for (const label pointi : mesh_.faces()[facei])
                {
                    bndPoint_.set(pointi);
                }
                for (const label edgei : mesh_.faceEdges()[facei])
                {
                    bndEdge_.set(edgei);
                }
            }
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

    // Reading porosity properties from constant directory
    IOdictionary porosityProperties
    (
        IOobject
        (
            "porosityProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    porosityEnabled_ =
        porosityProperties.getOrDefault<bool>("porosityEnabled", false);

    if (porosityEnabled_)
    {
        FatalErrorInFunction
            << "pointIsoAdvection does not support porosity"
            << exit(FatalError);
    }

    if (porosityEnabled_)
    {
        porosityPtr_ = mesh_.getObjectPtr<volScalarField>("porosity");

        if (porosityPtr_)
        {
            auto limits = gMinMax(porosityPtr_->primitiveField());

            if (limits.min() <= 0 || limits.max() > 1 + SMALL)
            {
                FatalErrorInFunction
                    << "Porosity field has values <= 0 or > 1"
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Porosity enabled in constant/porosityProperties "
                << "but no porosity field is found in object registry."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::surfaceScalarField Foam::pointIsoAdvection::makeWetFraction
(
    const volScalarField& alpha1,
    const bool faceFractions
)
{
    const fvMesh& mesh = alpha1.mesh();
    if (faceFractions)
    {
        IOobject io
        (
            IOobject::groupName("alphaf", alpha1.group()),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (io.typeHeaderOk<surfaceScalarField>(true))
        {
            return surfaceScalarField(io, mesh);
        }
        FatalErrorInFunction
            << "faceFractions: the field " << io.name() << " is not in "
            << mesh.time().timeName() << "/ (initDavofFields writes it)"
            << exit(FatalError);
    }
    return surfaceScalarField
    (
        IOobject
        (
            "pointIsoWetFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
}


void Foam::pointIsoAdvection::defaultFaceFractions()
{
    const scalarField& a = alpha1In_;
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const scalarField& phiIn = phi_.primitiveField();
    scalarField& afIn = wetFraction_.primitiveFieldRef();

    // The target of a face no plane set: alpha of its upwind cell after the
    // bounding (the owner's on a boundary face), as davofSwept
    for (label facei = 0; facei < nInt; ++facei)
    {
        if (!afSet_.test(facei))
        {
            afIn[facei] = (phiIn[facei] >= 0) ? a[own[facei]] : a[nei[facei]];
        }
    }

    // A processor face is an internal face of the mesh: its upwind cell is the
    // local one for phi >= 0 (local side), else the remote one, whose new alpha
    // is the patch neighbour field of alpha (the coupled-patch exchange of
    // correctBoundaryConditions). Only the faces of non-coupled patches take
    // the owner's alpha.
    alpha1_.correctBoundaryConditions();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    surfaceScalarField::Boundary& afb = wetFraction_.boundaryFieldRef();
    DynamicLabelList wet;
    for (label facei = 0; facei < nInt; ++facei)
    {
        if (afIn[facei] != 0) wet.append(facei);
    }
    forAll(afb, patchi)
    {
        fvsPatchScalarField& pf = afb[patchi];
        if (pf.size() == 0) continue;
        const label start = pbm[patchi].start();
        const fvPatchScalarField& ap = alpha1_.boundaryField()[patchi];
        const bool coupled = ap.coupled();
        tmp<scalarField> tnbr;
        if (coupled)
        {
            tnbr = refCast<const coupledFvPatchScalarField>(ap).patchNeighbourField();
        }
        const scalarField& phib = phi_.boundaryField()[patchi];
        forAll(pf, i)
        {
            const label facei = start + i;
            if (!afSet_.test(facei))
            {
                pf[i] = (coupled && phib[i] < 0) ? tnbr()[i] : a[own[facei]];
            }
            if (pf[i] != 0) wet.append(facei);
        }
    }
    wetFaces_ = wet;
    if (davofPtr_)
    {
        davofPtr_->setWetFaces(wetFaces_);
    }
}


void Foam::pointIsoAdvection::setProcessorPatches()
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
         && !patches[patchi].empty()
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}


void Foam::pointIsoAdvection::extendMarkedCells
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

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
        }
    }
}


void Foam::pointIsoAdvection::polygonPoints
(
    const label celli,
    const vector& n0,
    const point& x0,
    const scalar dt,
    DynamicPointList& displaced,
    DynamicList<bool>& isBoundary,
    scalar& maxBndDisp
) const
{
    displaced.clear();
    isBoundary.clear();

    const pointField& pts = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const labelList& cEdges = mesh_.cellEdges()[celli];
    const labelList& cPoints = mesh_.cellPoints()[celli];

    // The vertices on the plane (the same zero as cutFaceAdvect's, 10 SMALL)
    for (const label pointi : cPoints)
    {
        scalar g = (pts[pointi] - x0) & n0;
        if (mag(g) < 10*SMALL)
        {
            displaced.append(pts[pointi] + dt*pointU_[pointi]);
            isBoundary.append(bndPoint_.test(pointi));
            if (bndPoint_.test(pointi))
            {
                maxBndDisp = max(maxBndDisp, dt*mag(pointU_[pointi]));
            }
        }
    }

    // The crossings of the edges, u linear along the edge (edgeLinear)
    for (const label edgei : cEdges)
    {
        const edge& e = edges[edgei];
        scalar ga = (pts[e.start()] - x0) & n0;
        scalar gb = (pts[e.end()] - x0) & n0;
        if (mag(ga) < 10*SMALL) ga = 0;
        if (mag(gb) < 10*SMALL) gb = 0;

        if ((ga < 0 && gb > 0) || (ga > 0 && gb < 0))
        {
            const scalar s = ga/(ga - gb);
            const point p
            (
                pts[e.start()] + s*(pts[e.end()] - pts[e.start()])
            );
            const vector u
            (
                (1 - s)*pointU_[e.start()] + s*pointU_[e.end()]
            );
            displaced.append(p + dt*u);
            isBoundary.append(bndEdge_.test(edgei));
            if (bndEdge_.test(edgei))
            {
                maxBndDisp = max(maxBndDisp, dt*mag(u));
            }
        }
    }
}


Foam::scalar Foam::pointIsoAdvection::wettedArea
(
    const label facei,
    const scalar magSf,
    const scalarList& g0,
    const scalarList& g1,
    const scalar tau
) const
{
    const face& f = mesh_.faces()[facei];
    const pointField& pts = mesh_.points();
    const label n = f.size();

    scalarList h(n);
    bool anyNeg = false;
    bool anyPos = false;
    forAll(h, i)
    {
        h[i] = (1 - tau)*g0[i] + tau*g1[i];
        if (h[i] < 0) anyNeg = true;
        if (h[i] > 0) anyPos = true;
    }

    // The liquid is where g_tau < 0
    if (!anyNeg) return 0;      // dry (or the face only touches the plane)
    if (!anyPos) return magSf;  // wet

    // Clip the face polygon to g_tau <= 0
    DynamicPointList poly(n + 2);
    for (label i = 0; i < n; ++i)
    {
        const label j = (i + 1) % n;
        const point& pi = pts[f[i]];
        const point& pj = pts[f[j]];

        if (h[i] <= 0)
        {
            poly.append(pi);
        }
        if ((h[i] < 0 && h[j] > 0) || (h[i] > 0 && h[j] < 0))
        {
            const scalar t = h[i]/(h[i] - h[j]);
            poly.append(pi + t*(pj - pi));
        }
    }

    if (poly.size() < 3) return 0;

    // The area of the clipped polygon: the fan around its vertex average
    point c(Zero);
    for (const point& p : poly) c += p;
    c /= poly.size();

    vector a(Zero);
    forAll(poly, k)
    {
        a += 0.5*((poly[k] - c) ^ (poly[(k + 1) % poly.size()] - c));
    }

    return mag(a);
}


Foam::scalar Foam::pointIsoAdvection::timeIntegratedFaceFlux
(
    const label facei,
    const planePair& pl,
    const scalar dt,
    const scalar phi,
    const scalar magSf,
    scalar& aEnd
) const
{
    const face& f = mesh_.faces()[facei];
    const pointField& pts = mesh_.points();
    const label n = f.size();

    // g_0 and g_1 at the face's vertices (the zero of cutFaceAdvect)
    scalarList g0(n);
    scalarList g1(n);
    forAll(f, i)
    {
        scalar a = (pts[f[i]] - pl.x0) & pl.n0;
        scalar b = (pts[f[i]] - pl.x1) & pl.n1;
        if (mag(a) < 10*SMALL) a = 0;
        if (mag(b) < 10*SMALL) b = 0;
        g0[i] = a;
        g1[i] = b;
    }

    // Break points: 0, 1 and the vertex arrival times inside (0, 1)
    DynamicList<scalar> br(n + 2);
    br.append(0);
    forAll(g0, i)
    {
        const scalar d = g0[i] - g1[i];
        if (d != 0)
        {
            const scalar tauV = g0[i]/d;
            if (tauV > 0 && tauV < 1)
            {
                br.append(tauV);
            }
        }
    }
    br.append(1);
    std::sort(br.begin(), br.end());

    // Two-point Gauss-Legendre on every sub-interval
    const scalar r3 = 1.0/std::sqrt(3.0);
    scalar integral = 0;
    for (label k = 0; k < br.size() - 1; ++k)
    {
        const scalar ta = br[k];
        const scalar tb = br[k + 1];
        if (tb <= ta) continue;

        const scalar mid = 0.5*(ta + tb);
        const scalar half = 0.5*(tb - ta);
        integral +=
            half
           *(
                wettedArea(facei, magSf, g0, g1, mid - half*r3)
              + wettedArea(facei, magSf, g0, g1, mid + half*r3)
            );
    }

    aEnd = wettedArea(facei, magSf, g0, g1, 1.0)/magSf;

    return phi*dt/magSf*integral;
}


void Foam::pointIsoAdvection::timeIntegratedFlux()
{
    addProfilingInFunction(geometricVoF);
    // Get time step
    const scalar dt = mesh_.time().deltaTValue();

    // The mesh-point velocities of the step, once, from the velocity the
    // step's phi belongs to (inverse-distance weights, boundary points
    // carry the boundary condition's value)
    {
        const tmp<pointVectorField> tUp =
            volPointInterpolation::New(mesh_).interpolate(U_);
        pointU_ = tUp().primitiveField();
    }

    // Create object for interpolating velocity to isoface centres (only for
    // the print-only comparison with isoAdvection's face flux)
    autoPtr<interpolationCellPoint<vector>> UInterp;
    if (compareIso_)
    {
        UInterp.reset(new interpolationCellPoint<vector>(U_));
    }

    // For each downwind face of each surface cell we find dVf
    label nSurfaceCells = 0;

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    clearIsoFaceData();

    // The end-of-step wetted fraction of the swept faces
    wetFraction_ = dimensionedScalar(dimless, Zero);
    afSet_.resize(mesh_.nFaces());
    afSet_.reset();

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();
    scalarField& wetIn = wetFraction_.primitiveFieldRef();

    // Get necessary mesh data
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();

    // Storage for isoFace points. Only used if writeIsoFacesToFile_
    DynamicList<List<point>> isoFacePts;
    const labelUList& interfaceLabels = surf_->interfaceLabels();

    // Calculating isoface normal velocity (only for the comparison)
    scalarField Un0(compareIso_ ? interfaceLabels.size() : 0);
    if (compareIso_)
    {
        forAll(Un0, i)
        {
            const label celli = interfaceLabels[i];
            const point x0(surf_->centre()[celli]);
            const vector n0(normalised(-surf_->normal()[celli]));
            Un0[i] = UInterp().interpolate(x0, celli) & n0;
        }
    }

    // Statistics of the fitted planes
    label nPlanes = 0;
    label nBoundaryCells = 0;
    label nOverThree = 0;
    scalar sumResidual = 0;
    scalar maxResidual = 0;
    scalar maxAngle = 0;
    scalar maxDiff = 0;
    label nCompared = 0;
    scalar maxBndDisp = 0;
    scalar maxBndDist = 0;

    DynamicPointList qPts(16);
    DynamicList<bool> qBnd(16);

    // Loop through cells
    forAll(interfaceLabels, i)
    {
        const label celli = interfaceLabels[i];
        if (mag(surf_->normal()[celli]) != 0)
        {
            // This is a surface cell, increment counter, append and mark cell
            nSurfaceCells++;
            surfCells_.append(celli);
            const point x0(surf_->centre()[celli]);
            const vector n0(normalised(-surf_->normal()[celli]));

            // The displaced polygon points and the fitted plane g_1
            polygonPoints(celli, n0, x0, dt, qPts, qBnd, maxBndDisp);

            if (qPts.size() < 3)
            {
                FatalErrorInFunction
                    << "pointIsoAdvection: cell " << celli
                    << " (alpha = " << alpha1In_[celli] << ") has "
                    << qPts.size() << " PLIC polygon points, 3 are needed"
                    << " for the plane fit"
                    << exit(FatalError);
            }

            if (printPoints_)
            {
                Pout<< "pointIsoPoints cell " << celli << " n_k " << n0.x()
                    << ' ' << n0.y() << ' ' << n0.z() << " :";
                forAll(qPts, k)
                {
                    char buf[160];
                    std::snprintf
                    (
                        buf, sizeof(buf), " %.17g %.17g %.17g %d",
                        qPts[k].x(), qPts[k].y(), qPts[k].z(),
                        int(qBnd[k])
                    );
                    Pout<< buf;
                }
                Pout<< endl;
            }

            const pointIso::fit fit =
                pointIso::fitPlane(qPts, qBnd, n0);

            planePair pl;
            pl.n0 = n0;
            pl.x0 = x0;
            pl.n1 = fit.n;
            pl.x1 = fit.x;

            {
                const scalar L = std::cbrt(mesh_.V()[celli]);
                forAll(qPts, k)
                {
                    if (qBnd[k])
                    {
                        maxBndDist =
                            max(maxBndDist, mag(fit.n & (qPts[k] - fit.x))/L);
                    }
                }
                const scalar res = fit.rms/L;
                ++nPlanes;
                sumResidual += res;
                maxResidual = max(maxResidual, res);
                maxAngle =
                    max
                    (
                        maxAngle,
                        std::acos(min(1.0, max(-1.0, fit.n & n0)))
                    );
                if (fit.nBoundary > 0) ++nBoundaryCells;
                if (fit.overThree) ++nOverThree;
            }

            DebugInfo
                << "\n------------ Cell " << celli << " with alpha1 = "
                << alpha1In_[celli] << " and 1-alpha1 = "
                << 1.0 - alpha1In_[celli] << " ------------"
                << endl;

            // Estimate time integrated flux through each downwind face
            // Note: looping over all cell faces - in reduced-D, some of
            //       these faces will be on empty patches
            const cell& celliFaces = cellFaces[celli];
            for (const label facei : celliFaces)
            {
                if (mesh_.isInternalFace(facei))
                {
                    bool isDownwindFace = false;

                    if (celli == own[facei])
                    {
                        if (phiIn[facei] >= 0)
                        {
                            isDownwindFace = true;
                        }
                    }
                    else
                    {
                        if (phiIn[facei] < 0)
                        {
                            isDownwindFace = true;
                        }
                    }

                    if (isDownwindFace)
                    {
                        scalar aEnd = 0;
                        dVfIn[facei] = timeIntegratedFaceFlux
                        (
                            facei,
                            pl,
                            dt,
                            phiIn[facei],
                            magSfIn[facei],
                            aEnd
                        );
                        wetIn[facei] = aEnd;
                        afSet_.set(facei);

                        if (compareIso_)
                        {
                            const scalar ref =
                                advectFace_.timeIntegratedFaceFlux
                                (
                                    facei,
                                    x0,
                                    n0,
                                    Un0[i],
                                    dt,
                                    phiIn[facei],
                                    magSfIn[facei]
                                );
                            maxDiff = max(maxDiff, mag(dVfIn[facei] - ref));
                            ++nCompared;
                        }
                    }

                }
                else
                {
                    bsFaces_.append(facei);
                    bsx0_.append(x0);
                    bsn0_.append(n0);
                    bsUn0_.append(compareIso_ ? Un0[i] : scalar(0));
                    bsPlanes_.append(pl);

                    // Note: we must not check if the face is on the
                    // processor patch here.
                }
            }
        }
    }

    // Get references to boundary fields
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();
    surfaceScalarField::Boundary& wetb = wetFraction_.boundaryFieldRef();

    // Loop through boundary surface faces
    forAll(bsFaces_, i)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[i];
        const label patchi = boundaryMesh.patchID(facei);

        if (!phib[patchi].empty())
        {
            const label patchFacei = boundaryMesh[patchi].whichFace(facei);

            const scalar phiP = phib[patchi][patchFacei];

            if (phiP >= 0)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                scalar aEnd = 0;
                dVfb[patchi][patchFacei] = timeIntegratedFaceFlux
                (
                    facei,
                    bsPlanes_[i],
                    dt,
                    phiP,
                    magSf,
                    aEnd
                );
                wetb[patchi][patchFacei] = aEnd;
                afSet_.set(facei);

                if (compareIso_)
                {
                    const scalar ref =
                        advectFace_.timeIntegratedFaceFlux
                        (
                            facei,
                            bsx0_[i],
                            bsn0_[i],
                            bsUn0_[i],
                            dt,
                            phiP,
                            magSf
                        );
                    maxDiff =
                        max(maxDiff, mag(dVfb[patchi][patchFacei] - ref));
                    ++nCompared;
                }

                // Handling upwind cyclic boundary patches
                const polyPatch& pp = boundaryMesh[patchi];
                const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);
                if (cpp)
                {
                    const label neiPatchID(cpp->neighbPolyPatchID());
                    dVfb[neiPatchID][patchFacei] = -dVfb[patchi][patchFacei];
                }

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }

    // Synchronize processor patches (with the face fractions of the key
    // faceFractions)
    syncProcPatches(dVf_, phi_, false, faceFractions_);

    writeIsoFaces(isoFacePts);

    DebugInfo << "Number of isoAdvector surface cells = "
        << returnReduce(nSurfaceCells, sumOp<label>()) << endl;

    // The report of the step: the fitted planes
    {
        reduce(nPlanes, sumOp<label>());
        reduce(nBoundaryCells, sumOp<label>());
        reduce(nOverThree, sumOp<label>());
        reduce(sumResidual, sumOp<scalar>());
        reduce(maxResidual, maxOp<scalar>());
        reduce(maxAngle, maxOp<scalar>());
        reduce(maxBndDisp, maxOp<scalar>());
        reduce(maxBndDist, maxOp<scalar>());

        Info<< "pointIsoAdvection: planes " << nPlanes
            << ", fit residual (rms/V^(1/3)) mean "
            << (nPlanes ? sumResidual/nPlanes : scalar(0))
            << " max " << maxResidual
            << ", max angle(n', n_k) "
            << maxAngle*180.0/Foam::constant::mathematical::pi
            << " deg, cells with boundary points " << nBoundaryCells
            << ", with more than three non-collinear boundary points "
            << nOverThree << ", boundary points: max |dx| = " << maxBndDisp
            << " m, max distance to the fitted plane/V^(1/3) = " << maxBndDist
            << endl;
    }

    if (compareIso_)
    {
        scalar phiDtMax = 0;
        forAll(phiIn, facei)
        {
            phiDtMax = max(phiDtMax, mag(phiIn[facei]));
        }
        forAll(phib, patchi)
        {
            const scalarField& pp = phib[patchi];
            forAll(pp, facei)
            {
                phiDtMax = max(phiDtMax, mag(pp[facei]));
            }
        }
        phiDtMax *= dt;

        reduce(maxDiff, maxOp<scalar>());
        reduce(phiDtMax, maxOp<scalar>());
        reduce(nCompared, sumOp<label>());

        Info<< "pointIsoAdvection compare: faces " << nCompared
            << ", max |dVf - dVf_isoAdvection| = " << maxDiff
            << ", max |phi dt| = " << phiDtMax
            << ", ratio = " << maxDiff/(phiDtMax + VSMALL) << endl;
    }
}


void Foam::pointIsoAdvection::setDownwindFaces
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
    for (const label facei: c)
    {
        // Get face and corresponding flux
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


Foam::scalar Foam::pointIsoAdvection::netFlux
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

    for (const label facei : c)
    {
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


Foam::DynamicList<Foam::label>  Foam::pointIsoAdvection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces,
    bool withAlphaf
)
{
    DynamicLabelList syncedFaces(0);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (Pstream::parRun())
    {
        DynamicList<label> neighProcs;
        PstreamBuffers pBufs;

        // Send
        for (const label patchi : procPatchLabels_)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);
            const label nbrProci = procPatch.neighbProcNo();

            // Neighbour connectivity
            neighProcs.push_uniq(nbrProci);

            const scalarField& pFlux = dVf.boundaryField()[patchi];
            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            UOPstream toNbr(nbrProci, pBufs);
            toNbr << surfCellFacesOnProcPatch << dVfPatch;
            if (withAlphaf)
            {
                const UIndirectList<scalar> afPatch
                (
                    wetFraction_.boundaryField()[patchi],
                    surfCellFacesOnProcPatch
                );
                toNbr << afPatch;
            }
        }

        // Limited to involved neighbour procs
        pBufs.finishedNeighbourSends(neighProcs);


        // Receive and combine
        for (const label patchi : procPatchLabels_)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);
            const label nbrProci = procPatch.neighbProcNo();

            List<label> faceIDs;
            List<scalar> nbrdVfs;
            List<scalar> nbrAf;

            {
                UIPstream fromNbr(nbrProci, pBufs);
                fromNbr >> faceIDs >> nbrdVfs;
                if (withAlphaf)
                {
                    fromNbr >> nbrAf;
                }
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
                if (withAlphaf)
                {
                    // the upwind side's face fraction, no sign change
                    wetFraction_.boundaryFieldRef()[patchi][facei] = nbrAf[i];
                    afSet_.set(procPatch.start() + facei);
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


void Foam::pointIsoAdvection::checkIfOnProcPatch(const label facei)
{
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID(facei);

        if (isA<processorPolyPatch>(pbm[patchi]) && !pbm[patchi].empty())
        {
            const label patchFacei = pbm[patchi].whichFace(facei);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}


void Foam::pointIsoAdvection::applyBruteForceBounding()
{
    addProfilingInFunction(geometricVoF);
    bool alpha1Changed = false;

    const scalar snapAlphaTol = dict_.getOrDefault<scalar>("snapTol", 0);
    if (snapAlphaTol > 0)
    {
        alpha1_ =
            alpha1_
           *pos0(alpha1_ - snapAlphaTol)
           *neg0(alpha1_ - (1.0 - snapAlphaTol))
          + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    if (dict_.getOrDefault("clip", true))
    {
        alpha1_.clamp_range(zero_one{});
        alpha1Changed = true;
    }

    if (alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}


void Foam::pointIsoAdvection::writeSurfaceCells() const
{
    if (!mesh_.time().writeTime()) return;

    if (dict_.getOrDefault("writeSurfCells", false))
    {
        cellSet cSet
        (
            IOobject
            (
                "surfCells",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            )
        );

        cSet.insert(surfCells_);

        cSet.write();
    }
}


void Foam::pointIsoAdvection::writeIsoFaces
(
    const UList<List<point>>& faces
) const
{
    if (!writeIsoFacesToFile_ || !mesh_.time().writeTime()) return;

    // Writing isofaces to obj file for inspection, e.g. in paraview
    const fileName outputFile
    (
        mesh_.time().globalPath()
      / "isoFaces"
      / word::printf("isoFaces_%012d.obj", mesh_.time().timeIndex())
    );

    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<List<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = faces;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            mkDir(outputFile.path());
            OBJstream os(outputFile);
            Info<< nl << "pointIsoAdvection: writing iso faces to file: "
                << os.name() << nl << endl;

            for (const auto& procFacePts : allProcFaces)
            {
                for (const auto& facePts : procFacePts)
                {
                    os.writeFace(facePts, false);
                }
            }
        }
    }
    else
    {
        mkDir(outputFile.path());
        OBJstream os(outputFile);
        Info<< nl << "pointIsoAdvection: writing iso faces to file: "
            << os.name() << nl << endl;

        for (const auto& facePts : faces)
        {
            os.writeFace(facePts, false);
        }
    }
}


// ************************************************************************* //
