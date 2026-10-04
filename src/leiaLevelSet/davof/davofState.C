/*---------------------------------------------------------------------------*\
    leia :: davof :: davofState (see davofState.H)
\*---------------------------------------------------------------------------*/

#include "davofState.H"
#include "davofSimplexGeometry.H"
#include "davofSphereGeometry.H"
#include "davofQuadraticFaceGeometry.H"
#include "levelSetPlaneReconstruction.H"
#include "syncTools.H"
#include "boundBox.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofState::davofState(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    alphaSource_
    (
        canonicalSource(dict.getOrDefault<word>("alphaSource", "detrixheAslam"))
    ),
    wispTol_(dict.getOrDefault<scalar>("wispTol", 1e-3)),
    alpha_
    (
        IOobject
        (
            "alpha.davof",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaf_(mesh.nFaces(), Zero),
    mTri_(mesh.nCells(), Zero),
    momentTri_(mesh.nCells(), Zero),
    hasMTri_(false),
    alphafOut_
    (
        IOobject
        (
            "alphaf.davof",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    m_
    (
        IOobject
        (
            "m.davof",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimArea, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    xS_
    (
        IOobject
        (
            "xS.davof",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    AS_
    (
        IOobject
        (
            "AS.davof",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimArea, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    mTet_(mesh.nCells(), Zero),
    maxAlphaClip_(0),
    pPlane_
    (
        IOobject("pPlane.davof", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    xPlane_
    (
        IOobject("xPlane.davof", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedVector(dimLength, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    APlane_
    (
        IOobject("APlane.davof", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimArea, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaPlane_
    (
        IOobject("alphaPlane.davof", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    maxVolDiffPlane_(0),
    maxAreaDiffPlane_(0),
    nFaceFallback_(0)
{
    if
    (
        alphaSource_ != "detrixheAslam"
     && alphaSource_ != "quadraticFaces"
     && alphaSource_ != "planePhaseIndicator"
     && alphaSource_ != "exactSphere"
    )
    {
        FatalIOErrorInFunction(dict_)
            << "Unknown davof alphaSource '" << alphaSource_
            << "'. Valid: detrixheAslam (alias linearInterpolant), "
            << "quadraticFaces, planePhaseIndicator, exactSphere."
            << exit(FatalIOError);
    }
}


Foam::word Foam::davofState::canonicalSource(const word& name)
{
    return (name == "linearInterpolant") ? word("detrixheAslam") : name;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::davofState::resetBulkCell(const label cellI, const scalar alphaValue)
{
    alpha_[cellI] = alphaValue;
    mTet_[cellI] = Zero;
    AS_[cellI] = 0;
    xS_[cellI] = mesh_.cellCentres()[cellI];
}


template<class PhiCell, class PhiFace, class PhiPoint>
void Foam::davofState::cellFromValues
(
    const label cellI,
    PhiCell&& phiCell,
    PhiFace&& phiFace,
    PhiPoint&& phiPoint
)
{
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const labelList& cFaces = mesh_.cells()[cellI];

    const point& xc = C[cellI];
    const scalar phic = phiCell(cellI);

    // One-signed cells carry no interface: alpha is 0 or 1, no tets.
    scalar phiMin = phic, phiMax = phic;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const scalar phif = phiFace(cellI, faceI);
        phiMin = min(phiMin, phif);
        phiMax = max(phiMax, phif);
        const face& f = faces[faceI];
        forAll(f, ip)
        {
            const scalar phip = phiPoint(cellI, f[ip]);
            phiMin = min(phiMin, phip);
            phiMax = max(phiMax, phip);
        }
    }
    if (phiMax <= 0) { resetBulkCell(cellI, 1); return; }
    if (phiMin >  0) { resetBulkCell(cellI, 0); return; }

    scalar liquidVol = 0, totalVol = 0, as = 0;
    vector mt(Zero), xs(Zero);

    // The Detrixhe-Aslam decomposition: one tet per face edge, (cell centre,
    // face centre, edge end points); the same loop order and volume formula
    // as detrixheAslamPhaseIndicator, so the plane source reproduces it bit
    // for bit.
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        const scalar phif = phiFace(cellI, faceI);

        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            const point& p0 = points[f[ip]];
            const point& p1 = points[ip1];
            const scalar phi0 = phiPoint(cellI, f[ip]);
            const scalar phi1 = phiPoint(cellI, ip1);

            const scalar vol = mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
            if (vol <= VSMALL) continue;

            liquidVol += davof::tetNegativeFraction(phic, phif, phi0, phi1)*vol;
            totalVol  += vol;

            const point x[4] = {xc, xf, p0, p1};
            const scalar ph[4] = {phic, phif, phi0, phi1};
            vector av;
            point cen;
            scalar a;
            if (davof::tetZeroSet(x, ph, av, cen, a) > 0)
            {
                mt += av;
                as += a;
                xs += a*cen;
            }
        }
    }

    const scalar raw =
        (totalVol > VSMALL) ? liquidVol/totalVol : (phic < 0 ? 1 : 0);
    const scalar clamped = min(max(raw, scalar(0)), scalar(1));
    maxAlphaClip_ = max(maxAlphaClip_, mag(raw - clamped));

    alpha_[cellI] = clamped;
    mTet_[cellI] = mt;
    AS_[cellI] = as;
    xS_[cellI] = (as > VSMALL) ? xs/as : xc;
}


template<class PhiFace, class PhiPoint>
Foam::scalar Foam::davofState::faceFraction
(
    const label faceI,
    PhiFace&& phiFace,
    PhiPoint&& phiPoint
) const
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];
    const point& xf = mesh_.faceCentres()[faceI];
    const scalar phif = phiFace(faceI);

    scalar liqA = 0, totA = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point& p0 = points[f[ip]];
        const point& p1 = points[ip1];
        const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
        if (triA <= VSMALL) continue;
        liqA +=
            davof::triNegativeFraction(phif, phiPoint(f[ip]), phiPoint(ip1))*triA;
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (phif < 0 ? scalar(1) : scalar(0));
}


template<class PhiFace, class PhiPoint>
void Foam::davofState::faceWettedTriangles
(
    const label faceI,
    PhiFace&& phiFace,
    PhiPoint&& phiPoint,
    vector& S,
    scalar& M
) const
{
    const pointField& points = mesh_.points();
    const face& f = mesh_.faces()[faceI];
    const point& xf = mesh_.faceCentres()[faceI];
    const scalar phif = phiFace(faceI);

    S = Zero;
    M = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point& p0 = points[f[ip]];
        const point& p1 = points[ip1];
        const vector St = 0.5*((p0 - xf) ^ (p1 - xf));
        if (mag(St) <= VSMALL) continue;
        const scalar phi0 = phiPoint(f[ip]);
        const scalar phi1 = phiPoint(ip1);
        const vector Swet = davof::triNegativeFraction(phif, phi0, phi1)*St;
        vector Spoly;
        point xwet;
        davof::triNegativePolygon(xf, p0, p1, phif, phi0, phi1, Spoly, xwet);
        S += Swet;
        M += xwet & Swet;
    }
}


void Foam::davofState::fillAlphafOut()
{
    scalarField& in = alphafOut_.primitiveFieldRef();
    forAll(in, faceI)
    {
        in[faceI] = alphaf_[faceI];
    }
    surfaceScalarField::Boundary& bf = alphafOut_.boundaryFieldRef();
    forAll(bf, patchI)
    {
        const label start = mesh_.boundaryMesh()[patchI].start();
        fvsPatchScalarField& pf = bf[patchI];
        forAll(pf, i)
        {
            pf[i] = alphaf_[start + i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::davofState::computeFromImplicitSurface(const implicitSurface& surface)
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const pointField& P = mesh_.points();

    scalarField psiC(C.size());
    forAll(C, i) psiC[i] = surface.value(C[i]);
    scalarField psiF(Cf.size());
    forAll(Cf, i) psiF[i] = surface.value(Cf[i]);
    scalarField psiP(P.size());
    forAll(P, i) psiP[i] = surface.value(P[i]);

    maxAlphaClip_ = 0;
    forAll(C, cellI)
    {
        cellFromValues
        (
            cellI,
            [&](const label c) { return psiC[c]; },
            [&](const label, const label f) { return psiF[f]; },
            [&](const label, const label p) { return psiP[p]; }
        );
    }
    forAll(alphaf_, faceI)
    {
        alphaf_[faceI] = faceFraction
        (
            faceI,
            [&](const label f) { return psiF[f]; },
            [&](const label p) { return psiP[p]; }
        );
    }
    // The area normal from the face triangulation, summed on the fly:
    // the same fan triangles as faceFraction, their wetted area vectors
    // into the owner (outward) and the neighbour (inward).
    {
        const labelUList& own = mesh_.faceOwner();
        const labelUList& nei = mesh_.faceNeighbour();
        const label nInt = mesh_.nInternalFaces();
        mTri_ = Zero;
        momentTri_ = Zero;
        forAll(alphaf_, faceI)
        {
            vector wet;
            scalar M;
            faceWettedTriangles
            (
                faceI,
                [&](const label f) { return psiF[f]; },
                [&](const label p) { return psiP[p]; },
                wet, M
            );
            mTri_[own[faceI]] -= wet;
            momentTri_[own[faceI]] += M - (C[own[faceI]] & wet);
            if (faceI < nInt)
            {
                mTri_[nei[faceI]] += wet;
                momentTri_[nei[faceI]] -= M - (C[nei[faceI]] & wet);
            }
        }
        hasMTri_ = true;
    }
    alpha_.correctBoundaryConditions();
    fillAlphafOut();
}


void Foam::davofState::computeFromQuadraticFaces(const implicitSurface& surface)
{
    // Cell fractions, tet reference and the Detrixhe-Aslam face fractions
    // (the fallback) first.
    computeFromImplicitSurface(surface);

    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // mTri_ holds the Detrixhe-Aslam triangle sum of
    // computeFromImplicitSurface; every face whose quadratic cut succeeds
    // replaces its DA triangle contribution by the quadratic one (both
    // from the same fan about x_f, nothing stored per face).
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const vectorField& C = mesh_.cellCentres();
    scalarField psiPq(P.size());
    forAll(P, i) psiPq[i] = surface.value(P[i]);

    label nFallback = 0;
    forAll(alphaf_, faceI)
    {
        const scalar A = mag(Sf[faceI]);
        if (A <= VSMALL) continue;
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        const scalar phif = surface.value(xf);

        scalar liqA = 0, totA = 0;
        vector wet = Zero;
        scalar Mq = 0;
        bool fb = false;
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            const point& p0 = P[f[ip]];
            const point& p1 = P[ip1];
            const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
            if (triA <= VSMALL) continue;
            // vertices xf, p0, p1; midpoints of (xf,p0), (p0,p1), (p1,xf)
            const scalar ph[6] =
            {
                phif,
                surface.value(p0),
                surface.value(p1),
                surface.value(0.5*(xf + p0)),
                surface.value(0.5*(p0 + p1)),
                surface.value(0.5*(p1 + xf))
            };
            bool tfb = false;
            const scalar fr = davof::triQuadraticNegativeFraction(xf, p0, p1, ph, tfb);
            if (tfb) { fb = true; break; }
            liqA += fr*triA;
            totA += triA;
            const vector St = 0.5*((p0 - xf) ^ (p1 - xf));
            wet += fr*St;
            // The wet centroid of the triangle from its linear cut (exact
            // for a plane; immaterial on a planar face).
            vector Sl;
            point xwet;
            davof::triNegativePolygon(xf, p0, p1, ph[0], ph[1], ph[2], Sl, xwet);
            Mq += xwet & (fr*St);
        }
        if (fb) { ++nFallback; continue; }
        if (totA > VSMALL)
        {
            alphaf_[faceI] = min(max(liqA/totA, scalar(0)), scalar(1));
            vector wetDA;
            scalar MDA;
            faceWettedTriangles
            (
                faceI,
                [&](const label fi) { return surface.value(Cf[fi]); },
                [&](const label p) { return psiPq[p]; },
                wetDA, MDA
            );
            mTri_[own[faceI]] += wetDA - wet;
            momentTri_[own[faceI]] +=
                (Mq - (C[own[faceI]] & wet)) - (MDA - (C[own[faceI]] & wetDA));
            if (faceI < nInt)
            {
                mTri_[nei[faceI]] -= wetDA - wet;
                momentTri_[nei[faceI]] -=
                    (Mq - (C[nei[faceI]] & wet)) - (MDA - (C[nei[faceI]] & wetDA));
            }
        }
    }
    nFaceFallback_ = returnReduce(nFallback, sumOp<label>());

    fillAlphafOut();
}


void Foam::davofState::computeFromExactSphere(const implicitSphere& sphere)
{
    // The linear-interpolant state first: it provides the tet zero-set
    // reference (mTet, AS, xS) and the alpha/alphaf fallbacks.
    computeFromImplicitSurface(sphere);

    const point c = sphere.center();
    const scalar R = sphere.radius();
    const pointField& P = mesh_.points();
    const faceList& faces = mesh_.faces();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& V = mesh_.V();

    // Faces: the planar polygon inside the ball. A face with every vertex
    // inside the ball is inside (convex faces), fraction exactly 1; a polygon
    // disjoint from the circle sums its sectors to a round-off zero. Both
    // are snapped (1e-12), so the bulk carries exactly equal face fractions
    // and areaNormal() zeroes it instead of counting round-off wisps.
    const scalar snapTol = 1e-12;
    forAll(alphaf_, faceI)
    {
        const scalar A = mag(Sf[faceI]);
        if (A <= VSMALL) { alphaf_[faceI] = 0; continue; }
        const face& f = faces[faceI];
        bool allInside = true;
        forAll(f, i)
        {
            if (magSqr(P[f[i]] - c) > R*R) { allInside = false; break; }
        }
        if (allInside) { alphaf_[faceI] = 1; continue; }
        scalar a = davof::faceBallArea(f, P, Cf[faceI], Sf[faceI]/A, c, R)/A;
        if (a < snapTol) a = 0;
        else if (a > 1 - snapTol) a = 1;
        alphaf_[faceI] = min(max(a, scalar(0)), scalar(1));
    }

    // Cells: axis-aligned boxes only.
    const labelListList& cellPoints = mesh_.cellPoints();
    forAll(alpha_, cellI)
    {
        const pointField cp(P, cellPoints[cellI]);
        const boundBox bb(cp, false);
        const vector span = bb.span();
        const scalar Vbb = span.x()*span.y()*span.z();
        if (mag(Vbb - V[cellI]) > 1e-8*V[cellI])
        {
            FatalErrorInFunction
                << "davof alphaSource exactSphere needs axis-aligned hexahedral "
                << "cells; cell " << cellI << " has volume " << V[cellI]
                << " but bounding-box volume " << Vbb << exit(FatalError);
        }
        // Quick decisions: the box entirely inside or entirely outside.
        scalar farthest = 0, nearest2 = 0;
        for (label i = 0; i < 8; ++i)
        {
            const point corner
            (
                (i & 1) ? bb.max().x() : bb.min().x(),
                (i & 2) ? bb.max().y() : bb.min().y(),
                (i & 4) ? bb.max().z() : bb.min().z()
            );
            farthest = max(farthest, mag(corner - c));
        }
        for (direction k = 0; k < 3; ++k)
        {
            const scalar q = min(max(c[k], bb.min()[k]), bb.max()[k]) - c[k];
            nearest2 += q*q;
        }
        if (farthest <= R)
        {
            alpha_[cellI] = 1;
        }
        else if (nearest2 >= R*R)
        {
            alpha_[cellI] = 0;
        }
        else
        {
            const scalar a = davof::boxBallVolume(bb.min(), bb.max(), c, R)/V[cellI];
            alpha_[cellI] = min(max(a, scalar(0)), scalar(1));
        }
    }
    maxAlphaClip_ = 0;
    // The exact fractions ARE the state: on these planar faces the scalar
    // identity -sum_f alpha_f S_f is the exact wetted area vector, so the
    // Detrixhe-Aslam fan sums of the first pass must not stand in for it
    // (they did from 2026-10-02 to 2026-10-04: the exact state's normal was
    // then first order; the sphere ladder of 2026-09-29 predates the fan).
    hasMTri_ = false;

    alpha_.correctBoundaryConditions();
    fillAlphafOut();
}


void Foam::davofState::computeFromPlanePhaseIndicator
(
    const volScalarField& psi,
    const implicitSurface* analytic
)
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const pointField& P = mesh_.points();
    const label nCells = mesh_.nCells();

    // Per-cell plane n_c . x + d_c (planeSet), bulk sign otherwise: the same
    // reconstruction as detrixheAslamPhaseIndicator / computeFaceAreaFractions.
    vectorField nCell(nCells, Zero);
    scalarField dCell(nCells, Zero);
    labelList planeSet(nCells, 0);
    scalarField signVal(nCells);
    forAll(signVal, c)
    {
        signVal[c] = analytic ? analytic->value(C[c]) : psi[c];
    }

    // Every cell (see the header): the same plane as detrixheAslamPhaseIndicator
    // builds when its band marks every cell, as it does inside leiaSetFields.
    const coupledFaceNeighbours coupledNei(mesh_, psi);
    forAll(planeSet, cellI)
    {
        vector nc(Zero);
        scalar d = 0;
        if (analytic)
        {
            nc = analytic->grad(C[cellI]);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                nc /= nmag;
                d = analytic->value(C[cellI]) - (nc & C[cellI]);
            }
        }
        else
        {
            const scalarList pc =
                leastSquaresPlaneCoeffs(mesh_, psi, cellI, coupledNei);
            nc = vector(pc[0], pc[1], pc[2]);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                d = pc[3]/nmag;
            }
        }
        const scalar nmag = mag(nc);
        if (nmag < SMALL) continue;
        nCell[cellI] = nc/nmag;
        dCell[cellI] = d;
        planeSet[cellI] = 1;
    }

    // Cells.
    maxAlphaClip_ = 0;
    forAll(alpha_, cellI)
    {
        if (!planeSet[cellI])
        {
            resetBulkCell(cellI, (signVal[cellI] < 0) ? 1 : 0);
            continue;
        }
        const vector n = nCell[cellI];
        const scalar d = dCell[cellI];
        cellFromValues
        (
            cellI,
            [&](const label c) { return (n & C[c]) + d; },
            [&](const label, const label f) { return (n & Cf[f]) + d; },
            [&](const label, const label p) { return (n & P[p]) + d; }
        );
    }

    // Faces: the fraction of face f under the plane of cell c (bulk sign when
    // c has no plane), central average over owner and neighbour.
    auto planeFrac = [&](const vector& n, const scalar d, const label faceI)
    {
        return faceFraction
        (
            faceI,
            [&](const label f) { return (n & Cf[f]) + d; },
            [&](const label p) { return (n & P[p]) + d; }
        );
    };
    auto cellFrac = [&](const label c, const label faceI) -> scalar
    {
        return planeSet[c]
          ? planeFrac(nCell[c], dCell[c], faceI)
          : ((signVal[c] < 0) ? scalar(1) : scalar(0));
    };

    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    for (label faceI = 0; faceI < nInt; ++faceI)
    {
        alphaf_[faceI] =
            0.5*(cellFrac(own[faceI], faceI) + cellFrac(nei[faceI], faceI));
    }

    // Coupled patches: the neighbour plane arrives by one matched exchange of
    // the owner-side plane data; every rank evaluates it on its own copy of
    // the face geometry. Non-coupled boundary faces: owner side only.
    const label nBnd = mesh_.nBoundaryFaces();
    vectorField nB(nBnd, Zero);
    scalarField dB(nBnd, Zero);
    labelList setB(nBnd, 0);
    scalarField sB(nBnd, Zero);
    for (label bf = 0; bf < nBnd; ++bf)
    {
        const label c = own[nInt + bf];
        nB[bf] = nCell[c];
        dB[bf] = dCell[c];
        setB[bf] = planeSet[c];
        sB[bf] = signVal[c];
    }
    syncTools::swapBoundaryFaceList(mesh_, nB);
    syncTools::swapBoundaryFaceList(mesh_, dB);
    syncTools::swapBoundaryFaceList(mesh_, setB);
    syncTools::swapBoundaryFaceList(mesh_, sB);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const bool coupled = pp.coupled();
        forAll(pp, i)
        {
            const label faceI = pp.start() + i;
            const scalar aOwn = cellFrac(own[faceI], faceI);
            if (coupled)
            {
                const label bf = faceI - nInt;
                const scalar aNei =
                    setB[bf]
                  ? planeFrac(nB[bf], dB[bf], faceI)
                  : ((sB[bf] < 0) ? scalar(1) : scalar(0));
                alphaf_[faceI] = 0.5*(aOwn + aNei);
            }
            else
            {
                alphaf_[faceI] = aOwn;
            }
        }
    }

    alpha_.correctBoundaryConditions();
    fillAlphafOut();
}


void Foam::davofState::areaNormal()
{
    vectorField& m = m_.primitiveFieldRef();
    m = Zero;
    const vectorField& Sf = mesh_.faceAreas();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    forAll(Sf, faceI)
    {
        // S_f is outward for the owner and inward for the neighbour.
        const vector contrib = alphaf_[faceI]*Sf[faceI];
        m[own[faceI]] -= contrib;
        if (faceI < nInt)
        {
            m[nei[faceI]] += contrib;
        }
    }
    // The face-triangulation sum, when the state was initialised from a
    // surface, is the area normal: identical on planar faces, exact where
    // the scalar identity is not (warped faces).
    if (hasMTri_)
    {
        m = mTri_;
    }
    // A cell whose face fractions are all equal holds no interface (all wet,
    // all dry): its sum_f S_f vanishes only to round-off, which would leave a
    // ~1e-19 m^2 "wisp" in every bulk cell. Set m_c = 0 exactly there.
    const cellList& cells = mesh_.cells();
    forAll(cells, c)
    {
        const labelList& cf = cells[c];
        scalar lo = GREAT, hi = -GREAT;
        forAll(cf, i)
        {
            lo = min(lo, alphaf_[cf[i]]);
            hi = max(hi, alphaf_[cf[i]]);
        }
        if (lo == hi && !(regenerated_.size() && regenerated_[c]))
        {
            m[c] = Zero;
            mTri_[c] = Zero;
        }
    }
    m_.correctBoundaryConditions();
}


void Foam::davofState::planePosition()
{
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();

    // sum_f alpha_f (x_f - x_c).S_f^out per cell: S_f is outward for the
    // owner and inward for the neighbour, as in areaNormal(). When the
    // state came from a surface, the moment of the fan triangles
    // (momentTri_) replaces it: identical on planar faces, exact on warped
    // ones, where (x - x_c).n_f is not constant over the face.
    scalarField faceMoment(mesh_.nCells(), Zero);
    if (hasMTri_)
    {
        faceMoment = momentTri_;
    }
    else
    {
        forAll(Sf, faceI)
        {
            const scalar af = alphaf_[faceI];
            if (af == 0) continue;
            faceMoment[own[faceI]] += af*((Cf[faceI] - C[own[faceI]]) & Sf[faceI]);
            if (faceI < nInt)
            {
                faceMoment[nei[faceI]] -= af*((Cf[faceI] - C[nei[faceI]]) & Sf[faceI]);
            }
        }
    }

    scalarField& p = pPlane_.primitiveFieldRef();
    vectorField& xp = xPlane_.primitiveFieldRef();
    scalarField& Ap = APlane_.primitiveFieldRef();
    scalarField& ap = alphaPlane_.primitiveFieldRef();
    maxVolDiffPlane_ = 0;
    maxAreaDiffPlane_ = 0;

    forAll(cells, c)
    {
        p[c] = 0;
        xp[c] = C[c];
        Ap[c] = 0;
        ap[c] = alpha_[c];
        if (!isInterfaceCell(c)) continue;

        const scalar mMag = mag(m_[c]);
        const vector n = m_[c]/mMag;
        p[c] = (3.0*alpha_[c]*V[c] - faceMoment[c])/mMag;

        // The plane through the cell on the Detrixhe-Aslam tets: the plane is
        // linear, so the tet fractions and zero sets are exact for it.
        const point& xc = C[c];
        const scalar phic = -p[c];
        auto phiAt = [&](const point& x) -> scalar
        {
            return ((x - xc) & n) - p[c];
        };
        scalar liquidVol = 0, totalVol = 0, as = 0;
        vector mt(Zero), xs(Zero);
        const labelList& cFaces = cells[c];
        forAll(cFaces, cf)
        {
            const label faceI = cFaces[cf];
            const face& f = faces[faceI];
            const point& xf = Cf[faceI];
            const scalar phif = phiAt(xf);
            forAll(f, ip)
            {
                const label ip1 = f.nextLabel(ip);
                const point& p0 = points[f[ip]];
                const point& p1 = points[ip1];
                const scalar vol = mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
                if (vol <= VSMALL) continue;
                const scalar phi0 = phiAt(p0);
                const scalar phi1 = phiAt(p1);
                liquidVol += davof::tetNegativeFraction(phic, phif, phi0, phi1)*vol;
                totalVol += vol;
                const point x[4] = {xc, xf, p0, p1};
                const scalar ph[4] = {phic, phif, phi0, phi1};
                vector av;
                point cen;
                scalar a;
                if (davof::tetZeroSet(x, ph, av, cen, a) > 0)
                {
                    mt += av;
                    as += a;
                    xs += a*cen;
                }
            }
        }
        ap[c] = (totalVol > VSMALL) ? liquidVol/totalVol : alpha_[c];
        Ap[c] = as;
        xp[c] = (as > VSMALL) ? xs/as : (xc + p[c]*n);
        maxVolDiffPlane_ = max(maxVolDiffPlane_, mag(ap[c] - alpha_[c]));
        maxAreaDiffPlane_ =
            max(maxAreaDiffPlane_, mag(mt - m_[c])/Foam::pow(V[c], 2.0/3.0));
    }
    maxVolDiffPlane_ = returnReduce(maxVolDiffPlane_, maxOp<scalar>());
    maxAreaDiffPlane_ = returnReduce(maxAreaDiffPlane_, maxOp<scalar>());

    pPlane_.correctBoundaryConditions();
    xPlane_.correctBoundaryConditions();
    APlane_.correctBoundaryConditions();
    alphaPlane_.correctBoundaryConditions();
}


void Foam::davofState::setRegenerated
(
    const boolList& isI,
    const scalarField& alphaf,
    const vectorField& mTri,
    const scalarField& momentTri
)
{
    regenerated_ = isI;
    alphaf_ = alphaf;
    mTri_ = mTri;
    momentTri_ = momentTri;
    hasMTri_ = true;
    fillAlphafOut();
    areaNormal();
    planePosition();
}


Foam::scalar Foam::davofState::maxConsistency() const
{
    const scalarField& V = mesh_.V();
    scalar worst = 0;
    forAll(mTet_, c)
    {
        worst = max(worst, mag(m_[c] - mTet_[c])/Foam::pow(V[c], 2.0/3.0));
    }
    return returnReduce(worst, maxOp<scalar>());
}


Foam::scalar Foam::davofState::wispThreshold(const label cellI) const
{
    return wispTol_*Foam::pow(mesh_.V()[cellI], 2.0/3.0);
}


bool Foam::davofState::isInterfaceCell(const label cellI) const
{
    return mag(m_[cellI]) > wispThreshold(cellI);
}


bool Foam::davofState::isWisp(const label cellI) const
{
    const scalar a = mag(m_[cellI]);
    return (a > 0) && (a <= wispThreshold(cellI));
}


void Foam::davofState::write()
{
    alpha_.write();
    alphafOut_.write();
    m_.write();
    xS_.write();
    AS_.write();
    pPlane_.write();
    xPlane_.write();
    APlane_.write();
    alphaPlane_.write();
}


// ************************************************************************* //
