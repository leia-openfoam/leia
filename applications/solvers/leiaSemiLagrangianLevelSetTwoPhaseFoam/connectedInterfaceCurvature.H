/*---------------------------------------------------------------------------*\
    Connected-interface curvature for two-dimensional level-set calculations.

    This is deliberately not another cell-local curvature extension.  It first
    turns one globally continuous, point-interpolated zero contour into a
    topological polyline:

      * compute each zero crossing on a shared mesh face exactly once;
      * give that same node to both adjacent cut-cell elements;
      * audit the resulting component count, closure and fit coverage.

    Curvature is then a property of an interface ELEMENT, not a cell.  A fixed
    arclength-neighbourhood on the connected chain is fitted in its local
    tangent/normal frame.  The element curvatures are regularised only along
    the chain with a one-dimensional Helmholtz filter.  On every closed
    component the additive null mode is fixed by the planar Gauss--Bonnet
    constraint integral(kappa ds) = +/-2*pi.

    Finally the element value is copied along the element's normal ray directly
    to the registered face field consumed by the CSF model.  A cell field is
    filled by the same rule for diagnostics only; the production force never
    interpolates that cell field back to faces.

    Serial/2-D first.  Parallel needs reconstruction-node exchange across
    processor patches; 3-D needs connected surface patches and a regularised
    shape-operator fit.
\*---------------------------------------------------------------------------*/

#ifndef connectedInterfaceCurvature_H
#define connectedInterfaceCurvature_H

#include "fvCFD.H"
#include "HashSet.H"
#include "levelSetPlaneReconstruction.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "footPointCurvature.H"
#include "levelSetImplicitSurfaces.H"

namespace Foam
{

//- Small dense SPD solve used by the rdfQuadratic estimator (the 5x5 normal
//  equations of the quadratic RDF fit).  In-place Cholesky; b overwritten with
//  the solution; returns false when the matrix is not positive definite
//  (degenerate stencil) so the caller can fall back.
inline bool rdfQuadraticCholeskySolve
(
    scalar* A,
    scalar* b,
    const label n
)
{
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= SMALL) { return false; }
        d = Foam::sqrt(d);
        A[j*n + j] = d;
        for (label i = j + 1; i < n; ++i)
        {
            scalar s = A[i*n + j];
            for (label k = 0; k < j; ++k) { s -= A[i*n + k]*A[j*n + k]; }
            A[i*n + j] = s/d;
        }
    }
    for (label i = 0; i < n; ++i)
    {
        scalar s = b[i];
        for (label k = 0; k < i; ++k) { s -= A[i*n + k]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label k = i + 1; k < n; ++k) { s -= A[k*n + i]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    return true;
}

//- Element-level output used by manufactured geometry gates.  The production
//  solver normally leaves this null; exposing the exact reconstructed elements
//  avoids testing a second, duplicated implementation in post-processing.
struct connectedInterfaceCurvatureDiagnostics
{
    label nElements{0};
    label nComponents{0};
    label nOpenEndpoints{0};
    label nFits{0};
    label nFallbacks{0};
    vectorField elementCentres;
    scalarField elementLengths;
    scalarField rawCurvature;
    scalarField regularizedCurvature;
    boolList fitOk;

    void clear()
    {
        nElements = 0;
        nComponents = 0;
        nOpenEndpoints = 0;
        nFits = 0;
        nFallbacks = 0;
        elementCentres.clear();
        elementLengths.clear();
        rawCurvature.clear();
        regularizedCurvature.clear();
        fitOk.clear();
    }
};

inline label computeConnectedInterfaceCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand,
    volScalarField& kappa,
    surfaceScalarField& kappaFace,
    connectedInterfaceCurvatureDiagnostics* diagnostics = nullptr
)
{
    if (diagnostics) { diagnostics->clear(); }
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "connectedInterface curvature is serial-only: shared interface "
            << "nodes are not yet exchanged across processor patches."
            << abort(FatalError);
    }

    const dictionary& extDict =
        static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("curvatureExtension");
    const label halfWidth = max
    (
        label(2), extDict.getOrDefault<label>("fitHalfWidth", 3)
    );
    const scalar lambda = max
    (
        scalar(0), extDict.getOrDefault<scalar>("tangentialRegularization", 1.0)
    );
    const label nRegularizationIter = max
    (
        label(0), extDict.getOrDefault<label>("regularizationIterations", 30)
    );
    const word regularizationOperator =
        extDict.getOrDefault<word>("regularizationOperator", "helmholtz");
    const label preserveModesThrough = max
    (
        label(1), extDict.getOrDefault<label>("preserveModesThrough", 4)
    );
    const word curvatureEstimator =
        extDict.getOrDefault<word>("estimator", "connectedFit");
    const bool analyticCurvatureOracle =
        curvatureEstimator == "analyticImplicitSurface";
    // rdfQuadratic: STATELESS per-step signed-distance function (RDF) to the
    // connected polyline -- a pure function of the current zero set, never
    // transported -- fitted by per-cell CPC quadratics; curvature from the
    // fits; face-native blended foot-point delivery.  Parameter-free: the
    // tangential regularizer is bypassed; only the Gauss--Bonnet additive-mode
    // constraint (an exact invariant, not a parameter) is applied.
    const bool rdfQuadraticEstimator = curvatureEstimator == "rdfQuadratic";
    const word geometrySource =
        extDict.getOrDefault<word>("geometrySource", "levelSetField");
    const bool analyticGeometryOracle =
        geometrySource == "analyticImplicitSurface";
    if
    (
        curvatureEstimator != "connectedFit"
     && !analyticCurvatureOracle
     && !rdfQuadraticEstimator
    )
    {
        FatalErrorInFunction
            << "Unknown connected-interface curvature estimator '"
            << curvatureEstimator << "'. Valid: connectedFit, "
            << "rdfQuadratic, analyticImplicitSurface." << abort(FatalError);
    }
    if (geometrySource != "levelSetField" && !analyticGeometryOracle)
    {
        FatalErrorInFunction
            << "Unknown connected-interface geometry source '"
            << geometrySource << "'. Valid: levelSetField, "
            << "analyticImplicitSurface." << abort(FatalError);
    }

    autoPtr<implicitSurface> analyticSurface;
    if (analyticCurvatureOracle || analyticGeometryOracle)
    {
        const dictionary& surfaceDict =
            static_cast<const fvSolution&>(mesh).subDict("levelSet")
            .subDict("implicitSurface");
        analyticSurface = implicitSurface::New
        (
            surfaceDict.get<word>("type"),
            surfaceDict
        );
    }

    vector ez(Zero);
    const Vector<label> gd = mesh.geometricD();
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (gd[cmpt] == -1) { ez[cmpt] = 1; }
    }
    if (mag(ez) < SMALL)
    {
        FatalErrorInFunction
            << "connectedInterface curvature currently requires a 2-D mesh "
            << "with one empty direction." << abort(FatalError);
    }

    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = mesh.faceCentres();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField hCell(mesh.nCells(), GREAT);
    forAll(own, facei)
    {
        const scalar d = mag(C[own[facei]] - C[nei[facei]]);
        hCell[own[facei]] = min(hCell[own[facei]], d);
        hCell[nei[facei]] = min(hCell[nei[facei]], d);
    }
    const scalar hGlobal = 1.0/gMax(mesh.deltaCoeffs().primitiveField());
    forAll(hCell, celli)
    {
        if (hCell[celli] == GREAT) { hCell[celli] = hGlobal; }
    }

    const coupledFaceNeighbours coupledNei(mesh, psi);

    DynamicList<label> elemCell;
    DynamicList<vector> elemNormal;
    DynamicList<vector> elemP0;
    DynamicList<vector> elemP1;
    DynamicList<label> elemFace0;
    DynamicList<label> elemFace1;
    label nBandCells = 0;
    label nPlaneOk = 0;
    label nTwoCrossings = 0;
    forAll(narrowBand, celli) { nBandCells += narrowBand[celli] == 1; }

    // A single continuous piecewise-linear zero contour.  Point interpolation
    // gives one value per mesh vertex; each internal-face crossing is therefore
    // computed ONCE and supplied to both adjacent cells.  The two crossings of
    // a cut cell are the endpoints of its interface element, so neighbouring
    // elements share an identical node before any curvature fit is attempted.
    scalarField psiAtPoint(points.size(), 0);
    if (analyticGeometryOracle)
    {
        forAll(points, pointi)
        {
            psiAtPoint[pointi] = analyticSurface->value(points[pointi]);
        }
    }
    else
    {
        const pointScalarField psiPoint
        (
            volPointInterpolation::New(mesh).interpolate(psi)
        );
        psiAtPoint = psiPoint.primitiveField();
    }
    vectorField crossing0(mesh.nCells(), Zero);
    vectorField crossing1(mesh.nCells(), Zero);
    labelList crossingFace0(mesh.nCells(), -1);
    labelList crossingFace1(mesh.nCells(), -1);
    labelList nCrossing(mesh.nCells(), 0);

    auto addCrossing = [&]
    (
        const label celli,
        const label facei,
        const point& crossing
    )
    {
        if (nCrossing[celli] == 0)
        {
            crossing0[celli] = crossing;
            crossingFace0[celli] = facei;
        }
        else if (nCrossing[celli] == 1)
        {
            crossing1[celli] = crossing;
            crossingFace1[celli] = facei;
        }
        ++nCrossing[celli];
    };

    forAll(own, facei)
    {
        const face& f = faces[facei];
        point q0(Zero), q1(Zero);
        label nHit = 0;
        forAll(f, ip)
        {
            const label p0 = f[ip];
            const label p1 = f[f.fcIndex(ip)];
            const scalar s0 = psiAtPoint[p0];
            const scalar s1 = psiAtPoint[p1];
            if (s0*s1 < 0)
            {
                const scalar t = s0/(s0 - s1);
                const point q = points[p0] + t*(points[p1] - points[p0]);
                if (nHit == 0) { q0 = q; }
                else if (nHit == 1) { q1 = q; }
                ++nHit;
            }
        }
        if (nHit != 2) { continue; }
        point crossing = 0.5*(q0 + q1);
        crossing -= ((crossing - Cf[facei]) & ez)*ez;
        addCrossing(own[facei], facei, crossing);
        addCrossing(nei[facei], facei, crossing);
    }

    forAll(nCrossing, celli)
    {
        if (nCrossing[celli] != 2) { continue; }
        ++nTwoCrossings;
        vector n(Zero);
        if (analyticGeometryOracle)
        {
            n = analyticSurface->grad
            (
                0.5*(crossing0[celli] + crossing1[celli])
            );
        }
        else
        {
            const scalarList pc =
                leastSquaresPlaneCoeffs(mesh, psi, celli, coupledNei);
            n = vector(pc[0], pc[1], pc[2]);
        }
        n -= (n & ez)*ez;
        const scalar nmag = mag(n);
        if (nmag < SMALL) { continue; }
        ++nPlaneOk;

        elemCell.append(celli);
        elemNormal.append(n/nmag);
        elemP0.append(crossing0[celli]);
        elemP1.append(crossing1[celli]);
        elemFace0.append(crossingFace0[celli]);
        elemFace1.append(crossingFace1[celli]);
    }

    const label nElem = elemCell.size();
    if (!nElem)
    {
        kappaFace = dimensionedScalar("zero", kappaFace.dimensions(), 0);
        Info<< "connectedInterface curvature: no elements (band="
            << nBandCells << ", plane=" << nPlaneOk
            << ", twoCrossings=" << nTwoCrossings << ")" << endl;
        return 0;
    }

    // Endpoint pairing.  The global mesh-face label is the topological node:
    // two independently reconstructed segments crossing one shared face are
    // made conforming by replacing their crossings with the common midpoint.
    const label nEnd = 2*nElem;
    labelList epElem(nEnd), epSide(nEnd), epFace(nEnd), epMate(nEnd, -1);
    vectorField epPoint(nEnd);
    for (label ei = 0; ei < nElem; ++ei)
    {
        epElem[2*ei] = epElem[2*ei + 1] = ei;
        epSide[2*ei] = 0;
        epSide[2*ei + 1] = 1;
        epFace[2*ei] = elemFace0[ei];
        epFace[2*ei + 1] = elemFace1[ei];
        epPoint[2*ei] = elemP0[ei];
        epPoint[2*ei + 1] = elemP1[ei];
    }

    // Assemble chains with a minimum-distance spanning rule.  Pairing every
    // shared-face crossing independently can close short loops before the whole
    // interface is connected when neighbouring cell planes choose different
    // exit faces.  Treat each segment as an initial component, repeatedly join
    // the closest two FREE endpoints belonging to different components, and
    // only then close each resulting chain.  Shared-face candidates receive no
    // penalty; a nonconforming cross-face join receives a one-cell penalty.
    labelList parent(nElem);
    forAll(parent, ei) { parent[ei] = ei; }
    auto root = [&](label ei) -> label
    {
        while (parent[ei] != ei) { ei = parent[ei]; }
        return ei;
    };

    for (label join = 0; join < nElem - 1; ++join)
    {
        label aBest = -1, bBest = -1;
        scalar scoreBest = GREAT;
        for (label a = 0; a < nEnd; ++a)
        {
            if (epMate[a] >= 0) { continue; }
            for (label b = a + 1; b < nEnd; ++b)
            {
                if
                (
                    epMate[b] >= 0 || epElem[a] == epElem[b]
                 || root(epElem[a]) == root(epElem[b])
                )
                {
                    continue;
                }
                const label ca = elemCell[epElem[a]];
                const label cb = elemCell[epElem[b]];
                const scalar h = min(hCell[ca], hCell[cb]);
                const scalar d2 = magSqr(epPoint[a] - epPoint[b]);
                const scalar facePenalty = epFace[a] == epFace[b] ? 0 : 1;
                const scalar score = d2/sqr(h + VSMALL) + facePenalty;
                if (score < scoreBest)
                {
                    scoreBest = score;
                    aBest = a;
                    bBest = b;
                }
            }
        }
        if (aBest < 0) { break; }
        epMate[aBest] = bBest;
        epMate[bBest] = aBest;
        const label ra = root(epElem[aBest]);
        const label rb = root(epElem[bBest]);
        parent[rb] = ra;
    }

    // Every assembled component is now one open chain (two free endpoints).
    // Close it only after all possible inter-segment joins have been made.
    for (label a = 0; a < nEnd; ++a)
    {
        if (epMate[a] >= 0) { continue; }
        label bBest = -1;
        scalar dBest = GREAT;
        for (label b = a + 1; b < nEnd; ++b)
        {
            if
            (
                epMate[b] >= 0 || epElem[a] == epElem[b]
             || root(epElem[a]) != root(epElem[b])
            )
            {
                continue;
            }
            const scalar d2 = magSqr(epPoint[a] - epPoint[b]);
            if (d2 < dBest) { dBest = d2; bBest = b; }
        }
        if (bBest >= 0)
        {
            epMate[a] = bBest;
            epMate[bBest] = a;
        }
    }

    // Make every paired crossing a single geometrical node.
    for (label a = 0; a < nEnd; ++a)
    {
        const label b = epMate[a];
        if (b > a)
        {
            vector p = 0.5*(epPoint[a] + epPoint[b]);
            epPoint[a] = p;
            epPoint[b] = p;
        }
    }

    vectorField elemCentre(nElem);
    scalarField elemLength(nElem, 0);
    for (label ei = 0; ei < nElem; ++ei)
    {
        elemP0[ei] = epPoint[2*ei];
        elemP1[ei] = epPoint[2*ei + 1];
        elemCentre[ei] = 0.5*(elemP0[ei] + elemP1[ei]);
        elemLength[ei] = mag(elemP1[ei] - elemP0[ei]);
    }

    // Component labels and closure audit of the explicit segment chain.
    labelList component(nElem, -1);
    DynamicList<bool> componentClosed;
    label nComponents = 0;
    for (label seed = 0; seed < nElem; ++seed)
    {
        if (component[seed] >= 0) { continue; }
        labelList queue(nElem, -1);
        label qHead = 0, qTail = 0;
        queue[qTail++] = seed;
        component[seed] = nComponents;
        bool closed = true;
        while (qHead < qTail)
        {
            const label ei = queue[qHead++];
            for (label side = 0; side < 2; ++side)
            {
                const label mate = epMate[2*ei + side];
                if (mate < 0)
                {
                    closed = false;
                    continue;
                }
                const label ej = epElem[mate];
                if (component[ej] < 0)
                {
                    component[ej] = nComponents;
                    queue[qTail++] = ej;
                }
            }
        }
        componentClosed.append(closed);
        ++nComponents;
    }

    // ---- rdfQuadratic machinery ------------------------------------------
    // A per-time-step signed-distance function d to the connected zero set,
    // computed at cell centres in a band around the interface, then fitted by
    // per-cell quadratics over the compact cell-point-cell stencil.  d is a
    // pure function of the current zero set: curvature evaluated from these
    // fits never reads the transported psi PROFILE (|grad d| ~ 1 by
    // construction, unconditionally well scaled), which severs the measured
    // profile-corruption -> curvature feedback.  Temporal continuity: the
    // geometry moves continuously under sub-cell interface motion, so d, the
    // fit coefficients, the feet and the delivered kappa are all continuous
    // (no crossing-set switching).  The fits also serve the face-native
    // delivery at the bottom of this function.
    //
    // The distance is measured to a C1 REPRESENTATION of the zero set, not to
    // the raw chords: per-element arclength-window parabola fits (the same
    // chain stencil, frame and fitParabola as the connectedFit estimator),
    // overlapping so their pointwise minimum is continuous.  The raw polyline
    // is only C0 -- its chord sagitta error h^2 kappa/8 oscillates at
    // wavelength ~h, exactly the cell-centre sampling spacing, and aliases
    // into the fitted Hessian as O(1) RELATIVE curvature noise (measured:
    // +-20% on an exact circle, catastrophic on transported zero sets).  The
    // parabola pieces carry O(h^3) geometry error and remove that noise
    // source at its origin.  This smooths GEOMETRY, never curvature: the
    // estimator stays parameter-free (no tangential filter; fitHalfWidth is a
    // stencil size, and Gauss--Bonnet below fixes the additive mode exactly).
    label rdfAx0 = -1, rdfAx1 = -1;
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (gd[cmpt] == 1)
        {
            if (rdfAx0 < 0) { rdfAx0 = cmpt; } else { rdfAx1 = cmpt; }
        }
    }
    auto rdfInPlane = [&](const vector& v) -> vector
    {
        return v - (v & ez)*ez;
    };
    scalarField rdfD;
    boolList rdfHasD;
    boolList rdfFitValid;
    scalarField rdfCoeffs;
    label rdfGeomFits = -1;

    if (rdfQuadraticEstimator)
    {
        // C1 interface geometry: per-element parabola fits of the chain
        // element centres in the local (tangent, normal) frame -- identical
        // stencil walk, frame construction and fitParabola call to the
        // connectedFit estimator.  Elements whose chain window is deficient
        // keep chord distances (their pieces are simply absent from the min).
        boolList elemGeomValid(nElem, false);
        vectorField elemThat(nElem, Zero);
        vectorField elemNbarF(nElem, Zero);
        scalarField elemFitA(nElem, 0);
        scalarField elemFitB(nElem, 0);
        scalarField elemFitC(nElem, 0);
        scalarField elemFitH(nElem, 0);
        rdfGeomFits = 0;
        for (label ei = 0; ei < nElem; ++ei)
        {
            DynamicList<label> stencil(2*halfWidth + 1);
            labelHashSet used;
            stencil.append(ei);
            used.insert(ei);

            label sideCount[2] = {0, 0};
            for (label initialSide = 0; initialSide < 2; ++initialSide)
            {
                label current = ei;
                label exitSide = initialSide;
                for (label step = 0; step < halfWidth; ++step)
                {
                    const label mate = epMate[2*current + exitSide];
                    if (mate < 0) { break; }
                    const label next = epElem[mate];
                    if (used.found(next)) { break; }
                    stencil.append(next);
                    used.insert(next);
                    ++sideCount[initialSide];
                    current = next;
                    exitSide = 1 - epSide[mate];
                }
            }
            if (stencil.size() < 5 || sideCount[0] < 2 || sideCount[1] < 2)
            {
                continue;
            }

            vector nbar(Zero);
            DynamicList<vector> fitPoints(stencil.size());
            for (const label ej : stencil)
            {
                vector nj = elemNormal[ej];
                if ((nj & elemNormal[ei]) < 0) { nj = -nj; }
                nbar += elemLength[ej]*nj;
                fitPoints.append(elemCentre[ej]);
            }
            nbar = rdfInPlane(nbar);
            if (mag(nbar) < SMALL) { continue; }
            nbar /= mag(nbar);
            vector that = ez ^ nbar;
            that /= mag(that) + VSMALL;

            const scalar hE = hCell[elemCell[ei]];
            scalar A = 0, B = 0, Cc = 0;
            if
            (
                !fitParabola
                (
                    fitPoints,
                    elemCentre[ei],
                    that,
                    nbar,
                    hE,
                    A,
                    B,
                    &Cc
                )
            )
            {
                continue;
            }
            elemThat[ei] = that;
            elemNbarF[ei] = nbar;
            elemFitA[ei] = A;
            elemFitB[ei] = B;
            elemFitC[ei] = Cc;
            elemFitH[ei] = hE;
            elemGeomValid[ei] = true;
            ++rdfGeomFits;
        }

        // Squared in-plane distance from q to element ei's parabola piece
        // zeta(tau) = h (A tau^2/2 + B tau + C), tau = t/h, restricted to the
        // element window |t| <= max(elemLength, 0.25 h) EXTENDED by one
        // element length so adjacent pieces overlap: coarse bracketing over
        // the window, then damped-Newton polish on the foot equation
        // f(t) = (t - tq) + (y(t) - zq) y'(t) = 0 with monotone-decrease
        // acceptance (robust where 1 + y'^2 + (y - zq) y'' is small).
        auto rdfCurveDist2 = [&](const label ei, const vector& q) -> scalar
        {
            const vector d = q - rdfInPlane(elemCentre[ei]);
            const scalar tq = d & elemThat[ei];
            const scalar zq = d & elemNbarF[ei];
            const scalar hE = elemFitH[ei];
            const scalar A = elemFitA[ei];
            const scalar B = elemFitB[ei];
            const scalar Cc = elemFitC[ei];
            const scalar tMax = max(elemLength[ei], 0.25*hE);

            auto d2Of = [&](const scalar t) -> scalar
            {
                const scalar u = t/hE;
                const scalar dy = hE*(0.5*A*u*u + B*u + Cc) - zq;
                const scalar dt = t - tq;
                return dt*dt + dy*dy;
            };

            scalar tBest = min(max(tq, -tMax), tMax);
            scalar d2Best = d2Of(tBest);
            for (label k = 0; k <= 4; ++k)
            {
                const scalar t = -tMax + 0.5*k*tMax;
                const scalar d2 = d2Of(t);
                if (d2 < d2Best) { d2Best = d2; tBest = t; }
            }
            for (label it = 0; it < 3; ++it)
            {
                const scalar u = tBest/hE;
                const scalar y = hE*(0.5*A*u*u + B*u + Cc);
                const scalar yp = A*u + B;
                const scalar ypp = A/hE;
                const scalar f = (tBest - tq) + (y - zq)*yp;
                const scalar fp = 1.0 + yp*yp + (y - zq)*ypp;
                if (mag(fp) < SMALL) { break; }
                const scalar tNew =
                    min(max(tBest - f/fp, -tMax), tMax);
                const scalar d2New = d2Of(tNew);
                if (d2New <= d2Best) { d2Best = d2New; tBest = tNew; }
                else { break; }
            }
            return d2Best;
        };

        // Distance region: cut cells grown 5 face-sweeps; fit region: grown 3
        // (every fit stencil then lies inside the distance region).
        boolList dRegion(mesh.nCells(), false);
        forAll(elemCell, ei) { dRegion[elemCell[ei]] = true; }
        boolList fitRegion;
        for (label sweep = 0; sweep < 5; ++sweep)
        {
            if (sweep == 3) { fitRegion = dRegion; }
            boolList grow(dRegion);
            forAll(own, facei)
            {
                if (dRegion[own[facei]] || dRegion[nei[facei]])
                {
                    grow[own[facei]] = true;
                    grow[nei[facei]] = true;
                }
            }
            dRegion = grow;
        }

        // Signed minimal distance to the C1 geometry; the SIGN comes from
        // psi (robust: only the zero set matters, not the profile).  Chord
        // distances bracket the search: the exact parabola-piece distance is
        // evaluated only for elements whose chord distance lies within one
        // cell size of the chord minimum (the chord and its piece differ by
        // at most the sagitta, far below h).
        rdfD.setSize(mesh.nCells(), GREAT);
        rdfHasD.setSize(mesh.nCells(), false);
        scalarField chordD2(nElem);
        forAll(dRegion, celli)
        {
            if (!dRegion[celli]) { continue; }
            const vector p = rdfInPlane(C[celli]);
            scalar dmin2 = GREAT;
            for (label ei = 0; ei < nElem; ++ei)
            {
                const vector a = rdfInPlane(elemP0[ei]);
                const vector ab = rdfInPlane(elemP1[ei]) - a;
                const scalar ab2 = magSqr(ab);
                scalar t = (ab2 > VSMALL) ? (((p - a) & ab)/ab2) : scalar(0);
                t = min(max(t, scalar(0)), scalar(1));
                chordD2[ei] = magSqr(p - (a + t*ab));
                dmin2 = min(dmin2, chordD2[ei]);
            }
            const scalar reach = Foam::sqrt(dmin2) + hCell[celli];
            const scalar reach2 = reach*reach;
            scalar dmin2C1 = GREAT;
            for (label ei = 0; ei < nElem; ++ei)
            {
                if (chordD2[ei] > reach2) { continue; }
                const scalar d2 =
                    elemGeomValid[ei]
                  ? rdfCurveDist2(ei, p)
                  : chordD2[ei];
                dmin2C1 = min(dmin2C1, d2);
            }
            rdfD[celli] =
                (psi[celli] >= 0 ? 1.0 : -1.0)*Foam::sqrt(dmin2C1);
            rdfHasD[celli] = true;
        }

        // IDW-weighted quadratic fit of d over the cell-point-cell stencil,
        // constant term pinned to the cell value (basis dx,dy,dx2,dxdy,dy2).
        rdfFitValid.setSize(mesh.nCells(), false);
        rdfCoeffs.setSize(5*mesh.nCells(), 0);
        const labelListList& cPoints = mesh.cellPoints();
        const labelListList& pCells = mesh.pointCells();
        forAll(fitRegion, celli)
        {
            if (!fitRegion[celli] || !rdfHasD[celli]) { continue; }
            labelHashSet sten;
            for (const label pointi : cPoints[celli])
            {
                for (const label cj : pCells[pointi])
                {
                    if (cj != celli && rdfHasD[cj]) { sten.insert(cj); }
                }
            }
            if (sten.size() < 5) { continue; }

            scalar A[25] = {0};
            scalar bvec[5] = {0};
            const vector pc = rdfInPlane(C[celli]);
            const scalar dc = rdfD[celli];
            for (const label cj : sten)
            {
                const vector dvec = rdfInPlane(C[cj]) - pc;
                const scalar dx = dvec[rdfAx0];
                const scalar dy = dvec[rdfAx1];
                const scalar w = 1.0/max(mag(dvec), SMALL);
                const scalar w2 = w*w;
                const scalar bb[5] = {dx, dy, dx*dx, dx*dy, dy*dy};
                const scalar ds = rdfD[cj] - dc;
                for (label k = 0; k < 5; ++k)
                {
                    bvec[k] += w2*bb[k]*ds;
                    for (label l = 0; l < 5; ++l)
                    {
                        A[5*k + l] += w2*bb[k]*bb[l];
                    }
                }
            }
            if (!rdfQuadraticCholeskySolve(A, bvec, 5)) { continue; }
            for (label k = 0; k < 5; ++k)
            {
                rdfCoeffs[5*celli + k] = bvec[k];
            }
            rdfFitValid[celli] = true;
        }
    }

    // Value, gradient and curvature of cell celli's RDF quadratic at x.
    auto rdfEval = [&]
    (
        const label celli,
        const vector& x,
        scalar& value,
        vector& g,
        scalar& kap
    ) -> bool
    {
        if (!rdfFitValid[celli]) { return false; }
        const scalar* cf = &rdfCoeffs[5*celli];
        const vector dvec = rdfInPlane(x - C[celli]);
        const scalar dx = dvec[rdfAx0];
        const scalar dy = dvec[rdfAx1];
        value =
            rdfD[celli]
          + cf[0]*dx + cf[1]*dy
          + cf[2]*dx*dx + cf[3]*dx*dy + cf[4]*dy*dy;
        const scalar gx = cf[0] + 2.0*cf[2]*dx + cf[3]*dy;
        const scalar gy = cf[1] + cf[3]*dx + 2.0*cf[4]*dy;
        const scalar hxx = 2.0*cf[2], hxy = cf[3], hyy = 2.0*cf[4];
        const scalar g2 = gx*gx + gy*gy;
        if (g2 < SMALL*SMALL) { return false; }
        const scalar gm = Foam::sqrt(g2);
        const scalar gHg =
            gx*(hxx*gx + hxy*gy) + gy*(hxy*gx + hyy*gy);
        kap = ((hxx + hyy)*g2 - gHg)/(g2*gm);
        g = Zero;
        g[rdfAx0] = gx;
        g[rdfAx1] = gy;
        return true;
    };

    // Curvature of cell celli's RDF quadratic at the FOOT of x on its zero
    // set (two gradient-projection steps; the RDF is near-linear so this
    // converges immediately).  Trust region: foot within 2.5 cell sizes.
    auto rdfFootKappa = [&]
    (
        const label celli,
        const point& x,
        scalar& kap
    ) -> bool
    {
        scalar val = 0;
        vector g(Zero);
        scalar kTmp = 0;
        vector xx = x;
        for (label it = 0; it < 2; ++it)
        {
            if (!rdfEval(celli, xx, val, g, kTmp)) { return false; }
            const scalar g2 = magSqr(g);
            if (g2 < SMALL*SMALL) { return false; }
            xx -= (val/g2)*g;
        }
        if (mag(rdfInPlane(xx - C[celli])) > 2.5*hCell[celli])
        {
            return false;
        }
        scalar valFoot = 0;
        vector gFoot(Zero);
        return rdfEval(celli, xx, valFoot, gFoot, kap);
    };

    scalarField rawKappa(nElem, 0);
    boolList fitOk(nElem, false);
    label nFallback = 0;

    if (analyticCurvatureOracle)
    {
        // Variable-curvature velocity oracle.  Geometry, topology, element
        // ownership and normal-ray face delivery are identical to the
        // candidate path; only the scalar value attached to each element is
        // replaced by the exact curvature of the configured analytic surface.
        // Do not regularise or Gauss--Bonnet-shift these values below.
        for (label ei = 0; ei < nElem; ++ei)
        {
            rawKappa[ei] = analyticSurface->curvature(elemCentre[ei]);
            fitOk[ei] = true;
        }
    }
    else if (rdfQuadraticEstimator)
    {
        // Element curvature from the RDF quadratic of the element's own cell,
        // evaluated at the element centre (already on the zero set, so the
        // foot projection is a near no-op there).  No tangential filter; the
        // Gauss--Bonnet constraint below fixes the additive mode exactly.
        for (label ei = 0; ei < nElem; ++ei)
        {
            scalar kap = 0;
            if (rdfFootKappa(elemCell[ei], elemCentre[ei], kap))
            {
                rawKappa[ei] = kap;
                fitOk[ei] = true;
            }
            else
            {
                ++nFallback;
                rawKappa[ei] = kappa[elemCell[ei]];
            }
        }
    }
    else for (label ei = 0; ei < nElem; ++ei)
    {
        DynamicList<label> stencil(2*halfWidth + 1);
        labelHashSet used;
        stencil.append(ei);
        used.insert(ei);

        label sideCount[2] = {0, 0};
        for (label initialSide = 0; initialSide < 2; ++initialSide)
        {
            label current = ei;
            label exitSide = initialSide;
            for (label step = 0; step < halfWidth; ++step)
            {
                const label mate = epMate[2*current + exitSide];
                if (mate < 0) { break; }
                const label next = epElem[mate];
                if (used.found(next)) { break; }
                stencil.append(next);
                used.insert(next);
                ++sideCount[initialSide];
                current = next;
                exitSide = 1 - epSide[mate];
            }
        }

        if (stencil.size() < 5 || sideCount[0] < 2 || sideCount[1] < 2)
        {
            ++nFallback;
            rawKappa[ei] = kappa[elemCell[ei]];
            continue;
        }

        vector nbar(Zero);
        DynamicList<vector> fitPoints(stencil.size());
        for (const label ej : stencil)
        {
            vector nj = elemNormal[ej];
            if ((nj & elemNormal[ei]) < 0) { nj = -nj; }
            nbar += elemLength[ej]*nj;
            fitPoints.append(elemCentre[ej]);
        }
        if (mag(nbar) < SMALL)
        {
            ++nFallback;
            rawKappa[ei] = kappa[elemCell[ei]];
            continue;
        }
        nbar /= mag(nbar);
        vector that = ez ^ nbar;
        that /= mag(that) + VSMALL;

        scalar A = 0, B = 0;
        if
        (
            !fitParabola
            (
                fitPoints,
                elemCentre[ei],
                that,
                nbar,
                hCell[elemCell[ei]],
                A,
                B
            )
        )
        {
            ++nFallback;
            rawKappa[ei] = kappa[elemCell[ei]];
            continue;
        }

        rawKappa[ei] =
            -(A/hCell[elemCell[ei]])/Foam::pow(1.0 + B*B, 1.5);
        fitOk[ei] = true;
        elemNormal[ei] = nbar;
    }

    // Tangential regularisation on the connected one-dimensional interface
    // graph.  "helmholtz" solves (I + lambda L_s)k=raw (strong but attenuates
    // physical low modes).  "biharmonic" repeatedly applies
    // (I-theta L_s^2), theta<=1/16: grid modes are damped strongly while the
    // low Fourier modes required by an oscillating droplet are nearly retained.
    scalarField kappaElem(rawKappa);
    if (analyticCurvatureOracle || rdfQuadraticEstimator)
    {
        // Exact analytic values are already the oracle; leave them unchanged.
        // rdfQuadratic is deliberately PARAMETER-FREE: no tangential filter;
        // the multi-cell RDF fit is itself the smoothing, and Gauss--Bonnet
        // below fixes the additive mode.
    }
    else if (regularizationOperator == "biharmonic")
    {
        const scalar theta = min(lambda, scalar(0.0625));
        for (label iter = 0; iter < nRegularizationIter && theta > 0; ++iter)
        {
            scalarField lap(nElem, 0);
            for (label ei = 0; ei < nElem; ++ei)
            {
                const label m0 = epMate[2*ei];
                const label m1 = epMate[2*ei + 1];
                if (m0 < 0 || m1 < 0 || !fitOk[ei]) { continue; }
                const label e0 = epElem[m0];
                const label e1 = epElem[m1];
                if (!fitOk[e0] || !fitOk[e1]) { continue; }
                lap[ei] = 2*kappaElem[ei] - kappaElem[e0] - kappaElem[e1];
            }
            scalarField next(kappaElem);
            for (label ei = 0; ei < nElem; ++ei)
            {
                const label m0 = epMate[2*ei];
                const label m1 = epMate[2*ei + 1];
                if (m0 < 0 || m1 < 0 || !fitOk[ei]) { continue; }
                const label e0 = epElem[m0];
                const label e1 = epElem[m1];
                if (!fitOk[e0] || !fitOk[e1]) { continue; }
                const scalar lap2 = 2*lap[ei] - lap[e0] - lap[e1];
                next[ei] = kappaElem[ei] - theta*lap2;
            }
            kappaElem = next;
        }
    }
    else if
    (
        regularizationOperator == "helmholtz"
     || regularizationOperator == "helmholtzPreserveModes"
    )
    {
        for (label iter = 0; iter < nRegularizationIter && lambda > 0; ++iter)
        {
            scalarField next(kappaElem);
            for (label ei = 0; ei < nElem; ++ei)
            {
                const label m0 = epMate[2*ei];
                const label m1 = epMate[2*ei + 1];
                if (m0 < 0 || m1 < 0 || !fitOk[ei]) { continue; }
                const label e0 = epElem[m0];
                const label e1 = epElem[m1];
                if (!fitOk[e0] || !fitOk[e1]) { continue; }
                next[ei] =
                    (rawKappa[ei] + lambda*(kappaElem[e0] + kappaElem[e1]))
                   /(1.0 + 2.0*lambda);
            }
            kappaElem = next;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown tangential regularizationOperator '"
            << regularizationOperator << "'. Valid: helmholtz, "
            << "helmholtzPreserveModes, biharmonic."
            << abort(FatalError);
    }

    if
    (
        !analyticCurvatureOracle
     && !rdfQuadraticEstimator
     && regularizationOperator == "helmholtzPreserveModes"
    )
    {
        // A scalar Helmholtz parameter cannot simultaneously remove the replay
        // defect and retain resolved droplet modes.  On each closed chain,
        // project the difference (raw-filtered) onto the arclength Fourier
        // space m=2..preserveModesThrough and add that projection back.  The
        // high graph modes remain strongly damped; the physically relevant
        // low modes are not attenuated by the regularizer.  Mode 0 is fixed by
        // Gauss--Bonnet below and mode 1 is deliberately not restored (it is a
        // translation/centering mode, not a shape oscillation).
        for (label comp = 0; comp < nComponents; ++comp)
        {
            if (!componentClosed[comp]) { continue; }

            label seed = -1;
            label nCompElem = 0;
            forAll(component, ei)
            {
                if (component[ei] == comp)
                {
                    if (seed < 0) { seed = ei; }
                    ++nCompElem;
                }
            }
            if (seed < 0 || nCompElem < 2*preserveModesThrough + 1)
            {
                continue;
            }

            DynamicList<label> ordered(nCompElem);
            label current = seed;
            label exitSide = 0;
            for (label step = 0; step < nCompElem; ++step)
            {
                ordered.append(current);
                const label mate = epMate[2*current + exitSide];
                if (mate < 0) { break; }
                const label next = epElem[mate];
                if (next == seed) { break; }
                current = next;
                exitSide = 1 - epSide[mate];
            }
            if (ordered.size() != nCompElem) { continue; }

            scalar totalLength = 0;
            for (const label ei : ordered) { totalLength += elemLength[ei]; }
            if (totalLength < VSMALL) { continue; }

            const label nBasis = 2*(preserveModesThrough - 1);
            List<scalarField> phi(nBasis, scalarField(nCompElem, 0));
            scalar arclength = 0;
            forAll(ordered, oi)
            {
                const label ei = ordered[oi];
                const scalar s = arclength + 0.5*elemLength[ei];
                const scalar angle = constant::mathematical::twoPi*s/totalLength;
                label basis = 0;
                for (label mode = 2; mode <= preserveModesThrough; ++mode)
                {
                    phi[basis++][oi] = Foam::cos(mode*angle);
                    phi[basis++][oi] = Foam::sin(mode*angle);
                }
                arclength += elemLength[ei];
            }

            List<scalarField> gram(nBasis, scalarField(nBasis, 0));
            scalarField rhs(nBasis, 0);
            forAll(ordered, oi)
            {
                const label ei = ordered[oi];
                const scalar weight = elemLength[ei];
                const scalar defect = rawKappa[ei] - kappaElem[ei];
                for (label p = 0; p < nBasis; ++p)
                {
                    rhs[p] += weight*phi[p][oi]*defect;
                    for (label q = 0; q < nBasis; ++q)
                    {
                        gram[p][q] += weight*phi[p][oi]*phi[q][oi];
                    }
                }
            }

            bool nonsingular = true;
            for (label pivot = 0; pivot < nBasis; ++pivot)
            {
                label pivotRow = pivot;
                scalar pivotMag = mag(gram[pivot][pivot]);
                for (label row = pivot + 1; row < nBasis; ++row)
                {
                    if (mag(gram[row][pivot]) > pivotMag)
                    {
                        pivotMag = mag(gram[row][pivot]);
                        pivotRow = row;
                    }
                }
                if (pivotMag < VSMALL*totalLength)
                {
                    nonsingular = false;
                    break;
                }
                if (pivotRow != pivot)
                {
                    Swap(gram[pivot], gram[pivotRow]);
                    Swap(rhs[pivot], rhs[pivotRow]);
                }
                for (label row = pivot + 1; row < nBasis; ++row)
                {
                    const scalar factor =
                        gram[row][pivot]/gram[pivot][pivot];
                    gram[row][pivot] = 0;
                    for (label col = pivot + 1; col < nBasis; ++col)
                    {
                        gram[row][col] -= factor*gram[pivot][col];
                    }
                    rhs[row] -= factor*rhs[pivot];
                }
            }
            if (!nonsingular) { continue; }

            scalarField coeff(nBasis, 0);
            for (label row = nBasis - 1; row >= 0; --row)
            {
                scalar value = rhs[row];
                for (label col = row + 1; col < nBasis; ++col)
                {
                    value -= gram[row][col]*coeff[col];
                }
                coeff[row] = value/gram[row][row];
            }
            forAll(ordered, oi)
            {
                scalar restored = 0;
                for (label p = 0; p < nBasis; ++p)
                {
                    restored += coeff[p]*phi[p][oi];
                }
                kappaElem[ordered[oi]] += restored;
            }
        }
    }

    // The connected reconstruction supplies a global invariant unavailable to
    // independent cell fits.  Preserve all regularised fluctuations but fix
    // the additive mode so each closed planar component has total turning
    // +/-2*pi (sign taken from the reconstructed outward-normal convention).
    scalarField componentShift(max(nComponents, label(1)), 0);
    if (!analyticCurvatureOracle)
    {
        for (label comp = 0; comp < nComponents; ++comp)
        {
            if (!componentClosed[comp]) { continue; }
            scalar length = 0, integral = 0;
            for (label ei = 0; ei < nElem; ++ei)
            {
                if (component[ei] != comp || !fitOk[ei]) { continue; }
                length += elemLength[ei];
                integral += elemLength[ei]*kappaElem[ei];
            }
            if (length < VSMALL) { continue; }
            const scalar target =
                (integral < 0 ? -1.0 : 1.0)
               *constant::mathematical::twoPi;
            const scalar shift = (target - integral)/length;
            // Stored so the rdfQuadratic face-native delivery below can apply
            // the SAME additive-mode correction to its face values.
            componentShift[comp] = shift;
            for (label ei = 0; ei < nElem; ++ei)
            {
                if (component[ei] == comp && fitOk[ei])
                {
                    kappaElem[ei] += shift;
                }
            }
        }
    }

    // Normal-ray parent lookup shared by cell diagnostics and the face-native
    // CSF field.  The empty-direction displacement is removed explicitly.
    auto nearestElement = [&]
    (
        const point& x,
        const scalar maxDistanceFactor
    ) -> label
    {
        scalar best = GREAT;
        label bestElem = -1;
        for (label ei = 0; ei < nElem; ++ei)
        {
            if (!fitOk[ei]) { continue; }
            vector d = x - elemCentre[ei];
            d -= (d & ez)*ez;
            const scalar d2 = magSqr(d);
            const scalar h = hCell[elemCell[ei]];
            if (d2 > sqr(maxDistanceFactor*h)) { continue; }
            const scalar dn = d & elemNormal[ei];
            const scalar dray2 = max(scalar(0), d2 - dn*dn);
            const scalar score = dray2 + 0.01*dn*dn;
            if (score < best)
            {
                best = score;
                bestElem = ei;
            }
        }
        return bestElem;
    };

    boolList inDelivery(mesh.nCells(), false);
    forAll(narrowBand, celli)
    {
        if (narrowBand[celli] == 1) { inDelivery[celli] = true; }
    }
    for (label sweep = 0; sweep < 2; ++sweep)
    {
        boolList grow(inDelivery);
        forAll(own, facei)
        {
            if (inDelivery[own[facei]] || inDelivery[nei[facei]])
            {
                grow[own[facei]] = true;
                grow[nei[facei]] = true;
            }
        }
        inDelivery = grow;
    }

    scalarField& kCell = kappa.primitiveFieldRef();
    forAll(inDelivery, celli)
    {
        if (!inDelivery[celli]) { continue; }
        const label ei = nearestElement(C[celli], 6.0);
        if (ei >= 0) { kCell[celli] = kappaElem[ei]; }
    }
    kappa.correctBoundaryConditions();

    scalarField allFaceKappa(mesh.nFaces(), 0);
    const labelList& faceOwn = mesh.faceOwner();
    const scalarField& faceW = mesh.weights().primitiveField();
    forAll(allFaceKappa, facei)
    {
        const label ei = nearestElement(Cf[facei], 6.0);
        if (ei < 0) { continue; }
        if (rdfQuadraticEstimator)
        {
            // FACE-NATIVE delivery: curvature is born AT the face -- blended
            // owner/neighbour foot-point evaluations of the RDF quadratics,
            // projected from THIS face centre onto the zero set.  The linear
            // geometric blend (never an owner/neighbour switch) keeps kappa_f
            // continuous under sub-cell interface motion; the per-component
            // Gauss--Bonnet shift keeps the delivered mean identical to the
            // element diagnostics.  Falls back to the nearest-element value
            // where a fit is unavailable.
            const scalar shift = componentShift[component[ei]];
            const label o = faceOwn[facei];
            scalar kO = 0, kN = 0;
            const bool okO = rdfFootKappa(o, Cf[facei], kO);
            bool okN = false;
            if (facei < mesh.nInternalFaces())
            {
                okN = rdfFootKappa(nei[facei], Cf[facei], kN);
            }
            if (okO && okN)
            {
                const scalar w = faceW[facei];
                allFaceKappa[facei] = w*kO + (1.0 - w)*kN + shift;
            }
            else if (okO) { allFaceKappa[facei] = kO + shift; }
            else if (okN) { allFaceKappa[facei] = kN + shift; }
            else { allFaceKappa[facei] = kappaElem[ei]; }
        }
        else
        {
            allFaceKappa[facei] = kappaElem[ei];
        }
    }
    scalarField& kfInternal = kappaFace.primitiveFieldRef();
    forAll(kfInternal, facei) { kfInternal[facei] = allFaceKappa[facei]; }
    surfaceScalarField::Boundary& kfbf = kappaFace.boundaryFieldRef();
    forAll(kfbf, patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        forAll(kfbf[patchi], i)
        {
            kfbf[patchi][i] = allFaceKappa[patch.start() + i];
        }
    }

    label nOpenEnds = 0;
    forAll(epMate, endpointi) { nOpenEnds += epMate[endpointi] < 0; }

    if (diagnostics)
    {
        diagnostics->nElements = nElem;
        diagnostics->nComponents = nComponents;
        diagnostics->nOpenEndpoints = nOpenEnds;
        diagnostics->nFits = nElem - nFallback;
        diagnostics->nFallbacks = nFallback;
        diagnostics->elementCentres = elemCentre;
        diagnostics->elementLengths = elemLength;
        diagnostics->rawCurvature = rawKappa;
        diagnostics->regularizedCurvature = kappaElem;
        diagnostics->fitOk = fitOk;
    }
    Info<< "connectedInterface curvature: " << nElem << " elements, "
        << nComponents << " component(s), " << nOpenEnds << " open endpoint(s), "
        << nElem - nFallback << " fits, " << nFallback << " fallback; "
        << "estimator=" << curvatureEstimator
        << ", geometrySource=" << geometrySource;
    if (rdfGeomFits >= 0)
    {
        Info<< ", C1 geometry pieces=" << rdfGeomFits << "/" << nElem;
    }
    Info
        << ", fitHalfWidth=" << halfWidth << ", operator="
        << regularizationOperator << ", parameter=" << lambda
        << ", iterations=" << nRegularizationIter << ", kappaElem in ["
        << gMin(kappaElem) << ", " << gMax(kappaElem) << "] 1/m"
        << " (band=" << nBandCells << ", plane=" << nPlaneOk
        << ", twoCrossings=" << nTwoCrossings << ")" << endl;

    return nFallback + nOpenEnds;
}

} // End namespace Foam

#endif

// ************************************************************************* //
