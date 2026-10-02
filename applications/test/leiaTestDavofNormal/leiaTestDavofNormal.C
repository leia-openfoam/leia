/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestDavofNormal

Description
    DAVOF static gate 1: the Gauss-identity interface normal on an implicit
    surface. The cell volume fractions alpha_c and the face liquid-area
    fractions alpha_f of the surface named in fvSolution levelSet.implicitSurface
    are built by libleiaDavof (fvSolution davof { alphaSource ...; wispTol ...; }),
    the area normal

        m_c = - sum_f alpha_f S_f^out

    is recovered per cell, and its direction is compared with the exact surface
    normal at the reconstructed patch centroid xS_c (NOT at the cell centre,
    which would add an O(h) floor):

        e_c = | m_c/|m_c| - n_exact(xS_c) |,   c in S = { |m_c| > wispTol V_c^(2/3) }

        L1 = sum e_c / N_S,  L2 = sqrt(sum e_c^2 / N_S),  Linf = max e_c

    (plain averages over the interface cells; area-weighted variants and the
    Linf over every cell with |m_c| > 0, i.e. including the wisps, are
    reported next to them). The total interface area sum_c |m_c| is compared
    with the exact area (4 pi R^2 for implicitSphere, a Gauss-Legendre
    quadrature of the parametrisation for signedDistanceEllipsoid).

    Curvature of the DAVOF state (fvSolution davof.curvature { models (...);
    <name> { type quadricFit; ... } }, libleiaDavof curvature/): every model is
    fitted to the ring of plane-polygon centroids and area normals, and scored
    against the exact geometry (implicitSphere, implicitPlane,
    signedDistanceEllipsoid) in the convention kappa = kappa_1 + kappa_2 =
    div(n), 2/R on a sphere:
        E_KAPPA_*       at the interface centroids, |kappa_c - kappa_exact(q)|
                        with q the closest surface point of xPlane_c [1/m];
                        E_K_L2 the same for the Gaussian curvature;
        E_KAPPA_CELL_*  the contour-referenced cell field on the force band
                        (parallel-surface forward map of the model at the cell
                        centre's offset) against the exact parallel-surface
                        curvature at the cell centre;
        E_KAPPA_FACE_*  interpolate(kappaCell) inverted at the face with the
                        level-set solver's parallelSurfaceInverse(kappa_f, d_f,
                        K_f) against the exact interface curvature at the foot
                        of the face centre; E_KAPPA_FACE_FOOT_* the models'
                        interface value at that foot (no parallel-surface
                        algebra), the comparison delivery.
    L1 mean, L2 rms, Linf max over the interface cells, the filled band cells
    and the active faces (internal faces and the owner side of coupled
    patches); KAPPA_REF_L2 is the rms of the exact total curvature over the
    interface cells (relative errors in the tables). The first model is the
    headline (its values also in the wide and the models CSV), every model has
    a row in leiaTestDavofCurvature.csv.

    Optional cross-check on the SAME alpha_c (fvSolution davof.crossCheck
    { set geometricVoF; models (plicRDF gradAlpha isoAlpha); <model> {...} }):
    OpenFOAM's libgeometricVoF reconstruction schemes are run on alpha.davof
    and scored with the same norms; their normal_ is the PLIC area vector
    pointing INTO alpha = 1, so it is flipped before the comparison, and their
    centre_ is the evaluation point of the exact normal.

    Outputs (case directory): leiaTestDavofNormal.csv (one row, the DAVOF
    result and the consistency diagnostics) and leiaTestDavofNormalModels.csv
    (tidy, one row per MODEL: davof, plicRDF, gradAlpha, isoAlpha),
    leiaTestDavofCurvature.csv (tidy, one row per curvature model); fields
    alpha.davof, alphaf.davof, m.davof, xS.davof, AS.davof, nDavof,
    eNormal.davof, interfaceMarker.davof (0 bulk, 1 interface, 2 wisp),
    eNormal.<model>, kappa.davof.<curvature model>, K.davof.<...>,
    kappa1/kappa2.davof.<...>, eKappa.davof and kappaCell.davof (headline).

    Parallel-safe (every norm is reduced). Options:
        -alphaName <word>   the leiaSetFields field to compare alpha.davof with
                            (default alpha.water; MAX_ALPHA_DIFF_DA is 0 for
                            alphaSource planePhaseIndicator with geometrySource
                            levelSetField, the same construction)
        -expectExact        exact-solution gate for an implicitPlane: exit
                            non-zero unless the normal error, the tet/face
                            consistency, the alpha difference against the
                            plane cut of src/vof/foamGeometry.H, the position,
                            the realizability diagnostics and (when curvature
                            models are configured) the curvature times h at the
                            centroids and the faces are at round-off.

    Run in a meshed, leiaSetFields-initialised case:
        blockMesh; leiaSetFields -alphaName alpha.water; leiaTestDavofNormal

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cpuTime.H"
#include "davofState.H"
#include "davofCurvature.H"
#include "davofCurvatureDelivery.H"
#include "davofQuadraticFaceGeometry.H"
#include "levelSetImplicitSurfaces.H"
#include "leiaVersionRegistry.H"
#include "reconstructionSchemes.H"
#include "foamGeometry.H"
#include "plicInterfaceSurface.H"
#include <functional>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct normResult
{
    word model;
    label nInterface = 0;
    label nWisp = 0;
    scalar A = 0;
    scalar L1 = 0, L2 = 0, Linf = 0;
    scalar L1aw = 0, L2aw = 0;
    scalar LinfAll = 0;
    // Position: the distance [m] of the reconstructed plane's polygon centroid
    // from the exact surface, mean / rms / max over the interface cells.
    scalar P1 = 0, P2 = 0, Pinf = 0;
    // DAVOF only: the same for the foot point x_c + p_c n_c of the plane.
    scalar Pfoot2 = 0, Pfootinf = 0;
    scalar cpu = 0;
};


//- Position norms of the points xpos over the interface cells (the same set
//  as evaluateNormals: |av| > wispTol V^(2/3)), dist the exact unsigned
//  distance to the surface. Fills r.P1, r.P2, r.Pinf.
void evaluatePositions
(
    normResult& r,
    const fvMesh& mesh,
    const vectorField& av,
    const vectorField& xpos,
    const scalar wispTol,
    const std::function<scalar(const point&)>& dist,
    scalarField* eOut
)
{
    const scalarField& V = mesh.V();
    scalar s1 = 0, s2 = 0, li = 0;
    label nI = 0;
    forAll(av, c)
    {
        const scalar a = mag(av[c]);
        if (a <= wispTol*Foam::pow(V[c], 2.0/3.0)) continue;
        const scalar e = dist(xpos[c]);
        ++nI;
        s1 += e;
        s2 += e*e;
        li = max(li, e);
        if (eOut) (*eOut)[c] = e;
    }
    nI = returnReduce(nI, sumOp<label>());
    s1 = returnReduce(s1, sumOp<scalar>());
    s2 = returnReduce(s2, sumOp<scalar>());
    r.Pinf = returnReduce(li, maxOp<scalar>());
    r.P1 = (nI > 0) ? s1/nI : 0;
    r.P2 = (nI > 0) ? Foam::sqrt(s2/nI) : 0;
}


//- Norms of the normal error of the per-cell area vectors av (sign*av points
//  from the liquid into the gas), scored at the points xc against the exact
//  normal of the surface. eOut (optional) receives e_c on interface cells.
normResult evaluateNormals
(
    const word& model,
    const fvMesh& mesh,
    const vectorField& av,
    const scalar sign,
    const vectorField& xc,
    const implicitSurface& surface,
    const scalar wispTol,
    scalarField* eOut,
    scalarField* markerOut
)
{
    const scalarField& V = mesh.V();
    scalar s1 = 0, s2 = 0, sw1 = 0, sw2 = 0, li = 0, liAll = 0;
    scalar sA = 0, sAI = 0;
    label nI = 0, nW = 0;

    forAll(av, c)
    {
        const scalar a = mag(av[c]);
        if (a <= 0)
        {
            if (markerOut) (*markerOut)[c] = 0;
            continue;
        }
        sA += a;
        const vector n = sign*av[c]/a;
        vector nEx = surface.grad(xc[c]);
        const scalar nExMag = mag(nEx);
        if (nExMag > VSMALL) nEx /= nExMag;
        const scalar e = mag(n - nEx);
        liAll = max(liAll, e);

        if (a > wispTol*Foam::pow(V[c], 2.0/3.0))
        {
            ++nI;
            s1 += e;
            s2 += e*e;
            li = max(li, e);
            sw1 += a*e;
            sw2 += a*e*e;
            sAI += a;
            if (eOut) (*eOut)[c] = e;
            if (markerOut) (*markerOut)[c] = 1;
        }
        else
        {
            ++nW;
            if (markerOut) (*markerOut)[c] = 2;
        }
    }

    normResult r;
    r.model = model;
    r.nInterface = returnReduce(nI, sumOp<label>());
    r.nWisp = returnReduce(nW, sumOp<label>());
    r.A = returnReduce(sA, sumOp<scalar>());
    s1 = returnReduce(s1, sumOp<scalar>());
    s2 = returnReduce(s2, sumOp<scalar>());
    sw1 = returnReduce(sw1, sumOp<scalar>());
    sw2 = returnReduce(sw2, sumOp<scalar>());
    sAI = returnReduce(sAI, sumOp<scalar>());
    r.Linf = returnReduce(li, maxOp<scalar>());
    r.LinfAll = returnReduce(liAll, maxOp<scalar>());
    if (r.nInterface > 0)
    {
        r.L1 = s1/r.nInterface;
        r.L2 = Foam::sqrt(s2/r.nInterface);
    }
    if (sAI > 0)
    {
        r.L1aw = sw1/sAI;
        r.L2aw = Foam::sqrt(sw2/sAI);
    }
    return r;
}


//- Surface area of the axis-aligned ellipsoid with half-axes a, b, c: a
//  composite 12-point Gauss-Legendre rule (32 x 64 panels) over the
//  parametrisation x = a sin t cos p, y = b sin t sin p, z = c cos t,
//  dA = sin t sqrt(b^2 c^2 sin^2 t cos^2 p + a^2 c^2 sin^2 t sin^2 p
//  + a^2 b^2 cos^2 t) dt dp; the integrand is smooth, the rule is exact to
//  round-off for this purpose.
scalar ellipsoidArea(const scalar a, const scalar b, const scalar c)
{
    const List<Pair<scalar>>& gl = Foam::davof::gaussLegendre01(12);
    const label nT = 32, nP = 64;
    const scalar pi = constant::mathematical::pi;
    scalar A = 0;
    for (label it = 0; it < nT; ++it)
    {
        forAll(gl, gt)
        {
            const scalar t = pi*(it + gl[gt].first())/nT;
            const scalar wt = pi/nT*gl[gt].second();
            const scalar st = Foam::sin(t), ct = Foam::cos(t);
            for (label ip = 0; ip < nP; ++ip)
            {
                forAll(gl, gp)
                {
                    const scalar p = 2*pi*(ip + gl[gp].first())/nP;
                    const scalar wp = 2*pi/nP*gl[gp].second();
                    const scalar sp = Foam::sin(p), cp = Foam::cos(p);
                    A += wt*wp*st*Foam::sqrt
                    (
                        sqr(b*c*st*cp) + sqr(a*c*st*sp) + sqr(a*b*ct)
                    );
                }
            }
        }
    }
    return A;
}


//- Exact geometry at a point x: the total curvature kappa_1 + kappa_2 and the
//  Gaussian curvature of the surface at the closest point of x, and the
//  signed distance d of x (positive on the gas side).
typedef std::function<void(const point&, scalar&, scalar&, scalar&)> exactGeomType;


struct curvResult
{
    word model;
    label nInterface = 0, nBand = 0, nUnfilled = 0, nActive = 0, nOneSided = 0;
    label nSkipped = 0;
    label nFallback = 0;
    scalar kappaRef2 = 0;                    // rms of the exact total curvature
    scalar S1 = 0, S2 = 0, Sinf = 0;         // at the interface centroids
    scalar K2 = 0;                           // Gaussian curvature, rms error
    scalar C1 = 0, C2 = 0, Cinf = 0;         // cell field on the force band
    scalar F1 = 0, F2 = 0, Finf = 0;         // faces, interpolate + inverse
    scalar G1 = 0, G2 = 0, Ginf = 0;         // faces, the models at the foot
    scalar cpu = 0;
};


//- Curvature norms of one model (see the Description).
curvResult evaluateCurvature
(
    const word& model,
    const fvMesh& mesh,
    const davofState& st,
    const davofCurvature& cm,
    const davofCurvatureDelivery& dl,
    const exactGeomType& exactGeom,
    scalarField* eOut,
    scalarField* kappaCellOut
)
{
    curvResult r;
    r.model = model;
    const vectorField& xPl = st.xPlane().primitiveField();
    const scalarField& kap = cm.kappa().primitiveField();
    const scalarField& Kg = cm.K().primitiveField();

    scalar s1 = 0, s2 = 0, si = 0, sK2 = 0, ref2 = 0;
    label nI = 0;
    forAll(kap, c)
    {
        if (!st.isInterfaceCell(c)) continue;
        scalar ke, Ke, de;
        exactGeom(xPl[c], ke, Ke, de);
        const scalar e = mag(kap[c] - ke);
        const scalar eK = mag(Kg[c] - Ke);
        ++nI;
        s1 += e;
        s2 += e*e;
        si = max(si, e);
        sK2 += eK*eK;
        ref2 += ke*ke;
        if (eOut) (*eOut)[c] = e;
    }
    nI = returnReduce(nI, sumOp<label>());
    s1 = returnReduce(s1, sumOp<scalar>());
    s2 = returnReduce(s2, sumOp<scalar>());
    sK2 = returnReduce(sK2, sumOp<scalar>());
    ref2 = returnReduce(ref2, sumOp<scalar>());
    r.nInterface = nI;
    r.Sinf = returnReduce(si, maxOp<scalar>());
    if (nI > 0)
    {
        r.S1 = s1/nI;
        r.S2 = Foam::sqrt(s2/nI);
        r.K2 = Foam::sqrt(sK2/nI);
        r.kappaRef2 = Foam::sqrt(ref2/nI);
    }

    const vectorField& C = mesh.cellCentres();
    scalar c1 = 0, c2 = 0, ci = 0;
    label nB = 0;
    forAll(C, c)
    {
        if (!dl.bandCell()[c] || !dl.filledCell()[c]) continue;
        scalar ke, Ke, de;
        exactGeom(C[c], ke, Ke, de);
        const scalar kpar = parallelSurfaceForward(ke, de, Ke);
        const scalar e = mag(dl.kappaCell()[c] - kpar);
        ++nB;
        c1 += e;
        c2 += e*e;
        ci = max(ci, e);
        if (kappaCellOut) (*kappaCellOut)[c] = dl.kappaCell()[c];
    }
    nB = returnReduce(nB, sumOp<label>());
    c1 = returnReduce(c1, sumOp<scalar>());
    c2 = returnReduce(c2, sumOp<scalar>());
    r.Cinf = returnReduce(ci, maxOp<scalar>());
    if (nB > 0)
    {
        r.C1 = c1/nB;
        r.C2 = Foam::sqrt(c2/nB);
    }

    const vectorField& Cf = mesh.faceCentres();
    scalar f1 = 0, f2 = 0, fi = 0, g1 = 0, g2 = 0, gi = 0;
    label nF = 0;
    forAll(Cf, f)
    {
        if (!dl.activeFace()[f]) continue;
        scalar ke, Ke, de;
        exactGeom(Cf[f], ke, Ke, de);
        const scalar e = mag(dl.kappaFaceInv()[f] - ke);
        const scalar g = mag(dl.kappaFaceFoot()[f] - ke);
        ++nF;
        f1 += e; f2 += e*e; fi = max(fi, e);
        g1 += g; g2 += g*g; gi = max(gi, g);
    }
    nF = returnReduce(nF, sumOp<label>());
    f1 = returnReduce(f1, sumOp<scalar>());
    f2 = returnReduce(f2, sumOp<scalar>());
    g1 = returnReduce(g1, sumOp<scalar>());
    g2 = returnReduce(g2, sumOp<scalar>());
    r.Finf = returnReduce(fi, maxOp<scalar>());
    r.Ginf = returnReduce(gi, maxOp<scalar>());
    if (nF > 0)
    {
        r.F1 = f1/nF; r.F2 = Foam::sqrt(f2/nF);
        r.G1 = g1/nF; r.G2 = Foam::sqrt(g2/nF);
    }

    r.nBand = dl.nBandCells();
    r.nUnfilled = dl.nUnfilledCells();
    r.nActive = dl.nActiveFaces();
    r.nOneSided = dl.nOneSidedFaces();
    r.nSkipped = dl.nSkippedFaces();
    r.nFallback = cm.nFallback();
    return r;
}


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "alphaName", "word",
        "Name of the leiaSetFields volume fraction to compare with "
        "(default alpha.water)."
    );
    argList::addBoolOption
    (
        "expectExact",
        "Exact-solution gate (implicitPlane): fail unless normal error, "
        "tet/face consistency and the plane-cut alpha difference are at round-off."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    leia::reportVersions(Info);
    leia::writeVersions(runTime);

    const word alphaName = args.getOrDefault<word>("alphaName", "alpha.water");
    const bool expectExact = args.found("expectExact");

    const fvSolution& fvSol = static_cast<const fvSolution&>(mesh);
    const dictionary& levelSetDict = fvSol.subDict("levelSet");
    const dictionary& surfDict = levelSetDict.subDict("implicitSurface");
    const word surfType = surfDict.get<word>("type");
    autoPtr<implicitSurface> surface = implicitSurface::New(surfType, surfDict);
    const dictionary& davofDict = fvSol.subDict("davof");

    // Exact area of the closed surface, where it is known.
    scalar radius = 0;
    scalar aExact = -1;
    if (surfType == "implicitSphere")
    {
        radius = surfDict.get<scalar>("radius");
        aExact = 4.0*constant::mathematical::pi*radius*radius;
    }
    else if (surfType == "signedDistanceEllipsoid")
    {
        // RADIUS (and R_OVER_H) = the smallest half-axis.
        const vector ax = surfDict.get<vector>("axes");
        radius = min(ax.x(), min(ax.y(), ax.z()));
        aExact = ellipsoidArea(ax.x(), ax.y(), ax.z());
    }

    // Exact unsigned distance of a point from the surface: closed form for
    // the sphere and the plane, |psi|/|grad psi| otherwise (exact for the
    // signed-distance surfaces, first order for an algebraic one).
    std::function<scalar(const point&)> exactDistance;
    if (surfType == "implicitSphere")
    {
        const vector centre = surfDict.get<vector>("center");
        const scalar R = radius;
        exactDistance = [=](const point& x) -> scalar
        {
            return mag(mag(x - centre) - R);
        };
    }
    else if (surfType == "implicitPlane")
    {
        const implicitPlane pl(surfDict);
        const vector n0 = pl.normal()/mag(pl.normal());
        const vector x0 = pl.position();
        exactDistance = [=](const point& x) -> scalar
        {
            return mag((x - x0) & n0);
        };
    }
    else
    {
        const implicitSurface& s = surface();
        exactDistance = [&s](const point& x) -> scalar
        {
            return mag(s.value(x))/max(mag(s.grad(x)), VSMALL);
        };
    }

    // Exact geometry for the curvature norms (kappa_1 + kappa_2 and K at the
    // closest surface point of x, signed distance of x). implicitSphere's own
    // curvature() returns 1/R, half the div(n) convention, so the sphere is
    // written out; signedDistanceEllipsoid::curvature() is the total curvature
    // at the closest point q = x - psi grad(psi) already, and its Gaussian
    // curvature is K = 1/(a^2 b^2 c^2 (X^2/a^4 + Y^2/b^4 + Z^2/c^4)^2).
    exactGeomType exactGeom;
    bool haveExactGeom = true;
    if (surfType == "implicitSphere")
    {
        const vector centre = surfDict.get<vector>("center");
        const scalar R = radius;
        exactGeom = [=](const point& x, scalar& kappa, scalar& K, scalar& d)
        {
            kappa = 2.0/R;
            K = 1.0/(R*R);
            d = mag(x - centre) - R;
        };
    }
    else if (surfType == "implicitPlane")
    {
        const implicitSurface& s = surface();
        exactGeom = [&s](const point& x, scalar& kappa, scalar& K, scalar& d)
        {
            kappa = 0;
            K = 0;
            d = s.value(x)/max(mag(s.grad(x)), VSMALL);
        };
    }
    else if (surfType == "signedDistanceEllipsoid")
    {
        const signedDistanceEllipsoid* ell =
            dynamic_cast<const signedDistanceEllipsoid*>(&surface());
        const vector centre = ell->center();
        const vector ax = ell->axes();
        const implicitSurface& s = surface();
        exactGeom = [&s, centre, ax]
        (
            const point& x, scalar& kappa, scalar& K, scalar& d
        )
        {
            d = s.value(x);
            const vector q = x - d*s.grad(x);
            kappa = s.curvature(x);
            const vector X = q - centre;
            const scalar a2 = sqr(ax.x()), b2 = sqr(ax.y()), c2 = sqr(ax.z());
            const scalar sInv =
                sqr(X.x())/sqr(a2) + sqr(X.y())/sqr(b2) + sqr(X.z())/sqr(c2);
            K = 1.0/(a2*b2*c2*sqr(sInv));
        };
    }
    else
    {
        haveExactGeom = false;
    }

    // ---- The DAVOF state -----------------------------------------------------
    davofState st(mesh, davofDict);
    cpuTime timer;
    autoPtr<volScalarField> psiPtr;
    if (st.alphaSource() == "detrixheAslam")
    {
        st.computeFromImplicitSurface(surface());
    }
    else if (st.alphaSource() == "quadraticFaces")
    {
        st.computeFromQuadraticFaces(surface());
    }
    else if (st.alphaSource() == "exactSphere")
    {
        const implicitSphere* sphere =
            dynamic_cast<const implicitSphere*>(&surface());
        if (!sphere)
        {
            FatalErrorInFunction
                << "davof alphaSource exactSphere needs levelSet.implicitSurface "
                << "of type implicitSphere, found " << surfType
                << exit(FatalError);
        }
        st.computeFromExactSphere(*sphere);
    }
    else
    {
        // The plane source needs the level set leiaSetFields wrote; the plane
        // is fitted in every cell (davofState.H explains why: that is what the
        // initialiser does, and MAX_ALPHA_DIFF_DA = 0 states it).
        psiPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "psi", runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::NO_WRITE
                ),
                mesh
            )
        );
        const word geometrySource =
            levelSetDict.subDict("phaseIndicator").getOrDefault<word>
            (
                "geometrySource", "levelSetField"
            );
        const implicitSurface* analytic =
            (geometrySource == "analyticImplicitSurface") ? &surface() : nullptr;
        st.computeFromPlanePhaseIndicator(psiPtr(), analytic);
    }
    st.areaNormal();
    st.planePosition();
    const scalar cpuDavof = timer.cpuTimeIncrement();

    // ---- Mesh size -----------------------------------------------------------
    const bool is2D = mesh.nGeometricD() < 3;
    const boundBox bb(mesh.points());
    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar Lz = bb.max().z() - bb.min().z();
    const label nCellsGlobal = returnReduce(mesh.nCells(), sumOp<label>());
    const scalar dx =
        (nCellsGlobal > 0)
      ? (is2D
         ? Foam::sqrt(Lx*Ly/nCellsGlobal)
         : Foam::cbrt(Lx*Ly*Lz/nCellsGlobal))
      : 0;

    // ---- DAVOF norms ---------------------------------------------------------
    volScalarField eDavof
    (
        IOobject("eNormal.davof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField marker
    (
        IOobject("interfaceMarker.davof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    volVectorField nDavof
    (
        IOobject("nDavof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedVector(dimless, Zero),
        zeroGradientFvPatchVectorField::typeName
    );

    volScalarField ePosDavof
    (
        IOobject("ePos.davof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    );

    List<normResult> results;
    {
        normResult r = evaluateNormals
        (
            "davof", mesh, st.m().primitiveField(), 1.0,
            st.xS().primitiveField(), surface(), st.wispTol(),
            &eDavof.primitiveFieldRef(), &marker.primitiveFieldRef()
        );
        r.cpu = cpuDavof;
        // The explicit plane position: its polygon centroid, and the foot
        // point x_c + p_c n_c of the cell centre on the plane.
        evaluatePositions
        (
            r, mesh, st.m().primitiveField(), st.xPlane().primitiveField(),
            st.wispTol(), exactDistance, &ePosDavof.primitiveFieldRef()
        );
        vectorField xFoot(mesh.cellCentres());
        forAll(nDavof, c)
        {
            if (st.isInterfaceCell(c))
            {
                nDavof[c] = st.m()[c]/mag(st.m()[c]);
                xFoot[c] += st.pPlane()[c]*nDavof[c];
            }
        }
        normResult rFoot;
        evaluatePositions
        (
            rFoot, mesh, st.m().primitiveField(), xFoot, st.wispTol(),
            exactDistance, nullptr
        );
        r.Pfoot2 = rFoot.P2;
        r.Pfootinf = rFoot.Pinf;
        results.append(r);
    }
    // A COPY: results grows below (the cross-check appends), which reallocates
    // the list and would leave a reference dangling.
    const normResult rd = results[0];

    // The scalar identity -sum_f alpha_f S_f as a diagnostic next to the
    // face-triangulation normal that the state uses when it was initialised
    // from a surface: identical on planar faces, 8 % apart on a plane
    // through vertex-perturbed hexahedra (2026-10-02). Computed here on the
    // fly from alpha_f and S_f; nothing is stored.
    normResult rAlpha;
    if (st.hasMTri())
    {
        vectorField mAlpha(mesh.nCells(), Zero);
        const vectorField& SfAll = mesh.faceAreas();
        const labelUList& ownAll = mesh.faceOwner();
        const labelUList& neiAll = mesh.faceNeighbour();
        const label nIntAll = mesh.nInternalFaces();
        forAll(SfAll, faceI)
        {
            const vector contrib = st.alphaf()[faceI]*SfAll[faceI];
            mAlpha[ownAll[faceI]] -= contrib;
            if (faceI < nIntAll) mAlpha[neiAll[faceI]] += contrib;
        }
        forAll(mAlpha, c)
        {
            if (!st.isInterfaceCell(c)) mAlpha[c] = Zero;
        }
        rAlpha = evaluateNormals
        (
            "davofAlphaSf", mesh, mAlpha, 1.0,
            st.xS().primitiveField(), surface(), st.wispTol(),
            nullptr, nullptr
        );
    }

    // Closure of the closed surface (round-off), tet/face consistency,
    // piecewise-linear patch area, alpha clipping.
    const vector sumM = gSum(st.m().primitiveField());
    const scalar sumMRel = (aExact > 0) ? mag(sumM)/aExact : mag(sumM)/max(rd.A, VSMALL);
    const scalar maxConsistency = st.maxConsistency();
    const scalar aPL = gSum(st.AS().primitiveField());
    const scalar eAreaRel = (aExact > 0) ? mag(rd.A - aExact)/aExact : -1;
    const scalar eAreaPLRel = (aExact > 0) ? mag(aPL - aExact)/aExact : -1;
    const scalar maxAlphaClip = returnReduce(st.maxAlphaClip(), maxOp<scalar>());
    // Realizability of the state by the DAVOF plane (planePosition()).
    const scalar maxVolDiffPlane = st.maxVolDiffPlane();
    const scalar maxAreaDiffPlane = st.maxAreaDiffPlane();

    // alpha.davof against the leiaSetFields field, when present.
    scalar maxAlphaDiffDA = -1;
    {
        IOobject io
        (
            alphaName, runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        );
        if (io.typeHeaderOk<volScalarField>(true))
        {
            volScalarField alphaRef(io, mesh);
            maxAlphaDiffDA = gMax
            (
                mag(st.alpha().primitiveField() - alphaRef.primitiveField())()
            );
        }
    }

    // alpha.davof against the exact plane cut (implicitPlane only).
    scalar maxAlphaDiffPlaneCut = -1;
    if (surfType == "implicitPlane")
    {
        const implicitPlane plane(surfDict);
        const scalarField& V = mesh.V();
        scalar worst = 0;
        forAll(V, c)
        {
            const scalar aCut =
                intersectCell<volumeArea>(c, mesh, plane).volume()/V[c];
            worst = max(worst, mag(st.alpha()[c] - aCut));
        }
        maxAlphaDiffPlaneCut = returnReduce(worst, maxOp<scalar>());
    }

    // ---- Curvature of the DAVOF state ----------------------------------------
    // fvSolution davof.curvature { models (...); <name> { type ...; } }: each
    // model is fitted to the ring of centroids and normals (libleiaDavof
    // curvature/), delivered to the force band and the active faces with the
    // parallel-surface closed form and scored against the exact geometry; the
    // first model is the headline.
    const dictionary& curvDict = davofDict.subOrEmptyDict("curvature");
    wordList curvModels(curvDict.getOrDefault<wordList>("models", wordList()));
    {
        wordList kept;
        for (const word& nm : curvModels) if (nm != "none") kept.append(nm);
        curvModels = kept;
    }
    if (curvModels.size() && !haveExactGeom)
    {
        FatalIOErrorInFunction(curvDict)
            << "davof curvature models need an exact geometry: implicitSphere, "
            << "implicitPlane or signedDistanceEllipsoid, found " << surfType
            << exit(FatalIOError);
    }
    volScalarField eKappa
    (
        IOobject("eKappa.davof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField kappaCellField
    (
        IOobject("kappaCell.davof", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar(dimless/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    List<curvResult> curvResults;
    PtrList<davofCurvature> curvPtrs;
    forAll(curvModels, mi)
    {
        const word& name = curvModels[mi];
        const dictionary& mDict = curvDict.subDict(name);
        curvPtrs.append(davofCurvature::New(name, mesh, st, mDict));
        davofCurvature& cm = curvPtrs.last();
        cm.compute();
        cpuTime dTimer;
        davofCurvatureDelivery dl(mesh, st, cm);
        curvResult cr = evaluateCurvature
        (
            name, mesh, st, cm, dl, exactGeom,
            (mi == 0) ? &eKappa.primitiveFieldRef() : nullptr,
            (mi == 0) ? &kappaCellField.primitiveFieldRef() : nullptr
        );
        cr.cpu = cm.cpuSeconds() + dTimer.cpuTimeIncrement();
        curvResults.append(cr);
        cm.write();
    }
    const bool haveCurv = curvResults.size() > 0;
    const curvResult ch = haveCurv ? curvResults[0] : curvResult();

    // ---- The PLIC surfaces as legacy VTK polydata ----------------------------
    // postProcessing/davofInterface/<time>/plic.<model>.vtk, one polygon per
    // interface cell with the cell index and the error fields attached
    // (libleiaDavofInterface; per rank in parallel).
    const fileName vtkDir =
        runTime.path()/"postProcessing"/"davofInterface"/runTime.timeName();
    {
        boolList isI(mesh.nCells(), false);
        forAll(isI, c) isI[c] = st.isInterfaceCell(c);
        plicInterfaceSurface s
        (
            mesh, st.m().primitiveField(), st.xPlane().primitiveField(), &isI
        );
        wordList vtkNames({"eNormal", "ePos"});
        List<const scalarField*> vtkFields
        ({
            &eDavof.primitiveField(), &ePosDavof.primitiveField()
        });
        if (haveCurv)
        {
            vtkNames.append("eKappa");
            vtkFields.append(&eKappa.primitiveField());
        }
        s.writeLegacyVTK(vtkDir/"plic.davof.vtk", vtkNames, vtkFields);
        Info<< "PLIC surface of davof: " << s.size() << " polygons -> "
            << vtkDir/"plic.davof.vtk" << endl;
    }

    // ---- Cross-check: OpenFOAM geometricVoF schemes on the same alpha -------
    const dictionary& crossDict = davofDict.subOrEmptyDict("crossCheck");
    const word crossSet = crossDict.getOrDefault<word>("set", "none");
    PtrList<volScalarField> eModels;
    PtrList<volScalarField> ePosModels;
    if (crossSet == "geometricVoF")
    {
        const wordList models
        (
            crossDict.getOrDefault<wordList>("models", wordList())
        );
        volVectorField U
        (
            IOobject("U.davofCross", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedVector(dimVelocity, Zero),
            zeroGradientFvPatchVectorField::typeName
        );
        surfaceScalarField phi
        (
            IOobject("phi.davofCross", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar(dimVolume/dimTime, Zero)
        );
        forAll(models, mi)
        {
            const word& model = models[mi];
            const dictionary& mDict = crossDict.subDict(model);
            Info<< nl << "geometricVoF cross-check: " << model << endl;
            cpuTime mTimer;
            // One model per scope: geometricVoF registers interfaceNormal /
            // interfaceCentre / RDF fields by the alpha group, so two models
            // cannot coexist in one registry.
            autoPtr<reconstructionSchemes> rs =
                reconstructionSchemes::New(st.alpha(), phi, U, mDict);
            rs->reconstruct(true);
            const scalar cpuModel = mTimer.cpuTimeIncrement();
            eModels.append
            (
                new volScalarField
                (
                    IOobject("eNormal." + model, runTime.timeName(), mesh,
                             IOobject::NO_READ, IOobject::AUTO_WRITE),
                    mesh, dimensionedScalar(dimless, Zero),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            normResult r = evaluateNormals
            (
                model, mesh, rs->normal().primitiveField(), -1.0,
                rs->centre().primitiveField(), surface(), st.wispTol(),
                &eModels.last().primitiveFieldRef(), nullptr
            );
            r.cpu = cpuModel;
            ePosModels.append
            (
                new volScalarField
                (
                    IOobject("ePos." + model, runTime.timeName(), mesh,
                             IOobject::NO_READ, IOobject::AUTO_WRITE),
                    mesh, dimensionedScalar(dimLength, Zero),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            // Their centre_ is the centroid of the PLIC polygon.
            evaluatePositions
            (
                r, mesh, rs->normal().primitiveField(),
                rs->centre().primitiveField(), st.wispTol(), exactDistance,
                &ePosModels.last().primitiveFieldRef()
            );
            results.append(r);
            {
                // normal_ points into alpha = 1: flip it out of the liquid.
                const vectorField nOut(-rs->normal().primitiveField());
                boolList isI(mesh.nCells(), false);
                forAll(isI, c) isI[c] = mag(nOut[c]) > 0;
                plicInterfaceSurface s
                (
                    mesh, nOut, rs->centre().primitiveField(), &isI
                );
                s.writeLegacyVTK
                (
                    vtkDir/("plic." + model + ".vtk"),
                    wordList({"eNormal", "ePos"}),
                    List<const scalarField*>
                    ({
                        &eModels.last().primitiveField(),
                        &ePosModels.last().primitiveField()
                    })
                );
                Info<< "PLIC surface of " << model << ": " << s.size()
                    << " polygons -> " << vtkDir/("plic." + model + ".vtk")
                    << endl;
            }
        }
    }
    else if (crossSet != "none")
    {
        FatalIOErrorInFunction(crossDict)
            << "Unknown davof crossCheck set '" << crossSet
            << "'. Valid: none, geometricVoF." << exit(FatalIOError);
    }

    // ---- Report --------------------------------------------------------------
    Info<< nl
        << "surface                 : " << surfType << nl
        << "alphaSource             : " << st.alphaSource() << nl
        << "cells / h               : " << nCellsGlobal << " / " << dx << nl
        << "R / h                   : " << ((dx > 0 && radius > 0) ? radius/dx : 0) << nl
        << "interface cells / wisps : " << rd.nInterface << " / " << rd.nWisp << nl
        << "area exact / DAVOF / PL : " << aExact << " / " << rd.A << " / " << aPL
        << "   (rel. error " << eAreaRel << " / " << eAreaPLRel << ")" << nl
        << "|sum m_c| / A_exact     : " << sumMRel << nl
        << "max |m_c - sum caps|/h^2: " << maxConsistency << nl
        << "max alpha clip          : " << maxAlphaClip << nl
        << "max |alpha - " << alphaName << "| : " << maxAlphaDiffDA << nl
        << "max |alpha - plane cut| : " << maxAlphaDiffPlaneCut << nl
        << "plane realizability     : max |alphaPlane - alpha| = " << maxVolDiffPlane
        << ", max |A n - m_c|/h^2 = " << maxAreaDiffPlane << nl
        << "quadratic-face fallbacks: " << st.nFaceFallback() << nl
        << "alpha_f S_f identity    : L2/Linf = "
        << (st.hasMTri() ? rAlpha.L2 : rd.L2) << " / "
        << (st.hasMTri() ? rAlpha.Linf : rd.Linf)
        << (st.hasMTri() ? "  (the davof row uses the face-triangulation sum)" : "  (= davof)") << nl
        << "position (DAVOF foot pt): L2/Linf = " << rd.Pfoot2 << " / " << rd.Pfootinf
        << " m" << nl;
    for (const normResult& r : results)
    {
        Info<< "  " << r.model << "  N_S = " << r.nInterface
            << "  wisps = " << r.nWisp
            << "  A = " << r.A
            << "  L1/L2/Linf = " << r.L1 << " / " << r.L2 << " / " << r.Linf
            << "  (area-weighted " << r.L1aw << " / " << r.L2aw
            << "; Linf incl. wisps " << r.LinfAll
            << "; position L1/L2/Linf = " << r.P1 << " / " << r.P2 << " / " << r.Pinf
            << " m; cpu " << r.cpu << " s)" << nl;
    }
    for (const curvResult& r : curvResults)
    {
        Info<< "  curvature " << r.model
            << "  N_S = " << r.nInterface << "  band = " << r.nBand
            << " (unfilled " << r.nUnfilled << ")  faces = " << r.nActive
            << " (one-sided " << r.nOneSided << ", skipped " << r.nSkipped
            << ")  fallbacks = " << r.nFallback
            << nl
            << "    kappa_ref L2 = " << r.kappaRef2
            << "  centroid L1/L2/Linf = " << r.S1 << " / " << r.S2 << " / " << r.Sinf
            << "  (rel. L2 " << ((r.kappaRef2 > 0) ? r.S2/r.kappaRef2 : 0)
            << "; K rms err " << r.K2 << ")" << nl
            << "    cell L1/L2/Linf = " << r.C1 << " / " << r.C2 << " / " << r.Cinf
            << "  face inverse L1/L2/Linf = " << r.F1 << " / " << r.F2 << " / " << r.Finf
            << "  face foot L1/L2/Linf = " << r.G1 << " / " << r.G2 << " / " << r.Ginf
            << "  (cpu " << r.cpu << " s)" << nl;
    }
    Info<< endl;

    if (Pstream::master())
    {
        const scalar rOverH = (dx > 0 && radius > 0) ? radius/dx : 0;
        OFstream os("leiaTestDavofNormal.csv");
        os.precision(12);
        os << "DELTA_X,N_CELLS_MESH,N_GEOMETRIC_D,ALPHA_SOURCE,WISP_TOL,RADIUS,"
              "R_OVER_H,N_INTERFACE,N_WISP,A_EXACT,A_DAVOF,E_AREA_REL,A_PL,"
              "E_AREA_PL_REL,SUM_M_REL,MAX_CONSISTENCY,MAX_ALPHA_CLIP,"
              "MAX_ALPHA_DIFF_DA,MAX_ALPHA_DIFF_PLANECUT,"
              "E_L1_N,E_L2_N,E_LINF_N,E_L1_N_AW,E_L2_N_AW,E_LINF_N_WISP,"
              "E_POS_L1,E_POS_L2,E_POS_LINF,E_POS_FOOT_L2,E_POS_FOOT_LINF,"
              "MAX_VOL_DIFF_PLANE,MAX_AREA_DIFF_PLANE,N_FACE_FALLBACK,"
              "CPU_SECONDS,"
              "E_KAPPA_L1,E_KAPPA_L2,E_KAPPA_LINF,KAPPA_REF_L2,E_KAPPA_FACE_L2,"
              "N_CURV_FALLBACK" << nl;
        // The headline curvature model's columns (blank without models).
        auto curvCols = [&](Ostream& o)
        {
            if (haveCurv)
            {
                o << ',' << ch.S1 << ',' << ch.S2 << ',' << ch.Sinf << ','
                  << ch.kappaRef2 << ',' << ch.F2 << ',' << ch.nFallback;
            }
            else
            {
                o << ",,,,,,";
            }
        };
        os << dx << ',' << nCellsGlobal << ',' << mesh.nGeometricD() << ','
           << st.alphaSource() << ',' << st.wispTol() << ',' << radius << ','
           << rOverH << ',' << rd.nInterface << ',' << rd.nWisp << ','
           << aExact << ',' << rd.A << ',' << eAreaRel << ',' << aPL << ','
           << eAreaPLRel << ',' << sumMRel << ',' << maxConsistency << ','
           << maxAlphaClip << ',' << maxAlphaDiffDA << ','
           << maxAlphaDiffPlaneCut << ','
           << rd.L1 << ',' << rd.L2 << ',' << rd.Linf << ','
           << rd.L1aw << ',' << rd.L2aw << ',' << rd.LinfAll << ','
           << rd.P1 << ',' << rd.P2 << ',' << rd.Pinf << ','
           << rd.Pfoot2 << ',' << rd.Pfootinf << ','
           << maxVolDiffPlane << ',' << maxAreaDiffPlane << ','
           << st.nFaceFallback() << ','
           << rd.cpu;
        curvCols(os);
        os << nl;

        OFstream osM("leiaTestDavofNormalModels.csv");
        osM.precision(12);
        osM << "MODEL,ALPHA_SOURCE,DELTA_X,N_CELLS_MESH,R_OVER_H,N_INTERFACE,"
               "N_WISP,A_EXACT,A_MODEL,E_AREA_REL,E_L1_N,E_L2_N,E_LINF_N,"
               "E_L1_N_AW,E_L2_N_AW,E_LINF_N_WISP,"
               "E_POS_L1,E_POS_L2,E_POS_LINF,CPU_SECONDS,"
               "E_KAPPA_L1,E_KAPPA_L2,E_KAPPA_LINF,KAPPA_REF_L2,E_KAPPA_FACE_L2,"
               "N_CURV_FALLBACK" << nl;
        for (const normResult& r : results)
        {
            const scalar eA = (aExact > 0) ? mag(r.A - aExact)/aExact : -1;
            osM << r.model << ',' << st.alphaSource() << ',' << dx << ','
                << nCellsGlobal << ',' << rOverH << ',' << r.nInterface << ','
                << r.nWisp << ',' << aExact << ',' << r.A << ',' << eA << ','
                << r.L1 << ',' << r.L2 << ',' << r.Linf << ','
                << r.L1aw << ',' << r.L2aw << ',' << r.LinfAll << ','
                << r.P1 << ',' << r.P2 << ',' << r.Pinf << ','
                << r.cpu;
            // The curvature belongs to the DAVOF state, not to a comparator.
            if (r.model == "davof") curvCols(osM); else osM << ",,,,,,";
            osM << nl;
        }

        OFstream osC("leiaTestDavofCurvature.csv");
        osC.precision(12);
        osC << "CURV_MODEL,ALPHA_SOURCE,DELTA_X,N_CELLS_MESH,R_OVER_H,"
               "N_INTERFACE,N_BAND_CELLS,N_UNFILLED_CELLS,N_ACTIVE_FACES,"
               "N_ONE_SIDED_FACES,N_SKIPPED_FACES,KAPPA_REF_L2,"
               "E_KAPPA_L1,E_KAPPA_L2,E_KAPPA_LINF,E_K_L2,"
               "E_KAPPA_CELL_L1,E_KAPPA_CELL_L2,E_KAPPA_CELL_LINF,"
               "E_KAPPA_FACE_L1,E_KAPPA_FACE_L2,E_KAPPA_FACE_LINF,"
               "E_KAPPA_FACE_FOOT_L1,E_KAPPA_FACE_FOOT_L2,E_KAPPA_FACE_FOOT_LINF,"
               "N_CURV_FALLBACK,CPU_SECONDS" << nl;
        for (const curvResult& r : curvResults)
        {
            osC << r.model << ',' << st.alphaSource() << ',' << dx << ','
                << nCellsGlobal << ',' << rOverH << ','
                << r.nInterface << ',' << r.nBand << ',' << r.nUnfilled << ','
                << r.nActive << ',' << r.nOneSided << ',' << r.nSkipped << ','
                << r.kappaRef2 << ','
                << r.S1 << ',' << r.S2 << ',' << r.Sinf << ',' << r.K2 << ','
                << r.C1 << ',' << r.C2 << ',' << r.Cinf << ','
                << r.F1 << ',' << r.F2 << ',' << r.Finf << ','
                << r.G1 << ',' << r.G2 << ',' << r.Ginf << ','
                << r.nFallback << ',' << r.cpu << nl;
        }
    }

    st.write();
    runTime.writeNow();

    if (expectExact)
    {
        // Round-off, not discretisation: the linear interpolant and the
        // plane-DA state are exact for a plane (measured ~1e-14 for the
        // interpolant, ~1e-13 for the least-squares plane fit).
        const scalar tol = 1e-12;
        // The position in units of h; the realizability diagnostics are
        // dimensionless already. The position p = (3 alpha V - M)/|m|
        // divides the difference of two O(V) numbers by the cap area, which
        // in a corner-sliver cell is ~1e-3 V^(2/3): its round-off floor is
        // ten times the others (measured 1.5e-12 and 2.2e-12 h on perturbed
        // hexahedra, 2026-10-02), hence the ten-fold allowance.
        const scalar posLinfH = (dx > 0) ? rd.Pinf/dx : rd.Pinf;
        const scalar posTol = 10*tol;
        // The curvature of a plane is zero: kappa h at round-off at the
        // centroids and at the faces (every configured model).
        scalar curvLinfH = 0;
        for (const curvResult& r : curvResults)
        {
            curvLinfH = max(curvLinfH, max(r.Sinf, r.Finf)*dx);
        }
        if
        (
            rd.Linf > tol
         || maxConsistency > tol
         || (maxAlphaDiffPlaneCut >= 0 && maxAlphaDiffPlaneCut > tol)
         || posLinfH > posTol
         || maxVolDiffPlane > tol
         || maxAreaDiffPlane > tol
         || curvLinfH > tol
         || rd.nInterface == 0
        )
        {
            FatalErrorInFunction
                << "-expectExact FAILED: E_LINF_N = " << rd.Linf
                << ", MAX_CONSISTENCY = " << maxConsistency
                << ", MAX_ALPHA_DIFF_PLANECUT = " << maxAlphaDiffPlaneCut
                << ", E_POS_LINF/h = " << posLinfH
                << ", MAX_VOL_DIFF_PLANE = " << maxVolDiffPlane
                << ", MAX_AREA_DIFF_PLANE = " << maxAreaDiffPlane
                << ", max E_KAPPA_LINF h = " << curvLinfH
                << ", N_INTERFACE = " << rd.nInterface
                << " (tolerance " << tol << ", position " << posTol << ")"
                << exit(FatalError);
        }
        Info<< "-expectExact PASSED (tolerance " << tol << ", position "
            << posTol << ")" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
