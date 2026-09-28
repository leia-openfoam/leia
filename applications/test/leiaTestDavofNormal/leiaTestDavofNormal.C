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
    with the exact area (4 pi R^2 for implicitSphere).

    Optional cross-check on the SAME alpha_c (fvSolution davof.crossCheck
    { set geometricVoF; models (plicRDF gradAlpha isoAlpha); <model> {...} }):
    OpenFOAM's libgeometricVoF reconstruction schemes are run on alpha.davof
    and scored with the same norms; their normal_ is the PLIC area vector
    pointing INTO alpha = 1, so it is flipped before the comparison, and their
    centre_ is the evaluation point of the exact normal.

    Outputs (case directory): leiaTestDavofNormal.csv (one row, the DAVOF
    result and the consistency diagnostics) and leiaTestDavofNormalModels.csv
    (tidy, one row per MODEL: davof, plicRDF, gradAlpha, isoAlpha); fields
    alpha.davof, alphaf.davof, m.davof, xS.davof, AS.davof, nDavof,
    eNormal.davof, interfaceMarker.davof (0 bulk, 1 interface, 2 wisp),
    eNormal.<model>.

    Parallel-safe (every norm is reduced). Options:
        -alphaName <word>   the leiaSetFields field to compare alpha.davof with
                            (default alpha.water; MAX_ALPHA_DIFF_DA is 0 for
                            alphaSource planePhaseIndicator with geometrySource
                            levelSetField, the same construction)
        -expectExact        exact-solution gate for an implicitPlane: exit
                            non-zero unless the normal error, the tet/face
                            consistency and the alpha difference against the
                            plane cut of src/vof/foamGeometry.H are at round-off.

    Run in a meshed, leiaSetFields-initialised case:
        blockMesh; leiaSetFields -alphaName alpha.water; leiaTestDavofNormal

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cpuTime.H"
#include "davofState.H"
#include "levelSetImplicitSurfaces.H"
#include "leiaVersionRegistry.H"
#include "reconstructionSchemes.H"
#include "foamGeometry.H"

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
    scalar cpu = 0;
};


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

    // ---- The DAVOF state -----------------------------------------------------
    davofState st(mesh, davofDict);
    cpuTime timer;
    autoPtr<volScalarField> psiPtr;
    if (st.alphaSource() == "linearInterpolant")
    {
        st.computeFromImplicitSurface(surface());
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

    List<normResult> results;
    {
        normResult r = evaluateNormals
        (
            "davof", mesh, st.m().primitiveField(), 1.0,
            st.xS().primitiveField(), surface(), st.wispTol(),
            &eDavof.primitiveFieldRef(), &marker.primitiveFieldRef()
        );
        r.cpu = cpuDavof;
        results.append(r);
        forAll(nDavof, c)
        {
            if (st.isInterfaceCell(c))
            {
                nDavof[c] = st.m()[c]/mag(st.m()[c]);
            }
        }
    }
    // A COPY: results grows below (the cross-check appends), which reallocates
    // the list and would leave a reference dangling.
    const normResult rd = results[0];

    // Closure of the closed surface (round-off), tet/face consistency,
    // piecewise-linear patch area, alpha clipping.
    const vector sumM = gSum(st.m().primitiveField());
    const scalar sumMRel = (aExact > 0) ? mag(sumM)/aExact : mag(sumM)/max(rd.A, VSMALL);
    const scalar maxConsistency = st.maxConsistency();
    const scalar aPL = gSum(st.AS().primitiveField());
    const scalar eAreaRel = (aExact > 0) ? mag(rd.A - aExact)/aExact : -1;
    const scalar eAreaPLRel = (aExact > 0) ? mag(aPL - aExact)/aExact : -1;
    const scalar maxAlphaClip = returnReduce(st.maxAlphaClip(), maxOp<scalar>());

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

    // ---- Cross-check: OpenFOAM geometricVoF schemes on the same alpha -------
    const dictionary& crossDict = davofDict.subOrEmptyDict("crossCheck");
    const word crossSet = crossDict.getOrDefault<word>("set", "none");
    PtrList<volScalarField> eModels;
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
            results.append(r);
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
        << "max |alpha - plane cut| : " << maxAlphaDiffPlaneCut << nl;
    for (const normResult& r : results)
    {
        Info<< "  " << r.model << "  N_S = " << r.nInterface
            << "  wisps = " << r.nWisp
            << "  A = " << r.A
            << "  L1/L2/Linf = " << r.L1 << " / " << r.L2 << " / " << r.Linf
            << "  (area-weighted " << r.L1aw << " / " << r.L2aw
            << "; Linf incl. wisps " << r.LinfAll
            << "; cpu " << r.cpu << " s)" << nl;
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
              "CPU_SECONDS" << nl;
        os << dx << ',' << nCellsGlobal << ',' << mesh.nGeometricD() << ','
           << st.alphaSource() << ',' << st.wispTol() << ',' << radius << ','
           << rOverH << ',' << rd.nInterface << ',' << rd.nWisp << ','
           << aExact << ',' << rd.A << ',' << eAreaRel << ',' << aPL << ','
           << eAreaPLRel << ',' << sumMRel << ',' << maxConsistency << ','
           << maxAlphaClip << ',' << maxAlphaDiffDA << ','
           << maxAlphaDiffPlaneCut << ','
           << rd.L1 << ',' << rd.L2 << ',' << rd.Linf << ','
           << rd.L1aw << ',' << rd.L2aw << ',' << rd.LinfAll << ','
           << rd.cpu << nl;

        OFstream osM("leiaTestDavofNormalModels.csv");
        osM.precision(12);
        osM << "MODEL,ALPHA_SOURCE,DELTA_X,N_CELLS_MESH,R_OVER_H,N_INTERFACE,"
               "N_WISP,A_EXACT,A_MODEL,E_AREA_REL,E_L1_N,E_L2_N,E_LINF_N,"
               "E_L1_N_AW,E_L2_N_AW,E_LINF_N_WISP,CPU_SECONDS" << nl;
        for (const normResult& r : results)
        {
            const scalar eA = (aExact > 0) ? mag(r.A - aExact)/aExact : -1;
            osM << r.model << ',' << st.alphaSource() << ',' << dx << ','
                << nCellsGlobal << ',' << rOverH << ',' << r.nInterface << ','
                << r.nWisp << ',' << aExact << ',' << r.A << ',' << eA << ','
                << r.L1 << ',' << r.L2 << ',' << r.Linf << ','
                << r.L1aw << ',' << r.L2aw << ',' << r.LinfAll << ','
                << r.cpu << nl;
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
        if
        (
            rd.Linf > tol
         || maxConsistency > tol
         || (maxAlphaDiffPlaneCut >= 0 && maxAlphaDiffPlaneCut > tol)
         || rd.nInterface == 0
        )
        {
            FatalErrorInFunction
                << "-expectExact FAILED: E_LINF_N = " << rd.Linf
                << ", MAX_CONSISTENCY = " << maxConsistency
                << ", MAX_ALPHA_DIFF_PLANECUT = " << maxAlphaDiffPlaneCut
                << ", N_INTERFACE = " << rd.nInterface
                << " (tolerance " << tol << ")" << exit(FatalError);
        }
        Info<< "-expectExact PASSED (tolerance " << tol << ")" << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
