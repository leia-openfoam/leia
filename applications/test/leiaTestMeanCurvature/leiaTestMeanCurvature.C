/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestMeanCurvature

Description
    Standalone accuracy test for the mean curvature evaluated symbolically from the
    semi-Lagrangian quadratic reconstruction (slAdvection::meanCurvature, i.e.
    kappa = div(grad psi/|grad psi|) from the fit's gradient and Hessian), on a
    circular interface whose exact curvature is known.

    For a circle/sphere of radius R initialised by leiaSetFields (signed distance),
    the exact interface curvature is kappa_exact = (d-1)/R (1/R in 2D, 2/R in 3D).
    The error of the reconstructed curvature is measured in the band of cells that
    carry the surface-tension force -- those adjacent to a face with a non-zero
    snGrad(alpha) -- since only there does the curvature enter the momentum balance.
    Reports the L1, L2 and Linf norms (volume-weighted) of kappa - kappa_exact for
      * div  : kappa = div(grad psi/|grad psi|)     (the method used by the solver);
      * lap  : kappa = Laplacian(psi) = tr(H)        (the signed-distance simplification),
    and writes leiaTestMeanCurvature.csv (one row) for the study aggregator.

    FACE-CENTERED section (leiaTestFaceCurvature.csv, tidy: one row per
    MODEL x FOOT_POINT): the quantity the CSF force actually applies is the FACE
    curvature kappa_f in G_sigma,f = sigma*kappa_f*snGrad(alpha)*|Sf|, on the
    ACTIVE faces (|snGrad(alpha)| > 0). For every curvature model the app
    assembles kappa_f exactly as the corresponding surfaceTensionForce path
    does (arithmetic fvc::interpolate for cell models, the Kang/GFM inverse-psi
    weighting, the connectedInterface native face field, the FVM div/trace face
    formulas) and measures |Sf|-weighted L1/L2/Linf plus the force-weighted
    (|snGrad(alpha)||Sf|) L2 of kappa_f - kappa_exact over the active faces.

    FOOT_POINT = 1 rows re-deliver the same kappa_f at the interface through
    the STABILIZED FOOT-POINT algorithm (slReconstruction::footPointDistance,
    stable-foot-point-3d.md): the signed face-centre offset d_f (owner/
    neighbour fits, inverse-|psi| Kang weighting) feeds the parallel-curve
    inverse kappa = kappa_d/(1 - d_f*kappa_d) (|1 - d*kappa| > 1/2 guarded).
    The correction is emitted only for CONTOUR-referenced models (their kappa_f
    is the curvature of the level contour through the face centre); interface-
    referenced models (Kang, foot extensions, foot-point HF, connected, iso,
    interface mean) evaluate at d ~ 0 by construction, so the correction does
    not apply to them and only the FOOT_POINT = 0 row is written. The
    stableFootPoint model is fully foot-point-native: the fit curvature of each
    adjacent cell evaluated AT the face centre, offset-inverted with that
    cell's own stabilized foot-point distance, then Kang-combined.

    2D only (footPointCurvature/connectedInterfaceCurvature requirement) and
    meant to run SERIAL (the connected reconstruction is serial; face rows are
    computed over internal faces).

    Run in a meshed, leiaSetFields-initialised case:
        blockMesh; leiaSetFields -alphaName alpha.water; leiaTestMeanCurvature

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvSolution.H"
#include "slAdvection.H"
#include "narrowBand.H"

// Foot-point height-function curvature (solver-side header; -I in Make/options).
#include "footPointCurvature.H"
// Iso-agnostic geometric (multi-cell circle-fit) curvature (library header).
#include "isoGeometricCurvature.H"
// Connected zero-surface curvature with native face delivery (solver-side).
#include "connectedInterfaceCurvature.H"
// CSF-support-weighted spatially constant curvature diagnostic (solver-side).
#include "interfaceMeanCurvature.H"
// The production face delivery's shared helpers (solver-side): the implicit
// Gaussian curvature of a cell's fit and the dimension-seamless
// parallel-surface inverse -- reused here so the gate measures EXACTLY the
// formula the solver applies.
#include "stabilizedFootPointFaceCurvature.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "alphaName", "word",
        "Name of the volume-fraction field for the snGrad band (default alpha.water)."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word alphaName = args.getOrDefault<word>("alphaName", "alpha.water");

    // The face-centered section covers INTERNAL faces only: in a decomposed
    // run the processor-seam faces (active faces of a centered droplet!)
    // would silently vanish from every norm while the reduce() scaffolding
    // still produced plausible global numbers. Enforce the serial contract
    // instead of documenting it away.
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "leiaTestMeanCurvature is a SERIAL test: the face-centered "
            << "curvature norms cover internal faces only, so a parallel run "
            << "drops the processor-seam part of the CSF force support "
            << "(and the connected-interface model is serial anyway)."
            << exit(FatalError);
    }

    Info<< "Reading level set psi and volume fraction " << alphaName << nl << endl;
    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField alpha
    (
        IOobject(alphaName, runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    // Exact interface curvature of the initialised circle/sphere.
    const fvSolution& fvSol(mesh);
    const dictionary& implSurf =
        fvSol.subDict("levelSet").subDict("implicitSurface");
    const scalar R = implSurf.get<scalar>("radius");
    const label nd = mesh.nGeometricD();
    const scalar kappaExact = scalar(nd - 1)/R;   // 1/R in 2D, 2/R in 3D
    Info<< "R = " << R << ", geometric dims = " << nd
        << ", kappa_exact = " << kappaExact << nl << endl;

    // The FOOT_POINT = 1 face rows apply the parallel-curve inverse to the
    // interpolated CELL curvature, which must therefore be the raw
    // contour-referenced value: with fvSolution levelSet.semiLagrangian
    // offsetCorrection != none the reconstruction ALREADY offset-corrects at
    // the cell level and the face correction would be applied twice --
    // silently wrong "with foot point" rows. Fail loud instead.
    {
        const word offsetCorrection =
            fvSol.subDict("levelSet").subOrEmptyDict("semiLagrangian")
                 .getOrDefault<word>("offsetCorrection", "none");
        if (offsetCorrection != "none")
        {
            FatalErrorInFunction
                << "leiaTestMeanCurvature requires offsetCorrection none "
                << "(found " << offsetCorrection << "): the face-centered "
                << "FOOT_POINT rows apply the parallel-curve inverse "
                << "themselves and would double-correct."
                << exit(FatalError);
        }
    }

    // Curvature from the reconstruction: full div(n) and the SDF simplification tr(H).
    autoPtr<slAdvection> slAdv = slAdvection::New(mesh);

    volScalarField kappaDiv
    (
        IOobject("kappaDiv", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    // AUTO_WRITE (not a plain rename-copy of kappaDiv, which would not be written)
    // so the no-extension / Laplacian curvature fields land on disk for the
    // pole-error field analysis.
    volScalarField kappaLap
    (
        IOobject("kappaLap", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    volScalarField kappaNoExt
    (
        IOobject("kappaNoExt", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );

    slAdv->meanCurvature(psi, kappaDiv);              // div(n), normal-extended (foot)
    slAdv->meanCurvatureLaplacian(psi, kappaLap);     // tr(H), SDF simplification
    slAdv->meanCurvatureNoExtension(psi, kappaNoExt); // div(n) at cell centre, no extension

    // Trust-region closest-point (KKT) Newton foot on the cell's own quadratic
    // (the solver's curvatureExtension closestPointNewton path).
    volScalarField kappaCP
    (
        IOobject("kappaCP", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    slAdv->meanCurvatureClosestPoint(psi, kappaCP);

    // 2D-only models (foot-point height function, iso-geometric,
    // connectedInterface) FatalError on 3D meshes -- gate them and skip their
    // rows/columns in 3D (sentinel-empty CSV cells; the fig script tolerates
    // missing rows).
    const bool is2D = (mesh.nGeometricD() == 2);

    // Foot-point HEIGHT-FUNCTION curvature (footPointCurvature.H): local parabola
    // fit of the linear-plane interface feet, delivered normal-constant by uniform
    // nearest-foot assignment. Needs the NarrowBand field (signChange model).
    autoPtr<narrowBand> nb = narrowBand::New(mesh, psi);
    nb->calc();
    volScalarField kappaFoot
    (
        IOobject("kappaFoot", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        kappaNoExt   // fallback = the symbolic cell-centre curvature
    );
    kappaFoot.rename("kappaFoot");
    const label nFootFallback =
        is2D
      ? computeFootPointCurvature
        (
            mesh, psi,
            mesh.lookupObject<volScalarField>("NarrowBand"),
            kappaFoot
        )
      : -1;

    // Face-average, cell-centred curvature: interpolate the normal-extended kappa to
    // the faces (as the CSF force does, sigma*interpolate(kappa)*snGrad(alpha)) and
    // average back to the cell centres -- the effective curvature the force applies.
    volScalarField kappaFaceAvg
    (
        IOobject("kappaFaceAvg", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        fvc::average(fvc::interpolate(kappaDiv))
    );

    // ISO-AGNOSTIC GEOMETRIC curvature (isoGeometricCurvature.H) from the psi=0 face-edge
    // crossings, two ways: the CONORMAL integral (kappa n)_c = sum_f L_f x n_f with the
    // owner/neighbour-mean normals (kappaIso, method 0) and the LOCAL tangent-frame
    // parabola height-function fit (kappaParab, method 1). psi is a signed distance
    // increasing outward -> sign +1; extend two layers to cover the snGrad(alpha) band.
    volScalarField kappaIso
    (
        IOobject("kappaIso", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    volScalarField kappaParab
    (
        IOobject("kappaParab", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    if (is2D)
    {
        boolList fl; label hole = 0;
        const label gCon = computeIsoGeometricCurvature
        (
            mesh, psi, /*isoValue*/ 0.0, /*sign*/ 1.0,
            /*nExtend*/ 2, /*nSmooth*/ 0, kappaIso, fl, hole, /*method*/ 0
        );
        Info<< "iso-CONORMAL curvature: " << gCon << " cells, " << hole << " skipped" << nl;
        const label gPar = computeIsoGeometricCurvature
        (
            mesh, psi, /*isoValue*/ 0.0, /*sign*/ 1.0,
            /*nExtend*/ 2, /*nSmooth*/ 0, kappaParab, fl, hole, /*method*/ 1
        );
        Info<< "iso-PARABOLA curvature: " << gPar << " cells, " << hole << " skipped" << nl;
        kappaIso.correctBoundaryConditions();
        kappaParab.correctBoundaryConditions();
    }

    // CONNECTED zero-surface curvature (connectedInterfaceCurvature.H): one
    // value per conforming interface element, tangential regularisation along
    // the chain (fvSolution levelSet.curvatureExtension controls; defaults
    // fitHalfWidth 3, helmholtz lambda 1, 30 iterations, estimator
    // connectedFit), Gauss-Bonnet additive-mode fix, and NORMAL-RAY delivery
    // directly to the face field -- the production faceInterpolation
    // connectedInterface path of reconstructedCurvature. Serial only.
    volScalarField kappaConn
    (
        IOobject("kappaConn", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        kappaNoExt   // fallback for deficient elements = the symbolic value
    );
    kappaConn.rename("kappaConn");
    surfaceScalarField kappaConnFace
    (
        IOobject("kappaConnFace", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0)
    );
    bool haveConnected = false;
    if (!Pstream::parRun() && is2D)
    {
        const label nConnDeficient = computeConnectedInterfaceCurvature
        (
            mesh, psi,
            mesh.lookupObject<volScalarField>("NarrowBand"),
            kappaConn,
            kappaConnFace
        );
        haveConnected = true;
        Info<< "connected-interface curvature: " << nConnDeficient
            << " deficient fits/open endpoints" << nl;
    }
    else
    {
        Info<< "connected-interface curvature: SKIPPED ("
            << (is2D ? "serial only" : "2D only") << ")" << nl;
    }

    // CSF-support-weighted spatially CONSTANT curvature (the interfaceMean
    // diagnostic): the |snGrad(alpha)||Sf|-weighted mean of the reconstructed
    // cell curvature, uniform over the whole domain.
    volScalarField kappaMean
    (
        IOobject("kappaMean", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        kappaNoExt
    );
    kappaMean.rename("kappaMean");
    const scalar kappaBar = applyInterfaceMeanCurvature(mesh, alpha, kappaMean);
    Info<< "interface-mean curvature kappaBar = " << kappaBar
        << " (exact " << kappaExact << ")" << nl << endl;

    // Band of cells that feed the CSF force: those touching an internal face with
    // a non-zero snGrad(alpha). Only there does the curvature enter the momentum
    // balance (sigma*interpolate(kappa)*snGrad(alpha)).
    const surfaceScalarField snAlpha(fvc::snGrad(alpha));
    const scalarField& snI = snAlpha.primitiveField();
    const scalar snMax = max(gMax(mag(snI)), VSMALL);
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    boolList band(mesh.nCells(), false);
    forAll(snI, f)
    {
        if (mag(snI[f]) > 1e-8*snMax)
        {
            band[own[f]] = true;
            band[nei[f]] = true;
        }
    }

    // Volume-weighted error norms over the band, for every estimator. A band cell
    // whose curvature was not computed (kappa = 0) contributes its full error, as
    // that is exactly what the momentum equation would see there.
    const scalarField& V = mesh.V();
    scalar sumV = 0, nBand = 0;
    forAll(band, c) { if (band[c]) { sumV += V[c]; nBand += 1; } }
    reduce(sumV, sumOp<scalar>()); reduce(nBand, sumOp<scalar>());

    auto bandError = [&](const volScalarField& kf,
                         scalar& L1, scalar& L2, scalar& Linf)
    {
        scalar s1 = 0, s2 = 0, li = 0;
        forAll(band, c)
        {
            if (!band[c]) { continue; }
            const scalar e = mag(kf[c] - kappaExact);
            s1 += V[c]*e; s2 += V[c]*e*e; li = max(li, e);
        }
        reduce(s1, sumOp<scalar>()); reduce(s2, sumOp<scalar>());
        reduce(li, maxOp<scalar>());
        L1 = (sumV > 0) ? s1/sumV : 0;
        L2 = (sumV > 0) ? Foam::sqrt(s2/sumV) : 0;
        Linf = li;
    };

    scalar L1d, L2d, linfd, L1l, L2l, linfl;
    scalar L1n, L2n, linfn, L1f, L2f, linff;
    scalar L1p, L2p, linfp, L1i, L2i, linfi, L1b, L2b, linfb;
    bandError(kappaDiv,     L1d, L2d, linfd);   // div(n), normal-extended
    bandError(kappaLap,     L1l, L2l, linfl);   // tr(H)
    bandError(kappaNoExt,   L1n, L2n, linfn);   // div(n), no extension
    bandError(kappaFaceAvg, L1f, L2f, linff);   // face-averaged extended kappa
    bandError(kappaFoot,    L1p, L2p, linfp);   // foot-point height function
    bandError(kappaIso,     L1i, L2i, linfi);   // iso-agnostic geometric CONORMAL
    bandError(kappaParab,   L1b, L2b, linfb);   // iso-agnostic geometric PARABOLA

    // VORTICAL-DRIVER metric: the balanced-force CSF sigma*interp(kappa)*snGrad(alpha)
    // is exactly balanced by pressure for face-UNIFORM kappa; its spurious-current
    // driver is the NON-GRADIENT part, whose 2D vorticity source is
    //     omega_dot ~ sigma*(grad kappa x grad alpha)/rho.
    // A smooth radial kappa offset (grad kappa || grad alpha) drives ~nothing (it only
    // shifts the Laplace jump); a grid-scale TANGENTIAL kappa oscillation
    // (grad kappa perp grad alpha) is a maximal vortical forcing. Report the band
    // max/L2 of |grad kappa x grad alpha| per estimator [1/m^3].
    const volVectorField gradAlphaV(fvc::grad(alpha));
    auto vortDriver = [&](const volScalarField& kf, scalar& VL2, scalar& Vinf)
    {
        const volVectorField gk(fvc::grad(kf));
        scalar s2 = 0, li = 0;
        forAll(band, c)
        {
            if (!band[c]) { continue; }
            const scalar w = mag(gk[c] ^ gradAlphaV[c]);
            s2 += V[c]*w*w; li = max(li, w);
        }
        reduce(s2, sumOp<scalar>()); reduce(li, maxOp<scalar>());
        VL2 = (sumV > 0) ? Foam::sqrt(s2/sumV) : 0;
        Vinf = li;
    };
    scalar Vd2, Vdi, Vp2, Vpi, Vi2, Vii, Vb2, Vbi;
    vortDriver(kappaDiv,   Vd2, Vdi);
    Info<< "vortical driver |grad(kappa) x grad(alpha)| (band L2 / Linf) [1/m^3]:" << nl
        << "  div(n) extended : " << Vd2 << " / " << Vdi << nl;
    if (is2D)
    {
        vortDriver(kappaFoot,  Vp2, Vpi);
        vortDriver(kappaIso,   Vi2, Vii);
        vortDriver(kappaParab, Vb2, Vbi);
        Info<< "  foot-point HF   : " << Vp2 << " / " << Vpi << nl
            << "  iso-CONORMAL    : " << Vi2 << " / " << Vii << nl
            << "  iso-PARABOLA    : " << Vb2 << " / " << Vbi << nl;
    }
    Info<< endl;

    // ======================= FACE-CENTERED CURVATURE =======================
    // The CSF force applies curvature at FACES:
    //     G_sigma,f = sigma * kappa_f * snGrad(alpha) * |Sf|,
    // so the force-relevant accuracy is that of kappa_f on the ACTIVE faces
    // (|snGrad(alpha)| > 0) -- measured here for every model's face delivery,
    // each with and without the stabilized foot-point interface referencing.

    const label nIntFaces = mesh.nInternalFaces();
    boolList activeFace(nIntFaces, false);
    label nActive = 0;
    forAll(activeFace, f)
    {
        if (mag(snI[f]) > 1e-8*snMax) { activeFace[f] = true; ++nActive; }
    }
    reduce(nActive, sumOp<label>());

    const vectorField& Cf = mesh.Cf().primitiveField();
    const scalarField& magSfIn = mesh.magSf().primitiveField();
    const scalarField& psiIn = psi.primitiveField();

    // Stabilized foot-point signed offset of every active face centre from the
    // interface (stable-foot-point-3d.md engine on the owner/neighbour fits,
    // slReconstruction::footPointDistance) PLUS the implicit-surface GAUSSIAN
    // curvature of the adjacent fits (fitGaussianCurvature) -- both combined
    // with the inverse-|psi| Kang weighting so the near-interface side's fit
    // dominates. K is exactly +0 on pseudo-2D fits, so the general
    // parallel-SURFACE inverse below is bit-identical to the measured 2D
    // parallel-curve form there; in 3D K is load-bearing (without it the
    // inverse converges to 2/(R - d) on a sphere -- first order).
    const slReconstruction& recon = slAdv->reconstruction();
    scalarField dFace(nIntFaces, 0.0);
    scalarField KFace(nIntFaces, 0.0);
    boolList dOk(nIntFaces, false);
    label nFootFail = 0, nGaussFallback = 0;
    forAll(activeFace, f)
    {
        if (!activeFace[f]) { continue; }

        const scalar wP = mag(psiIn[nei[f]]);   // weight of the OWNER values
        const scalar wN = mag(psiIn[own[f]]);

        bool okO = false, okN = false;
        const scalar dO = recon.footPointDistance(own[f], Cf[f], 0.0, okO);
        const scalar dN = recon.footPointDistance(nei[f], Cf[f], 0.0, okN);
        if (okO && okN)
        {
            dFace[f] = (dO*wP + dN*wN)/(wP + wN + VSMALL);
            dOk[f] = true;
        }
        else if (okO) { dFace[f] = dO; dOk[f] = true; }
        else if (okN) { dFace[f] = dN; dOk[f] = true; }
        else { ++nFootFail; }

        bool okKO = false, okKN = false;
        const scalar KO = fitGaussianCurvature(recon, own[f], okKO);
        const scalar KN = fitGaussianCurvature(recon, nei[f], okKN);
        if (okKO && okKN) { KFace[f] = (KO*wP + KN*wN)/(wP + wN + VSMALL); }
        else if (okKO)    { KFace[f] = KO; }
        else if (okKN)    { KFace[f] = KN; }
        else              { ++nGaussFallback; }   // K = 0: the 2D form
    }
    reduce(nFootFail, sumOp<label>());
    reduce(nGaussFallback, sumOp<label>());

    // Dimension-seamless parallel-SURFACE inverse at the face (the exact
    // formula the solver delivery applies -- parallelSurfaceInverse from
    // stabilizedFootPointFaceCurvature.H), guarded past half the local
    // curvature radius. Faces without a converged foot stay uncorrected.
    auto footCorrect = [&](const scalarField& kf) -> tmp<scalarField>
    {
        auto tres = tmp<scalarField>::New(kf);
        scalarField& res = tres.ref();
        forAll(activeFace, f)
        {
            if (!activeFace[f] || !dOk[f]) { continue; }
            res[f] = parallelSurfaceInverse(kf[f], dFace[f], KFace[f]);
        }
        return tres;
    };

    // CONTROL: the 2D scalar inverse WITHOUT the Gaussian term. In 2D it is
    // identical to footCorrect (K = +0 exactly); in 3D it must FAIL at first
    // order (2/(R - d) on the sphere) -- the gate asserts that failure.
    auto footCorrectScalar2D = [&](const scalarField& kf) -> tmp<scalarField>
    {
        auto tres = tmp<scalarField>::New(kf);
        scalarField& res = tres.ref();
        forAll(activeFace, f)
        {
            if (!activeFace[f] || !dOk[f]) { continue; }
            res[f] = parallelSurfaceInverse(kf[f], dFace[f], scalar(0));
        }
        return tres;
    };

    // Arithmetic (fvc::interpolate) face values of the cell models -- the
    // reconstructedCurvature faceInterpolation=arithmetic CSF path.
    const scalarField kfQ(fvc::interpolate(kappaNoExt)().primitiveField());
    const scalarField kfLap(fvc::interpolate(kappaLap)().primitiveField());
    const scalarField kfExt(fvc::interpolate(kappaDiv)().primitiveField());
    const scalarField kfCP(fvc::interpolate(kappaCP)().primitiveField());
    const scalarField kfFoot(fvc::interpolate(kappaFoot)().primitiveField());
    const scalarField kfIso(fvc::interpolate(kappaIso)().primitiveField());
    const scalarField kfParab(fvc::interpolate(kappaParab)().primitiveField());
    const scalarField kfMean(fvc::interpolate(kappaMean)().primitiveField());

    // Kang/GFM inverse-|psi| interpolation TO the psi=0 crossing -- the
    // reconstructedCurvature faceInterpolation=interfaceWeighted path.
    scalarField kfKang(nIntFaces, 0.0);
    {
        const scalarField& kIn = kappaNoExt.primitiveField();
        forAll(kfKang, f)
        {
            const scalar wP = mag(psiIn[nei[f]]);
            const scalar wN = mag(psiIn[own[f]]);
            kfKang[f] = (kIn[own[f]]*wP + kIn[nei[f]]*wN)/(wP + wN + VSMALL);
        }
    }

    // FVM face curvatures exactly as the divGrad*/traceGrad* force models
    // assemble them (epsilon-guarded normal fields, then interpolate).
    tmp<volVectorField> tgradPsi = fvc::grad(psi);
    const dimensionedScalar deltaNP("deltaN", tgradPsi().dimensions(), SMALL);
    const volVectorField nPsi("nPsi", tgradPsi()/(mag(tgradPsi()) + deltaNP));
    const scalarField kfFvmDivPsi
    (
        fvc::interpolate(fvc::div(nPsi))().primitiveField()
    );
    const scalarField kfFvmTracePsi
    (
        fvc::interpolate(tr(fvc::grad(nPsi)))().primitiveField()
    );

    // interFoam-style alpha curvature: smoothed alpha (the divGradAlpha model
    // default nSmooth = 2), mesh-scale deltaN, kappa_f = interpolate(-div(n)).
    volScalarField alphaS("alphaS", alpha);
    for (label k = 0; k < 2; ++k)
    {
        alphaS = fvc::average(fvc::interpolate(alphaS));
    }
    tmp<volVectorField> tgradA = fvc::grad(alphaS);
    const dimensionedScalar deltaNA
    (
        "deltaN",
        tgradA().dimensions(),
        1e-8/Foam::cbrt(average(mesh.V()).value())
    );
    const volVectorField nAlpha("nAlpha", tgradA()/(mag(tgradA()) + deltaNA));
    const scalarField kfFvmDivAlpha
    (
        fvc::interpolate(-fvc::div(nAlpha))().primitiveField()
    );

    // Fully foot-point-NATIVE face curvature. A quadratic fit has a CONSTANT
    // Hessian, so its curvature is trustworthy only where the fit is centred:
    // evaluating kappa at the face centre and correcting with the face
    // centre's offset leaves an O(h) referencing residual (measured p ~ 1.1).
    // Instead each adjacent cell contributes its CELL-CENTRE fit curvature,
    // offset-inverted to the interface with that cell's own stabilized
    // foot-point distance OF THE CELL CENTRE -- per-side interface curvature
    // BEFORE the face combination -- then the two sides are Kang-combined.
    scalarField kfStable(nIntFaces, 0.0);
    label nStableUnset = 0;
    {
        const vectorField& C = mesh.C().primitiveField();
        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }

            scalar kSide[2] = {0, 0};
            bool okSide[2] = {false, false};
            const label cells[2] = {own[f], nei[f]};
            for (label s = 0; s < 2; ++s)
            {
                const label c = cells[s];
                vector g(Zero);
                symmTensor H(Zero);
                if (recon.fitDerivatives(c, g, H) < 2) { continue; }

                const scalar gm = mag(g);
                if (gm < SMALL) { continue; }

                scalar kd = (tr(H)*gm*gm - (g & (H & g)))/(gm*gm*gm);
                const scalar Kc = (g & (cof(H) & g))/(gm*gm*gm*gm);
                bool okD = false;
                const scalar d = recon.footPointDistance(c, C[c], 0.0, okD);
                if (okD)
                {
                    kd = parallelSurfaceInverse(kd, d, Kc);
                }
                kSide[s] = kd;
                okSide[s] = true;
            }

            if (okSide[0] && okSide[1])
            {
                const scalar wP = mag(psiIn[nei[f]]);
                const scalar wN = mag(psiIn[own[f]]);
                kfStable[f] = (kSide[0]*wP + kSide[1]*wN)/(wP + wN + VSMALL);
            }
            else if (okSide[0]) { kfStable[f] = kSide[0]; }
            else if (okSide[1]) { kfStable[f] = kSide[1]; }
            else { ++nStableUnset; }
        }
    }
    reduce(nStableUnset, sumOp<label>());

    // CUT-CELL delivery (dimension-agnostic, 2D and 3D identical): ONE
    // interface curvature per cut cell (0 < alpha < 1), obtained by inverting
    // that cell's own centre curvature with its own foot-point distance, then
    // ASSIGNED to the cell's active faces -- the mean where both sides are cut
    // cells. All active faces of a cut cell then carry the same value, so the
    // curvature difference across the force support is identically zero.
    scalarField kfCutCell(kfQ);        // fallback: the interpolated value
    label nCutCells = 0, nNoCutSide = 0, nCutFootFail = 0;
    {
        const vectorField& C = mesh.C().primitiveField();
        const scalarField& aIn = alpha.primitiveField();
        scalarField kappaGamma(mesh.nCells(), 0.0);
        boolList isCut(mesh.nCells(), false);

        forAll(aIn, c)
        {
            if (aIn[c] <= SMALL || aIn[c] >= 1 - SMALL) { continue; }

            vector g(Zero);
            symmTensor H(Zero);
            if (recon.fitDerivatives(c, g, H) < 2) { continue; }
            const scalar gm = mag(g);
            if (gm < SMALL) { continue; }

            const scalar kc = (tr(H)*gm*gm - (g & (H & g)))/(gm*gm*gm);
            const scalar Kc = (g & (cof(H) & g))/(gm*gm*gm*gm);
            bool okD = false;
            const scalar dc = recon.footPointDistance(c, C[c], 0.0, okD);
            if (!okD) { ++nCutFootFail; continue; }

            kappaGamma[c] = parallelSurfaceInverse(kc, dc, Kc);
            isCut[c] = true;
            ++nCutCells;
        }

        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }
            const bool cO = isCut[own[f]], cN = isCut[nei[f]];
            if (cO && cN)
            {
                kfCutCell[f] = 0.5*(kappaGamma[own[f]] + kappaGamma[nei[f]]);
            }
            else if (cO) { kfCutCell[f] = kappaGamma[own[f]]; }
            else if (cN) { kfCutCell[f] = kappaGamma[nei[f]]; }
            else { ++nNoCutSide; }
        }
        reduce(nCutCells, sumOp<label>());
        reduce(nNoCutSide, sumOp<label>());
        reduce(nCutFootFail, sumOp<label>());
    }

    // CELL-MEAN delivery: the PER-FACE inversion (footCorrect, i.e. exactly the
    // production stabilizedFootPointFace values), averaged over the active faces
    // of each cut cell and assigned back to them. Structurally identical to
    // cutCellInverse -- one value per cut cell -- but assembled from n averaged
    // inversions instead of one concentrated inversion. Measured separately
    // because the gain, not the across-support structure, is what governs the
    // coupled growth rate (plan sec. 10), and averaging lowers the gain.
    const scalarField kfPerFace(footCorrect(kfQ)());
    scalarField kfCellMean(kfPerFace);
    label nMeanCells = 0, nMeanNoCutSide = 0;
    {
        const scalarField& aIn = alpha.primitiveField();
        scalarField fSum(mesh.nCells(), 0.0), fCnt(mesh.nCells(), 0.0);
        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }
            fSum[own[f]] += kfPerFace[f]; fCnt[own[f]] += 1.0;
            fSum[nei[f]] += kfPerFace[f]; fCnt[nei[f]] += 1.0;
        }

        scalarField cellMean(mesh.nCells(), 0.0);
        boolList isCutM(mesh.nCells(), false);
        forAll(aIn, c)
        {
            if (aIn[c] <= SMALL || aIn[c] >= 1 - SMALL) { continue; }
            if (fCnt[c] < 0.5) { continue; }
            cellMean[c] = fSum[c]/fCnt[c];
            isCutM[c] = true;
            ++nMeanCells;
        }

        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }
            const bool cO = isCutM[own[f]], cN = isCutM[nei[f]];
            if (cO && cN) { kfCellMean[f] = 0.5*(cellMean[own[f]] + cellMean[nei[f]]); }
            else if (cO)  { kfCellMean[f] = cellMean[own[f]]; }
            else if (cN)  { kfCellMean[f] = cellMean[nei[f]]; }
            else          { ++nMeanNoCutSide; }
        }
        reduce(nMeanCells, sumOp<label>());
        reduce(nMeanNoCutSide, sumOp<label>());
    }

    // |Sf|-weighted error norms over the active faces (plus the force-weighted
    // |snGrad(alpha)||Sf| L2 -- the norm in which the error enters G_sigma,f).
    // An active face without a computed value (kappa_f = 0) contributes its
    // full error: that is exactly what the momentum equation would see there.
    scalar sumAf = 0, sumWf = 0;
    forAll(activeFace, f)
    {
        if (!activeFace[f]) { continue; }
        sumAf += magSfIn[f];
        sumWf += mag(snI[f])*magSfIn[f];
    }
    reduce(sumAf, sumOp<scalar>()); reduce(sumWf, sumOp<scalar>());

    struct faceRow
    {
        word model;
        label footPoint;
        scalar L1, L2, Linf, L2w;
        label nZero;
    };
    DynamicList<faceRow> faceRows;

    auto addFaceRow =
        [&](const word& model, const label fp, const scalarField& kf)
    {
        scalar s1 = 0, s2 = 0, li = 0, sw2 = 0;
        label nZero = 0;
        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }
            if (kf[f] == 0) { ++nZero; }
            const scalar e = mag(kf[f] - kappaExact);
            s1 += magSfIn[f]*e;
            s2 += magSfIn[f]*e*e;
            li = max(li, e);
            sw2 += mag(snI[f])*magSfIn[f]*e*e;
        }
        reduce(s1, sumOp<scalar>()); reduce(s2, sumOp<scalar>());
        reduce(li, maxOp<scalar>()); reduce(sw2, sumOp<scalar>());
        reduce(nZero, sumOp<label>());
        faceRows.append
        ({
            model, fp,
            (sumAf > 0) ? s1/sumAf : 0,
            (sumAf > 0) ? Foam::sqrt(s2/sumAf) : 0,
            li,
            (sumWf > 0) ? Foam::sqrt(sw2/sumWf) : 0,
            nZero
        });
    };

    // Contour-referenced models: kappa_f is the curvature of the level contour
    // through the face centre -> the stabilized foot-point correction applies.
    addFaceRow("quadraticCellCentre", 0, kfQ);
    addFaceRow("quadraticCellCentre", 1, footCorrect(kfQ)());
    addFaceRow("trHessian", 0, kfLap);
    addFaceRow("trHessian", 1, footCorrect(kfLap)());
    addFaceRow("fvmDivGradPsi", 0, kfFvmDivPsi);
    addFaceRow("fvmDivGradPsi", 1, footCorrect(kfFvmDivPsi)());
    addFaceRow("fvmTraceGradGradPsi", 0, kfFvmTracePsi);
    addFaceRow("fvmTraceGradGradPsi", 1, footCorrect(kfFvmTracePsi)());
    addFaceRow("fvmDivGradAlpha", 0, kfFvmDivAlpha);
    addFaceRow("fvmDivGradAlpha", 1, footCorrect(kfFvmDivAlpha)());

    // CONTROL: the 2D scalar inverse without the Gaussian term -- identical to
    // the corrected row in 2D (K = +0 exactly), first-order-wrong in 3D
    // (2/(R - d) on the sphere). The 3D gate asserts this failure.
    addFaceRow("scalarInverse2D", 1, footCorrectScalar2D(kfQ)());

    // Interface-referenced models: the evaluation offset is ~0 by the model's
    // own construction, so the parallel-curve correction does not apply
    // (FOOT_POINT = 0 rows only).
    addFaceRow("kangQuadratic", 0, kfKang);
    addFaceRow("quadraticNewtonFoot", 0, kfExt);
    addFaceRow("quadraticClosestPoint", 0, kfCP);
    if (is2D)
    {
        addFaceRow("footPointHeightFunction", 0, kfFoot);
        addFaceRow("isoConormal", 0, kfIso);
        addFaceRow("isoParabola", 0, kfParab);
    }
    if (haveConnected)
    {
        addFaceRow("connectedInterface", 0, kappaConnFace.primitiveField());
    }
    addFaceRow("interfaceMean", 0, kfMean);

    // The fully foot-point-native delivery (defined only WITH the algorithm).
    addFaceRow("stableFootPoint", 1, kfStable);

    // One interface curvature per CUT CELL, assigned to its active faces.
    addFaceRow("cutCellInverse", 1, kfCutCell);
    addFaceRow("cellMeanInverse", 1, kfCellMean);

    Info<< "face-centered curvature, " << nActive << " active faces"
        << " (foot-point unset: " << nFootFail
        << ", Gaussian fallbacks: " << nGaussFallback
        << ", stableFootPoint unset: " << nStableUnset
        << "; cut cells: " << nCutCells
        << ", active faces with no cut side: " << nNoCutSide
        << ", cut-cell foot failures: " << nCutFootFail << "):" << nl;
    for (const auto& r : faceRows)
    {
        Info<< "  " << r.model << (r.footPoint ? " +footPoint" : "")
            << "  L1/L2/Linf = " << r.L1 << " / " << r.L2 << " / " << r.Linf
            << "  (forceW L2 = " << r.L2w
            << ", zero faces = " << r.nZero << ")" << nl;
    }
    Info<< endl;

    // Representative cell size h (informational; the study uses N_CELLS):
    // sqrt(area/n) on pseudo-2D meshes, cbrt(volume/n) in 3D -- getting this
    // wrong corrupts every order fitted against DELTA_X by a factor 1.5.
    const boundBox bb(mesh.points());
    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar Lz = bb.max().z() - bb.min().z();
    const scalar dx =
        (mesh.nCells() > 0)
      ? (is2D
         ? Foam::sqrt(Lx*Ly/mesh.nCells())
         : Foam::cbrt(Lx*Ly*Lz/mesh.nCells()))
      : 0;

    Info<< "band cells                    : " << label(nBand) << nl
        << "kappa_exact                   : " << kappaExact << nl
        << "div(n) extended  L1/L2/Linf   : " << L1d << " / " << L2d << " / " << linfd
        << "   (rel L2 = " << L2d/kappaExact << ")" << nl
        << "div(n) NO extension L1/L2/Linf: " << L1n << " / " << L2n << " / " << linfn
        << "   (rel L2 = " << L2n/kappaExact << ")" << nl
        << "face-avg extended  L1/L2/Linf : " << L1f << " / " << L2f << " / " << linff
        << "   (rel L2 = " << L2f/kappaExact << ")" << nl
        << "tr(H)             L1/L2/Linf  : " << L1l << " / " << L2l << " / " << linfl
        << "   (rel L2 = " << L2l/kappaExact << ")" << nl
        << "foot-point HF     L1/L2/Linf  : " << L1p << " / " << L2p << " / " << linfp
        << "   (rel L2 = " << L2p/kappaExact << ", fallbacks = " << nFootFallback
        << ")" << nl
        << "iso-CONORMAL     L1/L2/Linf   : " << L1i << " / " << L2i << " / " << linfi
        << "   (rel L2 = " << L2i/kappaExact << ")" << nl
        << "iso-PARABOLA     L1/L2/Linf   : " << L1b << " / " << L2b << " / " << linfb
        << "   (rel L2 = " << L2b/kappaExact << ")" << nl << endl;

    if (Pstream::master())
    {
        OFstream os("leiaTestMeanCurvature.csv");
        os.precision(10);
        os << "TIME,DELTA_X,N_BAND,KAPPA_EXACT,"
              "E_L1_DIV,E_L2_DIV,E_LINF_DIV,E_L1_LAP,E_L2_LAP,E_LINF_LAP,"
              "E_L1_NOEXT,E_L2_NOEXT,E_LINF_NOEXT,"
              "E_L1_FACEAVG,E_L2_FACEAVG,E_LINF_FACEAVG,"
              "E_L1_FOOT,E_L2_FOOT,E_LINF_FOOT,N_FOOT_FALLBACK,"
              "E_L1_ISO,E_L2_ISO,E_LINF_ISO,E_L1_PARAB,E_L2_PARAB,E_LINF_PARAB" << nl;
        os << runTime.value() << ',' << dx << ',' << label(nBand) << ','
           << kappaExact << ',' << L1d << ',' << L2d << ',' << linfd << ','
           << L1l << ',' << L2l << ',' << linfl << ','
           << L1n << ',' << L2n << ',' << linfn << ','
           << L1f << ',' << L2f << ',' << linff << ',';
        // 2D-only models: sentinel-empty cells in 3D (their fields hold the
        // fallback copy of kappaNoExt there -- writing those numbers would
        // silently mislabel them).
        if (is2D)
        {
            os << L1p << ',' << L2p << ',' << linfp << ',' << nFootFallback << ','
               << L1i << ',' << L2i << ',' << linfi << ','
               << L1b << ',' << L2b << ',' << linfb << nl;
        }
        else
        {
            os << ",,," << nFootFallback << ",,,,,," << nl;
        }
    }

    // Face-centered curvature errors, tidy (one row per MODEL x FOOT_POINT),
    // for make_face_curvature_fig.py.
    if (Pstream::master())
    {
        OFstream osF("leiaTestFaceCurvature.csv");
        osF.precision(10);
        osF << "MODEL,FOOT_POINT,N_ACTIVE_FACES,DELTA_X,KAPPA_EXACT,"
               "E_L1,E_L2,E_LINF,E_L2_FORCEW,N_ZERO_FACES,N_FOOT_FAIL" << nl;
        for (const auto& r : faceRows)
        {
            osF << r.model << ',' << r.footPoint << ',' << nActive << ','
                << dx << ',' << kappaExact << ','
                << r.L1 << ',' << r.L2 << ',' << r.Linf << ','
                << r.L2w << ',' << r.nZero << ',' << nFootFail << nl;
        }
    }

    runTime.writeNow();   // write kappaDiv/kappaNoExt/kappaFaceAvg/kappaLap fields

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
