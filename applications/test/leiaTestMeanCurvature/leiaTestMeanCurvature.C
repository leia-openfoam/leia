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
    const label nFootFallback = computeFootPointCurvature
    (
        mesh, psi,
        mesh.lookupObject<volScalarField>("NarrowBand"),
        kappaFoot
    );

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
    vortDriver(kappaFoot,  Vp2, Vpi);
    vortDriver(kappaIso,   Vi2, Vii);
    vortDriver(kappaParab, Vb2, Vbi);
    Info<< "vortical driver |grad(kappa) x grad(alpha)| (band L2 / Linf) [1/m^3]:" << nl
        << "  div(n) extended : " << Vd2 << " / " << Vdi << nl
        << "  foot-point HF   : " << Vp2 << " / " << Vpi << nl
        << "  iso-CONORMAL    : " << Vi2 << " / " << Vii << nl
        << "  iso-PARABOLA    : " << Vb2 << " / " << Vbi << nl << endl;

    // Representative cell size h for a 2D square domain (informational; the study
    // uses N_CELLS): h = sqrt(area / nCells2D).
    const boundBox bb(mesh.points());
    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar dx = (mesh.nCells() > 0) ? Foam::sqrt(Lx*Ly/mesh.nCells()) : 0;

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
           << L1f << ',' << L2f << ',' << linff << ','
           << L1p << ',' << L2p << ',' << linfp << ',' << nFootFallback << ','
           << L1i << ',' << L2i << ',' << linfi << ','
           << L1b << ',' << L2b << ',' << linfb << nl;
    }

    runTime.writeNow();   // write kappaDiv/kappaNoExt/kappaFaceAvg/kappaLap fields

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
