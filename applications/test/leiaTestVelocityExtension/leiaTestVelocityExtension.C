/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestVelocityExtension

Description
    Static (t = 0) verification of the runtime-selected velocityExtension
    model: initialize the case fields (psi/alpha from leiaSetFields, U/phi
    from the prescribed velocityModel at t = 0), apply the extension ONCE (no
    advection) and measure whether the extended velocity is constant along
    the discrete level-set normal,

        e = | nHat . grad(Uext) |,   nHat = grad(psi)/|grad(psi)|,

    computed with fvc::grad (the models' own discretization). Writes the
    fields eNormalUext and eNormalU (raw-velocity reference) at time 0 and a
    CSV `leiaTestVelocityExtension.csv` with volume-weighted L2 / Linf norms
    over the extension-active region |psi| <= nLayers*cellSize (inside the
    Uext -> U fade start), the raw-velocity reference norms, the ratio, and
    the L2 split by sign(psi) (one-sidedness diagnosis). If the extension
    does not converge on a static interface, it never will under advection.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "narrowBand.H"
#include "phaseIndicator.H"
#include "velocityModel.H"
#include "velocityExtension.H"
#include "levelSetImplicitSurfaces.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // ------------------------------------------------------------------ //
    // Field/model initialization: identical order to the solver's
    // createFields.H so the extension sees exactly the t = 0 configuration.
    // ------------------------------------------------------------------ //

    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phi", dimVolume/dimTime, 0)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimVelocity, vector(0, 0, 0))
    );

    // ------------------------------------------------------------------ //
    // The input field (step A proper, 2026-11-14). velocityExtension.testInput:
    //   distance        psi = d (as written by leiaSetFields; the default)
    //   scaledDistance  psi = d g(x), g = exp(ln(1.25) sin(2 pi x/L) sin(2 pi y/L)),
    //                   L = testInputWavelength (default 0.5): q = |grad psi|
    //                   between 0.8 and 1.25 near the interface, smooth
    //   quadratic       psi = d (d + 2R) (= r^2 - R^2 for implicitSphere only):
    //                   q = 2r, about 2R at the interface
    // All inputs have the SAME zero set, so the exact extension u(Pi(x)) and
    // alpha are unchanged; only what an operator that reads psi as a distance
    // sees changes. d is the analytic signed distance of levelSet.implicitSurface.
    // ------------------------------------------------------------------ //
    {
        const fvSolution& fvSolIn(mesh);
        const dictionary& lsIn = fvSolIn.subDict("levelSet");
        const dictionary veIn(lsIn.subOrEmptyDict("velocityExtension"));
        const word input = veIn.getOrDefault<word>("testInput", "distance");
        if (input != "distance")
        {
            const dictionary& sIn = lsIn.subDict("implicitSurface");
            autoPtr<implicitSurface> surfIn =
                implicitSurface::New(sIn.get<word>("type"), sIn);
            const scalar lam = veIn.getOrDefault<scalar>("testInputWavelength", 0.5);
            scalar R = 0;
            if (input == "quadratic")
            {
                if (sIn.get<word>("type") != "implicitSphere")
                {
                    FatalErrorInFunction
                        << "testInput quadratic needs implicitSphere" << exit(FatalError);
                }
                R = sIn.get<scalar>("radius");
            }
            else if (input != "scaledDistance")
            {
                FatalErrorInFunction
                    << "testInput " << input << " is not one of"
                    << " (distance scaledDistance quadratic)" << exit(FatalError);
            }
            auto transform = [&](const point& x) -> scalar
            {
                const scalar d = surfIn->value(x);
                if (input == "quadratic")
                {
                    return d*(d + 2*R);
                }
                const scalar w = constant::mathematical::twoPi/lam;
                return d*Foam::exp
                (
                    Foam::log(1.25)*Foam::sin(w*x.x())*Foam::sin(w*x.y())
                );
            };
            forAll(psi, c)
            {
                psi[c] = transform(mesh.C()[c]);
            }
            forAll(psi.boundaryField(), patchI)
            {
                fvPatchScalarField& pf = psi.boundaryFieldRef()[patchI];
                if (pf.coupled()) continue;
                const vectorField& Cf = mesh.boundary()[patchI].Cf();
                forAll(pf, i)
                {
                    pf[i] = transform(Cf[i]);
                }
            }
            psi.correctBoundaryConditions();
            Info<< "testInput " << input << ": psi replaced (same zero set)" << endl;
        }
    }

    autoPtr<narrowBand> nband = narrowBand::New(mesh, psi);
    nband->calc();

    autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh);
    phaseInd->calcPhaseIndicator(alpha, psi);

    // Prescribed verification velocity at t = 0 (oscillation factor = 1).
    autoPtr<velocityModel> velModel = velocityModel::New(mesh);
    velModel->setVolumetricFlux(phi);
    velModel->setVelocity(U);

    // ------------------------------------------------------------------ //
    // ONE application of the runtime-selected extension model.
    // ------------------------------------------------------------------ //
    autoPtr<velocityExtension> velExt = velocityExtension::New(mesh);
    velExt->correct();

    // ------------------------------------------------------------------ //
    // Metric: e = |nHat . grad(Uext)| and the raw-velocity reference.
    // ------------------------------------------------------------------ //
    tmp<volVectorField> tGradPsi = fvc::grad(psi);
    const volVectorField& gradPsi = tGradPsi();
    volVectorField nHat
    (
        "eNormalNHat",
        gradPsi/(mag(gradPsi)
      + dimensionedScalar("e", gradPsi.dimensions(), SMALL))
    );

    volScalarField eNormalUext
    (
        IOobject
        (
            "eNormalUext",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(nHat & fvc::grad(velExt->Uext()))
    );

    volScalarField eNormalU
    (
        IOobject
        (
            "eNormalU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(nHat & fvc::grad(U))
    );

    eNormalUext.write();
    eNormalU.write();
    psi.write();   // for the psi = 0 overlay in the atlas

    // ------------------------------------------------------------------ //
    // Norms over the extension-active region |psi| <= nLayers*cellSize
    // (inside the Uext -> U fade start), volume-weighted, parallel-safe.
    // ------------------------------------------------------------------ //
    const fvSolution& fvSol(mesh);
    // subOrEmptyDict returns a temporary -> copy by value.
    const dictionary velExtDict
    (
        fvSol.subDict("levelSet").subOrEmptyDict("velocityExtension")
    );
    const label nLayers = velExtDict.getOrDefault<label>("nLayers", 3);

    // Per-cell characteristic size = min face-neighbour centre distance
    // (1/deltaCoeffs): the in-plane h for 2D meshes of arbitrary thickness
    // (cbrt(V) = (h^2 t)^(1/3) would mis-scale the measured region).
    scalarField cellSize(mesh.nCells(), GREAT);
    {
        const surfaceScalarField& dc = mesh.deltaCoeffs();
        const labelUList& fown = mesh.owner();
        const labelUList& fnei = mesh.neighbour();
        forAll(fown, f)
        {
            const scalar d = 1.0/dc[f];
            cellSize[fown[f]] = Foam::min(cellSize[fown[f]], d);
            cellSize[fnei[f]] = Foam::min(cellSize[fnei[f]], d);
        }
        forAll(mesh.boundary(), patchI)
        {
            if (mesh.boundary()[patchI].coupled())
            {
                const labelUList& fc = mesh.boundary()[patchI].faceCells();
                const fvsPatchScalarField& pdc = dc.boundaryField()[patchI];
                forAll(fc, i)
                {
                    cellSize[fc[i]] = Foam::min(cellSize[fc[i]], 1.0/pdc[i]);
                }
            }
        }
    }

    scalar sumE2V = 0, sumRaw2V = 0, sumV = 0;
    scalar sumE2Vin = 0, sumVin = 0, sumE2Vout = 0, sumVout = 0;
    scalar maxE = 0, maxRaw = 0;
    forAll(psi, c)
    {
        const scalar L = nLayers*cellSize[c];
        if (Foam::mag(psi[c]) > L)
        {
            continue;
        }
        const scalar V = mesh.V()[c];
        const scalar e2 = Foam::sqr(eNormalUext[c]);
        sumE2V += e2*V;
        sumRaw2V += Foam::sqr(eNormalU[c])*V;
        sumV += V;
        maxE = Foam::max(maxE, eNormalUext[c]);
        maxRaw = Foam::max(maxRaw, eNormalU[c]);
        if (psi[c] < 0)
        {
            sumE2Vin += e2*V; sumVin += V;
        }
        else
        {
            sumE2Vout += e2*V; sumVout += V;
        }
    }
    reduce(sumE2V, sumOp<scalar>());
    reduce(sumRaw2V, sumOp<scalar>());
    reduce(sumV, sumOp<scalar>());
    reduce(sumE2Vin, sumOp<scalar>());
    reduce(sumVin, sumOp<scalar>());
    reduce(sumE2Vout, sumOp<scalar>());
    reduce(sumVout, sumOp<scalar>());
    reduce(maxE, maxOp<scalar>());
    reduce(maxRaw, maxOp<scalar>());

    const scalar l2 = Foam::sqrt(sumE2V/Foam::max(sumV, SMALL));
    const scalar l2raw = Foam::sqrt(sumRaw2V/Foam::max(sumV, SMALL));
    const scalar l2in = Foam::sqrt(sumE2Vin/Foam::max(sumVin, SMALL));
    const scalar l2out = Foam::sqrt(sumE2Vout/Foam::max(sumVout, SMALL));
    const scalar h = Foam::max(Foam::pow(mesh.deltaCoeffs(), -1)).value();

    Info<< "Static extension verification (|psi| <= " << nLayers
        << "*cellSize):" << nl
        << "  E_NORMAL_L2   = " << l2 << "  (raw " << l2raw << ")" << nl
        << "  E_NORMAL_LINF = " << maxE << "  (raw " << maxRaw << ")" << nl
        << "  RATIO_L2      = " << l2/Foam::max(l2raw, SMALL) << nl
        << "  L2 in/out     = " << l2in << " / " << l2out << endl;

    if (Pstream::master())
    {
        OFstream csv("leiaTestVelocityExtension.csv");
        csv << "DELTA_X,E_NORMAL_L2,E_NORMAL_LINF,E_NORMAL_RAW_L2,"
            << "E_NORMAL_RAW_LINF,RATIO_L2,E_NORMAL_L2_IN,E_NORMAL_L2_OUT\n";
        csv << h << "," << l2 << "," << maxE << "," << l2raw << ","
            << maxRaw << "," << l2/Foam::max(l2raw, SMALL) << ","
            << l2in << "," << l2out << "\n";
    }

    // ------------------------------------------------------------------ //
    // THE EXACT REFERENCE (2026-11-12, docs/plan-velocity-extension-
    // verification.md, step A). The metric above, |nHat . grad Uext|, is zero
    // for ANY field that is constant along the normals, a wrong one included.
    // Here the computed extension is compared with the exact one,
    //     ue*(x) = u(Pi(x)),   Pi(x) = x - d(x) grad d/|grad d|,
    // where d is the analytic signed distance of levelSet.implicitSurface
    // (implicitSphere, implicitPlane with gradient 1, signedDistanceEllipse,
    // signedDistanceEllipsoid: exact signed distances) and u the prescribed
    // velocity model evaluated AT THE FOOT POINT. Written to a separate CSV,
    // so the CSV above is unchanged. Norms: volume-weighted L1 and L2 (never
    // L_inf), over the cells with |d| <= k h (k = 1..4, h the cell size) and
    // over the fixed physical band |d| <= W (exactBandWidth, default 4 h_max);
    // the face-flux error over the faces with |d(x_f)| <= W, which is what the
    // transport sees; and the same norms for the raw velocity U as the
    // "no extension" reference.
    // ------------------------------------------------------------------ //
    const dictionary& lsDict = fvSol.subDict("levelSet");
    if (lsDict.found("implicitSurface"))
    {
        const dictionary& sDict = lsDict.subDict("implicitSurface");
        autoPtr<implicitSurface> surf =
            implicitSurface::New(sDict.get<word>("type"), sDict);
        const volVectorField& C = mesh.C();
        const volVectorField& Ue = velExt->Uext();
        const scalar hMax = gMax(cellSize);
        const scalar W =
            velExtDict.getOrDefault<scalar>("exactBandWidth", 4*hMax);

        auto footVelocity = [&](const vector& x, scalar& d) -> vector
        {
            d = surf->value(x);
            vector g = surf->grad(x);
            g /= Foam::max(Foam::mag(g), SMALL);
            return velModel->velocity(x - d*g);
        };

        // Shells |d| <= 1h..8h (8h since 2026-11-13: the fade of the default
        // nLayers 3 reaches 2L = 8h; plan section 4a, repair 2).
        const label nK = 8;
        scalarField e1(nK + 1, 0), e2(nK + 1, 0), r1(nK + 1, 0), r2(nK + 1, 0), vv(nK + 1, 0);
        forAll(C, c)
        {
            scalar d = 0;
            const vector ueStar = footVelocity(C[c], d);
            const scalar e = Foam::mag(Ue[c] - ueStar);
            const scalar r = Foam::mag(U[c] - ueStar);
            const scalar V = mesh.V()[c];
            for (label k = 0; k <= nK; ++k)
            {
                const bool in = (k < nK)
                    ? (Foam::mag(d) <= (k + 1)*cellSize[c])
                    : (Foam::mag(d) <= W);
                if (!in) continue;
                e1[k] += e*V; e2[k] += e*e*V; r1[k] += r*V; r2[k] += r*r*V; vv[k] += V;
            }
        }
        // face flux: internal faces only (a coupled face is counted by both
        // ranks; the internal faces carry the comparison in serial and np 4)
        const surfaceScalarField& phiE = velExt->phi();
        const surfaceVectorField& Cf = mesh.Cf();
        const surfaceScalarField& magSf = mesh.magSf();
        scalar f1 = 0, f2 = 0, fr1 = 0, fr2 = 0, fa = 0;
        // The flux error also per shell |d(x_f)| <= k h (k = 1..8), so that an
        // interpolation error inside the band can be told from a jump at the
        // band edge (2026-11-12: the first ladder could not tell them apart).
        scalarField g1(nK, 0), g2(nK, 0), ga(nK, 0);
        // The flux error against the reference WITH the intended fade:
        //     ue_ref(x) = w(a) u(Pi(x)) + (1 - w(a)) u(x),
        //     w(a) = 1 (a <= 0), (1 + cos(pi a))/2 (0 < a < 1), 0 (a >= 1),
        // a = (|d| - L)/L (fadePlacement legacy, extendThroughFade) or
        // a = (|d| - L/2)/(L/2) (insideBand), L = (nLayers + 1) h, over the
        // faces with |d(x_f)| <= 3L. This is the operator the interfaceExtension
        // base intends in fadeMode levelSet on a signed-distance input; a
        // velocity jump at L is a first-order error against it. Not meaningful
        // for fadeMode layer (its fade is per layer) or for a model outside
        // the interfaceExtension base.
        const label nLay = velExtDict.getOrDefault<label>("nLayers", 3);
        const word placement =
            velExtDict.getOrDefault<word>("fadePlacement", "legacy");
        scalar q1 = 0, q2 = 0, qa = 0;
        const labelUList& fOwn = mesh.owner();
        forAll(phiE, f)
        {
            scalar d = 0;
            const vector ueStar = footVelocity(Cf[f], d);
            const scalar A = magSf[f];
            const scalar ref = ueStar & mesh.Sf()[f];
            const scalar e = Foam::mag(phiE[f] - ref)/A;
            for (label k = 0; k < nK; ++k)
            {
                if (Foam::mag(d) <= (k + 1)*cellSize[fOwn[f]])
                {
                    g1[k] += e*A; g2[k] += e*e*A; ga[k] += A;
                }
            }
            {
                const scalar hf = cellSize[fOwn[f]];
                const scalar L = (nLay + 1)*hf;
                if (Foam::mag(d) <= 3*L)
                {
                    const scalar a = (placement == "insideBand")
                        ? (Foam::mag(d) - 0.5*L)/(0.5*L)
                        : (Foam::mag(d) - L)/L;
                    const scalar w = (a <= 0) ? 1.0
                        : (a >= 1) ? 0.0 : 0.5*(1.0 + Foam::cos(M_PI*a));
                    const vector uf = velModel->velocity(Cf[f]);
                    const scalar refF =
                        (w*ueStar + (1.0 - w)*uf) & mesh.Sf()[f];
                    const scalar q = Foam::mag(phiE[f] - refF)/A;
                    q1 += q*A; q2 += q*q*A; qa += A;
                }
            }
            if (Foam::mag(d) > W) continue;
            const scalar r = Foam::mag(phi[f] - ref)/A;
            f1 += e*A; f2 += e*e*A; fr1 += r*A; fr2 += r*r*A; fa += A;
        }
        reduce(g1, sumOp<scalarField>()); reduce(g2, sumOp<scalarField>()); reduce(ga, sumOp<scalarField>());
        reduce(e1, sumOp<scalarField>()); reduce(e2, sumOp<scalarField>());
        reduce(r1, sumOp<scalarField>()); reduce(r2, sumOp<scalarField>());
        reduce(vv, sumOp<scalarField>());
        reduce(f1, sumOp<scalar>()); reduce(f2, sumOp<scalar>());
        reduce(fr1, sumOp<scalar>()); reduce(fr2, sumOp<scalar>()); reduce(fa, sumOp<scalar>());
        reduce(q1, sumOp<scalar>()); reduce(q2, sumOp<scalar>()); reduce(qa, sumOp<scalar>());

        Info<< "Exact extension reference ue* = u(Pi(x)) (" << sDict.get<word>("type")
            << "), W = " << W << ":" << nl;
        for (label k = 0; k <= nK; ++k)
        {
            const scalar V = Foam::max(vv[k], SMALL);
            Info<< "  " << (k < nK ? word("|d|<=" + Foam::name(k + 1) + "h") : word("|d|<=W"))
                << "  L1 " << e1[k]/V << "  L2 " << Foam::sqrt(e2[k]/V)
                << "   (raw U: L1 " << r1[k]/V << "  L2 " << Foam::sqrt(r2[k]/V) << ")" << nl;
        }
        const scalar FA = Foam::max(fa, SMALL);
        Info<< "  face flux |d|<=W  L1 " << f1/FA << "  L2 " << Foam::sqrt(f2/FA)
            << "   (base flux: L1 " << fr1/FA << "  L2 " << Foam::sqrt(fr2/FA) << ")" << endl;

        if (Pstream::master())
        {
            OFstream csv("leiaTestVelocityExtensionExact.csv");
            csv << "DELTA_X,W";
            for (label k = 0; k <= nK; ++k)
            {
                const word t = (k < nK) ? word(Foam::name(k + 1) + "H") : word("W");
                csv << ",E_EXT_L1_" << t << ",E_EXT_L2_" << t
                    << ",E_RAW_L1_" << t << ",E_RAW_L2_" << t;
            }
            csv << ",E_FLUX_L1_W,E_FLUX_L2_W,E_FLUX_RAW_L1_W,E_FLUX_RAW_L2_W";
            for (label k = 0; k < nK; ++k)
            {
                csv << ",E_FLUX_L1_" << (k + 1) << "H,E_FLUX_L2_" << (k + 1) << "H";
            }
            csv << ",E_FADEFLUX_L1,E_FADEFLUX_L2";
            csv << "\n";
            csv << h << "," << W;
            for (label k = 0; k <= nK; ++k)
            {
                const scalar V = Foam::max(vv[k], SMALL);
                csv << "," << e1[k]/V << "," << Foam::sqrt(e2[k]/V)
                    << "," << r1[k]/V << "," << Foam::sqrt(r2[k]/V);
            }
            csv << "," << f1/FA << "," << Foam::sqrt(f2/FA)
                << "," << fr1/FA << "," << Foam::sqrt(fr2/FA);
            for (label k = 0; k < nK; ++k)
            {
                const scalar GA = Foam::max(ga[k], SMALL);
                csv << "," << g1[k]/GA << "," << Foam::sqrt(g2[k]/GA);
            }
            const scalar QA = Foam::max(qa, SMALL);
            csv << "," << q1/QA << "," << Foam::sqrt(q2/QA);
            csv << "\n";
        }

        // -------------------------------------------------------------- //
        // STEP A ROUND 2 (2026-11-18, plan section 5, "the model's own
        // extension"). Two additions, written to a THIRD CSV
        // (leiaTestVelocityExtensionLayers.csv) so the two CSVs above are
        // unchanged:
        //  1. the reference at the Newton point of the EXACT input field,
        //         Y(x) = x - psi grad(psi)/|grad(psi)|^2,
        //     ue_Y = u(Y), next to u(Pi(x)). The theorems of the Bothe-Maric-
        //     Soga manuscript hold for Y; Y = Pi only for a signed distance.
        //     psi is the analytic input of testInput (distance,
        //     scaledDistance, quadratic) with its analytic gradient.
        //     Also the blend reference of projectionBlend:
        //         ue_B = delta u(x) + (1 - delta) u(Y),
        //     delta(a) from w1/w2, a = |psi|/(|grad psi| h): the operator
        //     that projectionBlend intends (meaningful only for that model);
        //  2. every norm PER TOPOLOGICAL LAYER l = 0..8 (l = 0: the cells
        //     with a sign change of psi across a face, internal or coupled;
        //     l + 1: a face neighbour of layer l), and over l = 0..3, the
        //     layers that the compact face curvature reads (plan section 4c,
        //     A0.4). A distance shell |d| <= k h mixes layers; the curvature
        //     requirement is stated per layer. Face flux per face layer
        //     min(l_own, l_nei), internal faces only (as above).
        // -------------------------------------------------------------- //
        {
            const word input =
                velExtDict.getOrDefault<word>("testInput", "distance");
            const scalar lam =
                velExtDict.getOrDefault<scalar>("testInputWavelength", 0.5);
            const scalar Rq =
                (input == "quadratic") ? sDict.get<scalar>("radius") : 0;
            const scalar pbW1 = velExtDict.getOrDefault<scalar>("w1", 1.5);
            const scalar pbW2 = velExtDict.getOrDefault<scalar>("w2", 4.0);

            // Newton point of the exact input; dN = psi/|grad psi| (signed).
            auto newtonY = [&](const vector& x, scalar& dN) -> vector
            {
                const scalar d = surf->value(x);
                vector n = surf->grad(x);
                n /= Foam::max(Foam::mag(n), SMALL);
                scalar p = d;
                vector g = n;
                if (input == "quadratic")
                {
                    p = d*(d + 2*Rq);
                    g = 2*(d + Rq)*n;
                }
                else if (input == "scaledDistance")
                {
                    const scalar w = constant::mathematical::twoPi/lam;
                    const scalar a = Foam::log(1.25);
                    const scalar sx = Foam::sin(w*x.x()), cx = Foam::cos(w*x.x());
                    const scalar sy = Foam::sin(w*x.y()), cy = Foam::cos(w*x.y());
                    const scalar gg = Foam::exp(a*sx*sy);
                    const vector dg = gg*a*w*vector(cx*sy, sx*cy, 0);
                    p = d*gg;
                    g = gg*n + d*dg;
                }
                const scalar q2 = Foam::max(Foam::magSqr(g), VSMALL);
                dN = p/Foam::sqrt(q2);
                return x - p*g/q2;
            };

            // Topological layer, parallel: sign change across internal and
            // coupled faces, then face-neighbour passes with a swap of the
            // layer index over coupled patches.
            const label nL = 8;
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            labelList layer(mesh.nCells(), labelMax);
            forAll(nei, f)
            {
                if ((psi[own[f]] < 0) != (psi[nei[f]] < 0))
                {
                    layer[own[f]] = 0;
                    layer[nei[f]] = 0;
                }
            }
            {
                scalarField psiNb;
                syncTools::swapBoundaryCellList(mesh, psi.primitiveField(), psiNb);
                forAll(mesh.boundaryMesh(), patchI)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchI];
                    if (!pp.coupled()) continue;
                    forAll(pp, i)
                    {
                        const label c = pp.faceCells()[i];
                        const label bf = pp.start() + i - mesh.nInternalFaces();
                        if ((psi[c] < 0) != (psiNb[bf] < 0)) layer[c] = 0;
                    }
                }
            }
            for (label pass = 1; pass <= nL; ++pass)
            {
                labelList next(layer);
                forAll(nei, f)
                {
                    if (layer[own[f]] == pass - 1 && next[nei[f]] == labelMax)
                        next[nei[f]] = pass;
                    if (layer[nei[f]] == pass - 1 && next[own[f]] == labelMax)
                        next[own[f]] = pass;
                }
                labelList layerNb;
                syncTools::swapBoundaryCellList(mesh, layer, layerNb);
                forAll(mesh.boundaryMesh(), patchI)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchI];
                    if (!pp.coupled()) continue;
                    forAll(pp, i)
                    {
                        const label c = pp.faceCells()[i];
                        const label bf = pp.start() + i - mesh.nInternalFaces();
                        if (layerNb[bf] == pass - 1 && next[c] == labelMax)
                            next[c] = pass;
                    }
                }
                layer.transfer(next);
            }

            // Per layer l = 0..nL, and slot nL + 1 = layers 0..3.
            const label nS = nL + 2;
            scalarField nC(nS, 0), vol(nS, 0), x1(nS, 0), x2(nS, 0),
                y1(nS, 0), y2(nS, 0), b1(nS, 0), b2(nS, 0), ry1(nS, 0),
                ry2(nS, 0), yp1(nS, 0);
            forAll(C, c)
            {
                const label l = layer[c];
                if (l > nL) continue;
                scalar d = 0, dN = 0;
                const vector uPi = footVelocity(C[c], d);
                const vector Y = newtonY(C[c], dN);
                const vector uY = velModel->velocity(Y);
                const scalar hc = cellSize[c];
                const scalar a = Foam::mag(dN)/hc;
                scalar delta = 1;
                if (a <= pbW1) delta = 0;
                else if (a < pbW2)
                {
                    const scalar s = (a - pbW1)/(pbW2 - pbW1);
                    delta = 0.5*(1 - Foam::cos(constant::mathematical::pi*s));
                }
                const vector uB = delta*U[c] + (1 - delta)*uY;
                vector nc = surf->grad(C[c]);
                nc /= Foam::max(Foam::mag(nc), SMALL);
                const vector xPi = C[c] - d*nc;
                const scalar ex = Foam::mag(Ue[c] - uPi);
                const scalar ey = Foam::mag(Ue[c] - uY);
                const scalar eb = Foam::mag(Ue[c] - uB);
                const scalar er = Foam::mag(U[c] - uY);
                const scalar V = mesh.V()[c];
                const label slots[2] = {l, (l <= 3 ? nL + 1 : label(-1))};
                for (const label s : slots)
                {
                    if (s < 0) continue;
                    nC[s] += 1; vol[s] += V;
                    x1[s] += ex*V; x2[s] += ex*ex*V;
                    y1[s] += ey*V; y2[s] += ey*ey*V;
                    b1[s] += eb*V; b2[s] += eb*eb*V;
                    ry1[s] += er*V; ry2[s] += er*er*V;
                    yp1[s] += Foam::mag(Y - xPi)/hc*V;
                }
            }
            scalarField fy1(nS, 0), fy2(nS, 0), fx1(nS, 0), fa(nS, 0);
            forAll(nei, f)
            {
                const label l = Foam::min(layer[own[f]], layer[nei[f]]);
                if (l > nL) continue;
                scalar d = 0, dN = 0;
                const vector uPi = footVelocity(Cf[f], d);
                const vector uY = velModel->velocity(newtonY(Cf[f], dN));
                const scalar A = magSf[f];
                const scalar ey = Foam::mag(phiE[f] - (uY & mesh.Sf()[f]))/A;
                const scalar ex = Foam::mag(phiE[f] - (uPi & mesh.Sf()[f]))/A;
                const label slots[2] = {l, (l <= 3 ? nL + 1 : label(-1))};
                for (const label s : slots)
                {
                    if (s < 0) continue;
                    fy1[s] += ey*A; fy2[s] += ey*ey*A; fx1[s] += ex*A; fa[s] += A;
                }
            }
            for
            (
                scalarField* fld :
                {&nC, &vol, &x1, &x2, &y1, &y2, &b1, &b2, &ry1, &ry2, &yp1,
                 &fy1, &fy2, &fx1, &fa}
            )
            {
                reduce(*fld, sumOp<scalarField>());
            }

            Info<< "Per-layer errors (round 2): against u(Pi), u(Y) of the exact"
                << " input (" << input << "), the projectionBlend reference (w1 "
                << pbW1 << ", w2 " << pbW2 << "):" << nl;
            for (label s = 0; s < nS; ++s)
            {
                const scalar V = Foam::max(vol[s], SMALL);
                const scalar A = Foam::max(fa[s], SMALL);
                Info<< "  " << (s <= nL ? word("l=" + Foam::name(s)) : word("l<=3"))
                    << "  n " << label(nC[s])
                    << "  Pi L1 " << x1[s]/V << "  Y L1 " << y1[s]/V
                    << "  B L1 " << b1[s]/V << "  flux(Y) L1 " << fy1[s]/A
                    << "  raw(Y) L1 " << ry1[s]/V
                    << "  |Y-Pi|/h " << yp1[s]/V << nl;
            }
            Info<< endl;

            if (Pstream::master())
            {
                OFstream csv("leiaTestVelocityExtensionLayers.csv");
                csv << "DELTA_X";
                for (label s = 0; s < nS; ++s)
                {
                    const word t =
                        (s <= nL) ? word("L" + Foam::name(s)) : word("L0TO3");
                    csv << ",N_" << t
                        << ",E_EXT_L1_" << t << ",E_EXT_L2_" << t
                        << ",E_EXTY_L1_" << t << ",E_EXTY_L2_" << t
                        << ",E_EXTB_L1_" << t << ",E_EXTB_L2_" << t
                        << ",E_RAWY_L1_" << t << ",E_RAWY_L2_" << t
                        << ",D_YPI_L1_" << t
                        << ",E_FLUXY_L1_" << t << ",E_FLUXY_L2_" << t
                        << ",E_FLUX_L1_" << t;
                }
                csv << "\n" << h;
                for (label s = 0; s < nS; ++s)
                {
                    const scalar V = Foam::max(vol[s], SMALL);
                    const scalar A = Foam::max(fa[s], SMALL);
                    csv << "," << label(nC[s])
                        << "," << x1[s]/V << "," << Foam::sqrt(x2[s]/V)
                        << "," << y1[s]/V << "," << Foam::sqrt(y2[s]/V)
                        << "," << b1[s]/V << "," << Foam::sqrt(b2[s]/V)
                        << "," << ry1[s]/V << "," << Foam::sqrt(ry2[s]/V)
                        << "," << yp1[s]/V
                        << "," << fy1[s]/A << "," << Foam::sqrt(fy2[s]/A)
                        << "," << fx1[s]/A;
                }
                csv << "\n";
            }
        }
    }
    else
    {
        Info<< "No levelSet.implicitSurface: the exact extension reference is skipped." << endl;
    }

    Info<< nl;
    runTime.printExecutionTime(Info);
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
