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

        const label nK = 4;
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
        // The flux error also per shell |d(x_f)| <= k h (k = 1..4), so that an
        // interpolation error inside the band can be told from a jump at the
        // band edge (2026-11-12: the first ladder could not tell them apart).
        scalarField g1(nK, 0), g2(nK, 0), ga(nK, 0);
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
            csv << "\n";
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
