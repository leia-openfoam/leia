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

    Info<< nl;
    runTime.printExecutionTime(Info);
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
