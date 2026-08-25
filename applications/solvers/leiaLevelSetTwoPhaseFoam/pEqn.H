{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    // Diagnostic pressure-compatibility oracle (default OFF): replace the
    // physical face-varying momentum coefficient by the volume-weighted mean
    // of the cell rAU field. This isolates variable-coefficient rAUf weighting
    // from the paired snGrad/Laplacian operator and fvc::reconstruct. The cell
    // rAU used for velocity recovery remains physical. Production runs never
    // set this environment variable and retain the original path byte-for-byte.
    if (!Foam::getEnv("LEIA_DIAGNOSTIC_CONSTANT_RAUF").empty())
    {
        const scalarField& cellVolumes = mesh.V();
        const scalarField& cellRAU = rAU().primitiveField();
        scalar volumeSum = 0.0;
        scalar volumeRAUSum = 0.0;
        forAll(cellVolumes, celli)
        {
            volumeSum += cellVolumes[celli];
            volumeRAUSum += cellVolumes[celli]*cellRAU[celli];
        }
        reduce(volumeSum, sumOp<scalar>());
        reduce(volumeRAUSum, sumOp<scalar>());
        const dimensionedScalar meanRAU
        (
            "meanRAU",
            rAU().dimensions(),
            volumeRAUSum/max(volumeSum, VSMALL)
        );
        rAUf = meanRAU;
        Info<< "rAUf diagnostic: constant volume-weighted mean = "
            << meanRAU.value() << nl;
    }
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    // Balanced-force invariant: capillarity and pressure live first in the
    // same integrated scalar face-flux space. At equilibrium the pressure-
    // equation flux cancels G_sigma face-by-face; only their difference is
    // passed to the common fvc::reconstruct below.
    surfaceScalarField phig
    (
        (GSigma - ghf*fvc::snGrad(rho)*mesh.magSf())*rAUf
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // WP8.0 diagnostic: the pressure-solver performance of the LAST solve, so
    // the capillary flux residual below can be split into the part that is not
    // absorbable in principle (non-gradient component of G_sigma) and the part
    // that is merely not yet absorbed (linear-solver residual).
    scalar pInitialResidual = 0, pFinalResidual = 0;
    label pIterations = 0;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        {
            const auto pPerf =
                p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
            pInitialResidual = pPerf.initialResidual();
            pFinalResidual = pPerf.finalResidual();
            pIterations = pPerf.nIterations();
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            // Secondary capillary/pressure-coupling diagnostic using the ACTUAL
            // PIMPLE weights and reconstruction. This deliberately remains
            // downstream of the direct geometrical force-error test: it answers
            // how much velocity increment survives the rAU/rAUf pressure
            // correction, not whether the surface-tension force itself is
            // geometrically accurate.
            //
            // For a uniformly translating equilibrium, all non-capillary
            // predictor fluxes should already be divergence-free. The expression
            // below therefore isolates the useful coupled residual without moving
            // the equations to a translating reference frame:
            //   U_sigma*   = rAU*reconstruct(phig/rAUf)
            //   U_sigma,res = rAU*reconstruct((phig-pFlux)/rAUf).
            // There is deliberately no cell-source contribution: every
            // production model is constrained to the G_sigma face-flux API.
            if (pimple.finalInnerIter() && pimple.finalIter())
            {
                const volVectorField USigmaPredictor
                (
                    "USigmaPredictor",
                    rAU()*fvc::reconstruct(phig/rAUf)
                );
                const volVectorField USigmaPressureResidual
                (
                    "USigmaPressureResidual",
                    rAU()*fvc::reconstruct
                    (
                        (phig - p_rghEqn.flux())/rAUf
                    )
                );

                // Hand the non-gradient part to the metrics writer when the
                // solver provides a store for it (the SL solver does; the plain
                // level-set solver does not, hence the registry guard). This is
                // the ONLY part of the capillary force that survives projection
                // and therefore the only part that can do work on the flow.
                if (mesh.foundObject<volVectorField>("USigmaResStore"))
                {
                    mesh.lookupObjectRef<volVectorField>("USigmaResStore")
                        == USigmaPressureResidual;
                }

                const scalarField& V = mesh.V();
                const vectorField& uPred =
                    USigmaPredictor.primitiveField();
                const vectorField& uRes =
                    USigmaPressureResidual.primitiveField();
                scalar pred2 = 0.0;
                scalar res2 = 0.0;
                scalar resInf = 0.0;
                forAll(uPred, c)
                {
                    pred2 += V[c]*magSqr(uPred[c]);
                    res2 += V[c]*magSqr(uRes[c]);
                    resInf = max(resInf, mag(uRes[c]));
                }
                reduce(pred2, sumOp<scalar>());
                reduce(res2, sumOp<scalar>());
                reduce(resInf, maxOp<scalar>());

                Info<< "rAU/rAUf capillary-pressure velocity residual: max="
                    << resInf << " m/s, relative-L2="
                    << Foam::sqrt(res2/max(pred2, VSMALL)) << endl;

                // ---- WP8.0: the residual in FLUX space, before reconstruct --
                // The pressure solve can only absorb the part of G_sigma that
                // lies in the range of snGrad; what it leaves face-by-face,
                // R_f = phig - p_rghEqn.flux(), is what fvc::reconstruct turns
                // into spurious velocity. Measuring R_f on the ACTIVE faces
                // (|snGrad(alpha)| > 0, the CSF support) and comparing it with
                // phig itself separates three candidate mechanisms:
                //   (i)   R_f/phig ~ O(1)      -> genuinely non-absorbable
                //                                 (non-gradient) capillary flux
                //   (ii)  R_f ~ linear-solver residual scale -> the projection
                //                                 is merely unconverged
                //   (iii) R_f small but max|U_res| large -> fvc::reconstruct is
                //                                 the amplifier
                // Reported per final corrector to capillaryFluxResidual.csv.
                {
                    const surfaceScalarField Rf("Rf", phig - p_rghEqn.flux());
                    const surfaceScalarField snA(fvc::snGrad(alpha1));

                    const scalarField& RfIn = Rf.primitiveField();
                    const scalarField& phigIn = phig.primitiveField();
                    const scalarField& snAIn = snA.primitiveField();
                    const scalar snMax = max(gMax(mag(snAIn)), VSMALL);

                    scalar gMaxAct = 0, rMaxAct = 0, g2 = 0, r2 = 0;
                    label nAct = 0;
                    forAll(RfIn, f)
                    {
                        if (mag(snAIn[f]) <= 1e-8*snMax) { continue; }
                        ++nAct;
                        gMaxAct = max(gMaxAct, mag(phigIn[f]));
                        rMaxAct = max(rMaxAct, mag(RfIn[f]));
                        g2 += sqr(phigIn[f]);
                        r2 += sqr(RfIn[f]);
                    }
                    reduce(nAct, sumOp<label>());
                    reduce(gMaxAct, maxOp<scalar>());
                    reduce(rMaxAct, maxOp<scalar>());
                    reduce(g2, sumOp<scalar>());
                    reduce(r2, sumOp<scalar>());

                    const scalar fracMax = rMaxAct/max(gMaxAct, VSMALL);
                    const scalar fracL2 =
                        Foam::sqrt(r2/max(g2, VSMALL));

                    Info<< "capillary flux residual (active faces): "
                        << "max|Rf|/max|phig| = " << fracMax
                        << ", L2 ratio = " << fracL2
                        << ", p iters = " << pIterations
                        << " (init " << pInitialResidual
                        << " -> final " << pFinalResidual << ")" << endl;

                    static autoPtr<OFstream> capResidFile;
                    if (Pstream::master())
                    {
                        if (!capResidFile.valid())
                        {
                            capResidFile.reset
                            (
                                new OFstream("capillaryFluxResidual.csv")
                            );
                            capResidFile().precision(10);
                            capResidFile()
                                << "TIME,nActiveFaces,phigMaxActive,"
                                   "RfMaxActive,RfFracMax,RfFracL2,"
                                   "maxUresidual,relL2Uresidual,"
                                   "pInitialResidual,pFinalResidual,pIterations"
                                << nl;
                        }
                        capResidFile()
                            << runTime.value() << ',' << nAct << ','
                            << gMaxAct << ',' << rMaxAct << ','
                            << fracMax << ',' << fracL2 << ','
                            << resInf << ','
                            << Foam::sqrt(res2/max(pred2, VSMALL)) << ','
                            << pInitialResidual << ',' << pFinalResidual << ','
                            << pIterations << nl;
                    }
                }
            }
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
