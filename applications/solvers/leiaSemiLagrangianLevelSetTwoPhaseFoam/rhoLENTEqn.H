{
    // Auxiliary mass conservation equation from rhoLENT (Liu et al., JCP
    // 493:112426, 2023, Eq. 40 / Algorithm 1):
    //
    //   ddt(rho) + div(rhof^{n+1} phi^o) = 0.
    //
    // The geometrical face density rhof^{n+1} is fixed by the newly advected
    // interface. The volumetric flux phi^o is the same outer-iteration flux
    // that forms rhoPhi in div(rhoPhi,U). Reassembling this diagonal equation
    // on every PIMPLE outer iteration therefore makes the mass flux in the
    // auxiliary density and momentum equations identical. Do not bound rho:
    // clipping would destroy the discrete free-stream identity.
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho) + fvc::div(rhoPhi)
    );
    rhoEqn.solve();

    // OPTIONAL DENSITY BOUNDING (massFlux { boundRho true; }; default false, so
    // every existing study is byte-unchanged). rho1 and rho2 are CONSTANTS, so any
    // rho outside [rho2, rho1] is unphysical by construction -- and a negative one
    // flips the sign of 1/rho in the pressure Laplacian, costing definiteness and
    // ending the run in a GAMG blow-up rather than a result.
    //
    // The clip is applied HERE, before the residual below and before momentum and
    // pressure read rho, because those are what the unbounded value destroys.
    //
    // IT IS NOT FREE. Clipping injects a source S_clip, so the discrete identity
    // becomes ddt(rho) + div(rhoPhi) = S_clip instead of 0, and the momentum
    // equation inherits U*S_clip -- exactly the free-stream identity rhoLENT exists
    // to enforce. So the clip is MEASURED, never silent: rhoClipL1 (volume-mean
    // |clip|) and rhoClipFraction (fraction of cells clipped) reach the metrics CSV,
    // and the mass residual printed below is computed AFTER the clip, so it reports
    // the equation the momentum solve actually inherits.
    if (boundRhoLENT)
    {
        scalarField& rhoIn = rho.primitiveFieldRef();
        const scalarField& Vc = mesh.V();
        const scalar rMin = rho2.value();
        const scalar rMax = rho1.value();
        scalar clipVol = 0;
        label nClipped = 0;
        forAll(rhoIn, celli)
        {
            const scalar r = rhoIn[celli];
            const scalar rc = Foam::min(Foam::max(r, rMin), rMax);
            if (rc != r)
            {
                clipVol += mag(rc - r)*Vc[celli];
                ++nClipped;
            }
            rhoIn[celli] = rc;
        }
        // COLLECTIVE: every rank must reach these reductions (a Pstream::master()
        // guard around a reduction deadlocked four cluster arms for 76 minutes).
        reduce(clipVol, sumOp<scalar>());
        reduce(nClipped, sumOp<label>());
        const label nTot = returnReduce(mesh.nCells(), sumOp<label>());
        rhoClipL1 = clipVol/Foam::max(gSum(mesh.V()), VSMALL);
        rhoClipFraction = scalar(nClipped)/Foam::max(scalar(nTot), VSMALL);
    }

    rho.correctBoundaryConditions();

    const volScalarField rhoDdtDiagnostic
    (
        IOobject
        (
            "rhoDdtDiagnostic",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::ddt(rho)
    );
    const volScalarField rhoFluxDivDiagnostic
    (
        IOobject
        (
            "rhoFluxDivDiagnostic",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::div(rhoPhi)
    );
    const volScalarField rhoMassResidual
    (
        IOobject
        (
            "rhoMassResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rhoDdtDiagnostic + rhoFluxDivDiagnostic
    );
    rhoMassResidualL1 =
        fvc::domainIntegrate(mag(rhoMassResidual)).value()
       /Foam::max(gSum(mesh.V()), VSMALL);
    rhoMassResidualLinf =
        gMax(mag(rhoMassResidual.primitiveField()));
    const scalar rhoBalanceScale =
        fvc::domainIntegrate
        (
            mag(rhoDdtDiagnostic) + mag(rhoFluxDivDiagnostic)
        ).value()/Foam::max(gSum(mesh.V()), VSMALL);
    rhoMassResidualRelL1 =
        rhoMassResidualL1/Foam::max(rhoBalanceScale, VSMALL);

    Info<< "rhoLENT auxiliary density: min/max(rho) = "
        << gMin(rho.primitiveField()) << ' '
        << gMax(rho.primitiveField())
        << ", mass residual L1/Linf/relative-L1 = "
        << rhoMassResidualL1 << ' ' << rhoMassResidualLinf << ' '
        << rhoMassResidualRelL1 << nl << endl;
}
