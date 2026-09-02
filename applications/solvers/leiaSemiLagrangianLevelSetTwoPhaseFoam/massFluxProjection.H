{
    // MASS-FLUX PROJECTION (author's proposal, 2026-09-02).
    //
    // A geometric mass flux does not satisfy the discrete mass balance: rho is set
    // algebraically from alpha and rhoPhi is built from the reconstructed planes, so
    //
    //     R_c = ddt(rho)_c + div(rhoPhi)_c
    //
    // is whatever the geometry leaves it -- measured at 4.79e-02 RELATIVE for
    // geometricFaceDensity (config/rhoBoundGate2D). For a droplet translating at
    // uniform U0 the momentum equation inherits U0*R, which is not a gradient, so the
    // pressure cannot absorb it, and which vanishes identically at U0 = 0.
    //
    // rhoLENT removes R by SOLVING for an auxiliary rho -- giving up rho = rho_geom,
    // and that auxiliary rho is unbounded (it needs massFlux.boundRho, whose clip
    // breaks the very identity being enforced). This construction removes R the other
    // way: leave rho purely GEOMETRIC, hence bounded in [rho2, rho1] by construction,
    // and correct the FLUX instead.
    //
    //     want   rhoPhi* = rhoPhi - rhoPhi'   with   div(rhoPhi') = R
    //
    // One equation per cell for a face field is underdetermined, so take the curl-free
    // member -- the minimal-L2 correction, i.e. the smallest change to the flux that
    // fixes the balance:
    //
    //     rhoPhi' = snGrad(lambda) * magSf,     laplacian(lambda) = R
    //
    // lambda is a scalar POTENTIAL for the flux correction; it is to the mass flux
    // exactly what pressure is to the volumetric flux in the projection step, and it
    // has no physical identity of its own. Its dimensions follow from the equation:
    // [lambda] = kg/(m s).
    //
    // COMPATIBILITY. All-Neumann (the correction vanishes where the flux is already
    // right: walls with phi = 0, an inlet with prescribed U and known alpha), so
    // lambda is fixed only up to a constant -- hence the reference cell -- and the
    // problem is solvable only if the volume integral of R vanishes. IT DOES NOT HERE,
    // and the reason is diagnostic rather than incidental: the integral is the rate of
    // change of total mass, and the droplet INFLATES by +5 to +8 percent over the
    // horizon (config/translatingClearOutlet2D), a defect of the semi-Lagrangian
    // transport that no flux redistribution can repair. Subtracting the mean is the
    // standard remedy, and what is subtracted IS that global drift. So this projection
    // fixes the DISTRIBUTION of the imbalance and leaves its GLOBAL part alone -- two
    // different defects, cleanly separated, and both are reported below.
    volScalarField massFluxResidual
    (
        IOobject
        (
            "massFluxResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::ddt(rho) + fvc::div(rhoPhi)
    );

    const scalar totVproj = Foam::max(gSum(mesh.V()), VSMALL);
    const scalar residualBefore =
        fvc::domainIntegrate(mag(massFluxResidual)).value()/totVproj;
    // The global part a flux redistribution provably cannot remove: report it, do not
    // hide it. This is the interface transport's own mass drift.
    const scalar residualGlobal =
        fvc::domainIntegrate(massFluxResidual).value()/totVproj;

    massFluxResidual -= dimensionedScalar
    (
        "residualGlobal", massFluxResidual.dimensions(), residualGlobal
    );

    fvScalarMatrix lambdaEqn
    (
        fvm::laplacian(massFluxPotential) == massFluxResidual
    );
    lambdaEqn.setReference(0, 0);
    lambdaEqn.solve();

    const surfaceScalarField rhoPhiCorr
    (
        "rhoPhiCorr",
        fvc::snGrad(massFluxPotential)*mesh.magSf()
    );
    rhoPhi -= rhoPhiCorr;

    // Report what the projection cost and what it achieved, measured on the CORRECTED
    // flux -- never asserted from the fact that a solve happened.
    const volScalarField residualAfterField
    (
        IOobject
        (
            "massFluxResidualAfter",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::ddt(rho) + fvc::div(rhoPhi)
    );
    massProjResidualBefore = residualBefore;
    massProjResidualAfter =
        fvc::domainIntegrate(mag(residualAfterField)).value()/totVproj;
    massProjGlobalDrift = mag(residualGlobal);
    massProjCorrFraction =
        gSum(mag(rhoPhiCorr.primitiveField()))
       /Foam::max(gSum(mag(rhoPhi.primitiveField())), VSMALL);

    Info<< "mass-flux projection: residual L1 " << massProjResidualBefore
        << " -> " << massProjResidualAfter
        << ", global drift (irreducible) " << massProjGlobalDrift
        << ", |rhoPhi'|/|rhoPhi| = " << massProjCorrFraction << nl << endl;
}
