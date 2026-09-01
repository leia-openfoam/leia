{
    // DISCRETE MASS RESIDUAL FOR A MODEL THAT SOLVES NO DENSITY EQUATION.
    //
    // rhoLENT computes this inside rhoLENTEqn.H because it is the equation it
    // solves. interpolatedDensity and geometricFaceDensity solve nothing: they set
    // rho algebraically from alpha and build rhoPhi geometrically, so
    //
    //     R = ddt(rho) + div(rhoPhi)
    //
    // is whatever the geometry and the flux happen to leave it -- and it had never
    // been measured. R is the quantity the mass-momentum consistency argument turns
    // on: for a droplet translating at uniform U0 the momentum transient plus
    // convection collapses, for constant U, to U0 * R, a source the pressure cannot
    // absorb because it is not a gradient, and one that vanishes IDENTICALLY at
    // U0 = 0. That is the stationary/translating asymmetry stated in one line.
    //
    // Pure diagnostic: no field is modified, and it runs only under
    // massFlux { massResidualDiagnostic true; } (default false).
    //
    // ddt(rho) uses the ddt(rho) scheme, which is a DIFFERENT key from momentum's
    // ddt(rho,U). Where the two differ -- Euler here against backward there -- the
    // momentum equation does not inherit U0*R but U0*R plus the difference of the
    // two ddt operators applied to rho. The banner prints both so the pairing is
    // never in doubt.
    const volScalarField rhoDdtDiag
    (
        IOobject
        (
            "rhoDdtDiag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::ddt(rho)
    );
    const volScalarField rhoDivDiag
    (
        IOobject
        (
            "rhoDivDiag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::div(rhoPhi)
    );
    const volScalarField massResidual
    (
        IOobject
        (
            "massResidual",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rhoDdtDiag + rhoDivDiag
    );

    // COLLECTIVE reductions: every rank reaches them; only master writes.
    const scalar totV = Foam::max(gSum(mesh.V()), VSMALL);
    rhoMassResidualL1 =
        fvc::domainIntegrate(mag(massResidual)).value()/totV;
    rhoMassResidualLinf = gMax(mag(massResidual.primitiveField()));
    const scalar balanceScale =
        fvc::domainIntegrate(mag(rhoDdtDiag) + mag(rhoDivDiag)).value()/totV;
    rhoMassResidualRelL1 =
        rhoMassResidualL1/Foam::max(balanceScale, VSMALL);

    Info<< "mass residual (" << massFluxType
        << "): ddt(rho)+div(rhoPhi) L1/Linf/relative-L1 = "
        << rhoMassResidualL1 << ' ' << rhoMassResidualLinf << ' '
        << rhoMassResidualRelL1 << nl << endl;
}
