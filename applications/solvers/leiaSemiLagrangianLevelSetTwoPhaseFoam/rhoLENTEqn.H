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
