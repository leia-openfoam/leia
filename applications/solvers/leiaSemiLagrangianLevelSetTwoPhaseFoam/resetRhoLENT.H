if (useRhoLENT && !freezeRhoPhi)
{
    // rho^{n+1} from the auxiliary equation has now served its purpose in the
    // converged pressure-velocity coupling. Reset the current-time density to
    // the Detrixhe-Aslam cell indicator, as rhoLENT Algorithm 1 step 13 does.
    // On the next ++runTime this geometrical field becomes rho.oldTime(), while
    // the next auxiliary solve again produces the momentum-consistent rho.
    const volScalarField rhoGeometric
    (
        IOobject
        (
            "rhoGeometric",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        alpha1*rho1 + alpha2*rho2
    );

    rhoLENTResetL1 =
        fvc::domainIntegrate(mag(rho - rhoGeometric)).value()
       /Foam::max(gSum(mesh.V()), VSMALL);

    rho == rhoGeometric;
    rho.correctBoundaryConditions();

    // p is a derived diagnostic field; p_rgh remains the pressure unknown.
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

    Info<< "rhoLENT reset from Detrixhe-Aslam alpha: mean |rhoAux-rhoGeom| = "
        << rhoLENTResetL1 << nl << endl;
}
