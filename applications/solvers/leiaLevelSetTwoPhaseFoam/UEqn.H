    MRF.correctBoundaryVelocity(U);

    // Evaluate the integrated scalar capillary force flux once and reuse this
    // exact field in both the momentum predictor and pressure corrector.
    // G_sigma already includes |Sf|; no downstream code may reinterpret a
    // model result as a pointwise face force density.
    tmp<surfaceScalarField> tGSigma =
        fSigma->faceSurfaceTensionForceFlux();
    const surfaceScalarField& GSigma = tGSigma();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                GSigma
              - ghf*fvc::snGrad(rho)*mesh.magSf()
              - fvc::snGrad(p_rgh)*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
