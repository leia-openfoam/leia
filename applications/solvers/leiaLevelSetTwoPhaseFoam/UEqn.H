    MRF.correctBoundaryVelocity(U);

    // Evaluate the integrated scalar capillary force flux once and reuse this
    // exact field in both the momentum predictor and pressure corrector.
    // G_sigma already includes |Sf|; no downstream code may reinterpret a
    // model result as a pointwise face force density.
    tmp<surfaceScalarField> tGSigma =
        fSigma->faceSurfaceTensionForceFlux();
    const surfaceScalarField& GSigma = tGSigma();

    // ONE-FIELD DNS VISCOUS TERM, written out.
    //
    // This replaces turbulence->divDevRhoReff(rho, U). We only ever run DNS -- laminar
    // transport, Newtonian phases -- so the wrapper resolves to linearViscousStress with
    // nuEff = nu, and writing it explicitly exposes the momentum equation actually being
    // solved instead of hiding it behind a turbulence interface. What the wrapper does is
    //
    //     divDevRhoReff(rho,U) = - fvc::div(mu dev2(T(grad U))) - fvm::laplacian(mu, U)
    //
    // i.e. minus the divergence of mu (grad U + (grad U)^T - (2/3)(div U) I). Both parts
    // are kept: the implicit Laplacian AND the explicit TRANSPOSE part, which for a
    // variable-viscosity two-phase flow does not vanish -- it is proportional to the
    // viscosity gradient and therefore lives exactly at the interface.
    //
    // The face viscosity muf is the coefficient of the implicit Laplacian, so it is what
    // the discretisation actually sees. It is a surfaceScalarField constructed in
    // createFields.H and, when massFlux.viscosityFaceModel selects it, filled from the
    // SAME geometric alpha_f that builds rho_f for the mass flux. That remains
    // runtime-selectable rather than unconditional because it demonstrably CHANGES the
    // answer (MEASURED 2026-09-03: arithmetic mu_f from alpha_f is 2.08x worse in L1 than
    // the interpolated default), so it is a model choice, not a refactor.
    //
    // With the default viscosityFaceModel = interpolated, muf = fvc::interpolate(mu) and
    // this assembly is bit-identical to the wrapper it replaces, because
    // fvm::laplacian(volScalarField, U) interpolates its coefficient to the faces anyway.
    const surfaceTensorField gradUTf(dev2(T(fvc::interpolate(fvc::grad(U)))));

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      - fvm::laplacian(muf, U)
      - fvc::div(muf*(mesh.Sf() & gradUTf))
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
