{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Level-set transport in the ADVECTIVE (Hamilton-Jacobi) form:
    // (v.grad)psi = div(phi psi) - psi (div phi) -- the conservative form
    // alone equals it only for a discretely solenoidal flux, which the
    // intermediate PIMPLE fluxes are not. RHS: the SDPLS source
    // S = f_nl(psi^n) psi (vanishes on the zero level set; noSource = off).
    //
    // DEFECT CORRECTION: this equation is re-assembled and re-solved on
    // EVERY PIMPLE outer corrector (this file is included inside
    // pimple.loop() without a firstIter() guard) -- nOuterCorrectors is
    // the deferred-correction control for second-order (linearUpwind)
    // div(phi,psi) schemes, and psi is re-advected with each corrector's
    // updated flux.
    const volScalarField divPhi("divPhi", fvc::div(phi));

    fvScalarMatrix psiEqn
    (
        fvm::ddt(psi)
      + fvm::div(phi, psi)
      - fvm::Sp(divPhi, psi)
     ==
        source->fvmsdplsSource(psi, U)
    );

    psiEqn.solve();

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // Interface updates ONCE per time step (final outer corrector): the
    // band refresh, the (criterion-gated) redistancing and the volume
    // fraction. Re-running these every corrector would repeatedly displace
    // the interface within one step (redistancing) and recompute the
    // tet-fill indicator against a psi that is still iterating.
    if (pimple.finalIter())
    {
        // Fresh sign-change band BEFORE redistancing (plane-anchored
        // models and the signedDistanceBounds criterion read it).
        narrowBand->calc();

        if (redist->redistance(psi))
        {
            // Cut cells may realign sub-h: refresh the band.
            narrowBand->calc();
        }

        phaseInd->calcPhaseIndicator(alpha1, psi);

        alpha2 = 1.0 - alpha1;
        mixture.correct();

        scalar domainFraction =
            alpha1.weightedAverage(mesh.Vsc()).value();

        Info<< "Phase-1 volume fraction = "
            << domainFraction
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }
}
