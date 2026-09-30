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
    // The flux that transports psi: the velocity extension's (default `none`:
    // a copy of phi, byte-identical to the transport before 2026-09-26).
    velExt->correct();
    const surfaceScalarField& phiPsi = velExt->phi();
    const volScalarField divPhi("divPhi", fvc::div(phiPsi));
    source->setTransportFlux(phiPsi);

    fvScalarMatrix psiEqn
    (
        fvm::ddt(psi)
      + fvm::div(phiPsi, psi)
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

    // Interface updates. The band refresh, the (criterion-gated) redistancing
    // and the volume correction run ONCE per time step (final outer corrector):
    // re-running redistancing every corrector would repeatedly displace the
    // interface within one step.
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

        // Uniform-shift volume correction, AFTER redistancing and BEFORE the
        // phase indicator. After, because redistancing rebuilds psi from the
        // zero level set and is NOT volume preserving -- the shift must be
        // the last operation of the step that touches psi, or a redistance
        // event discards the eps just chosen. It cannot undo what
        // redistancing established: eps is spatially constant, so
        // grad(psi + eps) = grad(psi) and both redistancing triggers
        // (signedDistanceBounds, gradPsiThreshold) are functions of
        // |grad psi| only. Before the indicator, because alpha1 -- and
        // through it alpha2, the mixture, rho/mu in UEqn/pEqn, and the CSV
        // row -- must be built from the corrected field.
        //
        // Default noVolumeCorrection returns false without touching psi.
        if (volCorr->correct(psi, phaseInd()))
        {
            // The zero crossing moved sub-h: the sign-change band can gain
            // or lose a cell, exactly as after a redistance event above.
            narrowBand->calc();
        }
    }
    else if (psiOuterCorrectors)
    {
        // The band of the current psi iterate, for the phase indicator and
        // the face density below.
        narrowBand->calc();
    }

    // The phase indicator, the mixture and the face density follow the
    // interface: on the final outer corrector (psiOuterCorrectors no), or on
    // every outer corrector from the current psi iterate (yes).
    if (pimple.finalIter() || psiOuterCorrectors)
    {
        phaseInd->calcPhaseIndicator(alpha1, psi);

        alpha2 = 1.0 - alpha1;
        mixture.correct();

        // Face area fraction, face density rho_f and face viscosity from the
        // NEW interface; the geometric cell density for the mass-flux modes
        // that solve no density equation (shared with the SL solver).
        #include "updateFaceDensity.H"
    }

    if (pimple.finalIter())
    {
        scalar domainFraction =
            alpha1.weightedAverage(mesh.Vsc()).value();

        Info<< "Phase-1 volume fraction = "
            << domainFraction
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;
    }

    // The mass flux rhoPhi and, for rhoLENT, the auxiliary density, on EVERY
    // outer corrector: pEqn changes phi, so rhoPhi must be reassembled before
    // the next UEqn (shared with the SL solver).
    #include "updateMassFlux.H"

    // The density bound breaks the identity rhoLENT enforces, so it is
    // reported, never silent. This solver writes no droplet CSV yet, so the
    // two numbers go to the log.
    if (useRhoLENT && boundRhoLENT && !freezeRhoPhi)
    {
        Info<< "rhoLENT density bound: rhoClipL1 = " << rhoClipL1
            << ", rhoClipFraction = " << rhoClipFraction << nl << endl;
    }
}
