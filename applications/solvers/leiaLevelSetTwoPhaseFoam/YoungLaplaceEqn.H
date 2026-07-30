// Initialise pressure from the same integrated capillary face flux used by
// PIMPLE. This is the sole production surface-tension delivery path.
surfaceScalarField fSigmaInitialFlux
(
    "fSigmaInitialFlux",
    fSigma->faceSurfaceTensionForceFlux()
);
// snGrad-based fluxes are owner-oriented even though the generic model field is
// not necessarily tagged as such.
fSigmaInitialFlux.setOriented();

fvScalarMatrix YoungLaplaceEqn
(
    fvm::laplacian(p_rgh) == fvc::div(fSigmaInitialFlux)
);
YoungLaplaceEqn.solve();
