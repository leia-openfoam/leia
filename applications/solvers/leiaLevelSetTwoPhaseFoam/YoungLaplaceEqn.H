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

// REFERENCE THE SINGULAR SYSTEM. p_rgh is zeroGradient on every wall and
// needReference() is true, so this Laplacian has a constant null space. It used
// to be solved unreferenced: it converged only because a closed-domain
// face-flux divergence happens to sum to zero, i.e. the RHS was orthogonal to
// the null space by luck rather than by construction, and GAMG on a singular
// system is a known stall risk. pEqn.H references its own solve; so does this
// one now, with the same idiom (pin the reference cell to its current value, a
// pure shift).
YoungLaplaceEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

// SOLVE IT TIGHT. A bare solve() resolves to the p_rgh solver entry, whose
// relTol is 0.01 in every droplet case. Starting from p_rgh = 0 that stops after
// a 100x residual drop, leaving the initial pressure ~1% off -- so the first
// momentum predictor saw an unbalanced impulse of order 0.01*sigma*kappa, in the
// one step where the balance is supposed to be exact. p_rgh.select(true) picks
// the p_rghFinal entry (relTol 0), the same one pEqn uses on its final inner
// iteration.
//
// STILL DIFFERENT, and deliberately so: the operator here is
// fvm::laplacian(p_rgh) with unit coefficient, while every subsequent step
// projects with fvm::laplacian(rAUf, p_rgh) against a source that also carries
// rAUf. For a spatially constant sigma*kappa the two coincide and the balance is
// exact either way, which is why the constant-curvature control cannot see this.
// For varying curvature they differ, and rAUf = 1/A varies by roughly the
// density ratio across the interface. Making them agree needs 1/A, which does
// not exist before the first UEqn is assembled; a uniform stand-in changes
// nothing, since laplacian(c*p) = div(c*G) is the same equation as
// laplacian(p) = div(G) for constant c. Left as is, and recorded, because pEqn
// re-projects on the first step.
YoungLaplaceEqn.solve(mesh.solver(p_rgh.select(true)));
