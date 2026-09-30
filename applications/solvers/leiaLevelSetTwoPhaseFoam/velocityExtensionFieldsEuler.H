// The velocity-extension composition root of the Eulerian two-phase solver
// (2026-09-26, docs/plan-halo-limited-gradient-control.md step C3). Included only
// from leiaLevelSetTwoPhaseFoam.C: the shared createFields.H cannot hold it,
// because the semi-Lagrangian solver builds its own velocityExtension and a mesh
// holds exactly one (it registers Uext). alphaEqn.H calls velExt->correct() on
// every outer corrector and transports psi with velExt->phi() in material form.
// The default model `none` makes velExt->phi() a copy of phi named "phi", so the
// div(phi,psi) scheme applies and every earlier run is byte-identical.
autoPtr<velocityExtension> velExt(velocityExtension::New(mesh));
velExt->correct();
