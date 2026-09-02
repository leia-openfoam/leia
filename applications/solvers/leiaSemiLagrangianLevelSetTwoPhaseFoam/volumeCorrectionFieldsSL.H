/*---------------------------------------------------------------------------*\
    PER-TIME-STEP PHASE-VOLUME CORRECTION for the SEMI-LAGRANGIAN two-phase
    solver.
\*---------------------------------------------------------------------------*/
//
// WHY THIS EXISTS, and why the previous decision not to have it was wrong.
//
// applications/solvers/leiaLevelSetTwoPhaseFoam/volumeCorrectionFields.H states
// that the correction is constructed "for the EULERIAN two-phase solver only",
// because the semi-Lagrangian solver "runs slAlphaEqn.H and would never call
// correct()". That reasoning is retracted: it is exactly backwards.
//
// Liu, Tolle, Bothe & Maric, JCP 493 (2023) 112426 -- the rhoLENT paper this
// solver's mass-flux model is taken from -- says in section 3.4:
//
//   "Level Set / Front Tracking methods are Lagrangian / Eulerian methods that
//    kinematically evolve the fluid interface WITHOUT UTILIZING FLUXES through
//    control-volume boundaries and are therefore INHERENTLY NOT MASS(VOLUME)-
//    CONSERVATIVE ... we ensure volume conservation by extending/contracting the
//    Front with a modified iso-value."
//
// and closes it with eqs. (46)-(49): a uniform normal shift of the iso-surface,
// h = (V_target - V_ini)/int dS. rhoLENT does not TOLERATE the interface volume
// error, it REPAIRS it every step -- and the semi-Lagrangian solver is precisely
// the non-flux-based, non-conservative transport the paper is describing.
// pointValueScheme.H says "Non-conservative" in its own header.
//
// MEASURED CONSEQUENCE of leaving it out (2026-09-02, translating droplet):
// the droplet INFLATES by +5 to +8 percent over the horizon, and at matched
// times the gain does NOT converge under refinement. Three separate attempts to
// remove the resulting mass-flux inconsistency through the FACE-DENSITY
// construction all failed -- two time levels (2%), an SL-advected donor plane
// (1%), and a curl-free flux projection, which removes the residual but stops
// the droplet. They failed because the error is not in the face density; it is
// the transport's own volume drift, which is what this corrects.
//
// newtonShiftVolumeCorrection is the same construction as the paper's: a single
// scalar eps with psi <- psi + eps chosen so that sum_c alpha_c(psi+eps)|Omega_c|
// returns to the captured target. Distance-preserving, so grad(psi) is untouched.
//
// ORDERING is load-bearing, exactly as in the Eulerian twin: constructed AFTER
// createFields.H has run calcPhaseIndicator, so the target volume is captured
// from a REAL alpha1 rather than from the zero the field was created with.
//
// Default noVolumeCorrection is the base class and a true no-op, so every study
// that does not select it is byte-unchanged.
autoPtr<volumeCorrection> volCorr = volumeCorrection::New(mesh);
volCorr->setTargetVolume(alpha1);
