/*---------------------------------------------------------------------------*\
    Per-time-step phase-volume correction, constructed for the EULERIAN
    two-phase solver only.

Why this is its own header and not part of createFields.H
---------------------------------------------------------
createFields.H is SHARED: leiaSemiLagrangianLevelSetTwoPhaseFoam includes the
same file ("All momentum/pressure/mesh machinery (createFields.H, UEqn.H,
pEqn.H, ...)" -- its own header comment). Constructing a volumeCorrection there
therefore requires the TYPE to be declared at file scope in BOTH solver .C
files, and the semi-Lagrangian solver belongs to a parallel session: adding an
include to it is a merge conflict waiting to happen for no benefit, since that
solver runs slAlphaEqn.H and would never call correct().

So the construction lives here, included from leiaLevelSetTwoPhaseFoam.C
immediately AFTER createFields.H -- at which point `mesh` and `alpha1` both
exist. Net footprint on the semi-Lagrangian line: zero files, zero bytes.

Ordering, which is load-bearing
-------------------------------
Constructed AFTER createFields.H has run phaseInd->calcPhaseIndicator(alpha1,
psi), so the target volume is captured from a REAL alpha1 rather than from the
zero the field was created with. The target must be the SAME functional the
error metric reports -- gSum(alpha*V), advectionErrors.H -- or the correction
would faithfully conserve a quantity nobody measures.

Default is noVolumeCorrection (the base class itself, exactly as sdplsSource's
base IS noSource): an absent levelSet/volumeCorrection sub-dictionary leaves psi
untouched for the entire run, so every existing case is unchanged.
\*---------------------------------------------------------------------------*/

autoPtr<volumeCorrection> volCorr = volumeCorrection::New(mesh);
volCorr->setTargetVolume(alpha1);
