// Initialize fields for error analysis.
// phi0 (a full surfaceScalarField) is read ONLY by oscillateVelocity() on the
// reversed/oscillating path; the steady path never touches it, so only allocate
// it when oscillating (saves ~one surfaceScalarField on steady cases).
autoPtr<surfaceScalarField> phi0Ptr;
if (velocityModel->isOscillating())
{
    phi0Ptr.reset(new surfaceScalarField("phi0", phi));
}
volVectorField U0("U0", U);

// Previous-step velocity u^n reused by the semi-Lagrangian trajectory
// (updated each step; = U0 scaled by the oscillation factor at t^n).
volVectorField Uold("Uold", U);

// Use psiEnd if present as field, otherwise copy psi
volScalarField psi0
(
    IOobject
    (
        "psiEnd",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    psi
);

// Use alphaEnd if present as field, otherwise copy alpha
volScalarField alpha0
(
    IOobject
    (
        "alphaEnd",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha
);

// Open error file for output
OFstream errorFile("leiaSemiLagrangeLevelSetFoam.csv");

// CSV Header (same columns as leiaLevelSetFoam.csv for direct comparison)
errorFile << "TIME,"
    << "DELTA_X,"
    << "L_INF_E_PSI,"
    << "E_VOL_ALPHA,"
    << "E_VOL_ALPHA_REL,"
    << "E_GEOM_ALPHA,"
    << "E_GEOM_ALPHA_REL,"
    << "E_BOUND_ALPHA,"
    << "ELAPSED_CPU_TIME,"
    << "MAX_CFL,"
    << "MAX_ABS_SUM_VOL_FLUXES\n";

#include "CourantNo.H"

// Report initial errors
reportErrors(errorFile, psi, psi0, alpha, alpha0, phi, CoNum);
