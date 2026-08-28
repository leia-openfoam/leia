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

// Trace-velocity source (levelSet.semiLagrangian.traceVelocity):
//   cellCentred   -- (default) trace with the prescribed CELL velocity.
//   projectedFlux -- trace with fvc::reconstruct(phi), the cell field of the
//                    DISCRETELY DIVERGENCE-FREE face flux. This is what
//                    Basilisk-style advection transports with. The prescribed
//                    cell velocity is only analytically solenoidal; its
//                    discrete divergence is not zero, and a semi-Lagrangian
//                    map preserves the enclosed volume exactly only when the
//                    tracing field is discretely divergence-free. The reversed
//                    vortex, whose exact answer is the initial condition, is
//                    the sharpest test of the difference.
const word slTraceVelocity
(
    static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("semiLagrangian")
        .getOrDefault<word>("traceVelocity", "cellCentred")
);
if (slTraceVelocity != "cellCentred" && slTraceVelocity != "projectedFlux")
{
    FatalErrorInFunction
        << "levelSet.semiLagrangian.traceVelocity = '" << slTraceVelocity
        << "'. Valid: cellCentred, projectedFlux." << exit(FatalError);
}
const bool traceFromFlux = (slTraceVelocity == "projectedFlux");
autoPtr<volVectorField> UtracePtr;
autoPtr<volVectorField> UtraceOldPtr;
if (traceFromFlux)
{
    Info<< "SL trace velocity: reconstruct(phi) -- divergence-free face flux"
        << nl << endl;
    UtracePtr.reset(new volVectorField("Utrace", fvc::reconstruct(phi)));
    UtraceOldPtr.reset(new volVectorField("UtraceOld", UtracePtr()));
}

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
errorFile << Foam::advectionErrorsCsvHeader();

#include "CourantNo.H"

// Report initial errors
reportErrors(errorFile, psi, psi0, alpha, alpha0, phi, CoNum);
