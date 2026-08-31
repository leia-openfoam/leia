/*---------------------------------------------------------------------------*\
    Transport-side fields shared by the semi-Lagrangian velocity-extension and
    rhoLENT mass/momentum-consistency gates. Included AFTER createFields.H:
    U, phi, psi, alpha1, rho, rhoPhi and the phase densities must exist.
\*---------------------------------------------------------------------------*/

// Normal-constant extension used ONLY to transport psi. The physical U and phi
// used by momentum remain untouched. "none" is the identity model.
autoPtr<velocityExtension> slVelExt = velocityExtension::New(mesh);
slVelExt->correct();

// The SL trajectory needs two velocity levels. Uext is recomputed from the
// current physical U and current psi before advection; UextOld stores the
// previous pre-momentum extension. Persist it so a restarted second-order SL
// trajectory uses the same two velocity levels as an uninterrupted run.
volVectorField UextOld
(
    IOobject
    (
        "UextOld",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    slVelExt->Uext()
);

// TRACE VELOCITY SOURCE (levelSet.semiLagrangian.traceVelocity):
//   cellCentred   -- (default) the velocity-extension field Uext, i.e. the
//                    cell-centred momentum velocity. Historical behaviour,
//                    byte-identical.
//   projectedFlux -- the cell field RECONSTRUCTED FROM THE PROJECTED FACE FLUX,
//                    fvc::reconstruct(phi). This is the Basilisk property: the
//                    interface is advected by the DISCRETELY DIVERGENCE-FREE
//                    velocity, not by a cell field that only approximates one.
//                    For a semi-Lagrangian scheme the enclosed volume is
//                    preserved exactly when the flow map is volume-preserving,
//                    which requires a solenoidal tracing velocity; phi is
//                    discretely solenoidal after the pressure correction and
//                    the cell-centred U is not.
//   reconstructedU -- THE CONFOUND CONTROL for projectedFlux, and nothing else:
//                    reconstruct(linearInterpolate(Uext) & Sf). It shares the
//                    EXTENSION (Uext) with cellCentred and the RECONSTRUCT
//                    operator with projectedFlux, but its face field is the
//                    plain interpolated flux, which is NOT the projected one
//                    and therefore NOT discretely divergence-free.
//                    Why it exists: cellCentred and projectedFlux differ in
//                    three ways at once -- solenoidality, extension, and the
//                    reconstruct smoothing -- so a projectedFlux win cannot be
//                    attributed to solenoidality without this arm. If
//                    reconstructedU tracks cellCentred, solenoidality is the
//                    mechanism; if it tracks projectedFlux, the mechanism is
//                    the reconstruct smoothing and the divergence-free reading
//                    is wrong.
//   reconstructedMomentum -- the SECOND control, which completes a ladder in
//                    which every step changes exactly ONE thing:
//                      cellCentred           Uext,          no reconstruct
//                      reconstructedU        Uext,          reconstruct
//                      reconstructedMomentum U (raw),       reconstruct
//                      projectedFlux         phi (projected), reconstruct
//                    so cellCentred -> reconstructedU isolates the RECONSTRUCT
//                    SMOOTHING, reconstructedU -> reconstructedMomentum isolates
//                    the off-interface EXTENSION, and reconstructedMomentum ->
//                    projectedFlux isolates SOLENOIDALITY (the only difference
//                    being whether the face field is the pressure-projected flux
//                    or the plain interpolated one).
//                    MEASURED 2026-08-31 on config/traceAmplifierDt2D: at the
//                    base step reconstructedU already removes 70% of the
//                    cellCentred amplifier with no solenoidality whatsoever
//                    (dr 85.85 -> 26.05), so the divergence-free property is a
//                    minority share and this arm is what splits the other 70%
//                    between extension and smoothing.
//                    NOTE none of these carries a tunable coefficient: they are
//                    fixed discretisation operators, not filters.
//                    MEASURED 2026-08-31: with levelSet.velocityExtension.type
//                    = none -- the default, and what every stationary-droplet
//                    study uses -- the extension is the IDENTITY, so Uext == U
//                    bitwise and reconstructedMomentum is BYTE-IDENTICAL to
//                    reconstructedU (verified on 267 steps, np 4). The extension
//                    rung of the ladder is therefore a null step in that
//                    configuration and the whole non-solenoidal share belongs to
//                    the reconstruct operator. This arm only separates from
//                    reconstructedU when an extension model is actually active.
const word slTraceVelocity
(
    static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("semiLagrangian")
        .getOrDefault<word>("traceVelocity", "cellCentred")
);
if
(
    slTraceVelocity != "cellCentred"
 && slTraceVelocity != "projectedFlux"
 && slTraceVelocity != "reconstructedU"
 && slTraceVelocity != "reconstructedMomentum"
)
{
    FatalErrorInFunction
        << "levelSet.semiLagrangian.traceVelocity = '" << slTraceVelocity
        << "'. Valid: cellCentred, projectedFlux, reconstructedU,"
        << " reconstructedMomentum." << exit(FatalError);
}
const bool traceReconstructU = (slTraceVelocity == "reconstructedU");
const bool traceReconstructMom = (slTraceVelocity == "reconstructedMomentum");
if
(
    traceReconstructMom
 && static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("velocityExtension")
        .getOrDefault<word>("type", "none") == "none"
)
{
    // Not fatal -- it is a legitimate selection -- but running it as a study
    // ARM alongside reconstructedU wastes a job on a bit-identical duplicate.
    WarningInFunction
        << "traceVelocity reconstructedMomentum with velocityExtension type"
        << " none: the extension is the identity, so Uext == U and this is"
        << " BIT-IDENTICAL to reconstructedU. It separates from it only when"
        << " an extension model is active." << nl << endl;
}
// Both non-default sources trace a RECONSTRUCTED cell field and so share the
// Utrace/UtraceOld two-level machinery; they differ only in the face field.
const bool traceFromFlux =
    (slTraceVelocity == "projectedFlux") || traceReconstructU
 || traceReconstructMom;
if (traceFromFlux)
{
    Info<< "SL trace velocity: "
        << (traceReconstructU
              ? "reconstruct(interpolate(Uext)&Sf) -- NOT divergence-free"
                " (smoothing+extension control)"
              : traceReconstructMom
              ? "reconstruct(interpolate(U)&Sf) -- NOT divergence-free,"
                " NO extension (smoothing-only control)"
              : "reconstruct(phi) -- divergence-free face flux")
        << nl << endl;
}

// The single definition of the traced face->cell field, used at every refresh
// site so the two sources can never drift apart.
auto slTraceField = [&]() -> tmp<volVectorField>
{
    if (traceReconstructU)
    {
        return fvc::reconstruct(linearInterpolate(slVelExt->Uext()) & mesh.Sf());
    }
    if (traceReconstructMom)
    {
        // Raw momentum U, NOT the extension: identical operator chain to
        // reconstructedU, so the only difference between the two is whether the
        // field was extended off the interface.
        return fvc::reconstruct(linearInterpolate(U) & mesh.Sf());
    }
    return fvc::reconstruct(phi);
};

// Trace-velocity levels used when traceFromFlux; allocated only then. Utrace
// mirrors Uext's role (new level) and UtraceOld mirrors UextOld (u^n).
autoPtr<volVectorField> UtracePtr;
autoPtr<volVectorField> UtraceOldPtr;
if (traceFromFlux)
{
    UtracePtr.reset
    (
        new volVectorField("Utrace", slTraceField())
    );
    UtraceOldPtr.reset
    (
        new volVectorField("UtraceOld", UtracePtr())
    );
}

// Iterated interface-force coupling (PIMPLE.psiOuterCorrectors, default no):
// re-advect psi from psi^n on EVERY outer corrector with the CURRENT outer
// iterate of the velocity, re-running band -> phase indicator -> curvature ->
// face delivery per pass, so the capillary force converges to its t^{n+1}
// value within the step. The momentum LHS is fvm::ddt (implicit Euler, all at
// t^{n+1}); with the default `no` the sigma*kappa_f*snGrad(alpha) RHS is
// evaluated ONCE per step from the AB2-predicted trajectory and then FROZEN
// across the outer correctors -- the once-lagged explicit coupling MEASURED as
// the dt-proportional 90% of the parasitic growth rate at N=128 (r = r0 + c*dt,
// config/stationaryDropletDtSweep.yaml). With `yes`, pass 1 keeps the AB2
// predictor and every later pass advects psi^n by the actual momentum iterate
// (the kernel's documented (u^{n+1}, u^n) slots, no extrapolation), so at
// outer-loop convergence the step is backward-Euler consistent in ALL terms.
// Default no = advect once on the first outer corrector = byte-identical.
const Switch psiOuterCorrectors
(
    static_cast<const fvSolution&>(mesh).subOrEmptyDict("PIMPLE")
        .getOrDefault<Switch>("psiOuterCorrectors", false)
);

// psi^n snapshot for psiOuterCorrectors: the semi-Lagrangian update mutates
// psi in place (psi^n -> psi^{n+1}), so outer passes >= 2 must restore the
// old-time level before re-advecting -- otherwise the passes would compound.
// Registered NO_WRITE: an in-step working copy, never a result.
volScalarField psiN
(
    IOobject
    (
        "psiN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    psi
);

if (psiOuterCorrectors)
{
    Info<< "PIMPLE.psiOuterCorrectors ON: psi re-advected from psi^n on every"
        << " outer corrector with the current velocity iterate" << endl;
}

// Two-phase mass-flux models for div(rhoPhi,U):
//   interpolatedDensity  : geometric rho_c, interpolate(rho_c)*phi
//   geometricFaceDensity : geometric rho_c, geometric rho_f*phi
//   rhoLENT              : geometric rho_f*phi plus an auxiliary conservative
//                          rho equation on every PIMPLE outer iteration; rho_c
//                          is reset from Detrixhe-Aslam alpha after convergence.
const fvSolution& transportFvSolution = mesh;
const dictionary massFluxDict =
    transportFvSolution.subDict("levelSet").subOrEmptyDict("massFlux");
const word massFluxType =
    massFluxDict.getOrDefault<word>("type", "interpolatedDensity");

if
(
    massFluxType != "interpolatedDensity"
 && massFluxType != "geometricFaceDensity"
 && massFluxType != "rhoLENT"
)
{
    FatalIOErrorInFunction(massFluxDict)
        << "Unknown massFlux type '" << massFluxType << "'. Valid types are "
        << "interpolatedDensity, geometricFaceDensity and rhoLENT."
        << exit(FatalIOError);
}

// rhoLENT uses one geometrical area fraction for each face. A central average
// of the owner/neighbour planes is the symmetric approximation to the unique
// face fraction used in the paper. The older geometricFaceDensity experiment
// retains its upwind default.
const word faceInterpolation = massFluxDict.getOrDefault<word>
(
    "faceInterpolation",
    massFluxType == "rhoLENT" ? word("central") : word("upwind")
);
if (faceInterpolation != "upwind" && faceInterpolation != "central")
{
    FatalIOErrorInFunction(massFluxDict)
        << "Unknown massFlux.faceInterpolation '" << faceInterpolation
        << "'. Valid values are upwind and central."
        << exit(FatalIOError);
}
const bool massFluxUpwind = (faceInterpolation == "upwind");
const bool useRhoLENT = (massFluxType == "rhoLENT");
const bool freezeRhoPhi = !Foam::getEnv("SL_FREEZE_RHOPHI").empty();

// Retain the final outer-iteration diagnostics until the once-per-time-step
// CSV write. They are zero for comparison modes that do not solve rhoLENT.
scalar rhoMassResidualL1 = 0;
scalar rhoMassResidualLinf = 0;
scalar rhoMassResidualRelL1 = 0;
scalar rhoLENTResetL1 = 0;

Info<< "Two-phase mass flux (rhoPhi): " << massFluxType;
if (massFluxType != "interpolatedDensity")
{
    Info<< " (" << faceInterpolation << " geometric alpha_f)";
}
if (freezeRhoPhi)
{
    Info<< " [SL_FREEZE_RHOPHI: frozen diagnostic]";
}
Info<< nl << endl;

// Face liquid-area fraction and density at the NEW interface location. rhof is
// held geometrical during each time step while the PIMPLE volumetric flux is
// iterated, exactly as rho_f^{n+1} in rhoLENT Algorithm 1.
surfaceScalarField alphaf
(
    IOobject
    (
        "alphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaf", dimless, Zero)
);

surfaceScalarField rhof
(
    IOobject
    (
        "rhof",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rho)
);
