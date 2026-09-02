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

// WHICH RECONSTRUCTED PLANE GIVES THE FACE LIQUID-AREA FRACTION alpha_f.
//
// This is NOT an interpolation scheme and has nothing to do with the
// div(rhoPhi,U) discretisation (that is fvSchemes, MOMENTUM_DIV_SCHEME, and is
// separately Gauss upwind). Each internal face has TWO candidate area
// fractions -- one from the owner cell's linear least-squares plane and one
// from the neighbour's -- and this selects between them:
//
//   donorPlane      alpha_f from the UPWIND (donor) cell's plane: owner if
//                   phi >= 0 else neighbour. The phase composition of the fluid
//                   the flux actually carries through the face. DEFAULT.
//   averagedPlanes  alpha_f = 1/2 (a_owner + a_neighbour), the symmetric
//                   average of the two independent geometric fractions.
//
// rho_f = alpha_f rho1 + (1 - alpha_f) rho2, and rhoPhi = rho_f phi.
//
// RENAMED 2026-09-01. The values used to be `upwind` and `central`, which are
// OpenFOAM interpolation-scheme words and invited exactly the wrong reading --
// `central` is NOT OpenFOAM's `linear`, since nothing is interpolated between
// cell values; two independent GEOMETRIC area fractions are averaged. The old
// key `faceInterpolation` and the old values are still accepted, with a
// warning, so committed cases keep working.
//
// donorPlane is the default because it is what works: "only Gauss upwind worked
// in rhoLENT" (author), and every droplet case in cases/ was developed and
// validated with it.
word alphaFSource = massFluxDict.getOrDefault<word>
(
    "alphaFSource",
    massFluxDict.getOrDefault<word>("faceInterpolation", "donorPlane")
);
if (massFluxDict.found("faceInterpolation") && !massFluxDict.found("alphaFSource"))
{
    WarningInFunction
        << "massFlux.faceInterpolation is DEPRECATED and renamed to "
        << "massFlux.alphaFSource -- it selects which reconstructed plane gives"
        << " the face area fraction, not an interpolation scheme." << nl << endl;
}
if (alphaFSource == "upwind" || alphaFSource == "central")
{
    const word mapped = (alphaFSource == "upwind") ? "donorPlane" : "averagedPlanes";
    WarningInFunction
        << "massFlux alpha_f source '" << alphaFSource << "' is DEPRECATED; "
        << "use '" << mapped << "'." << nl << endl;
    alphaFSource = mapped;
}
if (alphaFSource != "donorPlane" && alphaFSource != "averagedPlanes")
{
    FatalIOErrorInFunction(massFluxDict)
        << "Unknown massFlux.alphaFSource '" << alphaFSource
        << "'. Valid values are donorPlane and averagedPlanes "
        << "(deprecated: upwind, central)."
        << exit(FatalIOError);
}
const bool massFluxUpwind = (alphaFSource == "donorPlane");
const bool useRhoLENT = (massFluxType == "rhoLENT");
const bool freezeRhoPhi = !Foam::getEnv("SL_FREEZE_RHOPHI").empty();

// Retain the final outer-iteration diagnostics until the once-per-time-step
// CSV write. They are zero for comparison modes that do not solve rhoLENT.
scalar rhoMassResidualL1 = 0;
scalar rhoMassResidualLinf = 0;
scalar rhoMassResidualRelL1 = 0;
scalar rhoLENTResetL1 = 0;

// --- density bounding (rhoLENT only) and the mass-residual diagnostic ---------
//
// boundRho: clip the auxiliary density to [rho2, rho1] right after the density
// solve, before momentum and pressure see it. rho1/rho2 are constants, so a rho
// outside that interval is unphysical, and a NEGATIVE one flips the sign of 1/rho
// in the pressure Laplacian and ends the run in a GAMG blow-up. The clip is not
// free -- it breaks the discrete identity ddt(rho) + div(rhoPhi) = 0 that rhoLENT
// exists to enforce -- so it is MEASURED, never silent: rhoClipL1 and
// rhoClipFraction go to the metrics CSV.
//
// massResidualDiagnostic: compute ddt(rho) + div(rhoPhi) for EVERY mass-flux model,
// not only rhoLENT. For a droplet translating at uniform U0 the spurious momentum
// source is U0 times this residual and it vanishes identically at U0 = 0 -- the
// stationary/translating asymmetry. geometricFaceDensity solves no density equation,
// so its residual has never been measured. Pure diagnostic; changes no field.
const bool boundRhoLENT =
    massFluxDict.getOrDefault<bool>("boundRho", false);
const bool massResidualDiagnostic =
    massFluxDict.getOrDefault<bool>("massResidualDiagnostic", false);

// TIME LEVEL OF THE GEOMETRIC FACE DENSITY (author's proposal, 2026-09-02).
//   new        rho_f from the NEW interface only -- what every run to date used.
//   trapezoid  1/2 (rho_f^n + rho_f^{n+1}), both from the geometric reconstruction.
//
// WHY IT IS THE RIGHT OBJECT. The density that belongs in the mass flux is the
// density of the fluid that actually CROSSES the face during the step -- a time
// average over the swept region -- not the plane value at t^{n+1}. The trapezoid of
// the two endpoint reconstructions is a second-order approximation to that swept
// average, and it costs one extra plane reconstruction per outer iteration and no
// solve.
//
// WHAT IT CANNOT DO, stated up front so the measurement is read honestly: it cannot
// drive the discrete mass residual to zero. The balance is ONE CONSTRAINT PER CELL
// while a two-level blend supplies one value per face from a POINTWISE rule; with
// more faces than cells a satisfying face field exists, but no local formula reaches
// it -- that needs a projection. So the pre-registered expectation is a substantial
// REDUCTION of ddt(rho)+div(rhoPhi), not machine zero, and the size of that reduction
// is what decides whether the projection is worth building.
//
// Default `new` = every existing study is byte-unchanged.
const word alphaFTimeLevel =
    massFluxDict.getOrDefault<word>("alphaFTimeLevel", "new");
if (alphaFTimeLevel != "new" && alphaFTimeLevel != "trapezoid")
{
    FatalIOErrorInFunction(massFluxDict)
        << "Unknown alphaFTimeLevel '" << alphaFTimeLevel
        << "'. Valid: new, trapezoid." << exit(FatalIOError);
}
const bool alphaFTrapezoid = (alphaFTimeLevel == "trapezoid");

// MASS-FLUX PROJECTION (author's proposal, 2026-09-02). Correct the geometric mass
// flux so it satisfies the discrete mass balance exactly, instead of solving for an
// auxiliary density as rhoLENT does. rho stays purely geometric -- bounded in
// [rho2, rho1] BY CONSTRUCTION, no auxiliary equation, no clip -- and the imbalance is
// carried off by a curl-free flux correction rhoPhi' = snGrad(lambda)*magSf, with
// laplacian(lambda) = ddt(rho) + div(rhoPhi). Costs one Poisson solve of exactly the
// pressure equation's character per outer iteration.
//
// Requires a `massFluxPotential` entry in the case's solvers dict. Default false, so
// every existing study is byte-unchanged and no case needs that entry unless it opts in.
const bool projectMassFlux =
    massFluxDict.getOrDefault<bool>("projectMassFlux", false);
if (projectMassFlux && useRhoLENT)
{
    Info<< "massFlux: projectMassFlux with rhoLENT -- the auxiliary density equation"
        << " already drives the residual to round-off, so the projection is"
        << " redundant here. Running it anyway, as selected." << nl << endl;
}
scalar massProjResidualBefore = 0;
scalar massProjResidualAfter = 0;
scalar massProjGlobalDrift = 0;
scalar massProjCorrFraction = 0;

// Potential for the mass-flux correction. lambda is to the mass flux what pressure is
// to the volumetric flux; it has no physical identity, and its dimensions follow from
// laplacian(lambda) = ddt(rho) + div(rhoPhi): kg/(m s). Constructed unconditionally
// (a zeroed volScalarField costs nothing) and solved only when the projection is on.
// zeroGradient everywhere: the correction vanishes where the flux is already right.
volScalarField massFluxPotential
(
    IOobject
    (
        "massFluxPotential",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "massFluxPotential",
        dimMass/(dimLength*dimTime),
        Zero
    ),
    zeroGradientFvPatchScalarField::typeName
);
if (boundRhoLENT && !useRhoLENT)
{
    Info<< "massFlux: boundRho requested but type is '" << massFluxType
        << "' -- there is no auxiliary density equation to bound; ignored."
        << nl << endl;
}
scalar rhoClipL1 = 0;
scalar rhoClipFraction = 0;

Info<< "Two-phase mass flux (rhoPhi): " << massFluxType;
if (massFluxType != "interpolatedDensity")
{
    Info<< " (alpha_f from " << alphaFSource << ")";
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

// Companion face fraction at the OLD interface, for alphaFTimeLevel trapezoid. It is
// constructed unconditionally (a surfaceScalarField of zeros costs nothing) and only
// FILLED when the trapezoid form is selected, so the `new` default touches nothing.
surfaceScalarField alphafOld
(
    IOobject
    (
        "alphafOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphafOld", dimless, Zero)
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
