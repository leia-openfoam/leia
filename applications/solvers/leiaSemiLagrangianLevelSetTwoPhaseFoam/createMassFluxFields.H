// THE TWO-PHASE MASS FLUX, shared by the semi-Lagrangian two-phase solver
// (createTransportFields.H) and the Eulerian two-phase solver
// (leiaLevelSetTwoPhaseFoam.C): the model selection, its switches and
// diagnostics, and the face fields alphaf, alphafOld and rhof. Moved here
// VERBATIM from createTransportFields.H on 2026-09-27, so that both solvers
// read the same dictionary with the same defaults.
//
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
if (alphaFSource != "donorPlane"
 && alphaFSource != "averagedPlanes"
 && alphaFSource != "donorPlaneAdvected")
{
    FatalIOErrorInFunction(massFluxDict)
        << "Unknown massFlux.alphaFSource '" << alphaFSource
        << "'. Valid values are donorPlane, donorPlaneAdvected and averagedPlanes "
        << "(deprecated: upwind, central)."
        << exit(FatalIOError);
}
// donorPlaneAdvected is donorPlane with the plane displaced by the SAME
// semi-Lagrangian step the cell values undergo, so it is a donor (upwind) choice too.
const bool massFluxUpwind =
    (alphaFSource == "donorPlane" || alphaFSource == "donorPlaneAdvected");

// SL-ADVECTED FACE FRACTION (author's proposal, 2026-09-02).
//
// The cell update is psi_c^{n+1} = psi^n(x_c - u dt), reconstructed from the quadratic
// fit. alpha_f, however, is built from a LINEAR plane fit of the ALREADY-UPDATED
// psi^{n+1} cut against the face: a different field, a different reconstruction, and no
// shared displacement. This applies the SAME operator to the face --
//
//     alpha_f^{n+1} = area fraction of { x in f : psi^n(x - u_f dt) < 0 }
//
// -- which for a plane is one dot product, d -> d - (n . u_f) dt, not a geometric
// construction. It is deliberately NOT swept-volume flux clipping: that is unsplit
// geometric VOF, a different and far more expensive method.
const bool alphaFAdvected = (alphaFSource == "donorPlaneAdvected");
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
