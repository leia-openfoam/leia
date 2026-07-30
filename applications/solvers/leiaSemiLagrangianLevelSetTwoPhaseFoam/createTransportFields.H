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
