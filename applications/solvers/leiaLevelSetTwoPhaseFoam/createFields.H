#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
bool viscosityFaceAlphaGeometric = true;
word viscosityFaceModelName("linear");

// FACE DYNAMIC VISCOSITY -- ONE token, SIX explicit values naming both halves.
//
//   alg_lin   alg_harm   alg_blend      alpha_f = <alpha_c>_f, linear interpolation
//   geo_lin   geo_harm   geo_blend      alpha_f = geometric face area fraction of the
//                                       reconstructed plane, the SAME one that builds rho_f
//
//   *_lin    mu_f   = alpha_f mu1 + (1-alpha_f) mu2
//   *_harm   1/mu_f = alpha_f/mu1 + (1-alpha_f)/mu2
//   *_blend  w mu_harm + (1-w) mu_lin,  w = |n.Sf|/|Sf|
//
// The grid is spelled out rather than expressed as two crossed tokens because an earlier
// single token with three values (interpolated / geometric / geometricHarmonic) silently
// mixed the two halves and left alg_harm untested -- so "harmonic is worst" had only ever
// been measured on the geometric alpha_f. Six names make the missing cell impossible to
// hide.
//
// alg_lin IS the framework default: the linear model commutes with interpolation, so
// <alpha mu1 + (1-alpha) mu2>_f = <alpha>_f mu1 + (1-<alpha>_f) mu2. That identity is what
// let the axes collapse unnoticed. It fails for *_harm and *_blend.
//
// *_blend carries NO free parameter. Across a layered medium the phases are in SERIES
// normal to the interface (harmonic) and in PARALLEL tangent to it (arithmetic), so the
// weight is the alignment of the interface normal with the FINITE-VOLUME FACE normal --
// geometry, not a tuned coefficient. Pure harmonic applies the series formula even on faces
// lying nearly IN the interface plane, where arithmetic is the correct average; that is a
// candidate reason for its poor showing.
{
    const dictionary& mfD =
        mesh.solutionDict().subOrEmptyDict("levelSet").subOrEmptyDict("massFlux");
    word v = mfD.getOrDefault<word>("viscosityFaceModel", "geo_lin");
    // legacy names from before the grid was spelled out
    if (v == "interpolated")           v = "alg_lin";
    else if (v == "geometric")         v = "geo_lin";
    else if (v == "geometricHarmonic") v = "geo_harm";

    const wordList valid
    ({"alg_lin","alg_harm","alg_blend","geo_lin","geo_harm","geo_blend"});
    if (!valid.found(v))
    {
        FatalErrorInFunction << "massFlux.viscosityFaceModel = " << v
            << " is not one of " << valid << exit(FatalError);
    }
    viscosityFaceAlphaGeometric = v.starts_with("geo_");
    viscosityFaceModelName =
        v.ends_with("_lin") ? word("linear")
      : v.ends_with("_harm") ? word("harmonic") : word("normalBlend");
    Info<< "face viscosity: " << v << "  (alpha_f "
        << (viscosityFaceAlphaGeometric ? "geometric" : "algebraic")
        << ", model " << viscosityFaceModelName << ')' << endl;
}

surfaceScalarField muf
(
    IOobject("muf", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(mixture.mu())
);

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

// Level Set Field 
Info<< "Reading field psi\n" << endl;
volScalarField psi 
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Leve Set re-distancing model
autoPtr<redistancer> redist = redistancer::New(mesh);

// SDPLS source (S = f_nl(psi^n)*psi, vanishes on the zero level set;
// default noSource = off): the psi-maintenance mechanism of the two-phase
// level-set transport (alphaEqn.H).
autoPtr<sdplsSource> source = sdplsSource::New(mesh);
source->update(psi, U);

// Level Set narrowBand model
autoPtr<narrowBand> narrowBand = narrowBand::New(mesh, psi); 
narrowBand->calc();
narrowBand->write();

// Level Set phase-indicator model
autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh); 

phaseInd->calcPhaseIndicator(alpha1, psi);
alpha1.write();

// Production curvature (SL reconstruction), constructed ONLY when the force
// below is reconstructedCurvature -- whose constructor looks `kappa` up, so
// these fields must be registered first.
#include "reconstructedCurvatureFields.H"

// Surface Tension Force Model 
autoPtr<surfaceTensionForce> fSigma = surfaceTensionForce::New(mesh); 
