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
// GEOMETRIC FACE DYNAMIC VISCOSITY, runtime-selectable and INERT by default.
//
// OpenFOAM's viscous term, turbulence->divDevRhoReff(rho,U), takes mu from the mixture as
// a CELL field, mu_c = alpha_c rho1 nu1 + (1-alpha_c) rho2 nu2, and fvm::laplacian then
// interpolates THAT linearly to the faces. So mu_f descends from a linearly interpolated
// cell alpha while the mass flux rho_f = alpha_f rho1 + (1-alpha_f) rho2 descends from the
// GEOMETRIC face area fraction of the reconstructed interface. The two face coefficients
// of the same momentum equation are therefore built from different discrete objects.
//
// With viscosityFaceModel = geometric the SL solver overwrites muf below with
//     muf = alphaf mu1 + (1-alphaf) mu2
// using exactly the alphaf that produced rho_f, and UEqn.H assembles the viscous term
// explicitly from it. With "interpolated" (the default) muf is never read and the
// standard divDevRhoReff path runs unchanged, so every existing study is bit-identical.
const word viscosityFaceModel =
    mesh.solutionDict().subOrEmptyDict("levelSet").subOrEmptyDict("massFlux")
        .getOrDefault<word>("viscosityFaceModel", "interpolated");
if (viscosityFaceModel != "interpolated"
 && viscosityFaceModel != "geometric"
 && viscosityFaceModel != "geometricHarmonic")
{
    FatalErrorInFunction
        << "massFlux.viscosityFaceModel = " << viscosityFaceModel
        << " is not one of: interpolated, geometric, geometricHarmonic"
        << exit(FatalError);
}
// geometric          mu_f = alpha_f mu1 + (1-alpha_f) mu2          (arithmetic)
// geometricHarmonic  1/mu_f = alpha_f/mu1 + (1-alpha_f)/mu2        (harmonic)
//
// Both use the geometric alpha_f. The distinction matters because rho_f and mu_f are not
// the same kind of coefficient: rho_f is fixed by CONSERVATION -- it must match the mass
// equation, so the sharp geometric value is required -- whereas mu_f is a DIFFUSION
// coefficient, and the face value that reproduces the exact flux of a layered medium
// normal to the interface is the HARMONIC mean, not the arithmetic one. MEASURED
// 2026-09-03: the arithmetic form made things 2.08x WORSE in L1 than the interpolated
// default, which is what sharpening the wrong average does.
const bool geometricFaceViscosity = (viscosityFaceModel != "interpolated");
const bool harmonicFaceViscosity  = (viscosityFaceModel == "geometricHarmonic");

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
