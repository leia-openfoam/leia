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

// FACE DYNAMIC VISCOSITY. Two ORTHOGONAL choices, which an earlier single token conflated.
//
//   viscosityFaceAlpha : where alpha_f comes from
//       algebraic   <alpha_c>_f, linear interpolation of the cell phase indicator
//       geometric   the face area fraction of the reconstructed interface plane -- the
//                   SAME alpha_f that builds rho_f for the mass flux
//
//   viscosityFaceModel : how mu_f is formed from alpha_f
//       linear      mu_f   = alpha_f mu1 + (1-alpha_f) mu2
//       harmonic    1/mu_f = alpha_f/mu1 + (1-alpha_f)/mu2
//       normalBlend w mu_harmonic + (1-w) mu_linear,  w = |n.Sf|/|Sf|
//
// WHY normalBlend HAS NO FREE PARAMETER. Across a layered medium the correct face value
// depends on the orientation of the layering relative to the flux: NORMAL to the interface
// the two phases are in SERIES, which is the harmonic mean; TANGENT to it they conduct in
// PARALLEL, which is the arithmetic mean. The weight is therefore the alignment of the
// interface normal with the face normal, a geometric quantity, not a tuned coefficient.
//
// The old single-valued token still works and maps onto the grid:
//   interpolated -> (algebraic, linear);  geometric -> (geometric, linear);
//   geometricHarmonic -> (geometric, harmonic).
{
    const dictionary& mfD =
        mesh.solutionDict().subOrEmptyDict("levelSet").subOrEmptyDict("massFlux");
    word aSrc = mfD.getOrDefault<word>("viscosityFaceAlpha", "none");
    word mMod = mfD.getOrDefault<word>("viscosityFaceModel", "none");
    if (aSrc == "none" || mMod == "none")
    {
        // legacy single-token form
        const word legacy = (mMod != "none") ? mMod : word("geometric");
        if (legacy == "interpolated")           { aSrc = "algebraic"; mMod = "linear"; }
        else if (legacy == "geometric")         { aSrc = "geometric"; mMod = "linear"; }
        else if (legacy == "geometricHarmonic") { aSrc = "geometric"; mMod = "harmonic"; }
        else if (aSrc == "none")                { aSrc = "geometric"; }
    }
    if (aSrc != "algebraic" && aSrc != "geometric")
    {
        FatalErrorInFunction << "massFlux.viscosityFaceAlpha = " << aSrc
            << " is not one of: algebraic, geometric" << exit(FatalError);
    }
    if (mMod != "linear" && mMod != "harmonic" && mMod != "normalBlend")
    {
        FatalErrorInFunction << "massFlux.viscosityFaceModel = " << mMod
            << " is not one of: linear, harmonic, normalBlend" << exit(FatalError);
    }
    viscosityFaceAlphaGeometric = (aSrc == "geometric");
    viscosityFaceModelName = mMod;
}
Info<< "face viscosity: alpha_f "
    << (viscosityFaceAlphaGeometric ? "geometric" : "algebraic")
    << ", model " << viscosityFaceModelName << endl;

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
