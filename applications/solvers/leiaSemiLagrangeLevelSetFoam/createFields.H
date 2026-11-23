// Signed distance / level set
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

// Phase-indicator
Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Volumetric flux (used only for the convective Courant number report /
// adjustTimeStep; the semi-Lagrangian update advects with the cell velocity U).
Info<< "Reading field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("phi", dimVolume / dimTime, 0)
);

// Velocity (the advecting field; semi-Lagrangian trajectories use it directly)
Info<< "Initializing field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector(0,0,0))
);

// Level Set narrowBand model
autoPtr<narrowBand> narrowBand = narrowBand::New(mesh, psi);
narrowBand->calc();
narrowBand->write();

// Level Set phase-indicator model
autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh);
phaseInd->calcPhaseIndicator(alpha, psi);
alpha.write();

// Advection velocity model (prescribed verification velocity)
autoPtr<velocityModel> velocityModel = velocityModel::New(mesh);

velocityModel->setVolumetricFlux(phi);
phi.write();

velocityModel->setVelocity(U);
U.write();

// Semi-Lagrangian advection engine (reconstruction selected from
// levelSet.semiLagrangian.reconstruction; default quadraticWeightedLeastSquares).
autoPtr<slAdvection> slAdv = slAdvection::New(mesh);

// ---- dynamic mesh (AMR line, STATUS 11.44) -------------------------------------
// Nothing below changes a run on a static mesh: `none` is the inert transfer, and
// the band and the reference profile are built for a dynamic mesh only.
const bool dynamicMesh = mesh.dynamic();
autoPtr<slMeshTransfer> meshTransfer = slMeshTransfer::New(mesh);
autoPtr<slRefinementBand> refinementBand;
// The exact reference of the error metrics on the PRESENT mesh: the profile of
// the implicit surface, evaluated at x for a flow that returns and at x - U T for
// `endReference translate` (the one-way translation). On a static mesh the
// reference fields psiEnd/alphaEnd of the pre-processing are used, as before.
autoPtr<implicitSurface> referenceSurface;
autoPtr<profile> referenceProfile;
vector referenceShift(Zero);
if (dynamicMesh)
{
    refinementBand.reset(new slRefinementBand(mesh));
    const dictionary& lsDict = mesh.solutionDict().subDict("levelSet");
    const dictionary& surfDict = lsDict.subDict("implicitSurface");
    referenceSurface = implicitSurface::New(surfDict.get<word>("type"), surfDict);
    referenceProfile = profile::New(mesh, referenceSurface.ref());
    const word endReference
    (
        mesh.solutionDict().getOrDefault<word>("endReference", "none")
    );
    if (endReference == "translate")
    {
        referenceShift =
            mesh.solutionDict().subDict("velocityModel").get<vector>("velocity")
           *runTime.endTime().value();
    }
    else if (endReference != "none")
    {
        FatalErrorInFunction
            << "endReference " << endReference << ": none | translate"
            << exit(FatalError);
    }
    Info<< "Dynamic mesh: meshTransfer " << meshTransfer->type()
        << ", reference shift " << referenceShift << nl << endl;
}
