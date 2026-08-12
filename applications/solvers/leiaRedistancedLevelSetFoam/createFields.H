// Signed distance
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

// Volumetric flux
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

// Velocity
Info<< "Initializing field U (visualization)\n" << endl;
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

// Leve Set re-distancing model
autoPtr<redistancer> redist = redistancer::New(mesh);

// Level Set narrowBand model
autoPtr<narrowBand> narrowBand = narrowBand::New(mesh, psi);
narrowBand->calc();
narrowBand->write();

// Level Set phase-indicator model
autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh);
phaseInd->calcPhaseIndicator(alpha, psi);
alpha.write();

// Advection velocity model
autoPtr<velocityModel> velocityModel = velocityModel::New(mesh);

velocityModel->setVolumetricFlux(phi);
if (args.found("fluxCorrection"))
    correctFlux(phi);
phi.write();

velocityModel->setVelocity(U);
U.write();

// SDPLS Source term. MUST come after setVelocity: the source's strain field
// R = n . grad(U) . n is evaluated by update(), so constructing it against the
// still-zero U wrote an identically-zero R into 0/ and made the t=0 strain
// diagnostic useless. Only t=0 was affected -- every later update() sees the
// real U -- so no published result changes.
autoPtr<sdplsSource> source = sdplsSource::New(mesh);
source->update(psi, U);
source->write();
