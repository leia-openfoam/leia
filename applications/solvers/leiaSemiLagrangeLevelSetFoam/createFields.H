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
