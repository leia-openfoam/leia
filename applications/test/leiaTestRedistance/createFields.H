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

// Phase indicator (volume fraction): pre/post-event volume and shape
// metrics of the redistance event.
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

// Narrow band: the plane-anchored redistancer models look up the
// registered "NarrowBand" field to select their anchor cells (and the
// band-based phase indicators read it at construction).
autoPtr<narrowBand> nBand = narrowBand::New(mesh, psi);
nBand->calc();

// Level Set phase-indicator model
autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh);

// Leve Set re-distancing model
autoPtr<redistancer> redist = redistancer::New(mesh);

// // Volumetric flux
// Info<< "Reading field phi\n" << endl;
// surfaceScalarField phi 
// (
//     IOobject
//     (
//         "phi",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar ("phi", dimVolume / dimTime, 0)
// );

// Velocity 
// volVectorField U 
// (
//     IOobject
//     (
//         "U",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("U", dimVelocity, vector(1,0,0))
// );

