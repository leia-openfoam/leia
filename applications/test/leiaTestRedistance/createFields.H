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

