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

volScalarField psi0("psi0", psi);

// Phase-indicator  
volScalarField alpha 
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedScalar("alpha", dimless, 0)
);

// Narrow band model (registers "NarrowBand": band-based phase indicators
// look it up at CONSTRUCTION, the plane-anchored redistancer models at use).
autoPtr<narrowBand> nBand = narrowBand::New(mesh, psi);
nBand->calc();

// RTS
//- redistancer selection
autoPtr<redistancer> redist = redistancer::New(mesh);

//- phase-indicator selection
autoPtr<phaseIndicator> phaseInd = phaseIndicator::New(mesh);
