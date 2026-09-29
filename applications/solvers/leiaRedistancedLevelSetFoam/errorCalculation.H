// Initialize fields for error analysis
surfaceScalarField phi0 ("phi0", phi);
volVectorField U0("U0", U);

// Use psiEnd if present as field, otherwise copy psi
volScalarField psi0
(
    IOobject
    (
        "psiEnd",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    psi
);

// Use alphaEnd if present as field, otherwise copy alpha
volScalarField alpha0
(
    IOobject
    (
        "alphaEnd",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha
);
// For example psiEnd, alphaEnd get written in the 3Dtranslation test case by the Allrun scripts.
// leiaSetFiels is called in a different context for creating psiEnd and alphaEnd. See 3Dtranslation/system/init_End/ for further info.


// Open error file for output (same columns as leiaLevelSetFoam.csv /
// leiaSemiLagrangeLevelSetFoam.csv for cross-solver comparison).

// Say which reference the error metrics use (2026-09-29): psiEnd/alphaEnd when the case has
// them, otherwise the initial fields, which is right only for a flow that returns. The silent
// fallback hid a reversed translation for four weeks (STATUS 11.19). Every rank computes the
// sum (a collective), the master prints; the CSVs do not change.
{
    const scalar refDiff = gSumMag(psi0.primitiveField() - psi.primitiveField());
    Info<< "Error reference: "
        << (refDiff > 0
            ? "psiEnd/alphaEnd (the exact end state of the case)"
            : "the initial fields (valid only for a flow that returns to its start)")
        << nl << endl;
}

OFstream errorFile("leiaRedistancedLevelSetFoam.csv");

// CSV Header
errorFile << Foam::advectionErrorsCsvHeader();

#include "CourantNo.H"

// Report initial errors
reportErrors(errorFile, psi, psi0, alpha, alpha0, phi, CoNum);
