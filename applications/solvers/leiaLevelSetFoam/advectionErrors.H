// The CSV schema this file writes, defined once so the header and the write
// block below cannot drift apart. Included HERE rather than from
// errorCalculation.H, which is included inside main().
#include "advectionErrorsCsv.H"

template<typename Ostream>
void reportErrors
(
    Ostream& os,
    const volScalarField& psi,
    const volScalarField& psi0,
    const volScalarField& alpha, 
    const volScalarField& alpha0,
    const surfaceScalarField& phi,
    scalar maxCoNum
)
{
    // Compute the level set advection error. 
    volScalarField ePhi = Foam::mag(psi - psi0);

    // TODO: Improve the error calculation - use the narrow band. 
    // Zero the level set advection error in the bulk.
    forAll(alpha, cellI)
    {
        // If cellI is a bulk cell 
        if ((alpha[cellI] == 0) || 
            (alpha[cellI] == 1)) 
        {
            // and the sign of the signed distance didn't switch,
            // after the advection reversal, so cellI is a bulk
            // cell at the final simulation time. 
            if (Foam::sign(psi[cellI] == Foam::sign(psi0[cellI]))) 
                ePhi[cellI] = 0; 
            // This avoids the problem of maximal errors reported
            // away from the interface.
        }
    }

    scalar lInfEpsi = Foam::gMax(ePhi);

    const fvMesh& mesh = alpha0.mesh(); 
    const auto& runTime = mesh.time();
    const auto& V = mesh.V();

    scalar alphaV0 = gSum((alpha0 * V)());

    // SIGNED first. A magnitude cannot distinguish systematic volume LOSS from
    // an error that oscillates in sign, and those are different failure modes
    // with different causes. Measured on the 2D steady vortex at t=T: the
    // sourceless baseline is negative at all seven resolutions (monotone loss,
    // -0.245 -> -0.0023) while SDPLS:R changes sign twice (-0.138, +0.0096,
    // +0.058, +0.049, +0.026, +0.0061, -0.00048) -- it GAINS volume over most
    // of the ladder. mag() of that produces exactly the non-monotone sequence
    // that reads as noise, including a dip to 0.0096 at N=45 that is a zero
    // crossing rather than a convergent point. The signed value is emitted
    // ALONGSIDE the magnitude, so every existing column keeps its meaning.
    scalar alphaVolErrSigned = gSum((alpha * V)()) - alphaV0;

    scalar alphaVolErrSignedRel = alphaVolErrSigned / alphaV0;

    scalar alphaVolErr = Foam::mag(alphaVolErrSigned);

    scalar alphaVolErrRel = alphaVolErr / alphaV0;

    scalar alphaGeomErr = 
        gSum((V * Foam::mag(alpha - alpha0))());

    scalar alphaGeomErrRel = alphaGeomErr / alphaV0;
    
    scalar alphaBoundErr = Foam::max
    (
        // Max overshoot
        // Max positive value over c for
        // alpha_c - 1 > 0, otherwise 0.
        Foam::max(Foam::gMax((alpha - 1)()), 0.),
        // Max undershoot
        // Max positive value over c for
        // -alpha_c > 0, otherwise 0.
        //
        // This branch was a copy of the overshoot branch -- it computed
        // gMax(alpha - 1) a second time, so alphaBoundErr was max(over, over)
        // and an UNDERSHOOT (alpha < 0) could never be reported by any study
        // using this solver. The comment above already said what it should be.
        Foam::max(Foam::gMax((-alpha)()), 0.)
    );


    // Unstructured FVM length approximation
    scalar deltaX = Foam::gMin(
        Foam::pow(
            mesh.deltaCoeffs().internalField(),
            -1
        )()
    );

    scalar maxSumVolFlux = gMax(mag(fvc::surfaceSum(phi))());

    // The main MPI process writes into the error file.
    if (Pstream::myProcNo() == 0)
    {
        os << runTime.timeOutputValue() << ","
            << deltaX << "," 
            << lInfEpsi << "," 
            << alphaVolErr << ","
            << alphaVolErrRel << ","
            << alphaVolErrSigned << ","
            << alphaVolErrSignedRel << ","
            << alphaGeomErr << ","
            << alphaGeomErrRel << ","
            << alphaBoundErr << "," 
            << runTime.elapsedCpuTime() << ","
            << runTime.elapsedClockTime() << "," 
            << maxCoNum << ","
            << maxSumVolFlux << "\n";
    }
}
