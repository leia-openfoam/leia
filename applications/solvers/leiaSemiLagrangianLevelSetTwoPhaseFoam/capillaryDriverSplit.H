// WP8.1 (docs/plan-curvature-stabilization.md): split the variation of the
// DELIVERED face curvature into the part that varies ACROSS the force support
// (along the interface normal) and the part that varies ALONG the interface.
//
// Only the variation of kappa_f produces spurious velocity: writing
// kappa_f = kappaBar + dKappa_f, the flux splits as
//
//   sigma kappa_f snGrad(alpha)|Sf| = snGrad(sigma kappaBar alpha)|Sf|
//                                   + sigma dKappa_f snGrad(alpha)|Sf|,
//
// and the first term is a face-normal difference of a cell field, which the
// pressure equation reproduces exactly. Of the second term, the component of
// the variation ALONG the interface normal is unphysical by construction (the
// exact curvature is constant along normals within the support), while the
// component ALONG the interface carries the physical tangential capillary
// response. Which of the two grows during the coupled run decides whether a
// delivery that is structurally free of across-support variation (one
// curvature per interface element, extended along its normal) can help.
//
// Method, per cell that owns at least two ACTIVE faces (|snGrad(alpha)| > 0):
//   1. kappa_f is recovered from the delivered flux itself,
//      Q_f = G_sigma,f/(snGrad(alpha)_f |Sf|) = sigma kappa_f,
//      so the measurement is valid for EVERY model that returns the integrated
//      face flux -- including the registered face-curvature path, where
//      kappa_f is not the interpolate of any cell field.
//   2. n_hat from the cell's linear least-squares interface plane (the same
//      reconstruction the phase indicator uses; halo-only, fit-free).
//   3. for every pair of active faces of that cell, the difference quotient
//         q = |Q_f1 - Q_f2| / |x_f1 - x_f2|
//      is assigned to the ACROSS-SUPPORT bin when the separation is aligned
//      with n_hat (|dx.n_hat|/|dx| > 0.7) and to the ALONG-INTERFACE bin when
//      it is orthogonal to it (< 0.3); pairs in between are not counted.
//   4. both are normalised by the mean |Q| over the active faces, so the
//      reported numbers are relative variation rates [1/m] and are directly
//      comparable with each other.
//
// Serial/parallel note: pairs are formed from a cell's own internal faces, so
// the measure is rank-local by construction; cells whose faces lie on a
// processor patch contribute fewer pairs (counted, not silently dropped).

#ifndef capillaryDriverSplit_H
#define capillaryDriverSplit_H

#include "levelSetPlaneReconstruction.H"

namespace Foam
{

struct capillaryDriverSplitNorms
{
    scalar acrossSupportL2 = 0;    // [1/m] relative variation rate
    scalar alongInterfaceL2 = 0;   // [1/m]
    scalar acrossSupportMax = 0;
    scalar alongInterfaceMax = 0;
    label nAcross = 0;             // pairs counted per bin
    label nAlong = 0;
};


inline capillaryDriverSplitNorms computeCapillaryDriverSplit
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& alpha,
    const surfaceScalarField& GSigma
)
{
    capillaryDriverSplitNorms out;

    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snIn = snA.primitiveField();
    const scalarField& gIn = GSigma.primitiveField();
    const scalarField& magSfIn = mesh.magSf().primitiveField();
    const vectorField& CfIn = mesh.Cf().primitiveField();
    const scalar snMax = max(gMax(mag(snIn)), VSMALL);

    // Q_f = sigma kappa_f on the active faces.
    const label nIf = mesh.nInternalFaces();
    scalarField Q(nIf, 0.0);
    boolList activeFace(nIf, false);
    scalar qSum = 0;
    label nAct = 0;
    for (label f = 0; f < nIf; ++f)
    {
        if (mag(snIn[f]) <= 1e-8*snMax) { continue; }
        activeFace[f] = true;
        Q[f] = gIn[f]/(snIn[f]*magSfIn[f]);
        qSum += mag(Q[f]);
        ++nAct;
    }
    scalar qSumG = qSum, nActG = nAct;
    reduce(qSumG, sumOp<scalar>());
    reduce(nActG, sumOp<scalar>());
    const scalar qBar = (nActG > 0) ? qSumG/nActG : 0;
    if (qBar <= VSMALL) { return out; }

    const coupledFaceNeighbours coupledNei(mesh, psi);
    scalar sAcross = 0, sAlong = 0;

    forAll(mesh.cells(), c)
    {
        const cell& cFaces = mesh.cells()[c];

        // the cell's own active internal faces
        DynamicList<label> act(6);
        forAll(cFaces, i)
        {
            const label f = cFaces[i];
            if (f < nIf && activeFace[f]) { act.append(f); }
        }
        if (act.size() < 2) { continue; }

        const scalarList pc = leastSquaresPlaneCoeffs(mesh, psi, c, coupledNei);
        vector nHat(pc[0], pc[1], pc[2]);
        const scalar nMag = mag(nHat);
        if (nMag < SMALL) { continue; }
        nHat /= nMag;

        for (label a = 0; a < act.size(); ++a)
        {
            for (label b = a + 1; b < act.size(); ++b)
            {
                const vector dx = CfIn[act[a]] - CfIn[act[b]];
                const scalar dmag = mag(dx);
                if (dmag < SMALL) { continue; }

                const scalar align = mag(nHat & dx)/dmag;
                const scalar q = mag(Q[act[a]] - Q[act[b]])/(dmag*qBar);

                if (align > 0.7)
                {
                    sAcross += q*q;
                    out.acrossSupportMax = max(out.acrossSupportMax, q);
                    ++out.nAcross;
                }
                else if (align < 0.3)
                {
                    sAlong += q*q;
                    out.alongInterfaceMax = max(out.alongInterfaceMax, q);
                    ++out.nAlong;
                }
            }
        }
    }

    reduce(sAcross, sumOp<scalar>());
    reduce(sAlong, sumOp<scalar>());
    reduce(out.nAcross, sumOp<label>());
    reduce(out.nAlong, sumOp<label>());
    reduce(out.acrossSupportMax, maxOp<scalar>());
    reduce(out.alongInterfaceMax, maxOp<scalar>());

    out.acrossSupportL2 =
        (out.nAcross > 0) ? Foam::sqrt(sAcross/out.nAcross) : 0;
    out.alongInterfaceL2 =
        (out.nAlong > 0) ? Foam::sqrt(sAlong/out.nAlong) : 0;

    return out;
}

} // End namespace Foam

#endif

// ************************************************************************* //
