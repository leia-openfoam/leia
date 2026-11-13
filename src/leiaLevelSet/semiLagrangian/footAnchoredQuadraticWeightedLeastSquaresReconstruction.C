#include "footAnchoredQuadraticWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(footAnchoredQuadraticWeightedLeastSquaresReconstruction, 0);

    addToRunTimeSelectionTable
    (
        slReconstruction,
        footAnchoredQuadraticWeightedLeastSquaresReconstruction,
        Mesh
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- In-place Cholesky solve of a small symmetric positive-definite system.
//  A duplicate of the file-static helper of
//  uncachedQuadraticWeightedLeastSquaresReconstruction.C, which has internal
//  linkage there. Kept identical so the two models differ ONLY in the weight
//  anchor; a difference in the linear algebra would confound the comparison.
bool choleskySolveLocal(Foam::scalar* A, Foam::scalar* b, const Foam::label n)
{
    for (Foam::label i = 0; i < n; ++i)
    {
        for (Foam::label j = 0; j <= i; ++j)
        {
            Foam::scalar sum = A[i*n + j];
            for (Foam::label k = 0; k < j; ++k)
            {
                sum -= A[i*n + k]*A[j*n + k];
            }
            if (i == j)
            {
                if (sum <= Foam::SMALL) { return false; }
                A[i*n + j] = Foam::sqrt(sum);
            }
            else
            {
                A[i*n + j] = sum/A[j*n + j];
            }
        }
    }
    for (Foam::label i = 0; i < n; ++i)
    {
        Foam::scalar sum = b[i];
        for (Foam::label k = 0; k < i; ++k) { sum -= A[i*n + k]*b[k]; }
        b[i] = sum/A[i*n + i];
    }
    for (Foam::label i = n - 1; i >= 0; --i)
    {
        Foam::scalar sum = b[i];
        for (Foam::label k = i + 1; k < n; ++k) { sum -= A[k*n + i]*b[k]; }
        b[i] = sum/A[i*n + i];
    }
    return true;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::footAnchoredQuadraticWeightedLeastSquaresReconstruction::
footAnchoredQuadraticWeightedLeastSquaresReconstruction(const fvMesh& mesh)
:
    uncachedQuadraticWeightedLeastSquaresReconstruction(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::footAnchoredQuadraticWeightedLeastSquaresReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];
    const label nc = ncoeff_[c];

    // Constant reconstruction: no fit to re-anchor.
    if (nc == 0) { return psiC; }

    const List<scalar>& s = stencilPsi_[c];
    const label nNbr = nNbr_[c];
    if (nNbr < nc)
    {
        return uncachedQuadraticWeightedLeastSquaresReconstruction::
            evaluateRaw(c, x);
    }

    const point xc = stencilC(c, 0);        // the polynomial stays anchored here

    // Normal equations M = A^T W^2 A, g = A^T W^2 ds, with the weights measured
    // from the EVALUATION point x rather than from xc. The basis argument stays
    // (x_i - xc): only the weight moves, so the fit is still of the differences
    // about the cell centre and R_c(xc) = psi_c is preserved.
    scalar A[81];
    scalar g[9];
    scalar brow[9];

    for (label k = 0; k < nc; ++k)
    {
        g[k] = 0;
        for (label l = 0; l < nc; ++l) { A[k*nc + l] = 0; }
    }

    for (label i = 0; i < nNbr; ++i)
    {
        const point xi = stencilC(c, i + 1);
        const vector d = xi - xc;                  // basis displacement
        const scalar w = 1.0/Foam::max(Foam::mag(xi - x), SMALL);   // FOOT-anchored
        const scalar w2 = w*w;
        basis(d, nc, brow);
        const scalar ds = s[i + 1] - psiC;
        for (label k = 0; k < nc; ++k)
        {
            const scalar wb = w2*brow[k];
            g[k] += wb*ds;
            for (label l = 0; l < nc; ++l) { A[k*nc + l] += wb*brow[l]; }
        }
    }

    if (ridgeEps_ > 0)
    {
        scalar tr = 0;
        for (label k = 0; k < nc; ++k) { tr += A[k*nc + k]; }
        const scalar r = ridgeEps_*tr/nc;
        for (label k = 0; k < nc; ++k) { A[k*nc + k] += r; }
    }

    if (!choleskySolveLocal(A, g, nc))
    {
        // A degenerate foot-anchored system must never be worse than production.
        return uncachedQuadraticWeightedLeastSquaresReconstruction::
            evaluateRaw(c, x);
    }

    const vector dx = x - xc;
    scalar b[9];
    basis(dx, nc, b);
    scalar val = psiC;
    for (label k = 0; k < nc; ++k) { val += g[k]*b[k]; }
    return val;
}


// ************************************************************************* //
