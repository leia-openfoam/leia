/*---------------------------------------------------------------------------*\
    leia :: davof :: davofQuadricCurvature (davofQuadricCurvature.H)
\*---------------------------------------------------------------------------*/

#include "davofQuadricCurvature.H"
#include "CPCCellToCellStencil.H"
#include "extendedCentredCellToCellStencil.H"
#include "calculatedFvPatchFields.H"
#include "cpuTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(davofQuadricCurvature, false);
    addToRunTimeSelectionTable
    (
        davofCurvature, davofQuadricCurvature, Dictionary
    );
}

namespace
{

// In-place Cholesky solve of the SPD system A x = b, A row-major n x n
// (n <= 6), b overwritten with x; false if A is not positive-definite. The
// routine of uncachedQuadraticWeightedLeastSquaresReconstruction.C.
static bool choleskySolve(Foam::scalar* A, Foam::scalar* b, const Foam::label n)
{
    using Foam::scalar;
    using Foam::label;
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= Foam::SMALL) { return false; }
        d = Foam::sqrt(d);
        A[j*n + j] = d;
        for (label i = j + 1; i < n; ++i)
        {
            scalar s = A[i*n + j];
            for (label k = 0; k < j; ++k) { s -= A[i*n + k]*A[j*n + k]; }
            A[i*n + j] = s/d;
        }
    }
    for (label i = 0; i < n; ++i)
    {
        scalar s = b[i];
        for (label k = 0; k < i; ++k) { s -= A[i*n + k]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label k = i + 1; k < n; ++k) { s -= A[k*n + i]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    return true;
}

} // anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofQuadricCurvature::davofQuadricCurvature
(
    const word& name,
    const fvMesh& mesh,
    const davofState& state,
    const dictionary& dict
)
:
    davofCurvature(name, mesh, state, dict),
    positionWeight_(dict.getOrDefault<scalar>("positionWeight", 0)),
    weighting_(dict.getOrDefault<word>("weighting", "area")),
    minNeighbours_(dict.getOrDefault<label>("minNeighbours", 5)),
    alignTol_(dict.getOrDefault<scalar>("alignTol", 0.2)),
    ridge_(dict.getOrDefault<scalar>("ridge", 0)),
    offsetBeta_(dict.getOrDefault<scalar>("offsetBeta", 0.25)),
    modelIndex_(mesh.nCells(), -1)
{
    if (weighting_ != "area" && weighting_ != "uniform")
    {
        FatalIOErrorInFunction(dict)
            << "quadricFit weighting '" << weighting_
            << "': valid are area, uniform." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::davofQuadricCurvature::frameOf
(
    const vector& e3,
    vector& e1,
    vector& e2
)
{
    // The coordinate axis least aligned with e3, made orthogonal to it.
    direction k = 0;
    for (direction i = 1; i < 3; ++i)
    {
        if (mag(e3[i]) < mag(e3[k])) k = i;
    }
    vector a(Zero);
    a[k] = 1;
    e1 = a - (a & e3)*e3;
    e1 /= max(mag(e1), VSMALL);
    e2 = e3 ^ e1;
}


void Foam::davofQuadricCurvature::curvatureOf
(
    const FixedList<scalar, 6>& C,
    const scalar hc,
    const scalar U,
    const scalar V,
    scalar& kappa,
    scalar& K
)
{
    const scalar fu = C[1] + C[3]*U + C[4]*V;
    const scalar fv = C[2] + C[4]*U + C[5]*V;
    const scalar fuu = C[3]/hc;
    const scalar fuv = C[4]/hc;
    const scalar fvv = C[5]/hc;
    const scalar g2 = 1 + fu*fu + fv*fv;
    kappa = -((1 + fv*fv)*fuu - 2*fu*fv*fuv + (1 + fu*fu)*fvv)/(g2*Foam::sqrt(g2));
    K = (fuu*fvv - fuv*fuv)/(g2*g2);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::davofQuadricCurvature::compute()
{
    cpuTime timer;
    const fvMesh& mesh = mesh_;
    const scalarField& V = mesh.V();
    const vectorField& m = state_.m().primitiveField();

    // The interface weight: |m_c| on interface cells, 0 elsewhere and on
    // non-coupled boundaries (their stencil entries are then skipped).
    volScalarField wI
    (
        IOobject
        (
            "davofCurvatureWeight", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh, dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );
    label nModels = 0;
    forAll(wI, c)
    {
        wI[c] = state_.isInterfaceCell(c) ? mag(m[c]) : 0;
        if (wI[c] > 0) modelIndex_[c] = nModels++;
        else modelIndex_[c] = -1;
    }
    origin_.setSize(nModels, Zero);
    e1_.setSize(nModels, Zero);
    e2_.setSize(nModels, Zero);
    e3_.setSize(nModels, Zero);
    hc_.setSize(nModels, 0);
    coeffs_.setSize(nModels);
    forAll(coeffs_, i) coeffs_[i] = FixedList<scalar, 6>(Zero);

    // The ring data, gathered once for the interface cells only (the lists
    // of the other cells are cleared before the map is built); remote cells
    // arrive through the mapDistribute.
    List<List<vector>> xSt, mSt;
    List<List<scalar>> wSt;
    {
        CPCCellToCellStencil cpc(mesh);
        forAll(cpc, c)
        {
            if (!(wI[c] > 0)) cpc[c].clear();
        }
        const extendedCentredCellToCellStencil so(cpc);
        so.collectData(state_.xPlane(), xSt);
        so.collectData(state_.m(), mSt);
        so.collectData(wI, wSt);
    }

    scalarField& kappaS = kappaS_.primitiveFieldRef();
    scalarField& KS = KS_.primitiveFieldRef();
    scalarField& k1 = kappa1S_.primitiveFieldRef();
    scalarField& k2 = kappa2S_.primitiveFieldRef();
    kappaS = 0; KS = 0; k1 = 0; k2 = 0;
    filled_ = false;
    nFallback_ = 0;

    const vectorField& xPl = state_.xPlane().primitiveField();

    forAll(wI, c)
    {
        if (!(wI[c] > 0)) continue;
        const label mi = modelIndex_[c];

        const scalar mc = mag(m[c]);
        const vector e3 = m[c]/mc;
        vector e1, e2;
        frameOf(e3, e1, e2);
        const scalar hc = Foam::cbrt(V[c]);
        const point o = xPl[c];

        origin_[mi] = o; e1_[mi] = e1; e2_[mi] = e2; e3_[mi] = e3; hc_[mi] = hc;
        FixedList<scalar, 6>& C = coeffs_[mi];
        C = FixedList<scalar, 6>(Zero);
        filled_[c] = true;

        // Normal equations A C = b over the 6 coefficients; the slope rows never
        // touch C0, so with positionWeight 0 the 5 x 5 block of C1..C5 is solved
        // alone and C0 follows from the positions.
        scalar A[36] = {0};
        scalar b[6] = {0};
        // For C0 with positionWeight 0
        scalar sw = 0, swW = 0;
        // Usable neighbours (other than c)
        label nUse = 0;

        const List<vector>& xs = xSt[c];
        const List<vector>& ms = mSt[c];
        const List<scalar>& ws = wSt[c];

        // Position rows kept for the C0 pass
        DynamicList<FixedList<scalar, 4>> posRows(xs.size());   // U, V, W, w

        forAll(xs, i)
        {
            if (!(ws[i] > 0)) continue;
            const scalar mj = mag(ms[i]);
            if (!(mj > 0)) continue;
            const vector nj = ms[i]/mj;
            const scalar a3 = nj & e3;
            if (a3 <= alignTol_) continue;
            const vector r = xs[i] - o;
            const scalar U = (r & e1)/hc;
            const scalar Vv = (r & e2)/hc;
            const scalar W = (r & e3)/hc;
            const bool self = (mag(r) < VSMALL && mag(nj - e3) < VSMALL);
            if (!self) ++nUse;

            const scalar w =
                (weighting_ == "area") ? Foam::sqrt(mj/mc) : scalar(1);
            const scalar w2 = w*w;
            const scalar s1 = -(nj & e1)/a3;
            const scalar s2 = -(nj & e2)/a3;

            // slope row 1: C1 + C3 U + C4 V = s1
            {
                const label idx[3] = {1, 3, 4};
                const scalar val[3] = {1, U, Vv};
                for (label p = 0; p < 3; ++p)
                {
                    b[idx[p]] += w2*val[p]*s1;
                    for (label q = 0; q < 3; ++q)
                    {
                        A[idx[p]*6 + idx[q]] += w2*val[p]*val[q];
                    }
                }
            }
            // slope row 2: C2 + C4 U + C5 V = s2
            {
                const label idx[3] = {2, 4, 5};
                const scalar val[3] = {1, U, Vv};
                for (label p = 0; p < 3; ++p)
                {
                    b[idx[p]] += w2*val[p]*s2;
                    for (label q = 0; q < 3; ++q)
                    {
                        A[idx[p]*6 + idx[q]] += w2*val[p]*val[q];
                    }
                }
            }
            // position row: C0 + C1 U + C2 V + C3 U^2/2 + C4 U V + C5 V^2/2 = W
            if (positionWeight_ > 0)
            {
                const scalar wp2 = sqr(positionWeight_*w);
                const scalar val[6] = {1, U, Vv, 0.5*U*U, U*Vv, 0.5*Vv*Vv};
                for (label p = 0; p < 6; ++p)
                {
                    b[p] += wp2*val[p]*W;
                    for (label q = 0; q < 6; ++q)
                    {
                        A[p*6 + q] += wp2*val[p]*val[q];
                    }
                }
            }
            FixedList<scalar, 4> pr;
            pr[0] = U; pr[1] = Vv; pr[2] = W; pr[3] = w2;
            posRows.append(pr);
        }

        if (nUse < minNeighbours_)
        {
            ++nFallback_;
            continue;   // the plane model with zero coefficients
        }

        bool ok = false;
        if (positionWeight_ > 0)
        {
            if (ridge_ > 0)
            {
                scalar dmax = 0;
                for (label p = 0; p < 6; ++p) dmax = max(dmax, A[p*6 + p]);
                for (label p = 0; p < 6; ++p) A[p*6 + p] += ridge_*dmax;
            }
            ok = choleskySolve(A, b, 6);
            if (ok)
            {
                for (label p = 0; p < 6; ++p) C[p] = b[p];
            }
        }
        else
        {
            // The 5 x 5 block of C1..C5
            scalar A5[25], b5[5];
            for (label p = 0; p < 5; ++p)
            {
                b5[p] = b[p + 1];
                for (label q = 0; q < 5; ++q) A5[p*5 + q] = A[(p + 1)*6 + (q + 1)];
            }
            if (ridge_ > 0)
            {
                scalar dmax = 0;
                for (label p = 0; p < 5; ++p) dmax = max(dmax, A5[p*5 + p]);
                for (label p = 0; p < 5; ++p) A5[p*5 + p] += ridge_*dmax;
            }
            ok = choleskySolve(A5, b5, 5);
            if (ok)
            {
                for (label p = 0; p < 5; ++p) C[p + 1] = b5[p];
                // C0: the weighted mean of the position residuals
                forAll(posRows, i)
                {
                    const FixedList<scalar, 4>& pr = posRows[i];
                    const scalar U = pr[0], Vv = pr[1], W = pr[2], w2 = pr[3];
                    const scalar rest =
                        C[1]*U + C[2]*Vv + 0.5*C[3]*U*U + C[4]*U*Vv + 0.5*C[5]*Vv*Vv;
                    sw += w2;
                    swW += w2*(W - rest);
                }
                C[0] = (sw > 0) ? swW/sw : 0;
            }
        }

        if (!ok)
        {
            C = FixedList<scalar, 6>(Zero);
            ++nFallback_;
            continue;
        }

        scalar kappa, K;
        curvatureOf(C, hc, 0, 0, kappa, K);
        kappaS[c] = kappa;
        KS[c] = K;
        const scalar disc = Foam::sqrt(max(kappa*kappa - 4*K, scalar(0)));
        k1[c] = 0.5*(kappa + disc);
        k2[c] = 0.5*(kappa - disc);
    }

    kappaS_.correctBoundaryConditions();
    KS_.correctBoundaryConditions();
    kappa1S_.correctBoundaryConditions();
    kappa2S_.correctBoundaryConditions();
    reduce(nFallback_, sumOp<label>());
    cpuSeconds_ = timer.cpuTimeIncrement();
}


bool Foam::davofQuadricCurvature::pack(const label c, davofModelPack& p) const
{
    const label mi = modelIndex_[c];
    if (mi < 0 || !filled_[c])
    {
        p = davofModelPack(Zero);
        return false;
    }
    p[0] = 1;
    for (direction i = 0; i < 3; ++i)
    {
        p[1 + i] = origin_[mi][i];
        p[4 + i] = e1_[mi][i];
        p[7 + i] = e2_[mi][i];
        p[10 + i] = e3_[mi][i];
    }
    p[13] = hc_[mi];
    for (label k = 0; k < 6; ++k) p[14 + k] = coeffs_[mi][k];
    return true;
}


bool Foam::davofQuadricCurvature::offsetPacked
(
    const davofModelPack& p,
    const point& x,
    scalar& d,
    scalar& kappa,
    scalar& K
) const
{
    if (!(p[0] > 0.5)) return false;

    const vector o(p[1], p[2], p[3]);
    const vector e1(p[4], p[5], p[6]);
    const vector e2(p[7], p[8], p[9]);
    const vector e3(p[10], p[11], p[12]);
    const scalar hc = p[13];
    FixedList<scalar, 6> C;
    for (label k = 0; k < 6; ++k) C[k] = p[14 + k];

    const vector r = x - o;
    const scalar U = (r & e1)/hc;
    const scalar V = (r & e2)/hc;
    const scalar W = (r & e3)/hc;

    const scalar F = C[0] + C[1]*U + C[2]*V + 0.5*C[3]*U*U + C[4]*U*V + 0.5*C[5]*V*V;
    const scalar FU = C[1] + C[3]*U + C[4]*V;
    const scalar FV = C[2] + C[4]*U + C[5]*V;
    const scalar psi = W - F;
    const scalar gm = Foam::sqrt(FU*FU + FV*FV + 1);
    const scalar nU = -FU/gm, nV = -FV/gm;
    // h_nn = n.H.n with H = -Hess(F) in the (U, V) block
    const scalar hnn = -(C[3]*nU*nU + 2*C[4]*nU*nV + C[5]*nV*nV);
    const scalar D = gm*gm - 2*psi*hnn;
    scalar s;
    if (D >= offsetBeta_*gm*gm)
    {
        s = 2*psi/(gm + Foam::sqrt(D));
    }
    else
    {
        s = psi/gm;
    }
    d = s*hc;
    // The ray foot: x - d n, i.e. (U, V) - s (nU, nV)
    curvatureOf(C, hc, U - s*nU, V - s*nV, kappa, K);
    return true;
}


// ************************************************************************* //
