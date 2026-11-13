/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

\*---------------------------------------------------------------------------*/

#include "ghostCellQuadraticWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "cellModel.H"
#include "labelIOList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        ghostCellQuadraticWeightedLeastSquaresReconstruction,
        0
    );
    addToRunTimeSelectionTable
    (
        slReconstruction,
        ghostCellQuadraticWeightedLeastSquaresReconstruction,
        Mesh
    );
}

namespace
{
// In-place Cholesky solve of the SPD system A x = b (row-major n x n, n <= 9).
// Returns false if A is not positive definite.
static bool choleskySolve(Foam::scalar* A, Foam::scalar* b, const Foam::label n)
{
    using Foam::scalar;
    using Foam::label;
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= 1e-10) { return false; }       // A is Jacobi-scaled: diag = 1
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

// Weighted, anchored least-squares solve for the first n coefficients.
// B: m rows of nFull basis values, w2: squared weights, r: value - centre value.
static bool solveFit
(
    const Foam::scalar* B,
    const Foam::scalar* w2,
    const Foam::scalar* r,
    const Foam::label m,
    const Foam::label nFull,
    const Foam::label n,
    Foam::scalar* x
)
{
    using Foam::scalar;
    using Foam::label;
    if (m < n) { return false; }
    scalar M[81];
    scalar s[9];
    for (label a = 0; a < n; ++a)
    {
        x[a] = 0;
        for (label b = 0; b < n; ++b) { M[a*n + b] = 0; }
    }
    for (label i = 0; i < m; ++i)
    {
        const scalar* bi = B + i*nFull;
        for (label a = 0; a < n; ++a)
        {
            x[a] += w2[i]*bi[a]*r[i];
            for (label b = 0; b <= a; ++b) { M[a*n + b] += w2[i]*bi[a]*bi[b]; }
        }
    }
    for (label a = 0; a < n; ++a)
    {
        if (M[a*n + a] <= Foam::VSMALL) { return false; }
        s[a] = 1.0/Foam::sqrt(M[a*n + a]);
    }
    for (label a = 0; a < n; ++a)
    {
        x[a] *= s[a];
        for (label b = 0; b <= a; ++b)
        {
            M[a*n + b] *= s[a]*s[b];
            M[b*n + a] = M[a*n + b];
        }
    }
    if (!choleskySolve(M, x, n)) { return false; }
    for (label a = 0; a < n; ++a) { x[a] *= s[a]; }
    return true;
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::basis
(
    const vector& d,
    scalar* b
) const
{
    const label nd = activeDirs_.size();
    label k = 0;
    for (label a = 0; a < nd; ++a) { b[k++] = d[activeDirs_[a]]; }
    for (label a = 0; a < nd; ++a)
    {
        b[k++] = 0.5*d[activeDirs_[a]]*d[activeDirs_[a]];
    }
    for (label a = 0; a < nd; ++a)
    {
        for (label c = a + 1; c < nd; ++c)
        {
            b[k++] = d[activeDirs_[a]]*d[activeDirs_[c]];
        }
    }
}


Foam::scalar Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::poly
(
    const tensor& t,
    const scalar psiC,
    const vector& d
) const
{
    scalar b[9];
    basis(d, b);
    scalar v = psiC;
    for (label k = 0; k < ncoeff_; ++k) { v += t[k]*b[k]; }
    return v;
}


void Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::build()
{
    const label nCells = mesh_.nCells();

    // -- the refinement levels that refineHexMesh writes
    labelIOList cellLevel
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        labelList(nCells, Zero)
    );
    labelIOList pointLevel
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh_.nPoints(), Zero)
    );
    if (cellLevel.size() != nCells || pointLevel.size() != mesh_.nPoints())
    {
        FatalErrorInFunction
            << "polyMesh/cellLevel or polyMesh/pointLevel does not have the"
            << " size of the mesh" << exit(FatalError);
    }
    cellLevel_ = cellLevel;
    maxLevel_ = returnReduce(cellLevel_.size() ? max(cellLevel_) : 0, maxOp<label>());

    // -- the parent centre point of every refined cell: among the points of the
    //    cell at the cell's own level, the one that is farthest from the corner
    //    that the cell has from its parent (the point of the lowest level)
    const pointField& pts = mesh_.points();
    const labelListList& cp = mesh_.cellPoints();
    parentPoint_.setSize(nCells);
    parentPoint_ = -1;
    forAll(cellLevel_, c)
    {
        const label L = cellLevel_[c];
        if (L == 0) { continue; }
        label p0 = -1;
        forAll(cp[c], i)
        {
            const label p = cp[c][i];
            if (pointLevel[p] < L && (p0 < 0 || pointLevel[p] < pointLevel[p0]))
            {
                p0 = p;
            }
        }
        if (p0 < 0) { continue; }
        scalar dMax = -1;
        forAll(cp[c], i)
        {
            const label p = cp[c][i];
            if (pointLevel[p] == L)
            {
                const scalar d = magSqr(pts[p] - pts[p0]);
                if (d > dMax) { dMax = d; parentPoint_[c] = p; }
            }
        }
    }
    parentCount_.setSize(mesh_.nPoints());
    parentCount_ = 0;
    forAll(parentPoint_, c)
    {
        if (parentPoint_[c] >= 0) { parentCount_[parentPoint_[c]] += 1; }
    }
    syncTools::syncPointList(mesh_, parentCount_, plusEqOp<scalar>(), scalar(0));

    // -- the lattice of every hexahedron: the mean of its four parallel edges
    lattice_.setSize(nCells);
    lattice_ = Zero;
    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellShapeList& shapes = mesh_.cellShapes();
    static const label e[3][4][2] =
    {
        {{0, 1}, {3, 2}, {4, 5}, {7, 6}},
        {{0, 3}, {1, 2}, {4, 7}, {5, 6}},
        {{0, 4}, {1, 5}, {2, 6}, {3, 7}}
    };
    forAll(shapes, c)
    {
        if (shapes[c].model() != hex) { continue; }
        const cellShape& sh = shapes[c];
        vector a[3];
        for (label k = 0; k < 3; ++k)
        {
            a[k] = Zero;
            for (label j = 0; j < 4; ++j)
            {
                a[k] += 0.25*(pts[sh[e[k][j][1]]] - pts[sh[e[k][j][0]]]);
            }
        }
        lattice_[c] = tensor(a[0], a[1], a[2]);
    }

    // -- the static fields of the stencil entries
    forAll(cellLevel_, c)
    {
        levelFld_.primitiveFieldRef()[c] = cellLevel_[c];
        parentCentreFld_.primitiveFieldRef()[c] =
            parentPoint_[c] >= 0 ? pts[parentPoint_[c]] : mesh_.C()[c];
    }
    levelFld_.correctBoundaryConditions();
    parentCentreFld_.correctBoundaryConditions();
    stencil_.collectData(levelFld_, stLevel_);
    stencil_.collectData(parentCentreFld_, stParentC_);

    label nNoParent = 0, nNoHex = 0;
    forAll(cellLevel_, c)
    {
        if (cellLevel_[c] > 0 && parentPoint_[c] < 0) { ++nNoParent; }
        if (magSqr(lattice_[c]) < VSMALL) { ++nNoHex; }
    }
    reduce(nNoParent, sumOp<label>());
    reduce(nNoHex, sumOp<label>());
    Info<< "ghostCellQuadraticWeightedLeastSquares: " << maxLevel_ + 1
        << " refinement level(s); refined cells without a parent point: "
        << nNoParent << "; cells that are not a hexahedron: " << nNoHex
        << endl;
    built_ = true;
}


void Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::fitCell
(
    const label c
)
{
    const label n = stencilSize(c);
    const List<scalar>& s = stencilPsi_[c];
    const point x0 = stencilC(c, 0);
    const scalar psiC = s[0];
    const label Lc = cellLevel_[c];
    const scalar tol2 = sqr(1e-3*radius_[c]);

    DynamicList<vector> X(2*n);
    DynamicList<scalar> V(2*n);
    DynamicList<label> coarse(n);
    DynamicList<vector> same(n);

    // -- data at the own level, and one datum per refined parent
    for (label i = 1; i < n; ++i)
    {
        const label r = raw(c, i);
        const point xi = stencilC(c, i);
        if (stencilIsBoundaryFace(c, i))
        {
            X.append(xi); V.append(s[i]);
            continue;
        }
        const label Li = label(stLevel_[c][r] + 0.5);
        if (Li == Lc)
        {
            X.append(xi); V.append(s[i]); same.append(xi);
        }
        else if (Li > Lc)
        {
            const point& xp = stParentC_[c][r];
            bool have = false;
            forAll(X, j)
            {
                if (magSqr(X[j] - xp) < tol2) { have = true; break; }
            }
            if (!have)
            {
                X.append(xp); V.append(stParentMean_[c][r]); ++nRestricted_;
            }
        }
        else
        {
            coarse.append(i);
        }
    }

    // -- ghost data in the coarser cells: the positions of the cell's own
    //    lattice that no datum occupies and that lie in a coarser stencil cell
    if (coarse.size() && magSqr(lattice_[c]) > VSMALL)
    {
        const vector a0 = lattice_[c].x();
        const vector a1 = lattice_[c].y();
        const vector a2 = lattice_[c].z();
        const scalar occ2 = 0.09*min(magSqr(a0), min(magSqr(a1), magSqr(a2)));
        const Vector<label>& sd = mesh_.solutionD();
        for (label i0 = -1; i0 <= 1; ++i0)
        for (label i1 = -1; i1 <= 1; ++i1)
        for (label i2 = -1; i2 <= 1; ++i2)
        {
            if (!i0 && !i1 && !i2) { continue; }
            const vector d = scalar(i0)*a0 + scalar(i1)*a1 + scalar(i2)*a2;
            bool inPlane = true;
            for (direction k = 0; k < 3; ++k)
            {
                if (sd[k] < 0 && mag(d[k]) > 1e-6*mag(d)) { inPlane = false; }
            }
            if (!inPlane) { continue; }
            const point q = x0 + d;
            bool occupied = false;
            forAll(X, j)
            {
                if (magSqr(X[j] - q) < occ2) { occupied = true; break; }
            }
            if (occupied) { continue; }
            label iA = -1;
            scalar dA = GREAT;
            forAll(coarse, j)
            {
                const scalar dj = magSqr(stencilC(c, coarse[j]) - q);
                if (dj < dA) { dA = dj; iA = coarse[j]; }
            }
            const point xA = stencilC(c, iA);
            if (dA >= magSqr(xA - x0)) { continue; }     // q is not in that cell
            X.append(q);
            V.append(poly(stCoeff_[c][raw(c, iA)], s[iA], q - xA));
            ++nGhost_;
        }
    }
    else if (coarse.size())
    {
        // not a hexahedron: keep the coarse cells as they are
        forAll(coarse, j)
        {
            X.append(stencilC(c, coarse[j])); V.append(s[coarse[j]]);
        }
    }

    // -- the anchored weighted least-squares fit, weights 1/distance
    const label m = X.size();
    scalarField B(m*ncoeff_), w2(m), rhs(m);
    forAll(X, i)
    {
        const vector d = X[i] - x0;
        basis(d, B.begin() + i*ncoeff_);
        w2[i] = 1.0/max(magSqr(d), VSMALL);
        rhs[i] = V[i] - psiC;
    }
    scalar x[9];
    tensor t(Zero);
    const label nd = activeDirs_.size();
    if (solveFit(B.begin(), w2.begin(), rhs.begin(), m, ncoeff_, ncoeff_, x))
    {
        for (label k = 0; k < ncoeff_; ++k) { t[k] = x[k]; }
    }
    else if (solveFit(B.begin(), w2.begin(), rhs.begin(), m, ncoeff_, nd, x))
    {
        for (label k = 0; k < nd; ++k) { t[k] = x[k]; }
        ++nLinear_;
    }
    else
    {
        ++nConstant_;
    }
    coeffFld_.primitiveFieldRef()[c] = t;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::
ghostCellQuadraticWeightedLeastSquaresReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    activeDirs_(),
    ncoeff_(0),
    cellLevel_(),
    maxLevel_(0),
    parentPoint_(),
    parentCount_(),
    lattice_(),
    levelFld_
    (
        IOobject
        (
            "slGhostLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    parentCentreFld_
    (
        IOobject
        (
            "slGhostParentCentre",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimLength, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    parentMeanFld_
    (
        IOobject
        (
            "slGhostParentMean",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    coeffFld_
    (
        IOobject
        (
            "slGhostCoeffs",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedTensor(dimless, Zero),
        zeroGradientFvPatchTensorField::typeName
    ),
    built_(false),
    nGhost_(0),
    nRestricted_(0),
    nLinear_(0),
    nConstant_(0),
    reported_(false)
{
    const Vector<label> gd = mesh_.geometricD();
    DynamicList<label> dirs(3);
    for (direction k = 0; k < 3; ++k)
    {
        if (gd[k] > 0) { dirs.append(k); }
    }
    activeDirs_ = dirs;
    const label nd = activeDirs_.size();
    ncoeff_ = 2*nd + nd*(nd - 1)/2;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::update
(
    const volScalarField& psiOld
)
{
    collectStencil(psiOld);
    if (!built_) { build(); }

    // -- restriction: the mean of the children of every parent
    if (maxLevel_ > 0)
    {
        scalarField sum(mesh_.nPoints(), Zero);
        forAll(parentPoint_, c)
        {
            if (parentPoint_[c] >= 0) { sum[parentPoint_[c]] += psiOld[c]; }
        }
        syncTools::syncPointList(mesh_, sum, plusEqOp<scalar>(), scalar(0));
        scalarField& pm = parentMeanFld_.primitiveFieldRef();
        forAll(pm, c)
        {
            const label p = parentPoint_[c];
            pm[c] = p >= 0 ? sum[p]/max(parentCount_[p], scalar(1)) : psiOld[c];
        }
        parentMeanFld_.correctBoundaryConditions();
        stencil_.collectData(parentMeanFld_, stParentMean_);
    }

    // -- the fits, level by level from the coarsest level; the coefficients of
    //    the levels done so far are exchanged before the next level
    nGhost_ = nRestricted_ = nLinear_ = nConstant_ = 0;
    for (label L = 0; L <= maxLevel_; ++L)
    {
        if (L > 0)
        {
            coeffFld_.correctBoundaryConditions();
            stencil_.collectData(coeffFld_, stCoeff_);
        }
        forAll(cellLevel_, c)
        {
            if (cellLevel_[c] == L) { fitCell(c); }
        }
    }

    if (!reported_)
    {
        label g = nGhost_, r = nRestricted_, l = nLinear_, k = nConstant_;
        reduce(g, sumOp<label>());
        reduce(r, sumOp<label>());
        reduce(l, sumOp<label>());
        reduce(k, sumOp<label>());
        Info<< "ghostCellQuadraticWeightedLeastSquares: ghost data " << g
            << ", restricted data " << r << ", linear fits " << l
            << ", constant fits " << k << endl;
        reported_ = true;
    }

    computeLimiters();
}


Foam::scalar
Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    return poly(coeffFld_[c], stencilPsi_[c][0], x - mesh_.C()[c]);
}


// ************************************************************************* //
