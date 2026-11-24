/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

\*---------------------------------------------------------------------------*/

#include "cubicLatticeQuadraticWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "syncTools.H"
#include "cellModel.H"
#include "labelIOList.H"
#include "DynamicList.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        cubicLatticeQuadraticWeightedLeastSquaresReconstruction,
        0
    );
    addToRunTimeSelectionTable
    (
        slReconstruction,
        cubicLatticeQuadraticWeightedLeastSquaresReconstruction,
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

void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::basis
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


Foam::scalar Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::poly
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


Foam::label
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::latticeSlot
(
    const vector& xi
)
{
    label idx = 0;
    for (direction k = 0; k < 3; ++k)
    {
        const scalar r = std::round(xi[k]);
        if (mag(xi[k] - r) > 0.05 || mag(r) > 1.5) { return -1; }
        idx = 3*idx + (label(r) + 1);
    }
    if (idx == 13) { return -1; }                     // the centre
    return idx < 13 ? idx : idx - 1;
}


void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
buildUnitOperator()
{
    // The anchored weighted least-squares fit on the unit lattice, weights
    // 1/distance^2, basis (xi_a, xi_a^2/2, xi_a xi_b) in lattice coordinates.
    // Column j of P_ is the solution for the right-hand side e_j.
    scalar B[26*9], w2[26], rhs[26], x[9];
    for (label i0 = -1; i0 <= 1; ++i0)
    for (label i1 = -1; i1 <= 1; ++i1)
    for (label i2 = -1; i2 <= 1; ++i2)
    {
        const vector xi(i0, i1, i2);
        const label k = latticeSlot(xi);
        if (k < 0) { continue; }
        scalar* b = B + k*9;
        b[0] = xi.x(); b[1] = xi.y(); b[2] = xi.z();
        b[3] = 0.5*xi.x()*xi.x(); b[4] = 0.5*xi.y()*xi.y(); b[5] = 0.5*xi.z()*xi.z();
        b[6] = xi.x()*xi.y(); b[7] = xi.x()*xi.z(); b[8] = xi.y()*xi.z();
        w2[k] = 1.0/magSqr(xi);
    }
    haveP_ = true;
    for (label j = 0; j < 26 && haveP_; ++j)
    {
        for (label k = 0; k < 26; ++k) { rhs[k] = (k == j ? 1.0 : 0.0); }
        if (!solveFit(B, w2, rhs, 26, 9, 9, x)) { haveP_ = false; break; }
        for (label a = 0; a < 9; ++a) { P_[a*26 + j] = x[a]; }
    }
    if (!haveP_)
    {
        FatalErrorInFunction
            << "the unit operator of the 26-point lattice is singular"
            << exit(FatalError);
    }
}


void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::fitUnit
(
    const label c,
    const scalar* rl
)
{
    // -- coefficients in lattice coordinates
    scalar q[9];
    for (label a = 0; a < 9; ++a)
    {
        const scalar* Pa = P_ + a*26;
        scalar v = 0;
        for (label k = 0; k < 26; ++k) { v += Pa[k]*rl[k]; }
        q[a] = v;
    }
    // -- to physical coordinates: xi = d . M with M = inv(lattice), so
    //    g = M . g_xi and H = M . H_xi . M^T
    const tensor M = inv(lattice_[c]);
    const vector gxi(q[0], q[1], q[2]);
    const symmTensor Hxi(q[3], q[6], q[7], q[4], q[8], q[5]);
    const vector g = M & gxi;
    const tensor H = (M & Hxi) & M.T();
    tensor t(Zero);
    t[0] = g.x(); t[1] = g.y(); t[2] = g.z();
    t[3] = H.xx(); t[4] = H.yy(); t[5] = H.zz();
    t[6] = 0.5*(H.xy() + H.yx());
    t[7] = 0.5*(H.xz() + H.zx());
    t[8] = 0.5*(H.yz() + H.zy());
    coeffFld_.primitiveFieldRef()[c] = t;
    order_[c] = 2;
}


void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::build()
{
    const label nCells = mesh_.nCells();

    // -- the refinement levels: the registered objects of a refining dynamic
    //    mesh if they exist, else the files that refineHexMesh writes
    const bool registered =
        mesh_.foundObject<labelIOList>("cellLevel")
     && mesh_.foundObject<labelIOList>("pointLevel");
    labelIOList cellLevelFile
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
    labelIOList pointLevelFile
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
    const labelList& cellLevel =
        registered
      ? static_cast<const labelList&>(mesh_.lookupObject<labelIOList>("cellLevel"))
      : static_cast<const labelList&>(cellLevelFile);
    const labelList& pointLevel =
        registered
      ? static_cast<const labelList&>(mesh_.lookupObject<labelIOList>("pointLevel"))
      : static_cast<const labelList&>(pointLevelFile);
    if (cellLevel.size() != nCells || pointLevel.size() != mesh_.nPoints())
    {
        FatalErrorInFunction
            << "cellLevel or pointLevel ("
            << (registered ? "registered objects" : "files in polyMesh")
            << ") does not have the size of the mesh" << exit(FatalError);
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

    // -- a cell with hanging nodes is not a hexahedron in the cell-shape sense: a
    //    cell next to a finer cell has split faces. Its lattice comes from its
    //    eight anchor points (the points of a level not above the cell's level,
    //    which are the corners of the original hexahedron): the three edges
    //    from one anchor to its three nearest anchors.
    label nAnchorLattice = 0;
    forAll(shapes, c)
    {
        if (shapes[c].model() == hex) { continue; }
        DynamicList<label> anchors(8);
        forAll(cp[c], i)
        {
            if (pointLevel[cp[c][i]] <= cellLevel_[c]) { anchors.append(cp[c][i]); }
        }
        if (anchors.size() != 8) { continue; }
        const point& p0 = pts[anchors[0]];
        label near[3] = {-1, -1, -1};
        scalar dist[3] = {GREAT, GREAT, GREAT};
        for (label i = 1; i < 8; ++i)
        {
            const scalar d = magSqr(pts[anchors[i]] - p0);
            for (label k = 0; k < 3; ++k)
            {
                if (d < dist[k])
                {
                    for (label m = 2; m > k; --m) { dist[m] = dist[m-1]; near[m] = near[m-1]; }
                    dist[k] = d; near[k] = anchors[i];
                    break;
                }
            }
        }
        const vector a0 = pts[near[0]] - p0;
        const vector a1 = pts[near[1]] - p0;
        const vector a2 = pts[near[2]] - p0;
        if (mag(a0 & (a1 ^ a2)) > 1e-6*Foam::pow(dist[2], 1.5))
        {
            lattice_[c] = tensor(a0, a1, a2);
            ++nAnchorLattice;
        }
    }
    reduce(nAnchorLattice, sumOp<label>());

    // -- the static fields of the stencil entries
    forAll(cellLevel_, c)
    {
        levelFld_.primitiveFieldRef()[c] = cellLevel_[c];
        parentCentreFld_.primitiveFieldRef()[c] =
            parentPoint_[c] >= 0 ? pts[parentPoint_[c]] : mesh_.C()[c];
    }
    // An INCOMPLETE parent: some of its children are refined again, so fewer
    // than 2^d cells of this level share the parent centre point. The mean of
    // the present children is the value of a linear field at the CENTROID of
    // their centres, not at the parent centre. The datum of such a parent is
    // placed at that centroid. MEASURED 2026-11-23 (leiaTestMeshTransfer, STATUS
    // 11.44 item 12): with the datum at the parent centre the fit of a linear
    // field had a gradient error of h/64 ... 1.6e-2 in the cells next to such a
    // parent, on a blastAMR mesh with three levels.
    label nIncomplete = 0;
    {
        // Sums relative to the parent centre point: the coordinates of a point
        // are equal on every rank, and the moment is free of cancellation.
        const scalar full = scalar(label(1) << activeDirs_.size());
        vectorField sumD(mesh_.nPoints(), Zero);
        tensorField sumDD(mesh_.nPoints(), Zero);
        forAll(parentPoint_, c)
        {
            const label pp = parentPoint_[c];
            if (pp >= 0)
            {
                const vector dk(mesh_.C()[c] - pts[pp]);
                sumD[pp] += dk;
                sumDD[pp] += dk*dk;
            }
        }
        syncTools::syncPointList(mesh_, sumD, plusEqOp<vector>(), vector::zero);
        syncTools::syncPointList(mesh_, sumDD, plusEqOp<tensor>(), tensor::zero);
        tensorField& Mf = parentMomentFld_.primitiveFieldRef();
        Mf = tensor::zero;
        forAll(parentPoint_, c)
        {
            const label pp = parentPoint_[c];
            if (pp < 0) { continue; }
            const scalar nk = max(parentCount_[pp], scalar(1));
            const vector dbar(sumD[pp]/nk);
            Mf[c] = sumDD[pp]/nk - dbar*dbar;
            if (parentCount_[pp] < full - 0.5)
            {
                parentCentreFld_.primitiveFieldRef()[c] = pts[pp] + dbar;
                ++nIncomplete;
            }
        }
        reduce(nIncomplete, sumOp<label>());
    }
    Info<< "cubicLatticeQuadraticWeightedLeastSquares: cells of an incomplete "
        << "parent (datum at the centroid of the present children): "
        << nIncomplete << endl;
    // A physical patch face has the level -1. This marks a boundary-face entry
    // of the stencil on every rank: the index test stencilIsBoundaryFace() finds
    // only the faces of the own rank, and a patch face of another rank looked
    // like a cell with the level of its cell (STATUS 11.36 item 6).
    forAll(levelFld_.boundaryField(), patchi)
    {
        if (!levelFld_.boundaryField()[patchi].coupled())
        {
            levelFld_.boundaryFieldRef()[patchi] == scalar(-1);
        }
    }
    levelFld_.correctBoundaryConditions();
    parentCentreFld_.correctBoundaryConditions();
    stencil_.collectData(levelFld_, stLevel_);
    stencil_.collectData(parentCentreFld_, stParentC_);
    parentMomentFld_.correctBoundaryConditions();
    stencil_.collectData(parentMomentFld_, stParentM_);

    // -- which cells are cubes, and the lattice positions of the stencil of a
    //    cell that has exactly the 26 cells of its own level around it
    cube_.setSize(nCells);
    slot_.setSize(nCells);
    order_.setSize(nCells);
    order_ = 0;
    label nCube = 0, nSlot = 0;
    const bool threeD = (activeDirs_.size() == 3);
    forAll(cellLevel_, c)
    {
        slot_[c].clear();
        const vector a0 = lattice_[c].x();
        const vector a1 = lattice_[c].y();
        const vector a2 = lattice_[c].z();
        const scalar h2 = magSqr(a0);
        cube_[c] =
            threeD
         && h2 > VSMALL
         && mag(magSqr(a1) - h2) < 1e-6*h2
         && mag(magSqr(a2) - h2) < 1e-6*h2
         && mag(a0 & a1) < 1e-6*h2
         && mag(a0 & a2) < 1e-6*h2
         && mag(a1 & a2) < 1e-6*h2;
        if (!cube_[c]) { continue; }
        ++nCube;
        const label n = stencilSize(c);
        if (n != 27) { continue; }
        const tensor invLat = inv(lattice_[c]);
        const point x0 = stencilC(c, 0);
        labelList sl(n - 1, label(-1));
        boolList used(26, false);
        bool ok = true;
        for (label i = 1; i < n && ok; ++i)
        {
            const label r = raw(c, i);
            if
            (
                stencilIsBoundaryFace(c, i)
             || stLevel_[c][r] < -0.5
             || label(stLevel_[c][r] + 0.5) != cellLevel_[c]
            )
            {
                ok = false;
                break;
            }
            const label k = latticeSlot((stencilC(c, i) - x0) & invLat);
            if (k < 0 || used[k]) { ok = false; break; }
            used[k] = true;
            sl[i - 1] = k;
        }
        if (ok) { slot_[c].transfer(sl); ++nSlot; }
    }
    reduce(nCube, sumOp<label>());
    reduce(nSlot, sumOp<label>());

    label nNoParent = 0, nNoHex = 0;
    forAll(cellLevel_, c)
    {
        if (cellLevel_[c] > 0 && parentPoint_[c] < 0) { ++nNoParent; }
        if (magSqr(lattice_[c]) < VSMALL) { ++nNoHex; }
    }
    reduce(nNoParent, sumOp<label>());
    reduce(nNoHex, sumOp<label>());
    Info<< "cubicLatticeQuadraticWeightedLeastSquares: " << maxLevel_ + 1
        << " refinement level(s); refined cells without a parent point: "
        << nNoParent << "; cells that are not a hexahedron: " << nNoHex
        << " (" << nAnchorLattice << " of the cells with hanging nodes have a"
        << " lattice from their anchor points)"
        << "; cubes: " << nCube << "; cells with the 26 cells of their own"
        << " level around them: " << nSlot << endl;
    built_ = true;
}


void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::fitCell
(
    const label c
)
{
    const label n = stencilSize(c);
    const List<scalar>& s = stencilPsi_[c];
    const scalar psiC = s[0];

    // -- the 26 cells of the own level: the unit operator, no search
    if (haveP_ && slot_[c].size())
    {
        scalar rl[26];
        const labelList& sl = slot_[c];
        for (label i = 1; i < n; ++i) { rl[sl[i - 1]] = s[i] - psiC; }
        fitUnit(c, rl);
        ++nUnitSlot_;
        return;
    }

    const point x0 = stencilC(c, 0);
    const label Lc = cellLevel_[c];
    const scalar tol2 = sqr(1e-3*radius_[c]);

    DynamicList<vector> X(2*n);
    DynamicList<scalar> V(2*n);
    DynamicList<label> coarse(n);
    DynamicList<vector> same(n);
    // restricted data: position in X and the second moment of the children
    DynamicList<label> restrIdx(n);
    DynamicList<tensor> restrM(n);

    // -- data at the own level, and one datum per refined parent
    for (label i = 1; i < n; ++i)
    {
        const label r = raw(c, i);
        const point xi = stencilC(c, i);
        if (stencilIsBoundaryFace(c, i) || stLevel_[c][r] < -0.5)
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
                restrIdx.append(X.size() - 1); restrM.append(stParentM_[c][r]);
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
        const tensor invLat = inv(lattice_[c]);
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
            // -- the coarser stencil cell that CONTAINS q. A hexRef8 cell that
            //    is n levels coarser is the own lattice scaled by 2^n, so q is
            //    in it if the offset from its centre, in own lattice
            //    coordinates, is below 2^n/2 in every direction. The offsets are
            //    odd multiples of 1/2: the test has a margin of half a cell and
            //    does not depend on round-off or on the order of the stencil.
            //    A nearest-centre test has exact ties for a q outside the
            //    domain next to a boundary patch, and the order of the stencil
            //    then decided it (STATUS 11.36 item 6).
            label iA = -1;
            forAll(coarse, j)
            {
                const label ij = coarse[j];
                const label nLev =
                    Lc - label(stLevel_[c][raw(c, ij)] + 0.5);
                const scalar half = 0.5*scalar(label(1) << nLev);
                const vector cc = (q - stencilC(c, ij)) & invLat;
                if (mag(cc.x()) < half && mag(cc.y()) < half && mag(cc.z()) < half)
                {
                    iA = ij;
                    break;
                }
            }
            if (iA < 0) { continue; }             // q is in no coarser stencil cell
            const point xA = stencilC(c, iA);
            X.append(q);
            V.append(poly(stCoeff_[c][raw(c, iA)], s[iA], q - xA));
            ++nGhost_;
        }
    }
    else if (coarse.size())
    {
        // no lattice: keep the coarse cells as they are (counted: the stencil
        // of such a cell mixes two cell sizes)
        ++nMixed_;
        forAll(coarse, j)
        {
            X.append(stencilC(c, coarse[j])); V.append(s[coarse[j]]);
        }
    }

    // -- a cube whose 26 lattice positions hold one datum each: the unit operator
    const label m = X.size();
    // (not with restricted data: their rows carry the second moment, see below)
    if (haveP_ && cube_[c] && m == 26 && restrIdx.empty())
    {
        const tensor invLat = inv(lattice_[c]);
        scalar rl[26];
        bool used[26];
        for (label k = 0; k < 26; ++k) { used[k] = false; }
        bool ok = true;
        forAll(X, i)
        {
            const label k = latticeSlot((X[i] - x0) & invLat);
            if (k < 0 || used[k]) { ok = false; break; }
            used[k] = true;
            rl[k] = V[i] - psiC;
        }
        if (ok)
        {
            fitUnit(c, rl);
            ++nUnitMatched_;
            return;
        }
    }

    // -- every other cell: the anchored weighted least-squares fit per cell,
    //    weights 1/distance
    ++nPerCell_;
    scalarField B(m*ncoeff_), w2(m), rhs(m);
    forAll(X, i)
    {
        const vector d = X[i] - x0;
        basis(d, B.begin() + i*ncoeff_);
        w2[i] = 1.0/max(magSqr(d), VSMALL);
        rhs[i] = V[i] - psiC;
    }
    const label nd = activeDirs_.size();
    // A restricted datum is the MEAN of the children: its row is the mean of the
    // basis over the child centres, basis(centroid) plus the second moment M in
    // the quadratic entries (0.5 M_aa and M_ab). The fit of a quadratic field is
    // exact with it (supervisor, 2026-11-23; leiaTestMeshTransfer).
    forAll(restrIdx, j)
    {
        scalar* row = B.begin() + restrIdx[j]*ncoeff_;
        const tensor& M = restrM[j];
        label k = nd;
        for (label a = 0; a < nd; ++a)
        {
            row[k++] += 0.5*M(activeDirs_[a], activeDirs_[a]);
        }
        for (label a = 0; a < nd; ++a)
        {
            for (label e = a + 1; e < nd; ++e)
            {
                row[k++] += M(activeDirs_[a], activeDirs_[e]);
            }
        }
    }
    scalar x[9];
    tensor t(Zero);
    if (solveFit(B.begin(), w2.begin(), rhs.begin(), m, ncoeff_, ncoeff_, x))
    {
        for (label k = 0; k < ncoeff_; ++k) { t[k] = x[k]; }
        order_[c] = 2;
    }
    else if (solveFit(B.begin(), w2.begin(), rhs.begin(), m, ncoeff_, nd, x))
    {
        for (label k = 0; k < nd; ++k) { t[k] = x[k]; }
        order_[c] = 1;
        ++nLinear_;
    }
    else
    {
        order_[c] = 0;
        ++nConstant_;
    }
    coeffFld_.primitiveFieldRef()[c] = t;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
cubicLatticeQuadraticWeightedLeastSquaresReconstruction
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
            "slCubicLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    ),
    parentCentreFld_
    (
        IOobject
        (
            "slCubicParentCentre",
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
            "slCubicParentMean",
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
    parentMomentFld_
    (
        IOobject
        (
            "slCubicParentMoment",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedTensor(dimArea, Zero),
        zeroGradientFvPatchTensorField::typeName
    ),
    coeffFld_
    (
        IOobject
        (
            "slCubicCoeffs",
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
    haveP_(false),
    cube_(),
    slot_(),
    order_(),
    nUnitSlot_(0),
    nUnitMatched_(0),
    nPerCell_(0),
    nMixed_(0),
    nGhost_(0),
    nRestricted_(0),
    nLinear_(0),
    nConstant_(0),
    reported_(false),
    curvatureNewtonIters_(slDict_.getOrDefault<label>("curvatureNewtonIters", 3)),
    closestPointIters_(slDict_.getOrDefault<label>("closestPointNewtonIters", 10)),
    offsetCorrection_(slDict_.getOrDefault<word>("offsetCorrection", "none")),
    fpTolRel_(slDict_.getOrDefault<scalar>("footPointTolRel", 1e-6)),
    fpAlphaMax_(slDict_.getOrDefault<scalar>("footPointAlphaMax", 20)),
    fpCycles_(slDict_.getOrDefault<label>("footPointCycles", 10)),
    fpSurfIters_(slDict_.getOrDefault<label>("footPointSurfIters", 20)),
    offsetBeta_(slDict_.getOrDefault<scalar>("offsetBeta", 0.25))
{
    if
    (
        fpTolRel_ <= 0 || fpAlphaMax_ <= 0 || fpCycles_ < 1 || fpSurfIters_ < 1
    )
    {
        FatalIOErrorInFunction(slDict_)
            << "footPoint controls must be positive" << exit(FatalIOError);
    }
    if
    (
        offsetCorrection_ != "none"
     && offsetCorrection_ != "psi"
     && offsetCorrection_ != "psiOverGradPsi"
     && offsetCorrection_ != "quadraticRoot"
    )
    {
        FatalIOErrorInFunction(slDict_)
            << "offsetCorrection must be none, psi, psiOverGradPsi or "
            << "quadraticRoot, got '" << offsetCorrection_ << "'"
            << exit(FatalIOError);
    }
    if (offsetBeta_ <= 0 || offsetBeta_ >= 1)
    {
        FatalIOErrorInFunction(slDict_)
            << "offsetBeta must lie in (0,1), got " << offsetBeta_
            << exit(FatalIOError);
    }

    const Vector<label> gd = mesh_.geometricD();
    DynamicList<label> dirs(3);
    for (direction k = 0; k < 3; ++k)
    {
        if (gd[k] > 0) { dirs.append(k); }
    }
    activeDirs_ = dirs;
    const label nd = activeDirs_.size();
    ncoeff_ = 2*nd + nd*(nd - 1)/2;
    if (nd == 3) { buildUnitOperator(); }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::update
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
    nUnitSlot_ = nUnitMatched_ = nPerCell_ = 0;
    nMixed_ = 0;
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
        Info<< "cubicLatticeQuadraticWeightedLeastSquares: ghost data " << g
            << ", restricted data " << r << ", linear fits " << l
            << ", constant fits " << k << endl;
        label us = nUnitSlot_, um = nUnitMatched_, pc = nPerCell_;
        reduce(us, sumOp<label>());
        reduce(um, sumOp<label>());
        reduce(pc, sumOp<label>());
        Info<< "cubicLatticeQuadraticWeightedLeastSquares: unit operator in "
            << us << " cells of one level and " << um
            << " cells with restricted or ghost data; fit per cell in " << pc
            << " cells; cells that keep coarser cells in their stencil: "
            << returnReduce(nMixed_, sumOp<label>()) << endl;
        reported_ = true;
    }

    computeLimiters();
}


Foam::scalar
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    return poly(coeffFld_[c], stencilPsi_[c][0], x - mesh_.C()[c]);
}


Foam::label
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::fitDerivatives
(
    const label c,
    vector& g,
    symmTensor& H
) const
{
    g = Zero;
    H = Zero;
    const label order = order_.size() ? order_[c] : 0;
    if (order == 0) { return 0; }
    const tensor& t = coeffFld_[c];
    const label nd = activeDirs_.size();
    for (label a = 0; a < nd; ++a) { g[activeDirs_[a]] = t[a]; }
    if (order == 1) { return 1; }
    // symmTensor components: xx xy xz yy yz zz
    static const direction diag[3] = {symmTensor::XX, symmTensor::YY, symmTensor::ZZ};
    static const direction off[3][3] =
    {
        {symmTensor::XX, symmTensor::XY, symmTensor::XZ},
        {symmTensor::XY, symmTensor::YY, symmTensor::YZ},
        {symmTensor::XZ, symmTensor::YZ, symmTensor::ZZ}
    };
    label k = nd;
    for (label a = 0; a < nd; ++a) { H[diag[activeDirs_[a]]] = t[k++]; }
    for (label a = 0; a < nd; ++a)
    {
        for (label b = a + 1; b < nd; ++b)
        {
            H[off[activeDirs_[a]][activeDirs_[b]]] = t[k++];
        }
    }
    return 2;
}

// * * * * * * * * * * * * * * Curvature functions * * * * * * * * * * * * * //
//
// Every function reads the polynomial of ONE cell:
//     P(d) = psi_c + g.d + 0.5 d.H.d,    grad P(d) = g + H.d,
// with g and H from fitDerivatives. A linear fit has H = 0.

Foam::scalar
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::curvatureOffset
(
    const scalar psiC,
    const vector& g,
    const scalar gm,
    const symmTensor& H
) const
{
    if (offsetCorrection_ == "psi")
    {
        return psiC;
    }
    else if (offsetCorrection_ == "psiOverGradPsi")
    {
        return psiC/gm;
    }
    else if (offsetCorrection_ == "quadraticRoot")
    {
        const vector n = g/gm;
        bool fb = false;
        return offsetDistance(psiC, gm, n & (H & n), fb);
    }
    return 0;
}


bool Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::footPoint
(
    const label c,
    const vector& dp,
    const scalar level,
    vector& dFoot,
    vector& gFoot,
    symmTensor& H
) const
{
    vector g(Zero);
    if (fitDerivatives(c, g, H) == 0) { return false; }

    const scalar psiC = stencilPsi_[c][0];
    const scalar radius = stencilRadius(c);
    if (radius < SMALL) { return false; }
    const scalar eps = fpTolRel_*radius;

    auto G = [&](const vector& d) -> scalar
    {
        return psiC - level + (g & d) + 0.5*(d & (H & d));
    };
    auto gradG = [&](const vector& d) -> vector
    {
        return g + (H & d);
    };

    // Newton along the local gradient onto {G = 0}.
    auto surfacepoint = [&](vector& d) -> bool
    {
        for (label k = 0; k < fpSurfIters_; ++k)
        {
            const vector gr(gradG(d));
            const scalar g2 = gr & gr;
            if (g2 < SMALL) { return false; }
            const vector step = -(G(d)/g2)*gr;
            d += step;
            if (Foam::mag(step) < eps) { return true; }
        }
        return false;
    };

    vector dpi = dp;
    if (!surfacepoint(dpi)) { return false; }

    bool converged = false;
    for (label cycle = 0; cycle < fpCycles_ && !converged; ++cycle)
    {
        // (1) foot point of the query on the tangent plane at dpi
        const vector gi(gradG(dpi));
        const scalar gi2 = gi & gi;
        if (gi2 < SMALL) { return false; }
        const scalar lam = ((dp - dpi) & gi)/gi2;
        vector dqi = dp - lam*gi;

        // (2) back to the surface
        vector dpn = dqi;
        if (!surfacepoint(dpn)) { return false; }

        // (3) parabola correction from the two step vectors
        const vector f1 = dqi - dpi;
        const vector f2 = dpn - dqi;
        if (Foam::mag(f1) > eps)
        {
            const vector rel = dp - dpi;
            const scalar a0 = rel & f1;
            const scalar a1 = 2.0*(f2 & rel) - magSqr(f1);
            const scalar a2 = -3.0*(f1 & f2);
            const scalar a3 = -2.0*magSqr(f2);
            const scalar den = a1 + 2.0*a2 + 3.0*a3;
            if (Foam::mag(den) > SMALL)
            {
                const scalar alpha = 1.0 - (a0 + a1 + a2 + a3)/den;
                if (alpha > 0 && alpha < fpAlphaMax_)
                {
                    dqi = dpi + alpha*f1 + sqr(alpha)*f2;
                    dpn = dqi;
                    if (!surfacepoint(dpn)) { return false; }
                }
            }
        }

        converged = (Foam::mag(dpn - dpi) < eps);
        dpi = dpn;
    }

    // Trust region: the fit holds only where the stencil has data.
    if (!converged || Foam::mag(dpi) > radius) { return false; }

    gFoot = gradG(dpi);
    if (Foam::mag(gFoot) < SMALL) { return false; }
    dFoot = dpi;
    return true;
}


void Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::meanCurvature
(
    volScalarField& kappa
) const
{
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        vector g(Zero);
        symmTensor H(Zero);
        if (fitDerivatives(c, g, H) < 2) { continue; }

        const scalar psiC = stencilPsi_[c][0];
        // Band gate with the first-order distance |psi|/|grad psi|.
        if (Foam::mag(psiC)/max(Foam::mag(g), SMALL) > 3.0*stencilRadius(c))
        {
            continue;
        }

        // Newton along the gradient from the cell centre to the zero set.
        vector delta(Zero);
        for (label it = 0; it < curvatureNewtonIters_; ++it)
        {
            const scalar p = psiC + (g & delta) + 0.5*(delta & (H & delta));
            const vector gr(g + (H & delta));
            const scalar g2 = gr & gr;
            if (g2 < SMALL) { break; }
            delta -= (p/g2)*gr;
        }

        const vector gf(g + (H & delta));
        const scalar gm = Foam::mag(gf);
        if (gm < SMALL) { continue; }

        const scalar kd = kappaOf(gf, H);
        k[c] =
            (offsetCorrection_ == "none")
          ? kd
          : offsetCorrected(kd, curvatureOffset(psiC, gf, gm, H));
    }

    kappa.correctBoundaryConditions();
}


void
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
meanCurvatureLaplacian
(
    volScalarField& kappa
) const
{
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        vector g(Zero);
        symmTensor H(Zero);
        if (fitDerivatives(c, g, H) < 2) { continue; }

        const scalar psiC = stencilPsi_[c][0];
        const scalar gm = max(Foam::mag(g), SMALL);
        if (Foam::mag(psiC)/gm > 3.0*stencilRadius(c)) { continue; }

        k[c] =
            (offsetCorrection_ == "none")
          ? tr(H)
          : offsetCorrected(tr(H), curvatureOffset(psiC, g, gm, H));
    }

    kappa.correctBoundaryConditions();
}


void
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
meanCurvatureNoExtension
(
    volScalarField& kappa
) const
{
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        vector g(Zero);
        symmTensor H(Zero);
        if (fitDerivatives(c, g, H) < 2) { continue; }

        const scalar psiC = stencilPsi_[c][0];
        const scalar gm = Foam::mag(g);
        if (gm < SMALL) { continue; }
        if (Foam::mag(psiC)/gm > 3.0*stencilRadius(c)) { continue; }

        const scalar kd = kappaOf(g, H);
        k[c] =
            (offsetCorrection_ == "none")
          ? kd
          : offsetCorrected(kd, curvatureOffset(psiC, g, gm, H));
    }

    kappa.correctBoundaryConditions();
}


void
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
meanCurvatureClosestPoint
(
    volScalarField& kappa
) const
{
    // The foot is the closest point: min |delta|^2 with P(delta) = 0. Newton
    // on delta = lambda grad P(delta), P(delta) = 0, with the Jacobian
    //     [ I - lambda H   -grad P ]
    //     [ grad P^T          0    ]
    // Every guard failure keeps the cell-centre value.
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        vector g0(Zero);
        symmTensor H(Zero);
        if (fitDerivatives(c, g0, H) < 2) { continue; }
        const tensor Ht(H);

        const scalar psiC = stencilPsi_[c][0];
        const scalar gm0 = Foam::mag(g0);
        if (gm0 < SMALL) { continue; }
        if (Foam::mag(psiC)/gm0 > 3.0*stencilRadius(c)) { continue; }

        k[c] = kappaOf(g0, H);

        const scalar rTrust = stencilRadius(c);
        const scalar tolP = 1e-8*rTrust;

        scalar lambda = -psiC/(gm0*gm0);
        vector delta = lambda*g0;

        bool ok = false;
        bool failed = false;
        for (label it = 0; it < closestPointIters_; ++it)
        {
            const scalar p = psiC + (g0 & delta) + 0.5*(delta & (H & delta));
            const vector gr(g0 + (H & delta));
            if (magSqr(gr) < SMALL) { failed = true; break; }

            const vector rV = delta - lambda*gr;
            if (Foam::mag(p) < tolP && Foam::mag(rV) < 1e-8*rTrust)
            {
                ok = true;
                break;
            }

            scalar A[4][5];
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    A[i][j] = (i == j ? 1.0 : 0.0) - lambda*Ht(i, j);
                }
                A[i][3] = -gr[i];
                A[i][4] = -rV[i];
                A[3][i] = gr[i];
            }
            A[3][3] = 0;
            A[3][4] = -p;

            for (label col = 0; col < 4 && !failed; ++col)
            {
                label piv = col;
                for (label r = col + 1; r < 4; ++r)
                {
                    if (Foam::mag(A[r][col]) > Foam::mag(A[piv][col]))
                    {
                        piv = r;
                    }
                }
                if (Foam::mag(A[piv][col]) < VSMALL) { failed = true; break; }
                if (piv != col)
                {
                    for (label j = col; j < 5; ++j)
                    {
                        const scalar t = A[col][j];
                        A[col][j] = A[piv][j];
                        A[piv][j] = t;
                    }
                }
                for (label r = col + 1; r < 4; ++r)
                {
                    const scalar f = A[r][col]/A[col][col];
                    for (label j = col; j < 5; ++j) { A[r][j] -= f*A[col][j]; }
                }
            }
            if (failed) { break; }

            scalar sol[4];
            for (label i = 3; i >= 0; --i)
            {
                scalar s = A[i][4];
                for (label j = i + 1; j < 4; ++j) { s -= A[i][j]*sol[j]; }
                sol[i] = s/A[i][i];
            }
            if
            (
                !std::isfinite(sol[0]) || !std::isfinite(sol[1])
             || !std::isfinite(sol[2]) || !std::isfinite(sol[3])
            )
            {
                failed = true;
                break;
            }

            delta += vector(sol[0], sol[1], sol[2]);
            lambda += sol[3];
            if (Foam::mag(delta) > rTrust) { failed = true; break; }
        }

        if (failed || !ok) { continue; }

        const vector gf(g0 + (H & delta));
        if (Foam::mag(gf) < SMALL) { continue; }
        const scalar kFoot = kappaOf(gf, H);
        if (std::isfinite(kFoot)) { k[c] = kFoot; }
    }

    kappa.correctBoundaryConditions();
}


Foam::scalar
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::signedOffset
(
    const label c,
    bool& fallback
) const
{
    vector g(Zero);
    symmTensor H(Zero);
    const label order = fitDerivatives(c, g, H);
    const scalar gm = Foam::mag(g);

    if (order == 0 || gm < SMALL)
    {
        fallback = true;
        return 0;
    }

    const vector n = g/gm;
    return offsetDistance(stencilPsi_[c][0], gm, n & (H & n), fallback);
}


Foam::scalar
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::footPointDistance
(
    const label c,
    const point& p,
    const scalar level,
    bool& ok
) const
{
    const vector dp(p - mesh_.C()[c]);
    vector dFoot(Zero), gFoot(Zero);
    symmTensor H(Zero);
    ok = footPoint(c, dp, level, dFoot, gFoot, H);
    if (!ok) { return 0; }

    // Signed distance; positive on the side of larger psi.
    return ((dp - dFoot) & gFoot)/Foam::mag(gFoot);
}


Foam::scalar
Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction::
curvatureAtFootPoint
(
    const label c,
    const point& p,
    bool& ok
) const
{
    const vector dp(p - mesh_.C()[c]);
    vector dFoot(Zero), gFoot(Zero);
    symmTensor H(Zero);
    ok = footPoint(c, dp, 0.0, dFoot, gFoot, H);
    if (!ok) { return 0; }

    return kappaOf(gFoot, H);
}


// ************************************************************************* //
