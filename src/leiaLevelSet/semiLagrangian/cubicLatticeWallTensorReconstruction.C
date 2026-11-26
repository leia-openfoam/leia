/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "cubicLatticeWallTensorReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "SquareMatrix.H"
#include "emptyPolyPatch.H"
#include "boundBox.H"
#include "calculatedFvPatchFields.H"

namespace Foam
{
    defineTypeNameAndDebug(cubicLatticeWallTensorReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        cubicLatticeWallTensorReconstruction,
        Mesh
    );
}

namespace
{
    // Lagrange weights of the three nodes n at x.
    inline void lagrange3(const Foam::scalar n[3], const Foam::scalar x, Foam::scalar w[3])
    {
        w[0] = (x - n[1])*(x - n[2])/((n[0] - n[1])*(n[0] - n[2]));
        w[1] = (x - n[0])*(x - n[2])/((n[1] - n[0])*(n[1] - n[2]));
        w[2] = (x - n[0])*(x - n[1])/((n[2] - n[0])*(n[2] - n[1]));
    }

    // Least-squares quadratic in two variables through the points (s, t, v),
    // evaluated at (s0, t0). Exact for a quadratic when the points determine it.
    Foam::scalar planeQuadratic
    (
        const Foam::List<Foam::FixedList<Foam::scalar, 3>>& p,
        const Foam::scalar s0,
        const Foam::scalar t0
    )
    {
        using namespace Foam;
        SquareMatrix<scalar> M(6, Zero);
        FixedList<scalar, 6> r(Zero);
        for (const FixedList<scalar, 3>& q : p)
        {
            const scalar b[6] = {1, q[0], q[1], q[0]*q[0], q[1]*q[1], q[0]*q[1]};
            for (label i = 0; i < 6; ++i)
            {
                r[i] += b[i]*q[2];
                for (label j = 0; j < 6; ++j) { M(i, j) += b[i]*b[j]; }
            }
        }
        // Gauss elimination with partial pivoting (6 x 6).
        for (label k = 0; k < 6; ++k)
        {
            label piv = k;
            for (label i = k + 1; i < 6; ++i)
            {
                if (mag(M(i, k)) > mag(M(piv, k))) { piv = i; }
            }
            if (piv != k)
            {
                for (label j = 0; j < 6; ++j) { std::swap(M(k, j), M(piv, j)); }
                std::swap(r[k], r[piv]);
            }
            for (label i = k + 1; i < 6; ++i)
            {
                const scalar f = M(i, k)/M(k, k);
                for (label j = k; j < 6; ++j) { M(i, j) -= f*M(k, j); }
                r[i] -= f*r[k];
            }
        }
        FixedList<scalar, 6> a(Zero);
        for (label k = 5; k >= 0; --k)
        {
            scalar s = r[k];
            for (label j = k + 1; j < 6; ++j) { s -= M(k, j)*a[j]; }
            a[k] = s/M(k, k);
        }
        return a[0] + a[1]*s0 + a[2]*t0 + a[3]*s0*s0 + a[4]*t0*t0 + a[5]*s0*t0;
    }

    // Tensor index of the nodes (i, j, k) of x, y, z.
    inline Foam::label tIdx(const Foam::label i, const Foam::label j, const Foam::label k)
    {
        return 9*i + 3*j + k;
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::cubicLatticeWallTensorReconstruction::kindOf
(
    const FixedList<label, 3>& w
)
{
    const label n = w[0] + w[1] + w[2];
    if (n == 0) { return 0; }
    if (n == 1) { return 1 + (w[0] ? 0 : (w[1] ? 1 : 2)); }
    if (n == 2) { return (w[0] && w[1]) ? 4 : ((w[0] && w[2]) ? 5 : 6); }
    return 7;
}


void Foam::cubicLatticeWallTensorReconstruction::buildWall()
{
    const label nCells = mesh_.nCells();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& C = mesh_.C().primitiveField();

    side_.setSize(nCells);
    wallFace_.setSize(nCells);
    for (label c = 0; c < nCells; ++c)
    {
        side_[c] = 0;
        wallFace_[c] = -1;
    }
    boolList bad(nCells, false);

    // -- the walls of every cell: physical patch faces normal to an axis
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() || isA<emptyPolyPatch>(pp)) { continue; }
        forAll(pp, i)
        {
            const label f = pp.start() + i;
            const label c = mesh_.faceOwner()[f];
            const vector n = mesh_.faceAreas()[f]/mag(mesh_.faceAreas()[f]);
            label a = 0;
            for (label k = 1; k < 3; ++k)
            {
                if (mag(n[k]) > mag(n[a])) { a = k; }
            }
            if (mag(n[a]) < 1 - 1e-9) { bad[c] = true; continue; }
            const label s = n[a] > 0 ? 1 : -1;
            if (side_[c][a] != 0 && side_[c][a] != s) { bad[c] = true; continue; }
            if (wallFace_[c][a] >= 0) { bad[c] = true; continue; }  // split wall face
            side_[c][a] = s;
            wallFace_[c][a] = f;
        }
    }

    // -- the cell size; a candidate is an axis-aligned cube with 6 faces and 8 points
    h_.setSize(nCells);
    h_ = 0;
    const labelListList& cp = mesh_.cellPoints();
    const pointField& pts = mesh_.points();
    boolList cand(nCells, false);
    // h of EVERY axis-aligned cube (the owners of the interior tensor points are
    // not wall cells; their size is compared below)
    forAll(cand, c)
    {
        if (mesh_.cells()[c].size() != 6 || cp[c].size() != 8) { continue; }
        boundBox bb(pts, cp[c], false);
        const vector sp = bb.span();
        if
        (
            mag(sp.y() - sp.x()) > 1e-6*sp.x()
         || mag(sp.z() - sp.x()) > 1e-6*sp.x()
        )
        {
            continue;
        }
        h_[c] = sp.x();
        if (bad[c]) { continue; }
        if (side_[c][0] == 0 && side_[c][1] == 0 && side_[c][2] == 0) { continue; }
        cand[c] = true;
    }

    // -- the raw stencil centres and the static data of the stencil entries
    {
        volVectorField Cf("wallTensorC", mesh_.C());
        stencil_.collectData(Cf, stC_);
    }
    volScalarField sideFld
    (
        IOobject("wallTensorSide", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh_, dimensionedScalar(dimless, -1), calculatedFvPatchScalarField::typeName
    );
    volScalarField hFld
    (
        IOobject("wallTensorH", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh_, dimensionedScalar(dimLength, 0), calculatedFvPatchScalarField::typeName
    );
    volScalarField candFld
    (
        IOobject("wallTensorCand", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh_, dimensionedScalar(dimless, 0), calculatedFvPatchScalarField::typeName
    );
    forAll(side_, c)
    {
        sideFld.primitiveFieldRef()[c] =
            (side_[c][0] + 1) + 3*(side_[c][1] + 1) + 9*(side_[c][2] + 1);
        hFld.primitiveFieldRef()[c] = h_[c];
    }
    sideFld.correctBoundaryConditions();
    hFld.correctBoundaryConditions();
    List<List<scalar>> stSide, stH;
    stencil_.collectData(sideFld, stSide);
    stencil_.collectData(hFld, stH);

    // -- the 27 points of every candidate and their owners
    tRaw_.setSize(nCells);
    tKind_.setSize(nCells);
    tensor_.setSize(nCells);
    tensor_ = false;
    boolList geomOk(nCells, false);
    label nCand = 0, nMissing = 0, nSize = 0, nSide = 0;
    forAll(cand, c) { if (cand[c]) { ++nCand; } }
    forAll(cand, c)
    {
        if (!cand[c]) { continue; }
        const scalar h = h_[c];
        bool ok = true;
        for (label i = 0; i < 3 && ok; ++i)
        for (label j = 0; j < 3 && ok; ++j)
        for (label k = 0; k < 3 && ok; ++k)
        {
            const label id[3] = {i, j, k};
            vector off(Zero);
            FixedList<label, 3> atWall(0);
            for (label a = 0; a < 3; ++a)
            {
                scalar n[3];
                nodes(c, a, n);
                const scalar v = n[id[a]];
                if (mag(mag(v) - 0.5) < 1e-12) { atWall[a] = 1; off[a] = 0; }
                else { off[a] = v; }
            }
            const point target = C[c] + h*off;
            label r = -1;
            forAll(stC_[c], m)
            {
                if (magSqr(stC_[c][m] - target) < 1e-12*h*h) { r = m; break; }
            }
            if (r < 0) { ok = false; ++nMissing; break; }
            // a cell entry of the same size (a face entry has the level-free
            // code -1 in sideFld only on the own rank; the size test excludes
            // faces everywhere: a face carries h = 0)
            if (mag(stH[c][r] - h) > 1e-6*h) { ok = false; ++nSize; break; }
            // the owner of a boundary point lies on the same walls
            if (atWall[0] + atWall[1] + atWall[2] > 0)
            {
                const label code = label(stSide[c][r] + 0.5);
                const label os[3] = {code % 3 - 1, (code/3) % 3 - 1, code/9 - 1};
                for (label a = 0; a < 3; ++a)
                {
                    if (atWall[a] && os[a] != side_[c][a]) { ok = false; }
                }
                if (!ok) { ++nSide; break; }
            }
            tRaw_[c][tIdx(i, j, k)] = r;
            tKind_[c][tIdx(i, j, k)] = kindOf(atWall);
        }
        geomOk[c] = ok;
    }

    // -- an owner of a boundary point must itself be a tensor cell (it advances
    //    the value); prune until no candidate reads a non-tensor owner
    tensor_ = geomOk;
    for (label iter = 0; iter < 20; ++iter)
    {
        forAll(tensor_, c) { candFld.primitiveFieldRef()[c] = tensor_[c] ? 1 : 0; }
        candFld.correctBoundaryConditions();
        List<List<scalar>> stCand;
        stencil_.collectData(candFld, stCand);
        label nPruned = 0;
        forAll(tensor_, c)
        {
            if (!tensor_[c]) { continue; }
            for (label t = 0; t < 27; ++t)
            {
                if (tKind_[c][t] > 0 && stCand[c][tRaw_[c][t]] < 0.5)
                {
                    tensor_[c] = false; ++nPruned; break;
                }
            }
        }
        if (returnReduce(nPruned, sumOp<label>()) == 0) { break; }
    }

    // -- the state fields (index 1..7)
    state_.setSize(8);
    for (label k = 1; k < 8; ++k)
    {
        state_.set
        (
            k,
            new volScalarField
            (
                IOobject("wallTensorState" + Foam::name(k), mesh_.time().timeName(),
                         mesh_, IOobject::NO_READ, IOobject::NO_WRITE,
                         IOobject::NO_REGISTER),
                mesh_, dimensionedScalar(dimLength, 0),
                calculatedFvPatchScalarField::typeName
            )
        );
    }
    stState_.setSize(8);
    tVal_.setSize(nCells);

    label nWall = 0, nTensor = 0;
    forAll(side_, c)
    {
        if (side_[c][0] || side_[c][1] || side_[c][2] || bad[c]) { ++nWall; }
        if (tensor_[c]) { ++nTensor; }
    }
    reduce(nWall, sumOp<label>());
    reduce(nTensor, sumOp<label>());
    Info<< "cubicLatticeWallTensor: candidates (axis-aligned cubes at a wall) "
        << returnReduce(nCand, sumOp<label>()) << "; rejected: tensor point not in the"
        << " stencil " << returnReduce(nMissing, sumOp<label>()) << ", other cell size "
        << returnReduce(nSize, sumOp<label>()) << ", owner on another wall "
        << returnReduce(nSide, sumOp<label>()) << endl;
    Info<< "cubicLatticeWallTensor: " << nWall << " wall cells, " << nTensor
        << " with the tensor interpolant and advected boundary points, "
        << nWall - nTensor << " with the fit of the base class" << endl;
    wallBuilt_ = true;
}


Foam::scalar Foam::cubicLatticeWallTensorReconstruction::tensorValue
(
    const label c,
    const vector& xi
) const
{
    scalar w[3][3];
    for (label a = 0; a < 3; ++a)
    {
        scalar n[3];
        nodes(c, a, n);
        lagrange3(n, xi[a], w[a]);
    }
    scalar v = 0;
    for (label i = 0; i < 3; ++i)
    for (label j = 0; j < 3; ++j)
    for (label k = 0; k < 3; ++k)
    {
        v += w[0][i]*w[1][j]*w[2][k]*tVal_[c][tIdx(i, j, k)];
    }
    return v;
}


void Foam::cubicLatticeWallTensorReconstruction::initialiseState
(
    const volScalarField& psi
)
{
    // Faces: the patch values of psi.
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(tensor_, c)
    {
        if (!tensor_[c]) { continue; }
        for (label a = 0; a < 3; ++a)
        {
            const label f = wallFace_[c][a];
            if (side_[c][a] == 0 || f < 0) { continue; }
            const label pi = pbm.whichPatch(f);
            state_[1 + a].primitiveFieldRef()[c] =
                psi.boundaryField()[pi][f - pbm[pi].start()];
        }
    }
    // Edges and corners: the plane quadratics through the boundary data.
    List<FixedList<scalar, 3>> dummy;
    boolList all(mesh_.nCells(), true);
    // kinds 4..6, then 7 (inflowPointValue reads the collected state)
    for (label pass = 0; pass < 2; ++pass)
    {
        for (label k = 1; k < 8; ++k)
        {
            state_[k].correctBoundaryConditions();
            stencil_.collectData(state_[k], stState_[k]);
        }
        stencil_.collectData(psi, stPsi_);
        forAll(tensor_, c)
        {
            if (!tensor_[c]) { continue; }
            const label nw =
                (side_[c][0] != 0) + (side_[c][1] != 0) + (side_[c][2] != 0);
            if (pass == 0 && nw >= 2)
            {
                for (label k = 4; k <= 6; ++k)
                {
                    const label a = (k == 6) ? 1 : 0;
                    const label b = (k == 4) ? 1 : 2;
                    if (side_[c][a] && side_[c][b])
                    {
                        state_[k].primitiveFieldRef()[c] =
                            inflowPointValue(c, k, psi, dummy);
                    }
                }
            }
            if (pass == 1 && nw == 3)
            {
                state_[7].primitiveFieldRef()[c] = inflowPointValue(c, 7, psi, dummy);
            }
        }
    }
    stateSet_ = true;
}


Foam::scalar Foam::cubicLatticeWallTensorReconstruction::inflowPointValue
(
    const label c,
    const label k,
    const volScalarField& psi,
    const List<FixedList<scalar, 3>>&
) const
{
    // The data are read from the collected stencil lists stPsi_ (cell values)
    // and stState_ (boundary values); the caller collects them first.
    // Axes of the plane: an edge (a, b) uses the plane of a and b at the node 1
    // of the third axis; a corner uses the wall plane of axis 0 (axes 1 and 2).
    label a, b, e;
    if (k == 4) { a = 0; b = 1; e = 2; }
    else if (k == 5) { a = 0; b = 2; e = 1; }
    else if (k == 6) { a = 1; b = 2; e = 0; }
    else { a = 1; b = 2; e = 0; }    // corner: the plane of the wall of axis 0
    scalar na[3], nb[3];
    nodes(c, a, na);
    nodes(c, b, nb);

    // value of the tensor point with the node indices (ia, ib) in the plane and
    // node ie on the third axis, read from the collected data
    auto valueAt = [&](const label ia, const label ib, const label ie) -> scalar
    {
        label id[3];
        id[a] = ia; id[b] = ib; id[e] = ie;
        const label t = tIdx(id[0], id[1], id[2]);
        const label kd = tKind_[c][t];
        const label r = tRaw_[c][t];
        return kd == 0 ? stPsi_[c][r] : stState_[kd][c][r];
    };

    List<FixedList<scalar, 3>> pts;
    // the third axis node: 1 (offset 0) for an edge; for a corner node 0 (the wall)
    const label ie = (k == 7) ? 0 : 1;
    for (label ia = 0; ia < 3; ++ia)
    for (label ib = 0; ib < 3; ++ib)
    {
        if (ia == 0 && ib == 0) { continue; }   // the point itself
        FixedList<scalar, 3> q;
        q[0] = na[ia]; q[1] = nb[ib]; q[2] = valueAt(ia, ib, ie);
        pts.append(q);
    }
    return planeQuadratic(pts, na[0], nb[0]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cubicLatticeWallTensorReconstruction::cubicLatticeWallTensorReconstruction
(
    const fvMesh& mesh
)
:
    cubicLatticeQuadraticWeightedLeastSquaresReconstruction(mesh),
    side_(), tensor_(), h_(), wallFace_(), tRaw_(), tKind_(), tVal_(),
    state_(), stState_(), stPsi_(), stC_(),
    wallBuilt_(false), stateSet_(false), wallReported_(false),
    nOutside_(0), nInflowPoints_(0), nAdvectedPoints_(0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cubicLatticeWallTensorReconstruction::update
(
    const volScalarField& psiOld
)
{
    cubicLatticeQuadraticWeightedLeastSquaresReconstruction::update(psiOld);
    if (!wallBuilt_) { buildWall(); }
    if (!stateSet_) { initialiseState(psiOld); }

    for (label k = 1; k < 8; ++k)
    {
        state_[k].correctBoundaryConditions();
        stencil_.collectData(state_[k], stState_[k]);
    }
    stencil_.collectData(psiOld, stPsi_);
    forAll(tensor_, c)
    {
        if (!tensor_[c]) { continue; }
        for (label t = 0; t < 27; ++t)
        {
            const label kd = tKind_[c][t];
            const label r = tRaw_[c][t];
            tVal_[c][t] = kd == 0 ? stPsi_[c][r] : stState_[kd][c][r];
        }
    }
    nOutside_ = 0;
}


Foam::scalar Foam::cubicLatticeWallTensorReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    if (!wallBuilt_ || !tensor_[c])
    {
        return cubicLatticeQuadraticWeightedLeastSquaresReconstruction::evaluateRaw(c, x);
    }
    const vector xi((x - mesh_.C()[c])/h_[c]);
    return tensorValue(c, xi);
}


void Foam::cubicLatticeWallTensorReconstruction::advanceBoundary
(
    const volVectorField& uNew,
    const volVectorField& uOld,
    const scalar dt,
    volScalarField& psi
)
{
    if (!wallBuilt_) { return; }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    // The velocity gradient of a wall cell from the CELL values of its tensor
    // stencil only: a one-sided difference normal to a wall (node 0 and the node
    // one cell inward), a central difference along it. Exact for a linear
    // velocity. The Gauss gradient is not used: the trace velocity
    // reconstruct(phi) has the adjacent cell value on a wall face, which biases
    // it by O(|grad u|) in every wall cell (MEASURED 2026-11-26: a boundary-point
    // error of 0.2 h after one rotation, independent of dt).
    List<List<vector>> stU;
    {
        volVectorField uc("wallTensorU", uNew);
        uc.correctBoundaryConditions();
        stencil_.collectData(uc, stU);
    }
    tensorField G(mesh_.nCells(), Zero);
    forAll(tensor_, c)
    {
        if (!tensor_[c]) { continue; }
        vector row[3];
        for (label a = 0; a < 3; ++a)
        {
            label lo[3] = {1, 1, 1}, hi[3] = {1, 1, 1};
            scalar span;
            if (side_[c][a] == 0) { lo[a] = 0; hi[a] = 2; span = 2; }
            else { lo[a] = 1; hi[a] = 2; span = -side_[c][a]; }   // node 2 is inward
            const vector uHi = stU[c][tRaw_[c][tIdx(hi[0], hi[1], hi[2])]];
            const vector uLo = stU[c][tRaw_[c][tIdx(lo[0], lo[1], lo[2])]];
            row[a] = (uHi - uLo)/(span*h_[c]);
        }
        G[c] = tensor(row[0], row[1], row[2]);
    }

    // departure point of a boundary point at displacement xp (cell sizes) of
    // cell c, in cell sizes; inside == false if it leaves through a wall of c
    auto departure = [&](const label c, const vector& xp, bool& inside) -> vector
    {
        const scalar h = h_[c];
        const vector dx(h*xp);
        const vector up(uNew[c] + (dx & G[c]));
        const vector acc((uNew[c] - uOld[c])/dt + (up & G[c]));
        const vector xd(dx - up*dt + 0.5*acc*dt*dt);
        const vector xi(xd/h);
        inside = true;
        for (label a = 0; a < 3; ++a)
        {
            if (side_[c][a] != 0 && xi[a]*side_[c][a] > 0.5 + 1e-9) { inside = false; }
        }
        return xi;
    };

    nInflowPoints_ = nAdvectedPoints_ = 0;
    // new values per cell and kind; NaN marks "from boundary data, later"
    const scalar unset = -GREAT;
    List<FixedList<scalar, 8>> nv(mesh_.nCells());
    forAll(tensor_, c)
    {
        nv[c] = unset;
        if (!tensor_[c]) { continue; }
        for (label k = 1; k < 8; ++k)
        {
            FixedList<label, 3> w(0);
            if (k <= 3) { w[k - 1] = 1; }
            else if (k == 4) { w[0] = w[1] = 1; }
            else if (k == 5) { w[0] = w[2] = 1; }
            else if (k == 6) { w[1] = w[2] = 1; }
            else { w[0] = w[1] = w[2] = 1; }
            bool owned = true;
            vector xp(Zero);
            for (label a = 0; a < 3; ++a)
            {
                if (w[a]) { if (side_[c][a] == 0) { owned = false; } xp[a] = 0.5*side_[c][a]; }
            }
            if (!owned) { continue; }
            bool inside;
            const vector xi = departure(c, xp, inside);
            if (inside)
            {
                nv[c][k] = tensorValue(c, xi);
                ++nAdvectedPoints_;
            }
            else if (k <= 3)
            {
                const label f = wallFace_[c][k - 1];
                const label pi = pbm.whichPatch(f);
                nv[c][k] = psi.boundaryField()[pi][f - pbm[pi].start()];
                ++nInflowPoints_;
            }
            else
            {
                ++nInflowPoints_;     // edge or corner from boundary data, below
            }
        }
    }

    // faces first, then edges, then corners: an inflow edge reads the new faces,
    // an inflow corner reads the new edges (two exchanges)
    for (label stage = 0; stage < 3; ++stage)
    {
        const label kLo = stage == 0 ? 1 : (stage == 1 ? 4 : 7);
        const label kHi = stage == 0 ? 3 : (stage == 1 ? 6 : 7);
        if (stage > 0)
        {
            for (label k = 1; k < 8; ++k)
            {
                state_[k].correctBoundaryConditions();
                stencil_.collectData(state_[k], stState_[k]);
            }
            stencil_.collectData(psi, stPsi_);
            List<FixedList<scalar, 3>> dummy;
            forAll(tensor_, c)
            {
                if (!tensor_[c]) { continue; }
                for (label k = kLo; k <= kHi; ++k)
                {
                    // owned and not yet set: from the boundary data at the new time
                    bool owned = true;
                    if (k == 4) { owned = side_[c][0] && side_[c][1]; }
                    else if (k == 5) { owned = side_[c][0] && side_[c][2]; }
                    else if (k == 6) { owned = side_[c][1] && side_[c][2]; }
                    else if (k == 7) { owned = side_[c][0] && side_[c][1] && side_[c][2]; }
                    if (owned && nv[c][k] == unset)
                    {
                        nv[c][k] = inflowPointValue(c, k, psi, dummy);
                    }
                }
            }
        }
        forAll(tensor_, c)
        {
            if (!tensor_[c]) { continue; }
            for (label k = kLo; k <= kHi; ++k)
            {
                if (nv[c][k] != unset) { state_[k].primitiveFieldRef()[c] = nv[c][k]; }
            }
        }
    }

    if (!wallReported_ || mesh_.time().writeTime())
    {
        label a = nAdvectedPoints_, i = nInflowPoints_;
        reduce(a, sumOp<label>());
        reduce(i, sumOp<label>());
        Info<< "cubicLatticeWallTensor: boundary points advected " << a
            << ", from boundary data " << i << endl;
        wallReported_ = true;
    }
}

// ************************************************************************* //
