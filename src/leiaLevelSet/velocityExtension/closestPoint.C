/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "closestPoint.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "DynamicList.H"
#include "FixedList.H"
#include "zoneDistribute.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(closestPoint, 0);
    addToRunTimeSelectionTable(velocityExtension, closestPoint, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::closestPoint::closestPoint(const fvMesh& mesh)
:
    steadyUpwind(mesh),
    cpTol_(velExtDict_.getOrDefault<scalar>("cpTol", 0.1)),
    cpMaxDist_(velExtDict_.getOrDefault<scalar>("cpMaxDist", 1.5)),
    cpHaloReach_(velExtDict_.getOrDefault<scalar>("cpHaloReach", 1.5)),
    cpKnownVicinity_(velExtDict_.getOrDefault<bool>("cpKnownVicinity", true)),
    cpVerifySearch_(velExtDict_.getOrDefault<bool>("cpVerifySearch", false)),
    cpWalkMaxSteps_(velExtDict_.getOrDefault<label>("cpWalkMaxSteps", 0)),
    cpSeedReach_(velExtDict_.getOrDefault<scalar>("cpSeedReach", 4.0)),
    footCell_(),
    footPoint_(),
    nLocalise_(0),
    nWalkHit_(0),
    nReseed_(0),
    nOctree_(0),
    nWalkStep_(0),
    nNonConvex_(0),
    nOctreeTot_(0),
    nLocaliseTot_(0),
    footPointSource_
    (
        velExtDict_.getOrDefault<word>("footPointSource", "fieldInterpolation")
    ),
    recon_(nullptr)
{
    if (cpWalkMaxSteps_ <= 0)
    {
        // The descent's own path cap already bounds every query point by
        // |x_new - C_c| <= cpMaxDist_*nLayers_*h, and an ACCEPTED walk step
        // moves at least one cell towards the point, so cpMaxDist_*nLayers_
        // steps suffice whenever the walk path is straight. 4x that covers
        // walks that have to work around a concave region (the paper's Fig. 3
        // walks a tetrahedral mesh around a spherical obstacle), and the
        // floor of 16 keeps the budget sane for the usual nLayers_ = 3.
        cpWalkMaxSteps_ =
            Foam::max(label(16), label(4*cpMaxDist_*nLayers_ + 8));
    }

    if (cpKnownVicinity_)
    {
        Info<< "closestPoint: LENT known-vicinity search ON"
            << " (walk budget " << cpWalkMaxSteps_ << " steps, cached-seed"
            << " reach " << cpSeedReach_ << " cells); octree = fallback only."
            << endl;
    }
    else
    {
        Info<< "closestPoint: known-vicinity search OFF -- every descent"
            << " iterate re-localised by the UNSEEDED octree"
            << " (polyMesh::findCell), the shipped path." << endl;
    }

    if (cpVerifySearch_)
    {
        WarningInFunction
            << "cpVerifySearch is ON: the walk AND polyMesh::findCell run on"
            << " every re-localisation and the run FatalErrors on the first"
            << " disagreement. This is SLOWER than the plain octree path."
            << " Regression instrument only -- switch it off for production."
            << endl;
    }

    if (footPointSource_ == "reconstructionFit")
    {
        // Standalone quadratic LSQ fit, the same object and the same
        // construction path the Eulerian two-phase solver already uses for the
        // production curvature. No semi-Lagrangian file is touched, and the
        // class lives in this library, so nothing in the build changes.
        recon_ = slReconstruction::New(mesh_);

        Info<< "closestPoint: foot points from the QUADRATIC LSQ FIT"
            << " (slReconstruction::evaluateRaw), not from the interpolated"
            << " raw psi." << endl;
    }
    else if (footPointSource_ != "fieldInterpolation")
    {
        FatalErrorInFunction
            << "Unknown footPointSource '" << footPointSource_ << "'." << nl
            << "Valid: fieldInterpolation (default) | reconstructionFit."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::closestPoint::knownVicinityWalk
(
    const point& p,
    const label seed
)
{
    if (seed < 0 || seed >= mesh_.nCells())
    {
        return -1;
    }

    const cellList& cs = mesh_.cells();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const label nInt = mesh_.nInternalFaces();
    const volVectorField& C = mesh_.C();

    label c = seed;

    // The walk's descent quantity: squared distance from p to the centre of
    // the current cell. Squared, not mag: the comparison is monotone in it
    // and this saves one sqrt per candidate neighbour, of which there are
    // (walk steps) x (faces per cell) -- the inner loop of the whole model.
    scalar dc = Foam::magSqr(C[c] - p);

    for (label step = 0; step < cpWalkMaxSteps_; ++step)
    {
        ++nWalkStep_;

        // TERMINATION: point-in-polyhedron on the candidate cell.
        //
        // THE TRAP THE LENT PAPER CALLS OUT EXPLICITLY (Sec. 2.1.1): "the
        // important difference between the cell and the polyhedron [...] is
        // the assumption that the polyhedron has only outward directed
        // normals, which is not true for the mesh cell within the OpenFOAM
        // framework [...] face owner is defined as the cell with the smaller
        // cell ID with respect to the face neighbor cell. Therefore, the face
        // area normal vectors of a cell are not oriented only outward or only
        // inward." Sf points OWNER -> NEIGHBOUR, so for cell c it is outward
        // only on the faces c OWNS; on the faces where c is the NEIGHBOUR it
        // points INTO c, and a textbook half-space test that assumes outward
        // normals mis-classifies exactly those faces -- typically half of
        // them, which makes the test reject almost every cell and would send
        // every single query to the octree fallback (the optimisation would
        // silently do nothing while still being "correct").
        //
        // AVOIDED by delegating to polyMesh::pointInCell with FACE_PLANES,
        // which IS that faceOwner-oriented test and nothing else -- see
        // primitiveMeshFindCell.C:
        //     const labelList& owner = this->faceOwner();
        //     vector normal = Sf[nFace];
        //     if (owner[nFace] != celli) { normal = -normal; }
        //     if ((normal & proj) > 0) { return false; }
        // Delegating rather than re-deriving the flip here means this test
        // cannot drift away from the orientation convention OpenFOAM itself
        // uses in findCell.
        //
        // TWO-STAGE, and the second stage is what buys bit-identity.
        // FACE_PLANES is ~nFaces dot products and rejects essentially every
        // non-containing cell for the price of one; it is the cheap probe run
        // at every walk step. Acceptance is then confirmed with the DEFAULT
        // CELL_TETS decomposition, which is bit-exactly the predicate the
        // octree applies to its own candidates (polyMesh::cellTree() builds
        // treeDataCell(false, mesh, CELL_TETS) and treeDataCell::contains()
        // forwards to mesh.pointInCell(sample, celli, CELL_TETS)). Because
        // the walk only ever returns a cell that satisfies the octree's own
        // containment predicate, and defers to the octree otherwise, the cell
        // index handed back to the descent is the same one findCell returns.
        if (mesh_.pointInCell(p, c, polyMesh::FACE_PLANES))
        {
            if (mesh_.pointInCell(p, c))    // CELL_TETS: the octree's test
            {
                return c;
            }

            // Planes accept, tets reject. Only reachable for a NON-CONVEX
            // cell (p sits in the re-entrant region that is outside the cell
            // but inside every face plane) or for p within roundoff of a
            // face. The closest-centre rule has no better candidate to offer
            // from here, so hand the query to the octree and keep the answer
            // exactly findCell's. Counted separately: on hex/tet meshes this
            // must stay at zero, and a non-zero count on a polyhedral mesh is
            // the signal that the walk is paying for the fallback often.
            ++nNonConvex_;
            return -1;
        }

        // STEP (Fig. 2): to the face-neighbour whose CENTRE is closest to p.
        // Only a STRICT decrease of magSqr(C - p) is accepted, so the descent
        // quantity is strictly monotone along the walk: the path can never
        // revisit a cell and the loop always terminates -- inside a cell, at
        // a local minimum, or on the step budget. No visited-cell set is
        // needed, which is what keeps the step cost at nFaces dot products.
        label cBest = -1;
        for (const label f : cs[c])
        {
            if (f >= nInt)
            {
                // Boundary or COUPLED (processor / cyclic) face. The cell on
                // the far side is not addressable on this rank -- for a
                // processor patch it lives in another rank's mesh, and for a
                // cyclic patch its centre is the UNTRANSFORMED one, so the
                // distance comparison against p would be meaningless. The
                // walk simply cannot step through it; it stalls, the octree
                // answers, and for a point that genuinely left the rank that
                // answer is -1 -- exactly the condition the existing
                // halo/Taylor branch below already keys on.
                continue;
            }

            const label cn = (own[f] == c) ? nei[f] : own[f];
            const scalar d = Foam::magSqr(C[cn] - p);
            if (d < dc)
            {
                dc = d;
                cBest = cn;
            }
        }

        if (cBest < 0)
        {
            // Local minimum of the centre distance and p is not inside: the
            // classic known-vicinity failure (concave domain, obstacle
            // between seed and point -- the paper's Fig. 3 -- or p outside
            // the local mesh altogether).
            return -1;
        }

        c = cBest;
    }

    return -1;   // step budget exhausted
}


Foam::label Foam::closestPoint::findCellSeeded
(
    const point& p,
    const label seed,
    const label seed2
)
{
    ++nLocalise_;
    ++nLocaliseTot_;

    if (!cpKnownVicinity_)
    {
        ++nOctree_;
        ++nOctreeTot_;
        return mesh_.findCell(p);
    }

    label c = knownVicinityWalk(p, seed);

    if (c < 0 && seed2 >= 0 && seed2 != seed)
    {
        // Second chance from the band cell itself before paying for the
        // octree. This is the case where the CACHED seed has gone stale (the
        // interface moved away from where it was last step, or the walk from
        // it ran into a concave region) while the band cell c -- which is by
        // construction within cpMaxDist_*nLayers_ cells of its own foot point
        // -- is still a good seed. Cheap: a failed walk costs at most
        // cpWalkMaxSteps_ x nFaces dot products, against an octree descent
        // over the whole mesh.
        c = knownVicinityWalk(p, seed2);
        if (c >= 0)
        {
            ++nReseed_;
        }
    }

    if (c < 0)
    {
        // The paper's contract: "the octree-based search is always used as a
        // fall-back solution", which is what makes the combination converge
        // unconditionally for any point inside the flow domain.
        ++nOctree_;
        ++nOctreeTot_;
        return mesh_.findCell(p);
    }

    ++nWalkHit_;

    if (cpVerifySearch_)
    {
        const label cRef = mesh_.findCell(p);
        if (c != cRef)
        {
            FatalErrorInFunction
                << "known-vicinity search DISAGREES with polyMesh::findCell"
                << " at p = " << p << nl
                << "    walk   -> cell " << c
                << " (centre " << mesh_.C()[c] << ")" << nl
                << "    octree -> cell " << cRef << nl
                << "The walk accepts a cell only when"
                << " polyMesh::pointInCell(p, c, CELL_TETS) is true, which is"
                << " the octree's own containment predicate, so this can only"
                << " happen for a point lying exactly on a face shared by two"
                << " cells that both contain it."
                << exit(FatalError);
        }
    }

    return c;
}


void Foam::closestPoint::extend()
{
    interpolationCellPoint<scalar> psiInterp(psi_);
    interpolationCellPoint<vector> Uinterp(U_);

    // Rebuild the local fits against the CURRENT psi. The descent below then
    // reads psi from a smooth polynomial surrogate rather than from the
    // cell-point interpolation of the raw field.
    if (recon_)
    {
        recon_->update(psi_);
    }
    const auto& C = mesh_.C();

    // M_c (Eqs. (15)/(16)) allocated on first use and RETAINED across time
    // steps -- that persistence is the whole optimisation: it is what turns
    // the octree into the preprocessing step of Sec. 2.1.1 instead of a
    // per-iterate cost. Sized to the MESH, not to the band, because band
    // membership changes every step as the interface moves; a cell that
    // leaves the band keeps its entry and is seeded from it when it comes
    // back.
    //
    // The size test is also the topology-change guard: a refined,
    // redistributed or reread mesh reallocates the map and every seed
    // reverts to the band cell itself. Even a stale index that slipped past
    // this test cannot corrupt anything -- see knownVicinityWalk(): a seed is
    // only a starting point, and containment is PROVEN before any cell is
    // accepted, so a bad seed costs one octree fallback and nothing else.
    if (footCell_.size() != mesh_.nCells())
    {
        footCell_ = labelList(mesh_.nCells(), -1);
        footPoint_ = pointField(mesh_.nCells(), point::zero);
    }

    nLocalise_ = 0;
    nWalkHit_ = 0;
    nReseed_ = 0;
    nOctree_ = 0;
    nWalkStep_ = 0;
    nNonConvex_ = 0;

    // Build the face-diagonal (tet) decomposition ONCE, HERE, where every
    // rank is. polyMesh::tetBasePtIs() is lazily constructed by
    // polyMeshTetDecomposition::findFaceBasePts(), which ends in
    // syncTools::swapBoundaryFacePositions() and syncTools::syncFaceList():
    // both COLLECTIVE. polyMesh::findCell() knows this -- it calls
    // (void)tetBasePtIs() before it can early-return on an empty rank -- but
    // that only protects the ranks that actually reach findCell. The descent
    // loop below is rank-local: a rank whose band is empty this step calls
    // neither the CELL_TETS acceptance test of the walk nor the octree
    // fallback, so lazily constructing the decomposition down there would
    // have some ranks entering a collective that the others never enter.
    // Hoisting it to a point every rank passes every call removes that
    // hazard; when the decomposition already exists this is a pointer test.
    if (Pstream::parRun())
    {
        (void)mesh_.tetBasePtIs();
    }

    // Parallel finishing step (zoneDistribute; Scheufler & Roenby JCP 383
    // (2019), production-tested in geometricVoF/plicRDF): fetch ONE
    // cell-point-cell halo layer of remote (C, psi, grad psi, U, grad U)
    // around the band. A descent step that leaves the processor-local mesh
    // (findCell < 0) continues on first-order Taylor data anchored at the
    // nearest halo cell centre -- the same order as the local cell-point
    // interpolation -- instead of failing. Only excursions deeper than the
    // one-layer halo, physical-boundary overshoots and cyclic transforms
    // (whose untransformed centres misrepresent distances, so the reach cap
    // rejects them) still fall back to the pinned steadyUpwind fill.
    const bool haveHalo = Pstream::parRun();
    zoneDistribute* zd = nullptr;
    const labelListList* stencilPtr = nullptr;
    Map<vector> mapC;
    Map<scalar> mapPsi;
    Map<vector> mapGPsi;
    Map<vector> mapU;
    Map<tensor> mapGU;

    if (haveHalo)
    {
        boolList zone(mesh_.nCells(), false);
        forAll(band_, c)
        {
            zone[c] = band_[c] > 0.5;
        }

        zd = &zoneDistribute::New(mesh_);
        zd->setUpCommforZone(zone, true);

        // The same gradient flavours the family already uses: grad(psi)
        // resolves the fvSchemes 'grad(psi)' entry (as computeNormals);
        // grad(U) resolves the gradSchemes default.
        tmp<volVectorField> tgPsi = fvc::grad(psi_);
        tmp<volTensorField> tgU = fvc::grad(U_);

        mapC = zd->getDatafromOtherProc(zone, C);
        mapPsi = zd->getDatafromOtherProc(zone, psi_);
        mapGPsi = zd->getDatafromOtherProc(zone, tgPsi());
        mapU = zd->getDatafromOtherProc(zone, U_);
        mapGU = zd->getDatafromOtherProc(zone, tgU());
        stencilPtr = &zd->getStencil();
    }

    // Const views: Map::operator[] const FatalErrors on a missing key
    // instead of silently inserting.
    const Map<vector>& ccC = mapC;
    const Map<scalar>& ccPsi = mapPsi;
    const Map<vector>& ccGPsi = mapGPsi;
    const Map<vector>& ccU = mapU;
    const Map<tensor>& ccGU = mapGU;

    // Nearest REMOTE halo cell (global index) to point p, searched in the
    // stencil of the last local cell the walk visited. -1 if none lies
    // within cpHaloReach_*cellSize: one halo layer cannot represent points
    // beyond the immediate neighbour cells -- deeper walks must fall back.
    const auto nearestHalo = [&](const point& p, const label cLast) -> label
    {
        if (!haveHalo || cLast < 0)
        {
            return -1;
        }
        const labelList& st = (*stencilPtr)[cLast];
        const globalIndex& gn = zd->globalNumbering();
        label best = -1;
        scalar bestD = cpHaloReach_*cellSize_[cLast];
        for (const label gbl : st)
        {
            if (gn.isLocal(gbl))
            {
                continue;   // local containment already failed for p
            }
            const auto fnd = ccC.cfind(gbl);
            if (!fnd.good())
            {
                continue;
            }
            const scalar d = Foam::mag(p - *fnd);
            if (d < bestD)
            {
                bestD = d;
                best = gbl;
            }
        }
        return best;
    };

    // Every band cell gets the interface velocity at ITS closest interface
    // point: Newton descent along the LOCAL normal (re-evaluated at the cell
    // currently containing the iterate -- across the band the frozen starting
    // normal would be too crude), with acceptance guards. Successful cells
    // become Dirichlet data; failures (skeleton region, path cap, halo-reach
    // exceedance) are filled afterwards by the steadyUpwind solve pinned at
    // every successful cell.
    DynamicList<label> fixedCells(mesh_.nCells()/4);
    DynamicList<vector> fixedVals(mesh_.nCells()/4);
    label nFailed = 0;
    label nHalo = 0;

    forAll(band_, c)
    {
        if (band_[c] < 0.5)
        {
            continue;
        }

        const scalar hc = cellSize_[c];
        const scalar maxDist = cpMaxDist_*nLayers_*hc;
        // Taylor-vs-interpolation operator jumps at a processor crossing are
        // O(h^2); allow them without tripping the divergence guard.
        const scalar crossSlack = 1e-3*hc;

        point x = C[c];
        label cx = c;        // containing local cell (-1: halo/Taylor mode)
        label jg = -1;       // halo data cell (global index) when cx == -1
        label cLast = c;     // last local cell: stencil access + reach scale
        scalar apsi = Foam::mag(psi_[c]);
        bool ok = true;
        bool usedHalo = false;

        // Parallel: cap each step at ~one cell so the walk TRAVERSES the
        // cells on its way (cLast tracks to the true exit cell, whose halo
        // stencil covers the first remote layer at the crossing point). A
        // full Newton jump from an outer band cell would anchor the halo
        // search at the distant STARTING cell and miss. The capped walk
        // needs a larger iteration budget (~nLayers steps to reach the
        // interface + Newton polish). Serial keeps the exact original
        // full-jump iteration (bit-identical regression).
        const label nIter = haveHalo ? nDescent_ + nLayers_ + 2 : nDescent_;

        for (label it = 0; it < nIter; ++it)
        {
            scalar psix;
            vector nx;
            if (cx >= 0)
            {
                // SINGLE-VARIABLE SWAP: only WHERE psi is read changes. The
                // normal stays the same discrete nHat_ and every guard,
                // step cap and iteration budget is untouched, so a difference
                // in the result is attributable to the foot-point DEFINITION
                // and to nothing else.
                psix =
                    recon_
                  ? recon_->evaluateRaw(cx, x)
                  : psiInterp.interpolate(x, cx);
                nx = nHat_[cx];
            }
            else
            {
                const vector& Cj = ccC[jg];
                const vector& gj = ccGPsi[jg];
                psix = ccPsi[jg] + (gj & (x - Cj));
                nx = gj/Foam::max(Foam::mag(gj), SMALL);
            }

            vector dx = -psix*nx;
            if (haveHalo)
            {
                const scalar maxStep =
                    0.9*cellSize_[cLast >= 0 ? cLast : c];
                const scalar m = Foam::mag(dx);
                if (m > maxStep)
                {
                    dx *= maxStep/m;
                }
            }
            const point xNew = x + dx;

            if (Foam::mag(xNew - C[c]) > maxDist)
            {
                ok = false;   // walked too far (skeleton / bad normal)
                break;
            }

            // Seed for the known-vicinity search (LENT Sec. 2.1.1).
            //
            //  * Iterates > 0: the cell holding the PREVIOUS iterate. In
            //    parallel the step is capped at 0.9 cells, so this seed is
            //    the containing cell or one of its face neighbours and the
            //    walk terminates in one or two steps. cLast is used when the
            //    previous iterate left the rank (cx < 0) -- it is the last
            //    LOCAL cell the walk visited, i.e. the exit cell at the
            //    processor boundary, and therefore the only local index that
            //    means anything for a point out there.
            //
            //  * Iterate 0: the FULL Newton jump from C[c], which in serial
            //    spans up to |psi(C_c)| ~ nLayers_ cells and is the one query
            //    the octree used to be genuinely needed for. Seeded from the
            //    CACHED foot cell M_c(c) of the previous time step -- Eq.
            //    (16), "M_c(T) is always kept up-to-date", re-seeded from its
            //    own previous value. Between two time steps the foot point
            //    moves only O(|U| dt) = O(CFL h), so the cached cell is one
            //    or two face crossings from the new one while c itself is
            //    O(nLayers_) crossings away.
            //
            // The cpSeedReach_ test is what makes the cache SELF-VALIDATING
            // without any state beyond the two arrays: it compares the new
            // query point against the cached foot POINT (coordinates, which
            // are rank-independent and always meaningful) and rejects the
            // hint when they have drifted apart -- first call, restart, a
            // cell that has just entered the band, a foot point that has
            // since crossed onto another rank. On rejection the seed is c,
            // and findCellSeeded() retries from c anyway before the octree.
            label seed = (cx >= 0) ? cx : cLast;
            if
            (
                cpKnownVicinity_
             && it == 0
             && footCell_[c] >= 0
             && Foam::mag(xNew - footPoint_[c]) < cpSeedReach_*hc
            )
            {
                seed = footCell_[c];
            }

            // Re-localise: prefer the local mesh (a walk may re-enter it),
            // then the one-layer halo, then give up (deep excursion,
            // physical-boundary overshoot, cyclic transform).
            //
            // findCellSeeded() returns EXACTLY what mesh_.findCell(xNew)
            // returns -- walk first, octree whenever the walk fails -- so
            // every guard, step cap, budget and branch below is untouched and
            // reads the same cell index it read before.
            const label cn = findCellSeeded(xNew, seed, c);
            label jn = -1;
            if (cn < 0)
            {
                jn = nearestHalo(xNew, cLast);
                if (jn < 0)
                {
                    ok = false;
                    break;
                }
            }

            scalar apsiNew;
            if (cn >= 0)
            {
                apsiNew = Foam::mag(psiInterp.interpolate(xNew, cn));
            }
            else
            {
                apsiNew = Foam::mag
                (
                    ccPsi[jn] + (ccGPsi[jn] & (xNew - ccC[jn]))
                );
            }

            const bool pureLocal = (cx >= 0 && cn >= 0);
            if (apsiNew > apsi + (pureLocal ? 0.0 : crossSlack))
            {
                ok = false;   // diverging descent
                break;
            }

            x = xNew;
            cx = cn;
            jg = jn;
            apsi = apsiNew;
            if (cn >= 0)
            {
                cLast = cn;
            }
            else
            {
                usedHalo = true;
            }
        }

        if (ok && apsi < cpTol_*hc)
        {
            if (cx >= 0)
            {
                Uext_[c] = Uinterp.interpolate(x, cx);
            }
            else
            {
                // First-order Taylor of U at the foot point, anchored at the
                // halo data cell: U_j + (x - C_j) . (grad U)_j -- the same
                // order as the local linear (cell-point) interpolation.
                Uext_[c] = ccU[jg] + ((x - ccC[jg]) & ccGU[jg]);
            }
            // Update M_c (Eq. (16)) for the next time step.
            //
            // The foot POINT is cached unconditionally: coordinates are
            // rank-INDEPENDENT, so they stay meaningful no matter which rank
            // the foot point ends up on, and they are what the cpSeedReach_
            // test above vets the cell hint against.
            //
            // The cell index is cached as a purely LOCAL SEED HINT, valid
            // only on this rank. It is never written into a field, never
            // halo-exchanged and never sent anywhere: a local cell index that
            // crosses a processor boundary denotes a DIFFERENT cell on the
            // receiving rank, which is why the obvious "store the foot cell
            // in a volScalarField so it exchanges" design is wrong. Here
            // footCell_[c] is always an index THIS rank computed for ITS OWN
            // band cell c, so it is only ever dereferenced where it is valid.
            //
            // When the foot point ended on the halo/Taylor path (cx < 0)
            // there is no local containing cell, and cLast -- the exit cell
            // at the processor boundary -- is cached instead: next step the
            // walk starts there and marches back out towards the same region,
            // stalls at the same coupled face and hands over to the octree,
            // which is exactly the shipped behaviour for that branch.
            footPoint_[c] = x;
            footCell_[c] = (cx >= 0) ? cx : cLast;

            if (usedHalo)
            {
                ++nHalo;
            }
            fixedCells.append(c);
            fixedVals.append(Uext_[c]);
        }
        else
        {
            ++nFailed;
        }
    }

    // ONE collective for every counter this model reports. The three original
    // reductions plus the eight search counters would be eleven Allreduce
    // calls per time step on a routine whose work is otherwise entirely
    // rank-local -- packing them into a single contiguous FixedList makes it
    // one. The reduced values are bit-identical to the separate reductions
    // (integer sums over the same ranks), and in particular nFailed stays the
    // GLOBAL count, which is what makes the collective solveSteady() below
    // run on every rank whenever ANY rank has a failed cell.
    FixedList<label, 11> cnt;
    cnt[0] = fixedCells.size();
    cnt[1] = nHalo;
    cnt[2] = nFailed;
    cnt[3] = nLocalise_;
    cnt[4] = nWalkHit_;
    cnt[5] = nReseed_;
    cnt[6] = nOctree_;
    cnt[7] = nWalkStep_;
    cnt[8] = nNonConvex_;
    cnt[9] = nOctreeTot_;
    cnt[10] = nLocaliseTot_;
    reduce(cnt, sumOp<label>());

    const label nFixed = cnt[0];
    nHalo = cnt[1];
    nFailed = cnt[2];

    Info<< "closestPoint: " << nFixed << " foot-pointed band cells ("
        << nHalo << " via halo), " << nFailed << " fallback cells" << endl;

    // Search accounting, so the optimisation is MEASURED and not assumed.
    // The number that matters is the octree fraction: Sec. 2.1.1 claims the
    // octree degenerates to a preprocessing step, which here means the
    // cumulative octree count should stop growing after the first time step
    // apart from the parallel excursions that legitimately need it (those are
    // bounded by nHalo above plus the true off-rank misses). If the per-step
    // octree count stays comparable to the re-localisation count, the walk is
    // failing and the SEEDS -- not the search -- are what to look at.
    const label nLoc = cnt[3];
    const label nHit = cnt[4];
    const label nRes = cnt[5];
    const label nOct = cnt[6];
    const label nStp = cnt[7];
    const label nNcv = cnt[8];
    const label nOctT = cnt[9];
    const label nLocT = cnt[10];

    const scalar octPct = nLoc > 0 ? 100.0*nOct/nLoc : 0.0;
    const scalar stepAvg = nHit > 0 ? scalar(nStp)/nHit : 0.0;
    const scalar octPctTot = nLocT > 0 ? 100.0*nOctT/nLocT : 0.0;

    Info<< "closestPoint search: " << nLoc << " re-localisations, "
        << nHit << " known-vicinity hits (" << nRes << " after re-seeding, "
        << stepAvg << " steps/hit), " << nOct << " octree fallbacks ("
        << octPct << "%), " << nNcv << " plane/tet mismatches; cumulative "
        << nOctT << "/" << nLocT << " octree (" << octPctTot << "%)" << endl;

    // Fill failures (and keep the far field consistent) with the steady
    // upwind transport pinned at ALL successful closest-point cells.
    if (nFailed > 0)
    {
        solveSteady(Uext_, fixedCells, fixedVals);
    }
}

// ************************************************************************* //
