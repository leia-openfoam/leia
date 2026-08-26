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

Class
    Foam::closestPoint

Description
    Direct (geometric) closest-point velocity extension over the whole band:

        Uext(x_c) = U(x_fp(x_c)),   x_fp = foot point of x_c on Sigma,

    so (n . grad) Uext = 0 holds BY CONSTRUCTION (not asymptotically) inside
    the foot-point reach -- the general guarantee for unstructured FVM
    (closest point method: Ruuth & Merriman JCP 227 (2008); degree-p
    interpolation gives O(h^{p+1}): Macdonald & Ruuth SISC 31 (2009);
    front-tracking practice: Maric et al. Computers & Fluids 113 (2015)).

    Foot points by Newton descent along the local normal with cell
    re-localization (the interface-layer machinery generalized to every band
    cell) with acceptance guards: |psi(x_fp)| < cpTol*cellSize, non-decreasing
    |psi| abort, and path cap cpMaxDist*nLayers*cellSize.

    Parallel: a zoneDistribute cell-point-cell halo (Scheufler & Roenby
    JCP 383 (2019); production-tested in geometricVoF/plicRDF) carries one
    layer of remote (C, psi, grad psi, U, grad U) around the band. A descent
    step that leaves the processor-local mesh continues on first-order Taylor
    data anchored at the nearest halo cell centre (same order as the local
    cell-point interpolation) within cpHaloReach*cellSize of the exit cell.
    Only deeper excursions, cyclic transforms and the geometric guard
    failures fall back to the steadyUpwind solve pinned at all successful
    cells -- parallel-graceful by construction.

    SEARCH. Every descent iterate has to be re-localised (which cell contains
    it?). Done with polyMesh::findCell that is an UNSEEDED octree descent over
    the whole mesh, executed (band cells) x (nDescent, or nDescent + nLayers +
    2 in parallel) times per time step; MEASURED at 20.5x the wall clock of
    plain Euler advection on the 2D vortex at N = 512. Replaced by the LENT
    known-vicinity search (Maric et al., Computers & Fluids 113 (2015) 20-31,
    Sec. 2.1.1, Figs. 2 and 3): a Loehner-style walk on OpenFOAM's
    owner-neighbour addressing that starts from a SEED CELL, steps to the
    face-neighbour whose CENTRE is closest to the query point, and terminates
    on a point-in-polyhedron test. Seeds come from the cached foot-cell map
    M_c (Eqs. (15)/(16)), re-seeded from its own previous value. The octree is
    retained STRICTLY as the fallback, as in the paper: "octree-based search is
    executed for the whole computational domain only once, as a preprocessing
    step. All the other invocations, provided that the known-vicinity search
    fails in the first place, will start within an octree node defined by a
    previously cached cell."

SourceFiles
    closestPoint.C

\*---------------------------------------------------------------------------*/

#ifndef closestPoint_H
#define closestPoint_H

#include "steadyUpwind.H"
#include "slReconstruction.H"

namespace Foam
{

class closestPoint
:
    public steadyUpwind
{
    //- Foot-point acceptance tolerance: |psi(x_fp)| < cpTol_*cellSize.
    scalar cpTol_;

    //- Path cap: |x_fp - x_c| <= cpMaxDist_*nLayers*cellSize.
    scalar cpMaxDist_;

    //- Halo reach: Taylor data from a remote (one-layer) halo cell is valid
    //- within cpHaloReach_*cellSize of the exit cell; beyond -> fallback.
    scalar cpHaloReach_;

    //- Use the LENT known-vicinity walk to re-localise the descent iterates,
    //  with polyMesh::findCell (the octree) kept only as the fallback.
    //
    //  DEFAULT true. This is an EFFICIENCY switch, not a model switch: the
    //  walk is accepted only when polyMesh::pointInCell(p, c) -- with the
    //  DEFAULT CELL_TETS decomposition, which is bit-exactly the predicate
    //  the octree itself applies (polyMesh::cellTree() constructs
    //  treeDataCell(false, mesh, CELL_TETS), and treeDataCell::contains()
    //  forwards to mesh.pointInCell(sample, celli, decompMode_)) -- says the
    //  query point is inside that cell, and EVERY walk failure defers to
    //  findCell. A point strictly inside exactly one cell therefore gets the
    //  same cell index from both paths, so the arithmetic downstream is
    //  unchanged bit for bit.
    //
    //  Set false to reproduce the shipped octree path exactly, line for line:
    //  findCellSeeded() then degenerates to a bare mesh_.findCell(p) call and
    //  neither the cache nor the walk is consulted. Use it as the regression
    //  reference (see cpVerifySearch below).
    bool cpKnownVicinity_;

    //- Regression instrument: run the walk AND polyMesh::findCell on EVERY
    //  re-localisation and FatalError on the first disagreement. Strictly
    //  SLOWER than the octree path (it pays for both), so it is off by
    //  default -- it exists so bit-identity is MEASURED on the target case
    //  rather than argued for. The only construct that can make the two
    //  disagree is a query point lying exactly on a shared tet face, where
    //  both adjacent cells contain it and the two searches may pick
    //  different ones; cell-point interpolation of psi and U is continuous
    //  across that face (the cell-centre barycentric weight vanishes there,
    //  leaving only the shared face-centre and point values) so even then
    //  psix and Uext are unchanged -- only the cell-wise nHat_[cx] could
    //  jump, on a measure-zero set.
    bool cpVerifySearch_;

    //- Step budget of one known-vicinity walk; <= 0 auto-sizes it in the
    //  constructor. Bounds the worst case: a walk that has not terminated
    //  within the budget gives up and the octree answers, so the walk can
    //  never be more than a bounded constant slower than findCell.
    label cpWalkMaxSteps_;

    //- Reach, in cellSize_[c], within which the CACHED foot cell of band cell
    //  c is accepted as the seed for the first (largest) descent jump. The
    //  foot point moves by O(|U| dt) = O(CFL h) between time steps, so the
    //  cached cell is one or two face crossings from the new foot point,
    //  while the band cell c itself is O(nLayers) crossings away. Beyond this
    //  reach (first call, restart, a cell that has just entered the band, a
    //  topology change) the seed reverts to c.
    scalar cpSeedReach_;

    //- M_c, the cached foot-CELL map of Eq. (16), "always kept up-to-date"
    //  and re-seeded from its own previous value. Indexed by MESH cell (not
    //  by a band-local index) because cells enter and leave the band between
    //  time steps; -1 means "no hint". footCell_[c] is a RANK-LOCAL index
    //  that THIS rank computed for ITS OWN band cell c. It is never sent
    //  anywhere, so the classic error of exchanging a local cell index across
    //  a processor boundary -- where it denotes a DIFFERENT cell on the
    //  receiving rank -- cannot occur here. A stale entry (mesh change,
    //  restart, foot point that has since crossed onto another rank) cannot
    //  produce a wrong answer, only a slower one: it is a SEED, and the walk
    //  proves containment before accepting any cell.
    labelList footCell_;

    //- The cached foot POINT of band cell c. Coordinates are rank-INDEPENDENT
    //  and therefore always meaningful, which is exactly why the cache is
    //  anchored on them: footPoint_[c] is what validates footCell_[c] as a
    //  seed (the cpSeedReach_ test), and it is what remains usable when the
    //  foot point has crossed onto another rank and no local cell contains it.
    pointField footPoint_;

    //- Search statistics for the CURRENT call (reset at the top of extend()).
    label nLocalise_;
    label nWalkHit_;
    label nReseed_;
    label nOctree_;
    label nWalkStep_;
    label nNonConvex_;

    //- Cumulative octree calls / re-localisations over the whole run: the
    //  number that says whether the octree really became a preprocessing
    //  step (Sec. 2.1.1) or is still being hit every time step.
    label nOctreeTot_;
    label nLocaliseTot_;


    //- Where the descent reads psi when locating the foot point.
    //
    //  fieldInterpolation (DEFAULT, the shipped behaviour): the zero of the
    //  cell-point-interpolated RAW psi field, interpolationCellPoint<scalar>.
    //
    //  reconstructionFit: the zero of the QUADRATIC LEAST-SQUARES FIT that the
    //  semi-Lagrangian line already builds, via slReconstruction::evaluateRaw.
    //
    //  WHY THE OPTION EXISTS. On the 1D uniaxial-stretching gate, with ONE
    //  interface and therefore a perfectly smooth closest-point map (no medial
    //  axis anywhere near the band), this model does NOT reproduce the rigid
    //  translation it should: band mean |d(psi)/dx| = 0.9495 / 0.9421 / 0.8981
    //  at N = 32/64/128 against the exact 1, DEGRADING under refinement, and it
    //  misplaces the interface by 0.128 / 0.115 / 0.073 cells where advection
    //  with the true velocity is accurate to 4e-4 / 2e-4 / 1e-4 cells -- three
    //  orders of magnitude worse on the quantity that decides the solution.
    //
    //  Sweeping cpTol_ over FOUR DECADES (0.1 -> 0.0001 cells) moved that
    //  interface error by less than 1%: 0.1153, 0.1153, 0.1153, 0.1166 cells.
    //  The error is therefore NOT in converging to the foot point, it is in
    //  what the foot point IS. Tightening the tolerance only converges harder
    //  onto the zero of an interpolated raw psi, and if that zero is displaced
    //  from the interface no tolerance can help.
    //
    //  That is the one thing the static extension gate cannot see. It measures
    //  |nhat . grad(Uext)|, which is satisfied EXACTLY by any normal-constant
    //  field including one carrying a uniformly WRONG value -- a foot-point
    //  error shifts Uext by a constant along each normal and leaves the metric
    //  at zero. The gate tests constancy along normals, never Uext = U on the
    //  interface, which is why this model can measure O(h^1.98) statically and
    //  still damage advection.
    word footPointSource_;

    //- The quadratic LSQ fit, constructed ONLY when footPointSource_ is
    //  reconstructionFit. slReconstruction::New(mesh) is already proven
    //  constructible standalone outside the semi-Lagrangian solver -- the
    //  Eulerian two-phase solver does exactly this for the production
    //  curvature (leiaLevelSetTwoPhaseFoam/reconstructedCurvatureFields.H) --
    //  and it lives in the same library as this class, so no build change is
    //  needed. Not one semi-Lagrangian file is modified.
    autoPtr<slReconstruction> recon_;

    //- LENT known-vicinity walk (Sec. 2.1.1, Fig. 2): from cell `seed`, step
    //  to the face-neighbour whose centre is closest to p until p is inside.
    //  Returns the containing cell, or -1 when the walk fails (local minimum
    //  of the centre distance, coupled/boundary face in the way, step budget
    //  exhausted, non-convex cell) -- NEVER "p is outside the mesh", which is
    //  the octree's answer to give.
    label knownVicinityWalk(const point& p, const label seed);

    //- Re-localise p, returning EXACTLY what mesh_.findCell(p) returns:
    //  the walk from `seed`, retried from `seed2`, then the octree.
    label findCellSeeded(const point& p, const label seed, const label seed2);

protected:

    virtual void extend() override;

public:

    TypeName("closestPoint");

    closestPoint(const fvMesh& mesh);

    virtual ~closestPoint() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
