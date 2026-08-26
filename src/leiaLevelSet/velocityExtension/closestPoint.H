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
