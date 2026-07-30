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
    Foam::signedDistanceLinearWeightedLeastSquaresReconstruction

Description
    SIGNED-DISTANCE-PRESERVING linear weighted-least-squares reconstruction:
    the same per-cell IDW-weighted linear VALUE fit as
    linearWeightedLeastSquares,

        P(x) = psi_c + g . (x - x_c),

    but the departure-foot value is the fit's LOCAL SIGNED DISTANCE

        psi(x_foot) = P(x_foot) / |g|      (SMALL-guarded),

    i.e. the plane is renormalised to a UNIT gradient at evaluation. Two
    consequences, both by construction:

      1. |grad psi| is PINNED to ~1 through the transport step -- the runaway-
         gradient mode that destabilises the plain linear value fit
         (|grad psi| -> 1e21) and the gradedness erosion/steepening measured in
         every coupled two-phase failure (band min|grad psi| 0.998 -> 0.27 on
         the N=128 static SDF; 2932 -> 10150 on the moving algebraic ellipse)
         are projected out at the transport level.
      2. The plane's ZERO SET is untouched (P = 0 <=> P/|g| = 0), so the
         interface motion is identical to the plain linear fit's -- this is
         kinematic regularity MAINTENANCE, not a redistancing step.

    GUARDS (both measured necessary): the division psi = P/max(|g|, sdfGradMin)
    (dict entry, default 0.1) is applied ONLY in interface-BAND donor cells (a
    psi sign change within the stencil), where the plane is a valid single-
    interface model. Two failed designs inform this: (1) a GRADIENT trust
    window switches the division where val is LARGE -> O(1/gMin) value jumps
    ringing the SDF skeleton -> coupled run dead in ~10 steps; (2) ALWAYS
    dividing amplifies wherever the plane is an invalid interface model
    (skeleton, under-resolved filaments: stencil spans both interface sides,
    |g| -> 0 with val != 0) -> |grad psi| ~ 8e14 on the stretched 2D vortex.
    The band-edge switch cuts where val ~ 0, so the branches agree at the seam.

    NO CURVATURE: linear fit, zero Hessian (as linearWeightedLeastSquares);
    two-phase runs must take curvature from an FVM path or a force model that
    reads psi crossings directly (integralSurfaceTension).

    Selected with `reconstruction signedDistanceLinearWeightedLeastSquares`.

SourceFiles
    signedDistanceLinearWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef signedDistanceLinearWeightedLeastSquaresReconstruction_H
#define signedDistanceLinearWeightedLeastSquaresReconstruction_H

#include "linearWeightedLeastSquaresReconstruction.H"

namespace Foam
{

class signedDistanceLinearWeightedLeastSquaresReconstruction
:
    public linearWeightedLeastSquaresReconstruction
{
    //- Lower clamp on the fitted |g| in the signed-distance division (guards
    //  the SDF skeleton, where the fitted gradient degenerates; continuous).
    const scalar sdfGradMin_;

public:

    TypeName("signedDistanceLinearWeightedLeastSquares");

    explicit signedDistanceLinearWeightedLeastSquaresReconstruction
    (
        const fvMesh& mesh
    );

    virtual ~signedDistanceLinearWeightedLeastSquaresReconstruction() = default;

    //- Plane value renormalised to unit gradient (see class Description).
    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
