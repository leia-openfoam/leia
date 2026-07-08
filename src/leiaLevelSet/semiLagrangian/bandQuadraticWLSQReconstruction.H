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
    Foam::bandQuadraticWLSQReconstruction

Description
    Narrow-band semi-Lagrangian reconstruction: the O(h^2) weighted quadratic
    least-squares fit (quadraticWLSQReconstruction) is applied ONLY to cells in
    a band around the interface,

        band = { |psi_c| <= (nLayersBand + bandGuard) * h_c },

    and the far field is SLAVED to the band: after the band cells advect, psi
    outside the band is RE-EXTENDED as a signed distance from the band (a few
    Eikonal sweeps, band cells frozen). This is standard narrow-band level-set
    practice -- advect the band, re-extend the rest -- and it is the accepted
    mild reinitialization OUTSIDE the band (the interface band itself is never
    reinitialized). Independent far-field advection is NOT used: it is unstable
    (linear reconstruction blows up under the reversed vortex) and, as the
    interface sweeps the domain, contaminates the band. With the guard layer
    every metric-band cell's stencil lies inside the fitted band, so the band
    cells' reconstruction is identical to the full quadraticWLSQ model.

    Motivation: the per-cell QR pseudo-inverse of the full model is cached for
    EVERY cell (~O(N) matrices), which exhausts memory at high resolution under
    concurrency. Here the quadratic fit is done per step for band cells only,
    with NO persistent cache, so the memory is O(band) and high-resolution runs
    stay light. Same measured accuracy on the band, far less memory.

SourceFiles
    bandQuadraticWLSQReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef bandQuadraticWLSQReconstruction_H
#define bandQuadraticWLSQReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class bandQuadraticWLSQReconstruction
:
    public slReconstruction
{
    //- Active (non-empty) mesh directions (once).
    labelList activeDirs_;

    //- Quadratic basis size over the active directions (5 in 2D, 9 in 3D).
    label ncoeffFull_;

    //- Per-cell length scale = min face-neighbour centre distance (once).
    scalarField cellSize_;

    //- Band mask (recomputed each step as the interface moves).
    boolList band_;

    //- Per-cell fitted quadratic coefficients (only for band cells; else empty).
    List<scalarList> coeffs_;

    //- Per-cell fitted basis size (>0 for band cells; 0 for far-field / degen).
    labelList ncoeff_;

    //- Band half-width in cell layers and the guard added so that every
    //  metric-band cell's stencil stays inside the fitted band.
    label nLayersBand_;
    label bandGuard_;

    bool initialized_;

    //- Fill cellSize_ (min face-neighbour centre distance).
    void computeCellSize();

    //- Set activeDirs_/ncoeffFull_ from mesh.geometricD().
    void setupBasis();

    //- Quadratic polynomial basis over the active directions.
    void basis(const vector& d, const label ncoeff, scalarList& b) const;

    //- Weighted quadratic LSQ fit for cell c (per-step, no cache).
    label fitCell
    (
        const List<vector>& X,
        const List<scalar>& s,
        scalarList& coeffs
    ) const;

public:

    TypeName("bandQuadraticWLSQ");

    explicit bandQuadraticWLSQReconstruction(const fvMesh& mesh);

    virtual ~bandQuadraticWLSQReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;

    //- Re-extend psi outside the band as a signed distance from the band
    //  (Eikonal sweeps; band cells frozen). Slaves the far field to the band.
    virtual void postAdvect(volScalarField& psi) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
