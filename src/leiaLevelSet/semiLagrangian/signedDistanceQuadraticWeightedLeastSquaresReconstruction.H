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
    Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction

Description
    Signed-distance-preserving sibling of uncachedQuadraticWeightedLeastSquares.

    It runs the SAME IDW-weighted quadratic VALUE fit (identical stencil, basis,
    normal equations and Cholesky solve -- the slope is a byproduct of the fit,
    never an injected input), but instead of returning the reconstructed value
    P0 = P(x_d) it returns the signed distance from the foot to the reconstructed
    zero set,

        d1 = P0 / |grad P(x_d)|,                                   (reproject)

    which resets |grad psi| -> 1 at every foot every step. The exponential
    magnitude drift that a plain value-copy accumulates under sustained normal
    strain (D/Dt ln|grad psi| = -n.D.n) is thereby severed: the output carries
    only this step's strain and is reset again next step. The residual is removed
    by one normal-strain rescaling that maps the old-interface distance at the foot
    to the new-interface distance at the arrival cell,

        d  = d1 * (1 + eps_nn*dt),   eps_nn = n.D.n,               (rescale)
        D  = 0.5*(grad U + grad U^T)   at the arrival cell,

    upgrading the field from near-SDF to second-order SDF (sign of the factor is
    +, i.e. normal stretching pushes points farther from the interface). An
    optional conic term (n.H.n, from the fitted Hessian) is off by default.

    This is NOT reinitialisation: geometry is re-derived from psi each step; there
    is no PDE redistancing, no persistent triangulated surface, and (at CFL < 1)
    no extra communication beyond the one stencil-psi halo the base already does.

    Because a reprojected distance is NOT a stencil value, the Barth--Jespersen
    slope limiter and the corrector's stencil-range clamp -- both of which assume
    evaluateRaw reproduces psi_c inside the stencil range -- are incompatible with
    this model and are forced OFF in the constructor (limitSlope=false,
    clipToStencilBounds=false), so evaluate() passes evaluateRaw() through.

    U and grad U are looked up from the object registry inside update() (U == the
    new-time velocity at that point), so the Taylor departure-foot code in
    slAdvection is untouched. If U is absent (e.g. the unit-test harness) the
    reprojection still runs; only the strain rescaling is skipped.

    Selected with `reconstruction signedDistanceQuadraticWeightedLeastSquares`
    (short alias `sdQuadraticWLSQ`).

SourceFiles
    signedDistanceQuadraticWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef signedDistanceQuadraticWeightedLeastSquaresReconstruction_H
#define signedDistanceQuadraticWeightedLeastSquaresReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class signedDistanceQuadraticWeightedLeastSquaresReconstruction
:
    public slReconstruction
{
    // --- Quadratic value-fit machinery (identical to the uncached model) ---

    //- Active (non-empty) mesh directions (2 in 2D, 3 in 3D).
    labelList activeDirs_;

    //- Full quadratic coefficient count = 2*nd + nd*(nd-1)/2 (5 in 2D, 9 in 3D).
    label ncoeffFull_;

    //- Per-cell neighbour count (stencil size - 1) and its max.
    labelList nNbr_;
    label maxNbr_;

    //- Per-cell active coefficient count (ncoeffFull_ | nd linear | 0 constant).
    labelList ncoeff_;

    //- Fitted coefficients, one flat contiguous buffer [c*ncoeffFull_ + k].
    scalarField coeffsFlat_;

    //- One-off basis/fallback setup (mesh static within a run).
    bool built_;

    //- Optional Tikhonov ridge on the normal-equations diagonal (default 0).
    scalar ridgeEps_;

    // --- Signed-distance-preserving controls ---

    //- Apply the normal-strain rescaling d *= (1 + eps_nn*dt) (default true).
    Switch normalStrainRescale_;

    //- Newton iterations for the foot-point distance to the local quadric zero set,
    //  using the reconstruction's value and symbolic gradient g + H*delta each step
    //  (1 reproduces the linear P0/|grad P|; 2-3 use the Hessian). Default 3.
    label footpointNewtonIters_;

    //- Reproject only near the interface (the spec is band-only); far-field cells
    //  keep the plain value-copy. Reproject where |P0| < reprojBandRadii*stencilRadius.
    Switch reprojectBandOnly_;
    scalar reprojBandRadii_;

    //- Trust region on |grad P|: reproject only when the fitted gradient magnitude is
    //  in [gradLo, gradHi] (near the SDF value 1); otherwise the fit is untrustworthy
    //  for a P0/|grad P| division, so fall back to the raw value P0.
    scalar gradLo_;
    scalar gradHi_;

    // --- Per-step strain data, looked up from the registry in update() ---

    //- Cell velocity gradient grad U (via fvc::grad(U,"gradU")).
    tensorField gradU_;

    //- Current time step.
    scalar dt_;

    //- Whether U/grad U were available this step (else rescaling is skipped).
    bool haveStrain_;

    // --- Private members ---

    //- Constant-free quadratic basis at displacement d (identical to uncached).
    void basis(const vector& d, const label ncoeff, scalar* b) const;

    //- grad P at displacement d from the fitted coefficients (analytic).
    void gradFromCoeffs
    (
        const scalar* cf,
        const label nc,
        const vector& d,
        vector& g
    ) const;

    //- Active directions + ncoeffFull_ from mesh geometricD().
    void setupBasis();

    //- Set per-cell ncoeff_/nNbr_/maxNbr_ from the stencil sizes.
    void build();

public:

    TypeName("signedDistanceQuadraticWeightedLeastSquares");

    explicit signedDistanceQuadraticWeightedLeastSquaresReconstruction(const fvMesh& mesh);

    virtual ~signedDistanceQuadraticWeightedLeastSquaresReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
