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
    Foam::linearWeightedLeastSquaresReconstruction

Description
    Per-cell IDW-weighted least-squares fit of a LINEAR polynomial directly to
    the stencil psi VALUES, evaluated at the departure foot:

        psi(x) ~= psi_c + sum_a g_a * d_a,   d = x - x_c,
        g = argmin_g  sum_i w_i^2 ( psi_c + g.d_i - psi_i )^2,   w_i = 1/|d_i|.

    The linear-basis sibling of uncachedQuadraticWeightedLeastSquares (UQWLSR):
    the SAME per-cell weighted normal-equations assembly (M = A^T W^2 A,
    g = A^T W^2 ds) + in-place Cholesky (or householderQR), differing ONLY in
    the basis size: nd coefficients (2 in 2D, 3 in 3D) instead of the quadratic
    2*nd + nd(nd-1)/2. The slope is a byproduct of the VALUE fit, never an
    injected gradient field.

    Member of the WeightedLeastSquares (value-fitting) family. Unlike the Taylor
    family (linearTaylor, quadraticTaylor), which expand about the cell centre
    using a DIFFERENTIATED field (fvc::grad of psi), this fits VALUES over the
    stencil, so the reconstruction is bounded within the stencil's data range
    and never differentiates the field. It therefore does NOT feed the
    high-wavenumber steepening instability of the Taylor family (that
    instability is driven by fvc::grad amplifying 2h-mode noise; a value fit has
    no such amplification), and is STABLE BY CONSTRUCTION -- the same reason
    UQWLSR is stable, at linear (first) order and ~3-9x lower per-cell cost
    (2x2/3x3 normal equations vs 5x5/9x9).

    NO CURVATURE: a linear fit has an identically-zero Hessian, so this model
    does not override the base meanCurvature* methods (they stay NotImplemented,
    consistent with linearTaylor). Two-phase surface-tension coupling must take
    curvature from an FVM path (surfaceTensionForce divGradPsiSnGradAlpha with
    curvatureExtension fvm), never the reconstruction's Hessian.

    Selected with `reconstruction linearWeightedLeastSquares`. Composes with the
    default `correction direct` (the fit is exact per step; no iteration).

SourceFiles
    linearWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef linearWeightedLeastSquaresReconstruction_H
#define linearWeightedLeastSquaresReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class linearWeightedLeastSquaresReconstruction
:
    public slReconstruction
{
protected:
    // Protected (not private) so the signedDistanceLinearWeightedLeastSquares
    // subclass can reuse the fit and re-normalise the evaluation.

    //- Active (non-empty) mesh directions (2 in 2D, 3 in 3D).
    labelList activeDirs_;

    //- Linear coefficient count = nd (2 in 2D, 3 in 3D).
    label ncoeffFull_;

    //- Per-cell neighbour count (stencil size - 1) and its max (scratch sizing).
    labelList nNbr_;
    label maxNbr_;

    //- Per-cell active coefficient count (ncoeffFull_ = nd, or 0 constant).
    labelList ncoeff_;

    //- Fitted coefficients, ONE flat contiguous buffer [c*ncoeffFull_ + k],
    //  refreshed every step (no per-cell List handle / heap alloc).
    scalarField coeffsFlat_;

    //- One-off basis/fallback setup (mesh static within a run).
    bool built_;

    //- Optional Tikhonov ridge on the normal-equations diagonal (default 0):
    //  M(k,k) += ridgeEps * trace(M)/ncoeff. Conditions distorted-mesh solves.
    scalar ridgeEps_;

    //- Least-squares solver: "normalEquations" (default; A^T W^2 A + Cholesky)
    //  or "householderQR" (QR of the weighted design rows W A directly).
    const word fit_;

    //- Constant-free linear basis at displacement d (nd terms).
    void basis(const vector& d, const label ncoeff, scalar* b) const;

    //- Active directions + ncoeffFull_ = nd from mesh geometricD().
    void setupBasis();

    //- Set per-cell ncoeff_/nNbr_/maxNbr_ from the stencil sizes.
    void build();

public:

    TypeName("linearWeightedLeastSquares");

    explicit linearWeightedLeastSquaresReconstruction(const fvMesh& mesh);

    virtual ~linearWeightedLeastSquaresReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;

    //- No meanCurvature* overrides: a linear fit has zero Hessian, so curvature
    //  is inherited as NotImplemented (see class Description).
};

} // End namespace Foam

#endif

// ************************************************************************* //
