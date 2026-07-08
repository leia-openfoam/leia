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
    Foam::quadraticWLSQReconstruction

Description
    Production reconstruction: one weighted least-squares fit of the full
    quadratic polynomial per cell over its cell-point-cell stencil.

    Centred at x_c and fitting (psi_n - psi_c) eliminates the constant, so the
    cell value is reproduced exactly (evaluate(c, x_c) == psi_c). The basis over
    the non-empty directions is, in 3D,
        { dx, dy, dz, 1/2 dx^2, 1/2 dy^2, 1/2 dz^2, dxdy, dxdz, dydz }  (9)
    and in 2D the empty direction is dropped
        { dx, dy, 1/2 dx^2, 1/2 dy^2, dxdy }                            (5)
    (dropping it keeps the normal-equation system full rank on one-cell-thick
    meshes). Inverse-distance weights w_n = 1/|x_n - x_c|.

    The stencil geometry is fixed while the mesh is static, so the per-cell
    Moore-Penrose pseudo-inverse (MatrixTools::pinv, QR-based) that maps the
    weighted neighbour values to the polynomial coefficients is built ONCE and
    cached; the per-step cost is one halo exchange (stencil psiOld) plus one
    small matrix-vector product per cell. Cells with too few neighbours fall
    back per-cell to the linear fit (count reported).

    Performance: the cached map is stored per cell as ONE flat, row-major
    (ncoeff x Nneighbours) array with the inverse-distance weight folded into
    each column (MW[k,i] = pinv(k,i) w_i). The per-step coefficient update is
    then a contiguous pointer mat-vec coeffs_k = sum_i MW[k,i] (psi_i - psi_c),
    with no per-cell allocation and no operator() index arithmetic, so the inner
    loop auto-vectorises. hpctoolkit attributed ~28% of the solver to this loop
    before folding/flattening. See [[semi-lagrangian-solver]].

SourceFiles
    quadraticWLSQReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticWLSQReconstruction_H
#define quadraticWLSQReconstruction_H

#include "slReconstruction.H"
#include "scalarMatrices.H"

namespace Foam
{

class quadraticWLSQReconstruction
:
    public slReconstruction
{
protected:

    //- Active (non-empty) mesh directions, e.g. {0,1} on a 2D xy mesh.
    labelList activeDirs_;

    //- Full quadratic basis size over the active directions (5 in 2D, 9 in 3D).
    label ncoeffFull_;

    //- Cached per-cell map from neighbour differences to coefficients, stored
    //  flat row-major (ncoeff_[c] x nNbr_[c]) with the IDW weight folded in:
    //  MW_[c][k*nNbr + i] = pinv(k,i) * w_i. coeffs_k = sum_i MW[k,i](psi_i-psi_c).
    List<scalarList> MW_;

    //- Per-cell neighbour count (stencil size minus self); indexes MW_ rows.
    labelList nNbr_;

    //- Largest nNbr_ over all cells (size of the reused difference scratch).
    label maxNbr_;

    //- Per-cell fitted basis size (ncoeffFull_, or nd for linear fallback).
    labelList ncoeff_;

    //- Per-cell fitted coefficients (length ncoeff_[c]); refreshed each step.
    List<scalarList> coeffs_;

    bool built_;

    //- Quadratic (full) or linear basis of the displacement d over the active
    //  directions, written into the caller-provided buffer b (length ncoeff,
    //  no allocation).
    void basis(const vector& d, const label ncoeff, scalar* b) const;

    //- Build activeDirs_/ncoeffFull_ and the cached per-cell weighted maps
    //  from the collected stencil centres (called on the first update()).
    void build();

    //- Set activeDirs_/ncoeffFull_ from the mesh geometric directions (once).
    void setupBasis();

public:

    TypeName("quadraticWLSQ");

    explicit quadraticWLSQReconstruction(const fvMesh& mesh);

    virtual ~quadraticWLSQReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
