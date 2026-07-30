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
    Foam::quadraticTaylorReconstruction

Description
    Second-order (quadratic) Taylor reconstruction with derivatives from nested
    linear least squares:

        psi(x) ~= psi_c + g_c . d + 1/2 d^T H_c d,     d = x - x_c

    g_c = LSQ gradient of psiOld (fvSchemes key `gradPsi`);
    H_c = LSQ gradient of g_c, symmetrised (fvSchemes key `gradGradPsi`).
    Both passes are cell-point-cell least squares (linear-exact on arbitrary
    stencils), so the nested construction is exact for quadratic psi => global
    O(h^2). Nested Green-Gauss would NOT be consistent on skewed meshes and is
    deliberately avoided (the LSQ scheme is selected by name).

    Member of the TAYLOR family (linearTaylor, quadraticTaylor): it expands
    about the cell centre using DIFFERENTIATED FIELDS (fvc::grad of psi, then
    grad-of-grad). Differentiating the field amplifies high-wavenumber error
    (O(k) for the gradient, O(k^2) for the grad-of-grad Hessian), so this
    reconstruction is prone to a signed-distance steepening instability on
    under-resolved fronts, unlike the value-fitting WeightedLeastSquares family
    (quadraticWeightedLeastSquares / linearWeightedLeastSquares), which fit a
    polynomial to stencil VALUES and never differentiate the field.

    Selected with `reconstruction quadraticTaylor;`. (Formerly named `nestedLSQ`,
    for the nested LSQ gradient passes; renamed for the {linear,quadratic} x
    {Taylor,WLSQ} taxonomy. No back-compat alias -- update old dicts to the new name.)

SourceFiles
    quadraticTaylorReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticTaylorReconstruction_H
#define quadraticTaylorReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class quadraticTaylorReconstruction
:
    public slReconstruction
{
    //- Cell-centred LSQ gradient of psiOld (fvSchemes key "gradPsi").
    volVectorField gradPsi_;

    //- Symmetrised cell-centred LSQ Hessian of psiOld (key "gradGradPsi").
    volTensorField hessPsi_;

public:

    TypeName("quadraticTaylor");

    explicit quadraticTaylorReconstruction(const fvMesh& mesh);

    virtual ~quadraticTaylorReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
