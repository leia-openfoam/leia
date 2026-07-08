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
    Foam::linearTaylorReconstruction

Description
    Degenerate (linear) Taylor reconstruction: the O(h) baseline.

        psi(x) ~= psi_c + g_c . (x - x_c)

    g_c is the cell-point-cell least-squares gradient of psiOld (fvSchemes
    `gradPsi pointCellsLeastSquares`, linear-exact on arbitrary stencils). With
    the semi-Lagrangian foot this reproduces first-order upwind, so it VALIDATES
    the pipeline and exhibits the order jump against the quadratic models.

SourceFiles
    linearTaylorReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef linearTaylorReconstruction_H
#define linearTaylorReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class linearTaylorReconstruction
:
    public slReconstruction
{
    //- Cell-centred LSQ gradient of psiOld (fvSchemes key "gradPsi").
    volVectorField gradPsi_;

public:

    TypeName("linearTaylor");

    explicit linearTaylorReconstruction(const fvMesh& mesh);

    virtual ~linearTaylorReconstruction() = default;

    virtual void update(const volScalarField& psiOld) override;

    virtual scalar evaluateRaw(const label c, const point& x) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
