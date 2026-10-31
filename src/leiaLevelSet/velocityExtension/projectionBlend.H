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
    Foam::projectionBlend

Description
    The blended projection extension of D. Bothe (discussion of 2026-10):

        P(x)       = x - dist(x) n(x),   n = grad psi/|grad psi|
        Uext(x)    = delta(x) U(x) + (1 - delta(x)) U(P(x))

    dist is the pointwise distance estimate of profile/rescaledDistance.H:
    psi/|grad psi| (distance first) or the quadratic root (distance second),
    so P(x) is one Newton step to the zero set (first) or its second-order
    version. delta blends from 0 near the zero set to 1 far from it:

        delta = 0                                   for |dist| <= w1 h
              = (1 - cos(pi s))/2, s = (|dist| - w1 h)/((w2 - w1) h)
              = 1                                   for |dist| >= w2 h

    with h the cell size. U(P) is a first-order Taylor step from the cell
    centre, U(x_c) + grad U(x_c).(P - x_c): a compact stencil, no point
    search, decomposition-independent (grad U through its coupled patches).
    The flux is phiExt = interpolate(Uext) & Sf on internal and coupled
    faces; physical boundary faces keep the prescribed flux (delta = 1 there
    for any wall farther than w2 h from the zero set).

    Dictionary (fvSolution levelSet.velocityExtension):
        type          projectionBlend;
        levelSet      psi;
        distance      second;          // first | second
        w1            1.5;             // cells: delta = 0 inside
        w2            4;               // cells: delta = 1 outside
        gradScheme    rescaleGradPsi;  // fvSchemes gradSchemes
        hessianScheme rescaleHessian;
        gradUScheme   gradUSdpls;
        gradFloor     1e-3;            // |grad psi| below: Uext = U

SourceFiles
    projectionBlend.C

\*---------------------------------------------------------------------------*/

#ifndef projectionBlend_H
#define projectionBlend_H

#include "velocityExtension.H"

namespace Foam
{

class projectionBlend
:
    public velocityExtension
{
    // Private Data

        const volScalarField& psi_;
        label order_;
        scalar w1_;
        scalar w2_;
        word gradScheme_;
        word hessianScheme_;
        word gradUScheme_;
        scalar gradFloor_;

        //- Cell size per cell (geometric dimensions only)
        scalarField h_;


public:

    TypeName("projectionBlend");

    projectionBlend(const fvMesh& mesh);

    virtual ~projectionBlend() = default;

    virtual void correct();
};

} // End namespace Foam

#endif

// ************************************************************************* //
