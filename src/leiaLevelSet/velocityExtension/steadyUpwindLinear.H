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
    Foam::steadyUpwindLinear

Description
    Aslam's LINEAR extrapolation cascade (Aslam JCP 193 (2004); implicit /
    fast-sweeping form: Aslam, Luo & Zhao SISC 36 (2014); avoid the FMM
    variant of McCaslin et al. JCP 274 (2014) on unstructured meshes) built on
    the steadyUpwind operator:

      1.  V = n . grad(U)            (least-squares gradient, per cell)
      2.  extend V:                  sign(psi) (n . grad) V = 0
      3.  extend U with source:      sign(psi) (n . grad) Uext = sign(psi) V

    The extended field is then LINEAR along each normal ray (matches the
    normal derivative of the physical velocity at Sigma) instead of constant:
    O(h^2) extension quality by construction, the PDE counterpart of the
    closest-point model.

SourceFiles
    steadyUpwindLinear.C

\*---------------------------------------------------------------------------*/

#ifndef steadyUpwindLinear_H
#define steadyUpwindLinear_H

#include "steadyUpwind.H"

namespace Foam
{

class steadyUpwindLinear
:
    public steadyUpwind
{
protected:

    virtual void extend() override;

public:

    TypeName("steadyUpwindLinear");

    steadyUpwindLinear(const fvMesh& mesh);

    virtual ~steadyUpwindLinear() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
