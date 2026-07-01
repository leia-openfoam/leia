/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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
    Foam::pseudoTime

Description
    Velocity extension by the Adalsteinsson-Sethian normal-advection PDE
    sign(psi) (n.grad)Uext = 0, solved IMPLICITLY: an upwind convection matrix
    along w = sign(psi) n carries the interface velocity outward, with the
    interface seed layer (and the band exterior) imposed as Dirichlet
    constraints via fvMatrix::setValues. `nIterations` outer solves are used
    (the linear system converges quickly).

SourceFiles
    pseudoTime.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoTime_H
#define pseudoTime_H

#include "interfaceExtension.H"

namespace Foam
{

class pseudoTime
:
    public interfaceExtension
{
    //- Number of implicit pseudo-time steps.
    label nIterations_;

    //- Pseudo-time step as a fraction of the local cell size (CFL-like).
    scalar deltaTau_;

protected:

    virtual void extend() override;

public:

    TypeName("pseudoTime");

    pseudoTime(const fvMesh& mesh);

    virtual ~pseudoTime() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
