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
    Foam::anisotropicDiffusion

Description
    Velocity extension by anisotropic diffusion (Kummer): solves
    laplacian((n (x) n + eps I) Uext) = 0, diffusing Uext predominantly along
    the interface normal so that (n.grad)Uext -> 0. The interface seed layer is
    fixed as a Dirichlet constraint in the matrix (fvMatrix::setValues), and
    cells outside the band are pinned to the base velocity to confine the solve
    to the band. A small isotropic regularization eps keeps the system
    well-posed.

SourceFiles
    anisotropicDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicDiffusion_H
#define anisotropicDiffusion_H

#include "interfaceExtension.H"

namespace Foam
{

class anisotropicDiffusion
:
    public interfaceExtension
{
    //- Isotropic regularization of the normal-aligned diffusion tensor:
    //  the absolute value ('fixed') or the cell-size coefficient ('cellSize').
    scalar epsilon_;

    //- 'fixed' (legacy: eps = epsilon_, h-independent error floor) or
    //  'cellSize' (eps_c = epsilon_*cellSize_c, the min face-neighbour
    //  centre distance -> O(h)-consistent; NOT cbrt(V), which mis-scales
    //  on 2D meshes of arbitrary thickness).
    word epsilonModel_;

protected:

    virtual void extend() override;

public:

    TypeName("anisotropicDiffusion");

    anisotropicDiffusion(const fvMesh& mesh);

    virtual ~anisotropicDiffusion() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
