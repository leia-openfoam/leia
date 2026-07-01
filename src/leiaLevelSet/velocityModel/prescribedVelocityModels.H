/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

Description
    Concrete prescribed (verification) velocity models selected at runtime by
    velocityModel::New: shear2D, deformation3D, shear3D, translation, rotation,
    vortex2D, periodic2D. Velocity formulae are unchanged from the original
    advectionVerification.H implementation.

SourceFiles
    prescribedVelocityModels.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedVelocityModels_H
#define prescribedVelocityModels_H

#include "velocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- 2D shear (single-vortex) flow
class shear2D
:
    public velocityModel
{
public:
    TypeName("shear2D");
    explicit shear2D(const fvMesh& mesh) : velocityModel(mesh) {}
    virtual ~shear2D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 3D deformation flow (LeVeque)
class deformation3D
:
    public velocityModel
{
public:
    TypeName("deformation3D");
    explicit deformation3D(const fvMesh& mesh) : velocityModel(mesh) {}
    virtual ~deformation3D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 3D shear flow (Liovic et al. 2006)
class shear3D
:
    public velocityModel
{
    scalar R_, Umax_, x0_, y0_;
public:
    TypeName("shear3D");
    explicit shear3D(const fvMesh& mesh);
    virtual ~shear3D() = default;
    virtual vector velocity(const vector& p) const;
};


//- Uniform translation
class translation
:
    public velocityModel
{
    vector velocity_;
public:
    TypeName("translation");
    explicit translation(const fvMesh& mesh);
    virtual ~translation() = default;
    virtual vector velocity(const vector& p) const;
};


//- Rigid-body rotation
class rotation
:
    public velocityModel
{
    vector point_;
    vector omega_;
public:
    TypeName("rotation");
    explicit rotation(const fvMesh& mesh);
    virtual ~rotation() = default;
    virtual vector velocity(const vector& p) const;
};


//- 2D single vortex
class vortex2D
:
    public velocityModel
{
    scalar v0_;
public:
    TypeName("vortex2D");
    explicit vortex2D(const fvMesh& mesh);
    virtual ~vortex2D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 2D periodic shear
class periodic2D
:
    public velocityModel
{
    scalar v0_;
    scalar c1_;
    scalar c2_;
public:
    TypeName("periodic2D");
    explicit periodic2D(const fvMesh& mesh);
    virtual ~periodic2D() = default;
    virtual vector velocity(const vector& p) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
