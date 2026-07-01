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

\*---------------------------------------------------------------------------*/

#include "prescribedVelocityModels.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(shear2D, 0);
addToRunTimeSelectionTable(velocityModel, shear2D, Mesh);

defineTypeNameAndDebug(deformation3D, 0);
addToRunTimeSelectionTable(velocityModel, deformation3D, Mesh);

defineTypeNameAndDebug(shear3D, 0);
addToRunTimeSelectionTable(velocityModel, shear3D, Mesh);

defineTypeNameAndDebug(translation, 0);
addToRunTimeSelectionTable(velocityModel, translation, Mesh);

defineTypeNameAndDebug(rotation, 0);
addToRunTimeSelectionTable(velocityModel, rotation, Mesh);

defineTypeNameAndDebug(vortex2D, 0);
addToRunTimeSelectionTable(velocityModel, vortex2D, Mesh);

defineTypeNameAndDebug(periodic2D, 0);
addToRunTimeSelectionTable(velocityModel, periodic2D, Mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vector shear2D::velocity(const vector& p) const
{
    const scalar& x = p[0];
    const scalar& y = p[1];
    return vector
    (
         Foam::sin(2*M_PI*y)*Foam::sqr(Foam::sin(M_PI*x)),
        -Foam::sin(2*M_PI*x)*Foam::sqr(Foam::sin(M_PI*y)),
         0
    );
}


vector deformation3D::velocity(const vector& p) const
{
    const scalar& x = p[0];
    const scalar& y = p[1];
    const scalar& z = p[2];
    return vector
    (
         2*sin(2*M_PI*y)*pow(sin(M_PI*x),2)*sin(2*M_PI*z),
        -sin(2*M_PI*x)*pow(sin(M_PI*y),2)*sin(2*M_PI*z),
        -sin(2*M_PI*x)*sin(2*M_PI*y)*pow(sin(M_PI*z),2)
    );
}


shear3D::shear3D(const fvMesh& mesh)
:
    velocityModel(mesh),
    R_(velocityDict_.getOrDefault<scalar>("R", 0.5)),
    Umax_(velocityDict_.getOrDefault<scalar>("Umax", 1.)),
    x0_(velocityDict_.getOrDefault<scalar>("x0", 0.5)),
    y0_(velocityDict_.getOrDefault<scalar>("y0", 0.5))
{}

vector shear3D::velocity(const vector& p) const
{
    const scalar& x = p[0];
    const scalar& y = p[1];
    const scalar r = Foam::sqrt(Foam::sqr(x - x0_) + Foam::sqr(y - y0_));
    return vector
    (
         Foam::sin(2*M_PI*y)*Foam::sqr(Foam::sin(M_PI*x)),
        -Foam::sin(2*M_PI*x)*Foam::sqr(Foam::sin(M_PI*y)),
         Umax_ * Foam::sqr(1 - (r / R_))
    );
}


translation::translation(const fvMesh& mesh)
:
    velocityModel(mesh),
    velocity_(velocityDict_.get<vector>("velocity"))
{}

vector translation::velocity(const vector& /*p*/) const
{
    return velocity_;
}


rotation::rotation(const fvMesh& mesh)
:
    velocityModel(mesh),
    point_(velocityDict_.get<vector>("point")),
    omega_(velocityDict_.get<vector>("omega"))
{}

vector rotation::velocity(const vector& p) const
{
    return omega_ ^ (p - point_);
}


vortex2D::vortex2D(const fvMesh& mesh)
:
    velocityModel(mesh),
    v0_(velocityDict_.getOrDefault<scalar>("v0", -0.2))
{}

vector vortex2D::velocity(const vector& p) const
{
    return v0_ * Vector<scalar>
    (
        -sin(M_PI*p[0])*cos(M_PI*p[1]),
         cos(M_PI*p[0])*sin(M_PI*p[1]),
         0.0
    );
}


periodic2D::periodic2D(const fvMesh& mesh)
:
    velocityModel(mesh),
    v0_(velocityDict_.getOrDefault<scalar>("v0", -0.2)),
    c1_(0.1),
    c2_(-2)
{}

vector periodic2D::velocity(const vector& p) const
{
    return Vector<scalar>(v0_ + c1_*p[0] + c2_*p[1], -c1_*p[1], 0.0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
