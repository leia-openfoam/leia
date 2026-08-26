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

defineTypeNameAndDebug(uniaxialStrain, 0);
addToRunTimeSelectionTable(velocityModel, uniaxialStrain, Mesh);

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


uniaxialStrain::uniaxialStrain(const fvMesh& mesh)
:
    velocityModel(mesh),
    strainRate_(velocityDict_.getOrDefault<scalar>("strainRate", 1)),
    stagnationPoint_
    (
        velocityDict_.getOrDefault<vector>("stagnationPoint", vector::zero)
    )
{
    // The base class reads `oscillation` with default "on", which multiplies U
    // and phi by cos(pi t / tau) at every step. Uniaxial stretching is a
    // monotone test whose exact solutions (see prescribedVelocityModels.H)
    // assume a time-constant a; the cosine factor would silently invalidate
    // all of them, and the failure would look like a discretization error
    // rather than a set-up error. Re-read the entry here with default "off" so
    // the safe behaviour is the default one. An explicit `oscillation on;`
    // still selects the oscillating variant, and no other model is affected:
    // this assignment lives in this constructor only.
    // NOT getOrDefault<Switch>(key, "off"). Switch(const char*) is EXPLICIT, so
    // a string literal cannot convert to Switch directly; it decays to a
    // non-null const char* and converts to bool as TRUE. That default therefore
    // means ON -- the exact opposite of what the comment above requires, and it
    // would silently oscillate a case that omitted the key, destroying the
    // closed-form solution this model exists to provide. A bool literal is
    // parsed as written.
    isOscillating_ = velocityDict_.getOrDefault<Switch>("oscillation", false);

    Info<< "    uniaxialStrain: v(x) = " << strainRate_
        << " (x - " << stagnationPoint_.x() << ") e_x, div(v) = "
        << strainRate_ << " 1/s (NOT solenoidal, intrinsic to 1D);"
        << " oscillation " << isOscillating_ << endl;
}

vector uniaxialStrain::velocity(const vector& p) const
{
    // v_y = v_z = 0 exactly: on an `empty`-patch 1D mesh nothing else may be
    // non-zero, and d(v_x)/dx = strainRate_ is then the whole velocity
    // gradient, so nhat . grad(v) . nhat = strainRate_ exactly.
    return vector(strainRate_*(p.x() - stagnationPoint_.x()), 0, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
