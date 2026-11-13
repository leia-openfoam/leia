/*---------------------------------------------------------------------------*\
    leia :: davof :: davofNormalModels :: ringTilt
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::davofNormalModels::ringTilt

Description
    The ring area normal n_k (davofNormalModels::ringArea) rotated by a
    prescribed tilt theta_b = tiltC h / R, h = V_k^(1/3), R the radius of the
    exact sphere (key tiltRadius, default exactNormal/radius), a first-order
    angle error of known size. A test model of the normal-order study
    (leia davof-normal-order, 2026-11-13), not a reconstruction.

    The rotation is about the axis a = e x n_k / |e x n_k| (e the fixed
    direction tiltDir, default (0 0 1); where |e x n_k| < 1e-6 the axis
    e' x n_k with e' = (1 0 0) is taken), by the angle s theta_b:
        n' = n_k cos(s theta_b) + (a x n_k) sin(s theta_b),
    a unit-vector rotation (a . n_k = 0), so that n' turns away from e in
    the plane of e and n_k for s = +1. The magnitude |n_k| is kept.

    tiltMode biased: s = +1 in every cell, the same sense relative to the
    fixed direction e everywhere: a bias. tiltMode alternating: s = +1 or
    -1 by the parity of the cell's index sum floor(x/h + 1/2) + floor(y/h +
    1/2) + floor(z/h + 1/2) (a geometric parity, independent of the cell
    numbering and the decomposition): the same magnitude and order, no
    bias.

SourceFiles
    ringTiltNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringTilt_H
#define davofNormalModels_ringTilt_H

#include "ringAreaNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace davofNormalModels
{

class ringTilt
:
    public ringArea
{
    scalar tiltC_;
    scalar tiltR_;
    vector tiltDir_;
    bool alternating_;

public:

    //- Runtime type information
    TypeName("ringTilt");

    ringTilt
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );

    virtual ~ringTilt() = default;

    virtual void modelNormals(const scalar tState, vectorField& n);
};

} // End namespace davofNormalModels
} // End namespace Foam

#endif

// ************************************************************************* //
