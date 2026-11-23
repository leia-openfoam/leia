/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::quadraticSphere

Description
    The QUADRATIC level set of a sphere: with the signed distance d = |x - c| - R,

        psi = d (1 + d/(2 R)) = (|x - c|^2 - R^2)/(2 R).

    It has the zero set of the sphere and |grad psi| = 1 on it (|grad psi| = r/R
    elsewhere). It is a polynomial of degree 2 in x, also after a rigid motion.
    A quadratic reconstruction reproduces it exactly, on every mesh and at every
    refinement transition, so a translation or a rotation of it is an
    exact-solution gate of the whole semi-Lagrangian transport on a dynamic
    mesh: the error is the error of the departure point only (supervisor,
    2026-11-23, STATUS 11.44).

    The curvature and the phase indicator do not need the distance property.
    The metrics that measure | |grad psi| - 1 | in a band do: they are not a
    measure of error for this profile.

    Requires `implicitSurface { type implicitSphere; }`.
    Selected with `profile { type quadraticSphere; }` (token PROFILE).

SourceFiles
    quadraticSphere.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticSphere_H
#define quadraticSphere_H

#include "profile.H"

namespace Foam
{

class quadraticSphere
:
    public profile
{
    scalar radius_;

public:

    TypeName("quadraticSphere");

    quadraticSphere(const fvMesh& mesh, const implicitSurface& surface);

    virtual ~quadraticSphere() = default;

    virtual scalar value(const vector point) const;
};

} // End namespace Foam

#endif
