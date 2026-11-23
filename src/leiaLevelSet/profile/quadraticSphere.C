/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "quadraticSphere.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(quadraticSphere, false);
addToRunTimeSelectionTable(profile, quadraticSphere, Mesh);

quadraticSphere::quadraticSphere
(
    const fvMesh& mesh,
    const implicitSurface& surface
)
:
    profile(mesh, surface),
    radius_(refCast<const implicitSphere>(surface).radius())
{}

scalar quadraticSphere::value(const vector point) const
{
    const scalar d = surface_.value(point);
    return d*(1.0 + 0.5*d/radius_);
}

} // End namespace Foam
