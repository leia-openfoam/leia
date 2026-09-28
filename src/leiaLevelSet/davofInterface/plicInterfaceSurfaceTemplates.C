/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: plicInterfaceSurface templates
\*---------------------------------------------------------------------------*/

#include "plicInterfaceSurface.H"

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::plicInterfaceSurface::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords
) const
{
    return tmp<Field<Type>>::New(this->points().size(), Type(Zero));
}


// ************************************************************************* //
