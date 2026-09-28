/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: sampledDavofInterface templates
\*---------------------------------------------------------------------------*/

#include "sampledDavofInterface.H"

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledDavofInterface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    updateGeometry();
    return sampledSurface::sampleOnFaces
    (
        sampler,
        surface().meshCells(),
        surface(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledDavofInterface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // The polygons are disjoint: no point interpolation, the cell value
    // is what the surface carries (sample with interpolate false).
    updateGeometry();
    return tmp<Field<Type>>::New(points().size(), Type(Zero));
}


// ************************************************************************* //
