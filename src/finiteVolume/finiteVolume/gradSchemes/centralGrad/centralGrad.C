/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Julian Reitzel
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

#include "centralGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"
// #include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::centralGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    tmp<GradFieldType> tgGrad
    (
        new GradFieldType
        (
            IOobject
            (
                "grad(" + vsf.name() + ")",
                fileName(),
                vsf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vsf.mesh(),
            dimensioned<GradType>(vsf.dimensions()/dimLength, Foam::zero())
        )
    );
    GradFieldType& gGrad = tgGrad.ref();

    labelVector sizes = map_.sizes();


    if (sizes.x() != 1) // non-empty direction
    {
        forAll(gGrad, ID)
            {
                const labelVector ijk = map_(ID);
                const label i = ijk[0];
                const label j = ijk[1];
                const label k = ijk[2];

                if (i == 0) // left boundary
                {
                    gGrad[ID].x() = 
                        ( 
                            -   vsf[map_(i+2, j, k)] 
                            + 4*vsf[map_(i+1, j, k)]
                            - 3*vsf[map_(i  , j, k)]
                        )/(2*map_.h());
                }
                else if (i == sizes.x() - 1) // right boundary
                {
                    gGrad[ID].x() = 
                        ( 
                                vsf[map_(i-2, j, k)] 
                            - 4*vsf[map_(i-1, j, k)]
                            + 3*vsf[map_(i  , j, k)]
                        )/(2*map_.h());
                }
                else // (0 < i && i < sizes.x())
                {
                    gGrad[ID].x() = 
                        ( 
                                vsf[map_(i+1, j, k)] 
                            -   vsf[map_(i-1, j, k)]
                        )/(2*map_.h());
                }
            }
    }

    if (sizes.y() != 1) // non-empty direction
    {
        forAll(gGrad, ID)
            {
                const labelVector ijk = map_(ID);
                const label i = ijk[0];
                const label j = ijk[1];
                const label k = ijk[2];

                if (j == 0) // left boundary
                {
                    gGrad[ID].y() = 
                        ( 
                            -   vsf[map_(i, j+2, k)] 
                            + 4*vsf[map_(i, j+1, k)]
                            - 3*vsf[map_(i, j  , k)]
                        )/(2*map_.h());
                }
                else if (j == sizes.y() - 1) // right boundary
                {
                    gGrad[ID].y() = 
                        ( 
                                vsf[map_(i, j-2, k)] 
                            - 4*vsf[map_(i, j-1, k)]
                            + 3*vsf[map_(i, j  , k)]
                        )/(2*map_.h());
                }
                else // (0 < j && j < sizes.y())
                {
                    gGrad[ID].y() = 
                        ( 
                                vsf[map_(i, j+1, k)] 
                            -   vsf[map_(i, j-1, k)]
                        )/(2*map_.h());
                }
            }
    }

    if (sizes.z() != 1) // non-empty direction
    {
        forAll(gGrad, ID)
            {
                const labelVector ijk = map_(ID);
                const label i = ijk[0];
                const label j = ijk[1];
                const label k = ijk[2];

                if (k == 0) // left boundary
                {
                    gGrad[ID].z() = 
                        ( 
                            -   vsf[map_(i, j, k+2)] 
                            + 4*vsf[map_(i, j, k+1)]
                            - 3*vsf[map_(i, j, k  )]
                        )/(2*map_.h());
                }
                else if (k == sizes.z() - 1) // right boundary
                {
                    gGrad[ID].z() = 
                        ( 
                                vsf[map_(i, j, k-2)] 
                            - 4*vsf[map_(i, j, k-1)]
                            + 3*vsf[map_(i, j, k  )]
                        )/(2*map_.h());
                }
                else // (0 < i && i < sizes.x())
                {
                    gGrad[ID].z() = 
                        ( 
                                vsf[map_(i, j, k+1)] 
                            -   vsf[map_(i, j, k-1)]
                        )/(2*map_.h());
                }
            }
    }

    Field<GradType>& igGrad = gGrad;
    igGrad /= vsf.mesh().V();

     // TODO Implement coupled boundary case

    return tgGrad;
}




// ************************************************************************* //
