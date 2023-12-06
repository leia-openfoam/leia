/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::oneSidedGradientFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition applies a zero-gradient condition from the patch
    internal field onto the patch faces when \c evaluated but may also be
    assigned.  \c snGrad returns the patch gradient evaluated from the current
    internal and patch field values rather than returning zero.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            extrapolatedCalculated;
    }
    \endverbatim

SourceFiles
    oneSidedGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oneSidedGradientFvPatchField_H
#define oneSidedGradientFvPatchField_H

#include "calculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class oneSidedGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oneSidedGradientFvPatchField
:
    public calculatedFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("oneSidedGradient");


    // Constructors

        //- Construct from patch and internal field
        oneSidedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oneSidedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        oneSidedGradientFvPatchField
        (
            const oneSidedGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oneSidedGradientFvPatchField
        (
            const oneSidedGradientFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new oneSidedGradientFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oneSidedGradientFvPatchField
        (
            const oneSidedGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new oneSidedGradientFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueInternalCoeffs
        (
            const tmp<Field<scalar>>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueBoundaryCoeffs
        (
            const tmp<Field<scalar>>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


    private:

    // Member functions

        tmp<Field<Type>> values();

        bool isCartesianMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "oneSidedGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
