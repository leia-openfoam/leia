/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt 
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
    Foam::functionObjects::contactPoint

Group
    grpFieldFunctionObjects

Description
    - Works just for structured 2D meshes

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    contactPoint
    {
        type            contactPoint;
        libs            ( leiaFunctionObjects );
        writeControl    timeStep;
        field           gradPsiError;
    }

    \endverbatim



    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func contactPoint
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::contactPoint
    ...

SourceFiles
    contactPoint.C

\*---------------------------------------------------------------------------*/

#ifndef contactPoint_H
#define contactPoint_H

#include "fvCFD.H"
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class contactPoint Declaration
\*---------------------------------------------------------------------------*/

class contactPoint
:
    public fvMeshFunctionObject
{
    // Private Member Class
    class rootInterpolation
    {
        label size_;

        label label1_, label2_;
        
        scalar weight_;

        public:

            rootInterpolation() = default;

            void updateParameters(const scalarField& psi);

            scalar interpolate(const scalarField& vals);

            Pair<label> getLabels() const;

            scalar getWeight() const;

        private:

            Pair<label> findRootNeighbours(const scalarField& psi);

    };

    // Private Data Members

        word fieldName_;

        const volScalarField& field_;

        // wordList patchNames_;

        const word patchName_;

        const label patchID_;

        const fvPatch& patch_;

        const fvPatchField<scalar>& pfield_;

        uint component_;

        rootInterpolation rint_{};
        
        fileName filename_;
        OFstream outputFile_;




    // Private Member Functions

        bool is2D(const fvPatchField<scalar>& pfield);

        scalar calcContactPosition(const fvPatchField<scalar>& pfield);

        scalar calcContactAngle();

        scalar calcContactCurvature();

        scalar calcContactGradientNorm();

        template <template <class> class Field, class Type>
        tmp<Foam::scalarField> 
        getHalfField(const Field<Type>& pfield);

        // tmp<scalarField> getHalfField(const fvPatchField<scalar>& pfield);

public:

    //- Runtime type information
    TypeName("contactPoint");


    // Constructors

        //- Construct from Time and dictionary
        contactPoint
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        contactPoint(const contactPoint&) = delete;

        //- No copy assignment
        void operator=(const contactPoint&) = delete;


    //- Destructor
    virtual ~contactPoint() = default;


    // Public Member Functions

        virtual bool execute()
        {
            return true;
        }

        virtual bool write();
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
