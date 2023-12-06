/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "oneSidedGradientFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oneSidedGradientFvPatchField<Type>::
oneSidedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::oneSidedGradientFvPatchField<Type>::
oneSidedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchField<Type>(p, iF, dict, false)
{
    evaluate();
}


template<class Type>
Foam::oneSidedGradientFvPatchField<Type>::
oneSidedGradientFvPatchField
(
    const oneSidedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::oneSidedGradientFvPatchField<Type>::
oneSidedGradientFvPatchField
(
    const oneSidedGradientFvPatchField<Type>& ptf
)
:
    calculatedFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::oneSidedGradientFvPatchField<Type>::
oneSidedGradientFvPatchField
(
    const oneSidedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oneSidedGradientFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    calculatedFvPatchField<Type>::operator==(this->values());
    calculatedFvPatchField<Type>::evaluate();
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::oneSidedGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}
 
 
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::oneSidedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::oneSidedGradientFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}
 
 
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::oneSidedGradientFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::oneSidedGradientFvPatchField<Type>::values()
{
    // Global
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    // const surfaceVectorField deltas = mesh.delta()();
    // const surfaceScalarField magdeltas = mag(deltas);
    const scalar dx = 2*this->patch().deltaCoeffs()[0];


    // Patch
    tmp<Field<Type>> tfield(new Field<Type>(this->size()));
    Field<Type>& field = tfield.ref();
    // tmp<Field<Type>> tfieldPI = this->patchInternalField();
    // const Field<Type>& fieldPI = tfieldPI();

    // const labelUList& cellPILabels = this->faceCells();

    // const SubList<face> facesP = faces.slice(this->patch().start(), this->patch.size());

    // if 
    // (
    //         cellPILabels.size() != field.size() 
    //     || facesP.size() != field.size())
    // {
    //     FatalErrorInFunction
    //     << "NOT EQUAL SIZE" << nl
    //     << exit(FatalError);
    // }

    // 
    // cell2
    // _______ face2
    // cell1
    // _______ face1
    // cellPI: patchInternalCell
    //________ faceP
    forAll(field, ID)
    {
        // const face faceP = facesP[ID];
        const label facePLabel = this->patch().start() + ID;
        const label cellPILAbel = own[facePLabel];
        // const cell cellPI = cells[cellPILabels[ID]];

        const label face1Label = cells[cellPILAbel].opposingFaceLabel(facePLabel, faces);
        const label cell1Label = nei[face1Label];
        const label face2Label = cells[cell1Label].opposingFaceLabel(face1Label, faces);
        const label cell2Label = nei[face2Label];

        // const scalar f0 = fieldPI[ID];
        const Type f0 = this->internalField()[cellPILAbel];
        const Type f1 = this->internalField()[cell1Label];
        const Type f2 = this->internalField()[cell2Label];
        // const scalar dx = 1/mag(mesh.delta()()[face1Label]);
        // const scalar dx = mag(deltas[face1Label]);
        // const scalar dx = 1;
        const Type df = (-3*f0 + 4*f1 - f2)/(2*dx);
        field[ID] = f0 - df*dx/2;
    }

    return tfield;
}

// template<class Type>
// void Foam::oneSidedGradientFvPatchField<Type>::updateCoeffs()
// {


//     updated_ = true;
// }

template<class Type>
bool Foam::oneSidedGradientFvPatchField<Type>::isCartesianMesh()
{
    // const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    // const surfaceVectorField& Sf = mesh.Sf();
    // // const surfaceVectorField nf = mesh.Sf()/mesh.magSf();

    // forAll(Sf, ID)
    // {
    //     const vector n = Sf[ID]/mag(Sf[ID]);
    //     if (
    //             n != vector(1,0,0)
    //         ||  n != vector(0,1,0)
    //         ||  n != vector(0,0,1)
    //     )
    //     {
    //         return false;
    //         // FatalErrorInFunction
    //         // << "Mesh is not cartesian." << nl
    //         // << exit(FatalError);
    //     }
    // }
    return true;
}

// ************************************************************************* //
