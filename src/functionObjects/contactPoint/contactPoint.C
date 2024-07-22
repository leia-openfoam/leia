/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Julian Reitzel, TU Darmstadt 
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


#include "contactPoint.H"
#include "noBvGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "FieldOps.H"
#include "interpolationCellPoint.H"
#include <vector>

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * * //

template <class Type>
static Vector<double> findRootPosition(const fvPatchField<Type>& patchfield, int whichRoot = -1)
{
    const vectorField& Cf = patchfield.patch().Cf(); // face centres

    std::vector<Vector<double>> roots;
    for (int i=1; i < patchfield.size(); ++i)
    {
        if (sign(patchfield[i-1]) != sign(patchfield[i]))
        {
            Type y1 = patchfield[i-1];
            Type y2 = patchfield[i];
            Vector<double> x1 = Cf[i-1];
            Vector<double> x2 = Cf[i];
            Vector<double> rootPos = -y1*(x2 - x1)/(y2 - y1) + x1;
            roots.push_back(rootPos);
        }
    }
    if (sign(whichRoot) < 0)
    {
        whichRoot = roots.size() + whichRoot;
    }
    return roots.at(whichRoot);
}


template <class Type>
static Vector<Type> mulTenVec(const Tensor<Type>& ten, const Vector<Type>& vec)
{
    return Vector<Type>(
        ten.x() & vec,
        ten.y() & vec,
        ten.z() & vec
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(contactPoint, 0);
    addToRunTimeSelectionTable(functionObject, contactPoint, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::contactPoint::contactPoint
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject (name, runTime, dict),
    fieldName_(dict.get<word>("field")),
    field_(mesh_.lookupObject<volScalarField>(fieldName_)),
    patchName_(dict.get<word>("patch")),
    patchID_(mesh_.boundary().findPatchID(patchName_)),
    patch_(mesh_.boundary()[patchName_]),
    pfield_(patch_.patchField(field_)),
    component_(dict.getOrDefault<uint>("component",0)),
    filename_("contactPoint.csv"),
    outputFile_(filename_,  IOstreamOption(), IOstreamOption::appendType::APPEND)

{
    if ( !fileSize(filename_) && Pstream::myProcNo() == 0)
    {
            // CSV Header 
            outputFile_ << "TIME,"
                << "CONTACT_POINT,"
                << "CONTACT_ANGLE,"
                << "CONTACT_CURVATURE,"
                << "CONTACT_MAGGRAD"
                << endl;
    }
    write();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// bool Foam::functionObjects::contactPoint::read(const dictionary& dict)
// {
//     if (fvMeshFunctionObject::read(dict))
//     {
//         if (fieldName_.empty() || dict.found("field"))
//         {
//             dict.readEntry("field", fieldName_);
//         }
         
//         if (patchNames_.empty() || dict.found("field"))
//         {
//             dict.readEntry("patches", patchNames_);
//         }
 
//         return true;
//     }
 
//     return false;
// }

bool Foam::functionObjects::contactPoint::is2D(const fvPatchField<scalar>& pfield)
{
    return mesh_.nSolutionD() == 2;
}

template <template <class> class Field, class Type>
Foam::tmp<Foam::scalarField> 
Foam::functionObjects::contactPoint::getHalfField(const Field<Type>& pfield)

{
    // labelRange half(0, pfield.size()/2); // first half
    labelRange half(pfield.size()/2, pfield.size()); // second half

    tmp<scalarField> thpfield(new scalarField(pfield.slice(half)));
    return thpfield;
}

//- Calculates the position of the interface contact position with a wall patch
//
//  Assumptions:
//  - Mesh is 2D
//  - Case is symmetric -> just one half of the patch is used
scalar Foam::functionObjects::contactPoint::calcContactPosition(const fvPatchField<scalar>& pfield)
{
    if (!is2D(pfield))
    {
            FatalErrorInFunction
            << "Mesh is not 2D" << nl
            << exit(FatalError);
    }

    // tmp<scalarField> tvals = getHalfField(pfield);
    // const scalarField vals = tvals();
    const scalarField vals = getHalfField(pfield)();

    const scalarField pos = getHalfField(pfield.patch().patch().faceCentres().component(component_)());

    rint_.updateParameters(vals);
    return rint_.interpolate(pos);
}

//- Calculates the contact angle of the interface with a wall patch
//
//  contactAngle = 180° - arccos((\grad{psi} \cdot \Sf)/(|\grad{psi}| |\Sf|))
scalar Foam::functionObjects::contactPoint::calcContactAngle()
{
    const volVectorField gradPsi = fv::noBvGrad<scalar>(mesh_).grad(field_);
    volVectorField normal = gradPsi/mag(gradPsi);
    // normal.boundaryFieldRef().set(patchID_, fvPatchField<vector>::New("oneSidedGradient", patch_, normal));
    // normal.boundaryFieldRef().evaluate(); 

    const vectorField normal_p = patch_.patchInternalField(normal); // patchInternal
    // const vectorField& normal_p = normal.boundaryField()[patchID_]; // patch
    const vectorField& Sf_p = mesh_.Sf().boundaryField()[patchID_];

    const scalarField cosAlpha_p = (normal_p & Sf_p)/(mag(normal_p)*mag(Sf_p));
    
    return 180 - 180/constant::mathematical::pi*acos(rint_.interpolate(getHalfField(cosAlpha_p)));
}

scalar Foam::functionObjects::contactPoint::calcContactCurvature()
{
    const volVectorField gradPsi = fv::noBvGrad<scalar>(mesh_).grad(field_);
    volVectorField normal = gradPsi/mag(gradPsi);
    volVectorField tau1 = volVectorField("tau1", normal);
    volVectorField tau2 = volVectorField("tau2", normal);
    forAll(normal, ID)
    {
        tau1[ID] = normal[ID] ^ vector(0,0,1);
        tau2[ID] = tau1[ID] ^ normal[ID];
    }
    volTensorField gradNormal = fv::noBvGrad<vector>(mesh_).grad(normal);

    gradNormal.boundaryFieldRef().set(patchID_, fvPatchField<tensor>::New("zeroGradient", patch_, gradNormal));
    gradNormal.boundaryFieldRef().evaluate();

    volScalarField kappa = -tr(fv::noBvGrad<vector>(mesh_).grad(normal));
    forAll(kappa, ID)
    {
        kappa[ID] = (-1*mulTenVec(gradNormal[ID], tau1[ID]) & tau1[ID]) + (-1*mulTenVec(gradNormal[ID], tau2[ID]) & tau2[ID]);
    }
    kappa.boundaryFieldRef().set(patchID_, fvPatchField<scalar>::New("zeroGradient", patch_, kappa));
    kappa.boundaryFieldRef().evaluate();

    const Vector<scalar> root = findRootPosition(pfield_);
    const label rootCell = mesh_.findCell(root);

    tensor ten = interpolationCellPoint<tensor>(gradNormal).interpolate(root, rootCell);
    // tensor ten = gradNormal[762];
    std::string val_str = 
        "dnx_dy == dny_dx: " 
        + std::to_string(ten.xy()) + ", "
        + std::to_string(ten.yx()) + ";  "
        + "dnx_dz == dnz_dx: "
        + std::to_string(ten.xz()) + ", "
        + std::to_string(ten.zx()) + ";  "
        + "dnz_dy == dny_dz: "
        + std::to_string(ten.zy()) + ", "
        + std::to_string(ten.yz()) + ";  ";

    Info << "symmetry gradNormal: " << val_str << "\n";
    std::cout<< "\nsymmetry-error:\t" << ten.xy() - ten.yx() << "\n";
    return interpolationCellPoint<scalar>(kappa).interpolate(root, rootCell);

}

// scalar Foam::functionObjects::contactPoint::calcContactCurvature()
// {
//     const volVectorField gradPsi = fv::noBvGrad<scalar>(mesh_).grad(field_);
//     volVectorField normal = gradPsi/mag(gradPsi);
//     // normal.boundaryFieldRef().set(patchID_, fvPatchField<vector>::New("oneSidedGradient", patch_, normal));
//     // normal.boundaryFieldRef().evaluate(); 


//     // volScalarField kappa = -fvc::div(normal);
//     volScalarField kappa = -tr(fv::noBvGrad<vector>(mesh_).grad(normal));

//     // kappa.boundaryFieldRef().set(patchID_, fvPatchField<scalar>::New("oneSidedGradient", patch_, kappa));
//     // kappa.boundaryFieldRef().evaluate(); 

//     kappa.boundaryFieldRef().set(patchID_, fvPatchField<scalar>::New("zeroGradient", patch_, kappa));
//     kappa.boundaryFieldRef().evaluate(); 

//     const scalarField kappa_ph = getHalfField(kappa.boundaryField()[patchID_]); // patch
//     // const scalarField kappa_ph = getHalfField(patch_.patchInternalField(kappa)()); // patchInternal
//     return rint_.interpolate(kappa_ph);
// }

scalar Foam::functionObjects::contactPoint::calcContactGradientNorm()
{
    const volScalarField magGradPsi = mag(fvc::grad(field_));
    const scalarField magGradPsi_ph = getHalfField(patch_.patchInternalField(magGradPsi)());
    return rint_.interpolate(magGradPsi_ph);
}

bool Foam::functionObjects::contactPoint::write()
{

    scalar cpoint = calcContactPosition(pfield_);
    scalar cangle = calcContactAngle();
    scalar ccurvature = calcContactCurvature();
    scalar cmaggradient = calcContactGradientNorm();

    if (Pstream::myProcNo() == 0)
    {
        outputFile_ << mesh_.time().timeOutputValue() << ","
            << cpoint << ","
            << cangle << ","
            << ccurvature << ","
            << cmaggradient
            << endl;
    }
    return true;
}



// class rootInterpolation

// Foam::functionObjects::contactPoint::rootInterpolation::
// rootInterpolation()
// {}

scalar 
Foam::functionObjects::contactPoint::rootInterpolation::
interpolate(const scalarField& y)
{
    if (size_ != y.size())
    {
        FatalErrorInFunction
        << "Can not interpolate. scalarField has a different size." << nl
        << "size_: " << size_ << "\ty.size(): " << y.size()
        << exit(FatalError);
    }
    // y0 = w*y2 + (1 - w)y1
    return weight_*y[label2_] + (1 - weight_)*y[label1_];
}

Pair<label> Foam::functionObjects::contactPoint::rootInterpolation::getLabels() const
{
    return Pair<label>(label1_, label2_);
}

scalar Foam::functionObjects::contactPoint::rootInterpolation::getWeight() const
{
    return weight_;
}

//- Estimates the discrete root by returning its neighbours
//  Assumptions:
//  - field is 2D
//  - field is monotone -> just one root
Pair<label> 
Foam::functionObjects::contactPoint::rootInterpolation::
findRootNeighbours(const scalarField& psi)
{
    label label1 = 0, label2 = 0;
    Field<label> labels = labelRange(psi.size()).labels();
    Tuple2<scalar, label> minValLabel = Foam::FieldOps::findMinData(mag(psi)(), labels);

    label1 = minValLabel.second();
    scalar sign1 = sign(psi[label1]);

    if 
    (
            label1 - 1 >= 0 
        ||  sign1 != sign(psi[label1 - 1])
    )
    {
        label2 = label1 - 1;
    }
    else if
    (
            label1 + 1 >= 0 
        ||  sign1 != sign(psi[label1 + 1])
    )
    {
        label2 = label1 + 1;
    }
    else
    {
        FatalErrorInFunction
        << "No sign change in field" << nl
        << exit(FatalError);
    }
    return Pair<label>(label1, label2);
}

void 
Foam::functionObjects::contactPoint::rootInterpolation::
updateParameters(const scalarField& psi)
{
    if (size_ == 0)
    {
        size_ = psi.size();
    }
    Pair<label> labels = findRootNeighbours(psi);
    label1_ = labels.first();
    label2_ = labels.second();
    
    scalar psi1 = psi[labels.first()];
    scalar psi2 = psi[labels.second()];

    // weight_ w:  w*val2 + (1 - w)*val1 == 0
    weight_ = psi1/(psi1 - psi2);
}
// ************************************************************************* //
