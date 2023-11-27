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

\*---------------------------------------------------------------------------*/

#include "sdplsBeta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
 
namespace Foam
{

static scalar etaCutoff(const scalar val, const scalar b)
{
    if (val <= b)
    {
        return 1;
    }
    else if (val > b)
    {
        return 0;
    }
    else
    {
        return 0.5*cos((val-2*b)*constant::mathematical::pi) + 0.5;
    }
} 

static scalar psiCutoff(const scalar psi)
{
    scalar ib = 0.5;    // inner bound
    scalar ob = 1;      // outer bound
    if (mag(psi) < ib)
        return 1;
    else if (mag(psi) < ob)
        return 0.5*cos((psi - ib)*(constant::mathematical::pi)/(ob-ib)) + 0.5;
    else
        return 0;
} 


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sdplsBeta, false);
    addToRunTimeSelectionTable(sdplsSource, sdplsBeta, Dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::sdplsBeta::sdplsBeta(const dictionary& dict, const fvMesh& mesh)
    :
        sdplsSource(dict, mesh),
        alpha_(dict.getOrDefault<scalar>("alpha", -1)),
        beta_(dict.get<scalar>("beta")),
        etaCutoff_(dict.getOrDefault<bool>("etaCutoff", false)),
        psiCutoff_(dict.getOrDefault<bool>("psiCutoff", false))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<volScalarField> 
Foam::sdplsBeta::
nonLinearPart
    (
        const volScalarField& R, 
        const volScalarField& psi, 
        const volVectorField& U
    ) const
{
    const fvMesh& mesh = psi.mesh();

    tmp<volScalarField> tfield
    (
        new volScalarField
        (
            IOobject
            (
                word(),
                fileName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dimensioned<scalar>(dimless/dimTime, 1.0)
            *(beta_ - mag(grad(psi)))
        )
    );

    if (etaCutoff_ && alpha_ != -1)
    {
        const volScalarField maggradPsi = mag(grad(psi));
        forAll(tfield(), ID)
        {
            tfield.ref()[ID] *= etaCutoff(maggradPsi[ID], beta_ + alpha_);        
        }
        
    }
    if (psiCutoff_)
    {
        forAll(tfield(), ID)
        {
            tfield.ref()[ID] *= psiCutoff(psi[ID]);        
        }
        
    }

    return tfield;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
