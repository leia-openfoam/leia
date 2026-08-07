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

#include "discretization.H"
#include "addToRunTimeSelectionTable.H"
#include "sdplsSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(discretization, false);
    defineRunTimeSelectionTable(discretization, Dictionary);
    addToRunTimeSelectionTable(discretization, discretization, Dictionary);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<discretization> Foam::discretization::New(const word type)
{
    auto* ctorPtr = DictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
        << "Unknown discretization type \"" << type << "\"\n\n"
        << "Valid write types : "
        << flatOutput(DictionaryConstructorTablePtr_->sortedToc()) << nl
        << exit(FatalError);
    }
    Info << "Selecting SDPLS Source discretization: " << type << nl << endl;
    return autoPtr<discretization>(ctorPtr());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discretization::discretization()
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

Foam::tmp<fvScalarMatrix> 
Foam::discretization::
discretize(const volScalarField& nonLinearPart, const volScalarField& psi) const
{
    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            psi,
            psi.dimensions()*dimVol/dimTime
        )
    );

    fvScalarMatrix& fvm = tfvm.ref();
    const fvMesh& mesh = psi.mesh();

    // SIGN CONVENTION (this cost every published SDPLS number once).
    //
    // The caller writes  ddt + div - Sp(divPhi) == fvmsdplsSource(psi, U),
    // and fvMatrix::operator== is A - B with diag_ and source_ subtracted
    // termwise. Solving (A-B) gives
    //
    //     (A.M psi - A.source) - (B.M psi - B.source) = 0
    //  => V * Dpsi/Dt = B.diag * psi - B.source
    //
    // so the physical source represented by THIS matrix is
    // (diag*psi - source)/V, NOT (diag*psi + source): a term placed in
    // source_ enters the transport equation with a MINUS. OpenFOAM's own
    // fvm::Su(su, vf) encodes exactly this by doing `source_ -= V*su`.
    //
    // The intended model is S = Sc + Sp*psi^{n+1} (see the class docs), hence
    // the negation below. Without it, `explicit` solved Dpsi/Dt = -f_nl*psi --
    // DOUBLING the gradient drift instead of cancelling it -- and
    // `strictNegativeSpLinearImplicit` did so on its f_nl > 0 branch only,
    // which is why its band min |grad psi| collapsed to ~0.05 (mesh
    // independently) while the implicit branch correctly held the band max
    // near 1. `simpleLinearImplicit` (Sc = 0) was never affected.
    fvm.source()    = -mesh.V().field() * Sc(nonLinearPart, psi);
    fvm.diag()      = mesh.V().field() * Sp(nonLinearPart);
    return tfvm;
}

Foam::tmp<scalarField> 
Foam::discretization::
Sc(const volScalarField& nonLinearPart, const volScalarField& psi) const
{
    return tmp<scalarField>(new scalarField(nonLinearPart.size(), 0.0));
}

Foam::tmp<scalarField> 
Foam::discretization::Sp(const volScalarField& nonLinearPart) const
{
    return tmp<scalarField>(new scalarField(nonLinearPart.size(), 0.0));
}


bool Foam::discretization::iterative()
{
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
