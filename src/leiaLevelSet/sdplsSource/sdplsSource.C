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

#include "sdplsSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvSolution.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sdplsSource, false);
    defineRunTimeSelectionTable(sdplsSource, Dictionary);
    addToRunTimeSelectionTable(sdplsSource, sdplsSource, Dictionary);
} 

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<sdplsSource>
Foam::sdplsSource::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution (mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& sourceTermDict = 
        levelSetDict.optionalSubDict("sdplsSource");
    const word& type = sourceTermDict.getOrDefault<word>("type", "noSource");

    auto* ctorPtr = DictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "sdplsSource",
            type,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // A source with no discretization is a silent no-op: `discretization none`
    // resolves to the BASE discretization, whose Sc and Sp are both zero, so
    // the assembled matrix contributes nothing to psiEqn -- while the log still
    // reports "Selecting SDPLS source term type: R". Refuse the combination
    // rather than let a study measure an inactive source and call it SDPLS.
    const word& discretizationType =
        sourceTermDict.getOrDefault<word>("discretization", "none");

    Info << "Selecting SDPLS source term type: " << type << nl << endl;
    autoPtr<sdplsSource> model(ctorPtr(sourceTermDict, mesh));

    // ...EXCEPT for a source that assembles its own matrix and never calls the
    // discretization hierarchy (usesDiscretization() false: the divergence-form
    // family sdplsRdiv / sdplsRdivStrictSp). For it `discretization none` is NOT
    // a null study, and any other value is read by NOTHING. Until 2026-09-26
    // this test compared the type name with the literals "Rdiv" and
    // "RdivStrictSp", so every new self-assembling model had to edit this file.
    const bool divergenceFormSource = !model->usesDiscretization();

    if (divergenceFormSource && discretizationType != "none")
    {
        // Not fatal: every committed Rdiv study carries a discretization
        // keyword (the case templates always render one), and those runs must
        // stay runnable and behaviourally unchanged. But the keyword is inert
        // here, while aggregate.py and method_label.py copy SOURCE_SCHEME into
        // the curated CSV and the method label -- so without this line a table
        // row would claim a linearization the assembly never performed.
        WarningInFunction
            << "sdplsSource type '" << type << "' assembles its own matrix"
            << " (divergence form) and IGNORES"
            << " levelSet.sdplsSource.discretization, which is set to '"
            << discretizationType << "'." << nl
            << "The run is unaffected; the dictionary entry and any method"
            << " label derived from it are not descriptive of what was"
            << " assembled." << endl;
    }

    if (type != sdplsSource::typeName && !divergenceFormSource
     && discretizationType == "none")
    {
        FatalIOErrorInFunction(sourceTermDict)
            << "sdplsSource type '" << type << "' was selected with"
            << " discretization 'none'." << nl
            << "That combination assembles a zero matrix: the source term is"
            << " inactive and the run is indistinguishable from 'noSource'," << nl
            << "which is exactly the kind of silently-null study this guard"
            << " exists to prevent." << nl
            << "Set levelSet.sdplsSource.discretization to one of"
            << " simpleLinearImplicit | strictNegativeSpLinearImplicit |"
            << " explicit | exponential | exponentialImplicit," << nl
            << "or set type to noSource if an inactive source is intended."
            << exit(FatalIOError);
    }

    return model;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sdplsSource::sdplsSource(const dictionary& dict, const fvMesh& mesh)
    :
        discretization_
            (
                discretization::
                New(dict.getOrDefault<word>("discretization","none"))
            ),
        gradPsi_(gradPsi::New(dict.getOrDefault<word>("gradPsi","fvc"), mesh)),
        mollifier_(mollifier::New(dict.subOrEmptyDict("mollifier"))),
        R_(
            IOobject(
                "R",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ), 
            mesh, 
            dimensioned(dimless/dimTime, 0.0)
        ),
        nonLinearPart_(
            IOobject
            (
                "SDPLS_nonLinearPart",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 0.0)
        ),
        sourceDict_(dict)
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

Foam::tmp<fvScalarMatrix> 
Foam::sdplsSource::
fvmsdplsSource(const volScalarField& psi, const volVectorField& U)
{
    update(psi, U);
    tmp<fvScalarMatrix> tfvm;
    if(mollifier_->type() == "none")
    {
        tfvm = discretization_->discretize(nonLinearPart_, psi);
    }
    else
    {
        volScalarField mollified_nonLinearPart =
            nonLinearPart_*mollifier_->field(psi);
        tfvm = discretization_->discretize(mollified_nonLinearPart, psi);
    }

    return tfvm;
}

void
Foam::sdplsSource::update(const volScalarField& psi, const volVectorField& U)
{
    R_ = R(psi, U);
    nonLinearPart_ = nonLinearPart(R_, psi, U);
}


Foam::tmp<volScalarField> 
Foam::sdplsSource::R(const volScalarField& psi, const volVectorField& U) const
{
    volVectorField const grad_psi = grad(psi);
    // Dedicated fvSchemes keyword (was the unnamed "grad(U)" fallback):
    // R = (grad(U) & nHat) & nHat needs independent scheme control.
    volTensorField const grad_U = fvc::grad(U, "gradUSdpls").cref();
    dimensioned<scalar> const eps = dimensioned(grad_psi.dimensions(), SMALL);
    volVectorField const normal_interface = 
        (grad_psi/(mag(grad_psi) + eps)).cref();

    return (grad_U & normal_interface) & normal_interface;
}

Foam::tmp<volVectorField>
Foam::sdplsSource::grad(const volScalarField& psi) const
{
    return gradPsi_->grad(psi);
}

Foam::tmp<volScalarField> 
Foam::sdplsSource::nonLinearPart
    (
        const volScalarField& R,
        const volScalarField& psi,
        const volVectorField& U
    ) const
{
    const fvMesh& mesh = psi.mesh();
    return tmp<volScalarField>
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
            mesh,
            dimensionedScalar(nonLinearPart_.dimensions(), 0.0)
        )
    );
}

void Foam::sdplsSource::write() const
{
    R_.write();
    nonLinearPart_.write();
}

bool Foam::sdplsSource::iterative()
{
    return discretization_->iterative();
}

uint Foam::sdplsSource::maxIterations()
{
    if (iterative())
    {
        return sourceDict_.getOrDefault<uint>("maxIterations", 50);
    }
    else
    {
        return 1;
    }
}

const Foam::surfaceScalarField& Foam::sdplsSource::transportFlux() const
{
    if (!transportFlux_)
    {
        FatalErrorInFunction
            << "sdplsSource type '" << type() << "' needs the face flux that"
            << " transports psi, but the solver set none: the composition root"
            << " must call setTransportFlux() before fvmsdplsSource()."
            << exit(FatalError);
    }
    return *transportFlux_;
}


// * * * * * * * * * * * * * *  Global functions  * * * * * * * * * * * * * * //

Foam::tmp<fvScalarMatrix> 
Foam::fvm::sdplsSource(const volScalarField& psi, const volVectorField& U)
{
    return sdplsSource::New(psi.mesh()).ref().fvmsdplsSource(psi, U);
}

// ************************************************************************* //




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
