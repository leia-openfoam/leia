/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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

#include "pde3Redistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(pde3Redistancer, 0);
    addToRunTimeSelectionTable(redistancer, pde3Redistancer, Mesh);

} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pde3Redistancer::pde3Redistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    nCellsForDiffusion_(redistDict_.getOrDefault<label>("nCellsForDiffusion", 3)),
    deltaX_(Foam::max(1.0/mesh.deltaCoeffs()).value()),
    epsilon_(redistDict_.getOrDefault<dimensionedScalar>("epsilon", dimensionedScalar("epsilon", dimLength, 0.5*nCellsForDiffusion_*deltaX_))),
    deltaTau_(redistDict_.getOrDefault<scalar>("deltaT", 0.5*deltaX_)),
    nIter_(redistDict_.getOrDefault<label>("nIter", 10)),
    write_(redistDict_.getOrDefault<Switch>("write", false)),
    constrain_(redistDict_.getOrDefault<Switch>("constrainSystem", false))
{
    Info << "constrainSystem: " << constrain_ << nl
         << "Delta X: "  << deltaX_  << nl
         << "epsilon: "  << epsilon_ << nl 
         << "deltaT: " << deltaTau_  << nl << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pde3Redistancer::doRedistance(volScalarField& psi) 
{
   // Grab mesh
    const fvMesh& mesh = psi.mesh();

    // Grab time
    const Time& runTime = mesh.time();

    // Copy levelSet into another field for SS solution
    volScalarField restartPsi 
    (
        Foam::IOobject
        (
            "psiRestart",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false   // Dont add this to the registry
        ),
        mesh,
        psi.dimensions(),
        psi.boundaryField().types()
    );

    // Copy internalField and update BC
    restartPsi.internalFieldRef() = psi.internalField();
    restartPsi.correctBoundaryConditions();

    auto populateForwardBackwardDetivative = [] 
    (
        vectorField& forwardDiff,
        vectorField& backwardDiff,
        const volScalarField& psi
    ) -> void
    {
      // set reference to difference factors array
        const fvMesh& mesh = psi.mesh();  

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

        const surfaceVectorField deltas = mesh.delta();

        forAll(owner, faceI)
        {
            const vector& dPN = deltas[faceI];

            const vector direction = tensor::I & dPN;   

            label axisID {0};

            scalar directionValue{0};

            forAll(direction, i)
            {
                if (mag(direction[i]) > mag(directionValue))
                {
                    directionValue = direction[i];
                    axisID = i;
                }
            }

            if (directionValue > 0.0) // Points along axis
            {
                forwardDiff[owner[faceI]][axisID] = deltaCoeffs[faceI]*(psi[neighbour[faceI]] - psi[owner[faceI]]);
                backwardDiff[neighbour[faceI]][axisID] = deltaCoeffs[faceI]*(psi[neighbour[faceI]] - psi[owner[faceI]]);
            }
            else
            {
                forwardDiff[neighbour[faceI]][axisID] = deltaCoeffs[faceI]*(psi[owner[faceI]] - psi[neighbour[faceI]]);
                backwardDiff[owner[faceI]][axisID] = deltaCoeffs[faceI]*(psi[owner[faceI]] - psi[neighbour[faceI]] );
            }
        }

        forAll(psi.boundaryField(), patchI)
        {
            const fvPatchField<scalar>& psf = psi.boundaryField()[patchI];

            const Foam::labelUList& cellOwners = psf.patch().faceCells();

            const vectorField pDelta = psf.patch().delta();

            const scalarField& pDeltaCoeffs = deltaCoeffs.boundaryField()[patchI];


            forAll(cellOwners, faceI)
            {
                const vector& dPN = pDelta[faceI];

                const vector direction = tensor::I & dPN;  

                label axisID {0};

                scalar directionValue{0};

                forAll(direction, i)
                {
                    if (mag(direction[i]) > mag(directionValue))
                    {
                        directionValue = direction[i];
                        axisID = i;
                    }
                }

                const scalar faceValue = psi.boundaryField()[patchI][faceI];

                if (directionValue > 0.0) 
                {
                    forwardDiff[cellOwners[faceI]][axisID] = pDeltaCoeffs[faceI]*(faceValue - psi[cellOwners[faceI]]);
                }
                else
                {
                    backwardDiff[cellOwners[faceI]][axisID] = pDeltaCoeffs[faceI]*(psi[cellOwners[faceI]] - faceValue);
                }
            }
        }

        // Remove empty directions
        typename vector::labelType validComponents
        (
            mesh.template validComponents<vector>()
        );

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (validComponents[cmpt] == -1)
            {
                forwardDiff.replace
                (
                    cmpt,
                    0.0
                );

                backwardDiff.replace
                (
                    cmpt,
                    0.0
                );
            }
        }
    };


    auto HamiltonianGodunovGrad = [&populateForwardBackwardDetivative] 
    (
        const volScalarField& psi0,
        const volScalarField& psi
    ) -> tmp<volVectorField>
    {

        const fvMesh& mesh = psi.mesh();

        auto tgradPsi = volVectorField::New
        (
            "grad(" + psi.name() + ')',
            mesh,
            dimensionedVector("0", psi.dimensions()/dimLength, vector::zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        );

        auto& gradPsi = tgradPsi.ref();

        vectorField forwardDiff(gradPsi.size(), vector::zero);
        vectorField backwardDiff(gradPsi.size(), vector::zero);

        // Compute forward and backward differences
        populateForwardBackwardDetivative(forwardDiff, backwardDiff, psi);

        // Build the Gradient 
        forAll(gradPsi, cellI)
        {
            if (psi0[cellI] > 0)
            {
                for(label i=0; i<3; i++)
                {
                    gradPsi[cellI][i] = max(
                                                mag(max(backwardDiff[cellI][i], 0)),
                                                mag(min(forwardDiff[cellI][i], 0))
                                            );   
                }
            }
            else if (psi0[cellI] < 0)
            {
                for(label i=0; i<3; i++)
                {
                    gradPsi[cellI][i] = max(
                                                mag(min(backwardDiff[cellI][i], 0)),
                                                mag(max(forwardDiff[cellI][i], 0))
                                            );   
                }
            }
            else
            {
                gradPsi[cellI] = vector::zero;
            }
        }



        return tgradPsi;
    };
    

    auto HamiltonianGodunov = [&HamiltonianGodunovGrad]
    (
        const volScalarField& psi0,
        const volScalarField& psi
        ) -> tmp<volScalarField>
    {
        tmp<volScalarField> tvsf
        (
            new volScalarField
            (
                IOobject
                (
                    "magGrad(" + psi.name() + ')',
                    psi.mesh().time().timeName(),
                    psi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                psi.mesh(),
                dimensionedScalar("0", dimless, 0.0)
            )
        );

        volScalarField& sf = tvsf.ref();

        // Compute gradient
        tmp<volVectorField> tgradPsi = HamiltonianGodunovGrad(psi0, psi);
        const volVectorField& gradPsi = tgradPsi.ref();
     
        // Compute magnitude
        sf.primitiveFieldRef() = mag(gradPsi.primitiveField());

        return tvsf; 
    };

        
    // Compute sign of psi
    const volScalarField psiSign = psi/sqrt(   
                                              sqr(psi)*dimensionedScalar("1", dimless/dimArea, 1.0) 
                                            + sqr(epsilon_)*dimensionedScalar("1", dimless/dimArea, 1.0) 
                                           ); 

    for (int iter=0; iter<nIter_; iter++)
    {

        if (constrain_) // Force to keep original sign
        {
            forAll(psiSign, cellI)
            {
                const scalar curSign = restartPsi[cellI];
                if (curSign != psiSign[cellI])
                {
                    restartPsi[cellI] = sign(psiSign[cellI])*mag(restartPsi[cellI]);
                }
            }
        }

        
        tmp<volScalarField> tmagGradPsi = HamiltonianGodunov(psiSign, restartPsi);
        const volScalarField& magGradPsi = tmagGradPsi();

        const volScalarField rhs = psiSign*(1.0 - magGradPsi);

        // // Explicit euler
        restartPsi += deltaTau_*rhs;
        restartPsi.correctBoundaryConditions();
    }

    // Copy variable back to original lvlSet
    psi.internalFieldRef() = restartPsi.internalField();
    psi.correctBoundaryConditions();

    // Compute gradient for output
    tmp<volVectorField> tgradPsi = HamiltonianGodunovGrad(psi, psi);
    const volVectorField& gradPsi = tgradPsi.ref();

    if (write_)
    {
        psi.write();
        gradPsi.write();
        psiSign.write();


        volVectorField forwardDiff("forwardDiff", 0.0*gradPsi);
        volVectorField backwardDiff("backwardDiff", 0.0*gradPsi);

        populateForwardBackwardDetivative(forwardDiff.internalFieldRef(), backwardDiff.internalFieldRef(), psi);

        forwardDiff.write();
        backwardDiff.write();

    }
}


