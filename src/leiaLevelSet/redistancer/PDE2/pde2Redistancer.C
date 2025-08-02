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

#include "pde2Redistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(pde2Redistancer, 0);
    addToRunTimeSelectionTable(redistancer, pde2Redistancer, Mesh);

} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pde2Redistancer::pde2Redistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    nCellsForDiffusion_(redistDict_.getOrDefault<label>("nCellsForDiffusion", 3)),
    deltaX_(Foam::max(1.0/mesh.deltaCoeffs()).value()),
    epsilon_(redistDict_.getOrDefault<dimensionedScalar>("epsilon", dimensionedScalar("epsilon", dimLength, 0.5*nCellsForDiffusion_*deltaX_))),
    deltaTau_(redistDict_.getOrDefault<scalar>("deltaT", 0.5*deltaX_)),
    nIter_(redistDict_.getOrDefault<label>("nIter", 10)),
    maxCo_(redistDict_.getOrDefault<scalar>("maxCo", 1.0)),
    write_(redistDict_.getOrDefault<bool>("write", false))
{
    Info << "Delta X: "  << deltaX_  << nl
         << "epsilon: "  << epsilon_ << nl 
         << "deltaT: " << deltaTau_  << nl << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pde2Redistancer::doRedistance(volScalarField& psi) 
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

    
    dimensionedScalar dimSmall("0", dimless, SMALL);

    volVectorField gradPsi ("grad(psi)", fvc::grad(restartPsi, "grad(psi)"));
    volScalarField magGradPsi ("mag(grad(psi))", mag(gradPsi));
    volVectorField n ("levelSetK", gradPsi/(magGradPsi + dimSmall));


    // Sign of level set
    // Smooth variant
    const volScalarField psiSign = restartPsi/(
                                                Foam::sqrt(
                                                            sqr(restartPsi) 
                                                            + magGradPsi*sqr(epsilon_)
                                                          )
                                                + dimensionedScalar("0", dimLength, SMALL)
                                              ); 
    
    // // Sharp
    // const volScalarField psiSign = Foam::sign(psi); 

    surfaceScalarField w ("psiFlux",
                          (linearInterpolate(psiSign * n) & mesh.Sf())*dimensionedScalar("1", dimLength/dimTime, 1.0)
                          );


    scalarField sumW = fvc::surfaceSum(w)().internalField();

    scalarField rDeltaT = max(
                                1.0/(0.5*deltaX_),
                                sumW/(maxCo_*mesh.V().field())
                            );



    
    // Mark anchoring cells
    volScalarField anchoringCells 
    (
        Foam::IOobject
        (
            "anchoringCells",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false   // Dont add this to the registry
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

    {
        const scalarField& psiIf = psi.internalField();
        scalarField& anchoringCellsIf = anchoringCells.internalFieldRef();
        const labelUList& P = mesh.owner();
        const labelUList& N = mesh.neighbour();

        for (label fi=0; fi<P.size(); fi++)
        {
            const scalar tmp = psiIf[P[fi]]*psiIf[N[fi]];

            if (tmp < 0.0)
            {
                anchoringCellsIf[P[fi]] = 1.0;
                anchoringCellsIf[N[fi]] = 1.0;
            }
        }
    }

    const scalarField& anchoringCellsIf = anchoringCells.internalField();



    // Local Euler scheme lambda_fn 
    auto localEulerScheme = [] (
                                const volScalarField& vf,
                                const scalarField& rDeltaT
                            ) -> tmp<fvScalarMatrix>
    {
        tmp<fvScalarMatrix> tfvm
        (
            new fvScalarMatrix
            (
                vf,
                vf.dimensions()*dimVol/dimTime
            )
        );

        fvScalarMatrix& fvm = tfvm.ref();

        const fvMesh& mesh = vf.mesh();

        fvm.diag() = rDeltaT*mesh.Vsc();
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh.Vsc();
    
        return tfvm;
    };




    // Manual loop, re-writing the old time values for time-stepping

    for (int iter=0; iter<nIter_; iter++)
    {
        
        gradPsi = fvc::grad(restartPsi, "grad(psi)"); 

        magGradPsi = mag(gradPsi);

        n = gradPsi/(magGradPsi + dimSmall);

        // From paper
        w = (linearInterpolate(psiSign * n) & mesh.Sf())*dimensionedScalar("1", dimLength/dimTime, 1.0);

        // Alternative formulation
        // const surfaceScalarField w  = dimensionedScalar("1", dimLength/dimTime, 1.0)*linearInterpolate(psiSign)*fvc::snGrad(restartPsi)/(mag(fvc::snGrad(restartPsi) + dimSmall)) * mesh.magSf();



        fvScalarMatrix redistanceLevelSetEq
        (
            localEulerScheme(restartPsi, rDeltaT)     
            + fvm::div(w, restartPsi, "div(psiFlux, psi)")    
            + fvm::SuSp(-fvc::div(w), restartPsi)
            ==
            psiSign*dimensionedScalar("1", dimLength/dimTime, 1.0)
            
        );
            
        redistanceLevelSetEq.relax();
        
        redistanceLevelSetEq.solve();

        // Reset levelSet values
        forAll(anchoringCellsIf, cellI)
        {
            if (anchoringCellsIf[cellI] == 1.0)
            {
                restartPsi[cellI] = psi[cellI];
            }
        }

        restartPsi.correctBoundaryConditions();
    
        volScalarField& oldTime = restartPsi.oldTime();

        oldTime = restartPsi;        
        oldTime.correctBoundaryConditions();

        sumW = fvc::surfaceSum(w)().internalField();

        rDeltaT = max(
                        1.0/(0.5*deltaX_),
                        sumW/(maxCo_*mesh.V().field())
                    );
    }

    // Copy variable back
    psi.internalFieldRef() = restartPsi.internalField();
    psi.correctBoundaryConditions();

    if (write_)
    {
        gradPsi = fvc::grad(restartPsi, "grad(psi)"); 
        magGradPsi = mag(gradPsi);
        n = gradPsi/(magGradPsi + dimSmall);

        restartPsi.write();
        gradPsi.write();
        n.write();
        anchoringCells.write();
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
