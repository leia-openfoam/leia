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

#include "pde1Redistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(pde1Redistancer, 0);
    addToRunTimeSelectionTable(redistancer, pde1Redistancer, Mesh);

} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pde1Redistancer::pde1Redistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    nCellsForDiffusion_(redistDict_.getOrDefault<label>("nCellsForDiffusion", 3)),
    deltaX_(Foam::max(1.0/mesh.deltaCoeffs()).value()),
    epsilon_(redistDict_.getOrDefault<dimensionedScalar>("epsilon", dimensionedScalar("epsilon", dimLength, 0.5*nCellsForDiffusion_*deltaX_))),
    deltaTau_(redistDict_.getOrDefault<scalar>("deltaT", 0.5*deltaX_)),
    nIter_(redistDict_.getOrDefault<label>("nIter", 10)),
    write_(redistDict_.getOrDefault<bool>("write", false))
{
    Info << "Delta X: "  << deltaX_  << nl
         << "epsilon: "  << epsilon_ << nl 
         << "deltaT: " << deltaTau_  << nl << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pde1Redistancer::doRedistance(volScalarField& psi) 
{
    // Grab mesh
    const fvMesh& mesh = psi.mesh();

    // Grab time
    const Time& runTime = mesh.time();

    // Time& runTime = const_cast<Time&>(mesh.time());

    // // Store real time dT
    // const scalar deltaT = runTime.deltaT().value();

    // // Manually set deltaT
    // runTime.setDeltaT
    // (
    //     deltaTau_,
    //     false
    // );

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


    volVectorField gradPsi ("grad(psi)", fvc::grad(restartPsi, "grad(psi)"));
    volScalarField magGradPsi ("mag(grad(psi))", mag(gradPsi));

    // Sign of level set
    // Smooth variant
    const volScalarField psiSign = restartPsi/(
                                                Foam::sqrt(
                                                            sqr(restartPsi) 
                                                            + magGradPsi*sqr(epsilon_)
                                                          )
                                              ); 
    
    // // Sharp
    // const volScalarField psiSign = Foam::sign(psi); 

    // Euler scheme via lambda_fn 
    auto eulerScheme = [] (
                            const volScalarField& vf,
                            const scalar deltaT
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

        const scalar rDeltaT = 1.0/deltaT;

        fvm.diag() = rDeltaT*mesh.Vsc();

        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh.Vsc();
    
        return tfvm;
    };


        
    // Manual loop, re-writing the old time values for time-stepping

    for (int iter=0; iter<nIter_; iter++)
    {
        gradPsi = fvc::grad(restartPsi, "grad(psi)"); 
 
        magGradPsi = mag(gradPsi);

        fvScalarMatrix redistanceLevelSetEq
        (
            eulerScheme(restartPsi, deltaTau_)
            ==
            psiSign*(1.0 - magGradPsi)*dimensionedScalar("1", dimLength/dimTime, 1.0)
        );
        
        redistanceLevelSetEq.relax();
        
        redistanceLevelSetEq.solve();

    
        volScalarField& oldTime = restartPsi.oldTime();

        oldTime = restartPsi;        
        oldTime.correctBoundaryConditions();
    }

    // Copy variable back
    psi.internalFieldRef() = restartPsi.internalField();
    psi.correctBoundaryConditions();


    // // Restore simulation deltaT
    // runTime.setDeltaT
    // (
    //     deltaT, false
    // );

    if (write_)
    {
        gradPsi = fvc::grad(restartPsi, "grad(psi)"); 
        magGradPsi = mag(gradPsi);

        restartPsi.write();
        gradPsi.write();
        magGradPsi.write();
    }

}

// ************************************************************************* //
