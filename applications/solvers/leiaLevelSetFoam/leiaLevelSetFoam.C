/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt 
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

Application
    leiaLevelSetFoam 

Description
    \heading Solver details
    A level-set equation in conservative form   

    \f[
        \ddt{psi} + \div \left(\vec{v} \psi\right) = 0 
    \f]

    Where:
    \vartable
        psi       | Passive scalar
    \endvartable

    \heading Required fields
    \plaintable
        psi       | Passive scalar
        phi       | Volumetric Flux [m^3/s]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "advectionErrors.H"
#include "phaseIndicator.H"
#include "redistancer.H"
#include "narrowBand.H"
#include "sdplsSource.H"
#include "advectionVerification.H"
#include "subCycle.H"

// tmp
#include "fileName.H"
#include "uncollatedFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar maxDeltaT(surfaceScalarField phi, const dictionary& dict)
{
    scalar maxCo;
    if (dict.found("CFL"))
    {
        maxCo = dict.get<scalar>("CFL");
    }
    else
    {
        maxCo = dict.get<scalar>("maxCo");
    }
    scalar maxGrowFactor = 1.2;
    const fvMesh& mesh = phi.mesh();
    scalar deltaT = mesh.time().deltaT().value();
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    scalar deltaT_suggestion = maxCo / (0.5*gMax(sumPhi/mesh.V().field()));

    deltaT = min(deltaT * maxGrowFactor, deltaT_suggestion);

    return deltaT;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Level set equation solver."
    );
    
    argList::addBoolOption
    (
        "fluxCorrection", 
        "Use Helmholz decomposition to enforce a div-free volumetric flux."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // CFL based deltaT setting
    if (runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false))
    {
        runTime.setDeltaT(maxDeltaT(phi, runTime.controlDict()), false);
    }

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "errorCalculation.H"

    // TODO(TM): read from constant/functions/divDefCorr
    //           bool defCorr() member function using solver performance? 
    const label MAX_N_DEF_CORR = 20;

    while (runTime.run())
    {
        #include "CourantNo.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << endl;
        Info<< "deltaT = " << runTime.deltaT().value() << nl << endl;

        if (velocityModel->isOscillating())
        {
            velocityModel->oscillateVelocity(U, U0, phi, phi0, runTime);
        }

        // Applies the defect correction source for the upwind scheme.
        bool defCorr = true; 
        label nDefCorr = 0;
        while (defCorr && (nDefCorr < MAX_N_DEF_CORR))
    	{
            // TODO(TM): remove, ported to upwindSecondOrderDefCorr
            // Reset the error.
            psiErr == dimensionedScalar("psiErr", psi.dimensions(), 0);
            // Compute cell-centered gradient of psi.
            // TODO(TM): contraction for tensor fields or snGrad. 
            gradPsi = fvc::grad(psi); 
            // Get owner-neighbour addressing.
            const auto& own = mesh.owner();
            const auto& nei = mesh.neighbour();
            // Get cell centers.
            const auto& C = mesh.C();
            // Get face centers.
            const auto& Cf = mesh.Cf();

            // For internal faces
            forAll(own, faceI)
            {
                // If flux is positive, owner-cell is upwind. 
                if (phi[faceI] > 0)
                {
                    psiErr[faceI] = (gradPsi[own[faceI]] & (Cf[faceI] - C[own[faceI]]));
                }
                else // If flux is negative, neighbor-cell is upwind. 
                {
                    psiErr[faceI] = (gradPsi[nei[faceI]] & (Cf[faceI] - C[nei[faceI]])); 
                }
            }

            // Computing psiErr on coupled boundaries.
            auto& psiErrBdryField = psiErr.boundaryFieldRef(); 
            // Boundary data for psiErr calculation.
            const auto& gradPsiBdryField = gradPsi.boundaryField(); 
            const auto& cfBdryField = Cf.boundaryField(); 
            const auto& cBdryField = C.boundaryField(); 
            const auto& phiBdryField = phi.boundaryField();
            const auto& patches = mesh.boundary(); 
            const auto& faceOwner = mesh.faceOwner();

            // For all boundary patches (faces).
            forAll(psiErrBdryField, patchI)
            {
                const fvPatch& patch = patches[patchI];

                // Face centers patch field
                const auto& cfPatch = cfBdryField[patchI];

                // Face-centered error patch field for calculation 
                auto& psiErrPatch = psiErrBdryField[patchI];

                // Face-centered volumetric flux patch field
                const auto& phiPatchField = phiBdryField[patchI];

                // Compute psiErr on all outflow patches
                forAll(phiPatchField, faceI)
                {
                    const label faceG = faceI + patch.start(); // Global label. 
                    // If flux is positive, owner-cell is upwind. 
                    if (phiPatchField[faceI] > 0) // If flux is positive
                    {
                        psiErrPatch[faceI] = 
                        (   
                            gradPsi[faceOwner[faceG]] &  
                            (cfPatch[faceI] - C[faceOwner[faceG]])
                        );
                    }
                }

                if (isA<coupledFvPatch>(patch)) // coupled patch
                {
                    // Get gradPsi across coupled patch boundary 
                    const auto& gradPsiPatch = gradPsiBdryField[patchI];
                    auto gradPsiNeiTmp = gradPsiPatch.patchNeighbourField();
                    const auto& gradPsiNei = gradPsiNeiTmp();

                    // Get cell centers across coupled patch boundary 
                    const auto& cPatch = cBdryField[patchI];
                    auto cNeiTmp = cPatch.patchNeighbourField();
                    const auto& cNei = cNeiTmp();

                    // Compute psiErr on coupled patch.
                    forAll(phiPatchField, faceI)
                    {
                        if (phiPatchField[faceI] < 0) // If flux is negative 
                        {
                            // Coupled-patch neighbor is upwind
                            psiErrPatch[faceI] = 
                            (   
                                gradPsiNei[faceI] &  
                                (cfPatch[faceI] - cNei[faceI])
                            );
                        }
                    }
                }
            }
	    
            fvScalarMatrix psiEqn
            (
                fvm::ddt(psi) + fvm::div(phi, psi)
            	//== source->fvmsdplsSource(psi, U)
            );

            auto eqnSolverPerf = psiEqn.solve();
            // TODO(TM): obtain solver performance from the equation in fvOption.
            defCorr = (eqnSolverPerf.nIterations() > 0);  
            ++nDefCorr;

            psi.correctBoundaryConditions(); 
        }
        
        redist->redistance(psi);
        
        phaseInd->calcPhaseIndicator(alpha, psi);

        narrowBand->calc();

        reportErrors(
            errorFile, 
            psi, 
            psi0, 
            alpha, 
            alpha0, 
            phi,
            CoNum
        );

        runTime.write();
        runTime.printExecutionTime(Info);
        

        // CFL based deltaT setting
        if (runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false))
        {
            runTime.setDeltaT(maxDeltaT(phi, runTime.controlDict()), false);
        }

        // if last timestep would overshoot endTime, set deltaT
        if ((runTime.endTime() - runTime) < runTime.deltaT())
        {
            runTime.setDeltaT((runTime.endTime() - runTime), false);
        }

    }

    psi.write();
    alpha.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
