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
#include "velocityModel.H"
#include "prescribedVelocityModels.H"
#include "fluxCorrection.H"
#include "velocityExtension.H"

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
        // CFL from the flux that actually advects psi (velExt->phi() == phiExt;
        // == phi for the "none" model), not the base phi, so the extension
        // cannot silently exceed maxCo.
        runTime.setDeltaT(maxDeltaT(velExt->phi(), runTime.controlDict()), false);
    }

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "errorCalculation.H"

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

        // Correct the advecting velocity to be interface-normal-constant in a
        // narrow band (non-invasive: U/phi unchanged; advect with velExt->phi()).
        velExt->correct();

        // The level set obeys the ADVECTIVE (Hamilton-Jacobi) law Dpsi/Dt = 0,
        // not a conservation law. The conservative fvm::div(phiExt, psi) equals
        // (v.grad)psi only for solenoidal v; an interface-normal-constant
        // extension velocity is generically NOT solenoidal (div(v0*n) = v0*kappa),
        // so the conservative form hides a spurious compression source
        // psi*div(vExt). Assemble the advective derivative exactly,
        //   (v.grad)psi = div(phiExt psi) - psi (div phiExt),
        // (same idiom as pseudoTime's Uext transport). This also removes the
        // need to Helmholtz-project phiExt: a projection would destroy the
        // (n.grad)Uext = 0 property that preserves |grad psi| = 1 near Sigma
        // (Zhao et al. 1996; Adalsteinsson & Sethian 1999).
        const volScalarField divPhiExt("divPhiExt", fvc::div(velExt->phi()));

        // Defect-correction loop for the deferred second-order (linearUpwind)
        // spatial term: each pass re-assembles with the latest psi so the
        // explicit (linearUpwind - upwind) correction converges (~2-3 passes ->
        // formal 2nd order). Safe: fvm::ddt reuses psi.oldTime(), which is fixed
        // within the step (only updated at ++runTime), so re-solving here does
        // not corrupt the time derivative.
        const label nDefCorr =
            runTime.controlDict().getOrDefault<label>("nDefCorr", 2);
        for (label corr = 0; corr < nDefCorr; ++corr)
        {
            fvScalarMatrix psiEqn
            (
                fvm::ddt(psi)
                + fvm::div(velExt->phi(), psi)
                - fvm::Sp(divPhiExt, psi)
            ==
                source->fvmsdplsSource(psi, U)
            );

            psiEqn.solve();
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
            // CFL from the flux that actually advects psi (velExt->phi() == phiExt;
        // == phi for the "none" model), not the base phi, so the extension
        // cannot silently exceed maxCo.
        runTime.setDeltaT(maxDeltaT(velExt->phi(), runTime.controlDict()), false);
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
