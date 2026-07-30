/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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
    leiaRedistancedLevelSetFoam

Description
    Level-set advection with criterion-gated geometric redistancing.

    The advective (Hamilton-Jacobi) law Dpsi/Dt = 0 is assembled in the
    conservative-minus-compression form

    \f[
        \ddt{\psi} + \div(\vec{v} \psi) - \psi \div{\vec{v}} = 0
    \f]

    with the plain prescribed volumetric flux (optionally Helmholtz-projected
    to a solenoidal flux with -fluxCorrection). No velocity extension: the
    signed-distance property of psi is restored by the levelSet.redistancer
    model instead (e.g. planeFootWave -- the linear least-squares interface
    planes of the Detrixhe-Aslam phase indicator, propagated by MeshWave),
    gated by the levelSet.redistancer.trigger criterion.

    \heading Required fields
    \plaintable
        psi       | Level set (initial signed distance)
        alpha     | Phase indicator (volume fraction)
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
        "Level set equation solver with geometric redistancing."
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

        // The level set obeys the ADVECTIVE (Hamilton-Jacobi) law Dpsi/Dt = 0,
        // not a conservation law. The conservative fvm::div(phi, psi) equals
        // (v.grad)psi only for solenoidal v; prescribed verification
        // velocities are solenoidal analytically but not discretely.
        // Assemble the advective derivative exactly,
        //   (v.grad)psi = div(phi psi) - psi (div phi),
        // so the residual discrete divergence cannot act as a spurious
        // compression source (-fluxCorrection additionally projects phi).
        const volScalarField divPhi("divPhi", fvc::div(phi));

        // Defect-correction loop for the deferred second-order (linearUpwind)
        // spatial term: each pass re-assembles with the latest psi so the
        // explicit (linearUpwind - upwind) correction converges (~2-3 passes
        // -> formal 2nd order). fvm::ddt reuses psi.oldTime(), which is fixed
        // within the step, so re-solving does not corrupt the time derivative.
        const label nDefCorr =
            runTime.controlDict().getOrDefault<label>("nDefCorr", 2);
        for (label corr = 0; corr < nDefCorr; ++corr)
        {
            fvScalarMatrix psiEqn
            (
                fvm::ddt(psi)
                + fvm::div(phi, psi)
                - fvm::Sp(divPhi, psi)
            ==
                source->fvmsdplsSource(psi, U)
            );

            psiEqn.solve();
        }

        // Fresh sign-change band BEFORE redistancing: the plane-anchored
        // redistancer models read the registered "NarrowBand" field to
        // select their anchor cells.
        narrowBand->calc();

        phaseInd->calcPhaseIndicator(alpha, psi);

        const scalar Valpha0 = gSum(alpha.primitiveField()*mesh.V().field());

        if (redist->redistance(psi))
        {
            // A redistance event realigns cut-cell psi to the LLS planes:
            // the sign pattern in band cells may change sub-h -- refresh the
            // band and the volume fraction, and report the interface volume
            // the event shifted.
            narrowBand->calc();
            phaseInd->calcPhaseIndicator(alpha, psi);

            const scalar Valpha1 =
                gSum(alpha.primitiveField()*mesh.V().field());

            Info<< "redistance volume shift: " << Valpha0
                << " -> " << Valpha1
                << " (dV = " << Valpha1 - Valpha0
                << ", dV/V = " << (Valpha1 - Valpha0)/max(Valpha0, VSMALL)
                << ")" << endl;
        }

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
