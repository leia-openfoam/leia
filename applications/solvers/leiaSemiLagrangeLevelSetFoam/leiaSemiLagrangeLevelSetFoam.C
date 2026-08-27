/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    leiaSemiLagrangeLevelSetFoam

Description
    Cell-centred semi-Lagrangian level-set advection: psi is constant along
    characteristics, so psi^{n+1}(x_c) = psi^n(x_d) with x_d the Taylor backward
    foot and psi^n at x_d reconstructed by a runtime-selectable slReconstruction
    (linearTaylor | linearWeightedLeastSquares | quadraticTaylor |
    quadraticWeightedLeastSquares). No flux, no divergence, no linear solve; no
    reinitialization (the signed-distance property is not maintained).

    The parallel research line to the velocity-extension solver
    (leiaLevelSetFoam): same prescribed-velocity verification, same reversed
    vortex study, compared on the same error norms.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "advectionErrors.H"
#include "phaseIndicator.H"
#include "narrowBand.H"
#include "velocityModel.H"
#include "prescribedVelocityModels.H"
#include "slAdvection.H"

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
        "Semi-Lagrangian level set equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // CFL based deltaT setting (convective Co from the base flux; the SL foot
    // guard is the operational CFL <= 1 check).
    if (runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false))
    {
        runTime.setDeltaT(maxDeltaT(phi, runTime.controlDict()), false);
    }

    Info<< "\nCalculating semi-Lagrangian scalar transport\n" << endl;

    // Per-step metrics by default. 4044f47 gated them to writeTime() for speed
    // on fine/polyhedral meshes, but (a) this solver overrides deltaT AFTER
    // Time's write-time adjustment, so adjustableRunTime never lands on the
    // write instants, and (b) the endTime step is not a writeTime(), so the
    // reversal endpoint t = T -- the only instant the reversed benchmarks are
    // scored on -- was never logged and the aggregation scored t = T - dt: a
    // dt ~ h gap that manufactures a fake first order (measured 2026-08-27,
    // shape 9.18e-05 at T - dt vs 1.11e-06 at T, N = 256). The template knob
    // reportAtWriteTimesOnly (fvSolution/levelSet, REPORT_WRITE_ONLY) existed
    // but was never read; it now does what it says. The final step is logged
    // unconditionally in either mode.
    const Switch reportAtWriteTimesOnly
    (
        mesh.solutionDict().subOrEmptyDict("levelSet")
            .getOrDefault<Switch>("reportAtWriteTimesOnly", false)
    );

    #include "errorCalculation.H"

    while (runTime.run())
    {
        #include "CourantNo.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << endl;
        Info<< "deltaT = " << runTime.deltaT().value() << nl << endl;

        // u^{n+1}: rescale U (and phi) by cos(pi t^{n+1}/tau) for the reversed
        // (oscillating) test; identity for the steady case.
        if (velocityModel->isOscillating())
        {
            velocityModel->oscillateVelocity(U, U0, phi, phi0Ptr(), runTime);
        }

        // u^n: the previous-step velocity the Taylor trajectory needs. For the
        // oscillating field it is U0 scaled at t^n = t^{n+1} - dt; for the
        // steady field u^n = u^{n+1} = U0.
        if (velocityModel->isOscillating())
        {
            const scalar tn = runTime.value() - runTime.deltaT().value();
            Uold == U0*velocityModel->oscillationFactor(tn, velocityModel->tau());
        }
        else
        {
            Uold == U0;
        }

        // Semi-Lagrangian update: psi holds psi^n on entry, psi^{n+1} on exit.
        slAdv->advect(psi, U, Uold);

        // Diagnostics (phase indicator alpha, narrow band, error norms + CSV
        // row). The advection above uses neither alpha nor the narrow band, so
        // with reportAtWriteTimesOnly true they can be skipped between write
        // times (the per-cell-LLS + tet-fill indicator dominates the cost on
        // fine/polyhedral meshes) -- but the final step is ALWAYS logged: the
        // reversed benchmarks are scored exactly at t = endTime, and this
        // solver's own deltaT override keeps adjustableRunTime from ever
        // landing a writeTime() there.
        const bool finalStep =
            (runTime.endTime().value() - runTime.value())
          < 0.5*runTime.deltaTValue();
        if (!reportAtWriteTimesOnly || runTime.writeTime() || finalStep)
        {
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
        }

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
