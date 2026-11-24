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
#include "dynamicFvMesh.H"
#include "advectionErrors.H"
#include "phaseIndicator.H"
#include "narrowBand.H"
#include "velocityModel.H"
#include "prescribedVelocityModels.H"
#include "slAdvection.H"
#include "velocityExtension.H"
#include "slMeshTransfer.H"
#include "slRefinementBand.H"
#include "levelSetImplicitSurfaces.H"
#include "profile.H"

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


#include "leiaVersionRegistry.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Semi-Lagrangian level set equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    Foam::leia::reportVersions(Info);   // which leia libraries, from which commit
    Foam::leia::writeVersions(runTime); // the same lines in <case>/leia.version
    // The mesh is a dynamicFvMesh (composition root of the mesh-motion family, added
    // 2026-11-23 for the AMR line). Without constant/dynamicMeshDict, or with
    // `dynamicFvMesh staticFvMesh`, the mesh is static and mesh.update() does nothing:
    // every existing study is unchanged (gated by a bit-identity run, STATUS 11.44).
    #include "createDynamicFvMesh.H"

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

    if (dynamicMesh)
    {
        #include "initialAdaptation.H"
    }

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

        // Trace-velocity source. With projectedFlux the SAME two time levels
        // are built from the prescribed FLUX instead of the prescribed cell
        // velocity: phi is already rescaled to t^{n+1} above, and the old level
        // is the base flux rescaled to t^n (identity for the steady field).
        if (traceFromFlux && traceFluxExtension)
        {
            // Both levels with the psi^n geometry. For the reversed flow the
            // extension of u^n at the old level FIRST, from U and phi set to t^n
            // with the velocity model's own expression, then U and phi restored
            // bit for bit to t^{n+1} (the expression of oscillateVelocity) and
            // the extension of u^{n+1} at the new level. Two evaluations per
            // step, and the extension object ends the step at t^{n+1}.
            if (velocityModel->isOscillating())
            {
                const scalar tau = velocityModel->tau();
                const scalar fn = velocityModel->oscillationFactor
                (
                    runTime.value() - runTime.deltaT().value(), tau
                );
                const scalar fn1 =
                    velocityModel->oscillationFactor(runTime.timeOutputValue(), tau);
                phi == phi0Ptr()*fn;
                U == U0*fn;
                velExtPtr->correct();
                UtraceOldPtr() == fvc::reconstruct(velExtPtr->phi());
                phi == phi0Ptr()*fn1;
                U == U0*fn1;
                velExtPtr->correct();
                UtracePtr() == fvc::reconstruct(velExtPtr->phi());
            }
            else
            {
                velExtPtr->correct();
                UtracePtr() == fvc::reconstruct(velExtPtr->phi());
                UtraceOldPtr() == UtracePtr();
            }
        }
        else if (traceFromFlux)
        {
            UtracePtr() == fvc::reconstruct(phi);
            if (velocityModel->isOscillating())
            {
                const scalar tn = runTime.value() - runTime.deltaT().value();
                UtraceOldPtr() == fvc::reconstruct
                (
                    phi0Ptr()
                   *velocityModel->oscillationFactor(tn, velocityModel->tau())
                );
            }
            else
            {
                UtraceOldPtr() == UtracePtr();
            }
        }

        // Semi-Lagrangian update: psi holds psi^n on entry, psi^{n+1} on exit.
        if (traceFromFlux)
        {
            slAdv->advect(psi, UtracePtr(), UtraceOldPtr());
        }
        else
        {
            slAdv->advect(psi, U, Uold);
        }

        // Mesh adaptation for psi^{n+1}, at the time t^{n+1}, before anything of
        // the next step reads the mesh. The time and psi are at the same level
        // here, so a time-dependent patch value of psi (rigidMotionProfile) is
        // consistent with the cell values when the transfer evaluates it. Until
        // 2026-11-24 this block ran at the start of the step, after the time
        // increment, with psi^n; the sequence of advection and adaptation was
        // the same. A static mesh returns false, and none of the blocks runs.
        if (dynamicMesh)
        {
            refinementBand->calc(psi);     // the indicator that the mesh reads
            meshTransfer->store(psi);      // what the transfer of psi needs
        }
        mesh.update();
        if (dynamicMesh && mesh.topoChanging())
        {
            #include "rebuildAfterMeshChange.H"
        }

        // Narrow band FIRST, in every step, then the phase indicator: the order
        // of the two-phase SL solver (slAlphaEqn.H), whose interface step this
        // kinematic solver must reproduce. Both phase indicators read the
        // registered NarrowBand field to decide which cells get the geometric
        // volume fraction, so the band must be the band of psi^{n+1}. Until
        // 2026-10-01 this solver computed alpha first, with the band of psi^n,
        // and refreshed the band only on report steps: a cell that had just
        // entered the band got a sign-based 0/1 alpha for one step (MEASURED
        // 2026-09-30, STATUS 11.22: the relative volume error 6.5e-04 against
        // 9.2e-05 at step 53 of 2Dvortex N = 32). The default trajectory does
        // not read the band, so psi does not change; trajectoryVelocity
        // normalClosestPoint reads it in the next advection.
        narrowBand->calc();

        // Diagnostics (phase indicator alpha, error norms + CSV row). Alpha is
        // not used by the advection, so with reportAtWriteTimesOnly true it can
        // be skipped between write times (the per-cell-LLS + tet-fill indicator
        // dominates the cost on fine/polyhedral meshes) -- but the final step
        // is ALWAYS logged: the reversed benchmarks are scored exactly at
        // t = endTime, and this solver's own deltaT override keeps
        // adjustableRunTime from ever landing a writeTime() there.
        const bool finalStep =
            (runTime.endTime().value() - runTime.value())
          < 0.5*runTime.deltaTValue();
        if (!reportAtWriteTimesOnly || runTime.writeTime() || finalStep)
        {
            if (dynamicMesh)
            {
                #include "referenceOnPresentMesh.H"
            }
            phaseInd->calcPhaseIndicator(alpha, psi);

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
