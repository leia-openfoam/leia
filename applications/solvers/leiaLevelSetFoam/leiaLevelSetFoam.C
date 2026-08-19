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
    The UNIFIED level-set advection solver: the complete method is
    configured from the fvSolution levelSet dictionary --

    \verbatim
    levelSet
    {
        advection         { type eulerian; }   // eulerian | semiLagrangian
        velocityExtension { type none; ... }   // eulerian only
        sdplsSource       { type noSource; ... } // eulerian only
        redistancer       { type noRedistancing; ... }
        phaseIndicator    { type detrixheAslam; }
        narrowBand        { type signChange; }
        reportAtWriteTimesOnly false;          // gate alpha/error reporting
    }
    \endverbatim

    eulerian: implicit FV transport in the advective form
    ddt(psi) + div(v psi) - psi div(v) = S_SDPLS with a deferred-correction
    loop (controlDict nDefCorr) for second-order (linearUpwind) div schemes.
    semiLagrangian: characteristic update psi^{n+1}(x_c) = psi^n(x_d) with a
    runtime-selectable foot-point reconstruction.

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
#include "velocityModel.H"
#include "prescribedVelocityModels.H"
#include "fluxCorrection.H"
#include "levelSetAdvection.H"
#include "volumeCorrection.H"

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
        "Unified level set equation solver (dictionary-configured method)."
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

    // Per-time-step uniform-shift phase-volume correction. THE PRESCRIBED-
    // VELOCITY TRACK is this model's legitimate home: it restores total phase
    // volume by a spatially CONSTANT shift psi <- psi + eps, and because eps is
    // constant, grad(psi + eps) = grad(psi), so |grad psi| is untouched -- it
    // conserves volume without spending the signed-distance property SDPLS
    // exists to maintain. Target captured here, from the alpha createFields.H
    // has just built, so it is the SAME functional advectionErrors.H reports.
    //
    // Default noVolumeCorrection: an absent levelSet/volumeCorrection
    // sub-dictionary leaves psi untouched for the whole run, so every existing
    // kinematic study is bit-unchanged.
    autoPtr<volumeCorrection> volCorr = volumeCorrection::New(mesh);
    volCorr->setTargetVolume(alpha);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // CFL based deltaT setting from the flux that actually advects psi
    // (extension flux for eulerian-with-extension, base flux otherwise).
    if (runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false))
    {
        runTime.setDeltaT(
            maxDeltaT(adv->advectingFlux(), runTime.controlDict()), false);
    }

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "errorCalculation.H"

    while (runTime.run())
    {
        #include "CourantNo.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << endl;
        Info<< "deltaT = " << runTime.deltaT().value() << nl << endl;

        // u^{n+1}: rescale U (and phi) by cos(pi t^{n+1}/tau) for the
        // reversed (oscillating) tests; identity for steady velocities.
        if (velocityModel->isOscillating())
        {
            velocityModel->oscillateVelocity(U, U0, phi, phi0, runTime);
        }

        // psi^n -> psi^{n+1} by the dictionary-selected advection method.
        adv->advance(psi);

        // Fresh sign-change band BEFORE redistancing (the plane-anchored
        // redistancer models and the signedDistanceBounds criterion read
        // the registered "NarrowBand" field).
        narrowBand->calc();

        const bool reportNow =
            !reportAtWriteTimesOnly || runTime.writeTime();

        if (reportNow)
        {
            phaseInd->calcPhaseIndicator(alpha, psi);
        }

        const scalar Valpha0 =
            reportNow
          ? gSum(alpha.primitiveField()*mesh.V().field())
          : 0.0;

        if (redist->redistance(psi))
        {
            // A redistance event may realign cut cells sub-h: refresh the
            // band and the volume fraction, and report the interface volume
            // the event shifted.
            narrowBand->calc();
            phaseInd->calcPhaseIndicator(alpha, psi);

            if (reportNow)
            {
                const scalar Valpha1 =
                    gSum(alpha.primitiveField()*mesh.V().field());

                Info<< "redistance volume shift: " << Valpha0
                    << " -> " << Valpha1
                    << " (dV = " << Valpha1 - Valpha0
                    << ", dV/V = "
                    << (Valpha1 - Valpha0)/max(Valpha0, VSMALL)
                    << ")" << endl;
            }
        }

        // Uniform-shift volume correction, AFTER redistancing (which rebuilds
        // psi from the zero level set and is NOT volume preserving, so it
        // would discard an eps chosen before it) and BEFORE reportErrors, so
        // the volume AND shape errors describe the field the method actually
        // delivers. This ordering matters for honesty as much as correctness:
        // the correction trades local interface position for global volume by
        // construction, and measuring before it would hide the trade.
        if (volCorr->correct(psi, phaseInd()))
        {
            // The zero crossing moved sub-h: the sign-change band can gain or
            // lose a cell, exactly as after a redistance event.
            narrowBand->calc();

            if (reportNow)
            {
                phaseInd->calcPhaseIndicator(alpha, psi);
            }
        }

        if (reportNow)
        {
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
            runTime.setDeltaT(
                maxDeltaT(adv->advectingFlux(), runTime.controlDict()),
                false);
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
