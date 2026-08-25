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
    leiaSemiLagrangianLevelSetTwoPhaseFoam

Group
    grpMultiphaseSolvers

Description
    A solver for two incompressible, isothermal immiscible fluids using a
    SEMI-LAGRANGIAN level set. It is leiaLevelSetTwoPhaseFoam with
    the Eulerian psi transport swapped for the cell-centred semi-Lagrangian
    advection (slAdvection); the interface curvature that drives the CSF surface-
    tension force is evaluated SYMBOLICALLY from the quadratic reconstruction
    (fit gradient + Hessian), constant-extended along the normal, and consumed by
    the reconstructedCurvature surfaceTensionForce model. The level-set
    transport can use an interface-normal-constant velocity extension without
    changing the physical momentum velocity. The optional rhoLENT mass-flux
    mode solves an auxiliary conservative density equation with exactly the
    mass flux used by momentum, then resets cell density from the geometric
    Detrixhe-Aslam indicator after pressure-velocity convergence.

    All momentum/pressure/mesh machinery (createFields.H, UEqn.H, pEqn.H,
    YoungLaplaceEqn.H, error headers) is reused verbatim from the sibling solver
    via -I; only the level-set advection file (slAlphaEqn.H) and the extra fields
    (createSLFields.H) are local.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

// Leia Level Set Method
#include "phaseIndicator.H"
#include "redistancer.H"
#include "sdplsSource.H"     // reused sibling createFields.H constructs an sdplsSource
#include "surfaceTensionForce.H"
#include "narrowBand.H"

// Semi-Lagrangian level-set advection.
#include "slAdvection.H"
#include "velocityExtension.H"

// Geometric face liquid-area fractions -> consistent two-phase mass flux.
#include "faceAreaFraction.H"

// Foot-point height-function curvature (normal-constant, geometry-only).
#include "footPointCurvature.H"

// Connected zero-surface curvature (one value per interface element, tangential
// regularisation, normal extension directly to CSF faces).
#include "connectedInterfaceCurvature.H"

// Spatially constant, CSF-support-weighted curvature diagnostic.
#include "interfaceMeanCurvature.H"

// Stabilized foot-point re-referencing of the FACE curvature (the balanced-CSF
// delivery measured second-order on the face-centered curvature gate).
#include "stabilizedFootPointFaceCurvature.H"
#include "cellCentreInverseCurvature.H"

// Band renormalization: restore the parallel foliation of psi (psi <- psi/beta_Gamma)
// without moving the zero set. See plan sec. 14.2 for why the operand, not the
// operator, is the measured driver of the curvature-error growth.
#include "bandRenormalization.H"

#include "advectionErrors.H"

// WP0 band mode-spectrum diagnostic (normal-profile 2h/4h/8h amplitudes).
#include "bandModeSpectrum.H"

// WP8.1: across-support vs along-interface split of the delivered kappa_f.
#include "capillaryDriverSplit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A solver for two incompressible, isothermal immiscible fluids using a geometrically maintained semi-Lagrangian level set."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    // Semi-Lagrangian fields FIRST: kappa must be registered before createFields
    // constructs the reconstructedCurvature force (its ctor looks kappa up).
    #include "createSLFields.H"
    #include "createFields.H"
    #include "createTransportFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    // #include "CourantNo.H"
    const volScalarField& alpha = alpha1;
    #include "errorCalculation.H"

    // Parasitic-current + Laplace-jump metrics for the stationary-droplet study.
    #include "createDropletMetricsFile.H"

    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    // Mean curvature on the INITIALISED (t=0) interface: cell-centred, computed
    // symbolically from the quadratic reconstruction with NO foot-point normal
    // extension. Done here so (a) the initial Young-Laplace pressure balance below
    // uses the actual curvature (kappa was zero at construction), and (b) kappa is
    // stored at t=0 for inspection of the curvature on the clean signed-distance
    // field. It is only first-order accurate in the interface cells (the cell
    // centre is offset ~h/2 from the interface) -- that is accepted; no foot-point
    // correction is applied.
    if (fvmCurvature)
    {
        // Force model computes its own curvature; no symbolic fill (see
        // createSLFields.H) -- required for non-quadratic reconstructions.
    }
    else if (closestPointCurvatureExtension)
    {
        slAdv->meanCurvatureClosestPoint(psi, kappa);
    }
    else
    {
        slAdv->meanCurvatureNoExtension(psi, kappa);
    }

    // Same cell-centre inverse at t=0, so the initial Young-Laplace balance uses
    // the curvature the run-time force will use (the defect fixed in 7ad635d).
    if (cellCentreInverseExtension)
    {
        applyCellCentreInverseCurvature
            (
                mesh, slAdv->reconstruction(), kappa,
                cellCentreInverseGaussian
            );
    }

    // Optional harmonic (Laplace) curvature smoothing on the initialised field, so
    // the stored t=0 kappa and the initial Young-Laplace balance use the smoothed
    // field when enabled (no-op otherwise).
    #include "harmonicMeanCurvatureExtensionEqn.H"

    // Optional foot-point height-function curvature at t=0 (same path as the
    // time loop), so the initial Young-Laplace pressure is consistent with the
    // run-time force -- otherwise the first steps start from a polluted balance.
    if (footPointCurvatureExtension)
    {
        computeFootPointCurvature
        (
            mesh, psi,
            mesh.lookupObject<volScalarField>("NarrowBand"),
            kappa
        );
    }

    if (connectedInterfaceCurvatureExtension)
    {
        computeConnectedInterfaceCurvature
        (
            mesh, psi,
            mesh.lookupObject<volScalarField>("NarrowBand"),
            kappa,
            kappaInterfaceFace
        );
    }

    if (interfaceMeanCurvatureExtension)
    {
        applyInterfaceMeanCurvature(mesh, alpha1, kappa);
    }

    // Face delivery LAST (it consumes the final cell kappa), so the initial
    // Young-Laplace balance below uses the same corrected face curvature the
    // run-time force applies.
    if (stabilizedFootPointFaceExtension)
    {
        computeStabilizedFootPointFaceCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace
        );
    }
    else if (cutCellFootPointFaceExtension)
    {
        computeCutCellFootPointFaceCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace
        );
    }
    else if (cellMeanFootPointFaceExtension)
    {
        computeCellMeanFootPointFaceCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace
        );
    }
    else if (symmetricFaceMeanExtension)
    {
        computeSymmetricFaceMeanCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace,
            symmetricFaceMeanTheta
        );
    }
    else if (footPointEvaluatedFaceExtension)
    {
        computeFootPointEvaluatedFaceCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace
        );
    }
    else if (cellFootPointEvaluatedFaceExtension)
    {
        computeCellFootPointEvaluatedFaceCurvature
        (
            mesh, psi, alpha1, kappa,
            slAdv->reconstruction(),
            kappaStableFootFace
        );
    }

    // GUARD (added after footPointEvaluatedFace shipped WITHOUT its branch in
    // the chain above while having one in slAlphaEqn.H): a selected face
    // delivery that never ran here leaves kappaStableFootFace at its
    // construction value of ZERO, so YoungLaplaceEqn below initialises p_rgh
    // with NO capillary force at all and the first momentum step sees the full
    // unbalanced sigma*kappa. The run does not fail -- it silently starts from
    // a droplet with no Laplace jump, and the resulting impulse scales as
    // h^0.5 where a balanced start leaves an O(h^2) residual, which is exactly
    // large enough to corrupt a t_blow comparison between deliveries.
    if
    (
        (stabilizedFootPointFaceExtension || cutCellFootPointFaceExtension
      || cellMeanFootPointFaceExtension || symmetricFaceMeanExtension
      || footPointEvaluatedFaceExtension || cellFootPointEvaluatedFaceExtension)
     && gMax(mag(kappaStableFootFace.primitiveField())) <= VSMALL
    )
    {
        FatalErrorInFunction
            << "curvatureExtension " << curvatureExtensionType
            << " selected, but kappaStableFootFace is identically zero after "
            << "the t=0 face delivery: the initial Young-Laplace balance would "
            << "be solved with no surface tension. The delivery is missing its "
            << "branch in this dispatch chain."
            << exit(FatalError);
    }

    #include "YoungLaplaceEqn.H"
    kappa.write();
    kappaInterfaceFace.write();
    kappaStableFootFace.write();
    p_rgh.write();

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    // TODO(TM): recover the fluid interface on refined mesh.
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        // TODO(TM)
                        // Map the phase indicator.
                        // advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        // Compute cell-centered density.
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // TODO(TM): enable time-step sub-cycling.
            // #include "alphaControls.H"
            // #include "alphaEqnSubcycle.H"
            // Semi-Lagrangian level-set advection + symbolic curvature.
            #include "slAlphaEqn.H"

            // Update viscosity.
            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        // STORMER-VERLET second half drift. With the force centred at n+1/2
        // the level set left slAlphaEqn.H at n+1/2; carry it to n+1 with the
        // velocity the pressure solve just produced, and rebuild the interface
        // so alpha, rho, the curvature and the metrics all describe the end of
        // the step. Own scope: slAlphaEqn.H declares locals, and this is its
        // second textual inclusion in this function.
        if (midpointForceCentring)
        {
            slSecondHalfDrift = true;
            {
                #include "slAlphaEqn.H"
            }
            slSecondHalfDrift = false;
        }

        // The auxiliary rho is used only while pressure and velocity are
        // coupled; restore interface-consistent cell density for output and
        // for the next time level (rhoLENT Algorithm 1, step 13).
        #include "resetRhoLENT.H"

        reportErrors(
            errorFile,
            psi,
            psi0,
            alpha,
            alpha0,
            phi,
            CoNum
        );

        // Append the parasitic-current + Laplace-jump row for this step.
        #include "writeDropletMetrics.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
