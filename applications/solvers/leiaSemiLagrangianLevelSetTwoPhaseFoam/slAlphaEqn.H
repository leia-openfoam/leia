{
    // DEFAULT (psiOuterCorrectors no): advect the level set ONCE per time step,
    // on the first outer corrector: the semi-Lagrangian update mutates psi in
    // place (psi^n -> psi^{n+1}), so running it again on later PIMPLE iterations
    // would compound. The interface is then held fixed while the pressure-
    // velocity system is corrected (as interFoam holds alpha fixed across the
    // inner correctors) -- which FREEZES the capillary force at its once-per-
    // step AB2-predicted value against the implicit-Euler momentum LHS.
    //
    // ITERATED COUPLING (psiOuterCorrectors yes): every outer corrector
    // restores psi^n (from the psiN snapshot taken below) and re-advects it
    // with the CURRENT outer iterate of the velocity, then re-runs the whole
    // band -> indicator -> curvature -> face-delivery pipeline, so the force
    // the next UEqn sees converges to its t^{n+1} value. At outer-loop
    // convergence the step is backward-Euler consistent in ALL terms; the AB2
    // estimate is demoted to the pass-1 initial guess.
    if (pimple.firstIter() || psiOuterCorrectors || slSecondHalfDrift)
    {
        // Semi-Lagrangian advection: psi holds psi^n on entry, psi^{n+1} on exit.
        //
        // TIME LEVELS -- the coupled solver has only OLD-TIME velocity here.
        // slAdvection::advect(psi, Unew, Uold) documents Unew = u^{n+1}, and the
        // pointValue kernel therefore builds the arrival-time-centred foot
        //   x_d = x_c - Unew dt + 1/2 [ (Unew-Uold)/dt + (Unew.grad)Unew ] dt^2.
        // But this include runs BEFORE UEqn.H/pEqn.H, so the registered U is still
        // u^n (the velocity converged at the end of the previous step) and UextOld
        // is u^{n-1}. Feeding (u^n, u^{n-1}) into slots that require (u^{n+1}, u^n)
        // leaves an O(dt^2 d_t u) foot error -- the arrival-centred expansion needs
        // the NEW-time velocity in its linear term, and substituting the old one
        // adds +dt^2 d_t u that the temporal correction then reinforces instead of
        // cancelling. Per-step position error O(dt^2) instead of O(dt^3), i.e. the
        // transport degrades to FIRST order in time. It is invisible in every
        // steady or uniform test (d_t u = 0 makes the two forms identical) and only
        // bites once the interface velocity is genuinely unsteady.
        //
        // The departure-time-centred trace built from the two old levels is second
        // order and needs no u^{n+1} (the method is explicit in the current and the
        // previous velocity field by design):
        //   x_d = x_c - dt/2 (3 u^n - u^{n-1}) + 1/2 dt^2 (u^n.grad)u^n
        //       = x_c - dt u^{n+1/2} + 1/2 dt^2 (u.grad)u,
        // i.e. the LINEAR displacement must use the time-centred velocity, obtained
        // by Adams-Bashforth extrapolation from the two available levels; only the
        // temporal-derivative term changes sign relative to the arrival form, the
        // convective term keeps its plus sign under both centrings.
        //
        // Rather than branch the shared kernel (the prescribed-velocity solvers
        // genuinely do supply u^{n+1} after ++runTime, and the + sign is correct
        // for them), we hand the kernel the AB2 estimate of the arrival velocity,
        //   UextStar = u^n + (dt/dt0)(u^n - u^{n-1}) = u^{n+1} + O(dt^2),
        // and u^n as the old level. The unchanged kernel then reduces algebraically
        // to the expression above, with the dt/dt0 ratio making the temporal
        // difference quotient correct when the step size varies (the kernel divides
        // by the CURRENT step while the difference spans the PREVIOUS one). On the
        // first step UextOld == Uext, so the extrapolation is the identity -- the
        // usual one-off AB2 startup, inside the second-order budget.
        //
        // No mesh-motion flux subtraction: SL runs on a static mesh (cached
        // stencil).
        //
        // Diagnostic hook (default off): setting the environment variable
        // SL_FREEZE_INTERFACE=1 skips the advection so psi is held at its initial
        // value. On the stationary droplet this pins psi to the clean signed
        // distance and isolates the parasitic current driven by the fixed
        // curvature error from the transport-coupled feedback -- the controlled
        // experiment behind workflow/scripts/droplet_mechanism_probe.sh. Unset in
        // every production run, so the transport is byte-for-byte unchanged.
        // Compute an interface-normal-constant velocity from the CURRENT
        // physical U and pre-advection psi. The extension is used only in the
        // characteristic update; momentum continues to use U and phi. For
        // velocityExtension type none this is exactly the identity path.
        slVelExt->correct();

        if (slSecondHalfDrift)
        {
            // STORMER-VERLET, second half drift. psi is at n+1/2 and momentum
            // has just produced u^{n+1}; carry it the remaining dt/2 so the step
            // ends with psi^{n+1}, and so alpha, rho and the metrics below are at
            // the step end rather than at the midpoint. slVelExt->correct() above
            // re-read the registered U, which pEqn left at u^{n+1}, so the (new,
            // old) slots are (u^{n+1}, u^n) -- no extrapolation is needed here.
            if (Foam::getEnv("SL_FREEZE_INTERFACE").empty())
            {
                slAdv->setDtScale(0.5);
                slAdv->advect(psi, slVelExt->Uext(), UextOld);
                slAdv->setDtScale(1.0);
            }
        }
        else if (pimple.firstIter())
        {
            // Snapshot psi^n BEFORE the in-place advection: later outer passes
            // (psiOuterCorrectors) re-advect from this level.
            psiN == psi;

            if (Foam::getEnv("SL_FREEZE_INTERFACE").empty())
            {
                // dt/dt0: deltaT0 equals deltaT on a fixed step and is guarded
                // here because it is not guaranteed to be set on the very
                // first step.
                const scalar dtRatio =
                    runTime.deltaT0Value() > SMALL
                  ? runTime.deltaTValue()/runTime.deltaT0Value()
                  : 1.0;

                // Trace-velocity source: the extension field (default) or
                // the divergence-free reconstruct(phi) (Basilisk property).
                // Refresh the flux-based levels for this step before use.
                if (traceFromFlux)
                {
                    UtracePtr() == slTraceField();
                }
                const volVectorField& UtrNew =
                    traceFromFlux ? UtracePtr() : slVelExt->Uext();
                const volVectorField& UtrOld =
                    traceFromFlux ? UtraceOldPtr() : UextOld;

                const volVectorField UextStar
                (
                    "UextStar",
                    UtrNew + dtRatio*(UtrNew - UtrOld)
                );

                // midpoint centring drifts only HALF the step here; the other
                // half runs after the pressure solve. endStep keeps the full
                // step, so this line is byte-identical with the switch off.
                slAdv->setDtScale(midpointForceCentring ? 0.5 : 1.0);
                slAdv->advect(psi, UextStar, UtrNew);
                slAdv->setDtScale(1.0);
            }
            // Preserve the pre-momentum extension as the old velocity level.
            // Under psiOuterCorrectors this is also the OLD slot of every
            // later pass this step: after this line UextOld holds ext(u^n).
            // NOT updated on passes >= 2 -- the old level must stay u^n.
            UextOld == slVelExt->Uext();
            if (traceFromFlux)
            {
                // Same contract as UextOld: the OLD slot holds the t^n level
                // for the whole step and is never refreshed on passes >= 2.
                UtraceOldPtr() == UtracePtr();
            }
        }
        else
        {
            // psiOuterCorrectors pass >= 2: restore psi^n and re-advect it by
            // the CURRENT outer iterate of the velocity. slVelExt->correct()
            // above re-read the registered U, which the previous pass's
            // momentum+pressure solve left at u^{n+1,(k-1)} -- so the kernel's
            // documented (u^{n+1}, u^n) slots are finally both genuine: NO
            // extrapolation, arrival-centred foot with the actual iterate.
            psi == psiN;

            if (Foam::getEnv("SL_FREEZE_INTERFACE").empty())
            {
                if (traceFromFlux)
                {
                    // phi now carries the PROJECTED flux of the previous outer
                    // iterate, so this is the divergence-free u^{n+1,(k-1)}.
                    UtracePtr() == slTraceField();
                    slAdv->advect(psi, UtracePtr(), UtraceOldPtr());
                }
                else
                {
                    slAdv->advect(psi, slVelExt->Uext(), UextOld);
                }
            }
        }

        // Narrow band FIRST, then the phase indicator: the geometric indicator
        // reads the registered NarrowBand field to decide which cells get the
        // geometric volume fraction, so the band must be current with the freshly
        // advected psi (previously it lagged one step, giving cells that had just
        // entered the band a sign-based 0/1 alpha for one step).
        narrowBand->calc();

        // Optional grid-scale psi filter (default off; psiFilterEqn.H) --
        // the explicit delta-psi damping term of the stability inequality.
        #include "psiFilterEqn.H"

        // Optional redistancing (fvSolution levelSet.redistancer; default
        // noRedistancing = the original reinit-free pipeline, byte-identical).
        // The plane-anchored models (anchoredEikonal, planeFootWave) realign
        // psi to the SAME linear LSQ interface planes the phase indicator and
        // the geometric mass flux read, so redistancing is consistency with
        // the discrete interface, not displacement of it. With trigger
        // gradPsiThreshold events are rare (Abu-Al-Saud et al. 2018, Eq. 12) --
        // both papers show every-step reinit destroys the capillary
        // equilibrium. Anchors need the CURRENT band (hence after calc()); an
        // event may realign a cut-cell centre sub-cell, so the band is
        // refreshed afterwards, keeping the indicator consistent.
        if (redist->redistance(psi))
        {
            narrowBand->calc();
        }

        // Band renormalization (default off): psi <- psi/beta_Gamma on the band,
        // restoring the parallel foliation the parallel-surface inverse assumes.
        // Placed BEFORE the phase indicator and the curvature fill so both read
        // the renormalized field; the zero set is unchanged by construction, so
        // the band computed above stays valid and needs no recalc. The fit is
        // refreshed first (meanCurvatureNoExtension updates the reconstruction on
        // the current psi; the kappa it writes is overwritten by the fill below).
        if (bandRenormalization)
        {
            slAdv->meanCurvatureNoExtension(psi, kappa);
            renormalizeBandPsi
            (
                mesh, psi,
                mesh.lookupObject<volScalarField>("NarrowBand"),
                slAdv->reconstruction(),
                psiRenormThreshold,
                psiRenormSweeps,
                psiRenormRelaxation
            );
        }

        // PHASE-VOLUME CORRECTION (rhoLENT paper eqs. 46-49). Before the phase
        // indicator, because alpha1 -- and through it alpha2, the mixture, rho and
        // mu in UEqn/pEqn, the geometric face density and the CSV row -- must all be
        // built from the CORRECTED field. eps is spatially constant, so
        // grad(psi + eps) = grad(psi) and nothing the redistancer established is
        // undone. Default noVolumeCorrection returns false without touching psi.
        if (volCorr->correct(psi, phaseInd()))
        {
            // The zero crossing moved sub-h: the sign-change band can gain or lose
            // a cell, exactly as after a redistance event.
            narrowBand->calc();
        }

        phaseInd->calcPhaseIndicator(alpha1, psi);

        // Previous-step curvature, captured BEFORE any refill: needed by the
        // optional temporal under-relaxation of the foot-point path below.
        // NOTE under psiOuterCorrectors: on passes >= 2 this holds the previous
        // PASS's kappa, not the previous time step's -- the footPointCurvature
        // relax < 1 diagnostic then blends against the current-step iterate,
        // which is the consistent choice inside an iterated coupling.
        const scalarField kappaPrevStep(kappa.primitiveField());

        // Freeze-kappa diagnostic: skip every refill after the first step so the
        // force acts with the step-1 curvature for the whole run.
        if (!(freezeKappa && runTime.timeIndex() > 1))
        {

        // Mean curvature evaluated SYMBOLICALLY from the quadratic reconstruction
        // (fit gradient + Hessian): default at the CELL CENTRE (no normal
        // extension); with curvatureExtension type closestPointNewton the
        // evaluation point is the TRUE closest point of the cell's own
        // quadratic (trust-region KKT Newton -- NOT the gradient-projection
        // flow of the retired meanCurvature extension), giving the
        // normal-constant delivery in the force-support layer with per-cell
        // fallback to the centre value.
        if (fvmCurvature)
        {
            // The force model computes its own curvature (FVM operators);
            // the symbolic fill would FatalError on non-quadratic
            // reconstructions (linearTaylor). kappa stays 0.
        }
        else if (closestPointCurvatureExtension)
        {
            slAdv->meanCurvatureClosestPoint(psi, kappa);
        }
        else
        {
            slAdv->meanCurvatureNoExtension(psi, kappa);
        }

        // CELL-CENTRE parallel-surface inverse, in place on the cell field, so
        // everything downstream (interpolation to faces, the metrics, the
        // driver split) sees an interface-referenced curvature without any
        // face-side conversion.
        if (cellCentreInverseExtension)
        {
            applyCellCentreInverseCurvature
            (
                mesh, slAdv->reconstruction(), kappa,
                cellCentreInverseGaussian
            );
        }

        // Publish the fit normal alongside the curvature it belongs to: the
        // normalised gradient of the SAME per-cell quadratic, over the narrow
        // band (zero elsewhere). Consumed by name (`nHatFit`) by the
        // semiImplicitCapillaryForce fvOption; cost is one fitDerivatives
        // sweep of the band, negligible against the fit itself. Filled inside
        // the freeze-kappa guard so a frozen curvature keeps its matching
        // frozen normal.
        {
            const slReconstruction& siRecon = slAdv->reconstruction();
            const volScalarField& siBand =
                mesh.lookupObject<volScalarField>("NarrowBand");
            nHatFit == dimensionedVector("nHatFit", dimless, vector::zero);
            if (siRecon.type() == "uncachedQuadraticWeightedLeastSquares")
            {
                vectorField& nRef = nHatFit.primitiveFieldRef();
                vector g(Zero);
                symmTensor H(Zero);
                forAll(nRef, celli)
                {
                    if (siBand[celli] == 0) continue;
                    if
                    (
                        siRecon.fitDerivatives(celli, g, H) >= 1
                     && mag(g) > SMALL
                    )
                    {
                        nRef[celli] = g/mag(g);
                    }
                }
                nHatFit.correctBoundaryConditions();
            }
        }

        // Optional harmonic (Laplace) smoothing of the cell-centred curvature,
        // pinning the interface cells (no foot-point extension). No-op unless
        // enabled in fvSolution (levelSet.curvatureExtension.type harmonicLaplace).
        #include "harmonicMeanCurvatureExtensionEqn.H"

        // Optional foot-point height-function curvature: overwrites the band
        // (+2 rings) with the normal-constant height-function kappa; deficient
        // interface cells keep the symbolic value above (footPointCurvature.H).
        // Optional temporal under-relaxation (curvatureExtension.relax in (0,1]):
        // kappa^n = relax*kappa_new + (1-relax)*kappa^{n-1} -- filters the
        // per-step discrete switching of feet/parents (diagnostic for the
        // coupled instability trigger; 1 = no relaxation).
        if (footPointCurvatureExtension)
        {
            computeFootPointCurvature
            (
                mesh, psi,
                mesh.lookupObject<volScalarField>("NarrowBand"),
                kappa
            );
            const scalar krlx =
                static_cast<const fvSolution&>(mesh).subDict("levelSet")
                .subOrEmptyDict("curvatureExtension")
                .getOrDefault<scalar>("relax", 1.0);
            if (krlx < 1.0)
            {
                // TRUE temporal under-relaxation: blend against the previous
                // TIME STEP's delivered curvature (captured before the symbolic
                // refill above), low-pass filtering the per-step discrete
                // switching of feet/parents.
                kappa.primitiveFieldRef() =
                    krlx*kappa.primitiveField() + (1.0 - krlx)*kappaPrevStep;
                kappa.correctBoundaryConditions();
            }
        }

        // Connected reconstruction candidate: replace the independent cell
        // fits by one conforming zero-surface chain, regularise curvature only
        // along that chain, attach one value to each segment, and deliver it
        // directly along normals to the face-native CSF field.
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

        // Controlled pressure-null-space experiment: use the curvature computed
        // above only through its CSF-support-weighted interface mean.  The force
        // remains on the production Gsigma face-flux path; only the spatial
        // curvature fluctuations are removed.
        if (interfaceMeanCurvatureExtension)
        {
            applyInterfaceMeanCurvature(mesh, alpha1, kappa);
        }

        // Optional DIAGNOSTIC clamp (levelSet.curvatureClamp {min; max;}): caps the
        // grid-scale sign-indefinite kappa excursions to a physical range to test
        // whether they are load-bearing for the parasitic runaway. Off by default.
        if (kappaClampOn)
        {
            kappa.primitiveFieldRef() =
                min(max(kappa.primitiveField(), kappaClampMin), kappaClampMax);
            kappa.correctBoundaryConditions();
        }
        }   // end of the freeze-kappa guard around the curvature refill

        // Stabilized foot-point FACE delivery: re-reference the interpolated
        // face curvature to the interface (parallel-curve inverse at the
        // stabilized foot-point offset of each active face centre). Runs LAST
        // so it consumes the final cell kappa; the reconstruction's fits are
        // current (the symbolic fill above updated them for this psi).
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

        // Diagnostic hook (default off): SL_SMOOTH_RHO=1 replaces the sharp
        // geometric indicator by the erf-smoothed Heaviside of the level set,
        // H = 0.5*(1 - erf(psi/eps)) with eps = 1.5*h, exactly the mixture rule
        // of Abu-Al-Saud, Popinet & Tchelepi (2018, Eqs. 10-11), BEFORE the
        // mixture properties are formed. Probe for the CST force model: their
        // MAC reference never exposes pure light-phase density at faces that
        // receive the off-diagonal (corner) stress divergence one row away
        // from the alpha jump; our sharp indicator does (rho_f = rho_air at
        // those faces, 839x the blended inertia). Unset in production.
        if (!Foam::getEnv("SL_SMOOTH_RHO").empty())
        {
            const scalar h =
                1.0/gMax(mesh.deltaCoeffs().primitiveField());
            const scalar eps = 1.5*h;
            scalarField& aIn = alpha1.primitiveFieldRef();
            const scalarField& psiIn = psi.primitiveField();
            forAll(aIn, celli)
            {
                aIn[celli] = 0.5*(1.0 - Foam::erf(psiIn[celli]/eps));
            }
            alpha1.correctBoundaryConditions();
        }

        alpha2 = 1.0 - alpha1;
        mixture.correct();

        // Geometrical face density at the NEW interface location. In rhoLENT
        // this remains fixed during the outer iterations while phi^o changes;
        // the same rhof*phi^o then appears in the auxiliary density equation
        // and in div(rhoPhi,U). The older comparison modes also reuse rhof.
        if (!freezeRhoPhi && massFluxType != "interpolatedDensity")
        {
            const volScalarField& nb =
                mesh.lookupObject<volScalarField>("NarrowBand");
            if (alphaFAdvected)
            {
                // The face velocity that traces the plane is the SAME field the SL
                // scheme traces the cell feet with -- reconstruct(phi) for the
                // projectedFlux setting -- interpolated to faces. Using phi's own
                // reconstruction keeps the face and cell displacements from being two
                // different velocity fields, which is the whole point.
                const surfaceVectorField Uf
                (
                    "UfSLtrace", fvc::interpolate(fvc::reconstruct(phi))
                );
                computeFaceAreaFractions
                (
                    mesh, psi.oldTime(), nb, phi, alphaf, massFluxUpwind,
                    &Uf, runTime.deltaTValue()
                );
            }
            else
            {
                computeFaceAreaFractions
                (
                    mesh, psi, nb, phi, alphaf, massFluxUpwind
                );
            }
            if (alphaFTrapezoid)
            {
                // TIME-AVERAGED face density (author's proposal, 2026-09-02). The
                // density that crosses a face during the step is the average over the
                // SWEPT region, not the plane value at t^{n+1}; the trapezoid of the
                // two endpoint reconstructions is its second-order approximation.
                //
                // psi.oldTime() is the interface at t^n and is maintained by the
                // semi-Lagrangian scheme (pointValueScheme uses it directly). The
                // donor direction is taken from the CURRENT phi for both levels: phi
                // changes by O(dt) across the step, so using it for the old level is
                // a second-order-consistent choice and keeps the two reconstructions
                // on the same upwind side, which a mixed choice would not.
                computeFaceAreaFractions
                (
                    mesh, psi.oldTime(), nb, phi, alphafOld, massFluxUpwind
                );
                const surfaceScalarField alphafBar
                (
                    "alphafBar", 0.5*(alphaf + alphafOld)
                );
                rhof = alphafBar*rho1 + (scalar(1) - alphafBar)*rho2;
            }
            else
            {
                rhof = alphaf*rho1 + (scalar(1) - alphaf)*rho2;
            }
        }

        // The non-rhoLENT modes use the geometrical cell density immediately.
        // rhoLENT deliberately keeps rho.oldTime() as the previous geometrical
        // state and computes the current auxiliary rho below.
        if (!freezeRhoPhi && !useRhoLENT)
        {
            rho == alpha1*rho1 + alpha2*rho2;
            rho.correctBoundaryConditions();
        }
    }

    // Mass flux and auxiliary density are updated on EVERY PIMPLE outer
    // iteration. This is load-bearing: pEqn changes phi^o, so rhoPhi and the
    // rhoLENT continuity update must be reassembled before the next UEqn.
    if (!freezeRhoPhi)
    {
        if (massFluxType == "interpolatedDensity")
        {
            rhoPhi = fvc::interpolate(rho)*phi;
        }
        else
        {
            rhoPhi = rhof*phi;
        }

        if (projectMassFlux)
        {
            // Correct the flux BEFORE anything reads it -- the auxiliary density
            // equation below and div(rhoPhi,U) must both see the corrected flux, or
            // the projection fixes a balance nobody uses.
            #include "massFluxProjection.H"
        }

        if (useRhoLENT)
        {
            #include "rhoLENTEqn.H"
        }
        else if (massResidualDiagnostic)
        {
            // Same residual, for a model that solves no density equation. Pure
            // diagnostic; inert unless massResidualDiagnostic is switched on.
            #include "massResidualDiag.H"
        }
    }
}

scalar domainFraction = 0;
domainFraction = alpha1.weightedAverage(mesh.Vsc()).value();

Info<< "Phase-1 volume fraction = "
    << domainFraction
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
