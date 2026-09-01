/*---------------------------------------------------------------------------*\
    ONE CONSOLIDATED STARTUP BANNER: every runtime-selected algorithm and the
    numbers that decide the discretisation, printed from the dictionaries the
    solver ACTUALLY READS -- never from a study token or a config file.

    WHY THIS EXISTS. The method is assembled from a dozen independent runtime
    selections (transport, curvature, delivery, mass flux, filters, loop
    control), and until now each printed at most a line of its own, several
    printed nothing, and no single place said what the run was. Three defects
    in this campaign would have been visible immediately in a banner like this:

      * five cases whose fvSolution had NO traceVelocity entry, so every arm of
        a trace sweep silently ran the solver default and the arms differed in
        nothing;
      * a DUPLICATE faceInterpolation key, where OpenFOAM took the last one and
        the run used `upwind` while the materialized dictionary said `central`
        -- the configuration recording one thing and the solver doing another;
      * PHASE_INDICATOR left unpinned, so convergence orders were about to be
        fitted across a mixture of two interface reconstructions.

    Read it as the run's own account of itself: if a line here disagrees with
    the study config, the LINE is what happened.
\*---------------------------------------------------------------------------*/

{
    const fvSolution& fvSol = static_cast<const fvSolution&>(mesh);
    const dictionary& lsDict = fvSol.subDict("levelSet");

    auto sub = [&](const word& k) -> dictionary
    {
        return lsDict.subOrEmptyDict(k);
    };
    auto get = [&](const dictionary& d, const word& k, const word& dflt) -> word
    {
        return d.getOrDefault<word>(k, dflt);
    };
    // Pad by hand: <iomanip> is not pulled in here and Info is an OSstream,
    // not a std::ostream, so std::setw would not apply to it anyway.
    auto row = [&](const char* k, const std::string& v)
    {
        std::string key(k);
        if (key.size() < 26) key.append(26 - key.size(), ' ');
        Info<< "  " << key.c_str() << v.c_str() << nl;
    };

    const dictionary slD   = sub("semiLagrangian");
    const dictionary mfD   = sub("massFlux");
    const dictionary veD   = sub("velocityExtension");
    const dictionary nbD   = sub("narrowBand");
    const dictionary piD   = sub("phaseIndicator");
    const dictionary stfD  = sub("surfaceTensionForce");
    const dictionary ceD   = sub("curvatureExtension");
    const dictionary gfD   = sub("geometryFit");
    const dictionary pfD   = sub("psiFilter");
    const dictionary& pimple = fvSol.subDict("PIMPLE");

    const dimensionedScalar& rho1d = mixture.rho1();
    const dimensionedScalar& rho2d = mixture.rho2();
    const scalar r1 = rho1d.value(), r2 = rho2d.value();

    Info<< nl
        << "===================== LEIA METHOD CONFIGURATION ====================" << nl
        << "  (printed from the dictionaries in use; if this disagrees with the" << nl
        << "   study config, THIS is what the run did)" << nl << nl
        << "-- transport (level set) --" << nl;
    row("advection",        get(sub("advection"), "type", "semiLagrangian"));
    row("SL reconstruction", get(slD, "reconstruction", "<unset>"));
    row("SL correction",     get(slD, "correction", "<unset>"));
    row("SL footIntegrator", get(slD, "footIntegrator", "taylor (default)"));
    row("SL traceVelocity",  get(slD, "traceVelocity", "cellCentred (default)"));
    row("SL CFLmax",         slD.getOrDefault<scalar>("CFLmax", -1) < 0
                             ? std::string("<unset>")
                             : Foam::name(slD.get<scalar>("CFLmax")));
    row("velocityExtension", get(veD, "type", "none (default)"));
    row("narrowBand",        get(nbD, "type", "<unset>"));
    row("phaseIndicator",    get(piD, "type", "<unset>"));

    Info<< nl << "-- interface force / curvature --" << nl;
    const word stfType = get(stfD, "type", "<unset>");
    row("surfaceTensionForce", stfType);
    // The FACE curvature the CSF force actually applies is a separate runtime
    // choice from the cell model, and it was invisible here until 2026-09-01:
    // `model` is the force's own assembly (plain interpolation of the cell field),
    // `registered` takes a solver-registered face field, named on the next row.
    row("  faceCurvatureSource", get(stfD, "faceCurvatureSource", "model"));
    if (get(stfD, "faceCurvatureSource", "model") == "registered")
    {
        row("  faceCurvatureField",
            get(stfD, "faceCurvatureField", "<unset>"));
    }
    // A prescribed-curvature run is the well-balanced gate, not production: say so
    // loudly, because every error norm in it is measured against an exact kappa.
    if (stfType.find("onstantCurvature") != std::string::npos)
    {
        row("  PRESCRIBED curvature",
            Foam::name(stfD.getOrDefault<scalar>("curvature", 0)));
    }
    row("curvatureExtension",  get(ceD, "type", "none (default)"));
    if (get(ceD, "type", "none") == "cellCentreInverse")
    {
        row("  Gaussian-curvature K",
            ceD.getOrDefault<bool>("gaussianCurvature", true)
              ? "on (K-aware, 3D-active)" : "off (K == 0, the 2D algebra)");
    }
    row("geometryFit",         get(gfD, "type", "none (default)"));
    row("psiFilter",           get(pfD, "type", "none (default)"));
    if (get(pfD, "type", "none") != "none")
    {
        row("  psiFilter theta", Foam::name(pfD.getOrDefault<scalar>("theta", 0)));
    }

    Info<< nl << "-- two-phase mass/momentum coupling --" << nl;
    row("massFlux type",       get(mfD, "type", "interpolatedDensity (default)"));
    // Which reconstructed plane gives the face area fraction -- NOT an
    // interpolation scheme, and unrelated to div(rhoPhi,U).
    row("massFlux alpha_f source",
        get(mfD, "alphaFSource",
            get(mfD, "faceInterpolation", "donorPlane (default)")));
    row("rho1 / rho2",  Foam::name(r1) + " / " + Foam::name(r2));
    row("DENSITY RATIO", Foam::name(r1/Foam::max(r2, VSMALL)));
    // sigma is read straight from transportProperties: the mixture class
    // exposes sigmaK (sigma*kappa), not sigma itself.
    row("sigma", Foam::name(
        dimensionedScalar("sigma", dimMass/sqr(dimTime),
            mesh.lookupObject<IOdictionary>("transportProperties")).value()));

    // div(rhoPhi,U): the momentum face interpolation, distinct from the alpha_f
    // source above and repeatedly confused with it. A scheme entry is a TOKEN
    // STREAM, not a word -- "Gauss upwind" is two tokens, and reading it with
    // getOrDefault<word> is a FOAM FATAL IO ERROR ("1 excess tokens in stream"),
    // which is exactly how this line first shipped. Print the tokens.
    {
        const dictionary& divD = mesh.schemesDict().subDict("divSchemes");
        std::string divRhoPhiU("<not set>");
        if (divD.found("div(rhoPhi,U)"))
        {
            OStringStream buf;
            const ITstream& is = divD.lookup("div(rhoPhi,U)");
            forAll(is, i)
            {
                if (i) buf << ' ';
                buf << is[i];
            }
            divRhoPhiU = buf.str();
        }
        row("div(rhoPhi,U)", divRhoPhiU);
    }

    Info<< nl << "-- time integration & loop control --" << nl;
    const dictionary& ddtD = mesh.schemesDict().subDict("ddtSchemes");
    row("momentum ddt", ddtD.getOrDefault<word>("default", "<unset>"));
    // ddt(rho) is a SEPARATE key from ddt(rho,U) and is read only by the rhoLENT
    // auxiliary density equation. No case wrote it until 2026-09-01, so it silently
    // inherited `default` -- BDF2 on an equation whose boundedness argument assumes
    // Euler. Print it whenever it is actually consumed.
    if (get(mfD, "type", "geometricFaceDensity") == "rhoLENT")
    {
        row("ddt(rho)  [rhoLENT]",
            ddtD.getOrDefault<word>
            (
                "ddt(rho)",
                ddtD.getOrDefault<word>("default", "<unset>")
            ));
    }
    row("deltaT",       Foam::name(runTime.deltaTValue()));
    row("adjustTimeStep", runTime.controlDict()
                            .getOrDefault<bool>("adjustTimeStep", false)
                            ? "yes" : "no");
    row("nOuterCorrectors",
        Foam::name(pimple.getOrDefault<label>("nOuterCorrectors", 1)));
    row("nCorrectors",
        Foam::name(pimple.getOrDefault<label>("nCorrectors", 1)));
    row("psiOuterCorrectors",
        pimple.getOrDefault<Switch>("psiOuterCorrectors", false) ? "yes" : "no");
    row("residualControl",
        pimple.found("residualControl") ? "PRESENT" : "absent");

    Info<< nl << "-- mesh --" << nl;
    row("cells (global)", Foam::name(returnReduce(mesh.nCells(), sumOp<label>())));
    row("MPI ranks",      Foam::name(Pstream::nProcs()));
    Info<< "===================================================================" << nl
        << endl;
}
