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
    row("surfaceTensionForce", get(stfD, "type", "<unset>"));
    row("curvatureExtension",  get(ceD, "type", "none (default)"));
    row("geometryFit",         get(gfD, "type", "none (default)"));
    row("psiFilter",           get(pfD, "type", "none (default)"));
    if (get(pfD, "type", "none") != "none")
    {
        row("  psiFilter theta", Foam::name(pfD.getOrDefault<scalar>("theta", 0)));
    }

    Info<< nl << "-- two-phase mass/momentum coupling --" << nl;
    row("massFlux type",       get(mfD, "type", "interpolatedDensity (default)"));
    row("massFlux faceInterp", get(mfD, "faceInterpolation",
                                   get(mfD, "type", "") == "rhoLENT"
                                   ? "central (rhoLENT default)"
                                   : "upwind (default)"));
    row("rho1 / rho2",  Foam::name(r1) + " / " + Foam::name(r2));
    row("DENSITY RATIO", Foam::name(r1/Foam::max(r2, VSMALL)));
    // sigma is read straight from transportProperties: the mixture class
    // exposes sigmaK (sigma*kappa), not sigma itself.
    row("sigma", Foam::name(
        dimensionedScalar("sigma", dimMass/sqr(dimTime),
            mesh.lookupObject<IOdictionary>("transportProperties")).value()));

    Info<< nl << "-- time integration & loop control --" << nl;
    row("momentum ddt", mesh.schemesDict().subDict("ddtSchemes")
                            .getOrDefault<word>("default", "<unset>"));
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
