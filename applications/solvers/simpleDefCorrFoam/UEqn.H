    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        auto UEqnPerformance = solve(UEqn == -fvc::grad(p));
        Info << "Momentum predictor solved " << endl;
        correctDefect = 
	    (
	    	(UEqnPerformance.nIterations()[0] > 0) || 
	    	(UEqnPerformance.nIterations()[1] > 0) || 
	    	(UEqnPerformance.nIterations()[2] > 0) 
	    );
        fvOptions.correct(U);
    }
