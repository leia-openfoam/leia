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
	bool correctDefect = true;
	label nCorrectDefect = 0;

        while (correctDefect && simple.correctNonOrthogonal())
	{	
	    auto UEqnPerformance = solve(UEqn == -fvc::grad(p));
	    ++nCorrectDefect;
            correctDefect = 
	    (
	    	(UEqnPerformance.nIterations()[0] > 0) || 
	    	(UEqnPerformance.nIterations()[1] > 0) || 
	    	(UEqnPerformance.nIterations()[2] > 0) 
	    );
	}
        Info << "UEqn defect corrections = " << nCorrectDefect << endl; 

        fvOptions.correct(U);
   }
