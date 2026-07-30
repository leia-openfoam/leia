/*---------------------------------------------------------------------------*\
  Diagnostic interface-mean curvature delivery.

  Replace the reconstructed cell curvature by one spatially constant value,
  weighted by the actual CSF face support |snGrad(alpha)| |Sf|.  This is not a
  production curvature model: it is the controlled null-space experiment for
  the balanced-force chain.  With constant sigma it makes

      Gsigma = sigma*kappaBar*snGrad(alpha)*|Sf|

  an exact discrete pressure-gradient flux while retaining the mean value of
  the curvature computed from the transported interface.
\*---------------------------------------------------------------------------*/

#ifndef interfaceMeanCurvature_H
#define interfaceMeanCurvature_H

namespace Foam
{

inline scalar applyInterfaceMeanCurvature
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    volScalarField& kappa
)
{
    const surfaceScalarField support
    (
        "interfaceMeanCurvatureSupport",
        mag(fvc::snGrad(alpha))*mesh.magSf()
    );
    const surfaceScalarField kappaf
    (
        "interfaceMeanCurvatureFaceValue",
        fvc::interpolate(kappa)
    );

    const scalarField& w = support.primitiveField();
    const scalarField& kf = kappaf.primitiveField();
    scalar sumW = 0.0;
    scalar sumWKappa = 0.0;
    forAll(w, facei)
    {
        sumW += w[facei];
        sumWKappa += w[facei]*kf[facei];
    }
    reduce(sumW, sumOp<scalar>());
    reduce(sumWKappa, sumOp<scalar>());

    if (sumW <= VSMALL)
    {
        FatalErrorInFunction
            << "Cannot form interface-mean curvature: the internal-face "
            << "CSF support is zero." << abort(FatalError);
    }

    const scalar kappaBar = sumWKappa/sumW;
    kappa.primitiveFieldRef() = kappaBar;
    kappa.correctBoundaryConditions();

    Info<< "Interface-mean curvature diagnostic: kappaBar="
        << kappaBar << " 1/m, CSF support=" << sumW << " m" << endl;

    return kappaBar;
}

} // End namespace Foam

#endif

