/*---------------------------------------------------------------------------*\
  Manufactured Fourier-mode gate for connected 2-D interface curvature.

  The zero contour is

      r(theta) = R [1 + epsilon cos(m theta)]

  and its exact polar curvature is compared directly with the raw and
  tangentially regularised values attached to the reconstructed interface
  elements.  This tests the production connected-interface implementation,
  not a post-processing replica of it.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "narrowBand.H"
#include "connectedInterfaceCurvature.H"

using namespace Foam;

namespace
{

scalar wrappedPhaseDifference(const scalar a, const scalar b)
{
    scalar d = a - b;
    while (d > constant::mathematical::pi)
    {
        d -= constant::mathematical::twoPi;
    }
    while (d < -constant::mathematical::pi)
    {
        d += constant::mathematical::twoPi;
    }
    return d;
}

struct modeMetrics
{
    scalar mean{0};
    scalar amplitude{0};
    scalar phase{0};
};

modeMetrics projectMode
(
    const scalarField& values,
    const scalarField& theta,
    const scalarField& weights,
    const label mode
)
{
    modeMetrics result;
    scalar sumW = 0;
    forAll(values, i)
    {
        sumW += weights[i];
        result.mean += weights[i]*values[i];
    }
    result.mean /= sumW + VSMALL;

    scalar cosine = 0;
    scalar sine = 0;
    forAll(values, i)
    {
        const scalar fluctuation = values[i] - result.mean;
        cosine += weights[i]*fluctuation*Foam::cos(mode*theta[i]);
        sine += weights[i]*fluctuation*Foam::sin(mode*theta[i]);
    }
    cosine *= 2.0/(sumW + VSMALL);
    sine *= 2.0/(sumW + VSMALL);
    result.amplitude = Foam::sqrt(sqr(cosine) + sqr(sine));
    result.phase = Foam::atan2(sine, cosine);
    return result;
}

void errorNorms
(
    const scalarField& numerical,
    const scalarField& exact,
    const scalarField& weights,
    scalar& l2,
    scalar& linf,
    scalar& relativeL2,
    scalar& relativeLinf
)
{
    scalar sumW = 0;
    scalar sumE2 = 0;
    scalar sumExact2 = 0;
    scalar maxExact = 0;
    linf = 0;
    forAll(numerical, i)
    {
        const scalar error = mag(numerical[i] - exact[i]);
        sumW += weights[i];
        sumE2 += weights[i]*sqr(error);
        sumExact2 += weights[i]*sqr(exact[i]);
        linf = max(linf, error);
        maxExact = max(maxExact, mag(exact[i]));
    }
    l2 = Foam::sqrt(sumE2/(sumW + VSMALL));
    const scalar exactL2 = Foam::sqrt(sumExact2/(sumW + VSMALL));
    relativeL2 = l2/(exactL2 + VSMALL);
    relativeLinf = linf/(maxExact + VSMALL);
}

} // End anonymous namespace


int main(int argc, char* argv[])
{
    argList::addOption("mode", "label", "Manufactured azimuthal mode m.");
    argList::addOption
    (
        "epsilon", "scalar", "Dimensionless radial perturbation (default 0.02)."
    );
    argList::addOption("radius", "scalar", "Mean radius [m] (default 1e-3)."
    );
    argList::addOption
    (
        "center", "vector", "Interface centre (default (0.005 0.005 0))."
    );
    argList::addOption
    (
        "meshVariant", "word", "Mesh label written to the CSV."
    );
    argList::addOption
    (
        "output", "file", "Summary CSV path (default curvatureModeGate.csv)."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label mode = args.getOrDefault<label>("mode", 2);
    const scalar epsilon = args.getOrDefault<scalar>("epsilon", 0.02);
    const scalar radius = args.getOrDefault<scalar>("radius", 1e-3);
    const vector centre = args.getOrDefault<vector>
    (
        "center", vector(0.005, 0.005, 0)
    );
    const word meshVariant =
        args.getOrDefault<word>("meshVariant", "uniform");
    const fileName output =
        args.getOrDefault<fileName>("output", "curvatureModeGate.csv");

    if (mode < 2 || epsilon <= 0 || radius <= 0)
    {
        FatalErrorInFunction
            << "Require mode >= 2, epsilon > 0 and radius > 0."
            << abort(FatalError);
    }

    volScalarField psi
    (
        IOobject
        (
            "psi", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("psi", dimLength, Zero),
        "zeroGradient"
    );

    scalarField& psiInternal = psi.primitiveFieldRef();
    const vectorField& cellCentres = mesh.cellCentres();
    forAll(psiInternal, celli)
    {
        const vector d = cellCentres[celli] - centre;
        const scalar theta = Foam::atan2(d.y(), d.x());
        const scalar interfaceRadius =
            radius*(1.0 + epsilon*Foam::cos(mode*theta));
        psiInternal[celli] = Foam::sqrt(sqr(d.x()) + sqr(d.y()))
                           - interfaceRadius;
    }
    psi.correctBoundaryConditions();

    autoPtr<narrowBand> band = narrowBand::New(mesh, psi);
    band->calc();

    volScalarField kappa
    (
        IOobject
        (
            "kappa", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimless/dimLength, Zero),
        "zeroGradient"
    );
    surfaceScalarField kappaFace
    (
        IOobject
        (
            "kappaInterfaceFace", runTime.timeName(), mesh,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "kappaInterfaceFace", dimless/dimLength, Zero
        )
    );

    connectedInterfaceCurvatureDiagnostics diagnostics;
    computeConnectedInterfaceCurvature
    (
        mesh,
        psi,
        mesh.lookupObject<volScalarField>("NarrowBand"),
        kappa,
        kappaFace,
        &diagnostics
    );

    const bool topologyValid =
        diagnostics.nElements > 0
     && diagnostics.nComponents == 1
     && diagnostics.nOpenEndpoints == 0
     && diagnostics.nFits == diagnostics.nElements
     && diagnostics.nFallbacks == 0;

    const label nElem = diagnostics.nElements;
    scalarField theta(nElem, 0);
    scalarField exact(nElem, 0);
    forAll(theta, i)
    {
        const vector d = diagnostics.elementCentres[i] - centre;
        theta[i] = Foam::atan2(d.y(), d.x());
        const scalar mt = mode*theta[i];
        const scalar r = radius*(1.0 + epsilon*Foam::cos(mt));
        const scalar dr = -radius*epsilon*mode*Foam::sin(mt);
        const scalar d2r =
            -radius*epsilon*sqr(scalar(mode))*Foam::cos(mt);
        exact[i] =
            (sqr(r) + 2.0*sqr(dr) - r*d2r)
           /Foam::pow(sqr(r) + sqr(dr), 1.5);
    }

    const modeMetrics exactMode = projectMode
    (
        exact, theta, diagnostics.elementLengths, mode
    );
    const modeMetrics rawMode = projectMode
    (
        diagnostics.rawCurvature,
        theta,
        diagnostics.elementLengths,
        mode
    );
    const modeMetrics regularizedMode = projectMode
    (
        diagnostics.regularizedCurvature,
        theta,
        diagnostics.elementLengths,
        mode
    );

    scalar rawL2, rawLinf, rawRelL2, rawRelLinf;
    scalar regularizedL2, regularizedLinf, regularizedRelL2, regularizedRelLinf;
    errorNorms
    (
        diagnostics.rawCurvature,
        exact,
        diagnostics.elementLengths,
        rawL2,
        rawLinf,
        rawRelL2,
        rawRelLinf
    );
    errorNorms
    (
        diagnostics.regularizedCurvature,
        exact,
        diagnostics.elementLengths,
        regularizedL2,
        regularizedLinf,
        regularizedRelL2,
        regularizedRelLinf
    );

    const scalar rawTransfer =
        rawMode.amplitude/(exactMode.amplitude + VSMALL);
    const scalar regularizedTransfer =
        regularizedMode.amplitude/(exactMode.amplitude + VSMALL);
    const scalar rawPhaseError = mag
    (
        wrappedPhaseDifference(rawMode.phase, exactMode.phase)
    )*180.0/constant::mathematical::pi;
    const scalar regularizedPhaseError = mag
    (
        wrappedPhaseDifference(regularizedMode.phase, exactMode.phase)
    )*180.0/constant::mathematical::pi;

    const label resolution =
        label(Foam::sqrt(scalar(mesh.nCells())) + 0.5);
    const dictionary& extension =
        static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("curvatureExtension");
    const label fitHalfWidth =
        extension.getOrDefault<label>("fitHalfWidth", 3);
    const word regularizationOperator =
        extension.getOrDefault<word>("regularizationOperator", "helmholtz");
    const scalar regularization =
        extension.getOrDefault<scalar>("tangentialRegularization", 1.0);
    const label regularizationIterations =
        extension.getOrDefault<label>("regularizationIterations", 30);
    const label preserveModesThrough =
        extension.getOrDefault<label>("preserveModesThrough", 4);

    if (Pstream::master())
    {
        OFstream os(output);
        os.precision(12);
        os << "N_CELLS,mesh_variant,mode,epsilon,radius,fit_half_width,"
              "regularization_operator,tangential_regularization,"
              "regularization_iterations,preserve_modes_through,"
              "topology_valid,interface_elements,"
              "interface_components,open_endpoints,curvature_fits,"
              "curvature_fallbacks,exact_mode_amplitude,raw_mode_amplitude,"
              "regularized_mode_amplitude,raw_amplitude_transfer,"
              "regularized_amplitude_transfer,raw_phase_error_deg,"
              "regularized_phase_error_deg,raw_L2,raw_Linf,raw_relative_L2,"
              "raw_relative_Linf,regularized_L2,regularized_Linf,"
              "regularized_relative_L2,regularized_relative_Linf" << nl;
        os << resolution << ',' << meshVariant << ',' << mode << ','
           << epsilon << ',' << radius << ',' << fitHalfWidth << ','
           << regularizationOperator << ',' << regularization << ','
           << regularizationIterations << ',' << preserveModesThrough << ','
           << topologyValid << ','
           << diagnostics.nElements << ',' << diagnostics.nComponents << ','
           << diagnostics.nOpenEndpoints << ',' << diagnostics.nFits << ','
           << diagnostics.nFallbacks << ',' << exactMode.amplitude << ','
           << rawMode.amplitude << ',' << regularizedMode.amplitude << ','
           << rawTransfer << ',' << regularizedTransfer << ','
           << rawPhaseError << ',' << regularizedPhaseError << ','
           << rawL2 << ',' << rawLinf << ',' << rawRelL2 << ','
           << rawRelLinf << ',' << regularizedL2 << ','
           << regularizedLinf << ',' << regularizedRelL2 << ','
           << regularizedRelLinf << nl;

        OFstream elements("curvatureModeElements.csv");
        elements.precision(12);
        elements << "element,x,y,theta,length,kappa_exact,kappa_raw,"
                    "kappa_regularized,fit_ok" << nl;
        forAll(theta, i)
        {
            elements << i << ',' << diagnostics.elementCentres[i].x() << ','
                     << diagnostics.elementCentres[i].y() << ',' << theta[i]
                     << ',' << diagnostics.elementLengths[i] << ',' << exact[i]
                     << ',' << diagnostics.rawCurvature[i] << ','
                     << diagnostics.regularizedCurvature[i] << ','
                     << diagnostics.fitOk[i] << nl;
        }
    }

    Info<< "manufactured curvature mode m=" << mode
        << ": topology=" << (topologyValid ? "valid" : "INVALID")
        << ", raw transfer=" << rawTransfer
        << ", regularized transfer=" << regularizedTransfer
        << ", regularized phase error=" << regularizedPhaseError << " deg"
        << ", regularized relative L2=" << regularizedRelL2 << endl;

    runTime.writeNow();
    return topologyValid ? 0 : 2;
}

// ************************************************************************* //
