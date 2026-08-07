/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt 
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

\*---------------------------------------------------------------------------*/


#include "gradPsiErrorCSV.H"
#include "addToRunTimeSelectionTable.H"

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(gradPsiErrorCSV, 0);
    addToRunTimeSelectionTable(functionObject, gradPsiErrorCSV, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::gradPsiErrorCSV::gradPsiErrorCSV
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject (name, runTime, dict),
    field_(mesh_.lookupObject<volScalarField>(dict.get<word>("field"))),
    psi_(mesh_.lookupObject<volScalarField>("psi")),
    // gradPsi_(mesh_.lookupObject<volVectorField>("gradPsi")),// Does not work, but why?
    csvFile_(fileName("gradPsiError.csv"),  IOstreamOption(), IOstreamOption::appendType::APPEND)

{
    if ( !fileSize("gradPsiError.csv") && Pstream::myProcNo() == 0)
    {
            // CSV Header 
            csvFile_ << "TIME,"
                << "MAX_DELTA_X,"
                << "MEAN_DELTA_X,"
                << "E_MAX_GRAD_PSI,"
                << "E_MEAN_GRAD_PSI,"
                << "E_NARROW_MAX_GRAD_PSI,"
                << "E_NARROW_MEAN_GRAD_PSI,"
                << "E_L1_GRAD_PSI,"
                << "E_L2_GRAD_PSI,"
                << "E_NARROW_L1_GRAD_PSI,"
                << "E_NARROW_L2_GRAD_PSI,"
                << "MAX_MAG_GRAD_PSI,"
                << "MEAN_MAG_GRAD_PSI,"
                << "NARROW_MAX_MAG_GRAD_PSI,"
                << "NARROW_MEAN_MAG_GRAD_PSI,"
                // Conditioning diagnostic: the interface position error scales
                // like err_psi / |grad psi|, so a FLATTENING level set
                // (min |grad psi| -> 0 near the interface) amplifies psi noise
                // into O(1) position noise.
                << "MIN_MAG_GRAD_PSI,"
                << "NARROW_MIN_MAG_GRAD_PSI\n";
    }
    write();
}

bool Foam::functionObjects::gradPsiErrorCSV::write()
{

    const fvMesh& mesh = this->mesh_;
    const Time& runTime = mesh.time();

    const auto deltaX = pow(mesh.deltaCoeffs(),-1).cref();
    const auto max_deltaX = gMax(deltaX);
    const auto mean_deltaX = gAverage(deltaX);
    const auto gradPsiError = field_;
    const auto max_gradPsiError = gMax(gradPsiError);
    const auto mean_gradPsiError = gAverage(gradPsiError);

    // Empty-band sentinel. An empty narrow band means the phase has been
    // annihilated (no psi sign change anywhere) -- it does NOT mean the error
    // is zero. Reporting 0.0 made a dead interface read as a perfect one, and
    // the downstream ">0" filters in the plotting/order-fitting scripts then
    // dropped those rows silently, so a fitted order could be computed over an
    // unannounced subset of the ladder. NaN is both honest in the CSV and still
    // excluded by those same filters (NaN > 0 is false).
    const scalar undefined = std::numeric_limits<scalar>::quiet_NaN();

    scalar max_narrowGradPsiError = undefined;
    scalar mean_narrowGradPsiError = undefined;

    // Volume-weighted L1 and L2 (RMS) norms of e = ||grad psi| - 1|:
    //   L1 = sum(e V)/sum(V),  L2 = sqrt(sum(e^2 V)/sum(V)).
    const scalarField& cellV = mesh.V().field();
    const scalarField& eField = field_.primitiveField();
    const scalar sumV = gSum(cellV);
    const scalar L1_gradPsiError =
        gSum(eField*cellV)/Foam::max(sumV, VSMALL);
    const scalar L2_gradPsiError =
        Foam::sqrt(gSum(Foam::sqr(eField)*cellV)/Foam::max(sumV, VSMALL));
    scalar L1_narrowGradPsiError = undefined;
    scalar L2_narrowGradPsiError = undefined;

    // Same dedicated, unlimited keyword the gradPsiError field object uses --
    // the two must measure with the SAME operator, and it must not be the
    // advection scheme's limited `grad(psi)` (see gradPsiError.C).
    const volScalarField magGradPsi = mag(fvc::grad(psi_, "gradPsiMetric"));
    const auto max_magGradPsi = gMax(magGradPsi);
    const auto mean_magGradPsi = gAverage(magGradPsi);
    const auto min_magGradPsi = gMin(magGradPsi);
    scalar max_narrowMagGradPsi = undefined;
    scalar mean_narrowMagGradPsi = undefined;
    scalar min_narrowMagGradPsi = undefined;

    if (mesh.objectRegistry::found("NarrowBand"))
    {
        const auto narrowBand = mesh.lookupObject<volScalarField>("NarrowBand");
        List<scalar> narrowGradPsiError = subset(narrowBand, gradPsiError);
        // The band is EMPTY once the phase is fully annihilated (no psi sign
        // change): g{Max,Min} on an empty list return -/+VGREAT garbage.
        // Leave the NaN sentinel in place instead -- see `undefined` above.
        const label nBand =
            returnReduce(narrowGradPsiError.size(), sumOp<label>());
        if (nBand > 0)
        {
            max_narrowGradPsiError = gMax(narrowGradPsiError);
            mean_narrowGradPsiError = gAverage(narrowGradPsiError);

            List<scalar> narrowMagGradPsi = subset(narrowBand, magGradPsi);
            max_narrowMagGradPsi = gMax(narrowMagGradPsi);
            mean_narrowMagGradPsi = gAverage(narrowMagGradPsi);
            min_narrowMagGradPsi = gMin(narrowMagGradPsi);
        }

        // Volume-weighted L1/L2 of the error restricted to the narrow band.
        const scalarField& nbField = narrowBand.primitiveField();
        scalarField nbMask(eField.size(), 0.0);
        forAll(nbMask, c)
        {
            nbMask[c] = (nbField[c] > 0.5) ? 1.0 : 0.0;
        }
        const scalar nbV = gSum(nbMask*cellV);
        if (nbV > VSMALL)
        {
            L1_narrowGradPsiError = gSum(nbMask*eField*cellV)/nbV;
            L2_narrowGradPsiError =
                Foam::sqrt(gSum(nbMask*Foam::sqr(eField)*cellV)/nbV);
        }
    }

    if (Pstream::myProcNo() == 0)
    {
        csvFile_ << runTime.timeOutputValue() << ","
            << max_deltaX << ","
            << mean_deltaX << ","
            << max_gradPsiError << ","
            << mean_gradPsiError << ","
            << max_narrowGradPsiError << ","
            << mean_narrowGradPsiError << ","
            << L1_gradPsiError << ","
            << L2_gradPsiError << ","
            << L1_narrowGradPsiError << ","
            << L2_narrowGradPsiError << ","
            << max_magGradPsi << ","
            << mean_magGradPsi << ","
            << max_narrowMagGradPsi << ","
            << mean_narrowMagGradPsi << ","
            << min_magGradPsi << ","
            << min_narrowMagGradPsi << "\n";
    }
    return true;
}



// ************************************************************************* //
