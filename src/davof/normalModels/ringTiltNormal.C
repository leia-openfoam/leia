/*---------------------------------------------------------------------------*\
    leia :: davof :: davofNormalModels :: ringTilt (see ringTiltNormal.H)
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    GNU General Public License v3 or later, as ringTiltNormal.H.
\*---------------------------------------------------------------------------*/

#include "ringTiltNormal.H"
#include "davof.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringTilt, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringTilt, dictionary);
}
}


Foam::davofNormalModels::ringTilt::ringTilt
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    ringArea(dict, scheme),
    tiltC_(dict.get<scalar>("tiltC")),
    tiltR_
    (
        dict.found("tiltRadius")
      ? dict.get<scalar>("tiltRadius")
      : dict.subDict("exactNormal").get<scalar>("radius")
    ),
    tiltDir_(dict.getOrDefault<vector>("tiltDir", vector(0, 0, 1))),
    alternating_(dict.getOrDefault<word>("tiltMode", "biased") == "alternating")
{
    const word mode = dict.getOrDefault<word>("tiltMode", "biased");
    if (mode != "biased" && mode != "alternating")
    {
        FatalIOErrorInFunction(dict)
            << "tiltMode " << mode << ": biased | alternating"
            << exit(FatalIOError);
    }
    tiltDir_ /= mag(tiltDir_);
    Info<< "davof: normalModel ringTilt, tiltC " << tiltC_ << ", tiltRadius "
        << tiltR_ << ", tiltDir " << tiltDir_ << ", tiltMode " << mode
        << " (theta_b = tiltC h / R)" << endl;
}


void Foam::davofNormalModels::ringTilt::modelNormals
(
    const scalar tState,
    vectorField& n
)
{
    ringArea::modelNormals(tState, n);

    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();

    forAll(labels, i)
    {
        const label celli = labels[i];
        const scalar nMag = mag(n[i]);
        if (!(nMag > 0)) continue;
        const vector nu = n[i]/nMag;
        const scalar h = Foam::cbrt(V[celli]);

        scalar s = 1;
        if (alternating_)
        {
            const label par =
                label(std::floor(C[celli].x()/h + 0.5))
              + label(std::floor(C[celli].y()/h + 0.5))
              + label(std::floor(C[celli].z()/h + 0.5));
            s = (mag(par) % 2 == 0) ? 1 : -1;
        }
        const scalar th = s*tiltC_*h/tiltR_;

        vector a = tiltDir_ ^ nu;
        if (mag(a) < 1e-6)
        {
            a = vector(1, 0, 0) ^ nu;
        }
        a /= mag(a);
        n[i] = nMag*(nu*Foam::cos(th) + (a ^ nu)*Foam::sin(th));
    }
}


// ************************************************************************* //
