/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slRefinementBand.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"

Foam::slRefinementBand::slRefinementBand(const fvMesh& mesh)
:
    mesh_(mesh),
    nLayers_
    (
        mesh.solutionDict().subOrEmptyDict("levelSet")
            .subOrEmptyDict("semiLagrangian")
            .getOrDefault<label>("refinementBandLayers", 6)
    ),
    band_
    (
        IOobject
        (
            "amrBand",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{
    Info<< "slRefinementBand: " << nLayers_
        << " cell layers on each side of the interface" << endl;
}


void Foam::slRefinementBand::calc(const volScalarField& psi)
{
    const labelListList& cp = mesh_.cellPoints();
    const label nPoints = mesh_.nPoints();
    scalarField& b = band_.primitiveFieldRef();

    // -- start layer: a sign change among the cell-point-cell neighbours
    scalarField pMin(nPoints, GREAT), pMax(nPoints, -GREAT);
    forAll(cp, c)
    {
        forAll(cp[c], i)
        {
            const label p = cp[c][i];
            pMin[p] = min(pMin[p], psi[c]);
            pMax[p] = max(pMax[p], psi[c]);
        }
    }
    syncTools::syncPointList(mesh_, pMin, minEqOp<scalar>(), scalar(GREAT));
    syncTools::syncPointList(mesh_, pMax, maxEqOp<scalar>(), scalar(-GREAT));
    forAll(cp, c)
    {
        scalar lo = GREAT, hi = -GREAT;
        forAll(cp[c], i)
        {
            lo = min(lo, pMin[cp[c][i]]);
            hi = max(hi, pMax[cp[c][i]]);
        }
        b[c] = (lo <= 0 && hi >= 0) ? 1.0 : 0.0;
    }

    // -- the further layers
    scalarField pMark(nPoints);
    for (label layer = 1; layer < nLayers_; ++layer)
    {
        pMark = 0;
        forAll(cp, c)
        {
            if (b[c] > 0.5)
            {
                forAll(cp[c], i) { pMark[cp[c][i]] = 1.0; }
            }
        }
        syncTools::syncPointList(mesh_, pMark, maxEqOp<scalar>(), scalar(0));
        forAll(cp, c)
        {
            if (b[c] > 0.5) { continue; }
            forAll(cp[c], i)
            {
                if (pMark[cp[c][i]] > 0.5) { b[c] = 1.0; break; }
            }
        }
    }
    band_.correctBoundaryConditions();
}
