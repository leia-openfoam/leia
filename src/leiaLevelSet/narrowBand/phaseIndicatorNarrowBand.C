/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    Part of OpenFOAM; GPL v3 or later. See phaseIndicatorNarrowBand.H.
\*---------------------------------------------------------------------------*/

#include "phaseIndicatorNarrowBand.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvPatch.H"

namespace Foam
{

defineTypeNameAndDebug(phaseIndicatorNarrowBand, false);
addToRunTimeSelectionTable
(
    narrowBand, phaseIndicatorNarrowBand, Dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseIndicatorNarrowBand::phaseIndicatorNarrowBand
(
    const dictionary& dict,
    const volScalarField& psi
)
:
    neighbourNarrowBand(dict, psi),
    alphaName_(dict.getOrDefault<word>("alpha", "alpha")),
    alphaTol_(dict.getOrDefault<scalar>("alphaTol", 1e-8))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void phaseIndicatorNarrowBand::seed()
{
    volScalarField& NarrowBand = field();
    NarrowBand = dimensionedScalar(NarrowBand.dimensions(), 0.0);

    if (!mesh().foundObject<volScalarField>(alphaName_))
    {
        // The indicator is not registered yet. This happens on the very first
        // band calculation of a run, BEFORE the phase indicator has been
        // evaluated for the first time: createFields.H builds the band before
        // it builds alpha. Falling back to the psi sign change for that one
        // call is correct -- psi has just been initialised from the exact
        // implicit surface, so it IS a signed distance at that instant, which
        // is exactly the condition under which the psi-based seed is valid.
        // Every subsequent call finds alpha and uses it.
        WarningInFunction
            << "Phase indicator '" << alphaName_ << "' is not registered;"
            << " seeding this call from the psi sign change instead."
            << " This is expected once, before the first phase-indicator"
            << " evaluation." << endl;

        neighbourNarrowBand::seed();
        return;
    }

    // A COPY, synchronised. Differencing alpha across a coupled face against a
    // STALE far-side value makes the seed depend on the decomposition -- the
    // same defect narrowBand::psiSynced() exists to prevent for psi. alpha is
    // owned by the solver and other users depend on its current boundary
    // state, so it must not be mutated here.
    volScalarField alpha(mesh().lookupObject<volScalarField>(alphaName_));
    alpha.correctBoundaryConditions();

    const auto& own = mesh().owner();
    const auto& nei = mesh().neighbour();

    // (a) INTERNAL FACES: the discrete snGrad(alpha) != 0 test. A face across
    // which the geometric volume fraction jumps is cut by the interface, so
    // both cells sharing it are interface cells.
    forAll(own, faceI)
    {
        if (mag(alpha[own[faceI]] - alpha[nei[faceI]]) > alphaTol_)
        {
            NarrowBand[own[faceI]] = 1.0;
            NarrowBand[nei[faceI]] = 1.0;
        }
    }

    // (b) PARTIALLY FILLED CELLS. A cell clipped at a corner can have a
    // fractional alpha while every face neighbour agrees with it, so the jump
    // test in (a) alone would miss it -- the same class of omission the FIXME
    // in signChangeNarrowBand::calc() records for the psi sign test. This is
    // the direct geometric statement of a cut cell and costs one pass.
    forAll(NarrowBand, cellI)
    {
        const scalar a = alpha[cellI];
        if (a > alphaTol_ && a < 1.0 - alphaTol_)
        {
            NarrowBand[cellI] = 1.0;
        }
    }

    // (c) COUPLED PATCHES: the same jump test against the far-side cell value.
    // patchNeighbourField() returns the NEIGHBOUR CELL CENTRE value, which is
    // what must be differenced here -- not a face value. Writing or reading the
    // wrong quantity on a coupled patch is precisely the defect that was found
    // in velocityModel::setVelocity and that silently corrupted every parallel
    // SDPLS run.
    const auto& alphaBdry = alpha.boundaryField();
    const auto& patches = mesh().boundary();
    const auto& faceOwner = mesh().faceOwner();

    forAll(alphaBdry, patchI)
    {
        const fvPatch& patch = patches[patchI];
        if (!isA<coupledFvPatch>(patch))
        {
            continue;
        }

        const auto aNeiTmp = alphaBdry[patchI].patchNeighbourField();
        const auto& aNei = aNeiTmp();

        forAll(aNei, faceI)
        {
            const label faceJ = faceI + patch.start();
            const label ownCell = faceOwner[faceJ];

            if (mag(alpha[ownCell] - aNei[faceI]) > alphaTol_)
            {
                NarrowBand[ownCell] = 1.0;
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
