/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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

#include "redistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvcGrad.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(redistancer, false);
defineRunTimeSelectionTable(redistancer, Mesh);
addToRunTimeSelectionTable(redistancer, redistancer, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

redistancer::redistancer(const fvMesh& mesh)
    :
        fvSolution_(mesh),
        levelSetDict_(fvSolution_.subDict("levelSet")),
        redistDict_(levelSetDict_.subDict("redistancer")),
        redistanceInterval_(redistDict_.getOrDefault<label>("redistanceInterval", 1)),
        trigger_(redistDict_.getOrDefault<word>("trigger", "interval")),
        threshold_(redistDict_.getOrDefault<scalar>("threshold", -1.0)),
        gradBound_(redistDict_.getOrDefault<scalar>("gradBound", 2.0)),
        bandViolationFraction_(
            redistDict_.getOrDefault<scalar>("bandViolationFraction", 0.0)),
        criterionRegion_(
            redistDict_.getOrDefault<word>("criterionRegion", "band")),
        anchorLayers_(redistDict_.getOrDefault<label>("anchorLayers", 2))
{
    if (trigger_ != "interval" && trigger_ != "gradPsiThreshold"
     && trigger_ != "signedDistanceBounds")
    {
        FatalIOErrorInFunction(redistDict_)
            << "redistancer trigger must be 'interval', 'gradPsiThreshold' "
            << "or 'signedDistanceBounds', got '" << trigger_ << "'"
            << exit(FatalIOError);
    }
    if (criterionRegion_ != "band" && criterionRegion_ != "bulk")
    {
        FatalIOErrorInFunction(redistDict_)
            << "redistancer criterionRegion must be 'band' or 'bulk', got '"
            << criterionRegion_ << "'" << exit(FatalIOError);
    }
    if (criterionRegion_ == "bulk" && bandViolationFraction_ <= 0)
    {
        // A signed distance has |grad psi| = 0 on its skeleton (medial
        // axis): the strict min/max rule over the bulk fires forever there.
        FatalIOErrorInFunction(redistDict_)
            << "criterionRegion bulk requires bandViolationFraction > 0 "
            << "(the SDF skeleton has |grad psi| = 0 by construction; the "
            << "strict min/max rule cannot terminate there)."
            << exit(FatalIOError);
    }
    if (gradBound_ <= 1.0)
    {
        FatalIOErrorInFunction(redistDict_)
            << "gradBound must be > 1 (band |grad psi| bounds are "
            << "[1/gradBound, gradBound]), got " << gradBound_
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<Foam::redistancer> redistancer::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution (mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& redistDict = levelSetDict.subDict("redistancer");
    const word& modelType = redistDict.getOrDefault<word>("type", "noRedistancing");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = MeshConstructorTable(modelType);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "redistancer",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<redistancer>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar redistancer::resolveThreshold(const fvMesh& mesh) const
{
    if (threshold_ < 0)
    {
        // (h/L)^2: the paper's Delta^2 with Delta the nondimensional grid
        // spacing -- h = smallest cell-centre spacing, L = largest extent
        // of the mesh bounding box.
        const scalar h =
            1.0/gMax(mesh.deltaCoeffs().primitiveField());
        const scalar L = cmptMax(mesh.bounds().span());
        threshold_ = sqr(h/L);
        Info<< "redistancer: gradPsiThreshold auto threshold (h/L)^2 = "
            << threshold_ << " (h = " << h << ", L = " << L << ")" << endl;
    }
    return threshold_;
}


scalar redistancer::gradPsiBandL1(const volScalarField& psi) const
{
    const fvMesh& mesh = psi.mesh();

    const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());
    const scalar bandWidth = 6.0*h;

    const volScalarField magGradPsi(mag(fvc::grad(psi)));
    const scalarField& psiIn = psi.primitiveField();
    const scalarField& gIn = magGradPsi.primitiveField();
    const scalarField& V = mesh.V();

    scalar sumDev = 0;
    scalar sumV = 0;
    forAll(psiIn, celli)
    {
        if (mag(psiIn[celli]) < bandWidth)
        {
            sumDev += V[celli]*mag(gIn[celli] - 1.0);
            sumV += V[celli];
        }
    }
    // Global reductions: the fire/no-fire decision must be identical on
    // every rank (the anchor construction below contains matched MPI calls).
    reduce(sumDev, sumOp<scalar>());
    reduce(sumV, sumOp<scalar>());

    return (sumV > VSMALL) ? sumDev/sumV : 0.0;
}


void redistancer::topologicalMask
(
    const volScalarField& psi,
    const label nLayers,
    boolList& mask
) const
{
    const fvMesh& mesh = psi.mesh();

    if (!mesh.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "the topological anchor mask requires the registered "
            << "'NarrowBand' field: set levelSet.narrowBand.type to "
            << "signChange in fvSolution." << exit(FatalError);
    }
    const volScalarField& nb =
        mesh.lookupObject<volScalarField>("NarrowBand");

    // Marker field so each dilation sweep sees processor-halo neighbours
    // (constraint patches get their constraint type; physical patches
    // zeroGradient). A hole in a frozen guard layer at a processor boundary
    // would let the fill overwrite a cell that a serial run keeps.
    volScalarField m
    (
        IOobject
        (
            "redistAnchorMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        word("zeroGradient")
    );
    forAll(nb, celli)
    {
        m[celli] = (nb[celli] == 1) ? 1.0 : 0.0;
    }

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    for (label layer = 0; layer < nLayers; ++layer)
    {
        m.correctBoundaryConditions();
        scalarField add(mesh.nCells(), 0);
        forAll(own, facei)
        {
            if (m[own[facei]] > 0.5 || m[nei[facei]] > 0.5)
            {
                add[own[facei]] = 1;
                add[nei[facei]] = 1;
            }
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!p.coupled()) continue;
            const scalarField nbr
            (
                m.boundaryField()[patchi].patchNeighbourField()
            );
            const labelUList& fc = p.faceCells();
            forAll(fc, i)
            {
                if (nbr[i] > 0.5) add[fc[i]] = 1;
            }
        }
        forAll(add, celli)
        {
            if (add[celli] > 0.5) m[celli] = 1.0;
        }
    }

    mask.setSize(mesh.nCells());
    forAll(mask, celli)
    {
        mask[celli] = (m[celli] > 0.5);
    }
}


label redistancer::gradPsiBandBounds
(
    const volScalarField& psi,
    scalar& gMin,
    scalar& gMax,
    scalar& violFrac
) const
{
    const fvMesh& mesh = psi.mesh();

    if (!mesh.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "the signedDistanceBounds redistancing criterion requires "
            << "the registered 'NarrowBand' field: set "
            << "levelSet.narrowBand.type to signChange in fvSolution."
            << exit(FatalError);
    }
    const volScalarField& nb =
        mesh.lookupObject<volScalarField>("NarrowBand");

    boolList inBand(mesh.nCells(), false);
    if (criterionRegion_ == "bulk")
    {
        // Complement of the frozen anchor region + one junction layer:
        // anchors are untouched by construction and the anchor/fill junction
        // kink is a property of the Dirichlet coupling, not of degradation --
        // measuring either would retrigger forever.
        boolList anchorPlusJunction;
        topologicalMask(psi, anchorLayers_ + 1, anchorPlusJunction);
        forAll(inBand, celli)
        {
            inBand[celli] = !anchorPlusJunction[celli];
        }
    }
    else
    {
        // Topological band: sign-change cells + one dilation layer. (The
        // dilation is rank-local; a missing cross-processor neighbour cell
        // cannot change the O(1) firing behaviour of a global reduction.)
        const labelListList& cellCells = mesh.cellCells();
        forAll(nb, celli)
        {
            if (nb[celli] == 1)
            {
                inBand[celli] = true;
                for (const label nbrI : cellCells[celli])
                {
                    inBand[nbrI] = true;
                }
            }
        }
    }

    // Dedicated criterion gradient (fvSchemes "redistGradPsi"), decoupled
    // from the psi advection scheme's limiter.
    const volVectorField gradPsi(fvc::grad(psi, "redistGradPsi"));

    gMin = GREAT;
    gMax = -GREAT;
    label nBand = 0, nViol = 0;
    forAll(inBand, celli)
    {
        if (!inBand[celli]) continue;
        const scalar g = mag(gradPsi[celli]);
        gMin = min(gMin, g);
        gMax = max(gMax, g);
        ++nBand;
        if (g < 1.0/gradBound_ || g > gradBound_) ++nViol;
    }
    reduce(gMin, minOp<scalar>());
    reduce(gMax, maxOp<scalar>());
    reduce(nBand, sumOp<label>());
    reduce(nViol, sumOp<label>());

    violFrac = (nBand > 0) ? scalar(nViol)/nBand : 0.0;
    return nBand;
}


bool redistancer::redistance(volScalarField& psi)
{
    if (trigger_ == "signedDistanceBounds")
    {
        scalar gMin, gMax, violFrac;
        const label nBand = gradPsiBandBounds(psi, gMin, gMax, violFrac);

        if (nBand == 0)
        {
            WarningInFunction
                << "empty sign-change band: the interface is lost -- "
                << "nothing to redistance from (the criterion should have "
                << "fired earlier; check gradBound)." << endl;
            return false;
        }

        // Reported EVERY step: each run self-documents its margin to the
        // usability bounds.
        Info<< "redistancer: band |grad psi| in [" << gMin << ", " << gMax
            << "], bounds [" << 1.0/gradBound_ << ", " << gradBound_
            << "], violating fraction " << violFrac << endl;

        const bool fire =
            (bandViolationFraction_ > 0)
          ? (violFrac > bandViolationFraction_)
          : (gMin < 1.0/gradBound_ || gMax > gradBound_);

        if (!fire || !doRedistance(psi))
        {
            return false;
        }

        gradPsiBandBounds(psi, gMin, gMax, violFrac);
        Info<< "redistancer: EVENT (signedDistanceBounds) at t = "
            << psi.mesh().time().timeOutputValue()
            << ", post-event band |grad psi| in [" << gMin << ", " << gMax
            << "]" << endl;
        return true;
    }

    if (trigger_ == "interval")
    {
        if (psi.mesh().time().timeIndex() % redistanceInterval_ == 0)
        {
            return doRedistance(psi);
        }
        return false;
    }

    // trigger == "gradPsiThreshold" (Abu-Al-Saud et al. 2018, Eq. 12)
    const scalar threshold = resolveThreshold(psi.mesh());
    const scalar L1 = gradPsiBandL1(psi);

    // Reported EVERY evaluation: the achievable post-event floor of this
    // criterion (set by the anchor residual) must be calibrated against the
    // threshold from the run log.
    Info<< "redistancer: band L1(||grad psi|-1|) = " << L1
        << " (threshold " << threshold << ")" << endl;

    if (L1 <= threshold)
    {
        return false;
    }

    if (!doRedistance(psi))
    {
        return false;
    }

    const scalar L1post = gradPsiBandL1(psi);
    Info<< "redistancer: EVENT at t = "
        << psi.mesh().time().timeOutputValue()
        << ", band L1 " << L1 << " -> " << L1post << endl;

    if (L1post > L1)
    {
        // An event that fires must not injure: a fill that RAISES the band
        // deviation feeds back through this criterion (fires -> worse ->
        // fires again) and compounds. Loud, per-event, in the log.
        WarningInFunction
            << "redistance event INCREASED the band deviation ("
            << L1 << " -> " << L1post
            << "): the selected fill is injuring the field." << endl;
    }

    return true;
}

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
