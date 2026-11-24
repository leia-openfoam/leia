/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "layerBlendedReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

namespace Foam
{
    defineTypeNameAndDebug(layerBlendedReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        layerBlendedReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerBlendedReconstruction::layerBlendedReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    dict_(slDict_.subDict("layerBlended")),
    inner_(slReconstruction::New(mesh, dict_.get<word>("inner"))),
    outer_(slReconstruction::New(mesh, dict_.get<word>("outer"))),
    bandLayers_(dict_.get<label>("bandLayers")),
    rimLayers_(dict_.get<label>("rimLayers")),
    w_(),
    reported_(false)
{
    if (bandLayers_ < 0 || rimLayers_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "bandLayers must be >= 0 and rimLayers >= 1"
            << exit(FatalIOError);
    }
    if (inner_->type() == typeName || outer_->type() == typeName)
    {
        FatalIOErrorInFunction(dict_)
            << "inner and outer must not be " << typeName << exit(FatalIOError);
    }
    Info<< "layerBlended: inner " << inner_->type() << " for layers <= "
        << bandLayers_ << ", outer " << outer_->type() << " from layer "
        << bandLayers_ + rimLayers_
        << "; the slopeLimiter acts on the outer reconstruction only" << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layerBlendedReconstruction::calcWeights(const volScalarField& psi)
{
    const label nCells = mesh_.nCells();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    // layer 0: a sign change of psi across a face (coupled faces included)
    labelList layer(nCells, labelMax);
    forAll(nei, f)
    {
        if (psi[own[f]]*psi[nei[f]] <= 0)
        {
            layer[own[f]] = 0;
            layer[nei[f]] = 0;
        }
    }
    {
        scalarList psiCells(psi.primitiveField());
        scalarList psiNbr;
        syncTools::swapBoundaryCellList(mesh_, psiCells, psiNbr);
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        for (const polyPatch& pp : pbm)
        {
            if (!pp.coupled()) { continue; }
            const labelUList& fc = pp.faceCells();
            forAll(fc, i)
            {
                const label bf = pp.start() + i - mesh_.nInternalFaces();
                if (psi[fc[i]]*psiNbr[bf] <= 0) { layer[fc[i]] = 0; }
            }
        }
    }

    // layer l + 1: the cells that share a point with a cell of layer l
    const label nMax = bandLayers_ + rimLayers_;
    const labelListList& cp = mesh_.cellPoints();
    labelList mark(mesh_.nPoints());
    for (label l = 1; l <= nMax; ++l)
    {
        mark = 0;
        forAll(layer, c)
        {
            if (layer[c] == l - 1)
            {
                for (const label p : cp[c]) { mark[p] = 1; }
            }
        }
        syncTools::syncPointList(mesh_, mark, maxEqOp<label>(), label(0));
        forAll(layer, c)
        {
            if (layer[c] <= l) { continue; }
            for (const label p : cp[c])
            {
                if (mark[p]) { layer[c] = l; break; }
            }
        }
    }

    w_.setSize(nCells);
    label nIn = 0, nRim = 0;
    forAll(w_, c)
    {
        const label l = layer[c];
        if (l <= bandLayers_) { w_[c] = 1; ++nIn; }
        else if (l < nMax)
        {
            w_[c] = 1.0 - scalar(l - bandLayers_)/scalar(rimLayers_);
            ++nRim;
        }
        else { w_[c] = 0; }
    }

    if (!reported_)
    {
        reduce(nIn, sumOp<label>());
        reduce(nRim, sumOp<label>());
        Info<< "layerBlended: " << nIn << " cells in the H-band, " << nRim
            << " rim cells, "
            << returnReduce(nCells, sumOp<label>()) - nIn - nRim
            << " far-field cells" << endl;
        reported_ = true;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::layerBlendedReconstruction::update(const volScalarField& psiOld)
{
    // The base data that the scheme may read (stencil values, psiOldPtr_).
    collectStencil(psiOld);

    inner_->update(psiOld);
    outer_->update(psiOld);
    calcWeights(psiOld);

    // No limiter of this class: the limiter acts in outer_->evaluate.
    phi_.setSize(mesh_.nCells());
    phi_ = 1.0;
}


Foam::scalar Foam::layerBlendedReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar w = w_[c];
    if (w >= 1) { return inner_->evaluateRaw(c, x); }
    if (w <= 0) { return outer_->evaluate(c, x); }
    return w*inner_->evaluateRaw(c, x) + (1 - w)*outer_->evaluate(c, x);
}

// ************************************************************************* //
