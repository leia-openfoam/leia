/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::layerBlendedReconstruction

Description
    Two reconstructions, blended by the cell-layer distance to the interface:
    the INNER one (the quadratic fit) in the H-band, the OUTER one (a
    lower-order, bounded fit) in the far field, and a linear change of the
    weight over the rim layers between them.

        value(c, x) = w_c inner.evaluateRaw(c, x) + (1 - w_c) outer.evaluate(c, x)

    H-BAND. The band of cells that the curvature from the cell-centred
    quadratic fit needs. The layer l of a cell is counted in cells: l = 0 for a
    cell with a sign change of psi across one of its faces, and layer l + 1 is
    every cell that shares a POINT with a cell of layer l (the cell-point-cell
    stencil of the fit, so the fit of a cell of layer l reads layers <= l + 1).
        w = 1                       for l <= bandLayers
        w = 1 - (l - bandLayers)/rimLayers   in the rim
        w = 0                       for l >= bandLayers + rimLayers
    The two counts are cell counts, the same for every case. No physical length
    enters.

    LIMITER. The key slopeLimiter of the semiLagrangian dictionary acts on the
    OUTER reconstruction only (outer.evaluate). The inner reconstruction is
    read without a limiter (inner.evaluateRaw), and this class applies no
    limiter of its own.

    GEOMETRY. fitDerivatives and the curvature functions are those of the
    inner reconstruction.

    STATE 2026-11-24: a research candidate (supervisor: quadratic fit in the
    H-band, another scheme outside). Both reconstructions are updated in every
    cell, so this version saves no cost. HYPOTHESIS, not measured: the values
    that enter the H-band through the rim keep the H-band at the accuracy of
    the inner fit.

Usage
    \verbatim
    semiLagrangian
    {
        reconstruction  layerBlended;
        layerBlended
        {
            inner       cubicLatticeQuadraticWeightedLeastSquares;
            outer       linearWeightedLeastSquares;
            bandLayers  4;
            rimLayers   2;
        }
    }
    \endverbatim

SourceFiles
    layerBlendedReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef layerBlendedReconstruction_H
#define layerBlendedReconstruction_H

#include "slReconstruction.H"

namespace Foam
{

class layerBlendedReconstruction
:
    public slReconstruction
{
    const dictionary& dict_;
    autoPtr<slReconstruction> inner_;
    autoPtr<slReconstruction> outer_;
    label bandLayers_;
    label rimLayers_;

    //- Weight of the inner reconstruction in every cell (last update).
    scalarField w_;

    bool reported_;

    //- Fill w_ from the sign changes of psi.
    void calcWeights(const volScalarField& psi);

public:

    TypeName("layerBlended");

    explicit layerBlendedReconstruction(const fvMesh& mesh);

    virtual ~layerBlendedReconstruction() = default;

    virtual void update(const volScalarField& psiOld);

    virtual scalar evaluateRaw(const label c, const point& x) const;

    //- Weight of the inner reconstruction (for diagnostics).
    const scalarField& weights() const { return w_; }

    // Geometry: the inner reconstruction.

    virtual bool providesFitDerivatives() const override
    {
        return inner_->providesFitDerivatives();
    }
    virtual label fitDerivatives
    (
        const label c, vector& g, symmTensor& H
    ) const override
    {
        return inner_->fitDerivatives(c, g, H);
    }
    virtual void meanCurvature(volScalarField& k) const override
    {
        inner_->meanCurvature(k);
    }
    virtual void meanCurvatureLaplacian(volScalarField& k) const override
    {
        inner_->meanCurvatureLaplacian(k);
    }
    virtual void meanCurvatureNoExtension(volScalarField& k) const override
    {
        inner_->meanCurvatureNoExtension(k);
    }
    virtual void meanCurvatureClosestPoint(volScalarField& k) const override
    {
        inner_->meanCurvatureClosestPoint(k);
    }
    virtual scalar signedOffset(const label c, bool& fb) const override
    {
        return inner_->signedOffset(c, fb);
    }
    virtual scalar footPointDistance
    (
        const label c, const point& p, const scalar level, bool& ok
    ) const override
    {
        return inner_->footPointDistance(c, p, level, ok);
    }
    virtual scalar curvatureAtFootPoint
    (
        const label c, const point& p, bool& ok
    ) const override
    {
        return inner_->curvatureAtFootPoint(c, p, ok);
    }
};

} // End namespace Foam

#endif
