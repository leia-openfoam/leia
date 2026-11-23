/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slRefinementBand

Description
    The refinement indicator of the AMR line: the registered field `amrBand`,
    1 in the cells within nLayers cell layers of the interface and 0 elsewhere.
    A refining dynamic mesh reads it (blastAMR: errorEstimator fieldValue,
    fieldName amrBand; refine where the value is 1, unrefine where it is 0), so
    that the interface stays inside a uniform band of finest cells.

    Layers are counted through cell-point-cell neighbours on the present mesh,
    with point synchronisation across processor patches. The start layer is
    every cell that has a cell of the other sign of psi among its
    cell-point-cell neighbours, or psi = 0. A layer of a coarser cell counts as
    one layer, so the band is never narrower than nLayers finest cells.

    nLayers is a cell count that follows from the stencils (no tuning):
        3  the face curvature reads psi in the layers l <= 3,
      + 1  the quadratic fit of the last of those cells,
      + 1  so that no band cell reads a ghost value,
      + 1  travel of the interface between two adaptations at CFL <= 1,
      = 6  (levelSet.semiLagrangian.refinementBandLayers, default 6).
    The static curvature sweep of the banded meshes checks this count.

SourceFiles
    slRefinementBand.C

\*---------------------------------------------------------------------------*/

#ifndef slRefinementBand_H
#define slRefinementBand_H

#include "fvMesh.H"
#include "volFields.H"

namespace Foam
{

class slRefinementBand
{
    const fvMesh& mesh_;
    label nLayers_;
    volScalarField band_;

public:

    explicit slRefinementBand(const fvMesh& mesh);

    //- Set the indicator from psi. Call before mesh.update().
    void calc(const volScalarField& psi);

    label nLayers() const { return nLayers_; }

    const volScalarField& band() const { return band_; }
};

} // End namespace Foam

#endif
