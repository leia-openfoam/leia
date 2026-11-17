/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::ghostCellQuadraticWeightedLeastSquaresReconstruction

Description
    The anchored quadratic weighted least-squares fit with GHOST CELLS at a 2:1
    refinement transition. In the scheme names of
    docs/analysis-amplification-sl it is C-Cp/L2+g.

    WHY. The fit over the cell-point-cell stencil amplifies at a 2:1 transition,
    because the stencil of a cell there mixes two cell sizes. MEASURED with the
    matrix symbol (docs/analysis-amplification-sl/data/tables,
    scheme_matrix_hex21.csv and scheme_matrix_hex21_ghost.csv): growth factor - 1
    on the planar hexahedral transition is 2.0e-2 at Co 0.1 and 9.6e-2 at Co 0.5
    without ghost cells, and at most 4.4e-3 with them.

    WHAT. Every cell fits a stencil of data at its OWN refinement level.
      restriction   A stencil entry that is a FINER cell is replaced by one datum
                    per refined parent: the mean of the children of that parent,
                    at the parent's centre.
      prolongation  A stencil entry that is a COARSER cell is replaced by ghost
                    data: the fit of that coarse cell, read at the positions of
                    the cell's own lattice that lie in the coarse cell.
    Both are explicit. The cells are fitted level by level, from the coarsest
    level, and the coefficients of a level are exchanged before the next one.

    SCOPE. Hexahedral meshes refined by hexRef8 (refineHexMesh, static), read
    from polyMesh/cellLevel and polyMesh/pointLevel. No structured index is used:
    the parent of a cell is found through the parent's centre point, and the
    lattice of a cell is taken from the edges of its own hexahedron. Without the
    level files, or on a mesh with one level, every cell uses the plain fit over
    its stencil, which is the fit of uncachedQuadraticWeightedLeastSquares
    without its admissibility fallback.

    PARALLEL. The parent sums use syncTools on the parent centre point. The
    levels, parent centres, parent means and fit coefficients of the stencil
    entries come through the map of the stencil, as psi does.

    Selected with `reconstruction ghostCellQuadraticWeightedLeastSquares`.

SourceFiles
    ghostCellQuadraticWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef ghostCellQuadraticWeightedLeastSquaresReconstruction_H
#define ghostCellQuadraticWeightedLeastSquaresReconstruction_H

#include "slReconstruction.H"
#include "volFields.H"

namespace Foam
{

class ghostCellQuadraticWeightedLeastSquaresReconstruction
:
    public slReconstruction
{
    //- Active (non-empty) mesh directions and the coefficient count
    //  (5 in 2D, 9 in 3D).
    labelList activeDirs_;
    label ncoeff_;

    //- Refinement level of every cell, the largest level on all ranks, and the
    //  parent centre point of every cell (-1: no parent found).
    labelList cellLevel_;
    label maxLevel_;
    labelList parentPoint_;

    //- Number of children that share each parent centre point (all ranks).
    scalarField parentCount_;

    //- Lattice vectors of every cell (rows of the tensor); zero if the cell is
    //  not a hexahedron.
    tensorField lattice_;

    //- Fields that are exchanged through the stencil map.
    volScalarField levelFld_;
    volVectorField parentCentreFld_;
    volScalarField parentMeanFld_;
    volTensorField coeffFld_;

    //- The same quantities per stencil entry (raw stencil positions).
    List<List<scalar>> stLevel_;
    List<List<vector>> stParentC_;
    List<List<scalar>> stParentMean_;
    List<List<tensor>> stCoeff_;

    bool built_;

    //- Counters of the last update, for the one-off report.
    label nGhost_, nRestricted_, nLinear_, nConstant_;
    bool reported_;

    //- Raw position in stencil_[c] of the i-th kept entry.
    inline label raw(const label c, const label i) const
    {
        return excludeBoundaryFaces_ ? keep_[c][i] : i;
    }

    //- Basis b(d), without the constant.
    void basis(const vector& d, scalar* b) const;

    //- Value of the polynomial with coefficients t and centre value psiC at
    //  the displacement d from its centre.
    scalar poly(const tensor& t, const scalar psiC, const vector& d) const;

    //- Levels, parents and lattices. Called once.
    void build();

    //- Fit one cell from its data at its own level.
    void fitCell(const label c);

public:

    TypeName("ghostCellQuadraticWeightedLeastSquares");

    explicit ghostCellQuadraticWeightedLeastSquaresReconstruction
    (
        const fvMesh& mesh
    );

    virtual ~ghostCellQuadraticWeightedLeastSquaresReconstruction() = default;

    virtual void update(const volScalarField& psiOld);

    virtual scalar evaluateRaw(const label c, const point& x) const;
};

} // End namespace Foam

#endif
