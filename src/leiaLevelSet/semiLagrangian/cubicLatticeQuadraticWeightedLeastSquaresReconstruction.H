/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::cubicLatticeQuadraticWeightedLeastSquaresReconstruction

Description
    The anchored quadratic weighted least-squares fit with ghost cells at a 2:1
    refinement transition, for CUBIC cells, with ONE fit operator for all
    regular cells. It is the scheme of ghostCellQuadraticWeightedLeastSquares
    (C-Cp/L2+g); that class stays unchanged and is the reference.

    THE UNIT OPERATOR. A cubic cell whose 26 lattice positions all hold one
    datum (a cell centre of the own level, the mean of the children of a finer
    parent, or a ghost value read from the fit of a coarser cell) has a regular
    stencil. The weights are 1/distance^2, so the weighted least-squares
    solution in lattice coordinates xi = d . inv(lattice) is the same linear
    map for every such cell: a 9 x 26 matrix, computed once for the unit cube.
    The physical gradient and Hessian follow from the lattice of the cell:
        g = inv(L) . g_xi,     H = inv(L) . H_xi . inv(L)^T .
    This holds for a cube of any size and orientation, also at a transition.

    NOT REGULAR, fitted per cell as in the reference class: a cell with a
    physical patch face in its stencil (the datum is at the face centre, half a
    cell away), a cell with a ghost position that no coarser stencil cell
    contains, a cell that is not a cube, and every cell of a 2D case.

    LEVELS. The refinement levels are taken from the registered objects
    cellLevel and pointLevel (a dynamic mesh that refines keeps them in the
    registry and does not write them between write times); if they are not
    registered they are read from polyMesh, as refineHexMesh writes them. The
    object is built once for one mesh: after a topology change the owner
    constructs a new object.

    Also provides fitDerivatives (the gradient and the Hessian at the centre).

    Selected with `reconstruction cubicLatticeQuadraticWeightedLeastSquares`.

SourceFiles
    cubicLatticeQuadraticWeightedLeastSquaresReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef cubicLatticeQuadraticWeightedLeastSquaresReconstruction_H
#define cubicLatticeQuadraticWeightedLeastSquaresReconstruction_H

#include "slReconstruction.H"
#include "volFields.H"

namespace Foam
{

class cubicLatticeQuadraticWeightedLeastSquaresReconstruction
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

    //- The unit operator: 9 coefficients (in lattice coordinates) from the 26
    //  lattice data, row-major 9 x 26. Built in the constructor (3D only).
    scalar P_[9*26];
    bool haveP_;

    //- True if the cell is a cube (three equal, orthogonal lattice vectors).
    boolList cube_;

    //- For a cell whose stencil is exactly the 26 cells of its own level:
    //  the lattice position (0..25) of the kept stencil entries 1..n-1.
    //  Empty for every other cell.
    List<labelList> slot_;

    //- Fit order of every cell after the last update (2, 1 or 0).
    labelList order_;

    //- Counters of the last update: cells fitted with the unit operator
    //  through slot_, through the matched data, and per cell.
    label nUnitSlot_, nUnitMatched_, nPerCell_;

    //- Cells without a lattice that keep coarser cells in their stencil.
    label nMixed_;

    //- Lattice position 0..25 of the lattice coordinates xi, or -1 if xi is
    //  not within 0.05 of a lattice position other than the centre.
    static label latticeSlot(const vector& xi);

    //- Build P_.
    void buildUnitOperator();

    //- Store the coefficients of a regular cell from its 26 lattice data
    //  (value minus centre value, in lattice order).
    void fitUnit(const label c, const scalar* rl);

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

    TypeName("cubicLatticeQuadraticWeightedLeastSquares");

    explicit cubicLatticeQuadraticWeightedLeastSquaresReconstruction
    (
        const fvMesh& mesh
    );

    virtual ~cubicLatticeQuadraticWeightedLeastSquaresReconstruction() = default;

    virtual void update(const volScalarField& psiOld);

    virtual scalar evaluateRaw(const label c, const point& x) const;

    //- Gradient and Hessian of the fit at the cell centre; returns the order.
    virtual label fitDerivatives(const label c, vector& g, symmTensor& H) const;
};

} // End namespace Foam

#endif
