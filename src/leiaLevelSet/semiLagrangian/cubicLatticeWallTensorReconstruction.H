/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::cubicLatticeWallTensorReconstruction

Description
    The unit-cube reconstruction (cubicLatticeQuadraticWeightedLeastSquares)
    with the boundary treatment decided on 2026-11-26 (docs/plan-sl-amr-preprint.md
    section 6 item 8; symbol analysis docs/analysis-amplification-sl/sympy/
    wall_face_unknowns.py, candidate faceUnknownsTensorI2, and
    wall_edges_corners.py).

    BOUNDARY POINTS AS ADVECTED UNKNOWNS. A cell that touches a wall carries the
    psi value of every boundary point it owns: one face centre per wall it
    touches, one point on each boundary edge between two of its walls (at the
    height of the cell centre), and the corner if it touches three walls. Each
    value is moved by the semi-Lagrangian step from its own departure point,
    evaluated in the interpolant of its cell (advanceBoundary). Where the
    departure point lies outside the domain (inflow), the value comes from the
    boundary condition of psi at the new time: the patch value for a face; for
    an edge or a corner point the quadratic in the wall plane through the patch
    values and cell values next to it, which is exact for a quadratic field.

    THE WALL CELL. A cell that touches a wall is evaluated with the tensor-product
    quadratic interpolant on 27 points: the nodes {wall, 0, 1} (in cell sizes,
    mirrored at a high wall) along each axis on which it touches a wall, and
    {-1, 0, 1} along the other axes. The points are cell centres and boundary
    points of its own stencil. Every other cell keeps the fit of the base class.

    SCOPE (version 1). The treatment applies to a wall cell that is a cube with
    axes along x, y, z, whose wall faces are normal to those axes, and whose 27
    points exist at its own level. Every other wall cell keeps the base fit and
    is counted in the log. The unstructured formulation is a roadmap item
    (plan section 6 item 9). Static meshes only: after a topology change the
    owner constructs a new object, and the boundary values start again from
    the patch values.

    Requires the direct corrector: advanceBoundary uses the interpolants of
    psi^n that update() built.

    Selected with `reconstruction cubicLatticeWallTensor`.

SourceFiles
    cubicLatticeWallTensorReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef cubicLatticeWallTensorReconstruction_H
#define cubicLatticeWallTensorReconstruction_H

#include "cubicLatticeQuadraticWeightedLeastSquaresReconstruction.H"
#include "FixedList.H"

namespace Foam
{

class cubicLatticeWallTensorReconstruction
:
    public cubicLatticeQuadraticWeightedLeastSquaresReconstruction
{
    //- Kind of a boundary point by its set of wall axes: 0 cell centre,
    //  1..3 face on the wall of axis 0..2, 4..6 edge of the axis pairs
    //  (0 1), (0 2), (1 2), 7 corner.
    static label kindOf(const FixedList<label, 3>& atWall);

    //- Wall side of every cell along each axis: -1 low wall, +1 high wall, 0 none.
    List<FixedList<label, 3>> side_;

    //- True for a cell that uses the tensor interpolant.
    boolList tensor_;

    //- Cell size of every tensor cell.
    scalarField h_;

    //- Own boundary face of a tensor cell on each axis (mesh face index, -1 none).
    List<FixedList<label, 3>> wallFace_;

    //- For a tensor cell, per tensor point (27, index a*9 + b*3 + c over the
    //  nodes of x, y, z): the raw stencil index of the owning cell and the kind.
    List<FixedList<label, 27>> tRaw_;
    List<FixedList<label, 27>> tKind_;

    //- The values of the 27 points of every tensor cell, gathered in update().
    List<FixedList<scalar, 27>> tVal_;

    //- Boundary state: the value of each kind (1..7) owned by a cell.
    PtrList<volScalarField> state_;

    //- The state of the raw stencil entries, collected in update().
    List<List<List<scalar>>> stState_;

    //- psi of the raw stencil entries.
    List<List<scalar>> stPsi_;

    //- Raw stencil centres (built once).
    List<List<vector>> stC_;

    bool wallBuilt_;
    bool stateSet_;
    bool wallReported_;

    //- Counters of the last step.
    label nOutside_, nInflowPoints_, nAdvectedPoints_;

    //- Classify the walls and build the tensor maps. Called once.
    void buildWall();

    //- The node coordinates (cell sizes) along axis a of tensor cell c.
    inline void nodes(const label c, const label a, scalar n[3]) const
    {
        const label s = side_[c][a];
        if (s < 0) { n[0] = -0.5; n[1] = 0; n[2] = 1; }
        else if (s > 0) { n[0] = 0.5; n[1] = 0; n[2] = -1; }
        else { n[0] = -1; n[1] = 0; n[2] = 1; }
    }

    //- Tensor interpolant of cell c at the displacement xi (cell sizes).
    scalar tensorValue(const label c, const vector& xi) const;

    //- Value of a boundary point of kind k (4..7) of cell c from boundary data
    //  at the new time: the quadratic in a wall plane through the face values
    //  fv (per axis, of the cell and its neighbours) and the new cell values.
    scalar inflowPointValue
    (
        const label c,
        const label k,
        const volScalarField& psi,
        const List<FixedList<scalar, 3>>& faceNew
    ) const;

    //- Set the boundary state from the patch values of psi (the first call).
    void initialiseState(const volScalarField& psi);

public:

    TypeName("cubicLatticeWallTensor");

    explicit cubicLatticeWallTensorReconstruction(const fvMesh& mesh);

    virtual ~cubicLatticeWallTensorReconstruction() = default;

    virtual void update(const volScalarField& psiOld);

    virtual scalar evaluateRaw(const label c, const point& x) const;

    virtual void advanceBoundary
    (
        const volVectorField& uNew,
        const volVectorField& uOld,
        const scalar dt,
        volScalarField& psi
    ) override;
};

} // End namespace Foam

#endif
