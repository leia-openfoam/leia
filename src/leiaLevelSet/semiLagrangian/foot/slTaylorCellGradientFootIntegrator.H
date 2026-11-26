/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slTaylorCellGradientFootIntegrator

Description
    The taylor foot of the base class with the velocity gradient of its dt^2
    term from the CELL values only: the weighted least-squares gradient over the
    face neighbours of the cell (weights 1/|d|^2), with the neighbour cell of a
    processor face, and without the physical boundary faces. Exact for a linear
    velocity, also in a wall cell.

    Reason (2026-11-26, wallTensorRotationDt3D): the Gauss gradient reads the
    patch values of the trace velocity reconstruct(phi), which are the adjacent
    cell values, so it is wrong by O(|grad u|) in every wall cell; in the dt^2
    term that is an O(dt) error after a fixed time. HYPOTHESIS until measured.

    Selected with footIntegrator taylorCellGradient.

SourceFiles
    slTaylorCellGradientFootIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef slTaylorCellGradientFootIntegrator_H
#define slTaylorCellGradientFootIntegrator_H

#include "slFootIntegrator.H"

namespace Foam
{

class slTaylorCellGradientFootIntegrator
:
    public slFootIntegrator
{
public:

    TypeName("taylorCellGradient");

    explicit slTaylorCellGradientFootIntegrator(const fvMesh& mesh);

    virtual ~slTaylorCellGradientFootIntegrator() = default;

    //- The least-squares cell gradient of u (rows: d/dx_a), from cell values.
    static tmp<tensorField> cellGradient(const volVectorField& u);

    virtual void feet
    (
        const volVectorField& uNew,
        const volVectorField& uOld,
        const scalar dt,
        pointField& feet
    ) const;
};

} // End namespace Foam

#endif
