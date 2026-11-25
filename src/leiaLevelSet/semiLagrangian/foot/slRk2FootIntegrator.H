/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slRk2FootIntegrator

Description
    Midpoint RK2 foot with velocity SAMPLES only (Basilisk-style):
        u_h = (u_new + u_old)/2,  x_m = x_c - (dt/2) u_h(x_c),
        x_d = x_c - dt u_h(x_m),
    u_h(x_m) by interpolationCellPoint in the cell that contains x_m. Same
    formal order as taylor; reads no velocity gradient.
    Selected with footIntegrator rk2.

SourceFiles
    slRk2FootIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef slRk2FootIntegrator_H
#define slRk2FootIntegrator_H

#include "slFootIntegrator.H"

namespace Foam
{

class slRk2FootIntegrator
:
    public slFootIntegrator
{
public:

    TypeName("rk2");

    explicit slRk2FootIntegrator(const fvMesh& mesh);

    virtual ~slRk2FootIntegrator() = default;

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
