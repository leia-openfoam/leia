/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slFootIntegrator

Description
    The departure point (foot) x_d of the characteristic that arrives at a cell
    centre x_c after one step dt. Run-time-selectable family; the base class is
    the default model `taylor`, the arrival-centred Taylor foot:

        x_d = x_c - u_new dt + (dt^2/2) [ (u_new - u_old)/dt + (u_new . grad) u_new ]

    with grad(u_new) the cell gradient (scheme name gradU). Second order per
    step in dt (the per-step error is O(dt^3)).

    The caller gives two velocity levels: u_new and u_old (the kinematic
    solvers give u(t^{n+1}) and u(t^n); the two-phase solver gives the
    Adams-Bashforth estimate and u^n, see concepts/departure-foot-ab2-centring).

    Dictionary: levelSet { semiLagrangian { footIntegrator <type>; } },
    default taylor. Models: taylor, rk2, taylor3.

SourceFiles
    slFootIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef slFootIntegrator_H
#define slFootIntegrator_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class slFootIntegrator
{
protected:

    const fvMesh& mesh_;

public:

    TypeName("taylor");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slFootIntegrator,
        Mesh,
        (const fvMesh& mesh),
        (mesh)
    );

    explicit slFootIntegrator(const fvMesh& mesh);

    //- Select from levelSet.semiLagrangian.footIntegrator (default taylor).
    static autoPtr<slFootIntegrator> New(const fvMesh& mesh);

    virtual ~slFootIntegrator() = default;

    //- The foot of every cell centre.
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
