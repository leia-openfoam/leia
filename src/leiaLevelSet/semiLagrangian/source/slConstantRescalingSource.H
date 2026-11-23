/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slConstantRescalingSource

Description
    The rescaling of the level set by ONE CONSTANT (supervisor, 2026-11-23,
    STATUS 11.44). After the characteristic step the model computes

        G = mean |g| over the interface cells,

    with g the gradient of the semi-Lagrangian fit at the cell centre and the
    interface cells the cells of the registered field NarrowBand (narrowBand
    type signChange). The first call stores G as the reference G0. A call acts
    if |G/G0 - 1| > theta, and then multiplies psi in EVERY cell by G0/G.

    Properties.
    - The zero set and the sign of every cell value are unchanged, and so are
      the normal and the curvature of every level set: the factor is a constant.
    - The map adds no jump, because it has no region edge.
    - After an event G = G0, so the event clears its own trigger. The trigger of
      the model `rescaling` (maximum deviation over a region) is not cleared by
      its map; that map then ran in every step and destroyed the field (negative
      log 2026-11-23).
    - The reference is the initial mean gradient, not 1: the model does not ask
      for a signed distance.
    - It removes a uniform drift of |grad psi| only. It does not repair a
      gradient that varies along the interface.

    G0 is kept in the registry of the mesh (slConstantRescalingReference), so it
    survives the reconstruction of the semi-Lagrangian objects after a mesh
    change. It is not written: a restart takes a new reference.

    Selected with `semiLagrangian { source { type constantRescaling; } }`
    (token SL_SOURCE). Key: theta (token SL_RESCALE_THETA, default 0.1, one
    dimensionless number for every case).

SourceFiles
    slConstantRescalingSource.C

\*---------------------------------------------------------------------------*/

#ifndef slConstantRescalingSource_H
#define slConstantRescalingSource_H

#include "slSource.H"

namespace Foam
{

class slConstantRescalingSource
:
    public slSource
{
    const scalar theta_;

public:

    TypeName("constantRescaling");

    explicit slConstantRescalingSource(const fvMesh& mesh);

    virtual ~slConstantRescalingSource() = default;

    virtual bool active() const
    {
        return true;
    }

    virtual void apply
    (
        volScalarField& psi,
        const volVectorField& Utrajectory,
        slReconstruction& geometryFit,
        const scalar dt
    );
};

} // End namespace Foam

#endif
