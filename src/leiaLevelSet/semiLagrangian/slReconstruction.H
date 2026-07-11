/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::slReconstruction

Description
    Abstract base for runtime-selectable spatial reconstructions of the OLD
    level-set field at the semi-Lagrangian departure (foot) point. Selected by
    `levelSet { semiLagrangian { reconstruction ...; } }` in fvSolution.

    Contract: every reconstruction is centred on the ARRIVAL cell c and
    REPRODUCES the cell value there,

        evaluate(c, x_c) == psiOld[c]   (to machine precision),

    so that under the CFL <= 1 foot bound the evaluation at x_d is always
    interpolation inside c's point-neighbour (cell-point-cell) stencil, never
    extrapolation. The arrival-cell polynomial is built from a one-layer
    cell-point-cell (CPC) stencil whose remote (processor/cyclic) neighbours are
    gathered by OpenFOAM's extendedCentredCellToCellStencil map -- so evaluate()
    is only ever called for local arrival cells and needs no cell search.

    Base services: the cached CPC stencil, a per-step collect of psiOld (and a
    one-off collect of the cell centres) over that stencil, and the stencil
    min/max used by the optional monotonicity clip.

SourceFiles
    slReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef slReconstruction_H
#define slReconstruction_H

#include "fvMesh.H"
#include "fvSolution.H"
#include "volFields.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "extendedCentredCellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class slReconstruction Declaration
\*---------------------------------------------------------------------------*/

class slReconstruction
{
protected:

    const fvMesh& mesh_;
    const fvSolution& fvSolution_;
    const dictionary& levelSetDict_;
    //- Copy: subOrEmptyDict returns a temporary when the sub-dict is absent.
    dictionary slDict_;

    //- Optional monotonicity clip of the reconstructed value to the stencil
    //  bounds of psiOld (default false; quadratic fits can overshoot).
    Switch clipToStencilBounds_;

    //- Barth-Jespersen slope limiting of the reconstruction (default on): scale
    //  each cell's reconstruction increment by phi_[c] in [0,1] so it does not
    //  overshoot the stencil min/max at the neighbour centres -- prevents the
    //  runaway that a raw high-slope extrapolation causes, without a post-hoc
    //  value clamp. phi_ = 1 => unlimited.
    Switch limitSlope_;

    //- Per-cell Barth-Jespersen limiter (filled by computeLimiters()).
    scalarField phi_;

    //- Fill phi_ from the stencil bounds and the model's evaluateRaw() at the
    //  stencil neighbour centres. Call at the end of each model's update().
    void computeLimiters();

    //- Cached cell-point-cell stencil (arrival cell = stencil[c][0]); remote
    //  neighbours gathered via its mapDistribute in collectStencil().
    const extendedCentredCellToCellStencil& stencil_;

    //- Number of local cells. A compact stencil index g < nLocal_ is a LOCAL
    //  cell whose centre is read directly from mesh_.C() -- so the local stencil
    //  geometry is NOT duplicated per stencil entry (the dominant per-cell saving
    //  on the retained-centre models: >90% of interior cells store nothing here).
    label nLocal_;

    //- Compact-array TAIL: the [nLocal_ .. map.constructSize()) slice of the
    //  distributed cell centres -- boundary-face + processor-halo entries only,
    //  built once (geometry static). Indexed by (g - nLocal_) in stencilC().
    List<vector> centreTail_;
    bool haveCentres_;

    //- Per-cell interpolation radius (max stencil-centre distance), cached on the
    //  first collectStencil() via the stencilC() accessor.
    scalarField radius_;

    //- Per-cell stencil psiOld (compact, incl. remote) -- collected per step
    //  when a model or the clip needs it.
    List<List<scalar>> stencilPsi_;

    //- The field passed to the most recent update().
    const volScalarField* psiOldPtr_;

    //- Fill stencilPsi_ from psiOld (one paired halo exchange); build the centre
    //  tail + radius on the first call. Call once per step, uniformly on all ranks.
    void collectStencil(const volScalarField& psiOld);

    //- Build centreTail_ (+ set nLocal_): the boundary-face + processor-halo slice
    //  of the compact distributed cell centres. Called once (mesh static).
    void buildCentreTail();

    //- Number of entries in cell c's compact stencil (self + neighbours).
    inline label stencilSize(const label c) const
    {
        return stencil_.stencil()[c].size();
    }

    //- Centre of the i-th stencil entry of cell c. Local entries (compact index
    //  < nLocal_) are read straight from mesh_.C(); boundary/halo entries from
    //  centreTail_. Reproduces the former stencilC_[c][i] with no per-cell store.
    inline point stencilC(const label c, const label i) const
    {
        const label g = stencil_.stencil()[c][i];
        return (g < nLocal_) ? mesh_.C()[g] : centreTail_[g - nLocal_];
    }

    //- No-op retained for API compatibility: the per-cell stencil centres are no
    //  longer stored (stencilC() reads mesh_.C()/centreTail_ on demand), so there
    //  is nothing to free. Kept so existing callers compile unchanged.
    void releaseStencilCentres() {}


public:

    //- Runtime type information
    TypeName("slReconstruction");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slReconstruction,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );

    // Constructors

        //- Construct from mesh (reads the semiLagrangian sub-dict).
        explicit slReconstruction(const fvMesh& mesh);

    // Selectors

        //- Select from levelSet.semiLagrangian.reconstruction (default
        //  quadraticWeightedLeastSquares).
        static autoPtr<slReconstruction> New(const fvMesh& mesh);

    //- Destructor
    virtual ~slReconstruction() = default;

    // Member Functions

        Switch clipToStencilBounds() const { return clipToStencilBounds_; }

        //- Min/max of psiOld over the arrival cell's stencil (requires a prior
        //  collectStencil in update()).
        void stencilRange(const label c, scalar& lo, scalar& hi) const;

        //- Max distance from the arrival cell centre to any stencil-cell centre
        //  (the interpolation radius; the SL foot must stay within it). Requires
        //  a prior collectStencil in update().
        scalar stencilRadius(const label c) const;

        //- Prepare the reconstruction for the current step from psiOld
        //  (gradient/Hessian passes or the cached-pseudoinverse coefficients).
        //  Implementations MUST call computeLimiters() at the end.
        virtual void update(const volScalarField& psiOld) = 0;

        //- The UNLIMITED reconstruction at point x, using the arrival cell c's
        //  own polynomial. evaluateRaw(c, C[c]) == psiOld[c].
        virtual scalar evaluateRaw(const label c, const point& x) const = 0;

        //- The slope-limited reconstruction used by the advection:
        //  psi_c + phi_[c] * (evaluateRaw(c,x) - psi_c). evaluate(c,C[c]) == psi_c.
        scalar evaluate(const label c, const point& x) const;

        //- Optional post-advection fix-up of the new psi field, called once per
        //  step after the per-cell assignment (e.g. the band model's far-field
        //  signed-distance re-extension). Default: no-op.
        virtual void postAdvect(volScalarField& psi) {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
