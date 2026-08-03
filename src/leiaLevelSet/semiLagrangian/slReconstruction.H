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

    //- Legacy switch: limitSlope true == slopeLimiter barthJespersen (kept for
    //  backward compatibility with existing dictionaries).
    Switch limitSlope_;

    //- Slope-limiter family (default derived from limitSlope_ for compat):
    //    none            -- phi_ = 1, unlimited (default; bitwise unchanged).
    //    barthJespersen  -- hard min-ratio cap; MEASURED to collapse the
    //                       quadratic convergence order (3.0 -> 0.1) over the
    //                       wide point stencil.
    //    venkatakrishnan -- SMOOTH cap that deactivates in smooth regions
    //                       (epsilon^2 = (venkK*h)^3), so it bounds overshoot
    //                       near extrema WITHOUT collapsing the order.
    const word limiterType_;

    //- Venkatakrishnan threshold constant K (epsilon^2 = (venkK_*radius_[c])^3).
    //  Larger K => limiter deactivates over a wider smooth band (more accuracy,
    //  less bounding). Default 5.
    const scalar venkK_;

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

        //- Fill kappa (1/length) with the interface mean curvature
        //  kappa = div(grad psi / |grad psi|) evaluated SYMBOLICALLY from the
        //  reconstruction's fit gradient and Hessian (no FVM operators), constant-
        //  extended along the normal via the foot-point. Requires a prior update()
        //  on the current psi. Only the quadratic models provide this.
        virtual void meanCurvature(volScalarField& kappa) const
        {
            NotImplemented;
        }

        //- The SDF-assuming curvature kappa = Laplacian(psi) = tr(H), read from the
        //  fit Hessian. Equal to div(grad psi/|grad psi|) ONLY when psi is a signed
        //  distance (|grad psi| = 1); provided for comparison with meanCurvature.
        virtual void meanCurvatureLaplacian(volScalarField& kappa) const
        {
            NotImplemented;
        }

        //- As meanCurvature (same div(grad psi/|grad psi|) from the fit), but
        //  evaluated AT THE CELL CENTRE -- WITHOUT the foot-point projection, so it
        //  is NOT extended constant along the interface normal. Provided to compare
        //  curvature accuracy with and without the normal extension.
        virtual void meanCurvatureNoExtension(volScalarField& kappa) const
        {
            NotImplemented;
        }

        //- As meanCurvature -- symbolic kappa at the interface foot of the
        //  cell's own quadratic -- but the foot is the TRUE closest point,
        //  found by Newton on the KKT system of min |delta|^2 s.t.
        //  P(delta) = 0 (Jacobian uses the constant fit Hessian), with a
        //  trust region |delta| <= stencilRadius. Falls back per cell to the
        //  cell-centre (NoExtension) value on any guard failure, including
        //  feet beyond the trust region (outer band cells): only the
        //  force-support layer is foot-evaluated.
        virtual void meanCurvatureClosestPoint(volScalarField& kappa) const
        {
            NotImplemented;
        }

        //- Fitted gradient and Hessian of the local polynomial AT THE CELL
        //  CENTRE, read from the fit coefficients (no FVM re-differentiation).
        //  Returns the fit order actually available in this cell:
        //    2  full quadratic (g and H valid),
        //    1  linear fallback (g valid, H = 0),
        //    0  constant fallback (g = H = 0).
        //  Requires a prior update() on the current psi. Only the quadratic
        //  value-fit models provide this; used by the normal-projected SL
        //  scheme, which consults the fit for directions and geometry only.
        virtual label fitDerivatives
        (
            const label c,
            vector& g,
            symmTensor& H
        ) const
        {
            NotImplemented;
            return 0;
        }

        //- Signed geometric offset from the cell centre to the zero set of the
        //  cell's own fit along the normal ray: the stable quadratic root
        //      d_c = 2 psi_c / (|g| + sqrt(D)),  D = |g|^2 - 2 psi_c h_nn,
        //  with h_nn = n.H.n, accepted while D >= offsetBeta*|g|^2 and replaced
        //  by the first-order d_c = psi_c/|g| otherwise (normal-projected SL
        //  design note, Sec. 5). Exactly zero-preserving: psi_c = 0 => d_c = 0.
        //  fallback reports whether the first-order branch was taken. This is
        //  the SINGLE offset-distance definition shared by the nSL transport
        //  and the curvature offsetCorrection quadraticRoot mode, so transport
        //  and curvature never disagree about a cell's offset.
        virtual scalar signedOffset(const label c, bool& fallback) const
        {
            NotImplemented;
            return 0;
        }

        //- Signed geometric distance from the query point p to the iso-surface
        //  {R_c = level} of cell c's own local model, by the STABILIZED
        //  FOOT-POINT algorithm (docs/normal-projected-semi-lagrangian/
        //  stable-foot-point-3d.md): alternating tangent-plane/surfacepoint
        //  projections accelerated by the parabola correction built from the
        //  two step vectors of each cycle -- curvature-aware at the price of
        //  first derivatives only. All evaluations are polynomial arithmetic
        //  in the fit coefficients; no field access. The sign follows the
        //  model gradient at the foot (positive on the +psi side). ok reports
        //  convergence AND the foot staying inside the stencil trust region;
        //  on ok = false the returned value must not be used. This is the
        //  offset engine of nSL Variant (a).
        virtual scalar footPointDistance
        (
            const label c,
            const point& p,
            const scalar level,
            bool& ok
        ) const
        {
            NotImplemented;
            return 0;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
