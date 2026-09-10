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
    Foam::slCorrector

Description
    Runtime-selectable strategy for assembling psi^{n+1} from psi^n in the
    semi-Lagrangian step, GIVEN the (already computed, fixed) departure feet and
    a runtime-selected slReconstruction. This is the seam between the Taylor
    backtracking (which owns the feet, in slAdvection, and is never touched here)
    and the reconstruction of psi at those feet: it decides HOW the reconstruction
    is applied over the step.

    Selected with `levelSet { semiLagrangian { correction ...; } }` in fvSolution
    (default `direct`). Composes with any slReconstruction (`reconstruction ...`):

      * direct              -- one reconstruction pass (the plain scheme; what the
                               stable quadraticWeightedLeastSquares uses).
      * deferredCorrection  -- iterate the reconstruction, rebuilding it from the
                               current iterate each pass, so a lower-order gradient
                               used in the lift becomes self-consistent. (Research
                               path; see deferredCorrector.)

    Shared services for derived strategies: robustEvaluate() (reconstruction at a
    foot with the always-on value cap / optional monotone clip / non-finite reset)
    and footRadiusGuard() (the operational CFL <= 1 check). New strategies derive
    from slCorrector and register with addToRunTimeSelectionTable.

SourceFiles
    slCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef slCorrector_H
#define slCorrector_H

#include "fvMesh.H"
#include "volFields.H"
#include "pointField.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "slValueBound.H"

namespace Foam
{

class slReconstruction;

/*---------------------------------------------------------------------------*\
                        Struct slBoundTally Declaration
\*---------------------------------------------------------------------------*/

//- Per-step bound accounting. The derived corrector creates one, robustEvaluate
//  fills it cell by cell, and reportClipActivity reduces and writes it. Plain
//  per-cell data: nothing here is reduced until the whole cell loop is done, so
//  no collective is ever reached from inside a loop over cells.
struct slBoundTally
{
    //- Cells whose value the bound actually changed this step.
    boolList fired;

    //- Cells where the bound's own construction gave an EMPTY interval, so the
    //  recovery path supplied one. For the distance-cone bound this is a direct
    //  measurement of how far psi has drifted from a distance function.
    boolList inadmissible;

    //- The change the bound made, psi_bounded - psi_fit, in the units of psi
    //  (metres for a level set). Zero where the bound did not act. This is the
    //  interface damage, in the same quantity the falsified monotone clip was
    //  measured with.
    scalarField delta;

    //- Signed headroom, min(hi - v, v - lo)/max(|d|, SMALL), where the bound is
    //  enforced. NEGATIVE means the bound fired, POSITIVE is spare room, and
    //  zero means the bound was not enforced in that cell.
    scalarField slack;

    //- Non-finite reconstructions reset to the stencil mid this step.
    label nNonFinite = 0;

    explicit slBoundTally(const label nCells)
    :
        fired(nCells, false),
        inadmissible(nCells, false),
        delta(nCells, Zero),
        slack(nCells, Zero)
    {}

    //- Reuse across deferred-correction passes without reallocating.
    void reset()
    {
        fired = false;
        inadmissible = false;
        delta = Zero;
        slack = Zero;
        nNonFinite = 0;
    }
};


/*---------------------------------------------------------------------------*\
                         Class slCorrector Declaration
\*---------------------------------------------------------------------------*/

class slCorrector
{
protected:

    const fvMesh& mesh_;

    //- Copy of the semiLagrangian sub-dict (models read their own controls).
    dictionary dict_;

    //- The runtime-selected bound on the reconstructed value. Never null: the
    //  sentinel `valueBound fromClipSwitch` resolves to `none` when the legacy
    //  clipToStencilBounds Switch is false, so a case that predates the family
    //  gets the unbounded path bit for bit.
    autoPtr<slValueBound> bound_;

    //- Cumulative clip activity. A write-time SAMPLE of the per-step count is not
    //  enough and has already misled once: the inflowOnly gate reported 0 bounded
    //  cells at every write while its metric CSV left the baseline at step 20. The
    //  clip fires in bursts, so the total and the first firing are what matter.
    mutable label nClipFiredTotal_ = 0;
    mutable label nClipStepsFired_ = 0;
    mutable label firstFireIndex_ = -1;

    //- Cumulative cells whose bound interval was EMPTY. For the distance-cone
    //  bound a large value FALSIFIES the guard: the old field is not Lipschitz
    //  over those stencils, so the bound has no meaning there.
    mutable label nInadmissibleTotal_ = 0;

    //- Cells that have EVER been bounded, written as slClipFiredEver. The per-step
    //  field alone cannot find a firing that happened between two write times, and
    //  that is exactly the case that matters: with the extremum exemption the clip
    //  fires 3 cell-steps in the first few steps of a 1563-step run and never
    //  again, yet that one-off perturbation grows into a 1 % volume error.
    mutable boolList firedEver_;

    //- Reconstruct at a foot with the advection-level robustness safeguards:
    //  non-finite -> stencil mid, counted; then clip to the interval the
    //  runtime-selected slValueBound returns for this cell. Every bound returns
    //  the always-on runaway cap (stencil mid +- 10*range) for a cell it does not
    //  act on, so the value can never grow without limit whichever bound is
    //  selected. This is the ONE decision point for the bound, which is what
    //  makes the strategies comparable as a study axis.
    scalar robustEvaluate
    (
        const slReconstruction& recon,
        const label c,
        const point& foot,
        const bool eligible,
        slBoundTally& tally
    ) const;

    //- Per-cell mask of the cells where the quasi-monotone clip may act, built
    //  from psi^n once per step. Empty when the clip is off. The criterion is
    //  stated per cell over face neighbours, so it is compact-stencil,
    //  mesh-agnostic and decomposition-independent -- see slCorrector.C.
    void buildClipMask
    (
        const volScalarField& psi,
        const slReconstruction& recon,
        boolList& clipCell
    ) const;

    //- Accumulate the clip's activity EVERY step, and at write time report the
    //  per-step count, the cumulative total, the number of steps that fired and
    //  the first step that fired; then write slClipEligible / slClipFired so the
    //  firing cells can be LOCATED. A count says the clip acted; only the field
    //  says where, and where is what decides whether the region rule is right.
    void reportClipActivity
    (
        const slReconstruction& recon,
        const boolList& clipCell,
        const slBoundTally& tally
    ) const;

    //- Operational CFL <= 1 check: warn if any foot leaves its point-neighbour
    //  stencil hull (would make evaluate() extrapolate). Geometry only.
    void footRadiusGuard
    (
        const slReconstruction& recon,
        const pointField& feet
    ) const;

    //- Parallel-reduce and warn on non-finite reconstructions.
    void warnNonFinite(label nNonFinite) const;


public:

    //- Runtime type information
    TypeName("slCorrector");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slCorrector,
        Mesh,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );

    //- Construct from mesh and the semiLagrangian sub-dict.
    slCorrector(const fvMesh& mesh, const dictionary& dict);

    //- Select from semiLagrangian.correction (default "direct").
    static autoPtr<slCorrector> New(const fvMesh& mesh, const dictionary& dict);

    //- The selected bound on the reconstructed value.
    const slValueBound& bound() const
    {
        return bound_();
    }

    //- Destructor
    virtual ~slCorrector() = default;

    //- Assemble psi^{n+1} IN PLACE from psi^n (on entry) using the reconstruction
    //  evaluated at the fixed departure feet. May call recon.update() one or more
    //  times (rebuilding it from the current iterate).
    virtual void correct
    (
        volScalarField& psi,
        const pointField& feet,
        slReconstruction& recon
    ) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
