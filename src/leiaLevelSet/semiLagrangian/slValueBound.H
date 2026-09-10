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
    Foam::slValueBound

Description
    Runtime-selectable BOUND on the value the semi-Lagrangian reconstruction
    returns at a departure foot point. The bound decides one thing only: the
    interval [lo, hi] that psi^{n+1} in the arrival cell must lie in. It never
    touches the reconstruction, the feet or the field.

    Selected with `levelSet { semiLagrangian { valueBound ...; } }` in
    fvSolution:

      * none            -- no bound. Only the always-on runaway cap
                           (stencil mid +- 10*stencil range) that predates this
                           class. Reproduces `clipToStencilBounds false`.
      * stencilBounds   -- the quasi-monotone clip to the stencil [min, max],
                           with clipRegion and clipKeepExtrema. Reproduces
                           `clipToStencilBounds true`. FALSIFIED as a fix, kept
                           so its measurements stay reproducible.
      * lipschitzCone   -- the distance-cone bound from the Lipschitz property
                           of a signed distance function.

    The default is the sentinel word `fromClipSwitch`, which resolves to
    `stencilBounds` when the legacy Switch `clipToStencilBounds` is true and to
    `none` otherwise. Every case that predates this class therefore keeps its
    behaviour bit for bit, whether or not its template carries the token.

    WHY A FAMILY. The quadratic fit is not a convex combination of the stencil
    values, so it can create a new extremum, and the frozen-velocity power
    iteration measures a spectral radius above one on EVERY mesh tested
    (hexahedral production meshes included). The first bound tried, the
    quasi-monotone clip, was FALSIFIED: 59.2 % of the cells it must act on are
    themselves stencil extrema, because a signed distance field has genuine
    interior extrema on the medial axis. Several further candidates exist and
    they must be comparable as ONE study axis rather than as a sequence of
    branches inside robustEvaluate.

    SCOPE. The bound acts in slCorrector::robustEvaluate, so it reaches the
    `pointValue` scheme only -- which is the production default. fluxFormScheme
    calls slReconstruction::evaluate() directly and normalProjectedScheme never
    uses the corrector, so both stay unbounded.

SourceFiles
    slValueBound.C

\*---------------------------------------------------------------------------*/

#ifndef slValueBound_H
#define slValueBound_H

#include "fvMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class slReconstruction;

/*---------------------------------------------------------------------------*\
                        Struct slBoundResult Declaration
\*---------------------------------------------------------------------------*/

//- What a bound reports for one cell. Plain data, filled per cell, no
//  collectives: a bound is called once per cell per step and a reduction there
//  deadlocks as soon as the ranks hold unequal cell counts.
struct slBoundResult
{
    //- The interval the value must lie in.
    scalar lo = -GREAT;
    scalar hi = GREAT;

    //- True when [lo, hi] is a REAL bound. False when it is the loose runaway
    //  cap, so that a cap hit is never counted as a firing.
    bool enforced = false;

    //- True when the bound's own construction gave an EMPTY interval and the
    //  recovery path had to supply one.
    bool inadmissible = false;
};


/*---------------------------------------------------------------------------*\
                         Class slValueBound Declaration
\*---------------------------------------------------------------------------*/

class slValueBound
{
protected:

    const fvMesh& mesh_;

    //- Copy of the semiLagrangian sub-dict (models read their own controls).
    dictionary dict_;

    //- The runaway cap that predates this class: stencil mid +- 10*range. It is
    //  what every bound returns for a cell it does not act on, so the value can
    //  still not grow without limit. `enforced` stays false.
    inline void runawayCap
    (
        const scalar stencilLo,
        const scalar stencilHi,
        slBoundResult& r
    ) const
    {
        const scalar mid = 0.5*(stencilLo + stencilHi);
        const scalar cap = 10.0*Foam::max(stencilHi - stencilLo, SMALL);
        r.lo = mid - cap;
        r.hi = mid + cap;
        r.enforced = false;
    }


public:

    //- Runtime type information
    TypeName("slValueBound");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slValueBound,
        Mesh,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );

    //- Construct from mesh and the semiLagrangian sub-dict.
    slValueBound(const fvMesh& mesh, const dictionary& dict);

    //- Select from semiLagrangian.valueBound. The default sentinel
    //  `fromClipSwitch` resolves from the legacy clipToStencilBounds Switch,
    //  which the caller reads and passes in.
    static autoPtr<slValueBound> New
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const bool legacyClipSwitch
    );

    //- Destructor
    virtual ~slValueBound() = default;


    // Member Functions

        //- The interval the reconstructed value at `foot` must lie in.
        //  stencilLo/stencilHi are the arrival cell's stencil range, passed in
        //  so no bound recomputes it. `eligible` is the per-cell mask.
        virtual void interval
        (
            const slReconstruction& recon,
            const label c,
            const point& foot,
            const bool eligible,
            const scalar stencilLo,
            const scalar stencilHi,
            slBoundResult& r
        ) const = 0;

        //- True when this bound never enforces anything, so the corrector can
        //  skip its diagnostics and the spectrum probe can treat the transport
        //  operator as LINEAR.
        virtual bool inert() const
        {
            return false;
        }

        //- One line for the solver banner, naming the settings that matter.
        virtual void printBanner() const
        {}
};


} // End namespace Foam

#endif

// ************************************************************************* //
