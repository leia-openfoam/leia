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

namespace Foam
{

class slReconstruction;

/*---------------------------------------------------------------------------*\
                         Class slCorrector Declaration
\*---------------------------------------------------------------------------*/

class slCorrector
{
protected:

    const fvMesh& mesh_;

    //- Copy of the semiLagrangian sub-dict (models read their own controls).
    dictionary dict_;

    //- Reconstruct at a foot with the advection-level robustness safeguards:
    //  non-finite -> stencil mid (counts nNonFinite); optional tight monotone
    //  clip to the stencil bounds; else an always-on generous value cap
    //  (mid +- 10*range) that only fires on runaway growth.
    scalar robustEvaluate
    (
        const slReconstruction& recon,
        const label c,
        const point& foot,
        const bool clip,
        label& nNonFinite
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
