/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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
    Foam::redistancer

Description
    Abstract base class for implementing a redistancing/ reinitialization
    algorithmen for the Level-Set Method.

    Triggering (fvSolution levelSet.redistancer.trigger):
      - "interval" (default): redistance when
        timeIndex % redistanceInterval == 0 (original behaviour; every
        pre-existing case dictionary remains valid).
      - "gradPsiThreshold": the threshold criterion of Abu-Al-Saud, Popinet
        & Tchelepi (JCP 371:896, 2018, Eq. 12) -- redistance only when the
        volume-weighted L1 mean of ||grad psi| - 1| over the narrow band
        B = {|psi| < 6 h} exceeds a threshold. The threshold defaults to
        (h/L)^2 with h the minimum cell spacing and L the largest extent of
        the mesh bounding box (the paper states Delta^2 with Delta the
        NONDIMENSIONAL grid spacing; dividing by the domain extent
        reproduces that on a dimensional mesh). Both papers show that
        redistancing EVERY step destroys the discrete capillary equilibrium
        -- the criterion exists to keep events rare. The criterion value is
        reported every evaluation so that its floor (set by the anchor
        construction residual) can be calibrated against the threshold from
        a single run.

SourceFiles
    redistancer.C

\*---------------------------------------------------------------------------*/

#ifndef redistancer_H
#define redistancer_H

#include "dictionary.H"
#include "messageStream.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvSolution.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class redistancer Declaration
\*---------------------------------------------------------------------------*/

class redistancer
{
protected:

    const fvSolution& fvSolution_;
    const dictionary& levelSetDict_;
    const dictionary& redistDict_;

private:

    const label redistanceInterval_;

    //- Trigger mode: "interval" (original), "gradPsiThreshold" (documented
    //  negative result), or "signedDistanceBounds" (recommended)
    const word trigger_;

    //- Criterion threshold; < 0 selects the automatic (h/L)^2 default,
    //  resolved lazily on first use (needs the mesh).
    mutable scalar threshold_;

    //- signedDistanceBounds: redistance when the band |grad psi| leaves
    //  [1/g, g] -- the case- and resolution-independent usability bounds of
    //  the level set as a signed distance. |grad psi| -> 0: the zero
    //  crossing becomes unlocatable (interface loss); |grad psi| -> 1/h:
    //  the profile is a one-cell Heaviside jump (algebraic VOF). ln|grad
    //  psi| drifts at the interface-normal strain rate, so log-symmetric
    //  bounds fire equidistantly in strain-time for compression and
    //  expansion, and the event count is h-independent (the consistency
    //  requirement: O(h^2) displacement per step x O(1) steps).
    const scalar gradBound_;

    //- Optional robust variant: fire only when more than this fraction of
    //  band cells violates the bounds (0 = strict min/max).
    const scalar bandViolationFraction_;

    //- signedDistanceBounds criterion region:
    //    band (default) -- sign-change cells + one dilation layer (the
    //                      original behaviour; interface-quality criterion).
    //    bulk           -- the COMPLEMENT of the frozen anchor region
    //                      (sign-change band + anchorLayers dilations) plus
    //                      one junction layer. For the frozen-band bulk
    //                      re-extension: the anchors are untouched by
    //                      construction and the junction kink must not
    //                      retrigger, so neither is measured. NOTE: a signed
    //                      distance has |grad psi| = 0 on its skeleton
    //                      (medial axis), so the strict min/max firing rule
    //                      would fire forever there -- use the
    //                      bandViolationFraction robust variant with bulk.
    const word criterionRegion_;

protected:

    //- Anchor-region thickness in dilation layers beyond the sign-change
    //  band (shared key with the freezing fill models so the trigger and
    //  the fill agree on the region).
    const label anchorLayers_;

    //- mask = sign-change cells (registered NarrowBand == 1) dilated
    //  nLayers times, PARALLEL-CORRECT (halo exchange per sweep, so a
    //  guard layer never has a hole at a processor boundary -- the mask
    //  freezes VALUES, unlike the rank-local firing-criterion dilation).
    void topologicalMask
    (
        const volScalarField& psi,
        const label nLayers,
        boolList& mask
    ) const;

private:

    //- Volume-weighted L1 mean of ||grad psi| - 1| over {|psi| < 6h}
    scalar gradPsiBandL1(const volScalarField& psi) const;

    //- Extrema and bound-violating fraction of |grad psi| over the
    //  TOPOLOGICAL band (sign-change cells + one dilation layer -- NOT the
    //  value-based {|psi| < 6h}, which a degraded profile distorts).
    //  Gradient read through the dedicated "redistGradPsi" fvSchemes entry
    //  (decoupled from the advection limiter). Returns the global band
    //  cell count (0 = interface lost).
    label gradPsiBandBounds
    (
        const volScalarField& psi,
        scalar& gMin,
        scalar& gMax,
        scalar& violFrac
    ) const;

    //- Resolve the automatic threshold (h/L)^2 once
    scalar resolveThreshold(const fvMesh& mesh) const;

public:

    // Static Data Members

    TypeName("noRedistancing");


    declareRunTimeSelectionTable
    (
        autoPtr,
        redistancer,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );

    // Constructors

        //- Construct from components
        redistancer(const fvMesh& mesh);

    // Selectors

        //- Select default constructed
        static autoPtr<redistancer> New(const fvMesh& mesh);

    //- Destructor
    virtual ~redistancer() = default;

    // Member Functions

    //- Gate on the configured trigger and redistance if it fires.
    //  Returns true when a redistance event actually ran (callers use this
    //  to refresh sign-pattern-dependent state, e.g. the narrow band).
    bool redistance(volScalarField& psi);

    //- Perform the redistancing; returns whether psi was actually modified
    //  (the no-op base and models' early exits return false).
    virtual bool doRedistance(volScalarField& psi)
    {
        return false;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
