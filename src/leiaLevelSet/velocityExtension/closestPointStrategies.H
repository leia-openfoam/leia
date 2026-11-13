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

Description
    The two strategy families of the closestPoint velocity extension
    (plan docs/plan-velocity-extension-verification.md, section 4a, repair 3):

    closestPointStep (dictionary word cpStep) -- the descent step from the
    point x with the level-set value psi, the unit normal n and the gradient
    magnitude |grad psi| = q:
        distance (DEFAULT, the shipped step): dx = -psi n, which reads psi as
                 a signed distance;
        newton:  dx = -psi n/q = -psi grad(psi)/|grad psi|^2, the Newton step
                 of D. Bothe's notes (2026-10-29); exact for a linear psi of
                 any slope q.

    closestPointFill (dictionary word cpFill) -- the fill of the band cells
    whose descent failed:
        steadyUpwind (DEFAULT, the shipped fill): fill() returns false and the
                 model runs its steady upwind solve pinned at every successful
                 cell;
        footPointTaylor: a compact propagation of the foot-point data. Each
                 failed cell takes, from its face neighbours (processor
                 neighbours included), the foot point p nearest to its centre
                 C, with the normal n, the velocity u and the velocity gradient
                 G at p; Jacobi sweeps, so the result does not depend on the
                 cell order or on the decomposition. Then
                     Pi = p + (I - n n).(C - p),  Uext = u + (Pi - p).G,
                 second order in the tangential offset |Pi - p|.
    MEASURED 2026-11-13 (STATUS 11.31 item 8): with fadePlacement
    extendThroughFade, np 4 sends 7 % of the band to the fallback fill (serial:
    0 to 1 cell), because a descent that reaches 8 h leaves the one-layer halo;
    the steadyUpwind fill does not converge (STATUS 11.30 item 6).

SourceFiles
    closestPointStrategies.C

\*---------------------------------------------------------------------------*/

#ifndef closestPointStrategies_H
#define closestPointStrategies_H

#include "fvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class closestPointStep Declaration
\*---------------------------------------------------------------------------*/

class closestPointStep
{
public:

    TypeName("closestPointStep");

    declareRunTimeSelectionTable
    (
        autoPtr, closestPointStep, dictionary, (const dictionary& dict), (dict)
    );

    explicit closestPointStep(const dictionary&) {}

    //- Select by the word cpStep (default distance)
    static autoPtr<closestPointStep> New(const dictionary& dict);

    virtual ~closestPointStep() = default;

    //- True if step() reads the gradient magnitude (the caller computes it)
    virtual bool needsGradMag() const = 0;

    //- The descent step from psi, the unit normal n and q = |grad psi|
    virtual vector step
    (
        const scalar psi, const vector& n, const scalar q
    ) const = 0;
};


/*---------------------------------------------------------------------------*\
                      Class closestPointFill Declaration
\*---------------------------------------------------------------------------*/

class closestPointFill
{
public:

    TypeName("closestPointFill");

    declareRunTimeSelectionTable
    (
        autoPtr, closestPointFill, dictionary, (const dictionary& dict), (dict)
    );

    explicit closestPointFill(const dictionary&) {}

    //- Select by the word cpFill (default steadyUpwind)
    static autoPtr<closestPointFill> New(const dictionary& dict);

    virtual ~closestPointFill() = default;

    //- True if fill() reads the foot-point data (the caller records it)
    virtual bool needsFootData() const = 0;

    //- Fill Uext in the failed cells. COLLECTIVE: called on every rank when
    //  any rank has a failed cell. Returns false if the caller must fill.
    virtual bool fill
    (
        const fvMesh& mesh,
        volVectorField& Uext,
        const boolList& fixed,
        const boolList& failed,
        const pointField& footP,
        const vectorField& footN,
        const tensorField& footG
    ) const = 0;
};


// * * * * * * * * * * * * * * * Concrete strategies  * * * * * * * * * * * //

namespace closestPointSteps
{

class distance : public closestPointStep
{
public:
    TypeName("distance");
    explicit distance(const dictionary& dict) : closestPointStep(dict) {}
    bool needsGradMag() const override { return false; }
    vector step(const scalar psi, const vector& n, const scalar) const override
    {
        return -psi*n;
    }
};

class newton : public closestPointStep
{
    //- Floor on q = |grad psi| (a flat psi gives no direction anyway)
    scalar qFloor_;
public:
    TypeName("newton");
    explicit newton(const dictionary& dict);
    bool needsGradMag() const override { return true; }
    vector step
    (
        const scalar psi, const vector& n, const scalar q
    ) const override
    {
        return -psi*n/Foam::max(q, qFloor_);
    }
};

} // End namespace closestPointSteps


namespace closestPointFills
{

class steadyUpwind : public closestPointFill
{
public:
    TypeName("steadyUpwind");
    explicit steadyUpwind(const dictionary& dict) : closestPointFill(dict) {}
    bool needsFootData() const override { return false; }
    bool fill
    (
        const fvMesh&, volVectorField&, const boolList&, const boolList&,
        const pointField&, const vectorField&, const tensorField&
    ) const override
    {
        return false;
    }
};

class footPointTaylor : public closestPointFill
{
    //- Upper bound on the Jacobi sweeps
    label maxSweeps_;
public:
    TypeName("footPointTaylor");
    explicit footPointTaylor(const dictionary& dict);
    bool needsFootData() const override { return true; }
    bool fill
    (
        const fvMesh& mesh,
        volVectorField& Uext,
        const boolList& fixed,
        const boolList& failed,
        const pointField& footP,
        const vectorField& footN,
        const tensorField& footG
    ) const override;
};

} // End namespace closestPointFills

} // End namespace Foam

#endif

// ************************************************************************* //
