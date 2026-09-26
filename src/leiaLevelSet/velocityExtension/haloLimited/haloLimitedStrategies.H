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
    The strategy families of the halo-limited extension
    (Foam::velocityExtensions::haloLimited; halo-limited dossier, D. Bothe,
    2026-09-26, section "halo-limited directional velocity extension"):

    - extensionTravel    S_R(d): how far the sample point travels towards the
                         interface; odd, S_R(0) = 0, S_R'(0) = 1, |S_R| < R.
                         Its fraction reached c_R = S_R(d)/d (c_R(0) = 1).
    - extensionWeight    w(c_R): the weight of the sampled velocity.
    - extensionDirection (d, e) from psi and grad(psi): the nominal signed
                         distance and the direction away from the interface.
    - extensionSampler   the velocity at an arbitrary point Y near a cell, from
                         the cell's own local model of each velocity component.

    Each family is runtime-selectable from its sub-dictionary of
    levelSet.velocityExtension (\c travel, \c weight, \c direction,
    \c sampler); the entry \c type is required, and so is every coefficient.

SourceFiles
    haloLimitedStrategies.C

\*---------------------------------------------------------------------------*/

#ifndef haloLimitedStrategies_H
#define haloLimitedStrategies_H

#include "fvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class slReconstruction;

/*---------------------------------------------------------------------------*                        Class extensionTravel Declaration
\*---------------------------------------------------------------------------*/

class extensionTravel
{
public:

    TypeName("extensionTravel");

    declareRunTimeSelectionTable
    (
        autoPtr, extensionTravel, dictionary, (const dictionary& dict), (dict)
    );

    explicit extensionTravel(const dictionary&) {}

    static autoPtr<extensionTravel> New(const dictionary& dict);

    virtual ~extensionTravel() = default;

    //- The travelled distance S_R(d), |S_R(d)| < R
    virtual scalar S(const scalar d, const scalar R) const = 0;

    //- The fraction reached c_R(d) = S_R(d)/d, with c_R(0) = 1
    virtual scalar fraction(const scalar d, const scalar R) const = 0;
};


/*---------------------------------------------------------------------------*                        Class extensionWeight Declaration
\*---------------------------------------------------------------------------*/

class extensionWeight
{
public:

    TypeName("extensionWeight");

    declareRunTimeSelectionTable
    (
        autoPtr, extensionWeight, dictionary, (const dictionary& dict), (dict)
    );

    explicit extensionWeight(const dictionary&) {}

    static autoPtr<extensionWeight> New(const dictionary& dict);

    virtual ~extensionWeight() = default;

    //- The weight w in [0, 1] of the sampled velocity, from the fraction
    //  reached c_R in (0, 1]; w(1) = 1
    virtual scalar w(const scalar fractionReached) const = 0;
};


/*---------------------------------------------------------------------------*                      Class extensionDirection Declaration
\*---------------------------------------------------------------------------*/

class extensionDirection
{
public:

    TypeName("extensionDirection");

    declareRunTimeSelectionTable
    (
        autoPtr, extensionDirection, dictionary, (const dictionary& dict), (dict)
    );

    explicit extensionDirection(const dictionary&) {}

    static autoPtr<extensionDirection> New(const dictionary& dict);

    virtual ~extensionDirection() = default;

    //- The nominal signed distance d and the direction e (|e| <= 1, pointing
    //  from the interface to the point) from psi and g = grad(psi)
    virtual void evaluate
    (
        const scalar psi,
        const vector& g,
        scalar& d,
        vector& e
    ) const = 0;
};


/*---------------------------------------------------------------------------*                       Class extensionSampler Declaration
\*---------------------------------------------------------------------------*/

class extensionSampler
{
public:

    TypeName("extensionSampler");

    declareRunTimeSelectionTable
    (
        autoPtr,
        extensionSampler,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );

    extensionSampler(const fvMesh&, const dictionary&) {}

    static autoPtr<extensionSampler>
    New(const fvMesh& mesh, const dictionary& dict);

    virtual ~extensionSampler() = default;

    //- Build the local models of the scalar field f (one velocity component);
    //  f must stay alive until the last value() call
    virtual void update(const volScalarField& f) = 0;

    //- The local model of cell celli evaluated at the point x; equal to
    //  f[celli] at the cell centre
    virtual scalar value(const label celli, const point& x) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace extensionTravels
{

//- S_R(d) = d [1 + (d/R)^(2m)]^(-1/(2m)): S_R(d) = d - d^(2m+1)/(2m R^(2m)) + ...
//  near the interface, S_R -> R sign(d) far away.
class capped
:
    public extensionTravel
{
    //- The flatness exponent m >= 1
    const label m_;

public:

    TypeName("capped");

    explicit capped(const dictionary& dict);

    virtual scalar S(const scalar d, const scalar R) const;

    virtual scalar fraction(const scalar d, const scalar R) const;
};

} // End namespace extensionTravels


namespace extensionWeights
{

//- w = c_R^beta: the weight follows the fraction of the displacement reached,
//  so 1 - w = O((d/R)^(2m)) near the interface and w ~ (R/|d|)^beta far away.
class fractionReached
:
    public extensionWeight
{
    //- The exponent beta > 0
    const scalar beta_;

public:

    TypeName("fractionReached");

    explicit fractionReached(const dictionary& dict);

    virtual scalar w(const scalar c) const
    {
        return Foam::pow(c, beta_);
    }
};

} // End namespace extensionWeights


namespace extensionDirections
{

//- The pure level-set realization: d = psi/sqrt(Q2), e = g/sqrt(Q2), with
//  Q2 = q^2 + zeta(q^2) and the flat regularizer
//  zeta(z) = (1 - 4z)^3/4 for z < 1/4, zeta = 0 for z >= 1/4 (C2 at z = 1/4,
//  Q2 >= 0.154). zeta is zero for q >= 0.5, so d = psi/q and e = n exactly in
//  the admissible band and far beyond it (halo-limited dossier, flat
//  regularization of the algebraic point).
class levelSet
:
    public extensionDirection
{
public:

    TypeName("levelSet");

    explicit levelSet(const dictionary& dict);

    //- Q2 = q^2 + zeta(q^2)
    static scalar Q2(const scalar q2);

    virtual void evaluate
    (
        const scalar psi,
        const vector& g,
        scalar& d,
        vector& e
    ) const;
};

} // End namespace extensionDirections


namespace extensionSamplers
{

//- The cell's own weighted least-squares quadratic value fit over the
//  cell-point-cell stencil: the semi-Lagrangian reconstruction
//  uncachedQuadraticWeightedLeastSquares (same stencil, weights, linear
//  fallback and halo exchange), fitted to one velocity component at a time.
//  A quadratic model is required: a linear one turns the extension into a
//  weighted SDPLS source at first order (halo-limited dossier). The fit reads
//  its stencil options from levelSet.semiLagrangian.
class stencilFit
:
    public extensionSampler
{
    autoPtr<slReconstruction> fit_;

public:

    TypeName("stencilFit");

    stencilFit(const fvMesh& mesh, const dictionary& dict);

    virtual ~stencilFit();

    virtual void update(const volScalarField& f);

    virtual scalar value(const label celli, const point& x) const;
};

} // End namespace extensionSamplers

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
