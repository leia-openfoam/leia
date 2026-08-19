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

\*---------------------------------------------------------------------------*/

#include "exponentialDiscretization.H"
#include "addToRunTimeSelectionTable.H"
#include "sdplsSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialDiscretization, false);
    addToRunTimeSelectionTable
        (
            discretization,
            exponentialDiscretization,
            Dictionary
        );
}

// * * * * * * * * * * * *  File-local implementation  * * * * * * * * * * * //

namespace
{

// Overflow clamp on the ODE exponent z = f_nl*dt. See the class docs for the
// full argument; the two numbers that matter here are exp(30) = 1.07e13 (so
// the whole product chain V*psi*exp(z)/dt peaks at ~1e3, over 300 decades
// below the double-precision overflow at 1.8e308) and the MEASURED exponent
// of the coupled N=64 stationary droplet, z ~ 0.023, which is 1300x smaller.
// A run that reaches this clamp has an already-diverged velocity field; the
// clamp converts what would be inf/NaN -- silently poisoning alpha, U and p
// through the phase indicator -- into a large but FINITE source that fails
// visibly in the psi linear solver.
static const Foam::scalar zMax = 30.0;

// Below this |z| the first exponential-integrator weight is evaluated by its
// Maclaurin series instead of the quotient. The two branches are matched: at
// |z| = 1e-3 the series truncation after z^3/24 is z^4/120 = 8.3e-15 relative,
// while the quotient's cancellation error is eps/|z| = 2.2e-13 relative and
// only falls with growing |z|, so the worst relative error anywhere is
// 2.2e-13 -- better than two orders inside the 1e-10 tolerance of
// leiaTestSdplsSource.
static const Foam::scalar zSeries = 1.0e-3;

//- phi1(z) = (exp(z) - 1)/z, the first weight of the exponential integrator
//  family, with phi1(0) = 1.
//
//  Sc is built as psi^n * a * phi1(z) rather than as psi^n*(exp(z) - 1)/dt so
//  that the z -> 0 limit is EXACTLY a*psi^n in floating point: this member
//  then degenerates continuously to the tangent linearization that `explicit`
//  and `simpleLinearImplicit` apply, instead of losing digits to the
//  cancellation exp(z) - 1 (relative error ~ eps/|z|, i.e. 2.2e-12 already at
//  z = 1e-4). The affine unit test asserts the z = 0 row (alpha = 0) at 1e-10,
//  which the literal quotient would still pass but only by luck of exp(0)
//  being exactly 1; the series branch makes it structural.
inline Foam::scalar phi1(const Foam::scalar z)
{
    if (Foam::mag(z) < zSeries)
    {
        // Horner form of 1 + z/2 + z^2/6 + z^3/24.
        return 1.0 + z*(0.5 + z*((1.0/6.0) + z*(1.0/24.0)));
    }

    return (Foam::exp(z) - 1.0)/z;
}

} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponentialDiscretization::exponentialDiscretization()
    :
        discretization()
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

Foam::tmp<scalarField>
Foam::exponentialDiscretization::
Sc(const volScalarField& nonLinearPart, const volScalarField& psi) const
{
    // dt comes from the mesh's time registry: this linearization is the only
    // one in the hierarchy whose coefficient is not scale-free in time. The
    // exact factor exp(f_nl*dt) is a function of the STEP, not of f_nl alone,
    // which is precisely why the other three members can be written without dt
    // and this one cannot.
    const scalar dt = psi.mesh().time().deltaTValue();

    // psi^n is psi.oldTime(), NOT the psi argument.
    //
    // All three call sites (eulerianAdvection, leiaRedistancedLevelSetFoam,
    // alphaEqn.H) re-assemble this source inside an outer loop that has
    // already overwritten psi -- nDefCorr (default 3) or the PIMPLE outer
    // correctors. Feeding the latest iterate psi^k into a MULTIPLICATIVE
    // factor turns that loop into the Picard iteration
    //
    //     psi^{k+1} = psi^n + psi^k (e^z - 1),   fixed point psi^n/(2 - e^z),
    //
    // which has a pole at z = ln 2 = 0.693 and returns NEGATIVE psi beyond it;
    // at the default three passes and z = 0.7 it returns 4.08 psi^n against
    // the exact 2.01 psi^n. Reading oldTime() makes psi^n exp(z) a FIXED POINT
    // of the loop at the first pass, for any number of passes.
    //
    // This is safe with respect to OpenFOAM's lazy old-time storage: the
    // fvm::ddt(psi) term of the SAME matrix reads the same psi.oldTime(), so
    // whichever of the two is evaluated first triggers the identical creation
    // / storeOldTimes() refresh of psi_0, and the source and the time
    // derivative can never disagree about what psi^n is.
    const scalarField& psi0 = psi.oldTime().primitiveField();

    // OWNING copy, as everywhere in this hierarchy: a non-owning tmp bound to
    // a reference into a temporary field is the use-after-free that was live
    // in explicitDiscretization::Sc.
    tmp<scalarField> tSc(new scalarField(nonLinearPart.size(), 0.0));
    scalarField& ScField = tSc.ref();

    forAll(nonLinearPart, cellID)
    {
        const scalar a = nonLinearPart[cellID];       // f_nl, units 1/s
        const scalar z = a*dt;                        // dimensionless exponent

        if (Foam::mag(z) <= zMax)
        {
            // Sc = psi^n (e^z - 1)/dt, written as psi^n a phi1(z) so that the
            // z -> 0 limit is exactly a*psi^n and no division by dt round-trips
            // through a*dt/dt.
            ScField[cellID] = psi0[cellID]*a*phi1(z);
        }
        else
        {
            // Clamped branch: an already-diverged strain. Keep the source
            // finite and keep its SIGN, so the failure is a solver residual
            // rather than a NaN. exp(+30) = 1.07e13, exp(-30) = 9.4e-14 (the
            // negative clamp is inert to 1e-13 against the exact limit 0).
            const scalar zc = (z > 0 ? zMax : -zMax);
            ScField[cellID] = psi0[cellID]*(Foam::exp(zc) - 1.0)/dt;
        }
    }

    return tSc;
}

Foam::tmp<scalarField>
Foam::exponentialDiscretization::
Sp(const volScalarField& nonLinearPart) const
{
    // Sp = 0 BY CONSTRUCTION, not by omission.
    //
    // A diagonal coefficient can only realise the rational amplification
    // factor 1/(1 - Sp*dt); the exponential factor exp(f_nl*dt) is not of that
    // form for any Sp, so carrying any part of the term implicitly would
    // destroy the exactness this member exists to provide. The whole term goes
    // to Sc: an EXPLICIT delivery of an EXACT factor, trading the Patankar
    // Sp <= 0 diagonal support for e^{+z}e^{-z} = 1 on the sign-oscillating
    // parasitic strain. The M-matrix property of the psi system must then come
    // from the transport operator alone -- see the class documentation.
    return tmp<scalarField>(new scalarField(nonLinearPart.size(), 0.0));
}

// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
