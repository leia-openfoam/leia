/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    Part of OpenFOAM; GPL v3 or later. See exponentialImplicit.H.
\*---------------------------------------------------------------------------*/

#include "exponentialImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include <cmath>

namespace Foam
{

// Overflow clamp on the ODE exponent z = f_nl*dt, mirroring
// exponentialDiscretization so the two members cannot disagree about what a
// diverged velocity field does.
//
// ASYMMETRIC IN EFFECT, unlike the explicit member. For z -> +inf this
// coefficient SATURATES on its own: S_p dt = 1 - e^{-z} -> 1, so nothing
// overflows and the clamp is inert. It is only the z -> -inf side that needs
// it, where e^{-z} grows without bound. exp(30) = 1.07e13 leaves the assembled
// coefficient over 290 decades below double-precision overflow, and the
// MEASURED exponent on the coupled N=64 stationary droplet is z ~ 0.023 --
// 1300x smaller than the clamp. A run that reaches it has an already-diverged
// velocity field, and the clamp turns what would be inf/NaN silently poisoning
// alpha, U and p through the phase indicator into a large but FINITE
// coefficient that fails visibly in the psi linear solver.
static const scalar zMax = 30.0;

defineTypeNameAndDebug(exponentialImplicit, false);
addToRunTimeSelectionTable
(
    discretization, exponentialImplicit, Dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

exponentialImplicit::exponentialImplicit()
:
    discretization()
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

tmp<scalarField>
exponentialImplicit::Sc
(
    const volScalarField& nonLinearPart,
    const volScalarField& psi
) const
{
    // Identically zero: the entire source is carried by Sp.
    //
    // This is what removes the Picard hazard the explicit member has to work
    // around. With S_c = 0 no psi ITERATE is ever multiplied by a factor, so
    // re-assembling inside the nDefCorr / PIMPLE outer loop cannot build up the
    // fixed point psi^n/(2 - e^z) and its pole at z = ln 2. Nothing here reads
    // psi.oldTime(), and nothing needs to.
    //
    // OWNING copy, as everywhere in this hierarchy: a non-owning tmp bound to a
    // reference into a temporary field is the use-after-free that was once live
    // in explicitDiscretization::Sc.
    return tmp<scalarField>(new scalarField(nonLinearPart.size(), 0.0));
}


tmp<scalarField>
exponentialImplicit::Sp(const volScalarField& nonLinearPart) const
{
    // dt from the mesh's time registry: like the explicit exponential, and
    // unlike the other three members, this coefficient is a function of the
    // STEP and not of f_nl alone.
    const scalar dt = nonLinearPart.mesh().time().deltaTValue();

    tmp<scalarField> tSp(new scalarField(nonLinearPart.size(), 0.0));
    scalarField& SpField = tSp.ref();

    const scalarField& f = nonLinearPart.primitiveField();

    forAll(SpField, cellI)
    {
        // S_p = (1 - exp(-z))/dt with z = f_nl dt, so that
        //     1 - S_p dt = exp(-z)  and  psi^{n+1} = psi^n exp(z)
        // exactly when transport is absent.
        const scalar z = Foam::max(Foam::min(f[cellI]*dt, zMax), -zMax);

        // expm1(-z) = exp(-z) - 1 to full precision near z = 0, where the naive
        // form loses every significant digit to cancellation: at the measured
        // droplet exponent z ~ 0.023, 1 - exp(-z) evaluates a difference of two
        // numbers agreeing to two decimal places. -expm1(-z) is exact there and
        // reduces to f_nl dt as it must, so this member degenerates to
        // simpleLinearImplicit in the small-step limit rather than to noise.
        SpField[cellI] = -std::expm1(-z)/dt;
    }

    return tSp;
}

} // End namespace Foam

// ************************************************************************* //
