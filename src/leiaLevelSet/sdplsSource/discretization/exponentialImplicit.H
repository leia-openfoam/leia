/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    Part of OpenFOAM; GPL v3 or later.

Class
    Foam::exponentialImplicit

Description
    The EXACT integrating factor of the SDPLS source, delivered IMPLICITLY.

    The source is d(psi)/dt = f_nl(psi^n) psi, a linear ODE in psi whose exact
    solution over one step is the multiplicative factor psi^{n+1} = psi^n e^z
    with z = f_nl dt. Writing that as an implicit linearization,

        (psi^{n+1} - psi^n)/dt = S_p psi^{n+1}   and   psi^{n+1} = psi^n e^z
        =>  1 - S_p dt = e^{-z}
        =>  S_p = (1 - e^{-z})/dt ,   S_c = 0.

    NOTE the coefficient is NOT (e^z - 1)/dt. That is the EXPLICIT delivery
    (exponentialDiscretization); used implicitly it would give
    psi^{n+1} = psi^n/(2 - e^z), which has a pole at z = ln 2 and returns a
    NEGATIVE psi beyond it.

WHY THIS EXISTS: the explicit exponential costs shape and volume order
-----------------------------------------------------------------------
exponentialDiscretization reproduces exp(z) exactly, but it does so with
S_p = 0 -- the whole source is explicit -- so it gives up the implicit coupling
between the source and the transport operator. MEASURED on the 2D reversed
vortex, T = 2, seven rungs N = 32..256 (study sdplsExpSource2DvortexRev):

    arm                        gradBand(T/2)  shape(T)  volume(T)  volume(T/2)
    R / simpleLinearImplicit      +1.112       +1.883    +2.090     +3.294
    R / exponential               +1.109       +1.256    +1.718     +2.398

The gradient order is unchanged, confirming that a TEMPORAL correction cannot
repair a SPATIAL order loss -- but shape falls +1.883 -> +1.256 and volume
+2.090 -> +1.718. It is exactly the zero-contour metrics that record the loss,
which is what giving up the implicit coupling predicts.

This member keeps the exact factor AND the implicit coupling.

THREE PROPERTIES, all provable rather than hoped for
----------------------------------------------------
1. EXACT IN TIME. With transport absent the update is psi^n e^z identically, at
   any dt, for any sign of f_nl.

2. THE DIAGONAL NEVER CHANGES SIGN. The assembled diagonal carries
   1/dt - S_p = e^{-z}/dt, which is strictly positive for every finite z.
   simpleLinearImplicit carries 1/dt - f_nl, which goes NEGATIVE as soon as
   f_nl dt > 1 -- the |f_nl| dt < 1 restriction the SDPLS paper accepts. This
   member removes that restriction by construction, not by clipping.

3. NO PICARD POLE, and this is the practical difference from the explicit
   member. All three call sites re-assemble the source inside an outer loop
   (nDefCorr, default 3, or the PIMPLE correctors) that has already overwritten
   psi. exponentialDiscretization multiplies a psi ITERATE, so that loop becomes
   psi^{k+1} = psi^n + psi^k(e^z - 1), whose fixed point psi^n/(2 - e^z) has a
   pole at z = ln 2 = 0.693 -- it is safe only because it deliberately reads
   psi.oldTime(). Here S_c = 0 and S_p multiplies psi^{n+1} through the matrix,
   so no iterate is ever multiplied and the hazard cannot arise.

A strictly positive multiplicative factor also cannot flip the sign of psi in
any cell, so the discrete zero contour's topology is preserved at any dt -- the
sign-preservation property the exponential family exists for.

SourceFiles
    exponentialImplicit.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialImplicit_H
#define exponentialImplicit_H

#include "discretization.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class exponentialImplicit Declaration
\*---------------------------------------------------------------------------*/

class exponentialImplicit
:
    public discretization
{
public:

    TypeName("exponentialImplicit");

    exponentialImplicit();

    virtual ~exponentialImplicit() = default;

    //- Zero: the whole source is carried implicitly by Sp.
    virtual tmp<scalarField> Sc
    (
        const volScalarField& nonLinearPart,
        const volScalarField& psi
    ) const;

    //- (1 - exp(-f_nl dt))/dt
    virtual tmp<scalarField> Sp(const volScalarField& nonLinearPart) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
