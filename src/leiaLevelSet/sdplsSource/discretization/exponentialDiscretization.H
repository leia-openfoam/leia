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
    Foam::exponentialDiscretization

Description
    Exponential (integrating-factor) linearization of the SDPLS source: the
    EXACT one-step amplification factor, delivered explicitly.

    With the coefficient frozen at the old time level -- which is what every
    member of this hierarchy does, since \f$ f_{nl} \f$ is evaluated from
    \f$ \psi^n \f$ and \f$ \mathbf{v} \f$ before the matrix is assembled -- the
    SDPLS source is the LINEAR scalar ODE
    \f{equation}{
        \frac{d\psi}{dt} = f_{nl}(\psi^n)\,\psi,
        \qquad a \equiv f_{nl}(\psi^n),
    \f}
    whose exact solution over one step is the MULTIPLICATIVE integrating factor
    \f{equation}{
        \psi^{n+1} = \psi^n\,G(z), \quad G(z) = \exp(z),
        \quad z \equiv a\,\Delta t .
    \f}
    The three existing members are rational approximations of that same
    \f$ G \f$:
    \f{equation}{ \begin{aligned}
        \texttt{explicit}                       &:\; G \approx 1 + z \\
        \texttt{simpleLinearImplicit}           &:\; G \approx 1/(1 - z) \\
        \texttt{strictNegativeSpLinearImplicit} &:\; G \approx 1 + z \;(z > 0),
                                                   \;\; 1/(1 - z) \;(z < 0),
    \end{aligned} \f}
    each first order accurate in \f$ z \f$ and each carrying a \f$ |z| \f$
    restriction the exact factor does not have (\f$ 1 + z \f$ changes SIGN at
    \f$ z = -1 \f$; \f$ 1/(1 - z) \f$ has a POLE at \f$ z = 1 \f$, while
    \f$ \exp(z) > 0 \f$ and finite for every real \f$ z \f$).

    This class delivers \f$ G = \exp(z) \f$ itself. The linearization constants
    are
    \f{equation}{ \begin{aligned}
        S_{exponential,\,c} &= \psi^n_c\,
                               \frac{\exp(z_c) - 1}{\Delta t},
                               \qquad z_c = a_c\,\Delta t \\
        S_{exponential,\,p} &= 0,
    \end{aligned} \f}
    so the source term is
    \f{equation}{
        S_{exponential}(\psi^{n+1}) =
            \psi^n_c\,\frac{\exp(f_{nl}(\psi^n)\,\Delta t) - 1}{\Delta t},
    \f}
    and, paired with an implicit-Euler \f$ \partial_t \f$, the assembled step
    \f$ (\psi^{n+1} - \psi^n)/\Delta t = S_c \f$ returns
    \f$ \psi^{n+1} = \psi^n \exp(z) \f$ EXACTLY, at any \f$ \Delta t \f$.

    WHY, MEASURED. Under Navier--Stokes coupling the stationary droplet
    (cases/stationaryDroplet2D) develops a mesh-locked four-lobed (mode-4)
    capillary parasitic current. The normal strain the source reads from that
    current, \f$ a = \nhat\cdot\nabla\mathbf{v}\cdot\nhat \f$, OSCILLATES IN
    SIGN at the parasitic frequency with a near-zero mean, so what matters is
    not the per-step accuracy of \f$ G \f$ but how a \f$ +z \f$ step and a
    \f$ -z \f$ step COMPOSE:

    \verbatim
        exact           exp(+z) exp(-z)                    = 1
        Crank-Nicolson  [(1+z/2)/(1-z/2)][(1-z/2)/(1+z/2)] = 1
        strictNegSp     (1 + z) [1/(1 + z)]                = 1
        explicit        (1 + z)(1 - z)   = 1 - z^2   one-sided DECAY  of |psi|
        implicit Euler  1/[(1-z)(1+z)]   = 1/(1-z^2)
                                         = 1 + z^2 + O(z^4)
                                                     one-sided GROWTH of |psi|
    \endverbatim

    So \c simpleLinearImplicit RECTIFIES a zero-mean oscillating strain into
    one-sided growth of \f$ |\psi| \f$ at \f$ O(z^2) \f$ per pair, and
    \c explicit into one-sided decay; the sign-split
    \c strictNegativeSpLinearImplicit happens to be exactly neutral on a
    symmetric pair (its two branches are reciprocals), even though each half is
    only first order. The indictment therefore lands on the scheme the coupled
    droplet studies actually pin (config/sdplsDropletMechanism2D.yaml:
    \c SOURCE_SCHEME \c simpleLinearImplicit).

    The magnitude, from the coupled N=64 run: the amplified parasitic speed
    \f$ |\mathbf{u}| \sim 0.1\,\mathrm{m/s} \f$ over
    \f$ h = 1.56\times10^{-4}\,\mathrm{m} \f$ gives
    \f$ a \sim |\mathbf{u}|/h \sim 640\,\mathrm{s^{-1}} \f$, and the capillary
    step \f$ \Delta t \sim 3.6\times10^{-5}\,\mathrm{s} \f$ gives
    \f$ z \sim 0.023 \f$. About 1400 \f$ \pm \f$ pairs fit into the
    integration, so the accumulated spurious amplification is
    \f$ (1 + z^2)^{1400} = \exp(1400 \times 5.3\times10^{-4})
    = e^{0.74} \approx 2.1 \f$: a doubling of \f$ |\psi| \f$ in the band that
    no physical term put there. Refining \f$ \Delta t \f$ by 4x (11 315 steps)
    moved the MEASURED drain half-life from 0.046--0.049 to 0.052--0.060
    (+~20%), which is the direction and rough size that dropping the exponent
    4x predicts (\f$ e^{0.74} \to e^{0.19} \approx 1.2 \f$). The term is REAL
    but MODEST -- it is not the whole drain. \f$ \exp(z) \f$ removes it
    IDENTICALLY at any \f$ \Delta t \f$, because
    \f$ e^{+z}e^{-z} = 1 \f$ is an algebraic identity and not a truncation.

    WHAT THIS DOES NOT CLAIM. A strictly positive multiplicative factor cannot
    by itself move the zero level set: \f$ \psi G = 0 \iff \psi = 0 \f$ for
    \f$ G > 0 \f$. That is the SDPLS design property and it survives here
    unchanged, so the rectification does NOT displace the interface directly.
    It acts on \f$ \psi \f$ PERTURBATIONS in the band: for a spatially varying
    \f$ G \f$, \f$ \nabla(G\psi) = G\nabla\psi + \psi\nabla G \f$, whose second
    term vanishes ON \f$ \psi = 0 \f$ but NOT at the band cell centres, so
    \f$ |\nabla\psi| \f$ off the interface is inflated, the curvature
    \f$ \kappa = \nabla\cdot(\nabla\psi/|\nabla\psi|) \f$ evaluated at those
    cell centres changes, and the CSF force built from it changes the very
    parasitic current that produced \f$ a \f$. This is a LOOP-GAIN effect on
    the capillary feedback path, not an interface displacement, and it is not
    by itself an explanation of the measured volume drain (the conservation
    failure is addressed structurally by \c sdplsRdiv).

    THE TRADE, stated up front. The exponential factor is MULTIPLICATIVE and
    has no exact representation as a linear \c fvm::Sp coefficient: a diagonal
    contribution \f$ S_p \f$ can only produce the rational factor
    \f$ 1/(1 - S_p\Delta t) \f$. Reproducing \f$ \exp(z) \f$ exactly therefore
    forces the entire term onto \f$ S_c \f$, i.e. \f$ S_p = 0 \f$ and NO
    diagonal support: this member is an EXPLICIT delivery of an EXACT factor.
    It gives up the diagonal dominance that \c strictNegativeSpLinearImplicit
    exists to protect (Patankar's \f$ S_p \le 0 \f$ rule
    \cite patankar_numerical_1980), so the M-matrix property of the assembled
    \f$ \psi \f$ system must come from the transport operator alone. What it
    does NOT give up is boundedness of the factor itself: \f$ \exp(z) \f$ is
    positive and finite for every \f$ z \f$, so unlike \c explicit it cannot
    flip the sign of \f$ \psi \f$ at large \f$ |z| \f$ and unlike
    \c simpleLinearImplicit it has no pole.

    \f$ \psi^n \f$ IS \c psi.oldTime(), NOT THE ARGUMENT \c psi. Every call
    site re-assembles this source inside an outer loop that has ALREADY
    overwritten \c psi: the defect-correction loop of \c eulerianAdvection and
    \c leiaRedistancedLevelSetFoam (\c nDefCorr, default 3) and the PIMPLE
    outer-corrector loop of \c alphaEqn.H. Feeding the latest ITERATE into a
    multiplicative factor turns that loop into the Picard iteration
    \f$ \psi^{k+1} = \psi^n + \psi^k(e^z - 1) \f$, whose fixed point is
    \f$ \psi^n/(2 - e^z) \f$ -- a POLE at \f$ z = \ln 2 = 0.693 \f$ and
    NEGATIVE \f$ \psi \f$ beyond it. At the default three passes and
    \f$ z = 0.7 \f$ that iteration returns \f$ 4.08\,\psi^n \f$ against the
    exact \f$ 2.01\,\psi^n \f$: a factor-two error produced purely by the loop.
    Reading \c psi.oldTime() instead makes \f$ \psi^n e^z \f$ a FIXED POINT of
    the loop already at the first pass, for any number of passes. (For the
    record, the same iteration applied to \c explicit, \f$ S_c = a\psi^k \f$,
    converges to \f$ \psi^n/(1 - z) \f$: three passes of \c explicit deliver
    the truncated geometric series \f$ 1 + z + z^2 + z^3 \f$, not
    \f$ 1 + z \f$, which over-weights the quadratic term of \f$ \exp \f$ by 2x
    and the cubic by 6x.)

    Reading \c psi.oldTime() here is safe by construction: \c fvm::ddt(psi) in
    the SAME matrix reads the same field, so whichever of the two is evaluated
    first triggers the identical lazy creation / \c storeOldTimes() refresh of
    \c psi_0, and the source and the time derivative are consistent by
    construction rather than by evaluation order.

    THE EXACTNESS IS CONDITIONED ON \c ddt(psi) \c Euler. The identity
    \f$ \psi^{n+1} = \psi^n + \Delta t\,S_c = \psi^n e^z \f$ needs
    \c fvm::ddt to contribute \f$ (\psi^{n+1} - \psi^n)/\Delta t \f$.
    cases/stationaryDroplet2D renders \c ddt(psi) from the \c PSI_DDT study
    token, whose default (cases/default.parameter: \c EulerDdt) IS Euler, so
    the coupled droplet studies that motivate this class satisfy it. The
    verification cases (cases/2Dvortex, cases/3Dshear, cases/3Ddeformation) use
    \c default \c CrankNicolson \c 1 and do NOT: there \c fvm::ddt contributes
    \f$ 2(\psi^{n+1} - \psi^n)/\Delta t - (\partial_t\psi)^n \f$, so the step
    realises
    \f$ \psi^{n+1} = \psi^n + (\Delta t/2)(\partial_t\psi)^n
    + \psi^n(e^z - 1)/2 \f$ -- the exponential COEFFICIENT is still delivered,
    but the \f$ \pm z \f$ pair no longer cancels exactly. Select
    \c exponential together with \c ddt(psi) \c Euler when the exact-factor
    property is the point of the run.

    EVALUATION. \f$ S_c \f$ is formed as
    \f$ \psi^n a\,\varphi_1(z) \f$ with
    \f$ \varphi_1(z) = (e^z - 1)/z \f$, the first exponential-integrator
    weight, rather than as a literal \f$ (e^z - 1)/\Delta t \f$. The two are
    identical in exact arithmetic; the \f$ \varphi_1 \f$ form makes the
    \f$ z \to 0 \f$ limit exactly \f$ a\,\psi^n \f$ in floating point, so this
    member degenerates CONTINUOUSLY to the other three linearizations as
    \f$ z \to 0 \f$ instead of losing digits to the cancellation
    \f$ e^z - 1 \f$ (whose relative error is \f$ \sim \epsilon/|z| \f$).

    OVERFLOW CLAMP. \f$ z \f$ is clamped to \f$ |z| \le 30 \f$ before
    \f$ \exp \f$. \f$ e^{30} = 1.07\times10^{13} \f$; times a band value
    \f$ |\psi| \sim 10^{-3}\,\mathrm{m} \f$ over
    \f$ \Delta t \sim 3.6\times10^{-5}\,\mathrm{s} \f$ the worst
    \f$ S_c \f$ is \f$ \sim 3\times10^{14} \f$, and the matrix entry
    \f$ V S_c \f$ with \f$ V \sim h^3 \sim 4\times10^{-12}\,\mathrm{m^3} \f$ is
    \f$ \sim 10^{3} \f$ -- over 300 decades below the double-precision
    overflow at \f$ 1.8\times10^{308} \f$, so no product in the assembly chain
    can overflow. The clamp CANNOT be reached by a physical run: for
    \c type \c R the strain scales with the velocity gradient while
    \f$ \Delta t \f$ is Courant limited, so \f$ |z| = O(Co) = O(1) \f$, and the
    MEASURED worst case above is \f$ z \sim 0.023 \f$, 1300x below the clamp;
    for \c type \c beta, \f$ f_{nl} = 1\,[\mathrm{s^{-1}}]
    (\beta - |\nabla\psi|) \f$ with \f$ \beta = O(1) \f$ so
    \f$ |z| \ll 1 \f$. Reaching 30 at
    \f$ \Delta t = 3.6\times10^{-5}\,\mathrm{s} \f$ needs
    \f$ |a| \sim 8.3\times10^{5}\,\mathrm{s^{-1}} \f$, three orders beyond the
    parasitic strain -- only an ALREADY DIVERGED velocity field produces it,
    and the clamp then yields a large but FINITE source, so the failure
    surfaces as a \f$ \psi \f$ linear-solver residual instead of a NaN
    propagated silently into alpha, U and p. The negative clamp is inert to
    \f$ 10^{-13} \f$: \f$ e^{-30} = 9.4\times10^{-14} \f$ against the exact
    limit 0.

    COVERAGE. leiaTestSdplsSource carries \c exponential in its \c discs list
    with its OWN contract: at \f$ \psi = \psi^n \f$ the applied source is
    \f$ \psi(e^{f\Delta t} - 1)/\Delta t \f$, NOT \f$ f\psi \f$ (they differ at
    \f$ O(f^2\Delta t) \f$, i.e. \f$ 2.5\times10^{-3} \f$ at \f$ f = 0.7 \f$,
    \f$ \Delta t = 0.01 \f$ -- seven orders above the \f$ 10^{-10} \f$
    tolerance), so it is deliberately EXEMPT from the cross-discretization
    agreement check, and its one-step closed form is \f$ \psi e^{f\Delta t} \f$
    rather than \f$ \psi(1 + f\Delta t) \f$ or
    \f$ \psi/(1 - f\Delta t) \f$. It passes the first-order-consistency layer
    with zero deviation, being the exact solution.

    Selected as:
    \verbatim
        sdplsSource
        {
            type            R;
            discretization  exponential;
        }
    \endverbatim
    Not reachable unless a dictionary asks for it: \c sdplsSource defaults
    \c discretization to \c none, and every config/*.yaml pins
    \c SOURCE_SCHEME to one of the three pre-existing names, so every existing
    study is bit-for-bit unchanged by the presence of this member.

SourceFiles
    exponentialDiscretization.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialDiscretization_H
#define exponentialDiscretization_H

#include "discretization.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class exponentialDiscretization Declaration
\*---------------------------------------------------------------------------*/


class exponentialDiscretization
:
    public discretization
{
public:

    // Static Data Members

        // The dictionary name is the short "exponential"; the class carries
        // the "Discretization" suffix only for symmetry with
        // explicitDiscretization, which needs it because `explicit` is a C++
        // keyword.
        TypeName("exponential");

    // Constructors

        //- Default construct
        exponentialDiscretization();

    //- Destructor
    virtual ~exponentialDiscretization() = default;

protected:

    // Protected Member Functions

        //- Calculates the constant of the source term linearization.
        virtual tmp<scalarField> Sc
            (
                const volScalarField& nonLinearPart,
                const volScalarField& psi
            ) const;

        //- Calculates the proportional constant of the source term
        //  linearization.
        virtual tmp<scalarField> Sp(const volScalarField& nonLinearPart) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
