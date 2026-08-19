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
    Foam::sdplsRdivStrictSp

Description
    The conservative (divergence-form) SDPLS source of sdplsRdiv, with its
    ONE algebraic coefficient sign-split in the sense of Patankar (1980)
    \cite patankar_numerical_1980 so that the coefficient can never reduce the
    diagonal of the assembled level-set matrix.

    THIS IS A SINGLE-VARIABLE EXPERIMENT. The two divergence terms
    \f$ \nabla\cdot(\phi_w\psi) \f$ and \f$ \nabla\cdot(\phi\psi) \f$ are
    assembled by exactly the same calls, on exactly the same fluxes, with
    exactly the same fvSchemes entries (\c div(phiW,psi) and \c div(phi,psi))
    as sdplsRdiv. The normal \f$ \nhat \f$ is built from the same gradPsi
    strategy, lagged the same way, inside the same defect-correction loop. The
    strain \f$ a \f$ is the same \c strain() the base class refreshes in
    \c update(). Nothing is added, removed, damped, limited or filtered. The
    only difference between this class and sdplsRdiv is HOW the algebraic
    coefficient is delivered to the matrix: fully implicit there, sign-split
    here. sdplsRdiv itself is left untouched, so its committed measurements
    remain valid and the comparison of the two arms is a comparison of one
    thing.

    WHAT sdplsRdiv ASSEMBLES. With \f$ \mathbf{w} = (\nhat\cdot\mathbf{v})\nhat
    \f$, \f$ \phi_w = \mathrm{flux}(\mathbf{w}) \f$ and
    \f$ a = \nhat\cdot(\nabla\mathbf{v})\cdot\nhat \f$ the interfacial strain,
    sdplsRdiv::fvmsdplsSource returns the matrix

    \f{equation}{
        S_{\mathrm{matrix}} =
            \nabla\cdot(\phi_w\psi) - \nabla\cdot(\phi\psi)
            - \underbrace{(\nabla\cdot\mathbf{w} - a)}_{\textstyle -f}\,\psi ,
        \qquad
        f \;\equiv\; a - \nabla\cdot\mathbf{w}
        \quad [\mathrm{s}^{-1}] ,
    \f}

    i.e. the algebraic part it delivers is \f$ +f\,\psi^{n+1} \f$, entirely at
    the new time level, through \c fvm::Sp.

    WHICH SIGN OF f IS UNSAFE -- DERIVED, NOT ASSUMED. OpenFOAM's
    \c fvm::Sp(sp,psi) sets \c diag += V*sp, and this matrix is the RIGHT-hand
    side of the solver's equation (eulerianAdvection.C):

    \f{equation}{
        \underbrace{\mathrm{ddt}(\psi) + \nabla\cdot(\phi\psi)
        - (\nabla\cdot\phi)\,\psi}_{\textstyle A}
        \;=\; S_{\mathrm{matrix}} ,
    \f}

    and \c fvMatrix::operator== forms \f$ A - S_{\mathrm{matrix}} \f$ termwise,
    so the assembled system diagonal is
    \f$ \mathrm{diag}(A) - \mathrm{diag}(S_{\mathrm{matrix}}) \f$. With
    \f$ \mathrm{diag}(S_{\mathrm{matrix}}) = +V f \f$ from the \c fvm::Sp term,
    a POSITIVE \f$ f \f$ SUBTRACTS \f$ V f \f$ from the system diagonal and a
    negative \f$ f \f$ adds to it. In the language of the linearization
    \f$ S = S_c + S_p\,\psi^{n+1} \f$ used by the discretization hierarchy, the
    effective coefficient of sdplsRdiv is \f$ S_p = f = a -
    \nabla\cdot\mathbf{w} \f$, and Patankar's requirement is \f$ S_p \le 0 \f$.
    Written in the coefficient that appears literally in the sdplsRdiv source
    line, \f$ c = \nabla\cdot\mathbf{w} - a = -f \f$, the unsafe branch is
    \f$ c < 0 \f$.

    One cell, transport switched off, backward Euler, cell volume \f$ V \f$,
    makes the consequence explicit. Fully implicit (sdplsRdiv):

    \f{equation}{
        \left(\frac{1}{\Delta t} - f\right)\psi^{n+1}
        = \frac{\psi^{n}}{\Delta t} .
    \f}

    The diagonal \f$ 1/\Delta t - f \f$ passes through zero at
    \f$ f = 1/\Delta t \f$ and changes sign beyond it: the update is then
    amplifying and sign-reversing, and the linear solver still CONVERGES -- it
    is solving a well-posed but wrong-signed system. Sign-split (this class):

    \f{equation}{
        \left(\frac{1}{\Delta t} - f^{-}\right)\psi^{n+1}
        = \frac{\psi^{n}}{\Delta t} + f^{+}\psi^{n},
        \qquad f^{+} = \max(f,0) \ge 0, \quad f^{-} = \min(f,0) \le 0,
    \f}

    whose diagonal \f$ 1/\Delta t - f^{-} \ge 1/\Delta t \f$ for every cell,
    every step, every mesh. The positive branch is still delivered -- the
    source is not weakened -- but it is delivered to the right-hand side, where
    it cannot invert the diagonal.

    WHY THE FAILURE DOES NOT WASH OUT UNDER REFINEMENT (argument, not
    measurement). \f$ \nabla\cdot\mathbf{w} \f$ contains \f$ \nabla\nhat \f$,
    i.e. SECOND derivatives of \f$ \psi \f$; on a mesh whose discrete noise in
    a first derivative is \f$ O(\varepsilon/h) \f$, the noise carried in
    \f$ f \f$ is \f$ O(\varepsilon/h^{2}) \f$ relative to nothing that
    decreases. Meanwhile a fixed-CFL ladder shrinks \f$ \Delta t \propto h \f$,
    so the stability threshold \f$ 1/\Delta t \f$ grows only like \f$ 1/h \f$.
    The ratio \f$ f\,\Delta t \f$ therefore does not decrease under refinement
    and can grow -- which is consistent with a failure that gets WORSE on finer
    meshes rather than better.

    THE MEASUREMENT THIS CLASS EXISTS TO EXPLAIN, recorded verbatim.
    2D reversed vortex, \f$ T = 2 \f$, PRESCRIBED divergence-free velocity, no
    capillary force, no feedback of any kind: purely kinematic, so nothing here
    can be blamed on a coupling. Reversed-flow protocol: the gradient error is
    read at \f$ t = T/2 \f$ (maximal deformation, before reversal cancellation
    corrupts the \f$ t = T \f$ reading), the shape error at \f$ t = T \f$, the
    volume-conservation error at both. Band gradient error at \f$ t = T/2 \f$,
    \f$ N = 32, 45, 64, 90, 128, 181, 256 \f$:

    \verbatim
        noSource  1.222  1.529  1.637  1.666  1.675  1.663  1.673   order -0.116
        R         0.344  0.240  0.172  0.135  0.0958 0.0536 0.0309  order +1.112
        Rdiv      0.391  0.256  1.564  8.002  20.48  103.3  668.9   order -3.800
    \endverbatim

    Shape error order at \f$ t = T \f$: noSource +2.134, R +1.883, Rdiv -1.306.
    Volume-conservation error order: noSource +2.733, R +2.090, Rdiv -1.544.
    (The per-rung shape and volume VALUES are not reproduced here, only the
    fitted orders that were reported; the per-rung gradient values are the
    three rows above.) Rdiv's relative volume error reaches 1.36, i.e. 136 %,
    at \f$ N = 256 \f$.

    Read the three metrics together, as the standing rule requires: all three
    orders are negative for Rdiv, so this is not one metric misbehaving while
    the solution is otherwise sound -- the gradient, the shape and the enclosed
    volume all degrade together, and all three degrade FASTER on finer meshes.
    Rdiv is competitive with R only on the two coarsest rungs (0.391 vs 0.344
    at \f$ N = 32 \f$, 0.256 vs 0.240 at \f$ N = 45 \f$) and departs from
    \f$ N = 64 \f$ onwards, which is the signature of a threshold being crossed
    rather than of a consistently wrong operator.

    EVERY ONE OF THOSE RUNS COMPLETED. The \c solverFailed column is blank at
    all seven resolutions. These are converged linear solves producing
    increasingly wrong answers, not crashes -- exactly the failure mode a lost
    M-matrix property produces, and exactly the failure mode that a residual
    check or a "did it run" gate cannot detect.

    WHAT THIS CLASS ASSEMBLES.

    \f{equation}{
        S_{\mathrm{matrix}} =
            \nabla\cdot(\phi_w\psi) - \nabla\cdot(\phi\psi)
            + \max(f,0)\,\psi^{n} + \min(f,0)\,\psi^{n+1},
        \qquad f = a - \nabla\cdot\mathbf{w},
    \f}

    with \f$ \max(f,0)\psi^{n} \f$ placed through \c fvm::Su (right-hand side
    only) and \f$ \min(f,0)\psi^{n+1} \f$ through \c fvm::Sp (diagonal only).
    This is the decomposition of strictNegativeSpLinearImplicit, term for term,
    applied to the coefficient \f$ f \f$ of the divergence form instead of to
    \f$ f_{nl} \f$ of the algebraic form. At the continuum
    \f$ \max(f,0) + \min(f,0) = f \f$ identically, so the modified equation is
    the same equation; only the time level at which the positive branch is
    evaluated differs.

    WHICH \f$ \psi \f$ THE EXPLICIT BRANCH USES: \c psi.oldTime(), i.e.
    \f$ \psi^{n} \f$, NOT the current iterate. The reason is specific to this
    solver. The Eulerian advection loop re-assembles and re-solves \c nDefCorr
    times per step (default 3) with \c psi updated in place between passes, to
    converge the deferred \c linearUpwind correction. If the explicit branch
    were evaluated at the CURRENT \c psi, then at convergence of that loop
    \f$ \max(f,0)\psi = \max(f,0)\psi^{n+1} \f$ and the split would collapse
    back onto the fully implicit sdplsRdiv equation: the fixed point of the two
    schemes would be identical, the experiment would be guaranteed to reproduce
    Rdiv's answer whenever the loop converges, and a null result would be
    uninformative -- it could not distinguish "the unsigned coefficient was not
    the defect" from "the stabilisation was annihilated by the outer loop".
    With \c psi.oldTime() the positive branch stays at time level \f$ n \f$ for
    the whole step, so the split survives the loop and the scheme actually
    solved is the sign-split one. Note that
    strictNegativeSpLinearImplicit::Sc is written for \f$ \psi^{n} \f$ (its
    documented formula is \f$ \max(f_{nl},0)\psi^{n}_{c} \f$) and receives the
    current iterate only as a side effect of being called from inside the same
    defect-correction loop; \c psi.oldTime() is that class's stated intent, and
    with \c nDefCorr = 1 the two choices coincide exactly.

    THE PRICE OF THAT CHOICE, stated plainly: where \f$ f > 0 \f$ this variant
    is first-order-in-time explicit in the source coefficient, so it is not
    identical to Rdiv even when Rdiv behaves. On the two coarse rungs where
    Rdiv is still competitive, a small difference from Rdiv is expected and is
    not evidence of anything. The experiment is decided by the trend across the
    ladder, not by any single rung.

    THE DIAGNOSTIC THAT MAKES A NULL RESULT READABLE. Each assembly prints one
    line:

    \verbatim
        sdplsRdivStrictSp: f = a - div(w); cells with f > 0 (explicit branch
        active): <n>/<N>, max(f) dt = <x>, min(f) dt = <y>
    \endverbatim

    It exists because the split is a no-op wherever \f$ f \le 0 \f$: there
    \f$ \min(f,0) = f \f$, \f$ \max(f,0) = 0 \f$, and this class assembles
    exactly what sdplsRdiv assembles. If the run reproduces Rdiv's divergence
    while \c n is 0 or \c max(f)dt is far below 1, the conclusion is "the
    coefficient was already sign-safe on this problem, the hypothesis was never
    tested", which is a DIFFERENT statement from "the diagonal was not the
    defect" and implies a different next step. \f$ f\,\Delta t \f$ is the
    quantity to look at because the one-cell backward-Euler diagonal is
    \f$ 1/\Delta t - f \f$ per unit volume: \f$ \max(f)\Delta t \ge 1 \f$
    is a cell whose fully implicit diagonal has been driven through zero.
    Report this count alongside the three error metrics.

    WHAT THE EXPERIMENT DECIDES.

      - Converges kinematically -- band gradient error at \f$ t = T/2 \f$
        decreasing with a positive order comparable to R's +1.112, shape order
        at \f$ t = T \f$ positive, volume-conservation order positive: then the
        unsigned coefficient WAS the defect. The conservative form is sound and
        rescuable, and the sign-split delivery is what it needs; sdplsRdiv
        should then be replaced by this assembly rather than carried alongside
        it.

      - Still diverges -- any of the three orders still negative, or the
        gradient error still growing past \f$ N = 64 \f$: then the defect is
        NOT the diagonal, because after this change the diagonal is provably
        never reduced by the source. The remaining candidates are, in order of
        suspicion: (i) the \f$ \nabla\nhat \f$ content of
        \f$ \nabla\cdot\mathbf{w} \f$, whose second-derivative noise is the
        term that grows under refinement; (ii) the near-cancellation of
        \f$ \nabla\cdot(\phi_w\psi) - \nabla\cdot(\phi\psi) \f$, which is a
        difference of two large numbers whose exact cancellation is a continuum
        identity the discretization does not inherit -- the relative error of
        the difference is the relative error of the operands divided by the
        (small) difference. In that case the conservative form should be
        RETIRED, not carried as a half-tested option: two independent
        stabilisations of the same form, both failing, is a statement about the
        form.

    Either outcome is a result and must be reported as one. A negative result
    here closes a branch of the SDPLS conservation argument that has so far
    been argued structurally ("what enters through fluxes cannot drain") and
    never verified kinematically.

    MOLLIFIER: REFUSED, BY INHERITANCE. The constructor delegates to
    sdplsRdiv's, which raises a FatalError for any \c mollifier type other than
    \c none, because a cut-off multiplying a FLUX destroys the conservation the
    divergence form exists to provide. The refusal is therefore inherited, not
    replicated -- there is one copy of that rule, in sdplsRdiv.

    DIAGNOSTIC CAVEAT, inherited from sdplsRdiv and repeated here because it
    affects what the CSV columns mean: the registered field
    \c SDPLS_nonLinearPart is filled by \c update() with the ALGEBRAIC
    coefficient \f$ a \f$, not with \f$ f = a - \nabla\cdot\mathbf{w} \f$ and
    not with the fluxes. Any budget built from
    \c gSum(SDPLS_nonLinearPart*psi*V) is therefore NOT the source that ran.
    psiConservationCSV recognises this source type as non-algebraic and NaNs
    those columns.

    Selected as:
    \verbatim
        sdplsSource
        {
            type            RdivStrictSp;
            gradPsi         fvc;      // the normal's source, as for R and Rdiv
            mollifier       { type none; }   // anything else is a FatalError
            // discretization is IGNORED: the assembly is fixed by the form.
        }
    \endverbatim

SourceFiles
    sdplsRdivStrictSp.C

\*---------------------------------------------------------------------------*/

#ifndef sdplsRdivStrictSp_H
#define sdplsRdivStrictSp_H

#include "sdplsRdiv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sdplsRdivStrictSp Declaration
\*---------------------------------------------------------------------------*/

class sdplsRdivStrictSp
:
    public sdplsRdiv
{

public:

    // Static Data Members

        TypeName("RdivStrictSp");

    // Constructors

        //- Construct from components. Delegates to sdplsRdiv, which is what
        //  makes the mollifier refusal inherited rather than copied.
        sdplsRdivStrictSp(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~sdplsRdivStrictSp() = default;

    // Member Functions

        //- The conservative assembly with the algebraic coefficient
        //  f = a - div(w) sign-split: max(f,0) explicit against psi^n,
        //  min(f,0) implicit. The two divergence terms are identical to
        //  sdplsRdiv's. discretization() is unused, as for sdplsRdiv.
        virtual tmp<fvScalarMatrix>
        fvmsdplsSource(const volScalarField& psi, const volVectorField& U);
};

} // End namespace Foam

#endif

// ************************************************************************* //
