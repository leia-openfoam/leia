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
    Foam::newtonShiftVolumeCorrection

Description
    Distance-preserving GLOBAL volume correction: the uniform level set shift
    of Song et al., in which the whole field is displaced by a single scalar

    \f{equation}{
        \psi \;\leftarrow\; \psi + \varepsilon ,
    \f}

    with \f$ \varepsilon \in \mathbb{R} \f$ (one number for the entire mesh)
    chosen as the root of

    \f{equation}{
        f(\varepsilon) \;=\; V(\varepsilon) - V_{target} \;=\; 0,
        \qquad
        V(\varepsilon) \;=\; \sum_c \alpha_c(\psi + \varepsilon)\,|\Omega_c| ,
    \f}

    where \f$ \alpha_c \f$ is produced by the run-time selected
    Foam::phaseIndicator and \f$ V_{target} \f$ was captured once by
    Foam::volumeCorrection::setTargetVolume(). The root is found by
    Newton-Raphson

    \f{equation}{
        \varepsilon_{k+1} \;=\; \varepsilon_k
        - \frac{f(\varepsilon_k)}{f'(\varepsilon_k)},
        \qquad
        \varepsilon_0 = 0 ,
    \f}

    safeguarded by bisection on the bracket
    \f$ [-\varepsilon_{max}, +\varepsilon_{max}] \f$,
    \f$ \varepsilon_{max} = \f$ \c maxShiftCells \f$ \cdot h \f$.

    WHY IT PRESERVES THE SIGNED-DISTANCE PROPERTY
    \f$ \varepsilon \f$ is a CONSTANT field, so
    \f$ \nabla(\psi + \varepsilon) = \nabla \psi \f$ -- and not only in the
    continuum. The discrete Gauss gradient of a cell is
    \f$ (1/|\Omega_c|)\sum_f \mathbf{S}_f \psi_f \f$ and
    \f$ \sum_f \mathbf{S}_f = \mathbf{0} \f$ for a closed cell; the
    least-squares gradient is a linear combination of the differences
    \f$ \psi_N - \psi_P \f$. Both are invariant under a uniform shift
    BIT-FOR-BIT, provided the boundary values are shifted with the interior --
    which zeroGradient (the psi condition of every case in \c cases/) does on
    the following correctBoundaryConditions(), and which a \c fixedValue psi
    patch would NOT do. correct() therefore checks the psi boundary types once
    and warns by name if any patch fixes its value.

    Three consequences follow, and they are the whole argument for the method:
      - \f$ ||\nabla\psi| - 1| \f$, the metric this library uses to decide
        whether \f$ \psi \f$ is still a signed distance, is unchanged by the
        correction. The correction cannot buy volume by destroying the
        distance property, which is what a multiplicative or local
        (\f$ \varepsilon(\mathbf{x}) \f$) rescaling does.
      - Neither redistancing trigger can be provoked by it: both
        \c signedDistanceBounds (redistancer.C:336-374) and
        \c gradPsiThreshold (redistancer.C:385-398) are functions of
        \f$ |\nabla\psi| \f$ ONLY. The correction is therefore applied AFTER
        redistancing and cannot cause an extra redistancing event on the next
        step.
      - The shift restores the volume EXACTLY whatever \f$ |\nabla\psi| \f$
        happens to be, because \f$ V(\varepsilon) \f$ is MEASURED through the
        indicator rather than modelled. The interpretation of
        \f$ \varepsilon \f$ as an interface displacement normal to
        \f$ \Gamma \f$ additionally requires \f$ |\nabla\psi| \approx 1 \f$;
        that is what the SDPLS source and the redistancer maintain, and what
        the shift leaves untouched.

    KNOWN DRAWBACK, WHICH IS SEVERE
    One global \f$ \varepsilon \f$ conserves the TOTAL phase volume and
    nothing else.
      - It is ILL-POSED for several disconnected bodies. With bodies
        \f$ B_1,\dots,B_m \f$ the honest problem is m equations
        \f$ V_i(\varepsilon_i) = V_{i,target} \f$ in m unknowns; a single
        \f$ \varepsilon \f$ solves their SUM. A run in which one droplet is
        losing volume and another is gaining it reports zero total volume
        error while both bodies are wrong, and the shift actively makes one of
        them worse -- it inflates every body when the total is short,
        including the bodies that were already too large.
      - Even for ONE body it trades local interface position for a global
        number: the displacement \f$ \varepsilon \f$ is applied uniformly
        along \f$ \Gamma \f$, whereas the volume was lost where the flow
        strained the interface. The shape error \c E_GEOM_ALPHA
        \f$ = \sum_c |\alpha_c - \alpha^0_c| |\Omega_c| \f$
        (advectionErrors.H:65) can therefore INCREASE while the volume error
        goes to zero.
      - The documented remedy is connected-component labelling of
        \f$ \{\psi < 0\} \f$ with a per-body \f$ \varepsilon_i \f$ solved from
        a per-body target. It is NOT IMPLEMENTED HERE. Do not read a
        multi-body result from this class as a volume-conserving one.

    WHAT THIS CLASS IS NOT FOR
    It must NOT be used to "fix" the volume drain observed on the COUPLED
    Navier-Stokes droplet. That failure is a mode-4 parasitic current
    amplified 260x that physically squares and destroys the droplet; the
    volume loss is its SYMPTOM, not a bookkeeping error in the level set.
    Applying this correction there restores the volume NUMBER and leaves a
    square droplet of the correct area -- strictly worse than the visible
    failure, because a metric table would then show a converged volume column
    for a run whose interface has been destroyed by an instability. The
    legitimate use of this class is the PRESCRIBED-VELOCITY track (2D vortex,
    3D shear/deformation), where the velocity field is exact and the volume
    error is genuinely the transport scheme's bookkeeping error.

    REPORTING RULE
    Any study using this model must report the shape error
    (\c E_GEOM_ALPHA / \c E_GEOM_ALPHA_REL) ALONGSIDE the volume error, and
    the applied shift lastShift(). The correction trades shape for volume BY
    CONSTRUCTION, so a volume column alone cannot distinguish a scheme that
    conserves volume from a scheme whose volume was conserved for it.

    Usage
    \verbatim
        levelSet
        {
            volumeCorrection
            {
                type                newtonShift;

                maxIter             10;     // Newton iterations (default 10)
                relTol              1e-10;  // on |f|/V_target (default 1e-10)
                maxShiftCells       0.5;    // |eps| <= 0.5 h (default 0.5)
                derivativeStepCells 1e-3;   // delta = 1e-3 h (default 1e-3)

                // targetVolume    1.2e-3;  // optional, base-class entry:
                                            // overrides the captured target
            }
        }
    \endverbatim

See also
    Foam::volumeCorrection (the \c noVolumeCorrection default no-op)
    Foam::phaseIndicator (supplies \f$ \alpha(\psi) \f$, hence \f$ V \f$)

SourceFiles
    newtonShiftVolumeCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef newtonShiftVolumeCorrection_H
#define newtonShiftVolumeCorrection_H

#include "volumeCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class newtonShiftVolumeCorrection Declaration
\*---------------------------------------------------------------------------*/

class newtonShiftVolumeCorrection
:
    public volumeCorrection
{
private:

    // Private Data Members

        //- Maximum number of safeguarded Newton iterations. Each iteration
        //  costs THREE evaluations of V (two for the central-difference
        //  derivative, one for the new residual), and each evaluation is a
        //  full phaseIndicator pass over the narrow band plus one gSum. The
        //  default 10 therefore bounds the per-time-step cost at 31 indicator
        //  passes; quadratic convergence from eps = 0 reaches relTol in 2-4
        //  in practice, and the cap exists to stop a non-smooth
        //  V(eps) -- see maxShiftCells_ -- from spinning.
        label maxIter_;

        //- Convergence test on the RELATIVE volume residual
        //  |V(eps) - V_target| / V_target.
        //
        //  Default 1e-10, not the 1e-12 of the design note. Two floors sit
        //  above 1e-12 and neither is removable: gSum over N cells of doubles
        //  has a relative round-off floor of order sqrt(N)*2.2e-16, and the
        //  geometric indicators obtain alpha from a polyhedron-plane
        //  intersection (geometricPhaseIndicator.C:145-149) whose volume
        //  carries the clipping tolerance of that construction, which is well
        //  above machine epsilon. Requesting 1e-12 on the geometric family
        //  therefore exhausts maxIter_ and emits the non-convergence warning
        //  EVERY step, which trains the reader to ignore it. It is a
        //  dictionary entry: with the heaviside indicator, where V(eps) is an
        //  analytic C1 function of eps, 1e-12 is reachable and appropriate.
        scalar relTol_;

        //- Bracket half-width in cells: |eps| <= maxShiftCells_ * h.
        //
        //  A WELL-POSEDNESS bound, not a tuning knob. Two things break when
        //  the shift approaches h:
        //    (i) both geometric indicators classify NON-band cells by the
        //        SIGN of psi alone (geometricPhaseIndicator.C:101-107,
        //        detrixheAslamPhaseIndicator.C:152-159), so a bulk cell whose
        //        centre crosses zero flips alpha by a full 1 and V(eps) by a
        //        full |Omega_c|. V is then a staircase and the Newton
        //        derivative is meaningless across the step.
        //   (ii) the narrow band the indicator reads is a member reference
        //        bound at construction (geometricPhaseIndicator.C:58) and was
        //        built from the UNSHIFTED psi. Once the shifted zero crossing
        //        leaves that band the cut-cell reconstruction is no longer
        //        applied where the interface actually is.
        //  Both are sub-h effects at 0.5 h and dominate at h. A step needing
        //  more than half a cell of uniform displacement is not a bookkeeping
        //  error and must not be silently absorbed.
        scalar maxShiftCells_;

        //- Central-difference step in cells: delta = derivStepCells_ * h,
        //  f'(eps) ~ (V(eps+delta) - V(eps-delta)) / (2 delta).
        //
        //  Central rather than analytic BY DESIGN. The analytic derivative
        //  exists and is known for each indicator family:
        //    - geometric / detrixheAslam: the least-squares plane fit is
        //      EXACTLY translation-invariant in psi (the normal equations
        //      A p = b of levelSetPlaneReconstruction.C:92-155 have
        //      b -> b + eps A e_4 under psi -> psi + eps, hence
        //      p -> p + eps e_4: the unit normal n is unchanged and the
        //      offset d -> d + eps). V(eps) is thus the cut volume of the
        //      SAME planes translated by eps, and
        //      dV/deps = -sum_c A_c = -|Gamma_h|, minus the discrete
        //      interface area.
        //    - heaviside: alpha = 1 - H_e(psi) gives
        //      dV/deps = -sum_c delta_e(psi_c + eps)|Omega_c|, the smoothed
        //      Dirac of the reference, with the width e_h = (nCells/2) h
        //      already fixed by the indicator (heavisidePhaseIndicator.C:59).
        //  Coding either analytically would require a new virtual on
        //  phaseIndicator and would silently differentiate the WRONG function
        //  for any indicator added later. The central difference reproduces
        //  both to O(delta^2) through the same code path that evaluates the
        //  residual, so the derivative can never disagree with the function
        //  being zeroed.
        //
        //  Default 1e-3: delta must be small against the scale on which
        //  V(eps) curves (h) and large against the noise floor of the cut
        //  volume. With V = O(1) and a double-precision reduction the
        //  round-off contribution to f' is O(V*2.2e-16/delta), which at
        //  delta = 1e-3 h and h = 1/64 is O(1e-11) against a derivative of
        //  order |Gamma_h| = O(1).
        scalar derivStepCells_;

        //- Set on the first correct() call, after the psi boundary types
        //  have been inspected. A patch that FIXES psi breaks the
        //  distance-preservation argument (see correct()), and that must be
        //  said once per run -- not once per time step, which would bury it.
        bool psiBcChecked_;


    // Private Member Functions

        //- V(eps) = sum_c alpha_c(psi + eps)|Omega_c| through the SELECTED
        //  phase indicator. gSum: the returned value, and therefore every
        //  Newton iterate and every convergence test, is identical on every
        //  rank. That is mandatory -- phaseIndicator::calcPhaseIndicator
        //  performs MATCHED MPI (coupledFaceNeighbours,
        //  levelSetPlaneReconstruction.C:34-56), so a rank-divergent
        //  iteration count does not give a wrong answer, it deadlocks.
        scalar shiftedVolume
        (
            const volScalarField& psi,
            volScalarField& alphaTmp,
            volScalarField& psiShift,
            phaseIndicator& phaseInd,
            const scalar eps
        ) const;


public:

    // Static Data Members

        TypeName("newtonShift");


    // Constructors

        //- Construct from mesh (reads levelSet.volumeCorrection)
        newtonShiftVolumeCorrection(const fvMesh& mesh);


    //- Destructor
    virtual ~newtonShiftVolumeCorrection() = default;


    // Member Functions

        //- Solve V(eps) = V_target for the single scalar eps and apply
        //  psi += eps. Returns true if psi was modified.
        virtual bool correct
        (
            volScalarField& psi,
            phaseIndicator& phaseInd
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
