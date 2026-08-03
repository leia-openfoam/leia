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
    Foam::normalProjectedScheme

Description
    Normal-projected semi-Lagrangian level-set scheme (nSL) with per-step
    geometric renormalization of the distance property
    (docs/normal-projected-semi-lagrangian/normal-projected-semi-lagrangian.md).

    Replaces the baseline value transport psi^{n+1}(x_c) = psi^n(x_d) by the
    purely geometric update

        psi^{n+1}(x_c) = d_c + delta

    where, per cell and step, with g, H the gradient and Hessian of the local
    quadratic value fit of psi^n and n = g/|g|:

      - the departure foot is traced along the NORMAL-PROJECTED velocity
        u_n = (u.n) n. Because the tangential velocity satisfies
        u_t . grad(psi) = 0 identically, psi is exactly constant along
        trajectories of u_n for ANY field, drifted or not: nothing is lost by
        discarding u_t (design note Sec. 2);
      - d_c is the signed geometric offset from x_c to the zero set of the
        cell's own fit along the normal ray -- the stable quadratic root with
        a first-order fallback (slReconstruction::signedOffset; note Sec. 5).
        Exactly zero-preserving: the renormalization never moves the front;
      - delta = (x_d - x_c) . n is the normal displacement of the trace
        (projection definition; the arc-length variant of note Sec. 9.4 is
        deferred).

    No value of psi is ever interpolated at x_d: the reconstruction supplies
    directions, speeds and geometry only. The reconstruct-and-evaluate
    resampling of psi -- the operator with measured grid-scale (lambda ~ 2h)
    amplification > 1 in the baseline -- leaves the value path entirely, and
    |grad psi| is driven toward 1 in every active cell at every step without
    any global redistancing operator (which was measured to rebuild distance
    to advection noise).

    TIME LEVELS -- Unew is deliberately UNUSED. The trace is the
    departure-centred second-order form

        x_d = x_c - (dt/2) (3 u_n^n - u_n^{n-1})
                  + (dt^2/2) (u_n^n . grad) u_n^n

    built from the time-n velocity (the Uold argument at all three advect()
    call sites) and the STORED projected velocity of the previous step
    (uProjOld.nSL, projected with the PREVIOUS step's normals -- re-projecting
    the old velocity with the new normals mixes time levels of the geometry
    and is forbidden by the design note Sec. 4.2). This is the only centring
    whose normals exist at their own time level at every call site: an
    arrival-centred trace would need n^{n+1}, which does not exist before the
    advection. The gradient of u_n for the convective term is assembled
    ANALYTICALLY from grad(u) and the fit (Sec. 4.1), never by differencing
    the projected field.

    Cells whose fit is degenerate (constant fallback, or fitted |g| below
    minGradPsi) are FROZEN: psi^{n+1} = psi_c. This is exact for the clamped
    far-field plateau of the signedDistance profile (|grad psi| = 0 there by
    construction); frozen cells re-activate automatically when the interface
    approaches and their fitted |g| grows.

    The slCorrector is NOT used: there is no reconstruction evaluation at the
    foot to correct, and the stencil-bounds clip would be wrong for a
    renormalized value (which may legitimately leave the stencil range of a
    drifted psi). The foot-radius and non-finite guards are applied inline.

    Requires a reconstruction providing fitDerivatives()/signedOffset()
    (the quadratic value-fit models); selecting it with any other
    reconstruction fails with NotImplemented on the first step.

Usage
    fvSolution:
    \verbatim
    levelSet
    {
        semiLagrangian
        {
            scheme          normalProjected;
            reconstruction  uncachedQuadraticWeightedLeastSquares;
            minGradPsi      0.05;    // freeze threshold on the fitted |g|
            offsetBeta      0.25;    // root-acceptance threshold (Sec. 5.3)
        }
    }
    \endverbatim

SourceFiles
    normalProjectedScheme.C

\*---------------------------------------------------------------------------*/

#ifndef normalProjectedScheme_H
#define normalProjectedScheme_H

#include "slScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class normalProjectedScheme Declaration
\*---------------------------------------------------------------------------*/

class normalProjectedScheme
:
    public slScheme
{
    // Private Data

        //- Freeze threshold on the fitted gradient magnitude: cells with
        //  |g| < minGradPsi_ have no usable normal ray (far-field plateau of
        //  the clamped profile, level-set skeleton/ridge cells) and keep
        //  their value unchanged for the step.
        scalar minGradPsi_;

        //- Band gate in stencil radii: the quadratic-root conversion applies
        //  where |psi_c|/|g| <= bandRadii_*stencilRadius(c), i.e. ONLY where
        //  the local fit actually SAMPLES the zero set it converts to
        //  (default 1.0). Beyond one stencil radius the fit extrapolates to a
        //  zero set it never saw, and the root's psi_c*h_nn term amplifies
        //  fitted-Hessian noise linearly in the offset -- measured: ~30 h
        //  one-step errors at 45 h offset (domain corners), and a closed-loop
        //  runaway seeded at 2.6 h from a circle's cone tip when the gate
        //  admitted offsets up to 3 stencil radii (the tip kink pollutes the
        //  fitted h_nn of cells whose stencils reach it; the written error
        //  re-enters the next fit and compounds at ~1 h per step). Outside the gate
        //  the field is transported RAW, first order along the normal ray
        //  (psi_c + delta |g|): no division by the fitted |g| anywhere, so fit
        //  error enters only through the per-step increment (~1e-2 h), and a
        //  plateau cell (|g| ~ 0) freezes automatically. Only the band is
        //  renormalized -- and only the band feeds the physics. Any far-field
        //  RENORMALIZATION was measured to amplify fitted-|g| error linearly
        //  in the offset (boundary-ring cells with one-sided stencils: ~37 h
        //  one-step errors), so it is deliberately absent.
        scalar bandRadii_;

        //- In-band renormalization mode:
        //    geometric  write d_c + delta from the fit's signed offset (the
        //               design note's Sec. 5 mechanism);
        //    strain     transport the raw value along the ray and integrate
        //               the exact gradient-magnitude ODE instead,
        //               D|grad psi|/Dt = -|grad psi| eps_nn: the written value
        //               is (psi_c + delta |g|) (1 + eps_nn dt) with
        //               eps_nn = n.D.n from grad(u) (note Sec. 7.1). Fit noise
        //               enters this correction only MULTIPLIED BY dt;
        //    none       raw ray transport everywhere (no renormalization).
        //  MEASURED (sigma = 0 rigid translation, N = 64, 100 steps): the
        //  geometric write-back diverges -- band error 1.1 h and growing, for
        //  BOTH the quadratic-root and the first-order conversion -- because
        //  writing any fit-derived offset back into psi feeds neighbour noise
        //  into the next step's fit at order one per step, independent of dt.
        //  Raw transport is stable (0.012 h after the same 100 steps).
        word renormalization_;

        //- Offset engine of the geometric renormalization (Variant (a) vs (b)
        //  of the design note Sec. 7.2):
        //    rayRoot    the closed-form stable quadratic root along the
        //               normal ray (slReconstruction::signedOffset), and
        //               delta as the flat projection (x_d - x_c).n;
        //    footPoint  the stabilized foot-point algorithm on the local
        //               model (slReconstruction::footPointDistance,
        //               stable-foot-point-3d.md): d_c as the foot-point
        //               distance of x_c to {R_c = 0}, and the normal
        //               displacement as the foot-point distance of x_d to the
        //               cell's own iso-surface {R_c = psi_c} -- curvature-
        //               aware on both legs. Per-cell fallback to rayRoot when
        //               either search fails its guards (counted).
        //  Read only when renormalization is geometric.
        word offsetEngine_;

        //- True once uProjOld_ holds a valid previous-step projection --
        //  either read from disk on restart or written by a completed step.
        //  While false, the first step runs with u_n^{n-1} = u_n^n (the
        //  standard one-off Adams-Bashforth startup, within the second-order
        //  budget).
        bool haveHistory_;

        //- The previous step's normal-projected velocity, projected with the
        //  PREVIOUS step's normals (design note Sec. 4.2). Persisted with
        //  AUTO_WRITE/READ_IF_PRESENT so a restarted trace uses the same two
        //  velocity levels as an uninterrupted run (the UextOld precedent).
        volVectorField uProjOld_;


public:

    //- Runtime type information
    TypeName("normalProjected");


    // Constructors

        explicit normalProjectedScheme(const fvMesh& mesh);

    //- Destructor
    virtual ~normalProjectedScheme() = default;


    // Member Functions

        //- Advance psi^n -> psi^{n+1} by the normal-projected geometric
        //  update. Unew is unused (see Description); the corrector is unused.
        virtual void advance
        (
            volScalarField& psi,
            const volVectorField& Unew,
            const volVectorField& Uold,
            slReconstruction& recon,
            slCorrector& corrector
        ) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
