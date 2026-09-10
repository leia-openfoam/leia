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
    Foam::lipschitzConeValueBound

Description
    The distance-cone bound. A signed distance function is L-Lipschitz, so its
    exact value at the departure point x_d lies in the interval carried by EVERY
    stencil point j:

        psi_j - L|x_d - x_j|  <=  psi(x_d)  <=  psi_j + L|x_d - x_j|

    and the intersection over the arrival cell's stencil T_c is the bound:

        l_c = max_{j in T_c} ( psi_j - L|x_d - x_j| )
        u_c = min_{j in T_c} ( psi_j + L|x_d - x_j| )
        psi_c^{n+1} = clip(H_c, l_c, u_c)

    with H_c the quadratic candidate. It needs the stencil positions, the old
    values, the departure point and L. There is no interface projection, no
    closest-point search and no apex detector.

    WHY IT IS NOT THE MONOTONE CLIP. The bound admits a correctly deepening
    extremum, which is the case the stencil-range clip could not admit. One
    verified 1D example, from the external review of 2026-09-09 and reproduced
    independently: psi(x) = |x - 0.25| - 2 sampled at x = -1, 0, 1 gives
    psi = -0.75, -1.75, -1.25. The true value at x_d = 0.25 is -2, BELOW the
    sample minimum -1.75, so the monotone clip commits an error of 0.25. The cone
    interval is [-2, -1.5] and admits the true value at its lower end.

    WHAT IT GUARANTEES. One step, for data that satisfies the L-Lipschitz
    condition. Then the true departure value lies inside [l_c, u_c], and clipping
    onto an interval that contains the true value moves the candidate toward it:

        |clip(H_c, l_c, u_c) - psi_true|  <=  |H_c - psi_true|

    cell by cell, so the bound can only reduce the pointwise error. With the
    owner cell in the stencil it also gives |psi_c^{n+1} - psi_c^n| <= L|delta_c|.

    WHAT IT DOES NOT GUARANTEE, and the review is explicit about each: it does
    NOT make the new field globally L-Lipschitz, it does NOT make the composite
    map monotone or non-expansive, and the linear amplification diagnostic Lambda
    may NOT be set to one after a nonlinear clip. There is no order-of-accuracy
    claim. Repeated-update stability must be MEASURED on the mapping itself.

    WHERE L COMES FROM. Not from tuning.

      * lipschitzMode unity   -- L = lipschitzConstant, default exactly 1, the
                                 value for an exact signed distance field. Right
                                 wherever grad u = 0, which covers the
                                 frozen-velocity probe and the uniform-stream
                                 translating gate. The production candidate.
      * lipschitzMode stencil -- L = max(lipschitzConstant, largest measured
                                 owner-to-neighbour slope). The review warns that
                                 estimating L from the largest noisy slope "can
                                 make the guard ineffective", so this arm is a
                                 DIAGNOSTIC: it measures how loose a data-driven
                                 L becomes, not a fix.

    Any lipschitzConstant other than 1 is a tuned coefficient and voids the
    coefficient-free claim. It exists so the sensitivity can be measured.

    Advection does not preserve L = 1. The review gives the transport of the
    gradient and the resulting bound,

        D_t grad psi = -(grad u)^T grad psi,
        L(t) <= L(0) exp( int_0^t ||grad u(s)||_inf ds ),

    so a flow-map L is the defensible route under strain. It is NOT implemented
    here: it needs a ||grad u||_inf collective every step and a constant carried
    in time across restarts, and a collective inside this call path is the
    deadlock class that has already cost this campaign two runs. It is separate
    work, gated on the unity mode measuring something.

    ADMISSIBILITY. l_c > u_c is possible, and it is the guard's main risk and its
    best diagnostic. An empty interval requires
    psi_i - psi_k > L(r_i + r_k) >= L|x_i - x_k| by the triangle inequality
    through x_d, so it happens exactly when the OLD data already violates the
    Lipschitz condition between two stencil points. The crossing count therefore
    MEASURES how far the level set has drifted from a distance function. A large
    crossing fraction falsifies the guard.

    The review forbids the easy escapes -- "Quietly swapping the bounds, or
    inflating L without limit, removes their mathematical meaning" -- and asks
    for a documented recovery path:

      * onInadmissible cellOnly -- the owner-only interval
                                   [psi_c - L|delta_c|, psi_c + L|delta_c|],
                                   never empty, and it keeps the one bound that
                                   is proved. The default.
      * onInadmissible none     -- no bound in that cell, counted.

    BOUNDARY FACES. With the default `stencilBoundaryFaces include`, every
    physical boundary FACE is a stencil data point at about h/2 carrying the
    cell's own value through zeroGradient. That value is NOT a sample of the
    distance function -- the exact psi has a non-zero wall-normal gradient at the
    Popinet inlet -- so including it gives an interval as tight as
    [psi_c - L h/2, psi_c + L h/2] which can EXCLUDE the true value and clip a
    correct one. `coneBoundaryFaces` therefore defaults to false. It is exposed
    because the same mechanism is already on record as the frozen inlet layer
    (1 to 20 % of the exact transport).

    UNSTRUCTURED AND PARALLEL BY CONSTRUCTION. The bound is stated per cell over
    the arrival cell's own stencil, using positions and values that
    collectStencil already gathered across processor and cyclic patches. No
    collective, no extra halo exchange, no structured index space, and no
    dependence on the interface shape.

SourceFiles
    lipschitzConeValueBound.C

\*---------------------------------------------------------------------------*/

#ifndef lipschitzConeValueBound_H
#define lipschitzConeValueBound_H

#include "slValueBound.H"
#include "Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class lipschitzConeValueBound Declaration
\*---------------------------------------------------------------------------*/

class lipschitzConeValueBound
:
    public slValueBound
{
    // Private data

        //- unity | stencil. See the class description.
        const word lipschitzMode_;

        //- L for the unity mode, and the floor for the stencil mode. Default 1,
        //  the exact signed-distance value. Anything else is a tuned
        //  coefficient.
        const scalar L0_;

        //- cellOnly | none. What to do with an EMPTY interval.
        const word onInadmissible_;

        //- Whether physical boundary-face stencil entries carry the bound.
        //  Default false: a zeroGradient face value is not a sample of the
        //  distance function.
        const Switch useBoundaryFaces_;

        //- true when lipschitzMode_ == "stencil".
        const bool measureL_;

        //- true when onInadmissible_ == "cellOnly".
        const bool recoverCellOnly_;


public:

    //- Runtime type information
    TypeName("lipschitzCone");

    //- Construct from mesh and the semiLagrangian sub-dict.
    lipschitzConeValueBound(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~lipschitzConeValueBound() = default;


    // Member Functions

        virtual void interval
        (
            const slReconstruction& recon,
            const label c,
            const point& foot,
            const bool eligible,
            const scalar stencilLo,
            const scalar stencilHi,
            slBoundResult& r
        ) const override;

        virtual void printBanner() const override;
};


} // End namespace Foam

#endif

// ************************************************************************* //
