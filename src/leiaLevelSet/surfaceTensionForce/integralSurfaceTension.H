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
    Foam::integralSurfaceTension

Description
    Conservative, well-balanced INTEGRAL surface-tension model (CST) of
    Abu-Al-Saud, Popinet & Tchelepi (2018, J. Comput. Phys. 371:896-913),
    transcribed for the face-flux (interFoam-family) momentum balance.

    Surface tension enters as the divergence of a stress,
    \f[ -\nabla\cdot(p\,\mathbf I + \boldsymbol\sigma), \f]
    with the stress components built from the INTERFACE TANGENTS at the
    intersections of the psi=0 contour with the momentum-control-volume faces.
    The paper-specific transcription also contains gamma*kappa*s pressure-jump
    corrections at those intersections. Their curvature source can be selected
    as field, constant, or none for mechanism isolation; removing that
    correction is an ablation, not a complete curvature-free integral method.
    The tangent part itself uses no curvature field or interface thickness and
    bypasses the noise-amplifying second derivative. The reference
    implementation recovers Laplace equilibrium to machine precision on a
    reinitialisation-free level set (their Fig. 4).

    PROTOTYPE SCOPE: 2D, uniform, single-block hexahedral meshes (one empty
    direction), constant surface tension, serial. The structured (i,j)
    addressing is derived from the cell centres and guarded; unstructured and
    Marangoni generalisations are documented follow-ups.

SourceFiles
    integralSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef integralSurfaceTension_H
#define integralSurfaceTension_H

#include "surfaceTensionForce.H"
#include "fvSolution.H"
#include "volFieldsFwd.H"
#include "labelList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class integralSurfaceTension Declaration
\*---------------------------------------------------------------------------*/

class integralSurfaceTension
:
    public surfaceTensionForce
{
    protected:

        const fvSolution& fvSolutionDict_;
        const dictionary& levelSetDict_;
        const dictionary& surfTensionDict_;

        //- Level set (signed distance; psi < 0 inside phase 1).
        const volScalarField& psi_;

        //- Cell-centred mean curvature (filled by the solver; used only for
        //  the pressure-jump corrections at the interface crossings).
        const volScalarField& kappa_;

        //- Curvature used by the paper's pressure-jump correction: "field"
        //  linearly interpolates kappa_ to each crossing, "constant" uses the
        //  dictionary value, and "none" removes this correction so the force
        //  is assembled solely from interface-crossing tangents. The remaining
        //  integral-force machinery is shared.
        const word curvatureSource_;
        const scalar constantCurvature_;

        //- ISO-VALUE-AGNOSTIC normalisation (dict "normalizeGradient", default TRUE):
        //  divide the psi-derivative interface tangents by the local |grad psi| --
        //  SMOOTHED (3 Jacobi passes; raw central differences carry grid-scale noise
        //  that would re-enter the stress as a vortical driver) and interpolated
        //  consistently TO EACH CROSSING. With it the scheme accepts ANY smooth
        //  monotonic level set (epsilon-stabilised, like every normal/curvature model
        //  in the literature): measured stable with decaying currents on both the SDF
        //  circle and the strongly non-SDF algebraic ellipsoid (|grad psi| ~ 1000),
        //  no redistancing. false reproduces the historical SDF-only formulation.
        const bool normalizeGrad_;

        // Structured 2D addressing (built once; uniform-hex guarded)
        mutable bool addressingBuilt_;
        mutable label nx_, ny_;
        mutable scalar dx_;
        mutable labelList cellI_, cellJ_;     // cell -> (i,j)
        mutable labelList ijCell_;            // i + j*nx -> cell
        mutable label dirX_, dirY_;           // global components used as x,y

        void buildAddressing() const;

        //- cell label at (i,j), or -1 outside the block.
        label at(const label i, const label j) const
        {
            if (i < 0 || j < 0 || i >= nx_ || j >= ny_) { return -1; }
            return ijCell_[i + j*nx_];
        }

public:

    // Static Data Members
    TypeName ("integralSurfaceTension");

    // Constructors
    integralSurfaceTension(const fvMesh& mesh);

    //- Destructor
    virtual ~integralSurfaceTension() = default;

    // Member Functions
    virtual tmp<surfaceScalarField> calcFaceSurfaceTensionForceFlux() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
