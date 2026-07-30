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
    Foam::isoCurvatureSurfaceTension

Description
    ISO-VALUE-AGNOSTIC geometric surface-tension model (2D prototype).

    The interface is read from the ISO-VALUE crossings of ANY field on the mesh,
    so the model is independent of how the interface is encoded:
      - level set psi at iso 0,
      - volume fraction alpha at iso 0.5,
      - an algebraic implicit surface (e.g. ellipsoid) at iso 0.
    It only needs sign(field - isoValue) to change across a face.

    Curvature-normal (geometric conormal). For every internal face the interface
    crosses (sign(field-isoValue) differs between owner and neighbour), the
    interface-in-face segment vector L_f (in 2D the empty-direction segment
    h_z e_z) is crossed with the interface normal n_f at the crossing, where n_f
    is the OWNER<->NEIGHBOUR MEAN of nHat = grad(field)/|grad(field)| (a robust,
    iso-magnitude-independent unit normal, epsilon-guarded). Summed over a cell's
    crossed faces this is the discrete conormal line integral

        (kappa n)_c A_c = sum_f ( L_f x n_f ),

    the change of the interface tangent around the cut (the curvature-normal),
    which does NOT divide by |grad field| and is insensitive to whether the field
    is a signed distance. The scalar mean curvature is the projection
    kappa_c = ((kappa n)_c . nHat_c) / A_c, with A_c = |segment|_c h_z the
    interface area in the cell. kappa is written into the registered "kappa" field.

    The force is then the standard CSF face force (as reconstructedCurvature):

        f_{Sigma,f} = sigma_f interpolate(kappa)_f snGrad(alpha)_f    [N/m^3]

    Selected with `surfaceTensionForce { type isoCurvature; }`. Dict keys:
        field <name>    the field whose iso-value defines the interface (def "psi")
        isoValue <s>    the iso-value (def 0; use 0.5 for a VOF alpha)
        alpha <name>    phase indicator for the CSF weight snGrad(alpha)
        kappa <name>    curvature field to fill (def "kappa")
    Run with curvatureExtension { type fvm; } so the solver does not overwrite
    the kappa this model fills. 2D (one empty direction) prototype.

SourceFiles
    isoCurvatureSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef isoCurvatureSurfaceTension_H
#define isoCurvatureSurfaceTension_H

#include "surfaceTensionForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class isoCurvatureSurfaceTension Declaration
\*---------------------------------------------------------------------------*/

class isoCurvatureSurfaceTension
:
    public surfaceTensionForce
{
    protected:

        const fvSolution& fvSolutionDict_;
        const dictionary& levelSetDict_;
        const dictionary& surfTensionDict_;

        //- Phase indicator for the CSF force weight snGrad(alpha).
        const volScalarField& alpha_;

        //- Field whose iso-value defines the interface (psi, alpha, ...).
        const volScalarField& field_;

        //- Iso-value of field_ that IS the interface (0 for psi, 0.5 for alpha).
        const scalar isoValue_;

        //- Name of the curvature field to fill (looked up non-const each step).
        const word kappaName_;

        //- Overall curvature-sign factor (dict "sign", default +1; set -1 to flip).
        //  The geometric construction fixes a CONSISTENT sign around the interface;
        //  this one factor aligns it with the +1/R-for-a-drop convention.
        const scalar sgn_;

        //- Number of layers to EXTEND the interface-cell curvature outward into the
        //  neighbouring cells (dict "nExtend", default 2). Without it kappa is nonzero
        //  only in the one-cell interface band, so fvc::interpolate(kappa) HALVES it at
        //  the band-edge faces (where snGrad(alpha) is still nonzero) -> unbalanced CSF
        //  force -> growing parasitic currents. Extension fills those neighbour cells so
        //  the face-interpolated kappa is correct across the whole snGrad(alpha) support.
        const label nExtend_;

        //- Number of interface-restricted neighbour-averaging (Brackbill) smoothing
        //  passes applied to the geometric curvature (dict "nSmooth", default 0).
        const label nSmooth_;

        //- Time under-relaxation of kappa (dict "relax", default 1 = none): blend the
        //  new geometric kappa with the previous step's, kappa = relax*new +
        //  (1-relax)*old, to damp the per-step curvature jitter as the interface-cell
        //  set changes sub-cell (the driver of the localised force spikes).
        const scalar relax_;

        //- Curvature reduction (dict "method"): "conormal" (default, the integral
        //  identity (kappa n)_c A_c = sum_f L_f x n_f from the iso-points on the face
        //  edges) or "parabola" (local tangent-frame height-function fit). Stored as
        //  0 = conormal, 1 = parabola for computeIsoGeometricCurvature.
        const label method_;

public:

    // Static Data Members
    TypeName("isoCurvature");

    // Constructors
    isoCurvatureSurfaceTension(const fvMesh& mesh);

    //- Destructor
    virtual ~isoCurvatureSurfaceTension() = default;

    // Member Functions
    virtual tmp<surfaceScalarField> calcFaceSurfaceTensionForceFlux() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
