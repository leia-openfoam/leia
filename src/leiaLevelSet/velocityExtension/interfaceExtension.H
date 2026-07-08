/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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
    Foam::interfaceExtension

Description
    Intermediate (non-selectable) base implementing the shared backbone for
    interface-normal velocity extension:

      1. build a band (seed = the narrow-band interface layer, dilated nLayers);
      2. seed Uext at the interface layer with the true interface velocity:
         foot point xSigma = C - psi*n (one Newton step toward the interface),
         then cell-point (IDW) interpolation of U at xSigma;
      3. extend() -- pure virtual -- propagate Uext through the band so that
         (n.grad)Uext = 0, with the seed layer pinned as Dirichlet data;
      4. build the advection flux phiExt from Uext in the band.

    Concrete models implement only extend().

SourceFiles
    interfaceExtension.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceExtension_H
#define interfaceExtension_H

#include "velocityExtension.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interfaceExtension Declaration
\*---------------------------------------------------------------------------*/

class interfaceExtension
:
    public velocityExtension
{
protected:

    //- Level set and interface volume fraction (configurable registry keys).
    const volScalarField& psi_;
    const volScalarField& alpha_;

    //- Number of cell layers the interface layer is dilated by to form the
    //  application band (where Uext is used for the advecting flux).
    label nLayers_;

    //- Number of gradient-descent steps from cell centre to the interface.
    label nDescent_;

    //- Number of additional cell layers around the interface layer whose
    //  velocity is also foot-pointed (steepest descent + interpolation) and
    //  pinned as Dirichlet data in the extension solve. The PDE stencil of an
    //  interface-adjacent cell otherwise mixes solved and pinned values,
    //  producing kinks in Uext. 0 = pin the interface layer only (classic).
    //  Clamped to nLayers_ - 1 (anchors must stay inside the band).
    label nAnchorLayers_;

    //- Project phiExt to divergence-free after extension (sacrifices exact
    //  normal-constancy for volume conservation).
    bool projectFlux_;

    //- Safety cap on the extended velocity: |Uext| <= maxScale_ * max|U| in the
    //  band (guards against runaway feedback). Disabled if <= 0.
    scalar maxScale_;

    //- Extension band marker (AUTO_WRITE for visualization).
    volScalarField band_;

    //- Interface unit normal n = grad(psi)/|grad psi| (recomputed each step).
    volVectorField nHat_;

    //- Interface-layer (seed) cells where Uext is Dirichlet data.
    labelList seedCells_;

    //- Per-cell characteristic size: minimum face-neighbour centre distance
    //  (1/deltaCoeffs). Correct for 2D meshes of arbitrary thickness, where
    //  cbrt(V) = (h^2 t)^(1/3) != h mis-scales every h-proportional length
    //  (band/fade widths, tolerances, regularizations).
    scalarField cellSize_;


    // Backbone steps

        //- nHat_ = grad(psi)/|grad psi|.
        void computeNormals();

        //- Fill cellSize_ (min face-neighbour distance; once, mesh static).
        void computeCellSize();

        //- Fill band_ and seedCells_ (seed = narrowBand, dilated nLayers).
        void buildBand();

        //- Seed Uext_ at seedCells_ with the interpolated interface velocity.
        void buildSeedVelocity();

        //- Clip |Uext_| <= maxScale_*max|U_| in the band (runaway safety net).
        void clipExtended();

        //- Reset Uext_ = U_ outside the band; build phiExt_ from Uext_.
        void updateFlux();

        //- Model-specific extension of Uext_ through the band (seed pinned).
        virtual void extend() = 0;


public:

    // Constructors

        interfaceExtension(const fvMesh& mesh);

    //- Destructor
    virtual ~interfaceExtension() = default;

    // Member Functions

        //- Backbone: seed + extend + flux.
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
