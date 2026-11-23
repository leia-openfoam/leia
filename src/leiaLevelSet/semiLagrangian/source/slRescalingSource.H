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
    Foam::slRescalingSource

Description
    The algebraic rescaling of the level set (D. Bothe, notes of 2026-10) on
    the semi-Lagrangian line (work package 3 of the plan "SLLS on cubic
    meshes with dynamic AMR"). The map multiplies the cell value by a
    positive factor, so no cell value changes its sign:

      order first   psi_c <- psi_c/|g|
      order second  psi_c <- 2 psi_c/(|g| + sqrt(max(|g|^2 - 2 psi_c d_nn, 0)))
                    d_nn = n.H.n, n = g/|g|

    g and H are the gradient and the Hessian of the semi-Lagrangian fit of
    the arrival field at the cell centre (slReconstruction::fitDerivatives).
    The model does not differentiate psi with a finite-volume operator. The
    formula is rescaledDistance::estimate (libleiaCore), which the Eulerian
    rescalingRedistancer also calls.

    The reconstruction (the geometry fit of slAdvection) must implement
    fitDerivatives: quadraticWeightedLeastSquares or
    uncachedQuadraticWeightedLeastSquares. Every other reconstruction stops
    with NotImplemented in the base class.

    Per cell:
    - fit order 2: the selected order;
    - fit order 1 (the linear fallback, H = 0): the first-order map;
    - fit order 0, or |g| < gradFloor: the cell keeps its value and the
      model counts it as skipped.

    Region: the registered field NarrowBand of the narrowBand model
    (libleiaCore), dilated by bandLayers layers of face neighbours. Each
    layer is exchanged across the coupled patches before the next layer, so
    the region does not depend on the decomposition. The solver computes
    NarrowBand after each advection step. The seed is therefore the band of
    psi at the start of the step, not the band of the arrival field. With
    CFL <= 1 the sign change moves by at most one cell layer per step, and
    the dilation covers that layer. bandLayers is a count of cell layers
    with one default for every case: the compact curvature reads psi in the
    layers l <= 3 around the interface cells (CLAUDE.md "No free
    parameters", item 1).

    Trigger:
    - gradientDeviation (default): the map runs in a step only if
      max | |g| - 1 | > theta over the cells of the region with fit order
      >= 1. theta is one dimensionless number.
    - everyStep: the map runs in every step. A research instrument for the
      control arm. On the Eulerian line the every-step map failed
      (docs/knowledge-base/concepts/algebraic-rescaling-eulerian.md).
    maxEvents >= 0 limits the number of events (-1: no limit).

    All new values come from the arrival field (pass 1). Then the model
    assigns them (pass 2) and calls psi.correctBoundaryConditions().

    Every event prints one line: the time index, the measured maximum, the
    cells rescaled, the cells skipped and the cells with a negative
    discriminant. The destructor prints the totals. The reductions are
    collective and run on every rank.

    Not yet measured (STATUS 11.44).

    The keys are in levelSet.semiLagrangian.source, next to `type`, as for
    the model gradientControl:

    \verbatim
        semiLagrangian
        {
            source
            {
                type        rescaling;
                order       second;             // first | second
                trigger     gradientDeviation;  // gradientDeviation | everyStep
                theta       0.1;
                maxEvents   -1;                 // optional
                bandLayers  3;                  // optional
                gradFloor   1e-3;               // optional
            }
        }
    \endverbatim

SourceFiles
    slRescalingSource.C

\*---------------------------------------------------------------------------*/

#ifndef slRescalingSource_H
#define slRescalingSource_H

#include "slSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class slRescalingSource Declaration
\*---------------------------------------------------------------------------*/

class slRescalingSource
:
    public slSource
{
    // Private Data

        //- 1: psi/|g|; 2: the quadratic-root distance
        label order_;

        //- gradientDeviation | everyStep
        const word trigger_;

        //- The threshold of max | |g| - 1 | (trigger gradientDeviation)
        const scalar theta_;

        //- Largest number of events; -1 = no limit
        const label maxEvents_;

        //- Layers of face neighbours added to the registered NarrowBand
        const label bandLayers_;

        //- |g| below this value: the cell keeps its value
        const scalar gradFloor_;

        //- Totals since the start of the run (all ranks)
        label nEvents_;
        label nCalls_;
        label nRescaledTotal_;
        label nSkippedTotal_;


    // Private Member Functions

        //- The region: NarrowBand dilated by bandLayers_ face layers
        void regionMask(boolList& mask) const;


public:

    //- Runtime type information
    TypeName("rescaling");


    // Constructors

        explicit slRescalingSource(const fvMesh& mesh);


    //- Destructor: prints the totals
    virtual ~slRescalingSource();


    // Member Functions

        virtual bool active() const
        {
            return true;
        }

        //- The rescaling of psi in the region, if the trigger fires
        virtual void apply
        (
            volScalarField& psi,
            const volVectorField& Utrajectory,
            slReconstruction& geometryFit,
            const scalar dt
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
