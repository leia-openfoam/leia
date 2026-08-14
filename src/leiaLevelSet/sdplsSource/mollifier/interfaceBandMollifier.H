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
    Foam::interfaceBandMollifier

Description
    SDPLS cut-off built from the interface TOPOLOGY rather than from \f$|\psi|\f$.

    The cut-off of \cite fricke_locally_2022 eq. (24), implemented in
    Foam::mollifier1, thresholds \f$ |\psi| \f$ against physical lengths
    \f$ w_1, w_2 \f$. That is exact only while \f$ \psi \f$ is a signed distance,
    which is precisely the property SDPLS exists to maintain and does not fully
    achieve. Measured here against the unmollified control at \f$ N=128 \f$,
    fitted over the shared rungs 64/81/102/128:

    \verbatim
                      far field      band L2        volume        shape
      2D vortex       12x better   26% worse      equal         23% worse
      3D shear       8.6x better  3.4x worse   +0.23 -> -0.97  +0.92 -> -0.58
      3D deformation 1.5x better  7.8x worse   +2.36 -> +1.07  +0.99 -> +0.58
                                  order -0.54
    \endverbatim

    The failure tracks the deformation: on the LeVeque case the interface is
    drawn into sheets thinner than a cell, \f$ |\nabla\psi| \f$ departs far from
    unity, and a fixed \f$ |\psi| \f$ threshold stops marking a fixed distance --
    it clips the source out of the region where it is still needed, and the band
    error then DIVERGES.

    This cut-off never asks how far a cell is from the interface. It asks
    whether the phase indicator changes across a face:

      1. Seed. Mark every cell adjacent to a face with
         \f$ |\mathtt{snGrad}(\alpha)| > 0 \f$. \f$ \alpha \f$ is 0 or 1 in the
         bulk, so this is a purely topological marker of the interface, immune
         to whatever \f$ |\nabla\psi| \f$ is doing.
      2. Dilate \c nLayers times over cell neighbours, recording the layer index
         at which each cell was first reached. ACROSS PROCESSOR BOUNDARIES --
         see below.
      3. Taper. \f$ G = 1 \f$ on the seed layer and falls smoothly to 0 just
         outside the last layer.

    \b Parallel \b correctness \b is \b not \b optional \b here.
    Foam::neighbourNarrowBand dilates with \c mesh().cellCells(), which is local
    topology: a cell one layer across a processor boundary from a marked cell is
    never reached. The band is then torn at every processor boundary and the
    result depends on the DECOMPOSITION -- the same case on 4 and on 16 ranks
    would cut the source off in different places. This class exchanges the
    marker over coupled patches after every dilation pass, as
    Foam::signChangeNarrowBand already does for its seed.

    The taper keeps the two properties that make a cut-off safe for SDPLS, for
    the same reasons as eq. (24):
      - it is symmetric about the interface, because the seed straddles it and
        dilation is isotropic;
      - it is FLAT (identically unity) on the seed layer, so it contributes no
        variation of the effective source coefficient across the two cells
        straddling \f$ \psi = 0 \f$. That variation is the only thing that
        displaces the discrete zero crossing, so a taper with a kink at the
        interface would introduce the interface motion SDPLS is built to avoid.

    Selected at run time as \c band:
    \verbatim
        sdplsSource
        {
            type    R;
            mollifier
            {
                type      band;
                nLayers   3;      // dilation passes beyond the seed layer
            }
        }
    \endverbatim

See also
    Foam::mollifier1 (the |psi|-threshold cut-off of eq. 24)

SourceFiles
    interfaceBandMollifier.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceBandMollifier_H
#define interfaceBandMollifier_H

#include "mollifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class interfaceBandMollifier Declaration
\*---------------------------------------------------------------------------*/

class interfaceBandMollifier
:
    public mollifier
{
private:

    // Private Data Members

        //- Dilation passes beyond the seed layer. The cut-off is unity on the
        //  seed and reaches zero one layer past the last dilated layer, so the
        //  source acts over roughly nLayers cells either side of the interface.
        //  A CELL COUNT, deliberately: unlike a physical width it needs no
        //  assumption about |grad psi|, and it is what "a couple of cells
        //  either side" actually means.
        label nLayers_;

        //- Name of the phase indicator field to read the topology from.
        word alphaName_;

public:

    // Static Data Members

        TypeName("band");

    // Constructors

        //- Construct from fvSolution
        interfaceBandMollifier(const dictionary& dict);

    //- Destructor
    virtual ~interfaceBandMollifier() = default;

    // Member functions

        //- Cut-off field from the interface topology. Overrides the pointwise
        //  base implementation: this cut-off is not a function of psi[cellI].
        virtual tmp<volScalarField> field(const volScalarField& psi) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
