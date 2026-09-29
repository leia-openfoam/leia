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
    Foam::slGradientControlSource

Description
    The gradient-control source of the semi-Lagrangian line
    (docs/plan-halo-limited-gradient-control.md, D3). After the
    characteristic step it integrates \f$ D_t \psi = \psi F \f$ exactly for
    one step with F frozen:
    \f{equation}{
        \psi_c \leftarrow \psi_c \exp(\Delta t\, F_c),
        \qquad
        F = w(q)\, a + G(q, \sigma),
    \f}
    so no value changes sign and the zero set of psi does not move. The law
    and its strain weight are a Foam::gradientControlLaw
    (libleiaGradientControl) from the \c law sub-dictionary.

    The inputs are cell values of the arrival field psi^A:
    - \f$ q = |\mathbf{g}| \f$, \f$ \mathbf{n} = \mathbf{g}/q \f$ with
      \f$ \mathbf{g} \f$ = fvc::grad(psi, "gradPsiSource"): the case selects the
      gradient scheme, and every reconstruction works with it (the fit
      derivatives exist only in the uncached quadratic fit);
    - \f$ \mathbf{D} = \mathrm{symm}(\nabla\mathbf{u}) \f$ with
      fvc::grad(Utrajectory, "gradUSource") of the velocity that moved psi,
      only if the law reads \f$ \sigma = |\mathbf{D}|_F \f$ or the strain
      weight reads \f$ a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n} \f$.

    Band: F acts only where \f$ |\psi_c|/q_c \le \f$ bandCells \f$ h_c \f$,
    with \f$ h_c \f$ the smallest distance from the centre of cell c to the
    centre of a face neighbour (processor neighbours included), so the band
    does not depend on the decomposition. F = 0 exactly outside the band and
    where q = 0.

    |dt F| is clamped at 30. Every clamp is counted, and the field
    \c slSourceClamp marks the clamped cells of the last step. The field
    \c slSourceF holds the applied F. Both are written with the other fields.
    Per step the log line reports the band cells, the clamps and max |dt F|
    (collective reductions on every rank).

    \verbatim
        semiLagrangian
        {
            source
            {
                type        gradientControl;
                bandCells   3;
                law
                {
                    type    linearZ;
                    mu      10;
                    strainWeight { type none; }
                }
            }
        }
    \endverbatim

SourceFiles
    slGradientControlSource.C

\*---------------------------------------------------------------------------*/

#ifndef slGradientControlSource_H
#define slGradientControlSource_H

#include "slSource.H"
#include "gradientControlLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class slGradientControlSource Declaration
\*---------------------------------------------------------------------------*/

class slGradientControlSource
:
    public slSource
{
    // Private Data

        //- The band half width in cell sizes
        const scalar bandCells_;

        //- The law G and its strain weight w
        autoPtr<gradientControlLaw> law_;

        //- The cell size h_c: the smallest centre-to-centre distance to a
        //  face neighbour, processor neighbours included
        scalarField h_;

        //- The applied F [1/s] of the last step (0 outside the band)
        volScalarField F_;

        //- 1 in the cells where |dt F| was clamped in the last step
        volScalarField clamp_;

        //- The number of clamps since the start of the run (all ranks)
        label nClampTotal_;


    // Private Member Functions

        //- h_c from the face neighbours (coupled patches included)
        void computeCellSize();


public:

    //- The clamp of |dt F|
    static constexpr scalar maxExponent = 30;

    //- Runtime type information
    TypeName("gradientControl");


    // Constructors

        explicit slGradientControlSource(const fvMesh& mesh);


    //- Destructor
    virtual ~slGradientControlSource() = default;


    // Member Functions

        virtual bool active() const
        {
            return true;
        }

        //- psi <- psi exp(dt F) in the band; see the class description
        virtual void apply
        (
            volScalarField& psi,
            const volVectorField& Utrajectory,
            slReconstruction& geometryFit,
            const scalar dt
        );

        //- The law
        const gradientControlLaw& law() const
        {
            return *law_;
        }

        //- The cell size of the band criterion
        const scalarField& cellSize() const
        {
            return h_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
