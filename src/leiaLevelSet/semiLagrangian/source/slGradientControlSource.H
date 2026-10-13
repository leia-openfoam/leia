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

    qScheme (optional, default centred) selects how q is computed:
    - centred: \f$ q = |\mathbf{g}| \f$ as above. On a uniform hexahedral
      mesh this is the centred difference. It reverses the group velocity of
      the sawtooth mode: the source equation carries smooth errors away from
      the interface, but the centred q carries the sawtooth envelope onto it
      (zero-flow gate E1, STATUS 11.26). Kept as the inert default.
    - upwindFit: two passes over face neighbours. Pass 1 updates the
      semi-Lagrangian fit on the current psi and reads its gradient
      \f$ \mathbf{g}_u \f$ and Hessian \f$ \mathbf{H}_u \f$ at every cell
      centre (fitDerivatives; exchanged across coupled patches). Pass 2:
      \f[
          \mathbf{g}_c = \frac{\sum_u w_u\,(\mathbf{g}_u
              + \mathbf{H}_u (\mathbf{x}_c - \mathbf{x}_u))}{\sum_u w_u},
          \qquad w_u = s_c (\psi_c - \psi_u) > 0,
      \f]
      over the face neighbours u upwind of c (closer to the zero set, or
      across it), \f$ s_c = \mathrm{sign}(\psi_c) \f$, and
      \f$ q = |\mathbf{g}_c| \f$, \f$ \mathbf{n} = \mathbf{g}_c/q \f$.
      A cell without an upwind neighbour (a local extremum of
      \f$ s\psi \f$, or \f$ \psi_c = 0 \f$) gets F = 0. Compact (face
      neighbours of face neighbours), per cell, no structured index.

    anchorInterface (optional, default false): a research probe. With true, F = 0
    in every cell that has a face neighbour (processor neighbours included)
    of the opposite sign or psi_c = 0, so the source cannot move a zero
    crossing between two cells; only the semi-Lagrangian step moves them.

    strainNormal (optional, default psi): a research probe. With fixedX the
    strain weight reads a = e_x.D.e_x instead of n.D.n, a smooth F of the same
    size that does not depend on psi (STATUS 11.26, the tangential-transport
    hypothesis).

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
                qScheme     centred;    // centred | upwindFit
                anchorInterface false;  // true: F = 0 in the interface cells
                strainNormal psi;       // psi | fixedX (probe: a = D_xx)
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

        //- How q is computed: centred | upwindFit
        const word qScheme_;

        //- F = 0 in the cells with a face neighbour of the opposite sign
        const bool anchorInterface_;

        //- The direction of the strain term: psi (n = g/q) | fixedX (e_x)
        const word strainNormal_;

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

        //- qScheme upwindFit: the gradient g_c from the fits of the upwind
        //  face neighbours, and whether c has an upwind neighbour
        void upwindFitGradient
        (
            volScalarField& psi,
            slReconstruction& fit,
            vectorField& gUp,
            boolList& valid
        ) const;


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
