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
    Foam::rescalingRedistancer

Description
    Algebraic rescaling of the level set (D. Bothe, discussion of 2026-10,
    the notes "Dieter-Besprechung-LevelSet"): a pointwise map psi -> c psi
    with a POSITIVE factor c, so the sign of psi, the zero set and the normal
    on the zero set stay exactly as they are.

      order first   (0a)  r    = psi/|grad psi|
                          r = d - (B/(2A)) d^2 + O(d^3), A = |grad psi|,
                          B = d_nn psi on the zero set: a first-order
                          distance (exact value and slope on the zero set)
      order second  (0b)  dhat = 2 psi/(|grad psi|
                                 + sqrt((|grad psi|^2 - 2 psi d_nn psi)^+))
                          the stable root of psi - d d_n psi + d^2/2 d_nn psi
                          = 0, the Taylor expansion of psi(x - d n) = 0:
                          dhat = d + O(d^3), a second-order distance; exact
                          for psi = r^2 - R^2

    Both return d unchanged for psi = d (|grad psi| = 1, d_nn psi = 0).

    The UNSTRUCTURED-FVM version: grad psi = fvc::grad(psi, gradScheme), and
    the Hessian H = fvc::grad(grad psi, hessianScheme), the FVM gradient of
    the FVM gradient, d_nn psi = n.H.n with n = grad psi/|grad psi|. This
    Hessian is consistent on uniform hexahedra and is a wide, noise-amplifying
    operator: an error eps in psi gives an error of order eps/h^2 in H. The
    model exists also to measure that (author decision 2026-10-29).

    Dictionary (fvSolution levelSet.redistancer, next to the base keys
    trigger, redistanceInterval, ...):
        type            rescaling;
        order           second;          // first | second
        region          band;            // band | all
        nLayers         3;               // band = sign-change cells + nLayers
        maxEvents       -1;              // -1: no limit; 1: rescale once
        gradScheme      rescaleGradPsi;  // fvSchemes gradSchemes entry
        hessianScheme   rescaleHessian;  // fvSchemes gradSchemes entry
        positivePart    max;             // max | smooth
        smoothEps       0;               // smooth: (s + sqrt(s^2 + eps^2))/2
        gradFloor       1e-3;            // |grad psi| below it: cell unchanged
        magGrad         centred;         // centred | godunov (Rouy-Tourin)

    magGrad godunov: the one-sided |grad psi| does not vanish at a
    grid-scale extremum. Linearised about d, the centred first-order map is
    one step of the linearQ source with mu dt = 1, and it doubles the mode
    (-1)^(i+j) d per application (MEASURED 2026-10-29: applied every step it
    diverges, with Crank-Nicolson and with Euler).

    Every event logs the cells changed, the cells left by the gradient floor,
    the cells with a negative discriminant, and the largest |change|.

SourceFiles
    rescalingRedistancer.C

\*---------------------------------------------------------------------------*/

#ifndef rescalingRedistancer_H
#define rescalingRedistancer_H

#include "redistancer.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class rescalingRedistancer Declaration
\*---------------------------------------------------------------------------*/

class rescalingRedistancer
:
    public redistancer
{
    // Private Data

        //- 1: psi/|grad psi|; 2: the quadratic-root distance
        label order_;

        //- band | all
        word region_;

        //- Dilation layers of the band beyond the sign-change cells
        label nLayers_;

        //- Largest number of events; -1 = no limit
        label maxEvents_;

        //- Events so far
        label nEvents_;

        //- fvSchemes names of the gradient and of the Hessian
        word gradScheme_;
        word hessianScheme_;

        //- max | smooth, and the width of the smooth positive part
        word positivePart_;
        scalar smoothEps_;

        //- |grad psi| below this value: the cell is left unchanged
        scalar gradFloor_;

        //- centred (|fvc::grad psi|) | godunov (Rouy-Tourin one-sided)
        word magGrad_;


    // Private Member Functions

        //- Rouy-Tourin |grad psi|: per cell and axis the largest one-sided
        //  slope toward the zero set (the formula of
        //  pdeRedistancer::magGradGodunov, coupled patches included)
        tmp<scalarField> magGradGodunov(const volScalarField& psi) const;


public:

    //- Runtime type information
    TypeName("rescaling");


    // Constructors

        //- Construct from the mesh (reads fvSolution levelSet.redistancer)
        rescalingRedistancer(const fvMesh& mesh);


    //- Destructor
    virtual ~rescalingRedistancer() = default;


    // Member Functions

        //- Rescale psi in the region; false when the event limit is reached
        virtual bool doRedistance(volScalarField& psi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
