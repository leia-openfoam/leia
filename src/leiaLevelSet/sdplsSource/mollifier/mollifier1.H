/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt
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
    Foam::mollifier1

Description
    The SDPLS cut-off (mollifier) function \f$ G \f$ of Fricke, Marić, Vučković,
    Roisman and Bothe, "A locally signed-distance preserving level set method
    (SDPLS) for moving interfaces", arXiv:2208.01269
    \cite fricke_locally_2022, <b>equation (24)</b>:

    \f{equation}{
        G(x) =
        \begin{cases}
        1 & \text{if} \quad 0 \le x \le w_1,\\
        \exp \left( - \ln(10^3) \frac{(x-w_1)^2}{(w_2-w_1)^2} \right) & \text{if} \quad w_1 < x,\\
        G(-x) & \text{if} \quad x<0
        \end{cases}
    \f}

    It multiplies the SDPLS source coefficient, which is <b>equation (23)</b> of
    the same reference:

    \f{equation}{
        r^n_{ijk} = - (\nabla v) \frac{\nabla \phi}{|\nabla \phi| + \epsilon}
                      \cdot \frac{\nabla \phi}{|\nabla \phi| + \epsilon}
                      \, G(\phi^n_{ijk}) ,
    \f}

    so the source is forced to zero sufficiently far from the zero contour. The
    reference introduces it to keep a normal-field discontinuity at inflow
    boundaries -- produced by the homogeneous Neumann condition on \f$ \phi \f$
    -- from degenerating the regularity of \f$ r(x) \f$ (Sec. 3, item (ii)).

    There is a second, independent reason to use it, measured in this repository:
    without a cut-off the source acts wherever \f$ \phi \f$ is large and the
    normal strain does not cancel, and it amplifies \f$ \phi \f$ in proportion to
    itself. On the 2D steady vortex the global \f$ \max ||\nabla\phi|-1| \f$
    reaches 119.6 at \f$ N=256 \f$ and \e grows under refinement (fitted order
    -1.44) while the narrow-band error stays controlled.

    Two structural properties of \f$ G \f$ carry weight and must survive any
    replacement:
      - it is \b even in \f$ \phi \f$, and \b flat (identically unity) on
        \f$ [0, w_1] \f$. The discrete zero crossing of \f$ \phi \f$ lies
        between cell centres and is displaced only by a \e variation of the
        effective source coefficient across the two straddling cells. Because
        \f$ G \equiv 1 \f$ there, it contributes none. A taper with a kink at
        \f$ \phi = 0 \f$ (\c tanh of \f$ |\phi| \f$, say) would introduce exactly
        the interface motion the SDPLS source is constructed not to have.
      - \f$ w_1 \f$ and \f$ w_2 \f$ are \b physical lengths, not multiples of the
        cell size. Scaling them with \f$ h \f$ makes the modified PDE
        mesh-dependent, so a convergence order fitted across such a ladder is not
        the order of any one continuum problem. The reference uses
        \f$ w_1 = 0.05 \f$, \f$ w_2 = 0.15 \f$ for an initial radius
        \f$ R_0 = 0.3 \f$ (Sec. 4).

    Selected at run time by the name \c m1, not by the class name:
    \verbatim
        sdplsSource
        {
            type    R;
            mollifier
            {
                type    m1;
                w1      0.05;
                w2      0.15;   // must be > w1
            }
        }
    \endverbatim

See also
    Foam::mollifier (the \c none cut-off, which leaves the source unscaled)

SourceFiles
    mollifier1.C

\*---------------------------------------------------------------------------*/

#ifndef mollifier1_H
#define mollifier1_H

#include "mollifier.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class mollifier1 Declaration
\*---------------------------------------------------------------------------*/

class mollifier1
    :   public mollifier
{
private:

    // Private Data Members

        //- Half-width of the plateau: \f$ G \equiv 1 \f$ for
        //  \f$ |\phi| \le w_1 \f$. A physical length. Eq. (24) of
        //  \cite fricke_locally_2022.
        scalar w1_;

        //- Distance at which the cut-off has decayed to \f$ 10^{-3} \f$, which
        //  is what sets the Gaussian width \f$ (w_2 - w_1) \f$ in eq. (24).
        //  A physical length, and necessarily \f$ > w_1 \f$ -- the constructor
        //  rejects anything else. (The \f$ 10^{-3} \f$ is the VALUE \f$ G \f$
        //  attains at \f$ w_2 \f$, not a default for \f$ w_2 \f$ itself.)
        scalar w2_;
public:

    // Static Data Members
    
        TypeName("m1");

    // Constructors

        //- Construct from fvSolution
        mollifier1(const dictionary& dict);
    
    //- Destructor
    virtual ~mollifier1() = default;

protected:

    // Protected Member functions
        
        //- Evaluate the cut-off \f$ G \f$ of eq. (24) in
        //  \cite fricke_locally_2022 (arXiv:2208.01269) for one cell.
        //
        //  @param x the level set value \f$ \phi \f$ in that cell. Only
        //         \f$ |x| \f$ matters: \f$ G \f$ is even by construction, and
        //         the negative branch recurses on \f$ -x \f$.
        //  @return \f$ G(x) \in (0, 1] \f$, unity on the plateau
        //          \f$ |x| \le w_1 \f$.
        virtual double mollify(double x) const;


};

} // End namespace Foam

#endif

// ************************************************************************* //
