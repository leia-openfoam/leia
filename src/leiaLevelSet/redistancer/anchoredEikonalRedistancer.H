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
    Foam::anchoredEikonalRedistancer

Description
    Plane-anchored steady-Eikonal redistancer. The narrow-band (interface)
    cells and their first ring receive Dirichlet signed distances read from
    the NORMALIZED linear least-squares level-set planes (planeAnchors.H --
    the same planes the Detrixhe-Aslam indicator and the geometric mass flux
    consume, so the redistanced psi is realigned to the discrete interface
    the solver already sees). The remaining field solves the steady Eikonal
    equation |grad psi| = 1 in Tucker's regularized advection-diffusion
    linearization, iterated with the transport direction frozen from the
    previous iterate (the discretization of OpenFOAM's
    wallDist advectionDiffusion patch-distance method, applied here with
    fvMatrix::setValues internal Dirichlet anchors instead of wall
    boundaries):

        div(nHat y) - y div(nHat) - epsilon y laplacian(y) = 1,
        nHat = grad(y)/|grad(y)|,

    solved for the UNSIGNED distance y (initialized |psi|, anchored |dist|),
    then re-signed with sign(psi_old); anchored cells receive their signed
    plane distances exactly.

    fvSolution levelSet.redistancer entries (defaults = OpenFOAM's
    advectionDiffusion): epsilon 0.1; tolerance 1e-3; maxIter 10.
    Requires: an fvSolution solver entry matching "psiDist" (use the regex
    key "psiDist.*"; the matrix is ASYMMETRIC -- PBiCGStab/DILU or
    smoothSolver/GaussSeidel, never PCG/DIC), an fvSchemes entry
    div(yPhi,psiDist), and optionally relaxationFactors equations
    "psiDist.*" (without it the equation relaxation is a no-op).

SourceFiles
    anchoredEikonalRedistancer.C

\*---------------------------------------------------------------------------*/

#ifndef anchoredEikonalRedistancer_H
#define anchoredEikonalRedistancer_H

#include "redistancer.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class anchoredEikonalRedistancer Declaration
\*---------------------------------------------------------------------------*/

class anchoredEikonalRedistancer
:
    public redistancer
{
private:

    const fvMesh& mesh_;

    //- Diffusive regularization factor (Tucker; OF default 0.1)
    const scalar epsilon_;

    //- Iteration convergence tolerance on the initial residual
    const scalar tolerance_;

    //- Maximum number of frozen-direction iterations
    const label maxIter_;

public:

    TypeName("anchoredEikonal");

    anchoredEikonalRedistancer(const fvMesh& mesh);

    virtual ~anchoredEikonalRedistancer() = default;

    virtual bool doRedistance(volScalarField& psi);
};

} // End namespace Foam

#endif

// ************************************************************************* //
