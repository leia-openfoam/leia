/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

\*---------------------------------------------------------------------------*/

#include "fluxCorrection.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::correctFlux(surfaceScalarField& phi)
{
    const fvMesh& mesh = phi.mesh();
    const Time& runTime = phi.time();

    volScalarField Evol = fvc::div(phi);

    Info<< "gMax(fvc::div(phi)) = " << gMax(Evol) << endl;
    Info<< "gAverage(fvc::div(phi)) = " << gAverage(Evol) << endl;

    // Projection method enforcement of div(phi) = 0.
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("p", dimArea / dimTime, 0),
        "zeroGradient"
    );

    // https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
    // v = v_sol + grad p ; laplace p = div v ; v_sol = v - grad p
    fvScalarMatrix pEqn
    (
        fvm::laplacian(p) == fvc::div(phi)
    );
    // p has zeroGradient on every patch -> the Laplacian is pure-Neumann and
    // singular (defined only up to a constant). Pin the level at one cell so the
    // solver is non-singular; only grad(p) enters the projection, so the choice
    // of reference value is irrelevant. setReference is a no-op if a patch
    // already fixes p (it guards on p.needReference()). Without this: NaN.
    pEqn.setReference(0, scalar(0));
    pEqn.solve();
    phi = phi - pEqn.flux();
    Evol = fvc::div(phi);

    Info<< "Projection, gMax(fvc::div(phi)) = " << gMax(Evol) << endl;
    Info<< "Projection, gAverage(fvc::div(phi)) = " << gAverage(Evol) << endl;

    if (runTime.writeTime())
        Evol.write();
}

// ************************************************************************* //
