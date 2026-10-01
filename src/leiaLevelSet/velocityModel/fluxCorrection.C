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

void Foam::fluxDivergenceReport
(
    const surfaceScalarField& phi,
    const char* label
)
{
    const fvMesh& mesh = phi.mesh();
    const scalarField& V = mesh.V().field();

    // div(phi) = sum_f phi_f / V_c; surfaceSum(mag(phi)) = sum_f |phi_f|.
    const scalarField divPhi(mag(fvc::div(phi)().primitiveField()));
    const scalarField sumMagPhi(fvc::surfaceSum(mag(phi))().primitiveField());

    const scalar maxDiv = gMax(divPhi);
    const scalar meanDiv = gSum(divPhi*V)/gSum(V);
    const scalar maxFraction =
        gMax(divPhi*V/max(sumMagPhi, scalarField(sumMagPhi.size(), VSMALL)));

    Info<< label << ": max|div(phi)| = " << maxDiv
        << " 1/s, mean|div(phi)| = " << meanDiv
        << " 1/s, max |sum phi_f|/sum|phi_f| = " << maxFraction << endl;
}


void Foam::correctFlux(surfaceScalarField& phi)
{
    const fvMesh& mesh = phi.mesh();
    const Time& runTime = phi.time();

    fluxDivergenceReport(phi, "Flux projection, before");

    // A COMPATIBLE Neumann problem. With zeroGradient p on every patch the
    // projection cannot change a boundary flux, so the sum of div(phi) over
    // the domain stays equal to the net boundary flux. If that is not zero,
    // the pinned reference cell absorbs all of it: MEASURED 2026-10-01 on the
    // polyhedral 3D shear meshes, whose top and bottom patches carry 0.136
    // m3/s in and out but whose face-centre sums differ by 4.8e-05 m3/s
    // (cfMesh meshes the two patches differently): after the projection every
    // other cell was divergence-free and the reference cell kept
    // |div(phi)| = 31 1/s. Scale the outflow of the non-coupled patches so
    // that it equals the inflow, as adjustPhi does for a closed pressure
    // problem. A no-op where the net boundary flux already vanishes (every
    // uniform hex mesh measured).
    {
        scalar massIn = 0;
        scalar massOut = 0;
        forAll(phi.boundaryField(), patchi)
        {
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
            if (phip.coupled())
            {
                continue;
            }
            forAll(phip, facei)
            {
                if (phip[facei] < 0)
                {
                    massIn -= phip[facei];
                }
                else
                {
                    massOut += phip[facei];
                }
            }
        }
        reduce(massIn, sumOp<scalar>());
        reduce(massOut, sumOp<scalar>());

        const scalar totalFlux = VSMALL + gSum(mag(phi.primitiveField()));
        const scalar imbalance = massOut - massIn;

        if (mag(imbalance) > SMALL*totalFlux)
        {
            if (massOut > VSMALL)
            {
                const scalar corr = massIn/massOut;
                auto& bphi = phi.boundaryFieldRef();
                forAll(bphi, patchi)
                {
                    fvsPatchScalarField& phip = bphi[patchi];
                    if (phip.coupled())
                    {
                        continue;
                    }
                    forAll(phip, facei)
                    {
                        if (phip[facei] > 0)
                        {
                            phip[facei] *= corr;
                        }
                    }
                }
                Info<< "Flux projection: net boundary flux " << imbalance
                    << " m3/s balanced by scaling the outflow by " << corr
                    << endl;
            }
            else
            {
                WarningInFunction
                    << "Net boundary flux " << imbalance << " m3/s with no"
                    << " outflow to scale: the reference cell absorbs it."
                    << endl;
            }
        }
    }

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
    // singular (defined only up to a constant). Pin the level at ONE cell so
    // the solver is non-singular; only grad(p) enters the projection, so the
    // choice of reference value is irrelevant. setReference is a no-op if a
    // patch already fixes p (it guards on p.needReference()) and for a
    // negative cell index. Without a reference: NaN.
    //
    // setReference(celli) pins the LOCAL cell celli on every rank that passes
    // celli >= 0. Passing 0 on every rank (the code until 2026-10-01) pinned
    // one cell PER RANK in parallel: the N pinned equations cannot all hold
    // for a Neumann problem, so the solution put sources at N - 1 cells and
    // the projected flux kept a divergence there, different for every
    // decomposition. The master rank alone now pins its cell 0, which is the
    // serial choice in a serial run.
    const label refCell = (Pstream::master() ? 0 : -1);
    pEqn.setReference(refCell, scalar(0));
    pEqn.solve();
    phi = phi - pEqn.flux();

    fluxDivergenceReport(phi, "Flux projection, after");

    if (runTime.writeTime())
    {
        fvc::div(phi)().write();
    }
}

// ************************************************************************* //
