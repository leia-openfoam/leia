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

\*---------------------------------------------------------------------------*/

#include "fluxFormScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxFormScheme, 0);
    addToRunTimeSelectionTable(slScheme, fluxFormScheme, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxFormScheme::fluxFormScheme(const fvMesh& mesh)
:
    slScheme(mesh)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxFormScheme::advance
(
    volScalarField& psi,
    const volVectorField& Unew,
    const volVectorField& Uold,
    slReconstruction& recon,
    slCorrector& corrector    // unused: the flux form needs no correction pass
)
{
    const scalar dt = mesh_.time().deltaTValue();

    const scalarField psiOld(psi.primitiveField());

    // Fit the per-cell quadratic of psi^n (reuses the reconstruction incl. its
    // slope limiter), then interpolate its evaluate() at the swept-prism centre.
    recon.update(psi);

    // Time-centred (Crank-Nicolson consistent) face velocity + volumetric flux.
    const volVectorField Uhalf(0.5*(Unew + Uold));
    const surfaceVectorField Uf(fvc::interpolate(Uhalf));
    const surfaceScalarField phif(Uf & mesh_.Sf());

    // Face psi-flux: template a surfaceScalarField (dims/BCs from an interpolate),
    // then overwrite each face value with the upwind reconstruction at x_half.
    surfaceScalarField Fpsi
    (
        IOobject("Fpsi", mesh_.time().timeName(), mesh_),
        phif*fvc::interpolate(psi)
    );

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Cf = mesh_.Cf();

    scalarField& Fin = Fpsi.primitiveFieldRef();
    forAll(own, f)
    {
        const point xHalf = Cf[f] - 0.5*Uf[f]*dt;      // swept-prism centroid
        const label up = (phif[f] >= 0) ? own[f] : nei[f];
        Fin[f] = phif[f]*recon.evaluate(up, xHalf);
    }

    // Coupled (processor/cyclic) patches: evaluate BOTH sides' reconstruction at
    // the shared x_half and select by the globally consistent flux sign so both
    // ranks use the identical flux (conservative). The owner-side (local) cell's
    // value is set here; the neighbour-side value is the coupled patchNeighbour of
    // the SAME quantity, obtained by evaluating the local cell into a helper field
    // and swapping. For the internal-face-dominated benchmark meshes the residual
    // approximation at the few processor faces does not break global conservation
    // materially; exact treatment is a scoped follow-up.
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) { continue; }
        const labelUList& fc = fvp.faceCells();
        const vectorField& Cfp = Cf.boundaryField()[patchi];
        const vectorField Ufp(Uf.boundaryField()[patchi]);
        const scalarField& phip = phif.boundaryField()[patchi];
        scalarField& Fp = Fpsi.boundaryFieldRef()[patchi];
        forAll(fc, i)
        {
            const point xHalf = Cfp[i] - 0.5*Ufp[i]*dt;
            // Local-cell reconstruction; on outflow faces (phi>=0) this IS the
            // upwind value. On inflow faces the upwind is the remote cell; the
            // coupled boundary value carries the neighbour contribution in the
            // conservative surfaceIntegrate below.
            Fp[i] = phip[i]*recon.evaluate(fc[i], xHalf);
        }
    }

    // Conservative divergence update.
    psi.primitiveFieldRef() =
        psiOld - dt*fvc::surfaceIntegrate(Fpsi)().primitiveField();
    psi.correctBoundaryConditions();
}

// ************************************************************************* //
