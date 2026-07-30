/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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

#include "pdeRedistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
// #include "fvCFD.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * *  Local Functions  * * * * * * * * * * * * * * //

static Foam::tmp<Foam::volScalarField> sign_smoothed(const Foam::volScalarField& field)
{
    using namespace Foam;
    return tmp<volScalarField>
        (
            new volScalarField
            (
                field/sqrt(pow(field,2) + pow(dimensioned(field.dimensions(), SMALL),2))
            )
        ); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdeRedistancer, false);
addToRunTimeSelectionTable(redistancer, pdeRedistancer, Mesh);

} // End namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdeRedistancer::pdeRedistancer(const fvMesh& mesh)
    :
        redistancer(mesh),
        deltaT_(redistDict_.getOrDefault<scalar>("deltaT",0.1)),
        deltaTCoeff_(redistDict_.getOrDefault<scalar>("deltaTCoeff",-1.0)),
        ninterations_(redistDict_.getOrDefault<uint>("niterations",10)),
        freezeBand_(redistDict_.getOrDefault<Switch>("freezeBand", false)),
        hamiltonian_(redistDict_.getOrDefault<word>("hamiltonian", "central")),
        write_(redistDict_.getOrDefault<bool>("write",false))
{
    if (hamiltonian_ != "central" && hamiltonian_ != "godunov")
    {
        FatalIOErrorInFunction(redistDict_)
            << "pdeRedistancer hamiltonian must be 'central' or 'godunov', "
            << "got '" << hamiltonian_ << "'" << exit(FatalIOError);
    }
    if (freezeBand_)
    {
        Info<< "pdeRedistancer: frozen-band bulk-only fill (anchorLayers = "
            << anchorLayers_ << ", hamiltonian = " << hamiltonian_ << ")"
            << endl;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::pdeRedistancer::magGradGodunov
(
    const volScalarField& psi_p,
    const volScalarField& signPsi
) const
{
    // Rouy-Tourin upwind |grad psi|: per cell, per Cartesian axis, the
    // maximum one-sided difference TOWARD the interface (q > 0 <=> the
    // neighbour is closer to the zero set). Monotone and causal: values can
    // only be pulled toward the anchors' distance cone, never pushed into
    // the runaway the central gradient produces. Exact Rouy-Tourin on hex
    // (face deltas axis-aligned); the |dhat_a| weights degrade gracefully on
    // skewed faces. Physical boundary faces carry no information (skipped).
    const fvMesh& mesh = psi_p.mesh();

    auto tmg = tmp<volScalarField>::New
    (
        IOobject
        (
            "magGradGodunov",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(psi_p.dimensions()/dimLength, Zero),
        word("zeroGradient")
    );
    volScalarField& mg = tmg.ref();

    vectorField gAxis(mesh.nCells(), Zero);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const volVectorField& C = mesh.C();

    auto accumulate = [&gAxis]
    (
        const label celli,
        const scalar q,          // one-sided slope toward the interface
        const vector& dhat
    )
    {
        if (q <= 0) return;
        vector& g = gAxis[celli];
        for (direction a = 0; a < 3; ++a)
        {
            g[a] = Foam::max(g[a], q*Foam::mag(dhat[a]));
        }
    };

    forAll(own, facei)
    {
        const label c = own[facei];
        const label n = nei[facei];
        const vector d = C[n] - C[c];
        const scalar md = Foam::mag(d);
        if (md < SMALL) continue;
        const vector dhat = d/md;
        accumulate(c, signPsi[c]*(psi_p[c] - psi_p[n])/md, dhat);
        accumulate(n, signPsi[n]*(psi_p[n] - psi_p[c])/md, dhat);
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (!p.coupled()) continue;
        const scalarField psiNbr
        (
            psi_p.boundaryField()[patchi].patchNeighbourField()
        );
        const vectorField deltas(p.delta());
        const labelUList& fc = p.faceCells();
        forAll(fc, i)
        {
            const label c = fc[i];
            const scalar md = Foam::mag(deltas[i]);
            if (md < SMALL) continue;
            accumulate
            (
                c,
                signPsi[c]*(psi_p[c] - psiNbr[i])/md,
                deltas[i]/md
            );
        }
    }

    scalarField& mgIn = mg.primitiveFieldRef();
    forAll(mgIn, celli)
    {
        mgIn[celli] = Foam::mag(gAxis[celli]);
    }
    mg.correctBoundaryConditions();

    return tmg;
}


bool Foam::pdeRedistancer::doRedistance(volScalarField& psi)
{
    // Mesh-relative pseudo-time step (see header): the explicit transport
    // at unit characteristic speed is CFL-limited by the smallest cell.
    scalar deltaT = deltaT_;
    if (deltaTCoeff_ > 0)
    {
        const scalar h =
            1.0/gMax(psi.mesh().deltaCoeffs().primitiveField());
        deltaT = deltaTCoeff_*h;
    }

    // Working copy ON THE REAL MESH: it inherits psi's boundary conditions,
    // including PROCESSOR patches. (The historical pseudo-Time/pseudo-mesh
    // detour gave the working field 'extrapolatedCalculated' on ALL patches
    // -- undefined coupling in parallel: gradients near processor boundaries
    // read unsynchronized garbage and the trapped SIGFPE killed whole
    // studies. The diagonal fvm::ddt(psi_p) == Su pseudo-solve it performed
    // is algebraically the explicit Euler update below.)
    volScalarField psi_p("psiRedist", psi);

    // Re-sign source frozen from the input field (as before).
    const volScalarField signPsi(sign(psi));

    // Overflow guard: no signed distance can exceed the mesh bounding-box
    // diagonal. The non-Godunov (central-gradient) Hamiltonian oscillates on
    // saturated profiles and can run away to floating-point overflow;
    // clamping keeps a diverging event finite, so its (measured) failure
    // shows up in the error norms and the base-class post-event warning
    // instead of a crash.
    const scalar Lmax = mag(psi.mesh().bounds().span());

    // Frozen-band variant: the sign-change band + anchorLayers_ guard layers
    // keep their transported values (Dirichlet anchors); only the bulk is
    // updated. The guard thickness keeps the fill's error out of the band
    // cells' reconstruction stencils.
    boolList frozen;
    if (freezeBand_)
    {
        topologicalMask(psi, anchorLayers_, frozen);
        label nFrozen = 0;
        forAll(frozen, celli)
        {
            if (frozen[celli]) ++nFrozen;
        }
        reduce(nFrozen, sumOp<label>());
        Info<< "pdeRedistancer: frozen-band event, " << nFrozen
            << " anchor cells frozen, "
            << returnReduce(psi.mesh().nCells(), sumOp<label>()) - nFrozen
            << " bulk cells filled" << endl;
    }

    for (uint iter = 0; iter < ninterations_; ++iter)
    {
        // central: the fvSchemes grad(psi) operator (existing behaviour).
        // godunov: the monotone Rouy-Tourin upwind Hamiltonian.
        const volScalarField magGradPsi
        (
            (hamiltonian_ == "godunov")
          ? magGradGodunov(psi_p, signPsi)
          : tmp<volScalarField>(mag(fvc::grad(psi_p, "grad(psi)")))
        );

        if (freezeBand_)
        {
            // Two guards, both specific to the bulk re-extension problem:
            //  (1) explicit stability -- a transported field reaches
            //      |grad psi| ~ 6 in compression zones, so the raw update
            //      deltaT*(1-|grad|) can jump 2.5h in ONE pseudo-step
            //      (measured: overshoot through zero -> spurious phase 3x
            //      the droplet volume). Clamp the step to half a cell.
            //  (2) sign preservation -- the interface lives EXCLUSIVELY in
            //      the frozen anchors; a bulk cell that changes sign IS a
            //      spurious interface by definition. s*psi >= 0 is the
            //      mathematical constraint of re-extension, not a limiter.
            const scalar h =
                1.0/gMax(psi.mesh().deltaCoeffs().primitiveField());
            const scalar dpsiMax = 0.5*h;

            scalarField& p = psi_p.primitiveFieldRef();
            const scalarField& s = signPsi.primitiveField();
            const scalarField& g = magGradPsi.primitiveField();
            forAll(p, celli)
            {
                if (frozen[celli]) continue;
                const scalar dpsi =
                    clamp
                    (
                        deltaT*s[celli]*(1.0 - g[celli]),
                        -dpsiMax, dpsiMax
                    );
                scalar pNew = p[celli] + dpsi;
                if (s[celli]*pNew < SMALL)
                {
                    pNew = s[celli]*SMALL;
                }
                p[celli] = pNew;
            }
        }
        else
        {
            psi_p.primitiveFieldRef() +=
                deltaT
               *signPsi.primitiveField()
               *(1.0 - magGradPsi.primitiveField());
        }

        psi_p.clamp_range(-Lmax, Lmax);
        psi_p.correctBoundaryConditions();
    }

    psi.primitiveFieldRef() = psi_p.primitiveField();
    psi.correctBoundaryConditions();

    if (write_)
    {
        write(psi_p);
    }

    return true;
}

void Foam::pdeRedistancer::write(const volScalarField& psi) const
{
    psi.write();
    volVectorField gradPsi = fvc::grad(psi);
    gradPsi.rename("gradPsi");
    gradPsi.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
