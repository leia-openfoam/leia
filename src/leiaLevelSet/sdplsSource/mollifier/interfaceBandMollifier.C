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

#include "interfaceBandMollifier.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSnGrad.H"
#include "coupledFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interfaceBandMollifier, false);
addToRunTimeSelectionTable(mollifier, interfaceBandMollifier, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interfaceBandMollifier::interfaceBandMollifier(const dictionary& dict)
:
    mollifier(dict),
    nLayers_(dict.getOrDefault<label>("nLayers", 3)),
    alphaName_(dict.getOrDefault<word>("alpha", "alpha"))
{
    if (nLayers_ < 0)
    {
        FatalErrorInFunction
            << "SDPLS band mollifier requires nLayers >= 0, got " << nLayers_
            << "." << nl
            << "nLayers is the number of dilation passes beyond the seed layer;"
            << " 0 leaves the source acting on the interface cells only."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * //

tmp<volScalarField> interfaceBandMollifier::field
(
    const volScalarField& psi
) const
{
    const fvMesh& mesh = psi.mesh();

    tmp<volScalarField> tG
    (
        new volScalarField
        (
            IOobject
            (
                word(),
                fileName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 0.0)
        )
    );
    volScalarField& G = tG.ref();

    // The phase indicator carries the topology. Fall back to psi's own sign
    // change if it is not registered -- that is still topological, and it keeps
    // this cut-off usable in a solver that never builds alpha.
    const bool haveAlpha = mesh.foundObject<volScalarField>(alphaName_);

    // `layer` holds the dilation pass at which a cell was first reached;
    // -1 means "not in the band".
    labelField layer(mesh.nCells(), -1);

    // ---- 1. SEED: cells adjacent to a face across which the phase changes.
    //
    // alpha is 0 or 1 in the bulk, so |snGrad(alpha)| > 0 marks exactly the
    // faces the interface passes between. This is topology: it does not ask how
    // far anything is from the interface, and so does not care what
    // |grad psi| is doing -- which is the whole point, since the |psi|-threshold
    // cut-off fails precisely where |grad psi| leaves 1.
    {
        const volScalarField& marker =
            haveAlpha ? mesh.lookupObject<volScalarField>(alphaName_) : psi;

        const surfaceScalarField sg(mag(fvc::snGrad(marker)));

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        forAll(sg, faceI)
        {
            if (sg[faceI] > SMALL)
            {
                layer[own[faceI]] = 0;
                layer[nei[faceI]] = 0;
            }
        }

        // Boundary faces, including COUPLED ones: a face on a processor patch
        // can straddle the interface just as an internal face can, and missing
        // it would leave a hole in the seed exactly at the boundary.
        forAll(sg.boundaryField(), patchI)
        {
            const fvsPatchScalarField& psg = sg.boundaryField()[patchI];
            const labelUList& fc = mesh.boundary()[patchI].faceCells();
            forAll(psg, i)
            {
                if (psg[i] > SMALL)
                {
                    layer[fc[i]] = 0;
                }
            }
        }
    }

    // ---- 2. DILATE, and do it across processor boundaries.
    //
    // Foam::neighbourNarrowBand dilates with mesh().cellCells(), which is LOCAL
    // topology: a cell one layer across a processor boundary from a marked cell
    // is never reached, the band is torn at every processor boundary, and the
    // result depends on the decomposition. Here the marker is carried in a
    // volScalarField so that evaluating its coupled patches performs the halo
    // exchange, and the spread is done face-by-face including coupled patches.
    for (label pass = 1; pass <= nLayers_; ++pass)
    {
        // Current membership as a field, so coupled patches can be evaluated.
        volScalarField inBand
        (
            IOobject
            (
                word(), fileName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false
            ),
            mesh,
            dimensionedScalar(dimless, 0.0),
            "zeroGradient"
        );
        forAll(layer, cellI)
        {
            inBand[cellI] = (layer[cellI] >= 0) ? 1.0 : 0.0;
        }
        inBand.correctBoundaryConditions();

        labelField next(layer);

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        forAll(own, faceI)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            if (layer[o] >= 0 && next[n] < 0) next[n] = pass;
            if (layer[n] >= 0 && next[o] < 0) next[o] = pass;
        }

        // THE PARALLEL STEP. patchNeighbourField() is the membership on the
        // OTHER side of a processor patch; without this the dilation simply
        // stops at the boundary.
        forAll(inBand.boundaryField(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            if (!isA<coupledFvPatch>(patch)) continue;

            const scalarField neiField
            (
                inBand.boundaryField()[patchI].patchNeighbourField()
            );
            const labelUList& fc = patch.faceCells();
            forAll(neiField, i)
            {
                if (neiField[i] > 0.5 && next[fc[i]] < 0)
                {
                    next[fc[i]] = pass;
                }
            }
        }

        layer = next;
    }

    // ---- 3. TAPER: unity on the seed, smoothly to zero one layer past the
    // last. A raised cosine is C^1 at both ends, and -- this is the part that
    // matters for SDPLS -- it is FLAT on the seed layer, so the cut-off adds no
    // variation of the effective source coefficient across the two cells
    // straddling psi = 0. That variation is the only thing that displaces the
    // discrete zero crossing, so a taper with a kink at the interface would
    // introduce the very interface motion the source is constructed to avoid.
    const scalar denom = scalar(nLayers_) + 1.0;
    forAll(G, cellI)
    {
        const label l = layer[cellI];
        if (l < 0)
        {
            G[cellI] = 0.0;
        }
        else if (l == 0)
        {
            G[cellI] = 1.0;
        }
        else
        {
            G[cellI] = 0.5*(1.0 + Foam::cos(constant::mathematical::pi*l/denom));
        }
    }
    G.correctBoundaryConditions();

    return tG;
}

} // End namespace Foam

// ************************************************************************* //
