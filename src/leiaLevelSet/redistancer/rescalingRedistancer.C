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

#include "rescalingRedistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "rescaledDistance.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rescalingRedistancer, false);
    addToRunTimeSelectionTable(redistancer, rescalingRedistancer, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rescalingRedistancer::rescalingRedistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    order_(0),
    region_(redistDict_.getOrDefault<word>("region", "band")),
    nLayers_(redistDict_.getOrDefault<label>("nLayers", 3)),
    maxEvents_(redistDict_.getOrDefault<label>("maxEvents", -1)),
    nEvents_(0),
    gradScheme_(redistDict_.getOrDefault<word>("gradScheme", "rescaleGradPsi")),
    hessianScheme_
    (
        redistDict_.getOrDefault<word>("hessianScheme", "rescaleHessian")
    ),
    positivePart_(redistDict_.getOrDefault<word>("positivePart", "max")),
    smoothEps_(redistDict_.getOrDefault<scalar>("smoothEps", 0)),
    gradFloor_(redistDict_.getOrDefault<scalar>("gradFloor", 1e-3)),
    magGrad_(redistDict_.getOrDefault<word>("magGrad", "centred"))
{
    if (magGrad_ != "centred" && magGrad_ != "godunov")
    {
        FatalIOErrorInFunction(redistDict_)
            << "rescaling: magGrad must be 'centred' or 'godunov', got '"
            << magGrad_ << "'" << exit(FatalIOError);
    }
    const word order(redistDict_.getOrDefault<word>("order", "second"));
    if (order == "first")
    {
        order_ = 1;
    }
    else if (order == "second")
    {
        order_ = 2;
    }
    else
    {
        FatalIOErrorInFunction(redistDict_)
            << "rescaling: order must be 'first' or 'second', got '"
            << order << "'" << exit(FatalIOError);
    }
    if (region_ != "band" && region_ != "all")
    {
        FatalIOErrorInFunction(redistDict_)
            << "rescaling: region must be 'band' or 'all', got '"
            << region_ << "'" << exit(FatalIOError);
    }
    if (positivePart_ != "max" && positivePart_ != "smooth")
    {
        FatalIOErrorInFunction(redistDict_)
            << "rescaling: positivePart must be 'max' or 'smooth', got '"
            << positivePart_ << "'" << exit(FatalIOError);
    }
    Info<< "rescaling redistancer: order " << order << ", region "
        << region_ << " (nLayers " << nLayers_ << "), maxEvents "
        << maxEvents_ << ", grad " << gradScheme_ << ", Hessian "
        << hessianScheme_ << ", positive part " << positivePart_
        << ", gradFloor " << gradFloor_ << ", |grad psi| " << magGrad_ << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::rescalingRedistancer::magGradGodunov(const volScalarField& psi) const
{
    const fvMesh& mesh = psi.mesh();
    vectorField gAxis(mesh.nCells(), Zero);
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const volVectorField& C = mesh.C();

    auto accumulate = [&gAxis](const label celli, const scalar q, const vector& dhat)
    {
        if (q <= 0) return;
        vector& g = gAxis[celli];
        for (direction a = 0; a < 3; ++a)
        {
            g[a] = Foam::max(g[a], q*Foam::mag(dhat[a]));
        }
    };
    auto sgn = [](const scalar s) { return s > 0 ? 1.0 : (s < 0 ? -1.0 : 0.0); };

    forAll(own, facei)
    {
        const label c = own[facei];
        const label n = nei[facei];
        const vector d = C[n] - C[c];
        const scalar md = Foam::mag(d);
        if (md < SMALL) continue;
        const vector dhat = d/md;
        accumulate(c, sgn(psi[c])*(psi[c] - psi[n])/md, dhat);
        accumulate(n, sgn(psi[n])*(psi[n] - psi[c])/md, dhat);
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (!p.coupled()) continue;
        const scalarField psiNbr(psi.boundaryField()[patchi].patchNeighbourField());
        const vectorField deltas(p.delta());
        const labelUList& fc = p.faceCells();
        forAll(fc, i)
        {
            const scalar md = Foam::mag(deltas[i]);
            if (md < SMALL) continue;
            accumulate(fc[i], sgn(psi[fc[i]])*(psi[fc[i]] - psiNbr[i])/md, deltas[i]/md);
        }
    }
    auto tq = tmp<scalarField>::New(mesh.nCells());
    scalarField& q = tq.ref();
    forAll(q, celli)
    {
        q[celli] = Foam::mag(gAxis[celli]);
    }
    return tq;
}


bool Foam::rescalingRedistancer::doRedistance(volScalarField& psi)
{
    if (maxEvents_ >= 0 && nEvents_ >= maxEvents_)
    {
        return false;
    }

    const fvMesh& mesh = psi.mesh();

    // The region: the topological band (parallel-correct dilation of the
    // registered sign-change band) or every cell.
    boolList inRegion(mesh.nCells(), true);
    if (region_ == "band")
    {
        topologicalMask(psi, nLayers_, inRegion);
    }

    // FVM derivatives. The gradient's coupled patches are evaluated before
    // the Hessian reads them, so the Hessian is decomposition-independent.
    volVectorField gradPsi(fvc::grad(psi, gradScheme_));
    gradPsi.correctBoundaryConditions();

    tmp<volTensorField> tH;
    if (order_ == 2)
    {
        tH = fvc::grad(gradPsi, hessianScheme_);
    }

    psi.correctBoundaryConditions();   // coupled values for the one-sided slopes
    tmp<scalarField> tqG;
    if (magGrad_ == "godunov")
    {
        tqG = magGradGodunov(psi);
    }

    const scalarField psiOld(psi.primitiveField());
    scalarField& p = psi.primitiveFieldRef();

    label nChanged = 0, nFloor = 0, nNegDisc = 0;
    scalar maxChange = 0;

    forAll(p, celli)
    {
        if (!inRegion[celli])
        {
            continue;
        }

        const vector& g = gradPsi[celli];
        const scalar qc = mag(g);
        const scalar q = (magGrad_ == "godunov") ? tqG()[celli] : qc;
        const scalar ps = psiOld[celli];
        const scalar dnn =
            (order_ == 2 && qc > 0) ? ((g/qc) & (tH()[celli] & (g/qc))) : scalar(0);
        scalar newPsi = ps;
        bool negDisc = false;
        if
        (
            !rescaledDistance::estimate
            (
                order_, ps, q, dnn, gradFloor_, positivePart_ == "smooth",
                smoothEps_, newPsi, negDisc
            )
        )
        {
            ++nFloor;
            continue;
        }
        if (negDisc)
        {
            ++nNegDisc;
        }

        maxChange = max(maxChange, mag(newPsi - ps));
        p[celli] = newPsi;
        ++nChanged;
    }

    psi.correctBoundaryConditions();

    reduce(nChanged, sumOp<label>());
    reduce(nFloor, sumOp<label>());
    reduce(nNegDisc, sumOp<label>());
    reduce(maxChange, maxOp<scalar>());

    ++nEvents_;

    Info<< "rescaling: event " << nEvents_ << " at t = "
        << mesh.time().timeOutputValue() << ": " << nChanged
        << " cells rescaled, " << nFloor << " below gradFloor, " << nNegDisc
        << " with a negative discriminant, max |change| " << maxChange
        << endl;

    return true;
}


// ************************************************************************* //
