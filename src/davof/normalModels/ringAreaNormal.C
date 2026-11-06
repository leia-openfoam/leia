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

#include "ringAreaNormal.H"
#include "davof.H"
#include "zoneDistribute.H"
#include "syncTools.H"
#include "PstreamBuffers.H"
#include "HashSet.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace davofNormalModels
{
    defineTypeNameAndDebug(ringArea, 0);
    addToRunTimeSelectionTable(davofNormalModel, ringArea, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofNormalModels::ringArea::ringArea
(
    const dictionary& dict,
    const reconstruction::davof& scheme
)
:
    davofNormalModel(typeName, dict, scheme, "all"),
    mInterface_
    (
        IOobject
        (
            IOobject::groupName("davofRingArea", scheme.alpha1().group()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimArea, Zero)
    ),
    filled_(),
    stencil_(dict.getOrDefault<word>("ringStencil", "one")),
    kernel_(dict.getOrDefault<word>("ringKernel", "none")),
    kernelSupport_(dict.getOrDefault<scalar>("ringKernelSupport", 2)),
    angleFilter_(dict.getOrDefault<scalar>("ringAngleFilter", 0)),
    cosAngleFilter_(Foam::cos(degToRad(angleFilter_))),
    planeFit_(dict.getOrDefault<word>("ringPlaneFit", "none")),
    planeFitTol_(dict.getOrDefault<scalar>("ringPlaneFitTol", 1e-10)),
    nAngleDropped_(0),
    nAngleKept_(0),
    nRemoteRings_(0)
{
    if (stencil_ != "one" && stencil_ != "two")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringArea: unknown ringStencil '" << stencil_
            << "'. Valid: one, two." << exit(FatalIOError);
    }
    if (kernel_ != "none" && kernel_ != "gaussian" && kernel_ != "wendland")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringArea: unknown ringKernel '" << kernel_
            << "'. Valid: none, gaussian, wendland." << exit(FatalIOError);
    }
    if (planeFit_ != "none" && planeFit_ != "mean" && planeFit_ != "own")
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringArea: unknown ringPlaneFit '" << planeFit_
            << "'. Valid: none, mean, own." << exit(FatalIOError);
    }
    if (!(kernelSupport_ > 0))
    {
        FatalIOErrorInFunction(dict)
            << "normalModel ringArea: ringKernelSupport must be positive"
            << exit(FatalIOError);
    }
    if (!production())
    {
        Info<< "davof: ringArea stencil variant: ringStencil " << stencil_
            << ", ringKernel " << kernel_ << " (support " << kernelSupport_
            << " h), ringAngleFilter " << angleFilter_ << " deg, ringPlaneFit "
            << planeFit_ << " (tol " << planeFitTol_ << ")" << endl;
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::davofNormalModels::ringArea::fillInterfaceNormals()
{
    davofNormalModel::fillInterfaceNormals(scheme_, mInterface_, filled_);
}


bool Foam::davofNormalModels::ringArea::production() const
{
    return
        stencil_ == "one"
     && kernel_ == "none"
     && !(angleFilter_ > 0)
     && planeFit_ == "none";
}


Foam::scalar Foam::davofNormalModels::ringArea::kernelValue
(
    const scalar q
) const
{
    if (kernel_ == "gaussian")
    {
        return (q < kernelSupport_) ? Foam::exp(-q*q) : 0;
    }
    if (kernel_ == "wendland")
    {
        const scalar r = q/kernelSupport_;
        if (r >= 1) return 0;
        const scalar t = 1 - r;
        return t*t*t*t*(1 + 4*r);
    }
    return 1;
}


void Foam::davofNormalModels::ringArea::ringMembers
(
    labelListList& ringIdx,
    List<vectorField>& ringM
)
{
    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const label nCells = mesh_.nCells();
    zoneDistribute& exchange = zoneDistribute::New(mesh_);

    ringIdx.resize(labels.size());
    ringM.resize(labels.size());
    nRemoteRings_ = 0;

    if (stencil_ == "one")
    {
        const boolList& zone = scheme_.interfaceCell();
        exchange.setUpCommforZone(zone, true);
        const Map<vector> mapM
        (
            exchange.getDatafromOtherProc(zone, mInterface_)
        );
        const labelListList& stencil = exchange.getStencil();

        forAll(labels, i)
        {
            const labelList& st = stencil[labels[i]];
            ringIdx[i] = st;
            vectorField& mv = ringM[i];
            mv.resize(st.size());
            forAll(st, s)
            {
                mv[s] = exchange.getValue(mInterface_, mapM, st[s]);
            }
        }
        return;
    }

    // Two rings. The zone: the interface cells and their point neighbours,
    // marked through the mesh points so that a neighbour on another
    // processor is in that processor's zone (its stencil is computed there)
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();
    labelList pointMark(mesh_.nPoints(), 0);
    for (const label celli : labels)
    {
        for (const label pointi : cellPoints[celli])
        {
            pointMark[pointi] = 1;
        }
    }
    syncTools::syncPointList(mesh_, pointMark, maxEqOp<label>(), label(0));
    boolList zone(nCells, false);
    forAll(pointMark, pointi)
    {
        if (pointMark[pointi])
        {
            for (const label celli : pointCells[pointi])
            {
                zone[celli] = true;
            }
        }
    }

    exchange.setUpCommforZone(zone, true);
    const Map<vector> mapM(exchange.getDatafromOtherProc(zone, mInterface_));
    const labelListList& stencil = exchange.getStencil();
    const globalIndex& gi = exchange.globalNumbering();

    // The rings of the remote ring cells, requested from their owners: the
    // stencil of every requested cell (global indices) and the m values on
    // it, which the owner has in its own exchange map (the cell is in its
    // zone)
    Map<labelList> remoteIdx;
    Map<vectorField> remoteM;

    if (UPstream::parRun())
    {
        const label nProcs = UPstream::nProcs();
        const label myProc = UPstream::myProcNo();

        List<labelHashSet> neededSet(nProcs);
        for (const label celli : labels)
        {
            for (const label gblIdx : stencil[celli])
            {
                if (!gi.isLocal(gblIdx))
                {
                    neededSet[gi.whichProcID(gblIdx)].insert(gblIdx);
                }
            }
        }
        List<labelList> needed(nProcs);
        forAll(neededSet, proci)
        {
            needed[proci] = neededSet[proci].sortedToc();
            nRemoteRings_ += needed[proci].size();
        }

        // The requests
        PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);
        for (const int proci : pBufs.allProcs())
        {
            if (proci != myProc && !needed[proci].empty())
            {
                UOPstream toProc(proci, pBufs);
                toProc << needed[proci];
            }
        }
        pBufs.finishedSends();
        List<labelList> requested(nProcs);
        for (const int proci : pBufs.allProcs())
        {
            if (proci != myProc && pBufs.recvDataCount(proci))
            {
                UIPstream fromProc(proci, pBufs);
                fromProc >> requested[proci];
            }
        }

        // The replies: the stencil and its m values of every requested
        // local cell (a requested boundary face has no stencil: empty)
        PstreamBuffers pBufs2(UPstream::commsTypes::nonBlocking);
        for (const int proci : pBufs2.allProcs())
        {
            const labelList& req = requested[proci];
            if (proci == myProc || req.empty()) continue;

            List<labelList> st(req.size());
            List<vectorField> mv(req.size());
            forAll(req, r)
            {
                const label localIdx = gi.toLocal(req[r]);
                if (localIdx < nCells)
                {
                    st[r] = stencil[localIdx];
                    mv[r].resize(st[r].size());
                    forAll(st[r], s)
                    {
                        mv[r][s] =
                            exchange.getValue(mInterface_, mapM, st[r][s]);
                    }
                }
            }
            UOPstream toProc(proci, pBufs2);
            toProc << st << mv;
        }
        pBufs2.finishedSends();
        for (const int proci : pBufs2.allProcs())
        {
            if (proci != myProc && pBufs2.recvDataCount(proci))
            {
                UIPstream fromProc(proci, pBufs2);
                List<labelList> st;
                List<vectorField> mv;
                fromProc >> st >> mv;
                const labelList& sent = needed[proci];
                if (st.size() != sent.size())
                {
                    FatalErrorInFunction
                        << "ringArea two rings: processor " << proci
                        << " returned " << st.size() << " rings for "
                        << sent.size() << " requests" << abort(FatalError);
                }
                forAll(sent, r)
                {
                    remoteIdx.insert(sent[r], st[r]);
                    remoteM.insert(sent[r], mv[r]);
                }
            }
        }
    }

    // The union over the ring of the rings, every index once
    labelHashSet seen(256);
    DynamicList<label> idx(128);
    DynamicList<vector> mv(128);
    forAll(labels, i)
    {
        const label celli = labels[i];
        seen.clear();
        idx.clear();
        mv.clear();
        for (const label gblIdx : stencil[celli])
        {
            if (gi.isLocal(gblIdx))
            {
                const label localIdx = gi.toLocal(gblIdx);
                if (localIdx >= nCells) continue;   // a boundary face
                for (const label g : stencil[localIdx])
                {
                    if (seen.insert(g))
                    {
                        idx.append(g);
                        mv.append(exchange.getValue(mInterface_, mapM, g));
                    }
                }
            }
            else
            {
                const auto iter = remoteIdx.cfind(gblIdx);
                if (!iter.good()) continue;   // a remote boundary face
                const labelList& st = iter.val();
                const vectorField& sm = remoteM[gblIdx];
                forAll(st, s)
                {
                    if (seen.insert(st[s]))
                    {
                        idx.append(st[s]);
                        mv.append(sm[s]);
                    }
                }
            }
        }
        ringIdx[i] = idx;
        ringM[i] = mv;
    }
    nRemoteRings_ = returnReduce(nRemoteRings_, sumOp<label>());
}


void Foam::davofNormalModels::ringArea::variantNormals(vectorField& n)
{
    addProfilingInFunction(geometricVoF);
    fillInterfaceNormals();

    const DynamicField<label>& labels = scheme_.interfaceLabels();
    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.V();
    const vectorField& m = scheme_.areaNormal().primitiveField();

    const bool needCentres = (kernel_ != "none");
    const bool needPlanes = (planeFit_ != "none");
    if ((needCentres || needPlanes) && stencil_ == "two")
    {
        FatalErrorInFunction
            << "normalModel ringArea: ringKernel and ringPlaneFit are"
            << " implemented on one ring (ringStencil one); the two-ring"
            << " positions and planes are not exchanged" << abort(FatalError);
    }

    labelListList ringIdx;
    List<vectorField> ringM;
    ringMembers(ringIdx, ringM);

    // The exchange after ringMembers() is set up for the interface cells
    // (one ring): the cell centres and, for the plane fit, the positioned
    // planes' centroids and normals of the stencil's cells are read
    // through it
    zoneDistribute& exchange = zoneDistribute::New(mesh_);
    const volVectorField& Cfield = mesh_.C();
    const volVectorField& centre = scheme_.centre();
    const volVectorField& normal = scheme_.normal();
    const boolList& zone = scheme_.interfaceCell();

    Map<vector> mapC, mapX, mapN;
    if (needCentres)
    {
        mapC = exchange.getDatafromOtherProc(zone, Cfield);
    }
    if (needPlanes)
    {
        mapX = exchange.getDatafromOtherProc(zone, centre);
        mapN = exchange.getDatafromOtherProc(zone, normal);
    }

    n.resize(labels.size());
    nFallback_ = 0;
    nAngleDropped_ = 0;
    nAngleKept_ = 0;

    DynamicList<vector> xs(32);
    DynamicList<scalar> ws(32);

    forAll(labels, i)
    {
        const label celli = labels[i];
        const labelList& idx = ringIdx[i];
        const vectorField& mv = ringM[i];
        const scalar h = Foam::cbrt(V[celli]);

        // The (kernel-weighted) ring sum d_k
        vector d(Zero);
        scalarField kw(idx.size(), 1);
        forAll(idx, s)
        {
            const vector& mj = mv[s];
            if (!(mag(mj) > 0)) continue;   // not an interface cell
            if (needCentres)
            {
                const vector xj = exchange.getValue(Cfield, mapC, idx[s]);
                kw[s] = kernelValue(mag(xj - C[celli])/h);
            }
            d += kw[s]*mj;
        }

        // The angle filter on the sum
        if (angleFilter_ > 0 && mag(d) > 0)
        {
            const scalar dMag = mag(d);
            vector d2(Zero);
            forAll(idx, s)
            {
                const vector& mj = mv[s];
                const scalar mjMag = mag(mj);
                if (!(mjMag > 0)) continue;
                const scalar cosAngle = (mj & d)/(mjMag*dMag);
                if (cosAngle < cosAngleFilter_)
                {
                    ++nAngleDropped_;
                    continue;
                }
                d2 += kw[s]*mj;
            }
            if (mag(d2) > 0)
            {
                d = d2;
            }
            else
            {
                ++nAngleKept_;
            }
        }

        n[i] = d;

        if (!needPlanes) continue;

        // The weighted plane through the ring's patch centroids
        xs.clear();
        ws.clear();
        forAll(idx, s)
        {
            const vector& mj = mv[s];
            const scalar mjMag = mag(mj);
            if (!(mjMag > 0)) continue;
            const vector nj = exchange.getValue(normal, mapN, idx[s]);
            if (!(mag(nj) > 0)) continue;   // no positioned plane
            xs.append(exchange.getValue(centre, mapX, idx[s]));
            ws.append(kw[s]*mjMag);
        }

        if (xs.size() < 3)
        {
            ++nFallback_;
            continue;
        }

        // The centre c and the weights: the weighted centroid (mean), or
        // the cell's own patch centroid with the moment-corrected weights
        // (own)
        point c(Zero);
        scalar wSum = 0;
        if (planeFit_ == "mean")
        {
            forAll(xs, s)
            {
                c += ws[s]*xs[s];
                wSum += ws[s];
            }
            if (!(wSum > 0))
            {
                ++nFallback_;
                continue;
            }
            c /= wSum;
        }
        else
        {
            c = (mag(normal[celli]) > 0) ? centre[celli] : C[celli];

            // b = sum w (x - c), M = sum w (x - c)(x - c)^T in units of h;
            // lambda = -M^+ b (the pseudo-inverse over the eigenvalues above
            // tol: in-plane for coplanar centroids)
            vector b(Zero);
            symmTensor M(Zero);
            forAll(xs, s)
            {
                const vector r = (xs[s] - c)/h;
                b += ws[s]*r;
                M += ws[s]*sqr(r);
            }
            const vector mu = eigenValues(M);
            const tensor ev = eigenVectors(M, mu);
            const scalar muMax = max(mag(mu.x()), max(mag(mu.y()), mag(mu.z())));
            vector lambda(Zero);
            for (direction k = 0; k < 3; ++k)
            {
                if (mu[k] > planeFitTol_*muMax)
                {
                    const vector e = ev.row(k);
                    lambda -= ((e & b)/mu[k])*e;
                }
            }
            forAll(xs, s)
            {
                const vector r = (xs[s] - c)/h;
                ws[s] *= (1 + (lambda & r));
            }
        }

        // The weighted covariance about c (in units of h) and its smallest
        // eigenvalue's eigenvector
        symmTensor cov(Zero);
        forAll(xs, s)
        {
            const vector r = (xs[s] - c)/h;
            cov += ws[s]*sqr(r);
        }
        const vector ev = eigenValues(cov);
        const scalar trace = ev.x() + ev.y() + ev.z();
        label kMin = 0;
        for (direction k = 1; k < 3; ++k)
        {
            if (ev[k] < ev[kMin]) kMin = k;
        }
        // The second smallest: the in-plane spread must exist
        scalar second = GREAT;
        for (direction k = 0; k < 3; ++k)
        {
            if (k != label(kMin)) second = min(second, ev[k]);
        }
        if
        (
            !(trace > 0)
         || second <= planeFitTol_*trace
         || ev[kMin] < -planeFitTol_*trace
        )
        {
            ++nFallback_;
            continue;
        }
        const tensor evs = eigenVectors(cov, ev);
        vector nFit = evs.row(kMin);
        const scalar nFitMag = mag(nFit);
        if (!(nFitMag > 0) || nFitMag > ROOTVGREAT)
        {
            ++nFallback_;
            continue;
        }
        nFit /= nFitMag;

        // The sign from the ring sum (or m_k)
        const vector ref = (mag(d) > 0) ? d : m[celli];
        if ((nFit & ref) < 0)
        {
            nFit = -nFit;
        }
        // The magnitude of d_k kept so that the normal carries the ring's
        // area (the scheme normalises)
        n[i] = (mag(d) > 0) ? mag(d)*nFit : nFit;
    }

    nFallback_ = returnReduce(nFallback_, sumOp<label>());
    nAngleDropped_ = returnReduce(nAngleDropped_, sumOp<label>());
    nAngleKept_ = returnReduce(nAngleKept_, sumOp<label>());

    Info<< "davof: ringArea stencil " << stencil_ << ", kernel " << kernel_
        << ", angle " << angleFilter_ << ", planeFit " << planeFit_ << ": "
        << nAngleDropped_ << " dropped, " << nAngleKept_ << " kept d, "
        << nFallback_ << " fit fallbacks, " << nRemoteRings_
        << " remote rings" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::davofNormalModels::ringArea::modelNormals
(
    const scalar,
    vectorField& n
)
{
    addProfilingInFunction(geometricVoF);
    if (production())
    {
        // The one implementation shared with the noNormalFallback ringArea
        // rule of the base class
        ringAreaNormals(scheme_, mInterface_, filled_, n);
        return;
    }
    variantNormals(n);
}


// ************************************************************************* //
