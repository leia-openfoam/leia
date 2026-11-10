/*---------------------------------------------------------------------------*\
    leia :: davof :: isoAdvectorScheme (see isoAdvectorScheme.H)
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "isoAdvectorScheme.H"
#include "davof.H"
#include "cutFacePLIC.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoAdvectorScheme, 0);
    addToRunTimeSelectionTable
    (
        davofAdvectionScheme,
        isoAdvectorScheme,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoAdvectorScheme::isoAdvectorScheme
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    davofAdvectionScheme(typeName, alpha1, phi, U),
    iso_(alpha1, phi, U),
    davofPtr_(nullptr)
{
    Info<< "isoAdvector: OpenFOAM's isoAdvection (libgeometricVoF) with the"
        << " reconstruction scheme " << iso_.surf().type() << endl;
    if (isA<reconstruction::davof>(iso_.surf()))
    {
        davofPtr_ = &refCast<reconstruction::davof>(iso_.surf());
        Info<< "isoAdvector: the DAVOF state's face fractions "
            << davofPtr_->alphaf().name()
            << " are set after every step from the translated planes of the"
            << " step's reconstruction (the upwind cell's plane; the"
            << " face-fraction corrector " << davofPtr_->faceFractionCorrector()
            << ")" << endl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::isoAdvectorScheme::predictFaceFractions()
{
    const fvMesh& mesh = alpha1_.mesh();
    reconstruction::davof& dv = *davofPtr_;
    surfaceScalarField& af = dv.alphaf();
    // the planes and the interface cells of the step's reconstruction
    // (isoAdvection reconstructs at the start of advect; nothing has
    // reconstructed since)
    const vectorField& nrm = dv.normal().primitiveField();
    const vectorField& ctr = dv.centre().primitiveField();
    const boolList& ic = dv.interfaceCell();
    const vectorField& Uc = U_.primitiveField();
    const scalarField& aIn = alpha1_.primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const vectorField& Sf = mesh.faceAreas();
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();
    const label nBnd = mesh.nBoundaryFaces();
    const scalar dt = mesh.time().deltaTValue();
    if (ic.size() != mesh.nCells()) return;

    // the clip of a face by a plane translated by u dt (cutFacePLIC,
    // H. Scheufler, J. Roenby, J. Comput. Phys. 383 (2019) 1-23,
    // doi:10.1016/j.jcp.2019.01.009; the normal points into the liquid)
    cutFacePLIC cutF(mesh);
    // the UNIT normal (calcSubFace snaps |(x - base).n| < SMALL to zero,
    // davof.C inducedFaceFractions, 2026-11-10)
    auto clip = [&](const label facei, const vector& n, const point& c, const vector& u) -> scalar
    {
        cutF.calcSubFace(facei, n/mag(n), c + u*dt);
        return min(max(mag(cutF.subFaceArea())/mag(Sf[facei]), scalar(0)), scalar(1));
    };
    // the value of a face from the two cells (u the upwind one): u's
    // translated plane, else d's, else u's new alpha
    auto value = [&]
    (
        const label facei,
        const bool icU, const vector& nU, const point& cU, const vector& uU, const scalar aU,
        const bool icD, const vector& nD, const point& cD, const vector& uD
    ) -> scalar
    {
        if (icU && mag(nU) > 0) return clip(facei, nU, cU, uU);
        if (icD && mag(nD) > 0) return clip(facei, nD, cD, uD);
        return min(max(aU, scalar(0)), scalar(1));
    };

    label nPlane = 0, nAlpha = 0;
    scalarField& afIn = af.primitiveFieldRef();
    for (label facei = 0; facei < nInt; ++facei)
    {
        const label o = own[facei], n = nei[facei];
        if (!(ic[o] || ic[n])) continue;
        const bool oUp = (phiIn[facei] >= 0);
        const label u = oUp ? o : n, d = oUp ? n : o;
        afIn[facei] = value(facei, ic[u], nrm[u], ctr[u], Uc[u], aIn[u], ic[d], nrm[d], ctr[d], Uc[d]);
        if ((ic[u] && mag(nrm[u]) > 0) || (ic[d] && mag(nrm[d]) > 0)) ++nPlane; else ++nAlpha;
    }

    // the boundary faces: the remote cell's plane, velocity, alpha and flag
    // on the coupled patches by a swap of the owner-side values
    vectorField nB(nBnd), cB(nBnd), uB(nBnd);
    scalarField aB(nBnd), fB(nBnd);
    for (label bfi = 0; bfi < nBnd; ++bfi)
    {
        const label o = own[nInt + bfi];
        nB[bfi] = nrm[o]; cB[bfi] = ctr[o]; uB[bfi] = Uc[o]; aB[bfi] = aIn[o];
        fB[bfi] = ic[o] ? 1 : 0;
    }
    syncTools::swapBoundaryFaceList(mesh, nB);
    syncTools::swapBoundaryFaceList(mesh, cB);
    syncTools::swapBoundaryFaceList(mesh, uB);
    syncTools::swapBoundaryFaceList(mesh, aB);
    syncTools::swapBoundaryFaceList(mesh, fB);

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    scalarField both(nBnd, Zero);
    bitSet coupledSet(nBnd);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        fvsPatchScalarField& pf = af.boundaryFieldRef()[patchi];
        if (pf.size() != pp.size()) continue;
        const fvsPatchScalarField& pphi = phi_.boundaryField()[patchi];
        forAll(pp, i)
        {
            const label facei = pp.start() + i;
            const label bfi = facei - nInt;
            const label o = own[facei];
            if (pp.coupled())
            {
                const bool icR = (fB[bfi] > 0.5);
                if (!(ic[o] || icR)) continue;
                const bool oUp = (pphi[i] >= 0);
                const scalar v = oUp
                  ? value(facei, ic[o], nrm[o], ctr[o], Uc[o], aIn[o], icR, nB[bfi], cB[bfi], uB[bfi])
                  : value(facei, icR, nB[bfi], cB[bfi], uB[bfi], aB[bfi], ic[o], nrm[o], ctr[o], Uc[o]);
                both[bfi] = v;
                coupledSet.set(bfi);
            }
            else
            {
                if (!ic[o]) continue;
                pf[i] = (mag(nrm[o]) > 0)
                  ? clip(facei, nrm[o], ctr[o], Uc[o])
                  : min(max(aIn[o], scalar(0)), scalar(1));
            }
        }
    }
    if (Pstream::parRun() || coupledSet.any())
    {
        // one value per coupled face: the half-sum of the two sides (the
        // same plane clipped on both sides differs at round-off only)
        syncTools::syncBoundaryFaceList(mesh, both, plusEqOp<scalar>());
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            if (!pp.coupled()) continue;
            fvsPatchScalarField& pf = af.boundaryFieldRef()[patchi];
            if (pf.size() != pp.size()) continue;
            forAll(pp, i)
            {
                const label bfi = pp.start() + i - nInt;
                if (coupledSet.test(bfi)) pf[i] = 0.5*both[bfi];
            }
        }
    }
    Info<< "isoAdvector: face fractions from the translated planes on "
        << returnReduce(nPlane, sumOp<label>()) << " faces, from alpha on "
        << returnReduce(nAlpha, sumOp<label>()) << " faces of the interface cells"
        << endl;
}


void Foam::isoAdvectorScheme::afterAdvect()
{
    if (!davofPtr_) return;
    predictFaceFractions();
    davofPtr_->correctFaceFractions();
}


// ************************************************************************* //
