/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCurvatureDelivery (davofCurvatureDelivery.H)
\*---------------------------------------------------------------------------*/

#include "davofCurvatureDelivery.H"
#include "syncTools.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofCurvatureDelivery::davofCurvatureDelivery
(
    const fvMesh& mesh,
    const davofState& state,
    const davofCurvature& model
)
:
    mesh_(mesh),
    state_(state),
    model_(model),
    kappaCell_(mesh.nCells(), 0),
    KCell_(mesh.nCells(), 0),
    bandCell_(mesh.nCells(), false),
    filledCell_(mesh.nCells(), false),
    kappaFaceInv_(mesh.nFaces(), 0),
    kappaFaceFoot_(mesh.nFaces(), 0),
    dFace_(mesh.nFaces(), 0),
    KFace_(mesh.nFaces(), 0),
    activeFace_(mesh.nFaces(), false),
    oneSidedFace_(mesh.nFaces(), false),
    nBandCells_(0),
    nUnfilledCells_(0),
    nActiveFaces_(0),
    nOneSidedFaces_(0),
    nSkippedFaces_(0)
{
    const label nInt = mesh.nInternalFaces();
    const label nBnd = mesh.nBoundaryFaces();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = mesh.faceCentres();
    const scalarField& alpha = state.alpha().primitiveField();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const scalar aTol = 1e-12;

    // ---- Across coupled patches: the neighbour's alpha and model ---------------
    scalarField alphaNb(nBnd, 0);
    List<davofModelPack> packNb(nBnd, davofModelPack(Zero));
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        forAll(pp, i)
        {
            const label bf = pp.start() - nInt + i;
            const label c = pp.faceCells()[i];
            alphaNb[bf] = alpha[c];
            model.pack(c, packNb[bf]);
        }
    }
    syncTools::swapBoundaryFaceList(mesh, alphaNb);
    syncTools::swapBoundaryFaceList(mesh, packNb);

    // ---- Active faces and the force band ---------------------------------------
    // A face is active when alpha jumps across it AND a model exists on at
    // least one side: the faces of wisps (a jump of 1e-5 between a wisp and
    // the bulk, no model anywhere near) carry no deliverable curvature and are
    // excluded like the wisps themselves; they are counted (nSkippedFaces).
    boolList hasModel(mesh.nCells(), false);
    {
        davofModelPack tmp;
        forAll(hasModel, c) hasModel[c] = model.pack(c, tmp);
    }
    for (label f = 0; f < nInt; ++f)
    {
        if (mag(alpha[own[f]] - alpha[nei[f]]) > aTol)
        {
            if (!hasModel[own[f]] && !hasModel[nei[f]])
            {
                ++nSkippedFaces_;
                continue;
            }
            activeFace_[f] = true;
            bandCell_[own[f]] = true;
            bandCell_[nei[f]] = true;
        }
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled()) continue;
        const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(pp);
        forAll(pp, i)
        {
            const label bf = pp.start() - nInt + i;
            const label c = pp.faceCells()[i];
            if (mag(alpha[c] - alphaNb[bf]) > aTol)
            {
                if (!hasModel[c] && !(packNb[bf][0] > 0.5))
                {
                    if (cpp.owner()) ++nSkippedFaces_;
                    continue;
                }
                bandCell_[c] = true;
                // Scored once: on the owner side of the coupled pair.
                if (cpp.owner()) activeFace_[pp.start() + i] = true;
            }
        }
    }

    // ---- The cell field on the band --------------------------------------------
    const cellList& cells = mesh.cells();
    forAll(bandCell_, c)
    {
        if (!bandCell_[c]) continue;
        ++nBandCells_;
        scalar d, kappa, K;
        if (model.offset(c, C[c], d, kappa, K))
        {
            kappaCell_[c] = parallelSurfaceForward(kappa, d, K);
            KCell_[c] = parallelSurfaceGaussianForward(kappa, d, K);
            filledCell_[c] = true;
            continue;
        }
        // No own model: the face-adjacent models, Kang-weighted by 1/|d|.
        const scalar hc = Foam::cbrt(mesh.V()[c]);
        scalar sw = 0, sk = 0, sK = 0;
        const cell& cf = cells[c];
        forAll(cf, i)
        {
            const label f = cf[i];
            bool ok = false;
            if (f < nInt)
            {
                const label other = (own[f] == c) ? nei[f] : own[f];
                ok = model.offset(other, C[c], d, kappa, K);
            }
            else
            {
                const label patchi = pbm.whichPatch(f);
                if (patchi < 0 || !pbm[patchi].coupled()) continue;
                ok = model.offsetPacked(packNb[f - nInt], C[c], d, kappa, K);
            }
            if (!ok) continue;
            const scalar w = 1.0/max(mag(d), 1e-12*hc);
            sw += w;
            sk += w*parallelSurfaceForward(kappa, d, K);
            sK += w*parallelSurfaceGaussianForward(kappa, d, K);
        }
        if (sw > 0)
        {
            kappaCell_[c] = sk/sw;
            KCell_[c] = sK/sw;
            filledCell_[c] = true;
        }
        else
        {
            ++nUnfilledCells_;
        }
    }

    // The neighbour's cell value across coupled patches
    scalarField kappaCellNb(nBnd, 0);
    scalarField filledNb(nBnd, 0);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        forAll(pp, i)
        {
            const label bf = pp.start() - nInt + i;
            const label c = pp.faceCells()[i];
            kappaCellNb[bf] = kappaCell_[c];
            filledNb[bf] = filledCell_[c] ? 1 : 0;
        }
    }
    syncTools::swapBoundaryFaceList(mesh, kappaCellNb);
    syncTools::swapBoundaryFaceList(mesh, filledNb);

    // ---- The face fields on the active faces -----------------------------------
    const surfaceScalarField& wts = mesh.weights();

    auto deliver = [&]
    (
        const label f,
        const scalar wOwn,
        const bool okO, const scalar dO, const scalar kO, const scalar KO,
        const bool okN, const scalar dN, const scalar kN, const scalar KN,
        const bool filledO, const scalar kappaCellO,
        const bool filledN, const scalar kappaCellN
    )
    {
        ++nActiveFaces_;
        scalar dF = 0, KF = 0, kFoot = 0;
        if (okO && okN)
        {
            scalar wO = mag(dN), wN = mag(dO);
            if (wO + wN < VSMALL) { wO = wN = 0.5; }
            const scalar s = wO + wN;
            dF = (wO*dO + wN*dN)/s;
            KF =
            (
                wO*parallelSurfaceGaussianForward(kO, dO, KO)
              + wN*parallelSurfaceGaussianForward(kN, dN, KN)
            )/s;
            kFoot = (wO*kO + wN*kN)/s;
        }
        else if (okO || okN)
        {
            oneSidedFace_[f] = true;
            ++nOneSidedFaces_;
            const scalar dd = okO ? dO : dN;
            const scalar kk = okO ? kO : kN;
            const scalar KK = okO ? KO : KN;
            dF = dd;
            KF = parallelSurfaceGaussianForward(kk, dd, KK);
            kFoot = kk;
        }
        else
        {
            // No model on either side: nothing to deliver (left at zero).
            return;
        }
        scalar kInterp;
        if (filledO && filledN)
        {
            kInterp = wOwn*kappaCellO + (1 - wOwn)*kappaCellN;
        }
        else if (filledO || filledN)
        {
            kInterp = filledO ? kappaCellO : kappaCellN;
        }
        else
        {
            return;
        }
        dFace_[f] = dF;
        KFace_[f] = KF;
        kappaFaceInv_[f] = parallelSurfaceInverse(kInterp, dF, KF);
        kappaFaceFoot_[f] = kFoot;
    };

    for (label f = 0; f < nInt; ++f)
    {
        if (!activeFace_[f]) continue;
        scalar dO = 0, kO = 0, KO = 0, dN = 0, kN = 0, KN = 0;
        const bool okO = model.offset(own[f], Cf[f], dO, kO, KO);
        const bool okN = model.offset(nei[f], Cf[f], dN, kN, KN);
        deliver
        (
            f, wts[f],
            okO, dO, kO, KO, okN, dN, kN, KN,
            filledCell_[own[f]], kappaCell_[own[f]],
            filledCell_[nei[f]], kappaCell_[nei[f]]
        );
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!pp.coupled()) continue;
        const fvsPatchScalarField& wp = wts.boundaryField()[patchi];
        forAll(pp, i)
        {
            const label f = pp.start() + i;
            if (!activeFace_[f]) continue;
            const label bf = f - nInt;
            const label c = pp.faceCells()[i];
            scalar dO = 0, kO = 0, KO = 0, dN = 0, kN = 0, KN = 0;
            const bool okO = model.offset(c, Cf[f], dO, kO, KO);
            const bool okN = model.offsetPacked(packNb[bf], Cf[f], dN, kN, KN);
            deliver
            (
                f, wp[i],
                okO, dO, kO, KO, okN, dN, kN, KN,
                filledCell_[c], kappaCell_[c],
                filledNb[bf] > 0.5, kappaCellNb[bf]
            );
        }
    }

    reduce(nBandCells_, sumOp<label>());
    reduce(nUnfilledCells_, sumOp<label>());
    reduce(nActiveFaces_, sumOp<label>());
    reduce(nOneSidedFaces_, sumOp<label>());
    reduce(nSkippedFaces_, sumOp<label>());
}


// ************************************************************************* //
