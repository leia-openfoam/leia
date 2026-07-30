/*---------------------------------------------------------------------------*\
    Foot-point HEIGHT-FUNCTION curvature: an unstructured-compatible height
    function built from the LEVEL SET's linear least-squares planes, delivered
    NORMAL-CONSTANT for a balanced CSF force.

    Motivation (root-cause analysis of the parasitic-current instability): a
    cell-centred level-set curvature is profile-invariant, so it always reports
    the curvature of the ISOSURFACE through the cell centre -- parallel
    isosurfaces cannot have equal curvature at different normal offsets
    (kappa_d = kappa/(1 - d kappa)), so the face curvature varies along the
    interface normal across the snGrad(alpha) force support for EVERY interface
    shape: no discrete force-balanced state exists, and the permanent stirring
    pumps grid-scale psi-profile modes that the fit aliases into sign-indefinite
    curvature errors. Height functions equilibrate (Popinet 2009) because the
    column curvature is (i) a functional of interface geometry only, (ii) an
    integral (staircase-immune) estimate, and (iii) constant along the normal
    within the force support. This header reproduces all three properties
    without grid columns:

      1. FEET: each interface cell's linear LSQ plane (the same reconstruction
         the Detrixhe-Aslam indicator uses) is intersected at its zero set:
         p_i = x_i - (n^.x_i + d^) n^. The foot is a geometry-only observable
         of psi (a 2h profile mode of amplitude eps shifts the 5-point plane
         zero-set only ~eps*h/(5*pi)).
      2. LOCAL HEIGHT FUNCTION: the feet of the CPC-neighbour interface cells
         are expressed in the local (tangent, normal) frame and fitted with a
         parabola zeta(tau) = A tau^2/2 + B tau + C by least squares -- a height
         function in the local frame, mesh-topology-free.
         kappa_i = -(A/h_i)/(1+B^2)^{3/2}   [outward-gradient convention,
         kappa = div(grad psi/|grad psi|) = +1/R for the droplet].
      3. NORMAL-CONSTANT DELIVERY: every cell of the band + 2 rings inherits
         the kappa of its NEAREST foot (uniformly, interface cells included) --
         on a circle the nearest-foot assignment is exactly radial, so the two
         cells of every intra-column face carry the same kappa and the face
         force sigma*interp(kappa)*snGrad(alpha) is balanced along the normal;
         the only residual is the tangential fit variation, which the interface
         can relax away (the height-function equilibrium mechanism).

    NOTE this shares only the words "foot point" with the earlier (removed)
    Newton normal extension: there the cell's own contaminated quadratic was
    projected and re-evaluated (amplifying its error 40-85x); here the feet are
    geometric zero-set observables and kappa comes from a cross-cell fit.

    Serial first: parallel needs a feet/kappa band halo (zoneDistribute pattern
    of velocityExtension/closestPoint.C); a warning is emitted in parallel.
\*---------------------------------------------------------------------------*/

#ifndef footPointCurvature_H
#define footPointCurvature_H

#include "fvCFD.H"
#include "HashSet.H"
#include "levelSetPlaneReconstruction.H"

namespace Foam
{

//- Solve the symmetric 3x3 normal-equations system M x = b in place.
//  Returns false if (near-)singular by the same scale-free Hadamard-ratio
//  criterion as levelSetPlaneReconstruction.
inline bool solve3x3(scalar M[3][3], scalar b[3], scalar x[3])
{
    const scalar det =
        M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
      - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
      + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);

    const scalar had =
        Foam::sqrt
        (
            (sqr(M[0][0]) + sqr(M[0][1]) + sqr(M[0][2]))
           *(sqr(M[1][0]) + sqr(M[1][1]) + sqr(M[1][2]))
           *(sqr(M[2][0]) + sqr(M[2][1]) + sqr(M[2][2]))
        );

    if (had < VSMALL || mag(det)/had < 1e-10)
    {
        return false;
    }

    x[0] = (b[0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
          - M[0][1]*(b[1]*M[2][2] - M[1][2]*b[2])
          + M[0][2]*(b[1]*M[2][1] - M[1][1]*b[2]))/det;
    x[1] = (M[0][0]*(b[1]*M[2][2] - M[1][2]*b[2])
          - b[0]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*b[2] - b[1]*M[2][0]))/det;
    x[2] = (M[0][0]*(M[1][1]*b[2] - b[1]*M[2][1])
          - M[0][1]*(M[1][0]*b[2] - b[1]*M[2][0])
          + b[0]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))/det;
    return true;
}


//- Parabola LSQ fit zeta = A tau^2/2 + B tau + C (uniform weights) of the feet
//  (given in the frame (that, nbar) about origin, scaled by h). Returns false
//  if the system is singular. Outputs A, B (and C when Cptr is given).
inline bool fitParabola
(
    const DynamicList<vector>& feet,   // foot positions (absolute)
    const point& origin,
    const vector& that,
    const vector& nbar,
    const scalar h,
    scalar& A,
    scalar& B,
    scalar* Cptr = nullptr
)
{
    scalar M[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    scalar b[3] = {0,0,0};
    forAll(feet, k)
    {
        const vector d = feet[k] - origin;
        const scalar tau = (that & d)/h;
        const scalar zet = (nbar & d)/h;
        const scalar phi0 = 0.5*tau*tau, phi1 = tau, phi2 = 1.0;
        const scalar ph[3] = {phi0, phi1, phi2};
        for (label r = 0; r < 3; ++r)
        {
            b[r] += ph[r]*zet;
            for (label c = 0; c < 3; ++c) { M[r][c] += ph[r]*ph[c]; }
        }
    }
    scalar x[3];
    if (!solve3x3(M, b, x)) { return false; }
    A = x[0]; B = x[1];
    if (Cptr) { *Cptr = x[2]; }
    return true;
}


//- Fill kappa with the foot-point height-function curvature (see file header).
//  narrowBand: the signChange band field (1 = interface-straddling cell).
//  Returns the number of interface cells that kept the fallback (symbolic)
//  curvature because their fit was deficient.
inline label computeFootPointCurvature
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand,
    volScalarField& kappa
)
{
    if (Pstream::parRun())
    {
        WarningInFunction
            << "footPointHeightFunction curvature is serial-only in this "
            << "version (feet are not exchanged across processor patches); "
            << "results near processor boundaries are inconsistent." << endl;
    }

    const vectorField& C = mesh.cellCentres();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    // Empty direction (2D): geometricD() component == -1.
    vector ez(Zero);
    {
        const Vector<label> gd = mesh.geometricD();
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            if (gd[cmpt] == -1) { ez[cmpt] = 1; }
        }
        if (mag(ez) < SMALL)
        {
            FatalErrorInFunction
                << "footPointHeightFunction: 3D not implemented (needs a "
                << "paraboloid fit); use it on 2D (one empty direction) meshes."
                << abort(FatalError);
        }
    }

    // h_c: min face-neighbour centre distance (empty-direction-safe because
    // the empty direction has no internal faces).
    // Collected only for band cells (all we need).
    const coupledFaceNeighbours coupledNei(mesh, psi);

    // ---- Step 1: interface cells -> accepted feet -------------------------
    DynamicList<label>  fCell;     // owning interface cell
    DynamicList<vector> fFoot;     // foot position
    DynamicList<vector> fNormal;   // unit plane normal (outward grad psi)
    labelList footIndex(mesh.nCells(), -1);   // cell -> index in the lists

    scalarField hCell(mesh.nCells(), GREAT);
    forAll(own, facei)
    {
        const scalar d = mag(C[own[facei]] - C[nei[facei]]);
        hCell[own[facei]] = min(hCell[own[facei]], d);
        hCell[nei[facei]] = min(hCell[nei[facei]], d);
    }

    forAll(narrowBand, celli)
    {
        if (narrowBand[celli] != 1) { continue; }

        const scalarList pc =
            leastSquaresPlaneCoeffs(mesh, psi, celli, coupledNei);
        const vector n(pc[0], pc[1], pc[2]);
        const scalar nmag = mag(n);
        if (nmag < SMALL) { continue; }              // degenerate plane

        const vector nHat = n/nmag;
        const scalar dHat = pc[3]/nmag;
        const point  p = C[celli] - ((nHat & C[celli]) + dHat)*nHat;

        // Leverage guard: a sign-change cell's zero crossing must be near the
        // cell; a farther foot means a bad plane -- never let it into fits.
        if (mag(p - C[celli]) > 1.2*hCell[celli]) { continue; }

        footIndex[celli] = fCell.size();
        fCell.append(celli);
        fFoot.append(p);
        fNormal.append(nHat);
    }

    if (fCell.empty())
    {
        return 0;    // no interface: nothing to do (kappa keeps the fallback)
    }

    // CPC neighbourhood helper (cell -> cells sharing a point).
    const labelListList& cellPoints = mesh.cellPoints();
    const labelListList& pointCells = mesh.pointCells();

    auto gatherFeet = [&]
    (
        const label celli,
        const vector& nHatI,
        const label rings,               // 1 = CPC, 2 = CPC^2
        DynamicList<vector>& feet,
        DynamicList<vector>& normals
    )
    {
        feet.clear(); normals.clear();
        labelHashSet visited;
        labelHashSet current;
        current.insert(celli);
        for (label r = 0; r < rings; ++r)
        {
            labelHashSet next;
            for (const label c : current)
            {
                for (const label pt : cellPoints[c])
                {
                    for (const label cc : pointCells[pt])
                    {
                        if (!visited.found(cc)) { next.insert(cc); }
                    }
                }
                visited.insert(c);
            }
            for (const label cc : next) { visited.insert(cc); }
            current = next;
        }
        for (const label cc : visited)
        {
            const label fi = footIndex[cc];
            if (fi < 0) { continue; }
            // same-branch filter: reject feet whose normal disagrees > 45 deg
            if ((fNormal[fi] & nHatI) < 0.7071) { continue; }
            feet.append(fFoot[fi]);
            normals.append(fNormal[fi]);
        }
    };

    // ---- Steps 2-3: per-interface-cell parabola fit -----------------------
    scalarField kappaFoot(fCell.size(), 0);
    boolList    kappaOk(fCell.size(), false);
    label nFallback = 0;

    DynamicList<vector> feet(16), normals(16);
    forAll(fCell, fi)
    {
        const label celli = fCell[fi];
        const scalar h = hCell[celli];
        const vector& nHatI = fNormal[fi];
        const point&  pI = fFoot[fi];

        bool ok = false;
        scalar A = 0, B = 0;

        for (label rings = 1; rings <= 2 && !ok; ++rings)
        {
            gatherFeet(celli, nHatI, rings, feet, normals);

            if (feet.size() < 5) { continue; }

            // averaged frame normal (same-branch feet)
            vector nbar(Zero);
            forAll(normals, k) { nbar += normals[k]; }
            const scalar nbmag = mag(nbar);
            if (nbmag < SMALL) { continue; }
            nbar /= nbmag;
            vector that = ez ^ nbar;
            that /= mag(that);

            // sufficiency: two-sided coverage + tangential extent >= 2h
            scalar tmin = GREAT, tmax = -GREAT;
            label nneg = 0, npos = 0;
            forAll(feet, k)
            {
                const scalar tau = (that & (feet[k] - pI))/h;
                tmin = min(tmin, tau); tmax = max(tmax, tau);
                if (tau < -SMALL) { ++nneg; }
                else if (tau > SMALL) { ++npos; }
            }
            if (nneg < 2 || npos < 2 || (tmax - tmin) < 2.0) { continue; }

            if (!fitParabola(feet, pI, that, nbar, h, A, B)) { continue; }

            // strongly tilted frame: refit once in the rotated frame
            if (mag(B) > 0.7)
            {
                const scalar s = 1.0/Foam::sqrt(1.0 + B*B);
                const vector that2 = (that + B*nbar)*s;
                const vector nbar2 = (nbar - B*that)*s;
                scalar A2 = 0, B2 = 0;
                if (fitParabola(feet, pI, that2, nbar2, h, A2, B2))
                {
                    A = A2; B = B2;
                }
            }
            ok = true;
        }

        if (ok)
        {
            // Sign: zeta measured along +nbar ~ outward grad psi; for the
            // droplet the interface curves AWAY from outward (A < 0) and the
            // repo convention kappa = div(grad psi/|grad psi|) = +1/R.
            kappaFoot[fi] = -(A/hCell[celli])/Foam::pow(1.0 + B*B, 1.5);
            kappaOk[fi] = true;
        }
        else
        {
            ++nFallback;   // cell keeps the fallback (symbolic) kappa
        }
    }

    // ---- Step 4: normal-constant delivery (uniform nearest accepted foot) --
    // delivery set = band + 2 face-connectivity rings
    boolList inDelivery(mesh.nCells(), false);
    forAll(narrowBand, celli)
    {
        if (narrowBand[celli] == 1) { inDelivery[celli] = true; }
    }
    for (label sweep = 0; sweep < 2; ++sweep)
    {
        boolList grow(inDelivery);
        forAll(own, facei)
        {
            if (inDelivery[own[facei]] || inDelivery[nei[facei]])
            {
                grow[own[facei]] = true;
                grow[nei[facei]] = true;
            }
        }
        inDelivery = grow;
    }

    // accepted (fit-ok) feet only
    DynamicList<vector> okFoot(fCell.size());
    DynamicList<vector> okNormal(fCell.size());
    DynamicList<scalar> okKappa(fCell.size());
    DynamicList<scalar> okH(fCell.size());
    forAll(fCell, fi)
    {
        if (kappaOk[fi])
        {
            okFoot.append(fFoot[fi]);
            okNormal.append(fNormal[fi]);
            okKappa.append(kappaFoot[fi]);
            okH.append(hCell[fCell[fi]]);
        }
    }

    if (okFoot.size())
    {
        // COLUMN (ray-distance) parent assignment. Plain Euclidean nearest-foot
        // is provably radial on the CIRCLE, but under interface DEFORMATION its
        // Voronoi walls tilt and can cut through a normal column: two cells on
        // one normal ray get different parents -> a kappa jump INSIDE the force
        // support -> dipole force -> more deformation (measured: coupled runs
        // blow up in ~1 ms while the frozen-interface run is quiet). The column
        // criterion assigns each cell to the foot whose own NORMAL RAY passes
        // closest to the cell centre (minimal tangential ray distance),
        // preserving column integrity under deformation; the small normal-
        // distance term breaks ties smoothly and keeps parents local.
        scalarField& k = kappa.primitiveFieldRef();
        forAll(inDelivery, celli)
        {
            if (!inDelivery[celli]) { continue; }
            scalar smin = GREAT;
            label  jmin = -1;
            forAll(okFoot, j)
            {
                const vector d = C[celli] - okFoot[j];
                const scalar d2 = magSqr(d);
                if (d2 > sqr(6.0*okH[j])) { continue; }   // locality prefilter
                const scalar dn = d & okNormal[j];
                const scalar dray2 = d2 - dn*dn;          // tangential ray dist^2
                const scalar score = dray2 + 0.01*dn*dn;  // smooth tie-break
                if (score < smin) { smin = score; jmin = j; }
            }
            if (jmin < 0)
            {
                // fallback: plain nearest (should not trigger in the band)
                scalar dmin = GREAT;
                forAll(okFoot, j)
                {
                    const scalar d2 = magSqr(C[celli] - okFoot[j]);
                    if (d2 < dmin) { dmin = d2; jmin = j; }
                }
            }
            if (jmin >= 0) { k[celli] = okKappa[jmin]; }
        }
        kappa.correctBoundaryConditions();
    }

    Info<< "footPointHeightFunction curvature: " << fCell.size()
        << " feet (" << okFoot.size() << " fits ok, " << nFallback
        << " fallback); kappa in ["
        << gMin(kappa.primitiveField()) << ", "
        << gMax(kappa.primitiveField()) << "] 1/m" << endl;

    return nFallback;
}

} // End namespace Foam

#endif

// ************************************************************************* //
