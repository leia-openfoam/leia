/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Description
    ISO-VALUE-AGNOSTIC GEOMETRIC MEAN CURVATURE, shared by the isoCurvature surface-
    tension model and the leiaTestMeanCurvature accuracy harness. 3D algorithm that
    also runs in OpenFOAM pseudo-2D (one empty direction).

    Interface from the iso-value crossings of ANY field (only sign(field-isoValue) is
    used): the crossed internal faces are detected from the VERTEX (point-interpolated)
    field values -- a face is crossed if field-isoValue changes sign along one of its
    edges -- giving a clean interface cut.

    Two cell-centred curvature reductions (curvatureMethod):

    (0) CONORMAL / integral surface-tension identity (the default):
          (kappa n)_c A_c = sum_{f in cut faces of cell c} c_f,   c_f = n_f x L_f,
        where L_f is the ACTUAL interface-in-face segment vector (between the two iso-
        crossing points on the face's edges) and n_f the interface normal at the segment
        MIDPOINT (stabilised owner<->neighbour MEAN of the cell-centred grad-normal). The
        single per-face contribution c_f = n_f x L_f is ORIENTED to point OUT of the
        owner: sign so c_f . S_f >= 0 for the owner-oriented mesh face vector S_f.
        The face normal is not part of the conormal magnitude; it supplies the exact
        topological orientation, including corner cuts where x_mid-C_cell is ambiguous.
        The contribution is added to the owner with its outward sign and to the neighbour
        with the OPPOSITE sign (the neighbour's outward is reversed), so ONE contribution
        per face feeds both face-adjacent cells' curvature-normal sums.
        kappa_c = -(sum_f c_f).n_c / A_c,
        A_c = L_int h_z (L_int = in-plane chord). This is well conditioned: for a corner
        cut the two adjacent faces give nearly the SAME n_f x L_f, which sums tangentially
        and projects to ~0 on n_c (bounded), rather than blowing up.

    (1) PARABOLA / local height function (comparison): fit h(s)=a s^2+b s+c to the
        crossing points of the cell + its one-ring interface neighbours in the local
        tangent-normal frame; kappa = -2a/(1+b^2)^{3/2}. Local, shape-agnostic, but a
        fixed-stencil point fit has an O(1) noise floor (does not converge).

    Optionally EXTEND kappa nExtend layers and SMOOTH nSmooth Brackbill passes. Fills
    kappa in the interface band (0 elsewhere), sets `filled` and `nHole`, returns nGood.
    sgn aligns the sign with +1/R for a drop (psi outward normal -> +1).

\*---------------------------------------------------------------------------*/

#ifndef isoGeometricCurvature_H
#define isoGeometricCurvature_H

#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvcGrad.H"
#include "volPointInterpolation.H"
#include "HashSet.H"
#include "boundBox.H"

namespace Foam
{

//- Iso-agnostic geometric mean curvature into kappa. curvatureMethod: 0 = conormal
//  integral (sum_f L_f x n_f, orientation from the outward face normal), 1 = local
//  tangent-frame parabola fit. Returns nGood.
inline label computeIsoGeometricCurvature
(
    const fvMesh& mesh,
    const volScalarField& field,   // geometry field (psi, alpha, ...)
    const scalar isoValue,         // iso-value that IS the interface
    const scalar sgn,              // overall sign (+1 outward normal, -1 inward)
    const label nExtend,           // curvature band-extension layers
    const label nSmooth,           // Brackbill smoothing passes
    volScalarField& kappa,         // out: curvature field
    boolList& filled,              // out: cells carrying a valid curvature
    label& nHole,                  // out: fit-failed interface cells
    const label curvatureMethod = 0  // 0 = conormal integral, 1 = parabola fit
)
{
    // 1. Cell interface unit normal nHat = grad(field)/|grad(field)|, SMALL-guarded.
    volVectorField nHat("nHatIso", fvc::grad(field));
    {
        const dimensionedScalar deltaN("deltaN", nHat.dimensions(), SMALL);
        nHat /= (mag(nHat) + deltaN);
    }
    const vectorField& nHatI = nHat.primitiveField();

    // 2. Empty-direction unit vector e_z (used only to take the in-plane chord length).
    const Vector<label> gd = mesh.geometricD();
    vector ez(Zero);
    label nEmpty = 0;
    for (direction c = 0; c < 3; ++c) { if (gd[c] == -1) { ez[c] = 1; ++nEmpty; } }
    if (nEmpty != 1)
    {
        FatalErrorInFunction
            << "isoGeometricCurvature is a 2D routine and needs exactly one empty "
            << "direction; geometricD() reports " << nEmpty << exit(FatalError);
    }
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();                    // owner-outward face vectors

    // Point-interpolated field for robust vertex edge-crossing detection.
    const pointScalarField fieldP
    (
        volPointInterpolation::New(mesh).interpolate(field)
    );
    const scalarField& fp    = fieldP.primitiveField();
    const faceList&    faces  = mesh.faces();
    const pointField&  pts    = mesh.points();

    // 3. Single face sweep. For each crossed internal face build the interface-in-face
    //    SEGMENT L_f (between its two iso-crossing edge points) and its midpoint x_mid;
    //    n_f = stabilised owner<->neighbour MEAN normal; the single per-face contribution
    //    c_f = n_f x L_f is oriented to point OUT of the owner (co-normal) and added to
    //    the owner (+) and neighbour (-) curvature-normal sums knS. The midpoint is the
    //    per-cell iso-crossing point (xa,xb) for the interface chord / parabola stencil.
    const label ncell = mesh.nCells();
    vectorField xa(ncell, Zero), xb(ncell, Zero);   // the two iso-crossing midpoints
    vectorField knS(ncell, Zero);                   // sum of outward co-normal contribs
    labelField  ncross(ncell, 0);

    forAll(own, f)
    {
        // Collect the (two) iso-crossing points on this face's edges -> segment.
        const face& fc = faces[f];
        point q0(Zero), q1(Zero);
        label nq = 0;
        forAll(fc, ip)
        {
            const label i0 = fc[ip];
            const label i1 = fc[fc.fcIndex(ip)];
            const scalar g0 = fp[i0] - isoValue;
            const scalar g1 = fp[i1] - isoValue;
            if (g0*g1 < 0)
            {
                const scalar lam = g0/(g0 - g1);
                const point q = pts[i0] + lam*(pts[i1] - pts[i0]);
                if (nq == 0) { q0 = q; } else if (nq == 1) { q1 = q; }
                ++nq;
            }
        }
        if (nq != 2) { continue; }               // need a proper interface segment

        const label o = own[f], n = nei[f];
        const vector Lf = q1 - q0;               // interface-in-face segment vector
        const point  xm = 0.5*(q0 + q1);         // segment midpoint

        // Stabilised owner<->neighbour MEAN interface normal at the midpoint.
        vector nf = nHatI[o] + nHatI[n];
        nf /= (mag(nf) + SMALL);

        // Single per-face contribution, oriented to point OUT of the OWNER cell.
        // Use the actual owner-oriented face vector, not (crossing-cellCentre):
        // the latter is nearly orthogonal to the conormal at corner cuts and its
        // sign can flip under roundoff, making the two endpoint tractions add.
        vector cf = nf ^ Lf;                     // n_f x L_f  (a co-normal x |L_f|)
        if ((cf & Sf[f]) < 0) { cf = -cf; }      // make it point out of the owner
        knS[o] += cf;
        knS[n] -= cf;                            // neighbour outward is reversed

        if      (ncross[o] == 0) { xa[o] = xm; }
        else if (ncross[o] == 1) { xb[o] = xm; }
        ++ncross[o];
        if      (ncross[n] == 0) { xa[n] = xm; }
        else if (ncross[n] == 1) { xb[n] = xm; }
        ++ncross[n];
    }

    // In-plane axes (for the parabola method).
    label ax0 = -1, ax1 = -1;
    for (direction d = 0; d < 3; ++d)
    {
        if (gd[d] == 1) { if (ax0 < 0) { ax0 = d; } else { ax1 = d; } }
    }
    const labelListList& cc = mesh.cellCells();

    // Empty-direction thickness h_z, for the interface area A_c = L_int h_z.
    const boundBox bb(mesh.points());
    const scalar hz = max(cmptSum(cmptMultiply(bb.span(), ez)), SMALL);
    const scalarField& V = mesh.V();   // for the in-plane cell size h_ip = sqrt(V/h_z)

    scalarField& k = kappa.primitiveFieldRef();
    k = 0.0;
    filled.setSize(ncell); filled = false;
    label nGood = 0; nHole = 0;

    forAll(k, c)
    {
        if (ncross[c] == 0) { continue; }
        if (ncross[c] != 2) { ++nHole; continue; }
        // In-plane interface chord length L_int (drop the empty-direction component).
        const vector d = xb[c] - xa[c];
        const scalar ds = mag(d - (d & ez)*ez);
        if (ds < SMALL) { ++nHole; continue; }

        if (curvatureMethod == 0)
        {
            // CONORMAL: (kappa n)_c = (sum_f c_f)/A_c, A_c = L_int h_z; scalar curvature
            // kappa_c = -(sum_f c_f).n_c / A_c. The extra minus makes the outward-conormal
            // convention give +1/R for a convex drop (psi outward normal), so sgn=+1
            // matches the parabola. Through cuts (interface across OPPOSITE faces) are
            // exact; extracting a POINTWISE kappa needs /A_c, which is ill-posed for a
            // CORNER cut (adjacent faces, L_int -> 0 while the curvature-normal does not).
            // Skip those short-chord cells and fill them by extension from the exact
            // through-cells (the integrated curvature-normal/force sum_f c_f itself is
            // bounded regardless -- only the pointwise-kappa division is ill-conditioned).
            const scalar hip = Foam::sqrt(V[c]/hz);           // in-plane cell size
            if (ds < 0.5*hip) { ++nHole; continue; }          // ill-posed corner cut
            k[c] = -sgn*(knS[c] & nHatI[c])/(ds*hz);
            filled[c] = true;
            ++nGood;
        }
        else
        {
            // LOCAL PARABOLA height-function fit over the one-ring neighbourhood.
            labelHashSet sten;
            sten.insert(c);
            forAll(cc[c], i) { if (ncross[cc[c][i]] == 2) { sten.insert(cc[c][i]); } }

            const scalar nx = nHatI[c][ax0], ny = nHatI[c][ax1];
            const scalar tx = -ny, ty = nx;
            const point  xm = 0.5*(xa[c] + xb[c]);
            const scalar x0 = xm[ax0], y0 = xm[ax1];

            scalar sumr2 = 0; label np = 0;
            for (const label g : sten)
            {
                for (label e = 0; e < 2; ++e)
                {
                    const point& P = (e == 0) ? xa[g] : xb[g];
                    const scalar dx = P[ax0]-x0, dy = P[ax1]-y0;
                    sumr2 += dx*dx + dy*dy; ++np;
                }
            }
            if (np < 3) { ++nHole; continue; }
            const scalar L2 = sumr2/np;
            if (L2 < SMALL) { ++nHole; continue; }
            const scalar L = Foam::sqrt(L2);

            scalar S0=0,S1=0,S2=0,S3=0,S4=0,T0=0,T1=0,T2=0;
            for (const label g : sten)
            {
                for (label e = 0; e < 2; ++e)
                {
                    const point& P = (e == 0) ? xa[g] : xb[g];
                    const scalar dx = (P[ax0]-x0)/L, dy = (P[ax1]-y0)/L;
                    const scalar s  = dx*tx + dy*ty;
                    const scalar hh = dx*nx + dy*ny;
                    const scalar s2 = s*s;
                    S0 += 1; S1 += s; S2 += s2; S3 += s2*s; S4 += s2*s2;
                    T0 += hh; T1 += hh*s; T2 += hh*s2;
                }
            }
            const scalar det =
                S4*(S2*S0 - S1*S1) - S3*(S3*S0 - S1*S2) + S2*(S3*S1 - S2*S2);
            if (mag(det) < SMALL) { ++nHole; continue; }
            const scalar detA =
                T2*(S2*S0 - S1*S1) - S3*(T1*S0 - S1*T0) + S2*(T1*S1 - S2*T0);
            const scalar detB =
                S4*(T1*S0 - T0*S1) - T2*(S3*S0 - S1*S2) + S2*(S3*T0 - T1*S2);
            const scalar a = detA/det;
            const scalar b = detB/det;
            k[c] = sgn*(-2.0*a/L)/Foam::pow(1.0 + b*b, 1.5);
            filled[c] = true;
            ++nGood;
        }
    }

    // Band extension.
    for (label e = 0; e < nExtend; ++e)
    {
        scalarField ksum(ncell, 0.0), kw(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && !filled[n]) { ksum[n] += k[o]; kw[n] += 1.0; }
            if (filled[n] && !filled[o]) { ksum[o] += k[n]; kw[o] += 1.0; }
        }
        boolList newly(ncell, false);
        forAll(k, c)
        {
            if (!filled[c] && kw[c] > 0.5) { k[c] = ksum[c]/kw[c]; newly[c] = true; }
        }
        forAll(newly, c) { if (newly[c]) { filled[c] = true; } }
    }

    // Brackbill smoothing.
    for (label s = 0; s < nSmooth; ++s)
    {
        scalarField ksum(ncell, 0.0), kw(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && filled[n])
            {
                ksum[o] += k[n]; kw[o] += 1.0;
                ksum[n] += k[o]; kw[n] += 1.0;
            }
        }
        forAll(k, c)
        {
            if (filled[c] && kw[c] > 0.5) { k[c] = 0.5*k[c] + 0.5*(ksum[c]/kw[c]); }
        }
    }

    return nGood;
}


//- DIRECT conormal surface-tension body-force DENSITY (the well-conditioned use of the
//  integral formulation). Per crossed face the single outward-oriented contribution
//  c_f = (n_f x L_f)_outward is summed to the owner (+) and neighbour (-); the cell body
//  force per unit sigma is g_c = -sgn * (sum_f c_f)/V_c  [1/m^2] (surface-tension force
//  per volume = sigma*g_c [N/m^3]). This divides by the CELL VOLUME V_c (never zero), NOT
//  the interface area A_c (-> 0 at corner cuts), so unlike the pointwise curvature it is
//  well conditioned for every cut geometry. Optionally EXTEND nExtend layers (so the
//  face interpolation of the band force is not halved at the band edge) and SMOOTH
//  nSmooth Brackbill passes. Returns the number of interface cells filled.
inline label computeIsoConormalForceDensity
(
    const fvMesh& mesh,
    const volScalarField& field,
    const scalar isoValue,
    const scalar sgn,
    const label nExtend,
    const label nSmooth,
    volVectorField& g,             // out: -sgn*(sum_f c_f)/V_c ; body force = sigma*g
    boolList* filledOut = nullptr  // out (optional): cells carrying a force value
)
{
    volVectorField nHat("nHatIso", fvc::grad(field));
    {
        const dimensionedScalar deltaN("deltaN", nHat.dimensions(), SMALL);
        nHat /= (mag(nHat) + deltaN);
    }
    const vectorField& nHatI = nHat.primitiveField();

    const Vector<label> gd = mesh.geometricD();
    label nEmpty = 0;
    for (direction c = 0; c < 3; ++c) { if (gd[c] == -1) { ++nEmpty; } }
    if (nEmpty != 1)
    {
        FatalErrorInFunction
            << "isoGeometricCurvature is a 2D routine and needs exactly one empty "
            << "direction; geometricD() reports " << nEmpty << exit(FatalError);
    }
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();
    const scalarField& V  = mesh.V();

    const pointScalarField fieldP
    (
        volPointInterpolation::New(mesh).interpolate(field)
    );
    const scalarField& fp    = fieldP.primitiveField();
    const faceList&    faces  = mesh.faces();
    const pointField&  pts    = mesh.points();

    const label ncell = mesh.nCells();
    vectorField knS(ncell, Zero);   // sum of outward co-normal contributions per cell
    labelField  ncross(ncell, 0);

    forAll(own, f)
    {
        const face& fc = faces[f];
        point q0(Zero), q1(Zero);
        label nq = 0;
        forAll(fc, ip)
        {
            const label i0 = fc[ip];
            const label i1 = fc[fc.fcIndex(ip)];
            const scalar g0 = fp[i0] - isoValue;
            const scalar g1 = fp[i1] - isoValue;
            if (g0*g1 < 0)
            {
                const scalar lam = g0/(g0 - g1);
                const point q = pts[i0] + lam*(pts[i1] - pts[i0]);
                if (nq == 0) { q0 = q; } else if (nq == 1) { q1 = q; }
                ++nq;
            }
        }
        if (nq != 2) { continue; }

        const label o = own[f], n = nei[f];
        const vector Lf = q1 - q0;
        vector nf = nHatI[o] + nHatI[n];
        nf /= (mag(nf) + SMALL);
        vector cf = nf ^ Lf;
        if ((cf & Sf[f]) < 0) { cf = -cf; }          // out of the owner
        knS[o] += cf;
        knS[n] -= cf;
        ++ncross[o]; ++ncross[n];
    }

    vectorField& gI = g.primitiveFieldRef();
    gI = Zero;
    boolList filled(ncell, false);
    label nGood = 0;
    forAll(gI, c)
    {
        if (ncross[c] == 2)
        {
            // Keep only the NORMAL component of the conormal sum. Surface tension on a
            // smooth interface is purely normal (sigma*kappa*n); the TANGENTIAL part of
            // sum_f c_f is a spurious discrete artefact that is large exactly at corner
            // cuts (where the two adjacent-face contributions add tangentially) and would
            // drive parasitic currents ALONG the interface. Projecting on n_c gives the
            // correct through-cell force AND ~0 at corner cuts (their sum is tangential),
            // so it is well conditioned without any /A_c blow-up.
            const vector nc = nHatI[c];
            gI[c] = -sgn*(knS[c] & nc)*nc/V[c];
            filled[c] = true;
            ++nGood;
        }
    }

    // Vector band extension (avoid halving the band force under face interpolation).
    for (label e = 0; e < nExtend; ++e)
    {
        vectorField gsum(ncell, Zero); scalarField w(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && !filled[n]) { gsum[n] += gI[o]; w[n] += 1.0; }
            if (filled[n] && !filled[o]) { gsum[o] += gI[n]; w[o] += 1.0; }
        }
        boolList newly(ncell, false);
        forAll(gI, c)
        {
            if (!filled[c] && w[c] > 0.5) { gI[c] = gsum[c]/w[c]; newly[c] = true; }
        }
        forAll(newly, c) { if (newly[c]) { filled[c] = true; } }
    }

    // Vector Brackbill smoothing.
    for (label s = 0; s < nSmooth; ++s)
    {
        vectorField gsum(ncell, Zero); scalarField w(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && filled[n])
            {
                gsum[o] += gI[n]; w[o] += 1.0;
                gsum[n] += gI[o]; w[n] += 1.0;
            }
        }
        forAll(gI, c)
        {
            if (filled[c] && w[c] > 0.5) { gI[c] = 0.5*gI[c] + 0.5*(gsum[c]/w[c]); }
        }
    }

    g.correctBoundaryConditions();
    if (filledOut) { *filledOut = filled; }
    return nGood;
}


//- Curvature-free, pairwise-conservative integral surface-tension traction.
//  For every internal mesh face cut by the reconstructed interface, build the
//  interface-in-face segment L_f and the outward-owner conormal
//
//      c_f = (n_Gamma x L_f)_outward-owner .
//
//  The oriented face value t_f = -sgn*sigma*c_f [N] is stored ONCE. Therefore
//  fvc::surfaceIntegrate(t_f) adds +t_f/V_owner and -t_f/V_neighbour and gives
//  the cell-integrated surface-stress divergence without an explicit curvature
//  field, cell-force interpolation, or independently deposited owner/neighbour
//  sources. Closed interfaces away from boundaries have roundoff-zero net force.
//
//  Prototype scope matches the geometric routines above: serial, pseudo-2D,
//  internal faces. Boundary/coupled-face intersections are deliberately rejected
//  by omission until a single shared parallel reconstruction is implemented.
inline label computeIsoConormalFaceTraction
(
    const fvMesh& mesh,
    const volScalarField& field,
    const scalar isoValue,
    const scalar sgn,
    const scalar sigma,
    surfaceVectorField& traction
)
{
    volVectorField nHat("nHatIsoTraction", fvc::grad(field));
    {
        const dimensionedScalar deltaN("deltaN", nHat.dimensions(), SMALL);
        nHat /= (mag(nHat) + deltaN);
    }
    const vectorField& nHatI = nHat.primitiveField();

    const Vector<label> gd = mesh.geometricD();
    label nEmpty = 0;
    for (direction c = 0; c < 3; ++c)
    {
        if (gd[c] == -1) { ++nEmpty; }
    }
    if (nEmpty != 1)
    {
        FatalErrorInFunction
            << "computeIsoConormalFaceTraction is a 2D prototype and needs "
            << "exactly one empty direction; geometricD() reports " << nEmpty
            << exit(FatalError);
    }

    const pointScalarField fieldP
    (
        volPointInterpolation::New(mesh).interpolate(field)
    );
    const scalarField& fp = fieldP.primitiveField();
    const faceList& faces = mesh.faces();
    const pointField& pts = mesh.points();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    vectorField& tI = traction.primitiveFieldRef();
    tI = Zero;
    label nCut = 0;

    forAll(own, f)
    {
        const face& fc = faces[f];
        point q0(Zero), q1(Zero);
        label nq = 0;
        forAll(fc, ip)
        {
            const label i0 = fc[ip];
            const label i1 = fc[fc.fcIndex(ip)];
            const scalar g0 = fp[i0] - isoValue;
            const scalar g1 = fp[i1] - isoValue;
            if (g0*g1 < 0)
            {
                const scalar lambda = g0/(g0 - g1);
                const point q = pts[i0] + lambda*(pts[i1] - pts[i0]);
                if (nq == 0) { q0 = q; }
                else if (nq == 1) { q1 = q; }
                ++nq;
            }
        }
        if (nq != 2) { continue; }

        const label o = own[f];
        const label n = nei[f];
        const vector Lf = q1 - q0;
        vector nf = nHatI[o] + nHatI[n];
        nf /= (mag(nf) + SMALL);

        vector cf = nf ^ Lf;
        if ((cf & Sf[f]) < 0) { cf = -cf; }

        // surfaceIntegrate uses the owner-oriented value for the owner and its
        // negative for the neighbour, exactly matching the two outward conormals.
        tI[f] = -sgn*sigma*cf;
        ++nCut;
    }

    traction.correctBoundaryConditions();
    return nCut;
}


//- Scalar band extension + Brackbill smoothing (shared by the conormal-CSF curvature).
inline void extendSmoothCellScalar
(
    const fvMesh& mesh,
    boolList& filled,
    const label nExtend,
    const label nSmooth,
    scalarField& k
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const label ncell = mesh.nCells();

    for (label e = 0; e < nExtend; ++e)
    {
        scalarField ksum(ncell, 0.0), w(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && !filled[n]) { ksum[n] += k[o]; w[n] += 1.0; }
            if (filled[n] && !filled[o]) { ksum[o] += k[n]; w[o] += 1.0; }
        }
        boolList newly(ncell, false);
        forAll(k, c)
        {
            if (!filled[c] && w[c] > 0.5) { k[c] = ksum[c]/w[c]; newly[c] = true; }
        }
        forAll(newly, c) { if (newly[c]) { filled[c] = true; } }
    }

    for (label s = 0; s < nSmooth; ++s)
    {
        scalarField ksum(ncell, 0.0), w(ncell, 0.0);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (filled[o] && filled[n])
            {
                ksum[o] += k[n]; w[o] += 1.0;
                ksum[n] += k[o]; w[n] += 1.0;
            }
        }
        forAll(k, c)
        {
            if (filled[c] && w[c] > 0.5) { k[c] = 0.5*k[c] + 0.5*(ksum[c]/w[c]); }
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
