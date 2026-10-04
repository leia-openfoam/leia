/*---------------------------------------------------------------------------*\
    leia :: davof :: zero-step regeneration

    One application of the reconstruct-clip-recover loop that the DAVOF
    transport closes, without a time step: from a static state
    (alpha_k, {alpha_f}) each interface cell recovers its normal (the Gauss
    identity), its explicit plane and, from the ring of recovered normals,
    its curvature tensor; a clipping surface is built from these three and
    every face of the cell is cut again with it. The regenerated face
    fractions are compared with a reference AS FRACTIONS, and the normal and
    the position are recovered from them once more.

    Arms (the clipping surface of a cell, pre-registered predictions in
    config/davof/ellipsoidRegeneration3D.yaml):

        plane            the explicit plane (2)             q = 1 predicted
        paraboloidPlane  the paraboloid with the ring tensor, its vertex on
                         the plane polygon's centroid       q = 1 predicted
        paraboloidFaces  the same paraboloid, its offset along the normal
                         fixed by the cell's own face fractions (the total
                         wetted face area of the cell is reproduced; one
                         monotone scalar equation)          q = 2 predicted
        paraboloidVolume the same paraboloid, its offset fixed by the cell
                         fraction (the comparison arm; on the chord cell
                         fraction of the quadraticFaces state it lands on
                         the chord)                        q = 1 predicted

    Nothing here matches a volume on purpose except the comparison arm: the
    method's anchoring is the face fractions (2026-10-03).

    The regeneration is cell-local: a shared face is cut once per adjacent
    interface cell with that cell's surface; the two values are the cell's
    own views, their disagreement is reported as the C0 diagnostic. Serial
    studies only (the neighbour's surface across a processor boundary is not
    exchanged).
\*---------------------------------------------------------------------------*/

#ifndef davofRegeneration_H
#define davofRegeneration_H

#include "fvMesh.H"
#include "davofSimplexGeometry.H"
#include "davofQuadraticFaceGeometry.H"
#include "levelSetImplicitSurfaces.H"
#include <functional>

namespace Foam
{
namespace davof
{

//- The clipping surface of one cell: the plane {w = c} or the paraboloid
//  w = c + 1/2 (b11 u^2 + 2 b12 u v + b22 v^2) in the frame (o; e1, e2, e3),
//  u = (x - o).e1, v = (x - o).e2, w = (x - o).e3, e3 the DAVOF normal
//  pointing out of the liquid. value(x) < 0 in the liquid.
struct cellSurface
{
    bool paraboloid = false;
    point o = Zero;
    vector e1 = Zero, e2 = Zero, e3 = Zero;
    scalar b11 = 0, b12 = 0, b22 = 0, c = 0;

    scalar value(const point& x) const
    {
        const vector r = x - o;
        const scalar w = r & e3;
        if (!paraboloid) return w - c;
        const scalar u = r & e1, v = r & e2;
        return w - c - 0.5*(b11*u*u + 2*b12*u*v + b22*v*v);
    }
};


//- Liquid area of the triangle (a, b, c) under phi: the quadratic cut
//  (six values: vertices and edge midpoints; exact for a paraboloid
//  restricted to an edge); where it falls back (more than two boundary
//  crossings, an arc that is no graph over its chord, a vertex exactly on the
//  surface) the triangle is split into four and the cut retried, down to
//  `levels` more splits, then the linear cut (counted). The recursion keeps
//  the cut, and with it the anchoring equation of the paraboloid, continuous
//  in the offset.
template<class PhiAt>
inline scalar triLiquidArea
(
    const point& a, const point& b, const point& c,
    PhiAt&& phi,
    const label levels,
    label& nFallback
)
{
    const scalar ph[6] =
    {
        phi(a), phi(b), phi(c),
        phi(0.5*(a + b)), phi(0.5*(b + c)), phi(0.5*(c + a))
    };
    const scalar A = 0.5*mag((b - a) ^ (c - a));
    if (A <= VSMALL) return 0;
    bool fb = false;
    const scalar fr = triQuadraticNegativeFraction(a, b, c, ph, fb);
    if (!fb) return fr*A;
    if (levels <= 0)
    {
        ++nFallback;
        return triNegativeFraction(ph[0], ph[1], ph[2])*A;
    }
    const point mab = 0.5*(a + b), mbc = 0.5*(b + c), mca = 0.5*(c + a);
    return
        triLiquidArea(a, mab, mca, phi, levels - 1, nFallback)
      + triLiquidArea(mab, b, mbc, phi, levels - 1, nFallback)
      + triLiquidArea(mca, mbc, c, phi, levels - 1, nFallback)
      + triLiquidArea(mab, mbc, mca, phi, levels - 1, nFallback);
}


//- The cut of face faceI by phi(x): the fan about x_f, triLiquidArea on
//  every triangle. Returns the liquid area fraction; wet = sum fr S_t
//  (oriented like the face), M = sum x_wet,t . (fr S_t) with the wet
//  centroid of the linear cut (exact on a planar face).
template<class PhiAt>
inline scalar faceCut
(
    const fvMesh& mesh,
    const label faceI,
    PhiAt&& phi,
    vector& wet,
    scalar& M,
    label& nFallback
)
{
    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const point& xf = mesh.faceCentres()[faceI];
    const scalar phif = phi(xf);
    scalar liqA = 0, totA = 0;
    wet = Zero;
    M = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point& p0 = P[f[ip]];
        const point& p1 = P[ip1];
        const vector St = 0.5*((p0 - xf) ^ (p1 - xf));
        const scalar triA = mag(St);
        if (triA <= VSMALL) continue;
        const scalar fr = triLiquidArea(xf, p0, p1, phi, 3, nFallback)/triA;
        liqA += fr*triA;
        totA += triA;
        wet += fr*St;
        vector Sl;
        point xwet;
        triNegativePolygon(xf, p0, p1, phif, phi(p0), phi(p1), Sl, xwet);
        M += xwet & (fr*St);
    }
    return (totA > VSMALL) ? liqA/totA : (phif < 0 ? scalar(1) : scalar(0));
}


//- sum_f |S_f| alpha_f^phi over the faces of cellI (the anchoring quantity of
//  the paraboloidFaces arm).
template<class PhiAt>
inline scalar cellWettedArea(const fvMesh& mesh, const label cellI, PhiAt&& phi)
{
    const labelList& cFaces = mesh.cells()[cellI];
    const vectorField& Sf = mesh.faceAreas();
    scalar a = 0;
    label nfb = 0;
    vector wet;
    scalar M;
    forAll(cFaces, i)
    {
        const label f = cFaces[i];
        a += faceCut(mesh, f, phi, wet, M, nfb)*mag(Sf[f]);
    }
    return a;
}


//- The liquid volume of cellI under phi on the Detrixhe-Aslam tets
//  (x_c, x_f, p_i, p_{i+1}); linear in each tet, second order for a curved
//  surface, which is enough for the comparison arm paraboloidVolume.
template<class PhiAt>
inline scalar cellLiquidVolume(const fvMesh& mesh, const label cellI, PhiAt&& phi)
{
    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& Cf = mesh.faceCentres();
    const point& xc = mesh.cellCentres()[cellI];
    const labelList& cFaces = mesh.cells()[cellI];
    const scalar phic = phi(xc);
    scalar liquidVol = 0;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        const scalar phif = phi(xf);
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            const point& p0 = P[f[ip]];
            const point& p1 = P[ip1];
            const scalar vol = mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
            if (vol <= VSMALL) continue;
            liquidVol += tetNegativeFraction(phic, phif, phi(p0), phi(p1))*vol;
        }
    }
    return liquidVol;
}


//- The liquid volume of cellI under phi with every Detrixhe-Aslam tet split
//  `levels` times into eight sub-tets (the linear interpolant on each): the
//  chord defect of the tet quadrature shrinks by 4^levels, so this is the
//  paraboloid's own cut volume to the order the comparison needs.
template<class PhiAt>
inline scalar cellLiquidVolumeRefined
(
    const fvMesh& mesh,
    const label cellI,
    PhiAt&& phi,
    const label levels
)
{
    std::function<scalar(const point&, const point&, const point&, const point&, label)> rec =
        [&](const point& a, const point& b, const point& c, const point& d, label lev) -> scalar
        {
            const scalar vol = mag(((b - a) ^ (c - a)) & (d - a))/6.0;
            if (vol <= VSMALL) return 0;
            if (lev == 0)
            {
                return tetNegativeFraction(phi(a), phi(b), phi(c), phi(d))*vol;
            }
            // The eight sub-tets of the midpoint subdivision (four corner
            // tets and the central octahedron split along one diagonal).
            const point ab = 0.5*(a + b), ac = 0.5*(a + c), ad = 0.5*(a + d);
            const point bc = 0.5*(b + c), bd = 0.5*(b + d), cd = 0.5*(c + d);
            return
                rec(a, ab, ac, ad, lev - 1) + rec(ab, b, bc, bd, lev - 1)
              + rec(ac, bc, c, cd, lev - 1) + rec(ad, bd, cd, d, lev - 1)
              + rec(ab, ac, ad, bd, lev - 1) + rec(ab, ac, bc, bd, lev - 1)
              + rec(ac, ad, bd, cd, lev - 1) + rec(ac, bc, bd, cd, lev - 1);
        };

    const pointField& P = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& Cf = mesh.faceCentres();
    const point& xc = mesh.cellCentres()[cellI];
    const labelList& cFaces = mesh.cells()[cellI];
    scalar liquidVol = 0;
    forAll(cFaces, cf)
    {
        const label faceI = cFaces[cf];
        const face& f = faces[faceI];
        const point& xf = Cf[faceI];
        forAll(f, ip)
        {
            const label ip1 = f.nextLabel(ip);
            liquidVol += rec(xc, xf, P[f[ip]], P[ip1], levels);
        }
    }
    return liquidVol;
}


//- Solve g(c) = target for the offset c by bisection; g is increasing in c
//  (a larger c puts more of the cell into the liquid). The bracket
//  [lo, hi] is widened up to eight times if it does not bracket the root.
//  Returns |g(c) - target|; bracketed is false when no bracket was found
//  (c is then the midpoint of the last bracket).
template<class G>
inline scalar solveOffset
(
    G&& g,
    const scalar target,
    scalar lo,
    scalar hi,
    const scalar tol,
    scalar& c,
    bool& bracketed
)
{
    scalar glo = g(lo) - target, ghi = g(hi) - target;
    for (label k = 0; k < 8 && glo > 0; ++k)
    {
        lo -= (hi - lo);
        glo = g(lo) - target;
    }
    for (label k = 0; k < 8 && ghi < 0; ++k)
    {
        hi += (hi - lo);
        ghi = g(hi) - target;
    }
    bracketed = !(glo > 0 || ghi < 0);
    if (!bracketed)
    {
        c = 0.5*(lo + hi);
        return min(mag(glo), mag(ghi));
    }
    scalar mid = 0.5*(lo + hi), gm = 0;
    for (label it = 0; it < 400; ++it)
    {
        mid = 0.5*(lo + hi);
        gm = g(mid) - target;
        if (gm > 0) hi = mid; else lo = mid;
        if (hi - lo <= tol) break;
    }
    c = 0.5*(lo + hi);
    return mag(g(c) - target);
}


//- Reference face fraction of faceI for an exact surface s: the quadratic
//  cut on every fan triangle subdivided `levels` times (4^levels
//  sub-triangles; the reference error is 4^levels times below the
//  unsubdivided quadratic cut, which is itself O(h^2) as a fraction).
inline scalar refFaceFraction
(
    const fvMesh& mesh,
    const label faceI,
    const implicitSurface& s,
    const label levels,
    label& nFallback
)
{
    std::function<scalar(const point&, const point&, const point&, label)> rec =
        [&](const point& a, const point& b, const point& c, label lev) -> scalar
        {
            if (lev == 0)
            {
                const scalar ph[6] =
                {
                    s.value(a), s.value(b), s.value(c),
                    s.value(0.5*(a + b)), s.value(0.5*(b + c)), s.value(0.5*(c + a))
                };
                bool fb = false;
                scalar fr = triQuadraticNegativeFraction(a, b, c, ph, fb);
                if (fb)
                {
                    fr = triNegativeFraction(ph[0], ph[1], ph[2]);
                    ++nFallback;
                }
                return fr*0.5*mag((b - a) ^ (c - a));
            }
            const point mab = 0.5*(a + b), mbc = 0.5*(b + c), mca = 0.5*(c + a);
            return
                rec(a, mab, mca, lev - 1) + rec(mab, b, mbc, lev - 1)
              + rec(mca, mbc, c, lev - 1) + rec(mab, mbc, mca, lev - 1);
        };

    const pointField& P = mesh.points();
    const face& f = mesh.faces()[faceI];
    const point& xf = mesh.faceCentres()[faceI];
    scalar liqA = 0, totA = 0;
    forAll(f, ip)
    {
        const label ip1 = f.nextLabel(ip);
        const point& p0 = P[f[ip]];
        const point& p1 = P[ip1];
        const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
        if (triA <= VSMALL) continue;
        liqA += rec(xf, p0, p1, levels);
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (s.value(xf) < 0 ? scalar(1) : scalar(0));
}


//- Re-cut the faces of every cell with isI[c] with its surface surf[c]:
//  the cell's own view of its face fractions (alphafOwn for the faces it
//  owns, alphafNei for the faces it neighbours), its fan sums mTri (outward
//  area vector, the sign convention of davofState) and momentTri. Cells
//  with !isI keep the values passed in.
inline void regenerate
(
    const fvMesh& mesh,
    const boolList& isI,
    const List<cellSurface>& surf,
    scalarField& alphafOwn,
    scalarField& alphafNei,
    vectorField& mTri,
    scalarField& momentTri,
    label& nFallback
)
{
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const label nInt = mesh.nInternalFaces();
    const vectorField& C = mesh.cellCentres();
    forAll(isI, c)
    {
        if (isI[c])
        {
            mTri[c] = Zero;
            momentTri[c] = 0;
        }
    }
    forAll(own, faceI)
    {
        const label o = own[faceI];
        if (isI[o])
        {
            const cellSurface& S = surf[o];
            vector wet;
            scalar M;
            alphafOwn[faceI] = faceCut
            (
                mesh, faceI, [&](const point& x) { return S.value(x); },
                wet, M, nFallback
            );
            mTri[o] -= wet;
            momentTri[o] += M - (C[o] & wet);
        }
        if (faceI < nInt)
        {
            const label n = nei[faceI];
            if (isI[n])
            {
                const cellSurface& S = surf[n];
                vector wet;
                scalar M;
                alphafNei[faceI] = faceCut
                (
                    mesh, faceI, [&](const point& x) { return S.value(x); },
                    wet, M, nFallback
                );
                mTri[n] += wet;
                momentTri[n] -= M - (C[n] & wet);
            }
        }
    }
}

} // namespace davof
} // namespace Foam

#endif

// ************************************************************************* //
