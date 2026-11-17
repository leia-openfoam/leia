/*---------------------------------------------------------------------------*\
    leia :: pointIso :: plane fit of pointIsoAdvection
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

Description
    The plane fitted to the displaced points of a PLIC polygon in
    pointIsoAdvection (linear least squares, eigenvectors by Jacobi
    rotations; OpenFOAM's eigenVectors returned a wrong vector on a
    near-degenerate covariance in this project and is not used).

    - No boundary point: the plane through the mean c' of the points, the
      normal the eigenvector of the smallest eigenvalue of
      sum (q_i - c')(q_i - c')^T.
    - One boundary point b: the plane through b, the normal the smallest
      eigenvector of sum (q_i - b)(q_i - b)^T.
    - Two boundary points (or several collinear ones: the two farthest
      apart): the plane contains the line through them; the normal is
      perpendicular to the line and minimises the sum of squared distances in
      the two-dimensional subspace perpendicular to it.
    - Three non-collinear boundary points: the plane through them.
    - More than three non-collinear boundary points: counted (overThree) and
      the line through the two farthest apart is used.
    The orientation is n'.nk > 0. Collinear means the largest distance of a
    boundary point from the line through the two farthest apart is at most
    1e-12 of their distance (the round-off scale, no model parameter).

    References:
    \verbatim
        J. Roenby, H. Bredmose, H. Jasak,
        A computational method for sharp interface advection,
        Royal Society Open Science 3 (2016) 160405,
        doi:10.1098/rsos.160405, arXiv:1601.05392
        H. Scheufler, J. Roenby,
        Accurate and efficient surface reconstruction from volume fraction
        data on general meshes,
        Journal of Computational Physics 383 (2019) 1-23,
        doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef pointIsoPlaneFit_H
#define pointIsoPlaneFit_H

#include "vector.H"
#include "DynamicList.H"
#include <cmath>

namespace Foam
{
namespace pointIso
{

//- Cyclic Jacobi rotations of the symmetric N x N matrix A: on return
//  ev holds the eigenvalues and the columns of V the eigenvectors
template<int N>
inline void jacobiEigen(scalar (&A)[N][N], scalar (&V)[N][N], scalar (&ev)[N])
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < 100; ++sweep)
    {
        scalar off = 0;
        for (int p = 0; p < N - 1; ++p)
        {
            for (int q = p + 1; q < N; ++q)
            {
                off += A[p][q]*A[p][q];
            }
        }
        if (off == 0) break;

        for (int p = 0; p < N - 1; ++p)
        {
            for (int q = p + 1; q < N; ++q)
            {
                if (A[p][q] == 0) continue;

                const scalar theta = (A[q][q] - A[p][p])/(2*A[p][q]);
                const scalar t =
                    (theta >= 0 ? 1.0 : -1.0)
                   /(std::fabs(theta) + std::hypot(theta, 1.0));
                const scalar c = 1.0/std::sqrt(t*t + 1.0);
                const scalar s = t*c;
                const scalar Apq = A[p][q];

                A[p][p] -= t*Apq;
                A[q][q] += t*Apq;
                A[p][q] = 0;
                A[q][p] = 0;

                for (int r = 0; r < N; ++r)
                {
                    if (r == p || r == q) continue;
                    const scalar Arp = A[r][p];
                    const scalar Arq = A[r][q];
                    A[r][p] = c*Arp - s*Arq;
                    A[p][r] = A[r][p];
                    A[r][q] = c*Arq + s*Arp;
                    A[q][r] = A[r][q];
                }
                for (int r = 0; r < N; ++r)
                {
                    const scalar Vrp = V[r][p];
                    const scalar Vrq = V[r][q];
                    V[r][p] = c*Vrp - s*Vrq;
                    V[r][q] = s*Vrp + c*Vrq;
                }
            }
        }
    }

    for (int i = 0; i < N; ++i)
    {
        ev[i] = A[i][i];
    }
}


//- The result of a plane fit
struct fit
{
    vector n;           //- unit normal, n.nk > 0
    point x;            //- point of the plane
    scalar rms;         //- rms distance of the points to the plane
    label nBoundary;    //- number of boundary points
    bool overThree;     //- more than three non-collinear boundary points
};


//- The smallest eigenvector of sum (q_i - b)(q_i - b)^T
inline vector smallestAbout
(
    const UList<point>& q,
    const point& b
)
{
    scalar A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    scalar V[3][3];
    scalar ev[3];

    forAll(q, i)
    {
        const vector d(q[i] - b);
        for (int a = 0; a < 3; ++a)
        {
            for (int c = 0; c < 3; ++c)
            {
                A[a][c] += d[a]*d[c];
            }
        }
    }

    jacobiEigen<3>(A, V, ev);

    int k = 0;
    for (int i = 1; i < 3; ++i)
    {
        if (ev[i] < ev[k]) k = i;
    }
    return vector(V[0][k], V[1][k], V[2][k]);
}


//- Fit the plane to the points q; isB marks the boundary points; nk the
//  orientation (n'.nk > 0)
inline fit fitPlane
(
    const UList<point>& q,
    const UList<bool>& isB,
    const vector& nk
)
{
    fit res;
    res.overThree = false;

    DynamicList<label> bi;
    forAll(q, i)
    {
        if (isB[i]) bi.append(i);
    }
    const label nb = bi.size();
    res.nBoundary = nb;

    vector n(Zero);
    point x0(Zero);

    // The two boundary points farthest apart
    label ia = -1;
    label ib = -1;
    scalar dmax = 0;
    for (label i = 0; i < nb; ++i)
    {
        for (label j = i + 1; j < nb; ++j)
        {
            const scalar dd = magSqr(q[bi[i]] - q[bi[j]]);
            if (dd > dmax)
            {
                dmax = dd;
                ia = bi[i];
                ib = bi[j];
            }
        }
    }

    if (nb == 0)
    {
        point c(Zero);
        forAll(q, i) c += q[i];
        c /= scalar(q.size());
        n = smallestAbout(q, c);
        x0 = c;
    }
    else if (nb == 1 || dmax == 0)
    {
        x0 = q[bi[0]];
        n = smallestAbout(q, x0);
    }
    else
    {
        const point b1 = q[ia];
        vector d(q[ib] - b1);
        const scalar len = mag(d);
        d /= len;

        // The boundary point farthest from the line
        scalar maxDist = 0;
        label ic = -1;
        for (const label k : bi)
        {
            if (k == ia || k == ib) continue;
            const vector r(q[k] - b1);
            const scalar dist = mag(r - (r & d)*d);
            if (dist > maxDist)
            {
                maxDist = dist;
                ic = k;
            }
        }
        const bool collinear = (maxDist <= 1e-12*len);

        if (nb == 3 && !collinear)
        {
            n = (q[ib] - b1) ^ (q[ic] - b1);
            n /= mag(n);
            x0 = b1;
        }
        else
        {
            if (nb > 3 && !collinear) res.overThree = true;

            // A basis of the plane perpendicular to the line
            vector axis(1, 0, 0);
            if (mag(d.y()) < mag(d[0]) && mag(d.y()) <= mag(d.z()))
            {
                axis = vector(0, 1, 0);
            }
            else if (mag(d.z()) < mag(d[0]) && mag(d.z()) < mag(d.y()))
            {
                axis = vector(0, 0, 1);
            }
            vector e1(d ^ axis);
            e1 /= mag(e1);
            const vector e2(d ^ e1);

            scalar A[2][2] = {{0, 0}, {0, 0}};
            scalar V[2][2];
            scalar ev[2];
            forAll(q, i)
            {
                const vector r(q[i] - b1);
                const scalar u = r & e1;
                const scalar v = r & e2;
                A[0][0] += u*u;
                A[0][1] += u*v;
                A[1][0] += u*v;
                A[1][1] += v*v;
            }
            jacobiEigen<2>(A, V, ev);
            const int k = (ev[1] < ev[0]) ? 1 : 0;
            n = V[0][k]*e1 + V[1][k]*e2;
            x0 = b1;
        }
    }

    if ((n & nk) < 0) n = -n;

    scalar sum = 0;
    forAll(q, i)
    {
        const scalar dist = n & (q[i] - x0);
        sum += dist*dist;
    }

    res.n = n;
    res.x = x0;
    res.rms = std::sqrt(sum/scalar(q.size()));
    return res;
}

} // End namespace pointIso
} // End namespace Foam

#endif
