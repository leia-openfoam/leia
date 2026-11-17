/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, GPL-3.0-or-later.

Application
    davofTestPointIsoPlaneFit

Description
    Gate G2 of pointIsoAdvection (docs/davof/REPORT-swept.md "pointIsoAdvector"):
    the plane fit pointIso::fitPlane (pointIsoPlaneFit.H).
    (a) coplanar points, no boundary point: the recovered plane (normal and
        point) against the plane the points were drawn on, tolerance 1e-14;
    (b) one, two and three boundary points on a plane, the other points
        perturbed off the plane: the fitted plane passes through the
        boundary points to 1e-14 and its rms residual is not larger than
        that of any of 1000 random planes through the same boundary points
        (one boundary point: random normals; two: random planes about the
        line; three: the plane is unique and is compared with the plane
        through them).
    The orientation n'.nk > 0 is checked in every trial.

    Reference:
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

#include "pointIsoPlaneFit.H"
#include <random>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace Foam;

static std::mt19937_64 rng(20261116);

static double uni(double a, double b)
{
    return std::uniform_real_distribution<double>(a, b)(rng);
}

static vector randomUnit()
{
    std::normal_distribution<double> g(0, 1);
    vector v(g(rng), g(rng), g(rng));
    return v/mag(v);
}

static void basis(const vector& n, vector& e1, vector& e2)
{
    vector axis(1, 0, 0);
    if (mag(n.y()) < mag(n.x()) && mag(n.y()) <= mag(n.z()))
    {
        axis = vector(0, 1, 0);
    }
    else if (mag(n.z()) < mag(n.x()) && mag(n.z()) < mag(n.y()))
    {
        axis = vector(0, 0, 1);
    }
    e1 = n ^ axis;
    e1 /= mag(e1);
    e2 = n ^ e1;
}

static double rmsResidual
(
    const std::vector<point>& q,
    const vector& n,
    const point& x0
)
{
    double s = 0;
    for (const point& p : q)
    {
        const double d = n & (p - x0);
        s += d*d;
    }
    return std::sqrt(s/double(q.size()));
}

int main()
{
    bool pass = true;

    // (a) coplanar points
    {
        double maxN = 0, maxX = 0;
        int nWrongSide = 0;
        for (int t = 0; t < 1000; ++t)
        {
            const vector n = randomUnit();
            const point x(uni(-1, 1), uni(-1, 1), uni(-1, 1));
            vector e1, e2;
            basis(n, e1, e2);
            const double sgn = (t % 2) ? 1.0 : -1.0;

            List<point> q(8);
            List<bool> isB(8, false);
            std::vector<point> qv;
            forAll(q, i)
            {
                q[i] = x + uni(-1, 1)*e1 + uni(-1, 1)*e2;
                qv.push_back(q[i]);
            }
            const pointIso::fit f = pointIso::fitPlane(q, isB, sgn*n);
            if ((f.n & (sgn*n)) <= 0) ++nWrongSide;
            maxN = std::max(maxN, double(mag(f.n - sgn*n)));
            maxX = std::max(maxX, std::fabs(double(n & (f.x - x))));
        }
        const bool ok = (maxN <= 1e-14 && maxX <= 1e-14 && nWrongSide == 0);
        pass = pass && ok;
        std::printf
        (
            "G2a coplanar (1000 planes, 8 points): max |n' - n| = %.3e, "
            "max |n.(x' - x)| = %.3e, wrong side %d -> %s\n",
            maxN, maxX, nWrongSide, ok ? "PASS" : "FAIL"
        );
    }

    // (b) boundary points
    for (int nb = 1; nb <= 3; ++nb)
    {
        const int nTrials = 300;
        double maxThrough = 0;
        double maxExcess = 0;   // max of (rms_fit - min rms_random)
        int nViol = 0, nWrongSide = 0;
        double worstDist = 0, worstQual = 0, maxThroughWell = 0, maxScaled = 0;
        int nAbove = 0, nWell = 0;
        for (int t = 0; t < nTrials; ++t)
        {
            const vector n = randomUnit();
            const point x(uni(-1, 1), uni(-1, 1), uni(-1, 1));
            vector e1, e2;
            basis(n, e1, e2);

            const int m = 8;
            List<point> q(m);
            List<bool> isB(m, false);
            std::vector<point> qv;
            for (int i = 0; i < m; ++i)
            {
                q[i] = x + uni(-1, 1)*e1 + uni(-1, 1)*e2;
                if (i < nb)
                {
                    isB[i] = true;
                }
                else
                {
                    q[i] += 0.1*uni(-1, 1)*n;
                }
                qv.push_back(q[i]);
            }

            const pointIso::fit f = pointIso::fitPlane(q, isB, n);
            if ((f.n & n) <= 0) ++nWrongSide;
            if (nb == 3)
            {
                // the conditioning of the triple: 2|cross|/(longest edge)^2
                const double e01 = mag(q[1] - q[0]);
                const double e02 = mag(q[2] - q[0]);
                const double e12 = mag(q[2] - q[1]);
                const double lg = std::max(e01, std::max(e02, e12));
                const double qual = 2*mag((q[1] - q[0]) ^ (q[2] - q[0]))/(lg*lg);
                double dd = 0;
                for (int i = 0; i < 3; ++i)
                {
                    dd = std::max(dd, std::fabs(double(f.n & (q[i] - f.x))));
                }
                if (dd > worstDist) { worstDist = dd; worstQual = qual; }
                if (dd > 1e-14) ++nAbove;
                if (qual > 0.1)
                {
                    ++nWell;
                    maxThroughWell = std::max(maxThroughWell, dd);
                }
                // the distance scaled by the round-off amplification 1/quality
                maxScaled = std::max(maxScaled, dd*qual);
            }
            for (int i = 0; i < nb; ++i)
            {
                maxThrough =
                    std::max
                    (
                        maxThrough,
                        std::fabs(double(f.n & (q[i] - f.x)))
                    );
            }
            const double rFit = rmsResidual(qv, f.n, f.x);

            double rMin = 1e300;
            for (int r = 0; r < 1000; ++r)
            {
                if (nb == 1)
                {
                    const vector nr = randomUnit();
                    rMin = std::min(rMin, rmsResidual(qv, nr, q[0]));
                }
                else if (nb == 2)
                {
                    vector d(q[1] - q[0]);
                    d /= mag(d);
                    vector u1, u2;
                    basis(d, u1, u2);
                    const double phi = uni(0, 2*M_PI);
                    const vector nr = std::cos(phi)*u1 + std::sin(phi)*u2;
                    rMin = std::min(rMin, rmsResidual(qv, nr, q[0]));
                }
                else
                {
                    vector nr((q[1] - q[0]) ^ (q[2] - q[0]));
                    nr /= mag(nr);
                    rMin = std::min(rMin, rmsResidual(qv, nr, q[0]));
                }
            }
            maxExcess = std::max(maxExcess, rFit - rMin);
            if (rFit > rMin + 1e-14) ++nViol;
        }
        const bool ok = (maxThrough <= 1e-14 && nViol == 0 && nWrongSide == 0);
        if (nb == 3)
        {
            std::printf
            (
                "G2b diagnostic, triples: %d of %d trials above 1e-14; worst "
                "distance %.3e at conditioning %.3e; %d trials with "
                "conditioning > 0.1, their max distance %.3e; max of "
                "distance x conditioning %.3e\n",
                nAbove, nTrials, worstDist, worstQual, nWell, maxThroughWell,
                maxScaled
            );
        }
        pass = pass && ok;
        std::printf
        (
            "G2b %d boundary point(s) (%d trials, 1000 random planes each): "
            "max distance of a boundary point to the fitted plane = %.3e, "
            "max (rms_fit - min rms_random) = %.3e, trials with the fit "
            "worse than a random plane = %d, wrong side %d -> %s\n",
            nb, nTrials, maxThrough, maxExcess, nViol, nWrongSide,
            ok ? "PASS" : "FAIL"
        );
    }

    std::printf("G2 %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
