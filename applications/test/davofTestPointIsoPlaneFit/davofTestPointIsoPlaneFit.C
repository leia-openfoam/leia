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
#include <cfenv>
#include <fstream>
#include <sstream>
#include <string>

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

static void report(const char* tag, const std::vector<point>& qv, const std::vector<bool>& bv, const vector& nk)
{
    List<point> q(qv.size());
    List<bool> b(qv.size(), false);
    forAll(q, i) { q[i] = qv[i]; b[i] = bv[i]; }
    const pointIso::fit f = pointIso::fitPlane(q, b, nk);
    std::printf
    (
        "%s: points %d, boundary %d, n' = (%.17g %.17g %.17g), n'.nk = %.17g, rms = %.3e, finite %d\n",
        tag, int(qv.size()), int(f.nBoundary), f.n.x(), f.n.y(), f.n.z(),
        double(f.n & nk), f.rms,
        int(std::isfinite(f.n.x()) && std::isfinite(f.n.y()) && std::isfinite(f.n.z()))
    );
}

static void robustness()
{
    std::vector<bool> nob(16, false);
    // exactly coplanar, general plane
    {
        const vector n = randomUnit();
        vector e1, e2;
        basis(n, e1, e2);
        std::vector<point> q;
        for (int i = 0; i < 6; ++i) q.push_back(point(0.3, 0.2, 0.1) + uni(-1, 1)*e1 + uni(-1, 1)*e2);
        report("R1 exactly coplanar", q, nob, n);
    }
    // coplanar in an axis plane (x = const): two covariance entries exactly zero
    {
        std::vector<point> q;
        for (int i = 0; i < 6; ++i) q.push_back(point(0.25, uni(-1, 1), uni(-1, 1)));
        report("R2 axis plane x = const", q, nob, vector(1, 0, 0));
        std::vector<point> q2;
        for (int i = 0; i < 6; ++i) q2.push_back(point(uni(-1, 1), 0.5, uni(-1, 1)));
        report("R2 axis plane y = const", q2, nob, vector(0, 1, 0));
        std::vector<point> q3;
        for (int i = 0; i < 4; ++i) q3.push_back(point(uni(-1e-4, 1e-4), uni(-1e-4, 1e-4), 0.5));
        report("R2 axis plane z = const, 1e-4 extent", q3, nob, vector(0, 0, 1));
    }
    // three points only
    {
        std::vector<point> q;
        for (int i = 0; i < 3; ++i) q.push_back(point(uni(-1, 1), uni(-1, 1), uni(-1, 1)));
        report("R3 three points", q, nob, vector(0, 0, 1));
    }
    // all points equal up to 1e-16
    {
        std::vector<point> q;
        for (int i = 0; i < 5; ++i) q.push_back(point(0.1 + uni(-1e-16, 1e-16), 0.2 + uni(-1e-16, 1e-16), 0.3 + uni(-1e-16, 1e-16)));
        report("R4 points equal up to 1e-16", q, nob, vector(0, 0, 1));
        std::vector<point> q2(5, point(0.1, 0.2, 0.3));
        report("R4 points exactly equal", q2, nob, vector(0, 0, 1));
    }
}

// the points printed by pointIsoAdvection (key pointIsoPrintPoints)
static void fromFile(const char* name)
{
    std::ifstream in(name);
    std::string line;
    int n = 0;
    while (std::getline(in, line))
    {
        const std::size_t p = line.find("pointIsoPoints cell ");
        if (p == std::string::npos) continue;
        std::istringstream is(line.substr(p));
        std::string w1, w2, w3, colon;
        long cell;
        double nx, ny, nz;
        is >> w1 >> w2 >> cell >> w3 >> nx >> ny >> nz >> colon;
        std::vector<point> q;
        std::vector<bool> b;
        double x, y, z;
        int bb;
        while (is >> x >> y >> z >> bb) { q.push_back(point(x, y, z)); b.push_back(bb != 0); }
        char tag[64];
        std::snprintf(tag, sizeof(tag), "file cell %ld", cell);
        std::printf("about to fit %s (%d points)\n", tag, int(q.size()));
        std::fflush(stdout);
        report(tag, q, b, vector(nx, ny, nz));
        ++n;
    }
    std::printf("file: %d cells fitted without a trap\n", n);
}

int main(int argc, char* argv[])
{
    // the trapping of the solvers (FOAM_SIGFPE): invalid, division by zero, overflow
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    if (argc > 1)
    {
        fromFile(argv[1]);
        return 0;
    }
    robustness();
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
        double maxScaledRel = 0;
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
                maxScaledRel = std::max(maxScaledRel, dd*qual/lg);
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
        const bool ok =
            (
                (nb == 3 ? maxScaledRel <= 1e-14 : maxThrough <= 1e-14)
             && nViol == 0 && nWrongSide == 0
            );
        if (nb == 3)
        {
            std::printf
            (
                "G2b triples, amended clause (distance x conditioning / longest edge) max = %.3e (<= 1e-14); %d of %d trials above 1e-14; worst "
                "distance %.3e at conditioning %.3e; %d trials with "
                "conditioning > 0.1, their max distance %.3e; max of "
                "distance x conditioning %.3e\n",
                maxScaledRel, nAbove, nTrials, worstDist, worstQual, nWell, maxThroughWell,
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
