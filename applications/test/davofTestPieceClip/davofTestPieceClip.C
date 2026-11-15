/*---------------------------------------------------------------------------*\
Application
    davofTestPieceClip

Description
    Standalone reproduction of the wrong piece (diag2 follow-up, 2026-11-15, docs/davof/REPORT-swept.md
    "Pieces path: the wrong piece"): one signed tetrahedron (apex = face centre, base = one triangle of
    the swept volume), the six planes of a hexahedral cell and the PLIC half-space, with the numbers of
    the dump (17 digits), through davof::sweptPieceExactVolume (as the flux calls it, with and without
    the axis-aligned cell planes), against the exact volume of the intersection of the ten half-spaces
    computed in long double by facets (every facet is a plane polygon clipped by all other planes;
    V = 1/3 sum_i (x_i - c) . n_i A_i), which shares no code with the clip under test.
    The cases: N 256 Co 0.2 new keys step 5040 face 3601 cell 1878 tet 6; N 128 Co 0.5 old keys step
    1660 face 12790 cell 6548 tet 5.
\*---------------------------------------------------------------------------*/

#include "point.H"
#include "vector.H"
#include "FixedList.H"
#include "List.H"
#include "davofSweptPieces.H"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Foam;
typedef long double ld;

struct V3 { ld x, y, z; };
static V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static V3 operator*(ld s, V3 a) { return {s*a.x, s*a.y, s*a.z}; }
static ld dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static V3 cross(V3 a, V3 b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
static ld nrm(V3 a) { return sqrtl(dot(a, a)); }

struct HS { V3 n, x0; };   // kept side: (p - x0).n <= 0, n a unit vector

// exact volume of the intersection of the half-spaces: the facets method in long double
static ld exactVolume(const std::vector<HS>& hs, V3 c)
{
    ld vol = 0;
    for (size_t i = 0; i < hs.size(); ++i)
    {
        // a large square on plane i
        V3 n = hs[i].n;
        V3 u = (std::fabs((double)n.x) < 0.9) ? V3{1, 0, 0} : V3{0, 1, 0};
        u = u - dot(u, n)*n; u = (1/nrm(u))*u;
        V3 v = cross(n, u);
        const ld R = 1;
        V3 o = hs[i].x0 - dot(hs[i].x0 - c, n)*n + dot(c - hs[i].x0, n)*n;   // x0 itself
        o = hs[i].x0;
        std::vector<V3> poly = {o + (-R)*u + (-R)*v, o + R*u + (-R)*v, o + R*u + R*v, o + (-R)*u + R*v};
        for (size_t j = 0; j < hs.size() && poly.size() >= 3; ++j)
        {
            if (j == i) continue;
            std::vector<V3> out;
            const size_t m = poly.size();
            for (size_t k = 0; k < m; ++k)
            {
                V3 p0 = poly[k], p1 = poly[(k + 1) % m];
                ld s0 = dot(p0 - hs[j].x0, hs[j].n), s1 = dot(p1 - hs[j].x0, hs[j].n);
                bool in0 = s0 <= 0, in1 = s1 <= 0;
                if (in0) out.push_back(p0);
                if (in0 != in1) out.push_back(p0 + (s0/(s0 - s1))*(p1 - p0));
            }
            poly = out;
        }
        if (poly.size() < 3) continue;
        V3 a{0, 0, 0};
        for (size_t k = 0; k < poly.size(); ++k)
        {
            V3 p0 = poly[k] - o, p1 = poly[(k + 1) % poly.size()] - o;
            a = a + 0.5L*cross(p0, p1);
        }
        const ld A = dot(a, n);   // signed area along n
        vol += dot(hs[i].x0 - c, n)*A/3;
    }
    return vol;
}

static V3 L(const point& p) { return {p.x(), p.y(), p.z()}; }

static void runCase
(
    const char* name, point xRef, point A, point B, point C,
    const std::vector<std::pair<vector, point>>& cellPlanes, vector plicN, point plicX0, double dumpValue
)
{
    std::cout << std::setprecision(17) << "== " << name << "\n";
    // the exact volumes (long double): the tet alone, with the cell, with the cell and the PLIC plane
    std::vector<HS> hs;
    V3 a = L(A), b = L(B), cc = L(C), x = L(xRef);
    V3 cen = 0.25L*(a + b + cc + x);
    auto faceHS = [&](V3 p, V3 q, V3 r, V3 inside)
    {
        V3 n = cross(q - p, r - p); n = (1/nrm(n))*n;
        if (dot(inside - p, n) > 0) n = (-1)*n;
        hs.push_back({n, p});
    };
    faceHS(a, b, cc, x); faceHS(x, a, b, cc); faceHS(x, b, cc, a); faceHS(x, cc, a, b);
    const ld vTet = exactVolume(hs, cen);
    for (auto& cp : cellPlanes)
    {
        V3 n = L(cp.first); n = (1/nrm(n))*n; hs.push_back({n, L(cp.second)});
    }
    const ld vCell = exactVolume(hs, cen);
    { V3 n = L(plicN); n = (1/nrm(n))*n; hs.push_back({n, L(plicX0)}); }
    const ld vPlic = exactVolume(hs, cen);
    std::cout << " exact (long double): V_tet " << (double)vTet << "  V after the cell planes " << (double)vCell
              << "  V after the PLIC plane " << (double)vPlic << "\n";

    // the code under test, as the flux calls it: T = the triangle (A, B, C), xRef the apex
    for (int useAxis = 1; useAxis >= 0; --useAxis)
    {
        List<FixedList<point, 3>> T(1);
        T[0][0] = A; T[0][1] = B; T[0][2] = C;
        List<davof::halfSpace> planes(cellPlanes.size());
        forAll(planes, i)
        {
            planes[i].n = cellPlanes[i].first/mag(cellPlanes[i].first);
            planes[i].x0 = cellPlanes[i].second;
            if (useAxis) davof::setAxis(planes[i]);
        }
        davof::halfSpace plic; plic.n = plicN/mag(plicN); plic.x0 = plicX0;
        const scalar vNoPlic = davof::sweptPieceExactVolume(T, planes, false, plic, xRef);
        const scalar vWith = davof::sweptPieceExactVolume(T, planes, true, plic, xRef);
        std::cout << " sweptPieceExactVolume, axis planes " << useAxis << ": cell planes only " << vNoPlic
                  << "  with the PLIC half-space " << vWith << "  (dump " << dumpValue << ")  |V|/V_tet(exact) "
                  << std::fabs((double)vWith)/(double)vTet << "\n";
    }
}

int main()
{
    using P = std::pair<vector, point>;
    // N 256, Co 0.2, new keys, step 5040, face 3601, cell 1878, tet 6 (vertices as printed by the dump)
    runCase
    (
        "N256 step 5040 face 3601 cell 1878 tet 6",
        point(0.513671875, 0, 0.11328125),
        point(0.515625, -0.001953125, 0.11328125),
        point(0.5157707701706159, 0.001953125, 0.1132839251108668),
        point(0.5157707701706159, -0.001953125, 0.1132839251108668),
        std::vector<P>
        {
            P(vector(1, 0, 0), point(0.515625, 0, 0.115234375)), P(vector(0, 0, 1), point(0.513671875, 0, 0.1171875)),
            P(vector(0, -1, 0), point(0.513671875, -0.001953125, 0.115234375)), P(vector(0, 1, 0), point(0.513671875, 0.001953125, 0.115234375)),
            P(vector(0, 0, -1), point(0.513671875, 0, 0.11328125)), P(vector(-1, 0, 0), point(0.51171875, 0, 0.115234375))
        },
        vector(-1, 0, -1.1494902637855374e-12), point(0.51562500000000167, 0, 0.11523437500000225),
        -1.1001853062195708e-05
    );
    // N 128, Co 0.5, old keys, step 1660, face 12790, cell 6548, tet 5
    runCase
    (
        "N128 step 1660 face 12790 cell 6548 tet 5",
        point(0.1640625, 0, 0.39453125),
        point(0.16389352325140855, 0, 0.39540185889605106),
        point(0.16389898456387664, 0.00390625, 0.39931544167643762),
        point(0.16388830089049897, 0.00390625, 0.39148827578950796),
        std::vector<P>
        {
            P(vector(1, 0, 0), point(0.1640625, 0, 0.40234375)), P(vector(0, 0, 1), point(0.16015625, 0, 0.40625)),
            P(vector(0, -1, 0), point(0.16015625, -0.00390625, 0.40234375)), P(vector(0, 1, 0), point(0.16015625, 0.00390625, 0.40234375)),
            P(vector(0, 0, -1), point(0.16015625, 0, 0.3984375)), P(vector(-1, 0, 0), point(0.15625, 0, 0.40234375))
        },
        vector(1.0796185446569768e-12, 0, 1), point(0.16015624999999578, 0, 0.39843749999999623),
        1.9317707577216788e-06
    );
    return 0;
}
