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
    This file is part of OpenFOAM, licensed GPL-3.0-or-later.

Application
    davofTestCutVolume

Description
    The cell cut on warped faces (2026-11-14, docs/davof/REPORT-swept.md "The
    cell cut on warped faces"). For a list of cells and unit normals n the wet
    volume V(c) of the cell above the plane n.(x - C) = c (C the cell centre,
    the convention of cutCellPLIC: base = C + n c, the wet side is n.(x - base)
    > 0) is evaluated on a grid of c (nGrid points across the extent of the
    cell along n, plus the points L -+ 1e-9 extent at every vertex level L) by

        leia  : Foam::davofCutCellPLIC with the grazing area 10*SMALL
                (the copy in src/davof/reconstruction; davofCellCut routes to
                OpenFOAM's class by default)
        of    : Foam::cutCellPLIC (called as plicRDF calls it)
        iso   : Foam::cutCellIso (f = n.x at the points, cutValue = n.C + c)
        ref   : the reference volume of this application ONLY: each face is
                split once into the triangles (F, v_i, v_i+1) with F the face
                centre of the mesh (a triangular face is one triangle), each
                triangle with the cell centre mesh.C() forms a tetrahedron,
                each tetrahedron is clipped by the half-space in closed form,
                the signed volumes are summed
        refP0 : the same with F the average of the face points (the apex of
                primitiveMeshTools::makeFaceCentresAndAreas)

    Output (one line per quantity, comma separated; "T1".."T9" in the log):
    the jumps at the vertex levels, the non-monotone grid intervals, the
    deviation of the full-wet limit from V_k, the deviation of every method
    from ref, the root existence of alpha V_k, the bisection of ref, the cost,
    the tetrahedron statistics, the warp scatter (CSV) and its correlation.

    References:
    \verbatim
        T. Maric,
        Iterative Volume-of-Fluid interface positioning in general
        polyhedrons with Consecutive Cubic Spline interpolation,
        Journal of Computational Physics: X 11 (2021) 100093,
        doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227

        H. Scheufler, J. Roenby,
        Accurate and efficient surface reconstruction from volume fraction
        data on general meshes,
        Journal of Computational Physics 383 (2019) 1-23,
        doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellPLIC.H"
#include "cutCellIso.H"
#include "davofCutCellPLIC.H"
// the copy of the cell cut is compiled into this application (one translation unit, no leia library)
#include "davofCutCellPLIC.C"
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace Foam;

namespace
{

constexpr int nMeth = 5;
const char* methName[nMeth] = {"leia", "of", "iso", "ref", "refP0"};
constexpr int iRef = 3;
constexpr int iRefP0 = 4;

inline scalar tetVolAbs(const point& a, const point& b, const point& c, const point& d)
{
    return mag(((b - a) ^ (c - a)) & (d - a))/6.0;
}

// The volume of {s > 0} in the tetrahedron x[0..3], s the linear function at the vertices
scalar clipTet(const point x[4], const scalar s[4])
{
    int pos[4], neg[4], np = 0, nn = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (s[i] > 0) pos[np++] = i; else neg[nn++] = i;
    }
    if (np == 0) return 0;
    const scalar full = tetVolAbs(x[0], x[1], x[2], x[3]);
    if (np == 4) return full;
    auto cutPt = [&](int i, int j)
    {
        return point(x[i] + (s[i]/(s[i] - s[j]))*(x[j] - x[i]));
    };
    if (np == 1)
    {
        const int a = pos[0];
        return tetVolAbs(x[a], cutPt(a, neg[0]), cutPt(a, neg[1]), cutPt(a, neg[2]));
    }
    if (np == 3)
    {
        const int b = neg[0];
        return full - tetVolAbs(x[b], cutPt(b, pos[0]), cutPt(b, pos[1]), cutPt(b, pos[2]));
    }
    // two positive (a, b), two non-positive (c, d): a triangular prism
    // (a, c1, d1) - (b, c2, d2), split into three tetrahedra
    const int a = pos[0], b = pos[1], c = neg[0], d = neg[1];
    const point c1 = cutPt(a, c), d1 = cutPt(a, d), c2 = cutPt(b, c), d2 = cutPt(b, d);
    return tetVolAbs(x[a], c1, d1, x[b])
         + tetVolAbs(c1, d1, x[b], c2)
         + tetVolAbs(d1, x[b], c2, d2);
}

struct Tets
{
    std::vector<std::array<point, 4>> x;   // (C, v_i, v_i+1, F)
    std::vector<int> sgn;                  // sign of the signed volume
    scalar Vtri = 0;
    label nNeg = 0;
};

Tets buildTets(const fvMesh& mesh, const label celli, const bool apexAverage)
{
    Tets t;
    const point& C = mesh.C()[celli];
    const pointField& P = mesh.points();
    for (const label fi : mesh.cells()[celli])
    {
        const face& f = mesh.faces()[fi];
        const int s = (mesh.faceOwner()[fi] == celli) ? 1 : -1;
        const label nf = f.size();
        auto add = [&](const point& a, const point& b, const point& F)
        {
            const scalar sv = s*(((a - F) ^ (b - F)) & (F - C))/6.0;
            std::array<point, 4> tt{{C, a, b, F}};
            t.x.push_back(tt);
            t.sgn.push_back(sv < 0 ? -1 : 1);
            t.Vtri += sv;
            if (sv < 0) ++t.nNeg;
        };
        if (nf == 3)
        {
            add(P[f[1]], P[f[2]], P[f[0]]);
            continue;
        }
        point F = mesh.faceCentres()[fi];
        if (apexAverage)
        {
            F = Zero;
            for (const label pi : f) F += P[pi];
            F /= nf;
        }
        for (label i = 0; i < nf; ++i)
        {
            add(P[f[i]], P[f[(i + 1) % nf]], F);
        }
    }
    return t;
}

// the wet volume of the half-space n.(x - C) - c > 0
scalar volRef(const Tets& t, const vector& n, const point& C, const scalar c)
{
    scalar V = 0;
    const size_t nt = t.x.size();
    for (size_t k = 0; k < nt; ++k)
    {
        const point* x = t.x[k].data();
        scalar s[4];
        for (int i = 0; i < 4; ++i) s[i] = (n & (x[i] - C)) - c;
        V += t.sgn[k]*clipTet(x, s);
    }
    return V;
}

double quantile(std::vector<double> v, const double q)
{
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[size_t(q*(v.size() - 1))];
}

double pearson(const std::vector<double>& a, const std::vector<double>& b)
{
    const size_t n = a.size();
    if (n < 2) return 0;
    double ma = 0, mb = 0;
    for (size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= n; mb /= n;
    double sab = 0, saa = 0, sbb = 0;
    for (size_t i = 0; i < n; ++i)
    {
        sab += (a[i] - ma)*(b[i] - mb);
        saa += (a[i] - ma)*(a[i] - ma);
        sbb += (b[i] - mb)*(b[i] - mb);
    }
    return (saa > 0 && sbb > 0) ? sab/std::sqrt(saa*sbb) : 0;
}

struct GridPt { double c; int tag; int lev; };   // tag 0 plain, 1 left of level, 2 right of level

struct Acc
{
    std::vector<double> jumps;       // per vertex level
    std::vector<double> jumpsInner;  // per vertex level except the lowest and the highest
    std::vector<double> jumpsFirst;  // the lowest level (the gap between V_k and the cut's own full-wet limit)
    long nSt0 = 0;                   // (cell, normal) pairs for which the cut at c_min + d is a cut (status 0)
    std::vector<double> fullDev;     // per (cell, normal)
    std::vector<double> diffRef;     // per (cell, normal): max over the grid
    std::vector<double> cellJump;    // per selected cell: max jump
    std::vector<double> cellFull;    // per selected cell: max full-wet deviation
    long nonMono = 0, nonMonoJump = 0, nIntervals = 0, nUpJumpLevels = 0;
    double seconds = 0;
    long nEval = 0;
    // root existence for alpha V_k: [alpha][0 end gap, 1 inner jump, 2 other, 3 root found]
    long root[6][4] = {};
};

} // anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote("The cell cut on warped faces (davofCutCellPLIC, cutCellPLIC, cutCellIso, a reference volume)");
    argList::addOption("nCells", "label", "cells of the mesh (default 200)");
    argList::addOption("nRandom", "label", "random normals (default 100)");
    argList::addOption("nGrid", "label", "uniform grid points of c (default 201)");
    argList::addOption("tag", "word", "mesh label of the output lines");
    argList::addOption("extra", "file", "failure dump of the CCS round: its cells of this tag are included");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label nSelWanted = args.getOrDefault<label>("nCells", 200);
    const label nRand = args.getOrDefault<label>("nRandom", 100);
    const label nGrid = args.getOrDefault<label>("nGrid", 201);
    const word tag = args.getOrDefault<word>("tag", "mesh");
    const double dRel = 1e-9;      // the offset at a vertex level, relative to the extent
    const double levTol = 1e-8;    // two levels are one when they differ by less than this
    const double monoTol = 1e-13;  // relative to V_k

    // the cells: the extra ones of the dump, then evenly strided ones up to nSelWanted
    std::vector<label> cells;
    std::vector<double> extraVk;
    if (args.found("extra"))
    {
        std::ifstream in(args.get<fileName>("extra"));
        std::string line;
        const std::string key = "FAILDUMP," + std::string(tag) + ",cell ";
        while (std::getline(in, line))
        {
            if (line.compare(0, key.size(), key) != 0) continue;
            const label ci = std::stoi(line.substr(key.size()));
            double vk = -1;
            const size_t p = line.find(",V_k ");
            if (p != std::string::npos) vk = std::stod(line.substr(p + 5));
            if (ci < mesh.nCells() && std::find(cells.begin(), cells.end(), ci) == cells.end())
            {
                cells.push_back(ci);
                extraVk.push_back(vk);
            }
        }
    }
    const label nExtra = cells.size();
    {
        const label nStride = max(label(0), min(nSelWanted, mesh.nCells()) - nExtra);
        for (label is = 0; is < nStride; ++is)
        {
            const label ci = (nStride == 1) ? 0 : label((long long)(is)*mesh.nCells()/nStride);
            if (std::find(cells.begin(), cells.end(), ci) == cells.end()) cells.push_back(ci);
        }
    }
    const label nSel = cells.size();

    // the normals: 26 aligned, then random ones (the seed of the CCS gate)
    DynamicList<vector> normals;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <= 1; ++k)
            {
                if (i == 0 && j == 0 && k == 0) continue;
                vector v(i, j, k);
                v /= mag(v);
                normals.append(v);
            }
    {
        std::mt19937_64 rng(20261114);
        std::normal_distribution<double> gauss(0.0, 1.0);
        for (label q = 0; q < nRand; ++q)
        {
            vector v(gauss(rng), gauss(rng), gauss(rng));
            v /= mag(v);
            normals.append(v);
        }
    }

    const scalarField& V = mesh.V();
    const pointField& P = mesh.points();
    Info.stream().precision(10);

    // the cell data: tets, warp measure
    std::vector<Tets> tCf(nSel), tP0(nSel);
    std::vector<double> warp(nSel, 0);
    label nTetsMin = labelMax, nTetsMax = 0;
    double nTetsSum = 0;
    label nNegCf = 0, nNegP0 = 0, cellsNegCf = 0, cellsNegP0 = 0;
    std::vector<double> vtriDevCf(nSel), vtriDevP0(nSel);
    for (label ic = 0; ic < nSel; ++ic)
    {
        const label celli = cells[ic];
        tCf[ic] = buildTets(mesh, celli, false);
        tP0[ic] = buildTets(mesh, celli, true);
        const label nt = tCf[ic].x.size();
        nTetsMin = min(nTetsMin, nt);
        nTetsMax = max(nTetsMax, nt);
        nTetsSum += nt;
        nNegCf += tCf[ic].nNeg; if (tCf[ic].nNeg) ++cellsNegCf;
        nNegP0 += tP0[ic].nNeg; if (tP0[ic].nNeg) ++cellsNegP0;
        vtriDevCf[ic] = mag(tCf[ic].Vtri - V[celli])/V[celli];
        vtriDevP0[ic] = mag(tP0[ic].Vtri - V[celli])/V[celli];
        scalar w = 0;
        for (const label fi : mesh.cells()[celli])
        {
            const vector nf = normalised(mesh.faceAreas()[fi]);
            const scalar sc = Foam::sqrt(mag(mesh.faceAreas()[fi]));
            for (const label pi : mesh.faces()[fi])
            {
                w = max(w, mag((P[pi] - mesh.faceCentres()[fi]) & nf)/sc);
            }
        }
        warp[ic] = w;
    }

    Info<< "T0," << tag << ",nCellsMesh " << mesh.nCells() << ",nPoints " << mesh.nPoints()
        << ",nSelected " << nSel << ",nExtra " << nExtra << ",nNormals " << normals.size()
        << ",nGrid " << nGrid << ",warpMax " << *std::max_element(warp.begin(), warp.end())
        << ",warpMedian " << quantile(warp, 0.5)
        << ",cellsWarpAbove1e-12 " << std::count_if(warp.begin(), warp.end(), [](double w){ return w > 1e-12; })
        << nl;
    {
        label nBad = 0, nChk = 0;
        for (label ie = 0; ie < nExtra; ++ie)
        {
            if (extraVk[ie] > 0)
            {
                ++nChk;
                const double rel = mag(V[cells[ie]] - extraVk[ie])/extraVk[ie];
                if (rel > 1e-8) ++nBad;
            }
        }
        Info<< "EXTRACHK," << tag << ",cellsChecked " << nChk << ",V_k differs from the dump (rel > 1e-8) " << nBad << nl;
    }

    cutCellPLIC cutOF(mesh);
    davofCutCellPLIC cutLeia(mesh, 10*SMALL);
    scalarField fIso(mesh.nPoints(), 0.0);
    cutCellIso cutIso(mesh, fIso);

    const std::array<double, 6> alphas{{1e-2, 0.1, 0.5, 0.9, 0.99, 1.0 - 1e-4}};
    std::vector<Acc> acc(nMeth);
    for (Acc& a : acc)
    {
        a.cellJump.assign(nSel, 0);
        a.cellFull.assign(nSel, 0);
    }
    // bisection of ref: [alpha] residual max for T = alpha Vtri, iterations max, no root for alpha V_k > Vtri,
    // residual max for T = alpha V_k when solvable
    double bisRes[6] = {}, bisResVk[6] = {};
    int bisIt[6] = {};
    long bisNoRoot[6] = {}, bisSolvable[6] = {}, bisPairs = 0;
    long nPairs = 0;
    double tBuild = 0;

    for (const vector& n : normals)
    {
        for (label ip = 0; ip < mesh.nPoints(); ++ip) fIso[ip] = n & P[ip];

        for (label ic = 0; ic < nSel; ++ic)
        {
            const label celli = cells[ic];
            const scalar Vk = V[celli];
            const point& C = mesh.C()[celli];
            const scalar nC = n & C;

            // the vertex levels, relative to the cell centre
            std::vector<double> lv;
            for (const label pi : mesh.cellPoints(celli)) lv.push_back((P[pi] - C) & n);
            std::sort(lv.begin(), lv.end());
            const double cmin = lv.front(), cmax = lv.back();
            const double ext = cmax - cmin;
            const double dlt = dRel*ext;
            std::vector<double> lev;
            for (double l : lv)
            {
                if (lev.empty() || l - lev.back() > levTol*ext) lev.push_back(l);
            }

            std::vector<GridPt> g;
            for (size_t j = 0; j < lev.size(); ++j)
            {
                g.push_back({lev[j] - dlt, 1, int(j)});
                g.push_back({lev[j] + dlt, 2, int(j)});
            }
            for (label q = 0; q < nGrid; ++q)
            {
                const double c = cmin + ext*q/double(nGrid - 1);
                bool near = false;
                for (double l : lev) if (std::abs(c - l) < 3*dlt) { near = true; break; }
                if (!near) g.push_back({c, 0, -1});
            }
            std::sort(g.begin(), g.end(), [](const GridPt& a, const GridPt& b){ return a.c < b.c; });
            const size_t ng = g.size();
            std::vector<size_t> idxL(lev.size());
            for (size_t k = 0; k < ng; ++k) if (g[k].tag == 1) idxL[g[k].lev] = k;

            std::vector<std::vector<double>> Vm(nMeth, std::vector<double>(ng));
            for (int m = 0; m < nMeth; ++m)
            {
                const auto t0 = std::chrono::steady_clock::now();
                for (size_t k = 0; k < ng; ++k)
                {
                    const double c = g[k].c;
                    double v = 0;
                    if (m == 0 || m == 1 || m == 2)
                    {
                        label st;
                        if (m == 0) st = cutLeia.calcSubCell(celli, c, n);
                        else if (m == 1) st = cutOF.calcSubCell(celli, c, n);
                        else st = cutIso.calcSubCell(celli, nC + c);
                        if (st == 0)
                        {
                            v = (m == 0) ? cutLeia.subCellVolume()
                              : (m == 1) ? cutOF.subCellVolume() : cutIso.subCellVolume();
                        }
                        else v = (st < 0) ? Vk : 0;
                    }
                    else
                    {
                        v = volRef(m == iRef ? tCf[ic] : tP0[ic], n, C, c);
                    }
                    Vm[m][k] = v;
                }
                const auto t1 = std::chrono::steady_clock::now();
                acc[m].seconds += std::chrono::duration<double>(t1 - t0).count();
                acc[m].nEval += ng;
            }
            ++nPairs;

            for (int m = 0; m < nMeth; ++m)
            {
                Acc& a = acc[m];
                const std::vector<double>& v = Vm[m];
                // jumps at the vertex levels
                for (size_t j = 0; j < lev.size(); ++j)
                {
                    const double jmp = std::abs(v[idxL[j] + 1] - v[idxL[j]])/Vk;
                    a.jumps.push_back(jmp);
                    if (j > 0 && j + 1 < lev.size()) a.jumpsInner.push_back(jmp);
                    if (j == 0) a.jumpsFirst.push_back(jmp);
                    a.cellJump[ic] = std::max(a.cellJump[ic], jmp);
                    if (v[idxL[j] + 1] - v[idxL[j]] > monoTol*Vk) ++a.nUpJumpLevels;
                }
                // non-monotone intervals (the wet volume is non-increasing in c)
                for (size_t k = 0; k + 1 < ng; ++k)
                {
                    ++a.nIntervals;
                    if (v[k + 1] - v[k] > monoTol*Vk)
                    {
                        ++a.nonMono;
                        if (g[k].tag == 1 && g[k + 1].tag == 2 && g[k].lev == g[k + 1].lev) ++a.nonMonoJump;
                    }
                }
                // the full-wet limit
                double fd;
                if (m == iRef) fd = std::abs(tCf[ic].Vtri - Vk)/Vk;
                else if (m == iRefP0) fd = std::abs(tP0[ic].Vtri - Vk)/Vk;
                else
                {
                    const double c = cmin + dlt;
                    label st;
                    if (m == 0) st = cutLeia.calcSubCell(celli, c, n);
                    else if (m == 1) st = cutOF.calcSubCell(celli, c, n);
                    else st = cutIso.calcSubCell(celli, nC + c);
                    if (st == 0) ++a.nSt0;
                    double vf;
                    if (st == 0)
                    {
                        const double A = (m == 0) ? mag(cutLeia.faceArea())
                                       : (m == 1) ? mag(cutOF.faceArea()) : mag(cutIso.faceArea());
                        const double vv = (m == 0) ? cutLeia.subCellVolume()
                                        : (m == 1) ? cutOF.subCellVolume() : cutIso.subCellVolume();
                        vf = vv + A*dlt;
                    }
                    else vf = (st < 0) ? Vk : 0;
                    fd = std::abs(vf - Vk)/Vk;
                }
                a.fullDev.push_back(fd);
                a.cellFull[ic] = std::max(a.cellFull[ic], fd);
                // deviation from ref over the grid
                if (m != iRef)
                {
                    double dmax = 0;
                    for (size_t k = 0; k < ng; ++k) dmax = std::max(dmax, std::abs(v[k] - Vm[iRef][k])/Vk);
                    a.diffRef.push_back(dmax);
                }
                // root existence of alpha V_k
                for (int ia = 0; ia < 6; ++ia)
                {
                    const double T = alphas[ia]*Vk;
                    long contCross = 0;
                    int firstJumpLev = -1;
                    for (size_t k = 0; k + 1 < ng; ++k)
                    {
                        const bool cross = (v[k] - T)*(v[k + 1] - T) <= 0;
                        if (!cross) continue;
                        // a level pair is a jump when V changes by more than 1e-6 V_k across it (the continuous
                        // change across 2 d is below 1e-8 V_k); otherwise the interval is continuous
                        if (g[k].tag == 1 && g[k + 1].tag == 2 && g[k].lev == g[k + 1].lev
                            && std::abs(v[k + 1] - v[k]) > 1e-6*Vk)
                        {
                            if (firstJumpLev < 0) firstJumpLev = g[k].lev;
                        }
                        else ++contCross;
                    }
                    if (contCross > 0) ++a.root[ia][3];
                    else if (firstJumpLev == 0) ++a.root[ia][0];
                    else if (firstJumpLev > 0) ++a.root[ia][1];
                    else ++a.root[ia][2];
                }
            }

            // bisection of ref (apex F = face centre)
            {
                const auto t0 = std::chrono::steady_clock::now();
                ++bisPairs;
                for (int ia = 0; ia < 6; ++ia)
                {
                    for (int var = 0; var < 2; ++var)
                    {
                        const double T = alphas[ia]*(var == 0 ? tCf[ic].Vtri : Vk);
                        if (var == 1 && T > tCf[ic].Vtri) { ++bisNoRoot[ia]; continue; }
                        if (var == 1) ++bisSolvable[ia];
                        double lo = cmin - dlt, hi = cmax + dlt;   // V(lo) = Vtri >= T >= V(hi) = 0
                        double res = 0;
                        int it = 0;
                        for (; it < 200; ++it)
                        {
                            const double mid = 0.5*(lo + hi);
                            const double vm = volRef(tCf[ic], n, C, mid);
                            res = std::abs(vm - T)/Vk;
                            if (res <= 1e-13) break;
                            if (vm > T) lo = mid; else hi = mid;
                            if (hi - lo < 1e-16*std::max(1.0, std::abs(mid))) break;
                        }
                        if (var == 0)
                        {
                            bisRes[ia] = std::max(bisRes[ia], res);
                            bisIt[ia] = std::max(bisIt[ia], it);
                        }
                        else bisResVk[ia] = std::max(bisResVk[ia], res);
                    }
                }
                tBuild += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            }
        }
    }

    // ---------------------------------------------------------------- output
    Info<< "T1,header,tag,method,nPairs,nLevels,jumpMax,jumpMedian,jumpP99,nJumpAbove1e-6,nUpwardJumpLevels,"
        << "nIntervals,nNonMonotone,nNonMonotoneAtJump,fullDevMax,fullDevMedian,fullDevP99,nFullDevAbove1e-13,"
        << "diffRefMax,diffRefMedian,diffRefP99" << nl;
    for (int m = 0; m < nMeth; ++m)
    {
        const Acc& a = acc[m];
        const long nAbove = std::count_if(a.jumps.begin(), a.jumps.end(), [](double j){ return j > 1e-6; });
        const long nFull = std::count_if(a.fullDev.begin(), a.fullDev.end(), [](double j){ return j > 1e-13; });
        Info<< "T1," << tag << "," << methName[m] << "," << nPairs << "," << a.jumps.size() << ","
            << quantile(a.jumps, 1.0) << "," << quantile(a.jumps, 0.5) << "," << quantile(a.jumps, 0.99) << ","
            << nAbove << "," << a.nUpJumpLevels << "," << a.nIntervals << "," << a.nonMono << "," << a.nonMonoJump << ","
            << quantile(a.fullDev, 1.0) << "," << quantile(a.fullDev, 0.5) << "," << quantile(a.fullDev, 0.99) << ","
            << nFull << ","
            << (m == iRef ? 0.0 : quantile(a.diffRef, 1.0)) << ","
            << (m == iRef ? 0.0 : quantile(a.diffRef, 0.5)) << ","
            << (m == iRef ? 0.0 : quantile(a.diffRef, 0.99)) << nl;
    }

    Info<< "T1b,header,tag,method,jumpInnerMax,jumpInnerMedian,jumpInnerP99,nInnerAbove1e-6,nInner,jumpFirstMax,jumpFirstMedian,"
        << "jumpFirstP99,pairsWithCutAt_cmin+d" << nl;
    for (int m = 0; m < 3; ++m)
    {
        const Acc& a = acc[m];
        const long nIn = std::count_if(a.jumpsInner.begin(), a.jumpsInner.end(), [](double j){ return j > 1e-6; });
        Info<< "T1b," << tag << "," << methName[m] << "," << quantile(a.jumpsInner, 1.0) << "," << quantile(a.jumpsInner, 0.5)
            << "," << quantile(a.jumpsInner, 0.99) << "," << nIn << "," << a.jumpsInner.size() << ","
            << quantile(a.jumpsFirst, 1.0) << "," << quantile(a.jumpsFirst, 0.5) << "," << quantile(a.jumpsFirst, 0.99) << ","
            << a.nSt0 << nl;
    }

    Info<< "T2,header,tag,method,alpha,nPairs,noRootEndGap,noRootInnerJump,noRootOther,rootFound" << nl;
    for (int m = 0; m < nMeth; ++m)
    {
        for (int ia = 0; ia < 6; ++ia)
        {
            const Acc& a = acc[m];
            Info<< "T2," << tag << "," << methName[m] << "," << alphas[ia] << "," << nPairs << ","
                << a.root[ia][0] << "," << a.root[ia][1] << "," << a.root[ia][2] << "," << a.root[ia][3] << nl;
        }
    }

    Info<< "T3,header,tag,alpha,nPairs,bisectTargetAlphaVtri_residualMax,iterationsMax,"
        << "nTargetAlphaVk_aboveVtri_noRoot,nTargetAlphaVk_solvable,bisectTargetAlphaVk_residualMax" << nl;
    for (int ia = 0; ia < 6; ++ia)
    {
        Info<< "T3," << tag << "," << alphas[ia] << "," << bisPairs << "," << bisRes[ia] << "," << bisIt[ia] << ","
            << bisNoRoot[ia] << "," << bisSolvable[ia] << "," << bisResVk[ia] << nl;
    }

    Info<< "T4,header,tag,method,nEval,secondsTotal,nsPerEval" << nl;
    for (int m = 0; m < nMeth; ++m)
    {
        Info<< "T4," << tag << "," << methName[m] << "," << acc[m].nEval << "," << acc[m].seconds << ","
            << 1e9*acc[m].seconds/std::max(1L, acc[m].nEval) << nl;
    }
    Info<< "T4b," << tag << ",ref bisection seconds (not in T4) " << tBuild << nl;

    Info<< "T5," << tag << ",tetsPerCellMin " << nTetsMin << ",mean " << nTetsSum/nSel << ",max " << nTetsMax
        << ",negativeTets(face centre) " << nNegCf << ",cellsWithNegative " << cellsNegCf
        << ",negativeTets(point average) " << nNegP0 << ",cellsWithNegative " << cellsNegP0 << nl;
    Info<< "T6," << tag << ",|Vtri - V_k|/V_k (face centre apex): max " << quantile(vtriDevCf, 1.0)
        << ",median " << quantile(vtriDevCf, 0.5) << ",nAbove1e-13 "
        << std::count_if(vtriDevCf.begin(), vtriDevCf.end(), [](double d){ return d > 1e-13; })
        << ";  (point average apex): max " << quantile(vtriDevP0, 1.0) << ",median " << quantile(vtriDevP0, 0.5)
        << ",nAbove1e-13 " << std::count_if(vtriDevP0.begin(), vtriDevP0.end(), [](double d){ return d > 1e-13; }) << nl;

    // the scatter and the correlation of the largest jump of a cell against its warp
    {
        std::ofstream csv("scatter-" + std::string(tag) + ".csv");
        csv.precision(10);
        csv << "cell,warp,V_k,jumpLeia,jumpOF,jumpIso,jumpRef,jumpRefP0,fullDevLeia,fullDevOF,fullDevIso,VtriDevCf,VtriDevP0\n";
        for (label ic = 0; ic < nSel; ++ic)
        {
            csv << cells[ic] << "," << warp[ic] << "," << V[cells[ic]];
            for (int m = 0; m < nMeth; ++m) csv << "," << acc[m].cellJump[ic];
            for (int m = 0; m < 3; ++m) csv << "," << acc[m].cellFull[ic];
            csv << "," << vtriDevCf[ic] << "," << vtriDevP0[ic] << "\n";
        }
    }
    Info<< "T7,header,tag,method,pearsonJumpVsWarp,pearsonLog10,pearsonFullDevVsWarp,pearsonLog10Full" << nl;
    for (int m = 0; m < nMeth; ++m)
    {
        std::vector<double> lw(nSel), lj(nSel), lf(nSel);
        for (label ic = 0; ic < nSel; ++ic)
        {
            lw[ic] = std::log10(std::max(warp[ic], 1e-16));
            lj[ic] = std::log10(std::max(acc[m].cellJump[ic], 1e-16));
            lf[ic] = std::log10(std::max(acc[m].cellFull[ic], 1e-16));
        }
        Info<< "T7," << tag << "," << methName[m] << "," << pearson(warp, acc[m].cellJump) << ","
            << pearson(lw, lj) << "," << pearson(warp, acc[m].cellFull) << "," << pearson(lw, lf) << nl;
    }
    {
        std::vector<double> lw(nSel), lj(nSel);
        for (label ic = 0; ic < nSel; ++ic)
        {
            lw[ic] = std::log10(std::max(warp[ic], 1e-16));
            lj[ic] = std::log10(std::max(vtriDevCf[ic], 1e-16));
        }
        Info<< "T7b," << tag << ",VtriDevCf vs warp: pearson " << pearson(warp, vtriDevCf)
            << ",log10 " << pearson(lw, lj) << nl;
    }
    Info<< "TDONE," << tag << nl;
    return 0;
}

// ************************************************************************* //
