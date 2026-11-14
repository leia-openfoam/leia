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
    davofTestCcsPositioning

Description
    Gate A-unit of the CCS alignment with the published method (2026-11-14,
    docs/davof/REPORT-swept.md "CCS aligned with the published method").
    For a set of cells of the mesh of the case, 1000 random unit normals and
    the 26 axis, edge and corner aligned ones, and a ladder of volume
    fractions alpha, the plane is positioned by Foam::davof::ccsPositionPlane
    twice: the OLD setting (ccsEndInterval none, ccsTolerance absolute,
    tolerance 1e-14) and the NEW setting (ccsEndInterval cubic, ccsTolerance
    relative, tolerance 1e-12). The start is the linear guess
    p0 = pMin + (1 - alpha)(pMax - pMin) along the normal, the same for both.
    The volume error of the returned plane is measured by an independent cut:
    for alpha <= 1/2 the volume on the liquid side against alpha V_k, for
    alpha > 1/2 the volume of the complementary cut (normal and position
    mirrored) against (1 - alpha) V_k, so that the relative error
    |V_cut - alpha V_k| / (min(alpha, 1 - alpha) V_k) is not limited by the
    round-off of V_k - small.

    Reference:
    \verbatim
        T. Maric,
        Iterative Volume-of-Fluid interface positioning in general
        polyhedrons with Consecutive Cubic Spline interpolation,
        Journal of Computational Physics: X 11 (2021) 100093,
        doi:10.1016/j.jcpx.2021.100093, arXiv:2012.02227
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellPLIC.H"
#include "davofCcsPositioning.H"
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("Gate A-unit of the CCS alignment (davofCcsPositioning.H)");
    argList::addOption("nCells", "label", "cells of the mesh (default 200)");
    argList::addOption("nRandom", "label", "random normals (default 1000)");
    argList::addOption("tag", "word", "mesh label of the output lines");
    argList::addOption("maxIter", "label", "ccs maxIter (default 30)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label nSel = min(args.getOrDefault<label>("nCells", 200), mesh.nCells());
    const label nRand = args.getOrDefault<label>("nRandom", 1000);
    const word tag = args.getOrDefault<word>("tag", "mesh");
    const label maxIter = args.getOrDefault<label>("maxIter", 30);

    // the normals: 26 aligned, then the random ones (a fixed seed)
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

    // the ladder
    DynamicList<scalar> alphas;
    for (const double a : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 1e-20})
    {
        alphas.append(a);
    }
    for (const double e : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14})
    {
        alphas.append(1.0 - e);
    }

    std::vector<double> fullRels;
    label nDumpAll = 0;
    const label nDumpMax = 400;
    cutCellPLIC cutPos(mesh);
    cutCellPLIC cutChk(mesh);
    const scalarField& V = mesh.V();
    const pointField& P = mesh.points();

    Info.stream().precision(10);
    Info<< "RESULT,header,tag,arm,alpha,nCases,errMax,errMedian,errMean,nOver1e-10,errWMax,errWMedian,nOverW,nGross1e-3,errP99,"
        << "exit1,exit2,exit3,exit4,nNoPlane,cutsMean,cutsMax,endOk,endRej" << nl;

    for (label arm = 0; arm < 2; ++arm)
    {
        const bool endInterval = (arm == 1);
        const bool relativeTol = (arm == 1);
        const scalar tolRel = (arm == 1) ? 1e-12 : 1e-14;
        for (const scalar alpha : alphas)
        {
            std::vector<double> errs, errsW, errs4;
            label nDump = 0;
            label nFail = 0;
            label nExit[5] = {0, 0, 0, 0, 0};
            label nNoPlane = 0, cutsMax = 0, endOk = 0, endRej = 0;
            double cutsSum = 0;
            for (label is = 0; is < nSel; ++is)
            {
                const label celli = (nSel == 1) ? 0 : label((long long)(is)*mesh.nCells()/nSel);
                const scalar Vc = V[celli];
                const point& xc = mesh.C()[celli];
                for (const vector& n : normals)
                {
                    scalar pMin = GREAT, pMax = -GREAT;
                    for (const label pi : mesh.cellPoints(celli))
                    {
                        const scalar f = (P[pi] - xc) & n;
                        pMin = min(pMin, f);
                        pMax = max(pMax, f);
                    }
                    const scalar p0 = pMin + (1 - alpha)*(pMax - pMin);
                    Foam::davof::ccsRecord rec;
                    const label status = Foam::davof::ccsPositionPlane
                    (
                        mesh, cutPos, celli, alpha, n, p0, tolRel, maxIter, rec,
                        endInterval, relativeTol
                    );
                    const scalar p = rec.pFinal;
                    ++nExit[min(rec.exitCode, label(4))];
                    if (status != 0 && alpha >= 1e-16 && alpha <= 1 - 1e-14) ++nNoPlane;
                    cutsSum += rec.nEvals;
                    cutsMax = max(cutsMax, rec.nEvals);
                    if (rec.endState > 0) ++endOk;
                    if (rec.endState < 0) ++endRej;

                    // the independent volume error, as the gate states it:
                    // the cut volume V_cut on the liquid side of the returned
                    // plane (original frame) against alpha V_k; for alpha > 1/2
                    // also the volume W of the complementary cut (normal and
                    // position mirrored) against (1 - alpha) V_k. On a cell
                    // with planar faces V + W = V_k to round-off; on a warped
                    // face the two cuts differ by the fan triangulation.
                    double err, errW;
                    {
                        const label st = cutChk.calcSubCell(celli, p, n);
                        const scalar Vcut = (st == 0) ? cutChk.subCellVolume() : ((st < 0) ? Vc : 0);
                        err = mag(Vcut - alpha*Vc)/(min(alpha, 1 - alpha)*Vc);
                        errW = err;
                    }
                    if (alpha > 0.5)
                    {
                        const label st = cutChk.calcSubCell(celli, -p, -n);
                        const scalar W = (st == 0) ? cutChk.subCellVolume() : ((st < 0) ? Vc : 0);
                        errW = mag(W - (1 - alpha)*Vc)/((1 - alpha)*Vc);
                    }
                    errs.push_back(err);
                    errsW.push_back(errW);
                    if (rec.fullRel >= 0) fullRels.push_back(rec.fullRel);
                    if (rec.exitCode == 4) errs4.push_back(err);

                    // the located pre-existing defect: 20 gross failures of the
                    // OLD setting at alpha = 1e-2, dumped with the iteration history
                    if (arm == 0 && alpha == 1e-2 && err > 1e-3 && nFail < 20)
                    {
                        ++nFail;
                        DynamicList<scalar> hist, histBr;
                        Foam::davof::ccsRecord rec2;
                        rec2.hist = &hist;
                        rec2.histBr = &histBr;
                        Foam::davof::ccsPositionPlane
                        (
                            mesh, cutPos, celli, alpha, n, p0, tolRel, maxIter, rec2, false, false
                        );
                        // the non-planarity of the cell's faces: max over faces of
                        // max_i |(x_i - c_f).n_f| / sqrt(|S_f|)
                        scalar warp = 0;
                        for (const label fi : mesh.cells()[celli])
                        {
                            const vector nf = normalised(mesh.faceAreas()[fi]);
                            const scalar sc = Foam::sqrt(mag(mesh.faceAreas()[fi]));
                            for (const label pi : mesh.faces()[fi])
                            {
                                warp = max(warp, mag((P[pi] - mesh.faceCentres()[fi]) & nf)/sc);
                            }
                        }
                        // V_full of the cut at the far extreme (p_e) against V_k
                        scalar fullRel = -1;
                        {
                            const label st = cutChk.calcSubCell(celli, pMin + 1e-8*(pMax - pMin), n);
                            if (st == 0) fullRel = (cutChk.subCellVolume() - Vc)/Vc;
                        }
                        // the vertex levels, sorted
                        std::vector<double> lv;
                        for (const label pi : mesh.cellPoints(celli)) lv.push_back((P[pi] - xc) & n);
                        std::sort(lv.begin(), lv.end());
                        Info<< "FAILDUMP," << tag << ",cell " << celli << ",normal " << n
                            << ",aligned " << (label(&n - &normals[0]) < 26) << ",alpha " << alpha
                            << ",V_k " << Vc << ",warp/sqrt(S_f) " << warp
                            << ",(V(pMin+)-V_k)/V_k " << fullRel << ",exit " << rec2.exitCode
                            << ",pMin " << pMin << ",pMax " << pMax << ",p0 " << p0 << nl;
                        Info<< "  levels:";
                        for (double l : lv) Info<< " " << l;
                        Info<< nl << "  history (cut: p, r/V_k, A):";
                        for (label q = 0; q + 2 < hist.size(); q += 3)
                        {
                            Info<< " (" << hist[q] << ", " << hist[q+1] << ", " << hist[q+2] << ")";
                        }
                        Info<< nl << "  bracket (lo, hi) per iteration:";
                        for (label q = 0; q + 1 < histBr.size(); q += 2)
                        {
                            Info<< " (" << histBr[q] << ", " << histBr[q+1] << ")";
                        }
                        Info<< nl;
                        // V(p) on a grid of 41 points over [pMin, pMax]: monotone?
                        Info<< "  V(p)/V_k on a grid:";
                        for (label q = 0; q <= 40; ++q)
                        {
                            const scalar pq = pMin + (pMax - pMin)*q/40.0;
                            const label st = cutChk.calcSubCell(celli, pq, n);
                            const scalar Vq = (st == 0) ? cutChk.subCellVolume() : ((st < 0) ? Vc : 0);
                            Info<< " " << Vq/Vc;
                        }
                        Info<< nl;
                    }
                    if (err > 1e-10 && nDump < 3 && nDumpAll < nDumpMax)
                    {
                        ++nDump; ++nDumpAll;
                        Info<< "DUMP," << tag << "," << (arm ? "new" : "old") << ",alpha " << alpha
                            << ",cell " << celli << ",normal " << n << ",aligned "
                            << (label(&n - &normals[0]) < 26) << ",exit " << rec.exitCode
                            << ",cuts " << rec.nEvals << ",status " << status
                            << ",endState " << rec.endState << ",pMin " << pMin << ",pMax " << pMax
                            << ",pFinal " << p << ",rFinal " << rec.rFinal << ",err " << err
                            << ",errW " << errW << nl;
                    }
                }
            }
            std::sort(errs.begin(), errs.end());
            std::sort(errsW.begin(), errsW.end());
            const size_t m = errs.size();
            double sum = 0;
            size_t over = 0, overW = 0;
            for (double e : errs) { sum += e; if (e > 1e-10) ++over; }
            for (double e : errsW) { if (e > 1e-10) ++overW; }
            size_t gross = 0;
            for (double e : errs) { if (e > 1e-3) ++gross; }
            if (!errs4.empty())
            {
                std::sort(errs4.begin(), errs4.end());
                size_t small = 0;
                for (double e : errs4) if (e < 1e-6) ++small;
                Info<< "EXIT4," << tag << "," << (arm ? "new" : "old") << "," << alpha
                    << ",n " << errs4.size() << ",below1e-6 " << small
                    << ",median err " << errs4[errs4.size()/2] << ",max err " << errs4.back() << nl;
            }
            Info<< "RESULT," << tag << "," << (arm ? "new" : "old") << ","
                << alpha << "," << m << "," << errs.back() << ","
                << errs[m/2] << "," << sum/m << "," << over << ","
                << errsW.back() << "," << errsW[m/2] << "," << overW << "," << gross << ","
                << errs[size_t(0.99*(m - 1))] << ","
                << nExit[1] << "," << nExit[2] << "," << nExit[3] << "," << nExit[4] << ","
                << nNoPlane << "," << cutsSum/m << "," << cutsMax << ","
                << endOk << "," << endRej << nl;
        }
    }
    if (!fullRels.empty())
    {
        std::sort(fullRels.begin(), fullRels.end());
        size_t above = 0;
        for (double v : fullRels) if (v > 1e-13) ++above;
        Info<< "FULLREL," << tag << ",n " << fullRels.size() << ",max " << fullRels.back()
            << ",median " << fullRels[fullRels.size()/2] << ",p99 "
            << fullRels[size_t(0.99*(fullRels.size() - 1))] << ",above1e-13 " << above << nl;
    }
    return 0;
}

// ************************************************************************* //
