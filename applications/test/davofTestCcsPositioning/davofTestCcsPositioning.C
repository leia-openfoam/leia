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

    cutCellPLIC cutPos(mesh);
    cutCellPLIC cutChk(mesh);
    const scalarField& V = mesh.V();
    const pointField& P = mesh.points();

    Info<< "RESULT,header,tag,arm,alpha,nCases,errMax,errMedian,errMean,nOver1e-10,"
        << "exit1,exit2,exit3,exit4,nNoPlane,cutsMean,cutsMax,endOk,endRej" << nl;

    for (label arm = 0; arm < 2; ++arm)
    {
        const bool endInterval = (arm == 1);
        const bool relativeTol = (arm == 1);
        const scalar tolRel = (arm == 1) ? 1e-12 : 1e-14;
        for (const scalar alpha : alphas)
        {
            std::vector<double> errs;
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

                    // the independent volume error
                    double err;
                    if (alpha <= 0.5)
                    {
                        const label st = cutChk.calcSubCell(celli, p, n);
                        const scalar Vcut = (st == 0) ? cutChk.subCellVolume() : ((st < 0) ? Vc : 0);
                        err = mag(Vcut - alpha*Vc)/(alpha*Vc);
                    }
                    else
                    {
                        const label st = cutChk.calcSubCell(celli, -p, -n);
                        const scalar W = (st == 0) ? cutChk.subCellVolume() : ((st < 0) ? Vc : 0);
                        err = mag(W - (1 - alpha)*Vc)/((1 - alpha)*Vc);
                    }
                    errs.push_back(err);
                }
            }
            std::sort(errs.begin(), errs.end());
            const size_t m = errs.size();
            double sum = 0;
            size_t over = 0;
            for (double e : errs) { sum += e; if (e > 1e-10) ++over; }
            Info<< "RESULT," << tag << "," << (arm ? "new" : "old") << ","
                << alpha << "," << m << "," << errs.back() << ","
                << errs[m/2] << "," << sum/m << "," << over << ","
                << nExit[1] << "," << nExit[2] << "," << nExit[3] << "," << nExit[4] << ","
                << nNoPlane << "," << cutsSum/m << "," << cutsMax << ","
                << endOk << "," << endRej << nl;
        }
    }
    return 0;
}

// ************************************************************************* //
