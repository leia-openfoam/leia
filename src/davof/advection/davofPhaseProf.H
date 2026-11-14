/*---------------------------------------------------------------------------*\
Description
    Print-only wall-clock timers and counters of the swept-face pieces path
    (diag2, 2026-11-14; docs/davof/REPORT-swept.md "Pieces path: aborts and
    cost"). Behind the dictionary key phaseProfile (default false): with the
    key off every timer is a branch on a bool and nothing is printed. The
    timers read the steady clock only; no result of the method depends on
    them. One accumulator per rank and per step, reduced and printed by
    davofSweptAdvection::advect() (mean and maximum over the ranks).
\*---------------------------------------------------------------------------*/

#ifndef davofPhaseProf_H
#define davofPhaseProf_H

#include <array>
#include <chrono>
#include <vector>

namespace Foam
{
namespace davofProf
{

// wall-clock phases (seconds per step)
enum Phase : int
{
    // advect()
    pBetween,       // from the end of the previous advect() to this one's start:
                    // function objects, IO, Courant numbers, the solver's own work
    pReinit,        // reinitFaceFractions
    pRecon,         // surf_->reconstruct() in total
    pCurved,        // curvedFlux() in total
    pAlphaUpdate,   // alpha from the fluxes, mass sums (reductions)
    pLimit,         // limitFluxes (boundFlux, the faces' scale factors, syncs)
    pBrute,         // applyBruteForceBounding
    pDefaultAf,     // defaultAreaFractions, publishWetFaces
    pReport,        // reportCurved (the reductions of the counters and the print)
    pAdvect,        // advect() in total
    // reconstruct()
    pRcAreaNormal,  // start of the function, updateHasValue, m_k = sum alpha_f S_f
    pRcCells,       // the interface-cell rule (selection of the interface cells)
    pRcExplicit,   // face moment and the explicit-plane cut of every interface cell
    pRcSelect,      // the normal model (ring normal)
    pRcPositioning, // positionPlane (CCS) of every interface cell
    pRcTail,        // boundary updates, reductions and the report line
    // curvedFlux()
    pCfInit,        // counters, cell boxes, until buildCellSurfaces
    pCfSurfaces,    // buildCellSurfaces
    pCfVelSetup,    // point velocities, tet cache, octree
    pCfExch1,       // exchange of the cell surfaces across the processor patches
    pCfWisp,        // the surfaces of the wisps
    pCfCand,        // the candidate cells
    pCfExch2,       // exchange of the second ring
    pCfSweptExch,   // sweptExchange (the remote ring cells)
    pCfFaceLoop,    // the loops over the faces (processFace) in total
    pCfProcSync,    // syncProcPatches of the fluxes
    pCfTail,        // the surface defect diagnostic
    // per face (sums over the faces of the rank)
    pFFlux,         // sweptFaceCellsFlux in total
    pFTrace,        // sweptTriangles: the trace of the face = the signed tets
    pFPre,          // the ring, the planes and the pre-test
    pFVol,          // sweptPieceVolume: the clips of the tets by the cell planes
    pFArea,         // sweptPieceArea: the clips of the pre-image
    nPhase
};

static const char* const phaseName[nPhase] =
{
    "between", "reinit", "recon", "curved", "alphaUpdate", "limit", "brute", "defaultAf", "report", "advect",
    "rcAreaNormal", "rcCells", "rcExplicit", "rcSelect", "rcPositioning", "rcTail",
    "cfInit", "cfSurfaces", "cfVelSetup", "cfExch1", "cfWisp", "cfCand", "cfExch2", "cfSweptExch", "cfFaceLoop",
    "cfProcSync", "cfTail",
    "fFlux", "fTrace", "fPre", "fVol", "fArea"
};

// counts of the step (sums over the ranks)
enum Count : int
{
    cInterface,     // interface cells (set by the reconstruction)
    cExplicitCuts,  // cuts of the explicit planes (one per interface cell)
    cCcsCuts,       // cuts of the CCS positioning
    cFaces,         // faces processed by the flux (processFace)
    cFacesFlux,     // faces through sweptFaceCellsFlux
    cFacesSweep,    // faces with a swept volume
    cTets,          // triangles of the swept volumes (the signed tets)
    cRing,          // ring cells of the faces
    cTouched,       // ring cells the pre-test keeps (the cells whose planes clip the swept volume)
    cVolClips,      // plane clips of the swept volumes
    cAreaClips,     // plane clips of the pre-images
    cPieces,        // pieces produced
    cExch,          // exchange rounds issued by the davof code (explicit)
    cLimitCells,    // cells of the limiter's bounding set
    cTetsBuilt,     // signed tets (x_ref, T) built and clipped by the pieces (one per triangle and ring cell)
    cTetClips,      // tet-plane clips: one convex polyhedron clipped by one half-space
    cPolyClips,     // polygon clips of the tets' faces (the Sutherland-Hodgman steps of the tet-plane clips)
    cCovClips,      // clips of the triangle's own coverage polygon (the 2D clips of sweptPieceExactVolume)
    nCount
};

static const char* const countName[nCount] =
{
    "interface", "explicitCuts", "ccsCuts", "faces", "facesFlux", "facesSweep", "tets", "ring", "touched",
    "volClips", "areaClips", "pieces", "exch", "limitCells", "tetsBuilt", "tetClips", "polyClips", "covClips"
};

inline bool on = false;                 // phaseProfile
inline bool dump = false;               // alphaDump
inline double dumpRecord = 1e-8;        // the largest max(alpha) - 1 dumped so far
inline double lastEnd = -1;             // steady-clock seconds at the end of the previous advect()
inline std::array<double, nPhase> T{};
inline std::array<double, nCount> N{};
// alphaDump: per face of the step, the swept volume, the unscaled liquid volume, the flux scale
inline std::vector<double> faceVSwept, faceVLiq0, faceScale;

inline double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

//- A start time; 0 when the key is off
inline double start() { return on ? now() : 0.0; }

//- Add the time since t to phase k and restart t
inline void lap(double& t, const int k)
{
    if (!on) return;
    const double n = now();
    T[k] += n - t;
    t = n;
}

inline void add(const int k, const double dt) { if (on) T[k] += dt; }
inline void count(const int k, const double n = 1) { if (on) N[k] += n; }

inline void reset()
{
    T.fill(0.0);
    N.fill(0.0);
}

} // End namespace davofProf
} // End namespace Foam

#endif
