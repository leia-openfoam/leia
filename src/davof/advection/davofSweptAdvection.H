/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
    Modified code Copyright (C) 2018 Johan Roenby
    Modified code Copyright (C) 2019 DLR
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::davofSweptAdvection

Description
    The dual advection of the DAVOF state (the cell volume fractions alpha_c
    and the face area fractions alpha_f, reconstruction scheme davof) with
    ONE surface per interface cell: the volume flux of every face of the
    interface cells and their neighbours and the new face fraction come from
    the surface of the upwind cell (advectionScheme davofSwept; the alpha
    solver dictionary of system/fvSolution). Ported from the TwoPhaseFlow
    fork of Tomislav Maric (branch davof-swept, c5c00c5) into leia on
    2026-11-06 against OpenFOAM v2606's geometricVoF (docs/davof/PORT-DESIGN.md).

    The surface (key surface): plane (the plane in the frame of the used
    normal positioned by the volume anchor, planeSource anchor | exact),
    paraboloid (the six-coefficient paraboloid of every interface cell fitted
    to its own face fractions with the ring tensor as prior) or exactSphere
    (the exact sphere of the exactNormal dictionary, the control). The flux
    (key flux): sweptFaceCells (the semi-Lagrangian swept-face flux: the
    face's vertices and centre traced backward, the swept volume cut into the
    pieces of the cells it lies in, each wet, dry or clipped by its own
    cell's plane, the face fraction from the same clips; planes only;
    davofSweptAdvectionFlux.C, davofSweptPieces.H) or prism (the liquid
    volume of the prism between the face's pre-image and the face under the
    upwind surface, the paraboloid's flux; the volume rule of
    davofVolumeRule.H: kromerBothe | exactSphere | splitQuadrature). The
    fork's key names clippingSurface and fluxIntegration are accepted as
    aliases for one release (deprecated).

    The step (davofSweptAdvectionTemplates.C advect): the reconstruction,
    the cell surfaces (buildCellSurfaces, the volume anchor), the fluxes and
    the face update (curvedFlux), the explicit alpha update, the conservative
    bounding of isoAdvector (limitFluxes, boundFlux, extended by a second
    layer of cells; the non-conservative clipping and snapping of
    applyBruteForceBounding), the default face fractions of the faces no
    surface set (the upwind cell's new alpha). The bounding, the processor
    sync of the fluxes (syncProcPatches, with the face fractions) and the
    face-value helpers are those of OpenFOAM's isoAdvection, carried here
    as copies (OpenFOAM's class is closed), with their original copyright:

    \verbatim
        J. Roenby, H. Bredmose, H. Jasak, A computational method for sharp
        interface advection, R. Soc. Open Sci. 3 (2016) 160405,
        doi:10.1098/rsos.160405.
        H. Scheufler, J. Roenby, Accurate and efficient surface
        reconstruction from volume fraction data on general meshes,
        J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009.
    \endverbatim

    Original isoAdvection code supplied by Johan Roenby, DHI (2016),
    modified by Henning Scheufler, DLR (2019-2020). The dual advection, the
    swept-face flux, the volume rules and the anchor: Tomislav Maric, TU
    Darmstadt (2026); the citations of the algorithms stand at the functions.

SourceFiles
    davofSweptAdvection.C
    davofSweptAdvectionFlux.C
    davofSweptAdvectionSurfaces.C
    davofSweptAdvectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef davofSweptAdvection_H
#define davofSweptAdvection_H

#include "davofAdvectionScheme.H"
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "className.H"
#include "reconstructionSchemes.H"
#include "davofCurvedGeometry.H"
#include "davofVolumeRule.H"
#include "davofSweptRemote.H"
#include "Map.H"
#include "mapDistributeBase.H"
#include "interpolationCellPoint.H"
#include "pointFields.H"
#include "OFstream.H"
#include "bitSet.H"
#include "zeroField.H"
#include "boundBox.H"
#include "triFace.H"
namespace Foam { namespace davof { class velocityAt; } }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    class davof;
}



/*---------------------------------------------------------------------------*\
                        Class davofSweptAdvection Declaration
\*---------------------------------------------------------------------------*/

class davofSweptAdvection
    : public davofAdvectionScheme
{
    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;


    // Private data

        //- The reconstruction scheme (davof: the only one the dual advection
        //  accepts), constructed before everything else of this class
        autoPtr<reconstructionSchemes> surf_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- Face volumetric water transport
        surfaceScalarField dVf_;

        //- Time spent performing interface advection
        scalar advectionTime_;

        //- store timeIndex
        label timeIndex_;

        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for search of isoFace giving specified VOF value
            scalar isoFaceTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;


        // Cell and face cutting

            //- List of surface cells
            DynamicLabelList surfCells_;





        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


        // DAVOF face fractions (davof reconstruction): the moving plane's
        // submerged face area at the end of the step stored as alphaf.<group>

            //- Store the end-of-step face fractions (default: the
            //  reconstruction scheme is davof; dictionary key
            //  storeAreaFractions)
            bool storeAreaFractions_;

            //- The registered face fractions alphaf.<group>; nullptr when
            //  none is registered and storeAreaFractions_ is false
            surfaceScalarField* alphafPtr_;

            //- The face fractions created here when storeAreaFractions_
            //  and nothing was registered (initialised upwind)
            autoPtr<surfaceScalarField> alphafOwnPtr_;

            //- Faces set by a plane sweep in this step
            bitSet alphafSetFace_;

            //- The faces alphafSetFace_ marks, as a list (set by a plane
            //  or received in the sync), of this step and of the previous
            //  one; the faces whose sweep fell back to the default
            DynamicLabelList alphafSetFaces_;
            DynamicLabelList alphafPrevSetFaces_;
            DynamicLabelList alphafFallbackFaces_;

            //- The alpha the default face fractions were last evaluated
            //  with (a copy of alpha at the end of the previous step):
            //  every face no plane set holds alphaDefault_ of its upwind
            //  cell (boundary faces: its owner) as long as defaultsValid_
            bool defaultsValid_;
            scalarField alphaDefault_;

            //- The upwind direction of every internal face at the last
            //  evaluation (set: phi >= 0, the owner) and phi's event
            //  number then (phi written since: the directions are
            //  re-scanned)
            bitSet phiPos_;
            label phiEventNo_;

            //- The faces whose fraction is not 0 (kept at every write of
            //  a face value), and their list in ascending face order at
            //  the end of the step, handed to the davof scheme, which
            //  sums the area normal over them in that order
            bitSet wetFace_;
            labelList wetFaces_;

            //- The davof scheme when the reconstruction is one (nullptr
            //  otherwise)
            reconstruction::davof* davofPtr_;

            label nAlphafClipped_;


            //- Extend the conservative bounding (limitFluxes) by a second
            //  layer of cells when the face fractions are stored: the
            //  second-ring sweep sets plane fluxes on faces whose other
            //  cell lies beyond the first layer (dictionary key
            //  boundSecondLayer, default true)
            bool boundSecondLayer_;

            //- Cells beyond the first layer the extended bounding corrected
            //  in the step (summed over the bounding passes)
            label nBoundedBeyondRing_;










        // DAVOF curved clipping surface (davofSweptAdvectionSurfaces.C, 2026-10-23):
        // the dual advection with ONE surface per cell for the volume flux
        // and the face update, replacing the plane sweep when active

            //- none (the plane sweep) | exactSphere | paraboloid (key
            //  clippingSurface)
            word clippingSurface_;

            //- The face state on warped faces (key faceState, 2026-10-25):
            //  scalar (alphaf alone, the default) or vector (alphaf and
            //  the wetted area vector wettedAreaDavof.<group>, the sum over
            //  the fan triangles of the wetted area times the triangle
            //  normal, which carries the normal identity on non-planar
            //  faces; registered by the davof scheme, set here by the
            //  curved face update, the defaults alpha_up S_f and the
            //  processor sync with the sign flipped). Needs a clipping
            //  surface (the plane sweep produces no vector)
            word faceState_;
            surfaceVectorField* wetVecPtr_;

            //- The centroids the ring Monge fit of the curvature tensor
            //  reads (key tensorCentroid): cut (the alpha-positioned
            //  plane's cut-polygon centroid, the default) or anchored (the
            //  plane in the frame of the used normal re-anchored by the
            //  volume rule's cell volume, o + c e3: on a warped cell the
            //  cut polygon closes the wet part of a face by a chord, so
            //  the cut centroid is off the fan polyhedron's plane by
            //  O(warp) and the fitted tensor of an exact plane is not
            //  zero; counted as anchor-centroid failures where the root
            //  is not bracketed)
            word tensorCentroid_;
            label nCentroidFail_;

            //- The prism quadrature: relative tolerance of the flux volume
            //  and the maximum refinement depth (keys curvedTolRel,
            //  curvedMaxDepth)
            scalar curvedTolRel_;
            label curvedMaxDepth_;

            //- The six-coefficient fit: the prior weight of the ring
            //  tensor, the face-residual tolerance, the iteration cap
            //  (keys fitPriorWeight, fitTol, fitMaxIter)
            scalar fitPriorWeight_;
            scalar fitTol_;
            label fitMaxIter_;

            //- faces (the six coefficients from the face fractions, the
            //  surface interpolates them) or volume (the offset shifted
            //  after the fit so that the volume under the surface is
            //  alpha_k V, the feedback of the volume fraction into the
            //  surface; key fitAnchor)
            word fitAnchor_;

            //- six (all six coefficients from the face fractions, the
            //  Gauss-Newton), one (the offset alone from the total wetted
            //  area, the ring tensor and zero slopes: leia's paraboloidFaces)
            //  or none (the frame of the normal and the ring tensor only;
            //  the offset from the volume anchor); key fitDof
            word fitDof_;

            //- The volume anchor's tolerance (of the cell volume; key
            //  anchorTolRel) and the period of the surface-defect
            //  diagnostic in steps (key defectEvery)
            scalar anchorTolRel_;
            label defectEvery_;

            //- The quadrature diagnostic (key curvedDebug): the cell with
            //  the largest anchor residual before the anchor, printed with
            //  its pyramids against the linear tet cut
            bool curvedDebug_;

            //- The volume anchor's face-fraction diagnostics (key
            //  anchorDiagnostics, default off): the face residual of the
            //  surface before the anchor and the face shift after it, two
            //  face-fraction evaluations of every face of every interface
            //  cell per step (14 % of the advection step at N = 64,
            //  2026-10-25), printed in the Info line when on
            bool anchorDiagnostics_;

            //- The volume anchor's root finder (key anchorMethod): ccs (the
            //  default: Consecutive Cubic Spline positioning on the
            //  paraboloid, the root of the cubic Hermite interpolant through
            //  the bracket ends with the exact values and slopes dV/dc, the
            //  projected wetted patch area of the volume rule, re-bracketed
            //  and bisection-safeguarded; T. Maric, J. Comput. Phys.: X 11
            //  (2021) 100093, doi:10.1016/j.jcpx.2021.100093,
            //  arXiv:2012.02227), newton (safeguarded Newton on that slope),
            //  illinois (a bracket by doubling steps and the Illinois rule,
            //  the method until 2026-10-25) or compare (illinois and ccs,
            //  the run continues with ccs, the cells whose offsets differ
            //  counted and printed). ccs and newton take the edge of the
            //  flat root interval for a cell whose target is 0 or V (see
            //  davofSweptAdvectionSurfaces.C); illinois keeps its bracket end there.
            word anchorMethod_;

            //- anchorSkipConverged (2026-11-04, default false): with the
            //  reconstruction's positioning ccs, a plane whose cut residual
            //  |V_cut - alpha V| is already within the anchor's tolerance is
            //  kept as the reconstruction left it (no start evaluation, no
            //  anchor); counted as skipped
            bool anchorSkipConverged_;

            //- The volume flux of a face (key fluxIntegration): prism (the
            //  default, the liquid volume of the prism between the face's
            //  pre-image and the face), trapezoidal (phi_f dt times the
            //  mean of the face fractions at t^n and t^{n+1}) or simpson
            //  (with the half-step pre-image's fraction, phi_f dt
            //  (a^n + 4 a^{n+1/2} + a^{n+1})/6) or simpsonSplit (Simpson's
            //  rule on the sub-intervals between the passages of the trace
            //  over the face's vertices), simpsonSplitS (the same with every
            //  value under the surface S: the face's own fraction at
            //  sigma = 0 instead of the stored one, one-sided limits at the
            //  breakpoints) or gaussSplit (the three-point Gauss rule on
            //  those sub-intervals, no value at a breakpoint and none from
            //  the stored fraction), or sweptFace (the semi-Lagrangian
            //  swept-face flux: the face's vertices and centre traced
            //  backward, the swept volume's triangles clipped by the
            //  donor's PLANE, the flux and the new face fraction from the
            //  same clips; planes only, one donor cell; davofSweptAdvectionFlux.C,
            //  2026-10-28); see davofSweptAdvectionSurfaces.C
            word fluxIntegration_;

            //- fluxIntegration sweptFace: the tolerance of the plane test
            //  max|b_ij| h of the donor's surface (key sweptPlaneTol,
            //  default 1e-6; a surface above it is a FatalError) and the
            //  counters of the step: the faces computed, the sign flips of
            //  the triangle list, the faces whose pre-image leaves the
            //  donor, the clipped triangles, the largest closure residual
            scalar sweptPlaneTol_;
            label nSweptFaces_;
            label nSweptFlips_;
            label nSweptOutOfDonor_;
            label nSweptClipTris_;
            scalar sweptClosureMax_;

            //- The pre-image classification diagnostic (key sweptClassify,
            //  default false; davofSweptAdvectionFlux.C, 2026-10-29), for the
            //  prism and the swept flux alike: per face the class of the
            //  cells its pre-image points (vertices and centre) reach
            //  beyond the donor, 0 the donor only, 1 an interface cell (a
            //  cell with a surface), 2 a full cell, 3 an empty cell, 4 a
            //  cell outside the donor's vertex ring or outside the mesh
            //  (the priority 1 > 2 > 3 where points reach several), counted
            //  per step; the first cell that leaves [0, 1] by more than
            //  1e-10 after the conservative bounding is printed once with
            //  its faces' fluxes before and after the bounding and their
            //  classes (serial runs); the mass sum(alpha V) after the
            //  advection, the bounding and the clipping is printed every
            //  step relative to its value at t^n
            bool sweptClassify_;
            labelList sweptClass_;
            FixedList<label, 5> nSweptClass_;
            bool sweptEventReported_;
            scalarField dVfBeforeBound_;

            //- fluxIntegration sweptFaceCells (the piecewise clip,
            //  2026-10-29): the counters of the step: the pieces with area,
            //  the faces whose boundary is not covered by the donor's vertex
            //  ring (a part beyond one layer, skipped) and the largest
            //  uncovered share, the faces that fell back to the one-donor
            //  clip on a rank-deficient cap system, the partial cells
            //  without a surface met, the largest closure residual of the
            //  pieces
            label nSweptPieces_;
            label nSweptBeyond_;
            scalar sweptBeyondMax_;
            label nSweptRankFallback_;
            label nSweptPartialCells_;
            scalar sweptPieceClosureMax_;
            label nSweptDegenerate_;

            //- sweptFaceCells: the face update (key sweptCellsFaceUpdate):
            //  pieces (the default; the pre-image fan cut into the cells'
            //  pieces, each wet, dry or clipped by its cell's plane, the
            //  same clip as the flux) or donor (the pre-image fan clipped
            //  by the donor's plane alone, the face update of sweptFace;
            //  the flux stays piecewise). The ladder of 2026-10-29 found the
            //  piecewise face fractions worse than the prism's at N = 32
            //  while the volume was better: the key separates the two.
            word sweptCellsFaceUpdate_;

            //- sweptFaceCells, the face update: the rule for a part of the
            //  pre-image lying in a cell plane (key sweptCoplanar):
            //  donorSide (the default since 2026-10-29: the part belongs to
            //  the cell on the donor's side of the plane) or both (counted
            //  in both cells sharing the plane, the defect of the first
            //  ladder: every face parallel to the velocity got the sum of
            //  two fractions); the coplanar parts met in the step
            word sweptCoplanar_;
            label nSweptCoplanar_;

            //- sweptFaceCells, the pre-test of a ring cell before the clips
            //  (key sweptPreTest, 2026-11-04; davofSweptAdvectionFlux.C): none
            //  (the default: the code of 720cfc6), strict (a cell skipped
            //  when one of its planes has every vertex of the swept volume
            //  beyond 1e-10 L: bit for bit the unskipped result) or
            //  touching (also the cells the swept volume only touches;
            //  round-off in the face fractions). The ring cells skipped and
            //  clipped in the step.
            word sweptPreTest_;
            label nSweptPreSkipped_;
            label nSweptPreClipped_;

            //- sweptFaceCells, the state shortcut (key sweptStateShortcut,
            //  default false, 2026-11-04): a face whose touched ring cells
            //  are all full (alpha^n = 1, no surface) takes V_liq = |phi_f|
            //  dt and alpha_f = 1, all empty V_liq = 0 and alpha_f = 0,
            //  without a clip. The faces by path in the step, and the faces
            //  the wisp tolerance (1e-8) would take in addition.
            bool sweptStateShortcut_;

            //- sweptFaceCells, the check of the pre-test and the shortcut
            //  (key sweptPreTestCheck, default false, a debug key,
            //  2026-11-04): every skipped ring cell is clipped and a cell
            //  with a piece counted as a violation; a shortcut face runs
            //  the pieces and a sum differing by more than 1e-10 from the
            //  shortcut's is counted. The faces whose swept volume leaves
            //  the donor (no shortcut) are counted as well.
            bool sweptPreTestCheck_;
            label nSweptPreViol_;
            label nSweptPreSlab_;
            label nSweptShortViol_;
            label nSweptShortNotInside_;

            //- sweptFaceCells: the bounding boxes of the cells for the ring
            //  search of a face (2026-11-04), built once per step in
            //  curvedFlux (per face before: 26 boxes per face, hotspot 9 of
            //  the profile of 2026-11-04)
            List<boundBox> cellBoxes_;

            //- sweptFaceCells, the trace cache (key sweptTraceCache, default
            //  false, 2026-11-04): per mesh point the cell that held its
            //  midpoint and the velocity there (sweptTriangles); another
            //  face repeats only the cell search from its donor and takes
            //  the kept velocity when the same cell results (bit for bit
            //  the per-face trace); the hits and misses of the step
            bool sweptTraceCache_;
            labelList traceCell_;
            vectorField traceU_;
            bitSet traceDone_;
            label nSweptTraceHits_;
            label nSweptTraceMisses_;
            label nSweptShortFull_;
            label nSweptShortEmpty_;
            label nSweptShortPieces_;
            label nSweptShortFullTol_;
            label nSweptShortEmptyTol_;

            //- sweptFaceCells, the per-face-class statistics (key
            //  sweptClassStats, default false, 2026-11-05; profile round 3):
            //  the classes no-sweep (phi_f = 0: the face update alone, no
            //  volume clip), donor-only (a sweeping internal face whose
            //  pre-test touches the donor alone), multi-cell (two or more
            //  touched cells) and boundary (a patch face); per class and step the
            //  faces, the touched cells, the pieces, the half-space clips of
            //  the volume clip and of the face update, the triangles entered,
            //  and the time in sweptPieceVolume, sweptPieceArea, the trace
            //  (sweptTriangles) and the first pass (cell planes, state, the
            //  pre-test), by steady_clock around the calls; counters only
            bool sweptClassStats_;

            //- sweptFaceCells, the donor fast path (key sweptDonorFastPath,
            //  2026-11-05, B.1 of profile round 3): a swept volume whose
            //  vertices (x_f, the face points, the apex, the traced points)
            //  have no vertex outside any donor plane (signed distance <= 0,
            //  a vertex on a plane allowed) is clipped by the donor's PLIC
            //  plane alone, with no ring search, no pre-test of the ring and
            //  no cell-plane clips: those clips are no-ops (clipWith changes
            //  nothing without a vertex strictly outside; a triangle lying
            //  in a donor plane is kept uncut by the contact rule, the
            //  donor's centre being inside; a pre-image part lying in the
            //  face's own plane is kept by the coplanar rule on the donor's
            //  side), so the piece and the pre-image part are bit for bit
            //  the pieces path's; every other ring cell is separated from
            //  the swept volume by a donor plane (the pre-test skipped it)
            //  or, the cell across a face parallel to the velocity, adds
            //  nothing (its contact and its coplanar part are dropped). The
            //  faces taking it per step.
            bool sweptDonorFastPath_;
            label nSweptDonorFast_;

            //- sweptFaceCells, the face update from the volume clip (key
            //  sweptFaceUpdateFromVolume, default false, 2026-11-05, B.2 of
            //  profile round 3): the pre-image fan triangles are part of the
            //  swept volume's closed list, so their wetted parts from the
            //  volume clip (pieceResult::aPre, the same cell planes and PLIC
            //  plane) give the wetted area vector of a sweeping face without
            //  the separate sweptPieceArea clip. Round-off, not bit for bit:
            //  the triangles enter the clip in the opposite vertex order
            //  (outward), the crossing points are computed from the other
            //  end, and the volume clip's remnant filter (1e-14 L^2) applies.
            //  The faces without a flux, the dump and the debug face keep the
            //  separate clip, which stays selectable (the key off).
            bool sweptFaceUpdateFromVolume_;

            //- sweptFaceCells, the pre-test's coplanar exception refined (key
            //  sweptPreTestCoplanarSide, 2026-11-05, B.3 of profile round 3):
            //  a ring cell separated from the swept volume by one of its
            //  planes that holds a pre-image triangle was kept for the face
            //  update's coplanar rule; with the rule donorSide that part
            //  belongs to the cell on the donor's side of the plane, so a
            //  cell across the plane ((n . nRef) <= 0) adds nothing (its
            //  coplanar part is dropped, a contact is dropped) and is
            //  skipped: bit for bit. One third of the faces under a rotation
            //  carry no flux and touched eight cells, four of them across.
            bool sweptPreTestCoplanarSide_;

            //- Round 4 (2026-11-05), the keys default false:
            //  sweptNoFluxLayer: a face without a flux skips a ring cell
            //  across the face plane before its planes and state when the
            //  cell's own face in that plane separates it (the pre-test's
            //  test on that one plane, B.3 on: the same decision, the other
            //  five planes and the state spared); sweptAxisPlanes: the signed distance to an
            //  axis-aligned cell plane is the one coordinate difference
            //  (davofSweptPieces.H halfSpace::dist, bit for bit);
            //  sweptPointVelocityLazy: the point velocity (the
            //  volPointInterpolation of U) evaluated per point on first use,
            //  the interior points by the library's own inverse-distance
            //  weights repeated operation by operation, a point of a patch
            //  by the whole-mesh interpolation computed once (bit for bit);
            //  sweptTetCache: the trace's velocity lookup on the cell's
            //  tetrahedra cached per cell (davofVelocityAt.H, bit for bit)
            bool sweptNoFluxLayer_;
            bool sweptAxisPlanes_;
            bool sweptPointVelocityLazy_;
            bool sweptTetCache_;

            //- sweptSliverClosedForm (round 4, gain 5, default false,
            //  round-off): the piece of a touched ring cell without a PLIC
            //  plane that one plane alone cuts, by davof::sweptPieceOnePlane
            //  (davofSweptPieces.H); the pieces taken per step
            bool sweptSliverClosedForm_;
            label nSweptSliverClosed_;

            //- The deformation diagnosis (2026-11-06, leia REPORT-swept "The
            //  3D Enright deformation", all default false): sweptExactCaps,
            //  the exact piece volume by davof::sweptPieceExactVolume where
            //  the cap system is not solvable (control a, no rank fallback);
            //  sweptReinitFaceFractions, the face fractions set at the start
            //  of every step from the previous reconstruction's planes (the
            //  mean of the two clips on a face between two planes, one clip
            //  by one plane, the mean of the two alphas without a plane), the
            //  transported state discarded (controls c and d);
            //  sweptFragmentLog, the per-cell history of the cells that
            //  become face-fraction interface cells while alpha is within
            //  1e-6 of 0 or 1 (postProcessing/sweptFragmentLog/cells_<rank>.csv,
            //  the first 200 per rank; the per-step class counts in the log:
            //  a face of the cell computed by the rank fallback in the
            //  previous step or the one before, the cell in the second-layer
            //  bounding set, the cell a wisp without a neighbour surface)
            bool sweptExactCaps_;
            bool sweptReinitFaceFractions_;
            //- The per-step flux error against the exact flux of the
            //  translated sphere (key sweptExactFluxError, default false;
            //  2026-11-12, see exactFluxError)
            bool sweptExactFluxError_;
            //- Print-only diagnostics of the flux error (key
            //  sweptFluxErrorDiag, default false; 2026-11-14): the signed mean
            //  under three sign conventions, the per-cell net error e_k and
            //  the split by face orientation, one line per step after the
            //  exactFluxError line (the transport unchanged)
            bool sweptFluxErrorDiag_;
            bool sweptFragmentLog_;
            bitSet fragFallbackFaces_;
            bitSet fragFallbackFacesPrev_;
            bitSet fragBound2Cells_;
            bitSet fragBound2CellsPrev_;
            bitSet fragWispCells_;
            bitSet fragWispCellsPrev_;
            bitSet fragLogged_;
            label nFragLogged_;
            label nFragFallback_;
            label nFragBound2_;
            label nFragWisp_;
            label nFragNone_;
            autoPtr<OFstream> fragStream_;
            //- sweptExactCaps: the pieces that took the exact path per step,
            //  the pieces checked against the closure and the largest
            //  |V_exact - V_closure| over the donor's cell volume
            label nSweptExactUsed_;
            label nSweptExactChecked_;
            scalar sweptExactDiffMax_;
            //- sweptCoverageCheck (2026-11-06): per face the ring's unclipped
            //  piece volumes (no PLIC clip) summed by the closure path and by
            //  the exact path against the geometric swept volume V_g, for the
            //  faces the ring covers; the faces checked and the largest
            //  relative deviation of each path per step
            bool sweptCoverageCheck_;
            label nSweptCovChecked_;
            scalar sweptCovClosureMax_;
            scalar sweptCovExactMax_;
            //- sweptPieceDumpTol (2026-11-06, default -1 = off): with
            //  sweptExactCaps, a PLIC-clipped piece on which the closure and
            //  the exact volume differ by more than tol V_c is dumped to
            //  postProcessing/sweptPieceDump/<step>_<rank>_<n>.txt (the raw
            //  closed triangle list, the cell planes, the PLIC plane, both
            //  volumes); at most 200 per run
            scalar sweptPieceDumpTol_;
            label nSweptPieceDumps_;
            //- the first-layer bounding as a class of the fragment log: the
            //  faces whose dVf the bounding changed in the step (and the
            //  step before) with the summed |correction|, and the faces
            //  whose alphaf the bounding changed (expected none)
            bitSet fragBoundFaces_;
            bitSet fragBoundFacesPrev_;
            scalarField fragBoundMag_;
            scalarField fragBoundMagPrev_;
            label nFragAlphafChanged_;
            label nFragBoundFaces_;
            autoPtr<surfaceScalarField> fragAlphafCopy_;

            //- Decision 1 (Tomislav, 2026-11-07): the piece-volume algorithm,
            //  key sweptPieceVolume closure | closureThenTets | signedTets
            //  (oneDonorFallback: a deprecated alias of closure; the former
            //  sweptExactCaps true: closureThenTets, deprecated). closure = the
            //  cap closure for k <= 3 of full rank, the one-donor clip of the
            //  face where it does not solve (the production behaviour before
            //  2026-11-07, kept for reproducibility); closureThenTets = the
            //  signed tetrahedra (davof::sweptPieceExactVolume) where the
            //  closure does not solve, no one-donor clip anywhere;
            //  signedTets = every piece by the signed tetrahedra, no cap
            //  bookkeeping. The 2D face update is the same in all three.
            //  sweptPieceValidate: |V_closure - V_tets| on every piece where
            //  both apply (printed per step). The pieces per path per step.
            word sweptPieceVolumeWord_;
            label sweptPieceMode_;
            bool sweptPieceValidate_;
            label nSweptClosurePieces_;

            //- Decision 2 (Tomislav, 2026-11-07): the face fractions under
            //  the bounding, key sweptBoundedFaceFraction none | scale |
            //  recompute: none = the bounding changes dVf and leaves alphaf
            //  (before 2026-11-07); scale = wherever a bounding pass changes a
            //  face's dVf by the factor s = dVf_bounded / dVf, the stored face
            //  fraction (and the wetted area vector) is scaled by s, a removed
            //  flux (s = 0) dries the face (exact there, approximate for a
            //  partial change); recompute = the face fraction recomputed as the
            //  wetted area of the pre-image fan under the donor's plane
            //  positioned along its normal so that the swept liquid volume
            //  (the signed tetrahedra under that plane alone) equals
            //  |dVf_bounded| (davof::solveOffset), consistent by construction,
            //  scale where no plane or no bracket. Both sides of a processor
            //  face apply the same s. The counters per step; for recompute the
            //  closed swept list and the donor plane of every computed face
            //  are kept for the step.
            word sweptBoundedFaceFraction_;
            label nBoundedFaces_;
            label nBoundedZero_;
            label nBoundedPartial_;
            scalar maxOneMinusS_;
            label nRecomputed_;
            label nRecomputeNoPlane_;
            label nRecomputeNotBracketed_;
            label recomputeEvals_;
            label recomputeEvalsMax_;
            List<List<FixedList<point, 3>>> sweptFaceTris_;
            vectorField sweptFacePlaneN_;
            pointField sweptFacePlaneX0_;
            bitSet sweptFaceHasPlane_;
            mutable vectorField upLazy_;
            mutable bitSet upDone_;
            mutable bitSet anyPatchPoint_;
            mutable DynamicList<scalar> upWeights_;
            mutable label nUpLazy_;
            mutable List<List<triFace>> cellTets_;

            //- The point velocity of point pi (sweptPointVelocityLazy)
            void ensurePointVelocity(const label pi) const;

            //- Every point velocity by the whole-mesh interpolation (a patch
            //  point in serial; every rank up front in parallel: collective)
            void ensureAllPointVelocities() const;
            FixedList<label, 4> nClsFaces_;
            FixedList<label, 4> nClsTouched_;
            FixedList<label, 4> nClsPieces_;
            FixedList<label, 4> nClsVolClips_;
            FixedList<label, 4> nClsAreaClips_;
            FixedList<label, 4> nClsTris_;
            FixedList<scalar, 4> tClsVol_;
            FixedList<scalar, 4> tClsArea_;
            FixedList<scalar, 4> tClsTrace_;
            FixedList<scalar, 4> tClsPre_;

            //- sweptFaceCells, the cap set of a piece (key sweptCapSet,
            //  2026-11-02): final (the default: the planes that produced an
            //  edge of the final clipped polygons, a triangle lying in a
            //  cell plane with the donor outside that cell a contact; the
            //  rank fallback does not depend on the order of a cell's faces)
            //  or crossing (the rule before 2026-11-02: a plane counted when
            //  it crossed the polygon at its turn in the clipping order;
            //  kept for the comparison); davofSweptPieces.H
            word sweptCapSet_;

            //- The consistency of the swept flux under oblique travel
            //  (2026-11-02, davofSweptAdvectionFlux.C sweptConsistencyProbe).
            //  A region of a cell can lie in the swept volumes of two faces
            //  (under north-east travel a corner of C in those of C|N and
            //  N|NE); the cell that receives it balances only if both faces
            //  give that region the same liquid. Two ways break it: a face
            //  outside the candidate set keeps the algebraic upwind flux
            //  alpha_donor phi_f dt while its swept volume reaches a cell of
            //  another state, and a partial cell without a surface (a wisp)
            //  is wet or dry by alpha >= 1/2 in the pieces while the
            //  algebraic faces release alpha phi_f dt.
            //
            //  sweptConsistencyLog (default false; serial runs only): per
            //  step (i) the algebraic faces whose swept volume reaches a cell
            //  other than the donor in another state (interface, wisp, the
            //  opposite pure state) with the liquid mismatch of the consistent
            //  rule minus the algebraic flux, (ii) the wisp pieces of the
            //  computed faces with the mismatch of the 0/1 rule against
            //  alpha V_piece, (iii) per cell the net residue the mismatches
            //  leave; postProcessing/sweptConsistency/summary.dat, the
            //  per-step residues residue_<step>.csv and the field
            //  sweptResidue (cumulative, in units of the cell volume).
            //
            //  sweptConsistent (default false; sweptFaceCells): every face
            //  whose donor shares a vertex with an active cell (a cell with a
            //  surface, an interface cell, or a cell with sweptConsistentTol <
            //  alpha < 1 - sweptConsistentTol) is computed by the swept flux,
            //  also without a surface on either side, and a piece in a cell
            //  without a surface carries alpha_c V_piece (its area alpha_c
            //  a_piece in the face update), the value every face that sweeps
            //  it, computed or algebraic, gives it.
            bool sweptConsistencyLog_;
            bool sweptConsistent_;
            scalar sweptConsistentTol_;
            label nConsProbed_;
            label nConsReach_;
            label nConsReachInterface_;
            label nConsReachWisp_;
            label nConsReachOpposite_;
            label nConsProbeFallback_;
            scalar consAlgMisSum_;
            scalar consAlgMisAbs_;
            scalar consAlgMisMax_;
            label nConsWispPieces_;
            scalar consWispMisSum_;
            scalar consWispMisAbs_;
            scalar consWispMisMax_;
            //- Per computed face: its liquid volume minus the consistent
            //  rule's (zero unless the face met a cell without a surface)
            scalar sweptFaceConsDelta_;
            //- sweptConsistent: the candidates added, the faces computed
            //  without a surface, their rank fallbacks (the algebraic value)
            label nConsCandAdded_;
            label nConsNoSurface_;
            label nConsFallback_;
            //- sweptConsistent: the boundary faces without a surface left to
            //  the boundary condition's flux, as on the default path
            //  (2026-11-03)
            label nConsBoundaryKept_;
            //- The face computed without a surface (sweptConsistent)
            bool sweptNoDonorSurface_;

            //- The surface of a wisp's pieces in sweptFaceCells (key
            //  sweptWispSurface, 2026-11-02): none (a cell without a surface
            //  and with 1e-8 < alpha^n < 1 - 1e-8 takes its pieces wet for
            //  alpha^n >= 1/2 and dry otherwise) | neighbour (the plane of
            //  its closest face neighbour with a surface clips its pieces and
            //  its part of the face update: the second-ring rule of a face
            //  without a surface, applied to a ring cell). A plane cutting a
            //  cell in a small corner leaves |m_c| below the interface-cell
            //  threshold, so the cell has no plane, while its corner holds
            //  liquid: the 0/1 rule dropped it (the plane gate, 2026-11-02)
            //  | anchored (2026-11-03: the neighbour's normal, the plane
            //  positioned in the wisp so that the volume under it is
            //  alpha^n V, the PLIC of the wisp; neighbour, the neighbour's
            //  plane itself, gave a wisp of a curved interface more liquid
            //  than it holds, min alpha -1.7e-5 on the N = 64 rotation)
            word sweptWispSurface_;
            //- The band of the wisp class (key sweptWispTol, default 1e-8,
            //  the classification before the key): a ring cell without a
            //  surface is full for alpha^n >= 1 - tol, empty for alpha^n <=
            //  tol, a wisp between; with sweptWispSurface neighbour the wisps
            //  take the neighbour's plane. A plane cutting a corner of volume
            //  below 1e-8 V of a cell left that cell empty and the face
            //  fraction of a face of it short by 1e-5 (the oblique arm of the
            //  plane gate, 2026-11-02)
            scalar sweptWispTol_;
            List<davof::cellSurface> wispSurfaces_;
            bitSet hasWispSurface_;
            label nSweptWispClipped_;
            label nSweptWispNoNeighbour_;
            label nSweptWispAnchorFail_;
            //- matched (2026-11-03): the neighbour's plane where its volume in
            //  the wisp is alpha^n V within sweptWispMatchTol V (key, default
            //  1e-11, the anchor's tolerance anchorTolRel), else anchored
            label nSweptWispMatched_;
            scalar sweptWispMatchTol_;
            scalarField consResidueStep_;
            autoPtr<OFstream> consSummary_;
            autoPtr<volScalarField> consResiduePtr_;

            //- clippingSurface plane: where the plane comes from (key
            //  planeSource, 2026-11-01): anchor (the default: the frame of
            //  the used normal positioned by the volume anchor) or exact
            //  (the tangent plane of the exact sphere of the exactNormal
            //  dictionary at the point nearest to the cell centre: a
            //  diagnostic control without the anchor and without feedback
            //  from the transported state)
            word planeSource_;

            //- The one-step face dump of sweptFaceCells (key sweptFaceDump,
            //  default false): postProcessing/sweptFaceDump/faces.csv and
            //  pieces.csv at the first step (serial runs)
            bool sweptFaceDump_;
            autoPtr<OFstream> dumpFaces_;
            autoPtr<OFstream> dumpPieces_;

            //- The kinematic consistency of the swept fluxes (key
            //  sweptScaleToPhi, default false, 2026-10-29): the liquid volume
            //  of the swept volume scaled by |phi_f| dt / V_swept, so that a
            //  wholly wet swept volume gives the flux phi_f dt that the
            //  upwind flux of a full cell uses. Under a rotation the traced
            //  swept volume (straight edges between traced points) differs
            //  from phi_f dt at the order of the angular displacement per
            //  step squared, and a full cell next to an interface cell
            //  received 3e-7 more than it sent in the first step of the
            //  N = 64 rotation (the prism's volume d.S_f IS phi_f dt). The
            //  identity to round-off under a uniform translation.
            //  Since 2026-11-02 it is the mode scale of
            //  sweptVolumeCorrection below (true when that mode is in use).
            bool sweptScaleToPhi_;

            //- The volume correction of the swept volume (key
            //  sweptVolumeCorrection, 2026-11-02): none | scale (the liquid
            //  volume times V_f / V_g, sweptScaleToPhi; kept for comparison)
            //  | pyramid (the apex of the pre-image fan moved along the
            //  pre-image's area vector so that the closed triangle list holds
            //  V_f, the volume correction of the unsplit face-vertex flux
            //  calculation, UFVFC: T. Maric, H. Marschall, D. Bothe, An
            //  enhanced un-split face-vertex flux-based VoF method,
            //  J. Comput. Phys. 371 (2018) 967-993,
            //  doi:10.1016/j.jcp.2018.03.048, no arXiv version found; the
            //  acronym as in T. Maric, D. B. Kothe, D. Bothe, Unstructured
            //  un-split geometrical Volume-of-Fluid methods - A review,
            //  J. Comput. Phys. 420 (2020) 109695,
            //  doi:10.1016/j.jcp.2020.109695, arXiv:2007.02860). Without the
            //  key: scale where sweptScaleToPhi is true, none otherwise (the
            //  builds before it).
            word sweptVolumeCorrection_;
            bool sweptVolumeKeyGiven_;
            //- Per step: the largest |delta|/h of the pyramid, the largest
            //  |V_f - V_g|/|V_f| before the correction, the largest per-axis
            //  displacement of a traced vertex over h, the twisted faces
            //  (V_g below 0.1 of the unsigned swept volume), the faces whose
            //  pre-image has no area (no pyramid, counted)
            scalar sweptPyrDeltaMax_;
            scalar sweptVolDefectMax_;
            scalar sweptDispAxisMax_;
            label nSweptTwisted_;
            label nSweptPyrSkipped_;

            //- The swept fluxes in parallel (davofSweptRemote.H,
            //  2026-10-30): the cells of other ranks in the vertex rings of
            //  this step's donors (zoneDistribute's cell-point-cell stencil,
            //  fetched once per step in sweptExchange), per donor the indices
            //  of its remote ring cells (ascending global index); empty in
            //  serial, where the path is the local one unchanged
            List<davof::remoteCell> sweptRemote_;
            Map<labelList> sweptRemoteRing_;
            //- Per step: the remote cells fetched, the faces whose pieces
            //  or trace used a remote cell, the traces that fell back to
            //  the Euler step (no cell, local or remote, holds the
            //  midpoint), the wall times of the exchange, of curvedFlux and
            //  of the advect call (seconds, this rank)
            label nSweptRemoteCells_;
            label nSweptRemoteFaces_;
            label nSweptEulerTrace_;
            bool sweptFaceRemote_;
            //- The face took the rank fallback (the one-donor clip; the
            //  column fallback of sweptFaceLog, 2026-11-02)
            bool sweptFaceFallback_;
            scalar sweptExchangeTime_;
            scalar sweptCurvedTime_;
            //- The parts of the exchange: the wait at a barrier for the
            //  slowest rank, the stencil update, the map, the packing, the
            //  transfer, the unpacking (seconds, this rank)
            FixedList<scalar, 6> sweptExchangeParts_;
            //- The cells whose stencil zoneCPCStencil holds (the update is
            //  called only for new donors), and the map of the transfer,
            //  kept while no rank's set of remote cells changes
            bitSet sweptStencilDone_;
            autoPtr<mapDistributeBase> sweptMap_;
            labelList sweptPrevNeed_;
            labelList sweptPrevCompact_;
            //- The per-face log of the swept fluxes (key sweptFaceLog,
            //  default false): per step and rank the faces computed with the
            //  face centre, the face fraction, the flux in the face's
            //  canonical orientation and the remote flag, in
            //  postProcessing/sweptFaceLog/<step>_<rank>.csv (the
            //  decomposition comparison)
            bool sweptFaceLog_;
            //- The debug face of sweptFaceCells (key sweptDebugPoint, a
            //  point; absent: off): the face whose centre lies within 1e-9
            //  of it writes its pieces every step
            bool sweptDebug_;
            point sweptDebugPoint_;
            autoPtr<OFstream> sweptFaceLogFile_;
            scalar debugR0Max_;
            label debugCellR0_;
            davof::cellSurface debugSurface_;
            scalar debugAlpha_;
            label nAnchorFail_;
            label nAnchorLeaves_;

            //- The cell-volume evaluations of the anchor in the step
            label nAnchorEvals_;

            //- The cells anchored at the edge of the flat interval (a
            //  target of 0 or V) in the step
            label nAnchorEdge_;

            //- The anchor's start (2026-11-04): the cells with a start
            //  evaluation, the maximum and the sum of |r_start|/V, the cells
            //  whose start residual is already within the tolerance, the
            //  cells skipped (anchorSkipConverged), and the largest
            //  disagreement |V_KromerBothe - V_cut|/V between the volume
            //  rule and the reconstruction's cut of the same plane
            //  (positioning ccs; -1 otherwise)
            label nAnchorStartCells_;
            scalar anchorStartResMax_;
            scalar anchorStartResSum_;
            label nAnchorStartBelowTol_;
            label nAnchorSkipped_;
            scalar anchorKbVsCutMax_;

            //- anchorMethod compare: the cells printed in the step, the
            //  cells whose offsets differ, those of them at the edge
            label nCompare_;
            label nCompareDiff_;
            label nCompareDiffEdge_;
            label nCompareDiffResolved_;
            scalar anchorShiftMax_;
            scalar anchorShiftSum_;

            //- The refinement levels of the cell volume under the surface
            //  for the surface defect (key defectLevels; -1 switches the
            //  diagnostic off)
            label defectLevels_;

            //- Write the per-cell fit residual and surface defect fields
            //  (key writeCurvedFields)
            bool writeCurvedFields_;

            //- The surface of every cell that has one (interface cells)
            List<davof::cellSurface> cellSurfaces_;
            bitSet hasSurface_;

            //- The run-time selected rule for the cell volume, the prism
            //  volume and the face fraction under a surface (key
            //  volumeRule: splitQuadrature | kromerBothe | exactSphere)
            autoPtr<davofVolumeRule> volumeRulePtr_;
            label nRuleFallback_;

            //- Counters of the step
            label nCurvedFaces_;
            label nCurvedSecondRing_;
            label nCurvedBoundary_;
            label nCurvedSkipped_;
            label nCurvedFallback_;
            label nCurvedLeaves_;
            label nFitCells_;
            label nFitConverged_;
            label nFitStalled_;
            label nFitSingular_;
            label nFitNoBracket_;
            label nFitNoFrame_;
            scalar fitResMax_;
            scalar fitResSum_;
            scalar defectMax_;
            scalar defectSum_;
            label nDefectCells_;
            autoPtr<volScalarField> fitResidualPtr_;
            autoPtr<volScalarField> paraboloidDefectPtr_;


    // Private Member Functions

        //- No copy construct
        davofSweptAdvection(const davofSweptAdvection&) = delete;

        //- No copy assignment
        void operator=(const davofSweptAdvection&) = delete;


        // Advection functions

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(bitSet& markedCell) const;

            void setProcessorPatches();


            //- For a given cell return labels of faces fluxing out of this cell
            //  (based on sign of phi)
            void setDownwindFaces
            (
                const label celli,
                DynamicLabelList& downwindFaces
            ) const;

            // LimitFluxes
            template < class SpType, class SuType >
            void limitFluxes
            (
                const SpType& Sp,
                const SuType& Su
            );

            // Bound fluxes
            template < class SpType, class SuType >
            void boundFlux
            (
                const bitSet& nextToInterface,
                surfaceScalarField& dVfcorrectionValues,
                DynamicLabelList& correctedFaces,
                const SpType& Sp,
                const SuType& Su
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label celli
            ) const;


            //- Clear the surface-cell list of the step
            void clearIsoFaceData()
            {
                surfCells_.clear();
            }

        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value;
            //  alphaf (optional) is exchanged on the same faces, the upwind
            //  side's face fraction copied without the sign change
            DynamicList<label> syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi,
                bool returnSyncedFaces=false,
                surfaceScalarField* alphaf = nullptr,
                surfaceVectorField* wetVec = nullptr
            );


        // DAVOF face fractions

            //- Look up alphaf.<group> in the registry; create it (upwind)
            //  when storeAreaFractions_ and none is registered
            void initAreaFractions();

            //- The default of one face from the alpha field given: the
            //  upwind cell's value (phi >= 0: the owner), boundary faces
            //  the owner's; faces of empty patches carry no value
            void setDefaultAreaFraction
            (
                const label facei,
                const scalarField& alpha
            );

            //- After the bounding: every face not set by a plane in this
            //  step holds the upwind cell's new alpha. Evaluated on the
            //  faces whose stored value can differ from it only (the
            //  faces set in the previous step, the fallbacks of this one,
            //  the faces of the cells whose alpha changed bitwise against
            //  alphaDefault_, the faces whose upwind direction changed);
            //  on every face when defaultsValid_ is false
            void defaultAreaFractions();

            //- The list of the wet faces from wetFace_ (bitSet::toc(),
            //  ascending), handed to the davof scheme
            void publishWetFaces();






            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label facei);

            //- Apply the bounding based on user inputs
            void applyBruteForceBounding();

public:

        // DAVOF curved clipping surface (davofSweptAdvectionSurfaces.C)

            //- The surface of every interface cell for the state at tState
            //  (exactSphere: the translated sphere; paraboloid: the
            //  six-coefficient fit to the cell's face fractions)
            void buildCellSurfaces(const scalar tState);

            //- The volume flux and the face fraction of every face of the
            //  interface cells and their neighbours from the surface of the
            //  upwind cell (second ring: the closest neighbouring surface)
            void curvedFlux(const scalar dt);

            //- The deformation diagnosis (2026-11-06): the face fractions
            //  from the previous reconstruction's planes
            //  (sweptReinitFaceFractions) and the per-cell fragment log
            //  (sweptFragmentLog), both called from advect
            void reinitFaceFractions();
            void fragmentLog();

            //- The flux error of the step on the translated sphere
            //  (sweptExactFluxError, 2026-11-12): for every face of the
            //  interface cells, the computed volume flux dVf (after the flux,
            //  before the bounding) against the exact one, the liquid volume
            //  of the prism between f - U0 dt and f under the sphere of the
            //  exactNormal dictionary at t^n (the exact sphere-tetrahedron
            //  overlap, davofVolumeRules::exactSphere; S. Strobl, A. Formella,
            //  T. Poeschel, J. Comput. Phys. 311 (2016) 158-172,
            //  doi:10.1016/j.jcp.2016.02.003), exact for a uniform
            //  translation; |dVf - dVf_exact|/V of the upwind cell as the
            //  mean (L1) and maximum over the faces, one line per step.
            //  Skipped (one message) when the dictionary carries a rotation.
            void exactFluxError(const scalar dt);

            //- Decision 2: the face fraction of a face whose dVf the bounding
            //  changed from dVfOld to dVfNew (scale or recompute)
            void boundedFaceFraction
            (
                const label facei,
                const scalar dVfOld,
                const scalar dVfNew
            );

            //- tensorCentroid anchored: for every interface cell the plane
            //  in the frame of the used normal through the cut centroid
            //  re-anchored by the volume rule's cell volume (the fan
            //  polyhedron) to alpha_k V, the point o + c e3 handed to the
            //  scheme as the centroid of the ring Monge fit
            void anchorTensorCentroids();

            //- The Info line of the step
            void reportCurved() const;

            //- fluxIntegration sweptFace (davofSweptAdvectionFlux.C): the liquid
            //  volume vLiq of the volume swept by facei between its traced
            //  pre-image and itself under the PLANE S of the cell donor
            //  (positive; the caller gives it the sign of phiF), the
            //  pre-image's face fraction af and wetted area vector aWet in
            //  the face's orientation. sweeps false (the displacement
            //  parallel to the face): the fraction only, no volume, no
            //  orientation count. UInterp interpolates U in a cell, Up is
            //  the volPointInterpolation of U at the vertices.
            void sweptFaceFlux
            (
                const label facei,
                const label donor,
                const davof::cellSurface& S,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const davof::velocityAt& UInterp,
                const vectorField& Up,
                scalar& af,
                vector& aWet,
                scalar& vLiq,
                pointField* qOut = nullptr,
                point* xqOut = nullptr
            );

            //- fluxIntegration sweptFaceCells (davofSweptAdvectionFlux.C,
            //  davofSweptPieces.H): the same arguments and results as
            //  sweptFaceFlux, the swept volume cut into the pieces inside
            //  the donor and its vertex-ring cells, each wet in total (a
            //  full cell), dry (an empty cell) or clipped by that cell's
            //  own plane (an interface cell); the face fraction from the
            //  same clip of the pre-image fan
            void sweptFaceCellsFlux
            (
                const label facei,
                const label donor,
                const davof::cellSurface& S,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const davof::velocityAt& UInterp,
                const vectorField& Up,
                scalar& af,
                vector& aWet,
                scalar& vLiq,
                pointField* qOut = nullptr,
                point* xqOut = nullptr
            );

            //- The traced pre-image (q, xq) and the outward triangle list T
            //  of the swept volume of facei (shared by the two swept
            //  fluxes); preSign is the orientation of the stored pre-image
            //  triangles against the face's own order; counts the flips,
            //  the closure residual and the out-of-donor faces
            void sweptTriangles
            (
                const label facei,
                const label donor,
                const scalar dt,
                const scalar phiF,
                const bool sweeps,
                const davof::velocityAt& UInterp,
                const vectorField& Up,
                List<FixedList<point, 3>>& T,
                pointField& q,
                point& xq,
                label& preSign,
                scalar& vTotal
            );

            //- The swept fluxes in parallel: the remote cells of the vertex
            //  rings of the cells marked in donors (zoneDistribute's stencil,
            //  one mapDistributeBase transfer of the packed cells,
            //  davofSweptRemote.H) into sweptRemote_ and sweptRemoteRing_;
            //  collective, a no-op in serial
            void sweptExchange
            (
                const boolList& donors,
                const vectorField& Up
            );

            //- sweptConsistencyLog: the algebraic faces of the step (the
            //  internal faces whose flux curvedFlux did not set) probed for
            //  the consistent rule's liquid, the residues of the step and
            //  the summary row (serial)
            void sweptConsistencyProbe
            (
                const scalar dt,
                const davof::velocityAt& UInterp,
                const vectorField& Up
            );

            //- The packed cell of davofSweptRemote.H (the owner's side)
            void packSweptCell
            (
                const label celli,
                const vectorField& Up,
                DynamicList<scalar>& v
            ) const;

            //- One row of the face dump (sweptFaceDump)
            void writeFaceDump
            (
                const label facei,
                const label donor,
                const scalar phiF,
                const bool sweeps,
                const scalar afPieces,
                const scalar afDonor,
                const scalar afUsed,
                const bool rankFallback,
                const vector& Apre,
                const label nPieces
            );

            //- The class of the cells the pre-image points q (vertices)
            //  and xq (centre) of a face of the cell donor reach beyond
            //  the donor (sweptClassify_; the codes above)
            label classifyPreImage
            (
                const label donor,
                const pointField& q,
                const point& xq
            ) const;


public:

    //- Runtime type information
    TypeName("davofSwept");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        davofSweptAdvection
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~davofSweptAdvection() = default;


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        template < class SpType, class SuType >
        void advect(const SpType& Sp, const SuType& Su);

        //- advection of the interface
        virtual void advect()
        {
            zeroField Sp;
            zeroField Su;
            advect<zeroField,zeroField>(Sp,Su);
        }

        //- advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        )
        {
            advect<volScalarField::Internal,volScalarField::Internal>(Sp,Su);
        }

        //- The reconstruction of the scheme
        virtual reconstructionSchemes& surf()
        {
            return surf_();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "davofSweptAdvectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
