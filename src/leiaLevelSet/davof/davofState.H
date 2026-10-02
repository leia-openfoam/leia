/*---------------------------------------------------------------------------*\
    leia :: davof :: davofState

    The DAVOF (dual area/volume-of-fluid) state of a mesh: the cell volume
    fractions alpha_c, the face liquid-area fractions alpha_f (one value per
    face, all faces of the polyMesh including empty and coupled patches), and
    the Gauss-identity area normal recovered from them,

        m_c = - sum_f alpha_f S_f^out,    S_f^out = +S_f (owner), -S_f (neighbour),

    whose direction is the interface normal pointing out of the liquid and
    whose magnitude is the interface area in the cell.

    Two sources of (alpha_c, alpha_f):

      detrixheAslam       The exact signed distance of an implicit surface at
                          the cell centre, the face centres and the mesh points
                          in the Detrixhe-Aslam formulas: the tets (xc, xf, p0,
                          p1) of every cell for the volume, the triangle fans
                          (xf, p0, p1) with the SAME vertex values for the face
                          area. The formulas are exact for a linear function on
                          each simplex, so the state describes ONE piecewise-
                          linear C0 surface (the chord surface) and the identity
                          above is exact for it; the chords are O(h^2) off a
                          curved surface, which leaves the face fractions O(h)
                          as fractions and the normal first order (measured).
                          Alias: linearInterpolant (the name of the first runs).
      quadraticFaces      detrixheAslam for the cell fractions, and for the face
                          fractions the quadratic interpolant of the exact
                          signed distance at the vertices AND the edge midpoints
                          of every fan triangle, cut exactly
                          (davofQuadraticFaceGeometry.H): crossings O(h^3),
                          fractions O(h^2), a second-order normal from exact
                          signed distances alone, non-iteratively, for any
                          implicit surface. Triangles the routine cannot handle
                          keep the Detrixhe-Aslam value (nFaceFallback()).
      planePhaseIndicator leia's production construction: one plane per
                          narrow-band cell (least-squares fit of psi, or the
                          exact tangent plane at the cell centre when
                          levelSet.phaseIndicator.geometrySource is
                          analyticImplicitSurface), cell fractions as in
                          detrixheAslamPhaseIndicator, face fractions as the
                          central average of the owner and neighbour planes as
                          in computeFaceAreaFractions(upwind = false).

    Reference geometry kept next to the state (linearInterpolant: exact for the
    piecewise-linear surface): the tet zero-set patches' area vectors mTet_c
    (a consistency check against m_c), their area AS_c and their area-weighted
    centroid xS_c (where the exact normal is evaluated by the tests).

    The explicit plane position (planePosition()), the Gauss theorem for the
    position vector on the liquid sub-volume, exact on planar faces because
    (x - x_c).n_f is constant on a face:

        p_c = ( 3 alpha_c V_c - sum_f alpha_f (x_f - x_c).S_f^out ) / |m_c|,

    the plane {x : (x - x_c).n_c = p_c}, n_c = m_c/|m_c|. Together with the
    normal it recovers the PLIC plane of a cell from (alpha_c, alpha_f) alone,
    without iteration or neighbours; for a state no plane realizes, a volume
    defect dV moves it by 3 dV/|m_c| (three times a volume-matching shift).

      exactSphere         The exact state of an implicitSphere on an axis-aligned
                          hex mesh (davofSphereGeometry.H): face fractions from
                          the circle-polygon intersection, cell fractions from
                          an adaptive quadrature of the exact cross-section.

    Dictionary (fvSolution):
        davof
        {
            alphaSource  detrixheAslam;   // | quadraticFaces | planePhaseIndicator | exactSphere
            wispTol      1e-3;                // interface cell: |m_c| > wispTol*V_c^(2/3)
        }
\*---------------------------------------------------------------------------*/

#ifndef davofState_H
#define davofState_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "levelSetImplicitSurfaces.H"

namespace Foam
{

class davofState
{
    const fvMesh& mesh_;
    const dictionary& dict_;
    const word alphaSource_;
    const scalar wispTol_;

    volScalarField alpha_;          // alpha.davof
    scalarField alphaf_;            // one value per polyMesh face
    // The area normal from the barycentric face triangulation: for every
    // face the fan of triangles (x_f, p_i, p_{i+1}) over the face's own
    // points, each triangle's wetted fraction times its area vector, summed
    // into the owner (-) and neighbour (+) cells ON THE FLY during the
    // surface-based initialisations; nothing is stored per face. Equal to
    // -sum_f alpha_f S_f on planar faces; on a warped face the scalar
    // fraction cannot represent the wetted area vector (measured
    // 2026-10-02: 8 % normal error on a plane at 0.1 h vertex perturbation,
    // the triangle sum 1e-13). areaNormal() uses it as m_c when available.
    vectorField mTri_;
    // The face moment of the explicit position from the same fan, per
    // cell: sum_f sum_t (x_wet,t - x_c).S_wet,t^out with the wet centroid
    // x_wet,t of every triangle (the polygon centroid of the linear cut;
    // the quadratic path keeps its exact area and takes the centroid from
    // the linear cut of the same triangle, which is exact for a plane and
    // immaterial on a planar face, where (x - x_c).n is constant).
    scalarField momentTri_;
    bool hasMTri_;
    surfaceScalarField alphafOut_;  // alphaf.davof (for writing)
    volVectorField m_;              // m.davof, the DAVOF area normal [m^2]
    volVectorField xS_;             // xS.davof, patch centroid [m]
    volScalarField AS_;             // AS.davof, patch area [m^2]
    vectorField mTet_;              // sum of the tet zero-set area vectors
    scalar maxAlphaClip_;

    // The explicit plane position (planePosition()):
    volScalarField pPlane_;         // pPlane.davof, (x - x_c).n on the plane [m]
    volVectorField xPlane_;         // xPlane.davof, centroid of the plane-cell polygon
    volScalarField APlane_;         // APlane.davof, area of that polygon [m^2]
    volScalarField alphaPlane_;     // alphaPlane.davof, cell fraction cut by the plane
    scalar maxVolDiffPlane_;        // max_c |alphaPlane_c - alpha_c|
    scalar maxAreaDiffPlane_;       // max_c |A_polygon n - m_c| / V_c^(2/3)
    label nFaceFallback_;           // quadraticFaces: faces kept at the DA value

    //- Cell fractions and zero-set patches of cellI from vertex values given
    //  by phiCell(c), phiFace(c, f), phiPoint(c, p).
    template<class PhiCell, class PhiFace, class PhiPoint>
    void cellFromValues
    (
        const label cellI,
        PhiCell&& phiCell,
        PhiFace&& phiFace,
        PhiPoint&& phiPoint
    );

    //- Liquid area fraction of faceI from the triangle fan about the face
    //  centre with vertex values phiFace(f), phiPoint(p).
    template<class PhiFace, class PhiPoint>
    scalar faceFraction
    (
        const label faceI,
        PhiFace&& phiFace,
        PhiPoint&& phiPoint
    ) const;

    void resetBulkCell(const label cellI, const scalar alphaValue);

    void fillAlphafOut();

public:

    davofState(const fvMesh& mesh, const dictionary& davofDict);

    //- The canonical source name (the alias linearInterpolant -> detrixheAslam)
    static word canonicalSource(const word& name);

    //- alphaSource detrixheAslam (alias linearInterpolant)
    void computeFromImplicitSurface(const implicitSurface& surface);

    //- alphaSource quadraticFaces: detrixheAslam, then the face fractions from
    //  the exactly cut quadratic interpolant of the signed distance at the
    //  vertices and edge midpoints of every fan triangle. The tet reference
    //  (mTet, AS, xS) stays the chord surface's, so maxConsistency() is the
    //  O(h) chord-to-quadratic difference, not a gate, as for exactSphere.
    void computeFromQuadraticFaces(const implicitSurface& surface);

    //- Faces whose quadratic cut fell back to the Detrixhe-Aslam value
    label nFaceFallback() const { return nFaceFallback_; }

    //- alphaSource exactSphere: the exact face fractions (planar faces cut by
    //  the ball: davofSphereGeometry.H) and the exact cell fractions of
    //  axis-aligned hexahedral cells (adaptive quadrature of the exact
    //  cross-section; any other cell shape is a FatalError). The tet zero-set
    //  reference (mTet, AS, xS) stays the linear interpolant's, so
    //  maxConsistency() then measures the state difference exact - linear,
    //  O(h) by construction, not an implementation defect.
    void computeFromExactSphere(const implicitSphere& sphere);

    //- alphaSource planePhaseIndicator (analytic: the tangent-plane oracle,
    //  nullptr for the least-squares plane of psi). EVERY cell gets a plane,
    //  not only the narrow band: leiaSetFields builds its band in
    //  createFields.H before psi is set (a uniform field passes the
    //  sign-change test on every face), so its Detrixhe-Aslam initialisation
    //  runs on every cell, and a plane fitted 0.5-0.9 h from the interface
    //  still clips a corner of its cell. Reproducing that bit for bit is what
    //  MAX_ALPHA_DIFF_DA = 0 in leiaTestDavofNormal means; the planes of cells
    //  farther away never cut them.
    void computeFromPlanePhaseIndicator
    (
        const volScalarField& psi,
        const implicitSurface* analytic
    );

    //- The wetted part of face faceI from the fan triangles (x_f, p_i,
    //  p_{i+1}) over the face points under the linear interpolant: its
    //  area vector S = sum_t a_t S_t and its first moment M = sum_t
    //  x_wet,t . (a_t S_t), so that a cell with centre x_c sees the moment
    //  M - x_c.S. Computed where it is summed, never stored.
    template<class PhiFace, class PhiPoint>
    void faceWettedTriangles
    (
        const label faceI,
        PhiFace&& phiFace,
        PhiPoint&& phiPoint,
        vector& S,
        scalar& M
    ) const;

    //- m_c = -sum_f alpha_f S_f^out, or the face-triangulation sum mTri_
    //  when the state was initialised from a surface
    void areaNormal();

    //- max_c |m_c - mTet_c| / V_c^(2/3), reduced over ranks
    scalar maxConsistency() const;

    //- The explicit position of the DAVOF plane in every interface cell
    //  (after areaNormal()), the Gauss theorem for x on the liquid sub-volume:
    //
    //      p_c = ( 3 alpha_c V_c - sum_f alpha_f (x_f - x_c).S_f^out ) / |m_c|
    //
    //  with n_c = m_c/|m_c| pointing out of the liquid, the plane
    //  {x : (x - x_c).n_c = p_c}. Exact on planar faces ((x - x_c).n_f is
    //  constant there), no iteration, nothing beyond (alpha_c, alpha_f) and
    //  the mesh geometry. The plane is then cut through the cell on the same
    //  Detrixhe-Aslam tets (exact for a plane): xPlane_c, APlane_c and the
    //  cut volume fraction alphaPlane_c; the last two against alpha_c and
    //  m_c are the realizability diagnostics (round-off for a plane, O(h^2)
    //  for a curved surface whose state no single plane realizes).
    void planePosition();

    scalar maxVolDiffPlane() const { return maxVolDiffPlane_; }
    scalar maxAreaDiffPlane() const { return maxAreaDiffPlane_; }
    const volScalarField& pPlane() const { return pPlane_; }
    const volVectorField& xPlane() const { return xPlane_; }
    const volScalarField& APlane() const { return APlane_; }
    const volScalarField& alphaPlane() const { return alphaPlane_; }

    //- Interface cell: |m_c| > wispTol V_c^(2/3); wisp: 0 < |m_c| <= that.
    scalar wispThreshold(const label cellI) const;
    bool isInterfaceCell(const label cellI) const;
    bool isWisp(const label cellI) const;

    // Access
    const word& alphaSource() const { return alphaSource_; }
    scalar wispTol() const { return wispTol_; }
    volScalarField& alpha() { return alpha_; }
    const volScalarField& alpha() const { return alpha_; }
    const scalarField& alphaf() const { return alphaf_; }
    const volVectorField& m() const { return m_; }
    const volVectorField& xS() const { return xS_; }
    const volScalarField& AS() const { return AS_; }
    const vectorField& mTet() const { return mTet_; }
    const vectorField& mTri() const { return mTri_; }
    const scalarField& momentTri() const { return momentTri_; }
    bool hasMTri() const { return hasMTri_; }
    scalar maxAlphaClip() const { return maxAlphaClip_; }

    //- Write alpha.davof, alphaf.davof, m.davof, xS.davof, AS.davof
    void write();
};

} // End namespace Foam

#endif

// ************************************************************************* //
