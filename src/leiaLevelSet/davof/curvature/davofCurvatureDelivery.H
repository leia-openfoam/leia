/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCurvatureDelivery

    From a curvature model of the DAVOF state (kappa, K at the interface
    centroids) to the two fields a balanced-force CSF consumes, with the
    parallel-surface closed form the level-set solver uses, because the
    interface centroid is neither the cell centre nor a face centre.

    Total curvature kappa = kappa_1 + kappa_2 and Gaussian curvature K of the
    parallel surface at signed offset d (positive on the gas side) from an
    interface with kappa, K:

        kappa_i(d) = kappa_i/(1 + d kappa_i)
        kappa(d)   = (kappa + 2 d K)/(1 + d kappa + d^2 K)      (forward)
        K(d)       = K/(1 + d kappa + d^2 K)
        kappa      = (kappa_d - 2 K_d d)/(1 - d kappa_d + K_d d^2) (inverse,
                     the solver's parallelSurfaceInverse, trust-guarded)

    Sphere: 2/R -> 2/(R + d) -> 2/R.

    Cell field (contour-referenced, on the FORCE BAND = cells adjacent to a
    face with alpha_owner != alpha_neighbour): an interface cell takes its own
    model, d_c = offset(c, x_c), kappaCell = forward(kappa, d_c, K), KCell =
    K(d_c); a band cell without a model takes the Kang-weighted (1/|d|) forward
    values of its face-adjacent interface cells, across coupled patches through
    the packed models; no source: unfilled (excluded from the norms).

    Face field on the ACTIVE faces (alpha jumps across the face and a model
    exists on at least one side; internal faces and the owner side of coupled
    patches, so a decomposed run scores the same faces as the serial one; the
    faces of wisps, a jump with no model anywhere, are skipped and counted),
    two deliveries:
      inverse  kappa_f = parallelSurfaceInverse(interpolate(kappaCell)_f, d_f, K_f)
               with d_f and K_f (the parallel surface's Gaussian curvature at
               the face centre) Kang-combined from the owner's and neighbour's
               offset(x_f) (weights: the OTHER cell's |d|, as
               stabilizedFootPointFaceCurvature.H does with |psi|), one-sided
               when one side has no model -- the level-set solver's path;
      foot     kappa_f = the models' interface curvature at the ray foot of x_f,
               Kang-combined the same way -- the comparison row, no parallel-
               surface algebra.
\*---------------------------------------------------------------------------*/

#ifndef davofCurvatureDelivery_H
#define davofCurvatureDelivery_H

#include "davofCurvature.H"

namespace Foam
{

//- Total curvature of the parallel surface at offset d (guarded: past half the
//  local curvature radius the value is returned unchanged).
inline scalar parallelSurfaceForward
(
    const scalar kappa,
    const scalar d,
    const scalar K
)
{
    const scalar denom = 1.0 + d*kappa + K*d*d;
    return (mag(denom) > 0.5) ? (kappa + 2.0*K*d)/denom : kappa;
}

//- Gaussian curvature of the parallel surface at offset d.
inline scalar parallelSurfaceGaussianForward
(
    const scalar kappa,
    const scalar d,
    const scalar K
)
{
    const scalar denom = 1.0 + d*kappa + K*d*d;
    return (mag(denom) > 0.5) ? K/denom : K;
}

//- The parallel-surface inverse with trust guard; falls back to the
//  uncorrected value past half the local curvature radius. Verbatim from
//  applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/
//  stabilizedFootPointFaceCurvature.H (a solver header a library cannot
//  include); keep the two identical.
inline scalar parallelSurfaceInverse
(
    const scalar kf,
    const scalar d,
    const scalar K
)
{
    const scalar denom = 1.0 - d*kf + K*d*d;
    return (mag(denom) > 0.5) ? (kf - 2.0*K*d)/denom : kf;
}


class davofCurvatureDelivery
{
    const fvMesh& mesh_;
    const davofState& state_;
    const davofCurvature& model_;

    // Cells
    scalarField kappaCell_;
    scalarField KCell_;
    boolList bandCell_;
    boolList filledCell_;

    // Faces (indexed by the mesh face; only active faces are meaningful)
    scalarField kappaFaceInv_;
    scalarField kappaFaceFoot_;
    scalarField dFace_;
    scalarField KFace_;
    boolList activeFace_;
    boolList oneSidedFace_;

    label nBandCells_;
    label nUnfilledCells_;
    label nActiveFaces_;
    label nOneSidedFaces_;
    label nSkippedFaces_;     // alpha jumps, no model on either side (wisps)

public:

    davofCurvatureDelivery
    (
        const fvMesh& mesh,
        const davofState& state,
        const davofCurvature& model
    );

    const scalarField& kappaCell() const { return kappaCell_; }
    const scalarField& KCell() const { return KCell_; }
    const boolList& bandCell() const { return bandCell_; }
    const boolList& filledCell() const { return filledCell_; }

    const scalarField& kappaFaceInv() const { return kappaFaceInv_; }
    const scalarField& kappaFaceFoot() const { return kappaFaceFoot_; }
    const scalarField& dFace() const { return dFace_; }
    const scalarField& KFace() const { return KFace_; }
    const boolList& activeFace() const { return activeFace_; }
    const boolList& oneSidedFace() const { return oneSidedFace_; }

    label nBandCells() const { return nBandCells_; }
    label nUnfilledCells() const { return nUnfilledCells_; }
    label nActiveFaces() const { return nActiveFaces_; }
    label nOneSidedFaces() const { return nOneSidedFaces_; }
    label nSkippedFaces() const { return nSkippedFaces_; }
};

} // namespace Foam

#endif

// ************************************************************************* //
