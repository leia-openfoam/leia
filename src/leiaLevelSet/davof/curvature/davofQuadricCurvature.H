/*---------------------------------------------------------------------------*\
    leia :: davof :: davofQuadricCurvature   (type quadricFit)

    A quadric graph per interface cell fitted to the Hermite data of the
    point-neighbour ring: in the frame (o = xPlane_c, e1, e2, e3 = n_c), scaled
    by h_c = V_c^(1/3),

        F(U,V) = C0 + C1 U + C2 V + C3 U^2/2 + C4 U V + C5 V^2/2,

    every neighbouring interface cell j contributing two SLOPE rows from its
    normal, F_U(U_j,V_j) = s_j1, F_V(U_j,V_j) = s_j2 with
    s_j = -(n_j.e1, n_j.e2)/(n_j.e3), and one POSITION row F(U_j,V_j) = W_j
    from its centroid, weighted by positionWeight.

    positionWeight 0 (the default) lets the NORMALS decide the second
    derivatives (C1..C5 from the slope rows alone, C0 afterwards from the
    positions): the centroids carry a non-smooth O(h^2) error, and a fit in
    which they decide the curvature has an O(1) floor, the point-fit floor
    recorded in isoGeometricCurvature.H; the normals' O(h^2) error over a ring
    of about 2h gives O(h). positionWeight 1 is that naive Hermite fit, kept as
    the comparison arm.

    Curvature of the graph at (U,V) with f_u = C1 + C3 U + C4 V, f_v = C2 +
    C4 U + C5 V, f_uu = C3/h_c, f_uv = C4/h_c, f_vv = C5/h_c:

        kappa = -[(1+f_v^2) f_uu - 2 f_u f_v f_uv + (1+f_u^2) f_vv]/(1+f_u^2+f_v^2)^{3/2}
        K     = (f_uu f_vv - f_uv^2)/(1+f_u^2+f_v^2)^2

    (the minus sign is div(n) of the graph with n along +e3, out of the
    liquid: a droplet curves away from e3 and gets kappa = +2/R).

    Offset of a point x (cell or face centre): local (U,V,W), the level-set-
    like psi = W - F(U,V), its gradient g = (-F_U, -F_V, 1) and h_nn = n.H.n
    with H = -Hess(F); the STABLE quadratic root of the normal-projected SL
    design note (docs/normal-projected-semi-lagrangian, Sec. 5.2-5.3),
    d = 2 psi/(|g| + sqrt(D)), D = |g|^2 - 2 psi h_nn, accepted while
    D >= offsetBeta |g|^2, first order psi/|g| otherwise; times h_c. kappa and
    K are evaluated at the ray foot.

    Dictionary:
        <name>
        {
            type            quadricFit;
            positionWeight  0;        // 0: normals decide the curvature; 1: naive Hermite fit
            weighting       area;     // area: sqrt(|m_j|/|m_c|) | uniform
            minNeighbours   5;        // fewer usable neighbours: fallback (plane, kappa = K = 0)
            alignTol        0.2;      // neighbours with n_j.n_c <= alignTol are dropped
            ridge           0;        // Tikhonov ridge on the normal equations (relative)
            offsetBeta      0.25;     // acceptance of the quadratic offset root
        }

    Stencil: OpenFOAM's CPCCellToCellStencil (point neighbours, global cell
    indices) with the lists of the non-interface cells cleared before the
    extendedCentredCellToCellStencil builds its mapDistribute, so the ring data
    (centroid, area vector, weight) are gathered for the interface cells only:
    for every cell they would be 27 x 7 scalars x 4.1e6 cells = 6 GB at
    N = 160, which is what killed the first ladder run. Parallel-aware
    through the map.
\*---------------------------------------------------------------------------*/

#ifndef davofQuadricCurvature_H
#define davofQuadricCurvature_H

#include "davofCurvature.H"

namespace Foam
{

class davofQuadricCurvature
:
    public davofCurvature
{
    const scalar positionWeight_;
    const word weighting_;
    const label minNeighbours_;
    const scalar alignTol_;
    const scalar ridge_;
    const scalar offsetBeta_;

    // The local models, stored for the interface cells only (modelIndex_[c]
    // < 0 elsewhere): at N = 160 a per-cell store of 19 scalars is 0.6 GB.
    labelList modelIndex_;
    vectorField origin_;
    vectorField e1_;
    vectorField e2_;
    vectorField e3_;
    scalarField hc_;
    List<FixedList<scalar, 6>> coeffs_;

    //- kappa and K of the graph with coefficients C (scale hc) at (U, V)
    static void curvatureOf
    (
        const FixedList<scalar, 6>& C,
        const scalar hc,
        const scalar U,
        const scalar V,
        scalar& kappa,
        scalar& K
    );

    //- Orthonormal frame with e3 given
    static void frameOf(const vector& e3, vector& e1, vector& e2);

public:

    TypeName("quadricFit");

    davofQuadricCurvature
    (
        const word& name,
        const fvMesh& mesh,
        const davofState& state,
        const dictionary& dict
    );

    virtual ~davofQuadricCurvature() = default;

    virtual void compute();

    //- Pack layout: [0] filled, [1..3] origin, [4..6] e1, [7..9] e2,
    //  [10..12] e3, [13] hc, [14..19] C0..C5.
    virtual bool pack(const label c, davofModelPack& p) const;

    virtual bool offsetPacked
    (
        const davofModelPack& p,
        const point& x,
        scalar& d,
        scalar& kappa,
        scalar& K
    ) const;

    scalar positionWeight() const { return positionWeight_; }
};

} // namespace Foam

#endif

// ************************************************************************* //
