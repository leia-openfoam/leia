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
    Foam::davofNormalModels::ringArea

Description
    The area-weighted average of the ring's DAVOF normals,

        n_k = d_k / |d_k|,   d_k = sum_j m_j,

    the sum over the point-neighbour ring R(k) of the cell (the cells sharing
    at least one point with it, the cell included; zoneDistribute's stencil,
    across processor boundaries), where only the interface cells contribute
    (m_j = 0 otherwise). This is the DAVOF normal of the union of the ring's
    cells: exact for a plane on any mesh, no division by a small area (the
    division is by the ring's total patch area |d_k|, the mechanism of its
    stability against the own normal m_k/|m_k|; see the tilt diagnostic of
    2026-10-26). Default applyTo all. The computation is the base class's
    davofNormalModel::ringAreaNormals(), shared with the noNormalFallback
    ringArea rule of the other models.

    The stencil variants (2026-11-08, Tomislav's task: the stability and the
    accuracy of the ring normal's stencil). Every key is opt-in; at the
    defaults the model runs the production path above, bit for bit:

        ringStencil     one;       // one | two: the point-neighbour ring, or
                                   // the union of the rings of the ring's
                                   // cells (every cell within two point
                                   // adjacencies, each counted once; the
                                   // rings of ring cells on other processors
                                   // are requested from their owners)
        ringKernel      none;      // none | gaussian | wendland: the weight
                                   // w_j = |Sigma_j| K(|x_j - x_k| / h_k) with
                                   // the cell centres x and h_k = V_k^(1/3);
                                   // d_k = sum_j K_j m_j. gaussian: K(q) =
                                   // exp(-q^2) cut off at ringKernelSupport
                                   // (R. A. Gingold, J. J. Monaghan, Mon. Not.
                                   // R. Astron. Soc. 181 (1977) 375-389,
                                   // doi:10.1093/mnras/181.3.375); wendland:
                                   // the C2 function of minimal degree,
                                   // K(q) = (1 - q/s)^4 (1 + 4 q/s) for q < s
                                   // (H. Wendland, Adv. Comput. Math. 4 (1995)
                                   // 389-396, doi:10.1007/BF02123482)
        ringKernelSupport 2;       // s, in units of h_k
        ringAngleFilter 0;         // theta_max in degrees (0: off): a ring
                                   // cell whose own normal m_j deviates from
                                   // the ring sum d_k by more than theta_max
                                   // is dropped and the sum recomputed; a
                                   // cell whose filtered sum vanishes keeps
                                   // d_k (counted)
        ringPlaneFit    none;      // none | mean | own: the normal of the
                                   // weighted least-squares plane through the
                                   // ring's patch centroids x_j (the
                                   // positioned planes' centroids, the first
                                   // pass) with the area weights w_j =
                                   // |Sigma_j|: the eigenvector of the
                                   // smallest eigenvalue of the weighted
                                   // covariance sum_j w_j (x_j - c)(x_j - c)^T
                                   // (the principal-component plane, K.
                                   // Pearson, Philos. Mag. 2 (1901) 559-572,
                                   // doi:10.1080/14786440109462720; H. Hoppe,
                                   // T. DeRose, T. Duchamp, J. McDonald, W.
                                   // Stuetzle, Computer Graphics 26 (1992)
                                   // 71-78, doi:10.1145/133994.134011), the
                                   // sign from d_k. mean: c the weighted
                                   // centroid. own: c = x_k, the cell's own
                                   // patch centroid, with the weights
                                   // corrected to w'_j = w_j (1 + lambda .
                                   // (x_j - x_k)) so that the centring
                                   // condition sum_j w'_j (x_j - x_k) = 0
                                   // holds (the linear moment correction of
                                   // the reproducing kernel particle method,
                                   // W. K. Liu, S. Jun, Y. F. Zhang, Int. J.
                                   // Numer. Methods Fluids 20 (1995)
                                   // 1081-1106, doi:10.1002/fld.1650200824).
                                   // Fewer than three usable centroids, a
                                   // collinear set (second eigenvalue below
                                   // ringPlaneFitTol times the trace) or an
                                   // indefinite covariance: d_k, counted as
                                   // a fit fallback
        ringPlaneFitTol 1e-10;

    The variants compose (a kernel and an angle filter on one or two rings;
    the plane fit reads the ring of ringStencil); the studies of 2026-11-08
    test them one at a time.

SourceFiles
    ringAreaNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_ringArea_H
#define davofNormalModels_ringArea_H

#include "davofNormalModel.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                           Class ringArea Declaration
\*---------------------------------------------------------------------------*/

class ringArea
:
    public davofNormalModel
{
protected:

    // Protected Data

        //- m_c in the interface cells, zero elsewhere (the field exchanged
        //  across processors; its boundary values are zero, so boundary
        //  faces in the stencil contribute nothing)
        volVectorField mInterface_;

        //- The cells mInterface_ was last filled in (zeroed there before
        //  the next fill; zero elsewhere)
        labelList filled_;

        //- The stencil variants (2026-11-08); at the defaults the
        //  production path runs unchanged
        word stencil_;
        word kernel_;
        scalar kernelSupport_;
        scalar angleFilter_;
        scalar cosAngleFilter_;
        word planeFit_;
        scalar planeFitTol_;

        //- ringPlaneFitDebug N: the ring data, the covariance and its
        //  eigen-decomposition of the first N fitted cells of the first
        //  call printed (0: off)
        label planeFitDebug_;

        //- Counts of the last variant call (reduced): cells dropped by the
        //  angle filter, cells whose filtered sum vanished (kept d_k),
        //  remote ring cells whose rings were requested (two rings)
        label nAngleDropped_;
        label nAngleKept_;
        label nRemoteRings_;


    // Protected Member Functions

        //- Fill mInterface_ from the scheme's area normal
        void fillInterfaceNormals();

        //- Whether every stencil key is at its default
        bool production() const;

        //- The kernel value at q = r/h (1 for none)
        scalar kernelValue(const scalar q) const;

        //- The eigen-decomposition of a symmetric 3 x 3 matrix by cyclic
        //  Jacobi rotations (C. G. J. Jacobi, J. Reine Angew. Math. 30
        //  (1846) 51-94, doi:10.1515/crll.1846.30.51): eVals ascending,
        //  the rows of eVecs the unit eigenvectors in that order. Used in
        //  place of OpenFOAM's eigenValues()/eigenVectors() (the roots of
        //  the characteristic cubic), which on the covariance of a
        //  coplanar ring (eigenvalues 1e-24, 1.4e-8, 3.6e-8 in units of
        //  h^2 |Sigma|) returned a normal 45 degrees off without a
        //  warning (measured 2026-11-08 on the translated plane, L1(n)
        //  0.77); the Jacobi rotations give the plane's normal to 2e-15
        static void jacobiEigen
        (
            const symmTensor& T,
            vector& eVals,
            tensor& eVecs
        );

        //- The ring of every interface cell as a list of (global index,
        //  m) pairs: ringStencil one (the stencil) or two (the union of the
        //  stencils of the stencil's cells, across processors); the
        //  exchange of mInterface_ is set up here (zone: the interface
        //  cells, or the interface cells and their point neighbours).
        //  Returns the stencil's exchange map of mInterface_ for further
        //  reads through exchange.getValue (one ring: every ring index is
        //  in the stencil of the cell; two rings: the map covers the
        //  stencils of the zone's local cells, the remote rings' values are
        //  in ringM)
        void ringMembers
        (
            labelListList& ringIdx,
            List<vectorField>& ringM
        );

        //- The variant normals (any key off its default)
        void variantNormals(vectorField& n);


public:

    //- Runtime type information
    TypeName("ringArea");


    //- Construct from the davof dictionary and the scheme
    ringArea
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~ringArea() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);

        //- The plane fit reads the ring's positioned planes (their
        //  centroids) from the first pass
        virtual bool needsPositionedPlanes() const
        {
            return planeFit_ != "none";
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
