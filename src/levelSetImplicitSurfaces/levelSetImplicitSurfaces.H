/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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

Description
    Implicit surfaces with exact normal and curvature. 

\*---------------------------------------------------------------------------*/

#ifndef levelSetImplicitSurfaces_H
#define levelSetImplicitSurfaces_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "dictionary.H"

namespace Foam {

    class implicitSurface
    {
        public: 

            TypeName("implicitSurface");

            declareRunTimeSelectionTable
            (
                autoPtr,
                implicitSurface, 
                dictionary, 
                (
                    const dictionary& dict 
                ), 
                (dict)
            );

            static autoPtr<implicitSurface> New(
                const word& name, 
                const dictionary& dict 
            );

            implicitSurface() = default;

            implicitSurface(const dictionary& dict) {};

            virtual ~implicitSurface() = default;

            virtual scalar value(const vector&) const = 0;
            virtual vector grad(const vector&) const = 0;
            virtual scalar curvature(const vector&) const = 0;
    };

    class implicitPlane : public implicitSurface
    {
        vector position_; 
        vector normal_; 

        scalar gradient_;

        public:

            TypeName ("implicitPlane");

            implicitPlane() = default;

            implicitPlane(vector position, vector normal);

            implicitPlane(const dictionary& dict);

            virtual ~implicitPlane() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector position() const;

            vector normal() const;
    };

    class hesseNormalPlane : public implicitSurface
    {
        vector n_; 
        scalar d_; 
        
        public:

            TypeName ("hesseNormalPlane");
        
            hesseNormalPlane(vector n, scalar d);

            hesseNormalPlane(const dictionary&);

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            const vector& normal() const;
            
            scalar cutoff() const;
    };

    //- EXACT Euclidean signed distance to an infinite slab of half-width W,
    //  centred at c, with unit axis d; negative inside the slab:
    //
    //      psi(x) = |(x - c) . d| - W
    //
    //  Dictionary:
    //
    //      type       implicitSlab;
    //      centre     (0 0 0);
    //      direction  (1 0 0);
    //      halfWidth  0.25;
    //
    //  The two interfaces are the parallel planes (x - c).d = +W and
    //  (x - c).d = -W, so the slab thickness is 2W.
    //
    //  THIS IS EXACT, NOT AN APPROXIMATION. For a point between two parallel
    //  planes the closest boundary point is its orthogonal projection onto
    //  the nearer plane, and for a point outside onto the plane it lies
    //  beyond; in both cases the distance is | (x - c).d | - W in magnitude,
    //  so |grad psi| = 1 everywhere except on the single plane named below.
    //  A slab case therefore STARTS from a true signed-distance field, and any
    //  |grad psi| - 1 measured afterwards was produced by the run, not
    //  inherited from the initialisation. Contrast implicitEllipsoid, whose
    //  value() is the algebraic level set and is NOT a distance;
    //  signedDistanceEllipse exists for exactly that reason.
    //
    //  ITS MEDIAL AXIS IS THE PLANE (x - c).d = 0, midway between the two
    //  interfaces, and that plane is the entire reason this surface exists.
    //  There psi attains its minimum -W and is NOT differentiable: grad psi
    //  jumps from -d to +d, and the closest-point map jumps from the
    //  (x - c).d = +W interface to the (x - c).d = -W one. A velocity
    //  extension of closestPoint type, Uext(x) = v(foot point of x), therefore
    //  carries a JUMP of magnitude |v(x_R) - v(x_L)| across the medial axis
    //  wherever the medial axis lies inside the extension band, which happens
    //  exactly when
    //
    //      nLayers * h  >  W,
    //
    //  i.e. the controlling parameter of the diagnosis is W / (nLayers * h).
    //  This is what a single-interface case cannot expose: with one interface
    //  every band cell has the same foot point, closestPoint is exact, and
    //  |grad psi| = 1 is preserved with no source at all. A circle of radius
    //  0.15 with a 3-cell band never has its medial axis (its centre) inside
    //  the band either, which is why the existing static closestPoint gate
    //  reports O(h^1.98) and never saw this; a filament thinning towards one
    //  cell does.
    //
    //  ON THE MEDIAL AXIS itself grad() returns the one-sided limit +d, and
    //  curvature() returns 0 there as everywhere else (two parallel planes
    //  have zero mean curvature). Both are conventions for a measure-zero set
    //  on which the exact distance function is genuinely non-differentiable.
    class implicitSlab : public implicitSurface
    {
        vector centre_;
        vector direction_;
        scalar halfWidth_;

        //- Normalise and store the slab axis; fails loudly on a zero vector,
        //  which would otherwise make every value() equal to -halfWidth_ and
        //  produce an interface-free field that looks like a solver failure.
        void setDirection(const vector& d);

        public:

            TypeName ("implicitSlab");

            implicitSlab() = default;

            implicitSlab(vector centre, vector direction, scalar halfWidth);

            implicitSlab(const dictionary& dict);

            virtual ~implicitSlab() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector centre() const;

            vector direction() const;

            scalar halfWidth() const;
    };

    class implicitSphere : public implicitSurface
    {
        vector center_; 
        scalar radius_; 

        public:

            TypeName ("implicitSphere");

            implicitSphere() = default;

            implicitSphere(vector center, scalar radius);

            implicitSphere(const dictionary& dict);

            virtual ~implicitSphere() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector center() const;

            scalar radius() const;
    };

    class implicitSlottedSphere : public implicitSurface
    {
        /*
        Enright, D., Fedkiw, R. P., Ferziger, J., & Mitchell, I. (2002). A Hybrid
        Particle Level Set Method for Improved Interface Capturing. Journal of
        Computational Physics (Vol. 183). https://doi.org/10.1006/jcph.2002.7166
        */

        vector center_; 
        scalar radius_; 
        vector minSlotCorner_; 
        vector maxSlotCorner_;

        public:

            TypeName ("slottedSphere");

            implicitSlottedSphere() = default;

            implicitSlottedSphere(
                vector center, 
                scalar radius,
                vector minSlotCorner, 
                vector maxSlotCorner
            );

            implicitSlottedSphere(const dictionary& dict);

            virtual ~implicitSlottedSphere() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector center() const;

            scalar radius() const;

            vector minSlotCorner() const;

            vector maxSlotCorner() const;
    };   

    class implicitEllipsoid : public implicitSurface
    {
        vector center_; 
        vector axes_; 
        vector axesSqr_;

        void setAxesSqr(const vector& axes);

        public:

            TypeName ("implicitEllipsoid");

            implicitEllipsoid(vector center, vector axes);

            implicitEllipsoid(const dictionary& dict);

            virtual ~implicitEllipsoid() = default;

            virtual scalar value(const vector& x) const;
            
            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector center() const;

            vector axes() const;
    };

    //- Exact Euclidean signed distance to an axis-aligned TRIAXIAL ellipsoid
    //  (three distinct half-axes). The missing 3D twin of
    //  signedDistanceEllipse, and the reference the curvature campaign named:
    //  the only way to test orientation and varying curvature together in 3D.
    //
    //  closestPoint solves, in the first octant with p_i = |x_i - c_i|,
    //      F(t) = sum_i (a_i p_i / (t + a_i^2))^2 - 1 = 0
    //  (Eberly's formulation). F is strictly decreasing on (-a_min^2, inf) with
    //  F -> +inf at the left end and F -> -1 at t -> inf, so the root is UNIQUE
    //  and safeguarded bisection + Newton cannot converge to a wrong branch --
    //  unlike a naive Lagrange iteration, which has spurious stationary points
    //  for interior x.
    //
    //  value()     = +/- |x - q|, sign from the algebraic value (outside > 0);
    //  grad()      = unit vector sign * (x - q)/|x - q| (|grad| = 1 exactly);
    //  curvature() = the SURFACE's total curvature (kappa1 + kappa2, the
    //      div(n) convention that gives 2/R on a sphere) AT THE CLOSEST POINT
    //      q -- the same convention signedDistanceEllipse uses in 2D. It is
    //      evaluated by the machine-generated implicitEllipsoid::curvature
    //      formula at q, so no independently derived differential geometry can
    //      disagree with the algebraic surface.
    class signedDistanceEllipsoid : public implicitSurface
    {
        vector center_;
        vector axes_;
        implicitEllipsoid ellipsoid_;   // supplies value sign + curvature at q

        //- Closest point on the ellipsoid surface to x (world coordinates).
        vector closestPoint(const vector& x) const;

        public:

            TypeName ("signedDistanceEllipsoid");

            signedDistanceEllipsoid(vector center, vector axes);

            signedDistanceEllipsoid(const dictionary& dict);

            virtual ~signedDistanceEllipsoid() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector center() const { return center_; }

            vector axes() const { return axes_; }
    };

    //- Exact Euclidean signed distance to an axis-aligned 2-D ellipse.
    //  The z coordinate and axes.z are ignored, making this suitable for the
    //  one-cell-thick empty-patch droplet benchmarks.  This is deliberately a
    //  separate type: implicitEllipsoid retains its historical algebraic
    //  level-set value for backward compatibility.
    class signedDistanceEllipse : public implicitSurface
    {
        vector center_;
        vector axes_;

        scalar closestParameter(const vector& x) const;

        public:

            TypeName ("signedDistanceEllipse");

            signedDistanceEllipse(vector center, vector axes);

            signedDistanceEllipse(const dictionary& dict);

            virtual ~signedDistanceEllipse() = default;

            virtual scalar value(const vector& x) const;

            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector center() const;

            vector axes() const;
    };

    class implicitSinc : public implicitSurface
    {
        vector origin_; 
        scalar amplitude_; 
        scalar omega_;

        public:

            TypeName ("implicitSinc");

            implicitSinc(vector origin, scalar amplitude, scalar omega);

            implicitSinc(const dictionary& dict);

            virtual ~implicitSinc() = default;

            virtual scalar value(const vector& x) const;
            
            virtual vector grad(const vector& x) const;

            virtual scalar curvature(const vector& x) const;

            vector origin() const;

            scalar amplitude() const;

            scalar omega() const;
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
