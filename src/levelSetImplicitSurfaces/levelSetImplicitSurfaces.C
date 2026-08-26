/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 Tomislav Maric, TU Darmstadt
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

\*---------------------------------------------------------------------------*/

#include "doubleFloat.H"
#include "levelSetImplicitSurfaces.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>
#include "scalarField.H"
#include "mathematicalConstants.H"

namespace Foam {

// * * * * * * * * * * * * Class implicitSurface  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSurface, false);
defineRunTimeSelectionTable(implicitSurface, dictionary);

autoPtr<implicitSurface> implicitSurface::New 
(
    const word& modelType, 
    const dictionary& dict 
)
{
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "implicitSurface",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<implicitSurface>(ctorPtr(dict));
}

// * * * * * * * * * * * * Class implicitPlane  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitPlane, false);
addToRunTimeSelectionTable(implicitSurface, implicitPlane, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

implicitPlane::implicitPlane(vector position, vector normal)
: 
    position_(position), 
    normal_(normal)
{
    normal_ /= Foam::mag(normal_);
}

implicitPlane::implicitPlane(const dictionary& dict)
    :
        position_(dict.get<vector>("position")),
        normal_(dict.get<vector>("normal")),
        gradient_(dict.getOrDefault<scalar>("gradient", 1))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar implicitPlane::value(const vector& x) const
{
    return gradient_*Foam::dot(x - position_, normal_);
}

vector implicitPlane::grad(const vector& x) const
{
    return normal_*gradient_; 
}

vector implicitPlane::position() const
{
    return position_; 
}

vector implicitPlane::normal() const
{
    return normal_; 
}

scalar implicitPlane::curvature(const vector& x) const
{
    return 0;
}

// * * * * * * * * * * * * Class hesseNormalPlane  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hesseNormalPlane, false);
addToRunTimeSelectionTable(implicitSurface, hesseNormalPlane, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hesseNormalPlane::hesseNormalPlane(vector n, scalar d)
: 
    n_(n), 
    d_(d)
{}

hesseNormalPlane::hesseNormalPlane(const dictionary& dict)
    :
        n_(dict.get<vector>("n")),
        d_(dict.get<scalar>("d"))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar hesseNormalPlane::value(const vector& x) const
{
    return Foam::dot(x , n_) + d_;
}

vector hesseNormalPlane::grad(const vector& x) const
{
    return n_; 
}

const vector& hesseNormalPlane::normal() const
{
    return n_; 
}

scalar hesseNormalPlane::curvature(const vector& x) const
{
    return 0;
}

// * * * * * * * * * * * * Class implicitSlab  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSlab, false);
addToRunTimeSelectionTable(implicitSurface, implicitSlab, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void implicitSlab::setDirection(const vector& d)
{
    const scalar magD = Foam::mag(d);

    // A zero axis is not recoverable and must not be tolerated: it would give
    // (x - c).d = 0 for every x, hence psi = -halfWidth everywhere, a field
    // with no zero contour at all. Every downstream diagnosis (interface
    // position, |grad psi|, phase volume) would then report a degenerate
    // result that looks like a solver failure rather than a dictionary typo.
    if (magD < SMALL)
    {
        FatalErrorInFunction
            << "implicitSlab `direction` is the zero vector " << d << "." << nl
            << "The slab axis defines psi(x) = |(x - centre).direction|"
            << " - halfWidth; a zero axis leaves psi = -halfWidth everywhere"
            << " and the field has no interface." << nl
            << "Set `direction` to the (non-zero) slab normal, e.g. (1 0 0)."
            << exit(FatalError);
    }

    direction_ = d/magD;
}

implicitSlab::implicitSlab(vector centre, vector direction, scalar halfWidth)
:
    centre_(centre),
    halfWidth_(halfWidth)
{
    setDirection(direction);
}

implicitSlab::implicitSlab(const dictionary& dict)
:
    centre_(dict.get<vector>("centre")),
    halfWidth_(dict.get<scalar>("halfWidth"))
{
    setDirection(dict.get<vector>("direction"));
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar implicitSlab::value(const vector& x) const
{
    // Exact Euclidean signed distance, negative inside: the closest boundary
    // point of a slab is the orthogonal projection onto the nearer of the two
    // parallel planes, so no minimisation is needed (contrast
    // signedDistanceEllipse::closestParameter, which has to search).
    return Foam::mag(Foam::dot(x - centre_, direction_)) - halfWidth_;
}

vector implicitSlab::grad(const vector& x) const
{
    // grad(|s| - W) = sign(s) d with s = (x - c).d. |grad| = 1 exactly, which
    // is the property a slab case is initialised for. On the medial axis
    // s = 0 the exact distance function has a kink and is not
    // differentiable; the >= 0 branch returns the one-sided limit +d there.
    const scalar s = Foam::dot(x - centre_, direction_);

    return (s >= 0) ? direction_ : -direction_;
}

scalar implicitSlab::curvature(const vector& x) const
{
    // Both interfaces are planes: every principal curvature is zero, so the
    // mean curvature is zero everywhere, including (by the same one-sided
    // convention as grad()) on the medial axis.
    return 0;
}

vector implicitSlab::centre() const
{
    return centre_;
}

vector implicitSlab::direction() const
{
    return direction_;
}

scalar implicitSlab::halfWidth() const
{
    return halfWidth_;
}

// * * * * * * * * * * * * Class implicitSphere * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSphere, false);
addToRunTimeSelectionTable(implicitSurface, implicitSphere, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

implicitSphere::implicitSphere(vector center, scalar radius)
    : 
        center_(center), 
        radius_(radius)
{}

implicitSphere::implicitSphere(dictionary const& dict)
    : 
        center_(dict.get<vector>("center")),
        radius_(dict.get<scalar>("radius"))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar implicitSphere::value(const vector& x) const
{
    return Foam::mag(x - center_) - radius_; 
}

vector implicitSphere::grad(const vector& x) const
{
    scalar x0c0 = x[0] - center_[0];
    scalar x1c1 = x[1] - center_[1];
    scalar x2c2 = x[2] - center_[2];

    return vector(x0c0, x1c1, x2c2) / 
        sqrt(x0c0*x0c0 + x1c1*x1c1 + x2c2*x2c2);
}

vector implicitSphere::center() const
{
    return center_; 
}

scalar implicitSphere::radius() const
{
    return radius_; 
}

// NOTE the 2D (circle) convention: div(n) of a CIRCLE is 1/R; a 3D sphere's
// div(n) is 2/R. The only consumer today is the 2D-only connectedInterface
// analytic oracle, which expects exactly this value -- any future 3D consumer
// must use (nGeometricD - 1)/R instead of this member.
scalar implicitSphere::curvature(const vector& x) const
{
    return 1 / radius_;
}

// * * * * * * * * * * * * Class implicitSlottedSphere * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSlottedSphere, false);
addToRunTimeSelectionTable(implicitSurface, implicitSlottedSphere, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

implicitSlottedSphere::implicitSlottedSphere(
    vector center, 
    scalar radius, 
    vector minSlotCorner, 
    vector maxSlotCorner
)
    : 
        center_(center), 
        radius_(radius), 
        minSlotCorner_(minSlotCorner),
        maxSlotCorner_(maxSlotCorner)
{}

implicitSlottedSphere::implicitSlottedSphere(dictionary const& dict)
    : 
        /*
        Enright, D., Fedkiw, R. P., Ferziger, J., & Mitchell, I. (2002). A Hybrid
        Particle Level Set Method for Improved Interface Capturing. Journal of
        Computational Physics (Vol. 183). https://doi.org/10.1006/jcph.2002.7166
        */
        center_(dict.getOrDefault<vector>("center", vector(0.5, 0.5, 0.75))),
        radius_(dict.getOrDefault<scalar>("radius", 0.15)), 
        minSlotCorner_(dict.getOrDefault<vector>("minSlotCorner", vector(0.2, 0.46, 0.2))),
        maxSlotCorner_(dict.getOrDefault<vector>("maxSlotCorner", vector(0.8, 0.54, 0.735)))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar implicitSlottedSphere::value(const vector& p) const
{
    vector e1 (1.,0.,0.);
    vector e2 (0.,1.,0.);
    vector e3 (0.,0.,1.);

    return 
    -Foam::max(
        Foam::mag(p - center_) - radius_,
        Foam::min(
            Foam::min(
                Foam::min(
                    (e1 &  (p - minSlotCorner_)), 
                    (-e1 & (p - maxSlotCorner_)) 
                ),
                Foam::min(
                    (e2 &  (p - minSlotCorner_)), 
                    (-e2 & (p - maxSlotCorner_)) 
                )
            ),
            Foam::min(
                (e3 &  (p - minSlotCorner_)), 
                (-e3 & (p - maxSlotCorner_)) 
            )
        )
    );
}

vector implicitSlottedSphere::grad(const vector& x) const
{
    return vector(NAN, NAN, NAN);
}

vector implicitSlottedSphere::center() const
{
    return center_; 
}

scalar implicitSlottedSphere::radius() const
{
    return radius_; 
}

vector implicitSlottedSphere::minSlotCorner() const
{
    return minSlotCorner_; 
}

vector implicitSlottedSphere::maxSlotCorner() const
{
    return maxSlotCorner_; 
}

scalar implicitSlottedSphere::curvature(const vector& x) const
{
    return scalar(NAN); 
}

// * * * * * * * * * * * * Class implicitEllipsoid * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitEllipsoid, false);
addToRunTimeSelectionTable(implicitSurface, implicitEllipsoid, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

implicitEllipsoid::implicitEllipsoid(vector center, vector axes)
    : 
        center_(center), 
        axes_(axes)
{
    setAxesSqr(axes);
}

implicitEllipsoid::implicitEllipsoid(const dictionary& dict)
    :
        center_(dict.get<vector>("center")),
        axes_(dict.get<vector>("axes"))
{
    setAxesSqr(axes_);
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void implicitEllipsoid::setAxesSqr(const vector& axes)
{
    axesSqr_[0] = Foam::sqr(axes[0]);
    axesSqr_[1] = Foam::sqr(axes[1]);
    axesSqr_[2] = Foam::sqr(axes[2]);
}

scalar implicitEllipsoid::value(const vector& x) const
{
    return Foam::sqr(x[0] - center_[0]) / axesSqr_[0] + 
           Foam::sqr(x[1] - center_[1]) / axesSqr_[1] + 
           Foam::sqr(x[2] - center_[2]) / axesSqr_[2] - 1;
}

vector implicitEllipsoid::grad(const vector& x) const
{
    return 2*vector
    (
        (x[0] - center_[0])/axesSqr_[0], 
        (x[1] - center_[1])/axesSqr_[1], 
        (x[2] - center_[2])/axesSqr_[2]
    );
}

vector implicitEllipsoid::center() const
{
    return center_; 
}

vector implicitEllipsoid::axes() const
{
    return axes_; 
}

scalar implicitEllipsoid::curvature(const vector& x) const
{
    // Total curvature (div(gradF/|gradF|), the kappa1 + kappa2 convention that
    // gives 2/r on a sphere's level set) of the LEVEL SET of the quadratic
    // form F = sum_i (x_i - O_i)^2/A_i - 1 THROUGH the point x, from the
    // standard identity
    //     div(n) = (lapF - n . HessF . n)/|gradF|,   n = gradF/|gradF|,
    // exact at every x with gradF != 0.
    //
    // HISTORY. This replaces a machine-generated SymPy transcription that the
    // verification gate (leiaTestSignedDistanceEllipsoid) proved wrong: on a
    // degenerate sphere it returned 1/r (half the div convention) and at the
    // ends of a triaxial ellipsoid's axes values consistent with neither the
    // mean nor the total curvature -- symptomatic of squared semi-axes being
    // substituted where the generated symbols expected plain ones. The 2D
    // droplet cases use axes = (a, b, 1) with x_z - O_z = 0; the z-axis then
    // contributes nothing to gradF and only 2/A_z = 2 to lapF, which enters as
    // 2/|gradF| ~ O(h/R * kappa * (b/R)) -- negligible for every configured
    // case, so the identity is dimension-consistent under the axes trick.
    scalar g2 = 0, lap = 0, nHn = 0;
    for (direction i = 0; i < 3; ++i)
    {
        const scalar A = axesSqr_[i];
        const scalar gi = 2.0*(x[i] - center_[i])/A;
        g2 += sqr(gi);
        lap += 2.0/A;
        nHn += sqr(gi)*(2.0/A);
    }
    const scalar gmag = Foam::sqrt(max(g2, VSMALL));
    return (lap - nHn/max(g2, VSMALL))/gmag;
}

// * * * * * * * * Class signedDistanceEllipsoid  * * * * * * * * * * * //

defineTypeNameAndDebug(signedDistanceEllipsoid, false);
addToRunTimeSelectionTable
(
    implicitSurface,
    signedDistanceEllipsoid,
    dictionary
);

signedDistanceEllipsoid::signedDistanceEllipsoid(vector center, vector axes)
:
    center_(center),
    axes_(axes),
    ellipsoid_(center, axes)
{}

signedDistanceEllipsoid::signedDistanceEllipsoid(const dictionary& dict)
:
    center_(dict.get<vector>("center")),
    axes_(dict.get<vector>("axes")),
    ellipsoid_(center_, axes_)
{}

vector signedDistanceEllipsoid::closestPoint(const vector& x) const
{
    // First octant: p_i = |x_i - c_i|; remember the signs to map q back.
    vector p, sgn;
    for (direction i = 0; i < 3; ++i)
    {
        const scalar d = x[i] - center_[i];
        sgn[i] = (d < 0) ? -1.0 : 1.0;
        p[i] = mag(d);
    }

    // F(t) = sum (a_i p_i/(t + a_i^2))^2 - 1, strictly decreasing.
    auto F = [&](const scalar t) -> scalar
    {
        scalar f = -1.0;
        for (direction i = 0; i < 3; ++i)
        {
            const scalar den = t + sqr(axes_[i]);
            f += sqr(axes_[i]*p[i]/den);
        }
        return f;
    };

    // Bracket the unique root. Left end: just right of -a_j^2 for the
    // SMALLEST axis with p_j > 0 (there F -> +inf); if x is the exact centre
    // every p_i = 0 and the nearest point is the end of the shortest axis.
    scalar aMinSqr = GREAT;
    bool anyP = false;
    for (direction i = 0; i < 3; ++i)
    {
        if (p[i] > SMALL*axes_[i])
        {
            anyP = true;
            aMinSqr = min(aMinSqr, sqr(axes_[i]));
        }
    }
    if (!anyP)
    {
        direction j = 0;
        for (direction i = 1; i < 3; ++i)
        {
            if (axes_[i] < axes_[j]) j = i;
        }
        vector q = center_;
        q[j] += axes_[j];
        return q;
    }

    scalar lo = -aMinSqr*(1.0 - 1e-12);
    // ensure F(lo) > 0 (roundoff at razor-thin margins)
    for (label k = 0; k < 60 && F(lo) <= 0; ++k)
    {
        lo = -aMinSqr + 0.5*(lo + aMinSqr);   // move closer to the pole
    }
    scalar hi = max(max(axes_[0], max(axes_[1], axes_[2]))*mag(p), aMinSqr);
    for (label k = 0; k < 200 && F(hi) >= 0; ++k)
    {
        hi = 2.0*hi + aMinSqr;
    }

    // Safeguarded Newton (bisection fallback) on the bracketed unique root.
    scalar t = 0.5*(lo + hi);
    for (label it = 0; it < 100; ++it)
    {
        const scalar f = F(t);
        if (f > 0) { lo = t; } else { hi = t; }
        scalar df = 0;
        for (direction i = 0; i < 3; ++i)
        {
            const scalar den = t + sqr(axes_[i]);
            df -= 2.0*sqr(axes_[i]*p[i])/(den*den*den);
        }
        scalar tN = (mag(df) > VSMALL) ? t - f/df : t;
        if (tN <= lo || tN >= hi)
        {
            tN = 0.5*(lo + hi);
        }
        if (mag(tN - t) < 1e-14*(mag(t) + aMinSqr)) { t = tN; break; }
        t = tN;
    }

    vector q;
    for (direction i = 0; i < 3; ++i)
    {
        q[i] = center_[i] + sgn[i]*sqr(axes_[i])*p[i]/(t + sqr(axes_[i]));
    }
    return q;
}

scalar signedDistanceEllipsoid::value(const vector& x) const
{
    const vector q = closestPoint(x);
    const scalar d = mag(x - q);
    return (ellipsoid_.value(x) >= 0) ? d : -d;
}

vector signedDistanceEllipsoid::grad(const vector& x) const
{
    const vector q = closestPoint(x);
    const scalar d = mag(x - q);
    if (d < SMALL*min(axes_[0], min(axes_[1], axes_[2])))
    {
        // On the surface: the SDF gradient is the outward unit normal.
        const vector g = ellipsoid_.grad(q);
        return g/max(mag(g), VSMALL);
    }
    const scalar s = (ellipsoid_.value(x) >= 0) ? 1.0 : -1.0;
    return s*(x - q)/d;
}

scalar signedDistanceEllipsoid::curvature(const vector& x) const
{
    // Total curvature (kappa1 + kappa2) of the SURFACE at the closest point --
    // the interface value the parallel-surface-inverted delivery targets, in
    // the div(n) convention that gives (nd-1)/R on a sphere, and the same
    // at-closest-point convention as signedDistanceEllipse in 2D.
    //
    // Computed DIRECTLY from the implicit form F = sum x_i^2/A_i - 1
    // (A_i = a_i^2) via the standard identity, evaluated AT q on the surface:
    //     kappa1 + kappa2 = div(gradF/|gradF|)
    //                     = ( lapF - n . HessF . n ) / |gradF|,  n = gradF/|gradF|
    // with gradF = (2 q_i/A_i), HessF = diag(2/A_i), lapF = sum 2/A_i.
    // Verified limits: sphere -> 2/R; a-axis end -> a (1/b^2 + 1/c^2).
    //
    // Deliberately NOT delegated to implicitEllipsoid::curvature: the
    // verification gate (leiaTestSignedDistanceEllipsoid) shows that formula
    // returns 1/R on a degenerate sphere (half the div convention) and values
    // at triaxial axis ends consistent with neither H nor 2H -- its
    // machine-generated transcription is suspect and is kept untouched only
    // because recorded 2D ellipsoid-psi baselines were measured against it.
    const vector q = closestPoint(x);
    scalar g2 = 0, lap = 0, nHn = 0;
    for (direction i = 0; i < 3; ++i)
    {
        const scalar A = sqr(axes_[i]);
        const scalar gi = 2.0*(q[i] - center_[i])/A;
        g2 += sqr(gi);
        lap += 2.0/A;
        nHn += sqr(gi)*(2.0/A);
    }
    const scalar gmag = Foam::sqrt(max(g2, VSMALL));
    return (lap - nHn/g2)/gmag;
}

// * * * * * * * * * * Class signedDistanceEllipse * * * * * * * * * * * //

defineTypeNameAndDebug(signedDistanceEllipse, false);
addToRunTimeSelectionTable
(
    implicitSurface,
    signedDistanceEllipse,
    dictionary
);

signedDistanceEllipse::signedDistanceEllipse(vector center, vector axes)
:
    center_(center),
    axes_(axes)
{}

signedDistanceEllipse::signedDistanceEllipse(const dictionary& dict)
:
    center_(dict.get<vector>("center")),
    axes_(dict.get<vector>("axes"))
{}

scalar signedDistanceEllipse::closestParameter(const vector& x) const
{
    const scalar px = mag(x.x() - center_.x());
    const scalar py = mag(x.y() - center_.y());
    const scalar a = axes_.x();
    const scalar b = axes_.y();
    const scalar halfPi = constant::mathematical::pi/2;

    auto distanceSqr = [&](const scalar theta)
    {
        return
            sqr(a*Foam::cos(theta) - px)
          + sqr(b*Foam::sin(theta) - py);
    };

    // Bracket the global minimum first.  The distance to a point inside an
    // ellipse is not guaranteed to be unimodal over the complete quadrant;
    // this inexpensive scan avoids converging to the wrong stationary point.
    const label nBracket = 64;
    const scalar dTheta = halfPi/nBracket;
    label best = 0;
    scalar bestDistance = distanceSqr(0);
    for (label i = 1; i <= nBracket; ++i)
    {
        const scalar candidate = distanceSqr(i*dTheta);
        if (candidate < bestDistance)
        {
            best = i;
            bestDistance = candidate;
        }
    }

    scalar left = max(scalar(0), (best - 1)*dTheta);
    scalar right = min(halfPi, (best + 1)*dTheta);
    const scalar golden = 0.5*(Foam::sqrt(5.0) - 1);
    scalar x1 = right - golden*(right - left);
    scalar x2 = left + golden*(right - left);
    scalar f1 = distanceSqr(x1);
    scalar f2 = distanceSqr(x2);
    for (label iter = 0; iter < 48; ++iter)
    {
        if (f1 > f2)
        {
            left = x1;
            x1 = x2;
            f1 = f2;
            x2 = left + golden*(right - left);
            f2 = distanceSqr(x2);
        }
        else
        {
            right = x2;
            x2 = x1;
            f2 = f1;
            x1 = right - golden*(right - left);
            f1 = distanceSqr(x1);
        }
    }
    return 0.5*(left + right);
}

scalar signedDistanceEllipse::value(const vector& x) const
{
    const scalar theta = closestParameter(x);
    const scalar px = mag(x.x() - center_.x());
    const scalar py = mag(x.y() - center_.y());
    const scalar qx = axes_.x()*Foam::cos(theta);
    const scalar qy = axes_.y()*Foam::sin(theta);
    const scalar distance = Foam::sqrt(sqr(px - qx) + sqr(py - qy));
    const scalar algebraic =
        sqr(px/axes_.x()) + sqr(py/axes_.y()) - 1;
    return algebraic >= 0 ? distance : -distance;
}

vector signedDistanceEllipse::grad(const vector& x) const
{
    const scalar theta = closestParameter(x);
    vector normal
    (
        Foam::cos(theta)/axes_.x(),
        Foam::sin(theta)/axes_.y(),
        0
    );
    normal /= mag(normal) + VSMALL;
    normal.x() *= (x.x() < center_.x() ? -1 : 1);
    normal.y() *= (x.y() < center_.y() ? -1 : 1);
    return normal;
}

scalar signedDistanceEllipse::curvature(const vector& x) const
{
    const scalar theta = closestParameter(x);
    const scalar a = axes_.x();
    const scalar b = axes_.y();
    return a*b/Foam::pow
    (
        sqr(a*Foam::sin(theta)) + sqr(b*Foam::cos(theta)),
        1.5
    );
}

vector signedDistanceEllipse::center() const
{
    return center_;
}

vector signedDistanceEllipse::axes() const
{
    return axes_;
}

// * * * * * * * * * * * * Class implicitSinc * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSinc, false);
addToRunTimeSelectionTable(implicitSurface, implicitSinc, dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

implicitSinc::implicitSinc(vector origin, scalar amplitude, scalar omega)
    : 
        origin_(origin), 
        amplitude_(amplitude), 
        omega_(omega)
{}

implicitSinc::implicitSinc(const dictionary& dict)
    :
        origin_(dict.get<vector>("origin")),
        amplitude_(dict.get<scalar>("amplitude")),
        omega_(dict.get<scalar>("omega"))
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar implicitSinc::value(const vector& x) const
{
    double r = Foam::sqrt
    (
        Foam::sqr(x[0] - origin_[0]) + 
        Foam::sqr(x[1] - origin_[1]) 
    );

    if (r < std::numeric_limits<double>::min())
        return amplitude_; 
    else 
    {
        return x[2] - origin_[2] - amplitude_ * sin(omega_ * r) / (omega_ * r);
    }
}

vector implicitSinc::grad(const vector& x) const
{
    const scalar& A = amplitude_; 
    const scalar& O0 = origin_[0];
    const scalar& O1 = origin_[1];

    const scalar& x0 = x[0];
    const scalar& x1 = x[1];

    return vector // Expression calculated in sympy.
    (
        A*(O0 - x0)*(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)),
            A*(O1 - x1)*(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)),
            1
        );
}

vector implicitSinc::origin() const
{
    return origin_; 
}

scalar implicitSinc::amplitude() const
{
    return amplitude_; 
}

scalar implicitSinc::omega() const
{
    return omega_; 
}

scalar implicitSinc::curvature(const vector& x) const
{
    const scalar& x0 = x[0];
    const scalar& x1 = x[1];

    const scalar& O0 = origin_[0];
    const scalar& O1 = origin_[1];

    const scalar& A = amplitude_;
    const scalar& w = omega_;

    // Taken from sympy
    return 0.5*A*(O0 - x0)*(5*O0 - 5*x0)*(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 7.0/2.0)*sqrt(pow(A, 2)*pow(O0 - x0, 2)*pow(w
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)*pow(w
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + 1)) + 0.5*A*(O0 - x0)*(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))*(-1.0/2.0*pow(A, 2)*(-2*O0 + 2*x0)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) - 1.0/2.0*pow(A, 2)*pow(O0 - x0, 2)*(10*O0 - 10*x0)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 6)) - 1.0/2.0*pow(A, 2)*pow(O0 - x0, 2)*(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))*(-2*pow(w, 2)*(-O0 + x0)*(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-3*O0 + 3*x0)*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-O0 + x0)*(-pow(O0 - x0, 2)
           - pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           + 2*(2*O0 - 2*x0)*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) - 1.0/2.0*pow(A, 2)*(10*O0 - 10*x0)*pow(O1 - x1, 2)*pow(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 6)) - 1.0/2.0*pow(A, 2)*pow(O1 - x1, 2)*(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))*(-2*pow(w, 2)*(-O0 + x0)*(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           *sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-3*O0 + 3*x0)*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-O0 + x0)*(-pow(O0 - x0, 2) - pow(O1 - x1, 2))*cos(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)) + 2*(2*O0 - 2*x0)*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)))/(w
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)*pow(pow(A, 2)*pow(O0 - x0, 2)*pow(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5))
           + 1, 3.0/2.0)) + 0.5*A*(O0 - x0)*(-pow(w, 2)*(-O0 + x0)*(pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + w*(-3*O0 + 3*x0)*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + w*(-O0 + x0)*(-pow(O0 - x0, 2) - pow(O1 - x1, 2))
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)) + (2*O0 - 2*x0)
           *sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)
           *sqrt(pow(A, 2)*pow(O0 - x0, 2)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + 1)) + 0.5*A*(O1 - x1)*(5*O1 - 5*x1)*(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 7.0/2.0)*sqrt(pow(A, 2)*pow(O0 - x0, 2)*pow(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + 1))
           + 0.5*A*(O1 - x1)*(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))
           *(-1.0/2.0*pow(A, 2)*pow(O0 - x0, 2)*(10*O1 - 10*x1)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 6))
           - 1.0/2.0*pow(A, 2)*pow(O0 - x0, 2)*(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))))*(-2*pow(w, 2)*(-O1 + x1)*(pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))) + 2*w*(-3*O1 + 3*x1)*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))) + 2*w*(-O1 + x1)*(-pow(O0 - x0, 2) - pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)) + 2*(2*O1 - 2*x1)*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))))/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) - 1.0/2.0*pow(A, 2)*(-2*O1 + 2*x1)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) - 1.0/2.0*pow(A, 2)*pow(O1 - x1, 2)*(10*O1 - 10*x1)*pow(w
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 6)) - 1.0/2.0*pow(A, 2)*pow(O1 - x1, 2)*(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))*(-2*pow(w, 2)*(-O1 + x1)*(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-3*O1 + 3*x1)*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + 2*w*(-O1 + x1)*(-pow(O0 - x0, 2)
           - pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           + 2*(2*O1 - 2*x1)*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)))/(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)*pow(pow(A, 2)*pow(O0 - x0, 2)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)*pow(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) + 1, 3.0/2.0)) + 0.5*A*(O1 - x1)*(-pow(w, 2)*(-O1 + x1)*(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + w*(-3*O1 + 3*x1)*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) + w*(-O1 + x1)*(-pow(O0 - x0, 2)
           - pow(O1 - x1, 2))*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))/sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))
           + (2*O1 - 2*x1)*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5.0/2.0)*sqrt(pow(A, 2)*pow(O0 - x0, 2)*pow(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2)
           + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w
           *sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + 1))
           - 1.0*A*(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(w*pow(pow(O0 - x0, 2)
           + pow(O1 - x1, 2), 5.0/2.0)*sqrt(pow(A, 2)*pow(O0 - x0, 2)*pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)
           *cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2)
           + pow(O1 - x1, 2))), 2)/(pow(w, 2)*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + pow(A, 2)*pow(O1 - x1, 2)
           *pow(w*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2)))
           - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(w*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))), 2)/(pow(w, 2)
           *pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5)) + 1));
}

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
