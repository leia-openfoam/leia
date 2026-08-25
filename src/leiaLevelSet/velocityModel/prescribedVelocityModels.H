/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Concrete prescribed (verification) velocity models selected at runtime by
    velocityModel::New: shear2D, deformation3D, shear3D, translation, rotation,
    vortex2D, periodic2D, uniaxialStrain. Velocity formulae are unchanged from
    the original advectionVerification.H implementation, except uniaxialStrain,
    which is new and is the only NON-SOLENOIDAL model here (div v = strainRate).

SourceFiles
    prescribedVelocityModels.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedVelocityModels_H
#define prescribedVelocityModels_H

#include "velocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- 2D shear (single-vortex) flow
class shear2D
:
    public velocityModel
{
public:
    TypeName("shear2D");
    explicit shear2D(const fvMesh& mesh) : velocityModel(mesh) {}
    virtual ~shear2D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 3D deformation flow (LeVeque)
class deformation3D
:
    public velocityModel
{
public:
    TypeName("deformation3D");
    explicit deformation3D(const fvMesh& mesh) : velocityModel(mesh) {}
    virtual ~deformation3D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 3D shear flow (Liovic et al. 2006)
class shear3D
:
    public velocityModel
{
    scalar R_, Umax_, x0_, y0_;
public:
    TypeName("shear3D");
    explicit shear3D(const fvMesh& mesh);
    virtual ~shear3D() = default;
    virtual vector velocity(const vector& p) const;
};


//- Uniform translation
class translation
:
    public velocityModel
{
    vector velocity_;
public:
    TypeName("translation");
    explicit translation(const fvMesh& mesh);
    virtual ~translation() = default;
    virtual vector velocity(const vector& p) const;
};


//- Rigid-body rotation
class rotation
:
    public velocityModel
{
    vector point_;
    vector omega_;
public:
    TypeName("rotation");
    explicit rotation(const fvMesh& mesh);
    virtual ~rotation() = default;
    virtual vector velocity(const vector& p) const;
};


//- 2D single vortex
class vortex2D
:
    public velocityModel
{
    scalar v0_;
public:
    TypeName("vortex2D");
    explicit vortex2D(const fvMesh& mesh);
    virtual ~vortex2D() = default;
    virtual vector velocity(const vector& p) const;
};


//- 2D periodic shear
class periodic2D
:
    public velocityModel
{
    scalar v0_;
    scalar c1_;
    scalar c2_;
public:
    TypeName("periodic2D");
    explicit periodic2D(const fvMesh& mesh);
    virtual ~periodic2D() = default;
    virtual vector velocity(const vector& p) const;
};


//- Uniaxial stretching along x about a stagnation point x0:
//
//      v(x) = a (x - x0) . e_x  e_x,    a = strainRate [1/s], default 1,
//
//  so that d(v_x)/dx = a identically and v_y = v_z = 0.
//
//  Dictionary (the top-level fvSolution sub-dictionary `velocityModel`):
//
//      velocityModel
//      {
//          type            uniaxialStrain;
//          strainRate      1;          // a [1/s]
//          stagnationPoint (0 0 0);    // x0, v(x0) = 0
//          oscillation     off;        // see OSCILLATION below
//      }
//
//  WHY THIS MODEL EXISTS, AND WHY IT IS NOT SOLENOIDAL. Every other
//  prescribed field in this file is divergence-free. This one is not:
//
//      div v = a  !=  0.
//
//  That is intrinsic to one dimension, not a defect of the construction: the
//  only strain a purely 1D flow possesses is the normal strain, and there is
//  no tangential direction in which a compensating contraction could live, so
//  a non-trivial solenoidal 1D flow does not exist. Any case or document
//  built on this model must state div v = a explicitly.
//  The level-set equation is nevertheless assembled correctly, because it is
//  assembled in ADVECTIVE form,
//
//      fvm::ddt(psi) + fvm::div(phi, psi) - fvm::Sp(fvc::div(phi), psi) == S
//
//  (eulerianAdvection::advance, leiaLevelSetTwoPhaseFoam/alphaEqn.H), and
//  div(v psi) - psi div(v) = v . grad(psi) is an identity for ANY v, not only
//  a solenoidal one. So the equation actually solved is
//  d(psi)/dt + v . grad(psi) = S regardless of div v. It is the
//  -fvm::Sp(divPhi, psi) term that does that work here; with every other
//  (solenoidal) model in this file that term is zero and is never exercised.
//  Do NOT run a case built on this model with the solver option
//  -fluxCorrection: correctFlux() projects phi onto div(phi) = 0 and would
//  delete the entire test.
//
//  INTERFACE MOTION. Initialising the interface at x_i(0) = xi0, the
//  characteristics of dx/dt = a (x - x0) give, with and without any source,
//
//      x_i(t) = x0 + (xi0 - x0) e^{a t}.
//
//  EXACT SOLUTIONS of d(psi)/dt + v . grad(psi) = S for psi0(x) = x - xi0,
//  which is the exact 1D signed distance:
//
//  (1) S = 0 (sdplsSource type noSource). psi is constant along
//      characteristics, so
//
//          psi(x,t) = x0 + (x - x0) e^{-a t} - xi0,   d(psi)/dx = e^{-a t}.
//
//      The signed-distance property DECAYS EXPONENTIALLY. At a = 1 1/s and
//      t = 1 s, d(psi)/dx = e^{-1} = 0.3678794412, i.e.
//      |grad psi| - 1 = -0.6321205588. The zero contour is untouched: the
//      interface still sits exactly at x_i(t). Losing the distance property
//      does not move the interface, which is exactly why a zero-contour-only
//      (shape or volume) error metric cannot see this failure, and why the
//      gradient error must always be reported next to shape and volume.
//
//  (2) S = a_strain psi, the SDPLS source R. Then d(psi)/dt = a psi along
//      characteristics and
//
//          psi(x,t) = x - x_i(t) = x - x0 - (xi0 - x0) e^{a t},
//
//      hence d(psi)/dx = 1 EXACTLY, for all t. The source holds the
//      signed-distance property exactly. This is the verification target: the
//      measured band |grad psi| must stay at 1 to discretization error while
//      the sourceless arm decays as e^{-a t}, and BOTH arms must track the
//      same x_i(t).
//
//  THE STRAIN THE SOURCE MUST SUPPLY. With psi0 the exact signed distance,
//  nhat = grad(psi)/|grad psi| = e_x and grad(v) = a e_x (x) e_x, so
//
//      a_strain = nhat . grad(v) . nhat = a,
//
//  CONSTANT IN SPACE AND TIME and equal to the dictionary value exactly
//  (1 1/s by default). This is what makes the 1D uniaxial-strain case the
//  cleanest available gate for the SDPLS source: the coefficient the source
//  must reproduce is an exact constant, so any deviation measured in a run is
//  discretization error, never a modelling error.
//
//  UNTESTED PREDICTION for the conservative source form Rdiv. Its derivation
//  rests on the identity a_strain = -div_s(v), which its own header states
//  holds for INCOMPRESSIBLE flow. In 1D the interface is a point: there are
//  no tangential directions, so div_s(v) = 0 identically while the true
//  strain is a, and the identity fails. Concretely, w = (nhat . v) nhat = v
//  exactly here because v is entirely normal, so phiW = phi, the two fvm::div
//  terms cancel DISCRETELY as well as at the continuum level, and
//  div(w) - a_strain = div(v) - a = 0: Rdiv assembles the ZERO MATRIX. The
//  PREDICTION is therefore that Rdiv behaves identically to noSource, with
//  |grad psi| decaying as e^{-a t}. Recorded here as a prediction made before
//  the measurement exists; it is NOT an established result.
//
//  DISCRETE FLUX, AND EMPTY DIRECTIONS. velocityModel::setVolumetricFlux
//  evaluates velocity(Cf) & Sf at each face centroid; it does NOT interpolate
//  U to the faces. Because v is linear in x and the face centroid is the
//  area-weighted mean position of a planar face, that product is the EXACT
//  face integral of v . n dA. Summing over the faces of a cell therefore
//  gives fvc::div(phi) = a to round-off on any mesh of planar faces, so the
//  -fvm::Sp(divPhi, psi) coefficient carries no discretization error of its
//  own. On a mesh whose y and z directions are `empty`, Sf on the x faces is
//  exactly +/- A e_x and Sf elsewhere is orthogonal to e_x, and v_y = v_z = 0
//  contribute nothing anywhere: the purely x-directed field needs no special
//  handling for the empty directions, in either setVelocity or
//  setVolumetricFlux.
//
//  OSCILLATION. The velocityModel base class defaults `oscillation` to "on",
//  which rescales U and phi by cos(pi t / tau) at every step for the reversed
//  (there-and-back) tests. That factor would invalidate every formula above,
//  so this model flips the DEFAULT to "off" in its constructor. An explicit
//  `oscillation on;` in the dictionary is still honoured, and no other
//  model's behaviour changes.
class uniaxialStrain
:
    public velocityModel
{
    //- Strain rate a [1/s]. d(v_x)/dx = a identically, so div v = a.
    scalar strainRate_;

    //- Stagnation point x0, v(x0) = 0. Only its x component enters the
    //  velocity; y and z are carried so a case can record the full point.
    vector stagnationPoint_;

public:
    TypeName("uniaxialStrain");
    explicit uniaxialStrain(const fvMesh& mesh);
    virtual ~uniaxialStrain() = default;
    virtual vector velocity(const vector& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
