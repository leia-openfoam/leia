/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::rigidMotionProfileFvPatchScalarField

Description
    The exact level set of a rigidly moving implicit surface as the patch value.

    The value at the face centre x and the time t is profile(X), with X the
    position at t = 0 of the material point that is at x at the time t:
        translation   X = x - velocity t
        rotation      X = point + R(-|omega| t)(x - point), axis omega/|omega|
    The profile and the implicit surface are those of levelSet in fvSolution,
    the same objects that leiaSetFields uses.

    USE. Exact-solution gates of the transport (translation, rotation). The
    semi-Lagrangian fit reads the patch value as a datum; with
    `stencilBoundaryFaces inflowOnly` only the inflow faces give data.
    zeroGradient gives the wrong datum there for a field with a gradient
    normal to the patch (MEASURED 2026-11-23, STATUS 11.44).

    Only one of velocity and omega may be non-zero.

Usage
    \verbatim
    left
    {
        type      rigidMotionProfile;
        velocity  (1 0 0);        // default (0 0 0)
        omega     (0 0 0);        // default (0 0 0)
        point     (0 0 0);        // default (0 0 0)
    }
    \endverbatim

SourceFiles
    rigidMotionProfileFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rigidMotionProfileFvPatchScalarField_H
#define rigidMotionProfileFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "implicitSurfaces.H"
#include "profile.H"

namespace Foam
{

class rigidMotionProfileFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    vector velocity_;
    vector omega_;
    vector point_;

    //- Built at the first use.
    mutable autoPtr<implicitSurface> surface_;
    mutable autoPtr<profile> profile_;

    //- Position at t = 0 of the point that is at x at the time t.
    point origin(const point& x, const scalar t) const;

    //- The exact values at the face centres at the present time.
    tmp<scalarField> exact() const;

public:

    TypeName("rigidMotionProfile");

    rigidMotionProfileFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    rigidMotionProfileFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    rigidMotionProfileFvPatchScalarField
    (
        const rigidMotionProfileFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    rigidMotionProfileFvPatchScalarField
    (
        const rigidMotionProfileFvPatchScalarField&
    );

    rigidMotionProfileFvPatchScalarField
    (
        const rigidMotionProfileFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return fvPatchField<scalar>::Clone(*this);
    }

    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return fvPatchField<scalar>::Clone(*this, iF);
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
