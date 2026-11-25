/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "rigidMotionProfileFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvSolution.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::point Foam::rigidMotionProfileFvPatchScalarField::origin
(
    const point& x,
    const scalar t
) const
{
    const scalar w = mag(omega_);
    if (w < VSMALL) { return x - velocity_*t; }

    const vector k(omega_/w);
    const vector v(x - point_);
    const scalar th = -w*t;
    return
        point_ + v*Foam::cos(th) + (k ^ v)*Foam::sin(th)
      + k*(k & v)*(1 - Foam::cos(th));
}


Foam::tmp<Foam::scalarField>
Foam::rigidMotionProfileFvPatchScalarField::exact() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    if (!profile_)
    {
        const dictionary& lsDict = mesh.solutionDict().subDict("levelSet");
        const dictionary& surfDict = lsDict.subDict("implicitSurface");
        surface_ = implicitSurface::New(surfDict.get<word>("type"), surfDict);
        profile_ = profile::New(mesh, surface_.ref());
    }

    const scalar t = mesh.time().value();
    const vectorField& Cf = patch().Cf();
    auto tv = tmp<scalarField>::New(Cf.size());
    scalarField& v = tv.ref();
    forAll(v, i) { v[i] = profile_->value(origin(Cf[i], t)); }
    return tv;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidMotionProfileFvPatchScalarField::
rigidMotionProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    velocity_(Zero),
    omega_(Zero),
    point_(Zero)
{}


Foam::rigidMotionProfileFvPatchScalarField::
rigidMotionProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::LAZY_READ),
    velocity_(dict.getOrDefault<vector>("velocity", Zero)),
    omega_(dict.getOrDefault<vector>("omega", Zero)),
    point_(dict.getOrDefault<vector>("point", Zero))
{
    if (mag(velocity_) > VSMALL && mag(omega_) > VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "only one of velocity and omega may be non-zero"
            << exit(FatalIOError);
    }
    // A dictionary with a value (a written field, or a field that a mesh
    // redistribution sends between ranks) keeps that value: the mesh of a
    // redistribution step has no fvSolution of the case. A dictionary without a
    // value (the initial field of a case) gets the exact profile.
    if (!dict.found("value"))
    {
        fvPatchField<scalar>::operator=(exact());
    }
}


Foam::rigidMotionProfileFvPatchScalarField::
rigidMotionProfileFvPatchScalarField
(
    const rigidMotionProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    velocity_(ptf.velocity_),
    omega_(ptf.omega_),
    point_(ptf.point_)
{}


Foam::rigidMotionProfileFvPatchScalarField::
rigidMotionProfileFvPatchScalarField
(
    const rigidMotionProfileFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    velocity_(ptf.velocity_),
    omega_(ptf.omega_),
    point_(ptf.point_)
{}


Foam::rigidMotionProfileFvPatchScalarField::
rigidMotionProfileFvPatchScalarField
(
    const rigidMotionProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    velocity_(ptf.velocity_),
    omega_(ptf.omega_),
    point_(ptf.point_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::rigidMotionProfileFvPatchScalarField::updateCoeffs()
{
    if (updated()) { return; }
    fvPatchField<scalar>::operator=(exact());
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::rigidMotionProfileFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("velocity", velocity_);
    os.writeEntry("omega", omega_);
    os.writeEntry("point", point_);
    fvPatchField<scalar>::writeValueEntry(os);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rigidMotionProfileFvPatchScalarField
    );
}

// ************************************************************************* //
