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

\*---------------------------------------------------------------------------*/

#include "interfaceMarkers.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "interpolationCellPoint.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceMarkers, 0);
    addToRunTimeSelectionTable(functionObject, interfaceMarkers, dictionary);
}
}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

Foam::vector Foam::functionObjects::interfaceMarkers::velocityAt
(
    const point& x,
    const scalar t
) const
{
    const scalar f =
        velModel_->isOscillating()
      ? velModel_->oscillationFactor(t, velModel_->tau())
      : scalar(1);
    return f*velModel_->velocity(x);
}


void Foam::functionObjects::interfaceMarkers::advance(const scalar t1)
{
    const scalar T = t1 - tMarkers_;
    if (T <= 0)
    {
        return;
    }
    const scalar dt = T/nSubSteps_;
    const bool twoD = (mesh_.nGeometricD() == 2);
    forAll(markers_, l)
    {
        point& x = markers_[l];
        scalar t = tMarkers_;
        for (label s = 0; s < nSubSteps_; ++s)
        {
            const vector k1 = velocityAt(x, t);
            const vector k2 = velocityAt(x + 0.5*dt*k1, t + 0.5*dt);
            const vector k3 = velocityAt(x + 0.5*dt*k2, t + 0.5*dt);
            const vector k4 = velocityAt(x + dt*k3, t + dt);
            x += dt/6*(k1 + 2*k2 + 2*k3 + k4);
            t += dt;
        }
        if (twoD)
        {
            x.z() = mesh_.bounds().centre().z();
        }
    }
    tMarkers_ = t1;
}


void Foam::functionObjects::interfaceMarkers::sampleAndWrite()
{
    const volScalarField& psi = lookupObject<volScalarField>(fieldName_);
    const volVectorField gradPsi(fvc::grad(psi, gradScheme_));
    interpolationCellPoint<scalar> iPsi(psi);
    interpolationCellPoint<vector> iGrad(gradPsi);

    scalar nFound = 0, q1 = 0, q2 = 0, qs = 0, p1 = 0, p2 = 0;
    forAll(markers_, l)
    {
        const point& x = markers_[l];
        const label celli = mesh_.findCell(x);
        if (celli < 0)
        {
            continue;
        }
        const scalar qh = Foam::mag(iGrad.interpolate(x, celli));
        const scalar ph = iPsi.interpolate(x, celli);
        const scalar eq = qh - q0_[l];
        const scalar ep = ph/q0_[l];
        nFound += 1;
        q1 += Foam::mag(eq); q2 += eq*eq; qs += eq;
        p1 += Foam::mag(ep); p2 += ep*ep;
    }
    reduce(nFound, sumOp<scalar>());
    reduce(q1, sumOp<scalar>()); reduce(q2, sumOp<scalar>());
    reduce(qs, sumOp<scalar>());
    reduce(p1, sumOp<scalar>()); reduce(p2, sumOp<scalar>());

    const scalar n = Foam::max(nFound, scalar(1));
    Log << type() << ' ' << name() << ": t " << tMarkers_
        << " markers " << label(nFound) << '/' << markers_.size()
        << " E_Q L1 " << q1/n << " L2 " << Foam::sqrt(q2/n)
        << " E_POS L1 " << p1/n << " L2 " << Foam::sqrt(p2/n) << endl;

    if (Pstream::master())
    {
        if (!csv_)
        {
            // globalPath: the case directory also in parallel (path() is processor0)
            csv_.reset(new OFstream(time_.globalPath()/"interfaceMarkers.csv"));
            csv_()
                << "TIME,N_MARKERS,N_FOUND,E_Q_L1,E_Q_L2,E_Q_MEAN,"
                << "E_POS_L1,E_POS_L2,Q_EXACT_MIN,Q_EXACT_MAX\n";
        }
        csv_().precision(12);
        csv_()
            << tMarkers_ << ',' << markers_.size() << ',' << label(nFound)
            << ',' << q1/n << ',' << Foam::sqrt(q2/n) << ',' << qs/n
            << ',' << p1/n << ',' << Foam::sqrt(p2/n)
            << ',' << Foam::min(q0_) << ',' << Foam::max(q0_) << nl;
        csv_().flush();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceMarkers::interfaceMarkers
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_("psi"),
    nMarkers_(720),
    nSubSteps_(4),
    gradScheme_("gradPsiMetric"),
    tMarkers_(runTime.value())
{
    read(dict);

    velModel_ = velocityModel::New(mesh_);

    const fvSolution& fvSol(mesh_);
    const dictionary& sDict =
        fvSol.subDict("levelSet").subDict("implicitSurface");
    surface_ = implicitSurface::New(sDict.get<word>("type"), sDict);
    const vector c = sDict.get<vector>("center");
    const scalar R = sDict.get<scalar>("radius");

    markers_.setSize(nMarkers_);
    if (mesh_.nGeometricD() == 2)
    {
        if (mesh_.geometricD()[vector::Z] != -1)
        {
            FatalErrorInFunction
                << "interfaceMarkers: a 2D mesh must be empty in z"
                << exit(FatalError);
        }
        const scalar zc = mesh_.bounds().centre().z();
        forAll(markers_, l)
        {
            const scalar th =
                constant::mathematical::twoPi*(l + 0.5)/nMarkers_;
            markers_[l] =
                point(c.x() + R*Foam::cos(th), c.y() + R*Foam::sin(th), zc);
        }
    }
    else
    {
        // Fibonacci sphere: nearly uniform markers, no structure assumed.
        const scalar ga =
            constant::mathematical::pi*(3 - Foam::sqrt(scalar(5)));
        forAll(markers_, l)
        {
            const scalar y = 1 - 2*(l + 0.5)/nMarkers_;
            const scalar r = Foam::sqrt(Foam::max(1 - y*y, scalar(0)));
            const scalar ph = ga*l;
            markers_[l] =
                c + R*vector(r*Foam::cos(ph), y, r*Foam::sin(ph));
        }
    }

    q0_.setSize(nMarkers_);
    forAll(markers_, l)
    {
        q0_[l] = Foam::mag(surface_->grad(markers_[l]));
    }

    Info<< type() << ' ' << name << ": " << nMarkers_ << " markers on the "
        << sDict.get<word>("type") << " (R " << R << "), q0 in ["
        << Foam::min(q0_) << ", " << Foam::max(q0_) << "], RK4 with "
        << nSubSteps_ << " sub-steps per step, grad " << gradScheme_ << endl;

    // The row of the start time: the interpolation error of the initial field.
    if (foundObject<volScalarField>(fieldName_))
    {
        sampleAndWrite();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceMarkers::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    fieldName_ = dict.getOrDefault<word>("field", "psi");
    nMarkers_ = dict.getOrDefault<label>("nMarkers", 720);
    nSubSteps_ = dict.getOrDefault<label>("nSubSteps", 4);
    gradScheme_ = dict.getOrDefault<word>("gradScheme", "gradPsiMetric");
    return true;
}


bool Foam::functionObjects::interfaceMarkers::execute()
{
    advance(time_.value());
    sampleAndWrite();
    return true;
}


bool Foam::functionObjects::interfaceMarkers::write()
{
    return true;
}


// ************************************************************************* //
