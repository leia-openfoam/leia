/*---------------------------------------------------------------------------*\
    leia :: davofAdvectorFoam :: alphaSuSp.H: the zero sources Sp, Su, divU of the pure advection
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    From the isoAdvector benchmark solver advectorVoF (Johan Roenby, DHI,
    2016; TwoPhaseFlow apps/benchmark/advectorVoF): J. Roenby, H. Bredmose,
    H. Jasak, A computational method for sharp interface advection, R. Soc.
    Open Sci. 3 (2016) 160405, doi:10.1098/rsos.160405
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless/dimTime, Zero)
);

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phi + mesh.phi())
  : fvc::div(phi)
);
