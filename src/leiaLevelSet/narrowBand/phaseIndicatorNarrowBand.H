/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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

Class
    Foam::phaseIndicatorNarrowBand

Description
    A narrow band seeded from the GEOMETRIC PHASE INDICATOR alpha rather than
    from psi, then dilated by owner--neighbour connectivity.

    A cell enters the seed when the phase indicator JUMPS across one of its
    faces -- the discrete statement of snGrad(alpha) != 0 -- or when the cell is
    itself partially filled, alphaTol < alpha_c < 1 - alphaTol, which is the
    geometric definition of a cut cell and catches a cell clipped at a corner
    whose face neighbours happen to agree. The union of the two is taken.

WHY THIS EXISTS, and what is wrong with the psi-based band it replaces
----------------------------------------------------------------------
The extension band is currently
    band = { |psi_c| <= (nLayers + 1) * cellSize_c }
(interfaceExtension.C). That is a DISTANCE band only if |grad psi| = 1 -- which
is precisely the property the velocity extension exists to maintain. The band
definition therefore depends on the quantity the machinery using it is supposed
to preserve, and when that degrades the band is no longer the region it claims
to be. It is circular by construction.

This is measured, not a concern in principle. On cases/1Dstretch (uniaxial
stretching, exact solution known, psi flattens as e^{-a t} without a source),
sweeping the extension width with the MEASUREMENT window held fixed at
|psi| <= 3h gave band mean |d(psi)/dx| against the exact 1:

    nLayers      3          8         32        200
    |dpsi/dx|  0.942131   0.999111  1.000000  1.000000

The ENTIRE gradient error was the fade region at the band edge -- a boundary
that lives in psi-space and therefore MOVES outward in x as psi flattens.
Remove it and the model reproduces the rigid translation exactly.

alpha is a geometric volume fraction produced by the Detrixhe--Aslam phase
indicator from a reconstructed interface. It does not depend on |grad psi| at
all, so a band seeded from it stays the region it claims to be no matter what
the level set does between redistancings. That is the whole point.

PARALLEL
--------
Two failure modes are avoided, both of which have already bitten this code base:

  (1) alpha is SYNCHRONISED before any face difference is taken, on a COPY, in
      the same way narrowBand::psiSynced() handles psi. Differencing against a
      stale coupled-patch value makes the seed depend on the decomposition.

  (2) The dilation is inherited from neighbourNarrowBand and NOT reimplemented.
      It exchanges once per layer, because layer k+1 on this rank may only be
      reachable through a cell that layer k placed on the far side of a
      processor boundary. Reimplementing that loop here would let the two copies
      drift apart.

Both mean the band field must be identical for a given mesh regardless of how
many ranks it is decomposed onto, which is the verification below.

Usage
    \verbatim
    narrowBand
    {
        type      phaseIndicator;
        alpha     alpha;      // registered indicator field, default "alpha"
        alphaTol  1e-8;       // jump / partial-fill tolerance, default 1e-8
        ntimes    3;          // dilation layers (neighbourNarrowBand)
    }
    \endverbatim

SourceFiles
    phaseIndicatorNarrowBand.C

\*---------------------------------------------------------------------------*/

#ifndef phaseIndicatorNarrowBand_H
#define phaseIndicatorNarrowBand_H

#include "neighbourNarrowBand.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class phaseIndicatorNarrowBand Declaration
\*---------------------------------------------------------------------------*/

class phaseIndicatorNarrowBand
:
    public neighbourNarrowBand
{
    //- Name of the registered geometric phase indicator field.
    word alphaName_;

    //- Jump / partial-fill tolerance on alpha.
    scalar alphaTol_;

public:

    TypeName("phaseIndicator");

    phaseIndicatorNarrowBand
    (
        const dictionary& dict,
        const volScalarField& psi
    );

    virtual ~phaseIndicatorNarrowBand() = default;

    //- Seed from snGrad(alpha) != 0 and from partially filled cells.
    //  Overrides the psi sign change; the dilation in
    //  neighbourNarrowBand::calc() is inherited unchanged.
    virtual void seed();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
