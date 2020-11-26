/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "muImpactSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muImpactSmagorinsky, 0);
addToRunTimeSelectionTable(muLESModel, muImpactSmagorinsky, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muImpactSmagorinsky::muImpactSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muLESModel(word("muImpactSmagorinsky"), U, phi, transport, 
               RAvg, epsilonAvg),
    nuSgs_(U.mesh().lookupObject<const volScalarField>("nuSgs")),
    B_(U.mesh().lookupObject<const volSymmTensorField>("B"))
{
//    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> muImpactSmagorinsky::devBeff() const
{
    FatalErrorIn("muImpactSSmagorinsky::devBeff")
        << "This function is dummy and should never be called."
        << exit(FatalError);
    
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> muImpactSmagorinsky::divDevBeff(volVectorField& U) const
{
    FatalErrorIn("muImpactSSmagorinsky::divDevBeff")
        << "This function is dummy and should never be called."
        << exit(FatalError);

    return
    (
      - fvm::laplacian(nuEff(), U) - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<volScalarField> muImpactSmagorinsky::epsilon() const
{
    return (  
        nuSgs() * magSqr(
            dev(
                twoSymm(
                    fvc::grad(U()
                    )
                )
            )
        ) 
    );
}


void muImpactSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    // muLESModel::correct(gradU);
}


bool muImpactSmagorinsky::read()
{
    if (muLESModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
