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

#include "muImpactADMRT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muImpactADMRT, 0);
addToRunTimeSelectionTable(muLESModel, muImpactADMRT, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muImpactADMRT::muImpactADMRT
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muLESModel(word("muImpactADMRT"), U, phi, transport, 
               RAvg, epsilonAvg),
    epsilonADM_
    (
        IOobject
        (
            "epsilonADM",
            runTime_.timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        ("epsilonADM", sqr(dimLength)/pow3(dimTime), SMALL)
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        k0()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// dummy, will never be called
tmp<volSymmTensorField> muImpactADMRT::devBeff() const
{
    FatalErrorIn("muImpactADMRT::devBeff")
        << "This function is dummy and should never be called."
        << exit(FatalError);

    return -nuEff()*dev(twoSymm(fvc::grad(U()))); 
}


// dummy, will never be called
tmp<fvVectorMatrix> muImpactADMRT::divDevBeff(volVectorField& U) const
{
    FatalErrorIn("muImpactADMRT::divDevBeff")
        << "This function is dummy and should never be called."
        << exit(FatalError);
    
    return
        (
            - fvm::laplacian(nuEff(), U) - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
        );
}


void muImpactADMRT::correct(const tmp<volTensorField>& gradU)
{}


bool muImpactADMRT::read()
{
    if (muLESModel::read())
    {
        // ck_.readIfPresent(coeffDict());
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
