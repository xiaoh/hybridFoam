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

#include "muLESModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muLESModel, 0);
defineRunTimeSelectionTable(muLESModel, dictionary);
addToRunTimeSelectionTable(muTurbulenceModel, muLESModel, muTurbulenceModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void muLESModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

muLESModel::muLESModel
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muTurbulenceModel(U, phi, lamTransportModel,
                      RAvg, epsilonAvg),

    IOdictionary
    (
        IOobject
        (
            "LESProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    //    coeffDict_(subOrEmptyDict(type + "Coeffs")),
    coeffDict_(subDictPtr(type + "Coeffs")),

    k0_("k0", dimVelocity*dimVelocity, SMALL),
    delta_(LESdelta::New("delta", U.mesh(), *this))
{
    readIfPresent("k0", k0_);

    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<muLESModel> muLESModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the muTurbulenceModel is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "LESProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("muLESModel") >> modelName;
    }

    Info<< "Selecting LES turbulence model " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "muLESModel::New(const volVectorField& U, const "
            "surfaceScalarField& phi, transportModel&)"
        )   << "Unknown muLESModel type " << modelName
            << endl << endl
            << "Valid muLESModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<muLESModel>(cstrIter()(U, phi, transport,
                                          RAvg, epsilonAvg));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void muLESModel::correct(const tmp<volTensorField>&)
{
    muTurbulenceModel::correct();
    delta_().correct();
}


void muLESModel::correct()
{
    correct(fvc::grad(U_));
}


bool muLESModel::read()
{
    if (regIOobject::read())
    {
        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        delta_().read(*this);

        readIfPresent("k0", k0_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
