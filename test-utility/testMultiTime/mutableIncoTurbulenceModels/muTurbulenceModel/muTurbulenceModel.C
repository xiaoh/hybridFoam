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

#include "muTurbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(muTurbulenceModel, 0);
defineRunTimeSelectionTable(muTurbulenceModel, muTurbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muTurbulenceModel::muTurbulenceModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg,
    const volScalarField & lesFlag
)
:
    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    phi_(phi),
    transportModel_(lamTransportModel),
    RAvg_(RAvg),
    epsilonAvg_(epsilonAvg),
    lesFlag_(lesFlag)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<muTurbulenceModel> muTurbulenceModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg,
    const volScalarField & lesFlag
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
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("simulationType") >> modelName;
    }

    Info<< "Selecting turbulence model type " << modelName << endl;

    muTurbulenceModelConstructorTable::iterator cstrIter =
        muTurbulenceModelConstructorTablePtr_->find(modelName);

    if (cstrIter == muTurbulenceModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "muTurbulenceModel::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown muTurbulenceModel type " << modelName
            << endl << endl
            << "Valid muTurbulenceModel types are :" << endl
            << muTurbulenceModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<muTurbulenceModel>(cstrIter()(U, phi, transport,
                                                 RAvg, epsilonAvg, lesFlag));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void muTurbulenceModel::correct()
{
    transportModel_.correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
