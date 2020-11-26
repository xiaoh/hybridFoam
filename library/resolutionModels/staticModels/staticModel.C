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

#include "staticModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(staticModel, 0);
defineRunTimeSelectionTable(staticModel, dictionary);
addToRunTimeSelectionTable(resolutionModel, staticModel, resolutionModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void staticModel::printParams()
{
    if (printParams_)
    {
        Info<< type() << "Params" << paramDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

staticModel::staticModel
(
    const word& type,
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    resolutionModel(meshL, meshR),

    IOdictionary
    (
        IOobject
        (
            "resolutionProperties",
            meshL.time().constant(),
            meshL.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    paramDict_(subDictPtr(type + "Params")),
    printParams_(lookupOrDefault<Switch>("printParams", false))
{

}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<staticModel> staticModel::New
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the resolutionModel is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "resolutionProperties",
                meshL.time().constant(),
                meshL.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("staticModelType") >> modelName;
    }

    Info<< "Selecting static resolution model: " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "staticModel::New(const fvMesh& meshL, const fvMesh& meshR)"
        )   << "Unknown staticModel type " << modelName
            << endl << endl
            << "Valid staticModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<staticModel>(cstrIter()(meshL, meshR));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool staticModel::read()
{
    if (regIOobject::read())
    {
        if (const dictionary* dictPtr = subDictPtr(type() + "Params"))
        {
            paramDict_ <<= *dictPtr;
        }
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
