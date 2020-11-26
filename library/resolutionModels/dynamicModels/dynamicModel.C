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

#include "dynamicModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicModel, 0);
defineRunTimeSelectionTable(dynamicModel, dictionary);
addToRunTimeSelectionTable(resolutionModel, dynamicModel, resolutionModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void dynamicModel::printParams()
{
    if (printParams_)
    {
        Info<< type() << "Params" << paramDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicModel::dynamicModel
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
    printParams_(lookupOrDefault<Switch>("printParams", false)),
    printBounds_(lookupOrDefault<Switch>("printBounds", false)),
    flagInterpolate_(meshL),
    distToWall_(meshL),

    resIndicator_
    (
        IOobject
        (
            "resIndicator",
            runTime_.timeName(),
            meshL,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar("resIndicator", dimless, scalar(1))
    )
{

}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<dynamicModel> dynamicModel::New
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

        dict.lookup("dynamicModelType") >> modelName;
    }

    Info<< "Selecting dynamic resolution model: " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dynamicModel::New(const fvMesh& meshL, const fvMesh& meshR)"
        )   << "Unknown dynamicModel type " << modelName
            << endl << endl
            << "Valid dynamicModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<dynamicModel>(cstrIter()(meshL, meshR));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicModel::assertOrder(scalar& l, scalar& u)
{
    if (l > u)
    {
        scalar tmp = u;
        u = l;
        l = tmp;
    }
}

  
// Assert Flags are either 1 or 0
void dynamicModel::assertFlagInRange(volScalarField& lesFlagX)
{
    word name =  lesFlagX.name();
    scalar small(1e-4);
    forAll(lesFlagX, celli)
    {      
        scalar flag = lesFlagX[celli];
        if( ! (mag(flag) < small || mag(flag - 1.0) < small) )
            FatalErrorIn
                (
                    "dynamicModel::assertFlagInRange()"
                )
                << name << " should be either 0 or 1!" << endl
                    << "Debug Info: " 
                    << name << "[" << celli << "] =" << flag
                    << abort(FatalError);
    }
}


void dynamicModel::diagnosisFlagConnectivity()
    // Warning: No. of interFaces not accurate for parallel runs.
{
    // RANS/LES Interface distribution
    surfaceScalarField interfaceLR = flagInterpolate_.interpolate(lesFlagL_);
    surfaceScalarField interFaceFlag = 1.0 - Foam::mag(interfaceLR - 0.5) * 2.0;
    
    scalar nLesCell = gSum(lesFlagL_.internalField()); // Global No. LES Cells
    scalar nRansCell = nGlobal_ - nLesCell; // Global No. RANS Cells
    scalar nMinorCell = min(nLesCell, nRansCell);
    scalar interFacePerMC = gSum(interFaceFlag.internalField()) / max(nMinorCell, 1);
    
    // Distance of rans cells to walls
    scalarField RCellToWall = (distToWall_.y() * ( 1 - lesFlagL_));
    
    // Mean Distance (weighted by cell volume) of RANS cells away from the walls
    scalar meanRCWD = 
        gSum(RCellToWall * meshL_.V()) 
        / max(gSum(meshL_.V()*( 1 - lesFlagL_)), SMALL*10.0);
    
    word warn("APPROXIMATE");
    if (! Pstream::parRun()) warn="";
    
    Info << "Connectivity: " 
        << warn << " no. total interface = " << interFacePerMC * nMinorCell << ", "
        << "interface/minority cell = " << interFacePerMC << "; "
        << "Avg wallDist RANS zone = " << meanRCWD << endl;
}



bool dynamicModel::read()
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
