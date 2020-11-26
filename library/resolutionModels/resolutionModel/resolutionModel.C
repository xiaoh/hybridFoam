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

#include "resolutionModel.H"
#include "zeroGradientFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(resolutionModel, 0);
defineRunTimeSelectionTable(resolutionModel, resolutionModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

resolutionModel::resolutionModel
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    runTime_(meshL.time()),
    meshL_(meshL),
    meshR_(meshR),
    nLocal_(meshL.nCells()),
    nGlobal_( returnReduce(nLocal_, sumOp<label>()) ),

    lesFlagL_
    (
        IOobject
        (
            "lesFlagL",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar("lesFlagL", dimless, scalar(1)),
        zeroGradientFvPatchScalarField::typeName 
    ),

    lesFlagR_
    (
        IOobject
        (
            "lesFlagR",
            runTime_.timeName(),
            meshR,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshR,
        dimensionedScalar("lesFlagR", dimless, scalar(0)),
        zeroGradientFvPatchScalarField::typeName 
    )

{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<resolutionModel> resolutionModel::New
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

        dict.lookup("resolutionType") >> modelName;
    }

    Info<< "Selecting resolution model type: " << modelName << endl;

    resolutionModelConstructorTable::iterator cstrIter =
        resolutionModelConstructorTablePtr_->find(modelName);

    if (cstrIter == resolutionModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "resolutionModel::New(const fvMesh&, const fvMesh&) "
        )   << "Unknown resolutionModel type " << modelName
            << endl << endl
            << "Valid resolutionModel types are :" << endl
            << resolutionModelConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<resolutionModel>(cstrIter()(meshL, meshR));
}


// Interpolate lesFlag from RANS mesh to LES mesh
void resolutionModel::interpLesFlagR2L()
{
    meshToMesh R2LInterp(meshR_, meshL_, HashTable<word>(0), wordList());
    
    R2LInterp.interpolateInternalField
        (
            lesFlagL_, lesFlagR_, 
            // meshToMesh::CELL_POINT_INTERPOLATE
            meshToMesh::INTERPOLATE,
	    eqOp<scalar>()	    
        );
    
    lesFlagL_ = pos(lesFlagL_ - 0.5);
}


// Interpolate lesFlag from RANS mesh to LES mesh
void resolutionModel::interpLesFlagL2R()
{
    meshToMesh L2RInterp(meshL_, meshR_, HashTable<word>(0), wordList());
    
    L2RInterp.interpolateInternalField
        (
            lesFlagR_, lesFlagL_, 
            meshToMesh::CELL_POINT_INTERPOLATE,
	    eqOp<scalar>()	    
        );
    
    lesFlagR_ = pos(lesFlagR_ - 0.5);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void resolutionModel::diagnosisFlagRegion()
{
    scalar flagPercent = gAverage(lesFlagL_) * 100.0;
    scalarField flagedVol = lesFlagL_.internalField() * meshL_.V();
    scalar flagPercentVol = gSum(flagedVol) / gSum(meshL_.V()) * 100.0;
    
    Info << "LES flagged region: "
        << flagPercent << "% cells; " 
        << flagPercentVol << "% volume;" << endl
        << endl;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
