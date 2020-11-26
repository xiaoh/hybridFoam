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

#include "LKE.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorList.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace passiveModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LKE, 0);
addToRunTimeSelectionTable(passiveModel, LKE, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LKE::LKE
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    passiveModel(typeName, meshL, meshR),

    runTime_(meshL.time()),

    RResolvedAvg_
    (
        meshL.lookupObject<const volSymmTensorField>("RResolvedAvg")
    ),
    
    kSgsAvg_
    (
        meshL.lookupObject<const volScalarField>("kSgsAvg")
    ),

    UAvg_
    (
        meshL.lookupObject<const volVectorField>("UAvg")
    ),

    epsilonAvg_
    (
        meshL.lookupObject<const volScalarField>("epsilonAvg")
    ),

    
    turbLengthScale_
    (
        IOobject
        (
            "turbLengthScale",
            runTime_.timeName(),
            meshL_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar
        ("turbLengthScale", dimLength, scalar(0.0))
    ),
    
    deltaL_
    (
        IOobject
        (
            "deltaL",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        turbLengthScale_
    ),
    
    kTotal_(scalar(0.5) * tr(RResolvedAvg_) + kSgsAvg_),
    
    volLengthWeight_
    (
        paramDict_.lookupOrAddDefault<scalar>
        (
            "volLengthWeight", 1.0
        )
    ),

    nuSgs_
    (
        meshL.lookupObject<const volScalarField>("nuSgs")
    ),

    relaxDict_
    (
        meshL.lookupObject<const IOdictionary>("relaxParameters")
    ),

    transportDict_
    (
        meshL.lookupObject<const IOdictionary>("transportProperties")
    ),

    averagingTime_
    (
        dimensioned<scalar>::lookupOrDefault
        ("averagingTime", relaxDict_.subDictPtr("timeScales"), 
        runTime_.endTime().value() * scalar(0.1),
        dimTime)
    ),
    
    nu_(transportDict_.lookup("nu")),

    epsilonMMAvg_
    (
        IOobject
     (
         "epsilonMMAvg",
         meshL.time().timeName(),
         meshL,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
        meshL,
        dimensionedScalar
        ("epsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0))
    ),

    resIndicatorL_
    (
        IOobject
        (
            "resIndicatorL",
            runTime_.timeName(),
            meshL,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar("resiL", dimless, scalar(0))
    ),

    resIndicatorK_
    (
        IOobject
        (
            "resIndicatorK",
            runTime_.timeName(),
            meshL,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar("resiK", dimless, scalar(0))
    ),

    resIndicatorE_
    (
        IOobject
        (
            "resIndicatorE",
            runTime_.timeName(),
            meshL,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar("resiE", dimless, scalar(0))
    ),

    epsSmall_
    (
        "epsSmall", dimVelocity*dimVelocity/dimTime, SMALL*10.0
    )

{
    computeEquivDelta();

    printParams();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LKE::computeEquivDelta()
{
    Info  <<  "Initialize according to mesh." << nl << endl;
    
    typedef List<vectorList> vectorListList;
    
    pointField points = meshL_.points();
    labelListList cellptLab = meshL_.cellPoints();
    
    vectorListList* vertexPosPtr = new vectorListList(cellptLab.size());
    vectorListList & vertexPos = * vertexPosPtr;
    
    forAll(cellptLab, cellI)
    {
        vertexPos[cellI].setSize(cellptLab[cellI].size()); 
        forAll(cellptLab[cellI], pI)
        {
            label p = cellptLab[cellI][pI];
            vertexPos[cellI][pI] = points[p];
        }
    }
    
    scalar tmpDist(-1);
    forAll(deltaL_, cellI)
    {
        scalar maxDist = 0;
        forAll(vertexPos[cellI], pI) 
        {
            forAll(vertexPos[cellI], pJ)
            {
                if 
                    ( pI!=pJ && 
                    (tmpDist=mag(vertexPos[cellI][pI]-vertexPos[cellI][pJ])) > maxDist 
                    )
                    maxDist = tmpDist;
            }
        }
        
        deltaL_.internalField()[cellI] = maxDist;
    }
    
    // blending with sqrt cell volume
    scalar beta = Foam::min(1, Foam::max(volLengthWeight_, 0));
    deltaL_.internalField() *= (1-beta);
    deltaL_.internalField() += beta * pow(meshL_.V(), 0.333333);
}


void LKE::evaluateFlag()
{
    kTotal_ = scalar(0.5) * tr(RResolvedAvg_) + kSgsAvg_;
    
    turbLengthScale_ = 
        Foam::sqrt(pow3(kTotal_)) / Foam::max(epsilonAvg_, epsSmall_); 

    // turbLength criteria
    resIndicatorL_.internalField() = 
        turbLengthScale_.internalField() / deltaL_.internalField();
    
    // TKE criteria (Pope)
    resIndicatorK_.internalField() = 
        1.0 - kSgsAvg_.internalField() / kTotal_.internalField();

    // Dissipation criteria (Davidson)
    dimensionedScalar existWt = 1.0/(1.0 + (runTime_.deltaT()/averagingTime_) );
    dimensionedScalar instantWt  = 1 - existWt;
    
    volSymmTensorField DAvg = dev(symm(fvc::grad(UAvg_)));
    
    epsilonMMAvg_ = existWt * epsilonMMAvg_
        + 2.0 * (nu_ + nuSgs_) * (magSqr(DAvg) * instantWt);

    resIndicatorE_.internalField() = 1.0 - 
        epsilonMMAvg_.internalField() / Foam::max(epsSmall_.value(), epsilonAvg_.internalField());

    lesFlagR_ = 0;
    lesFlagL_ = 1;
    
}


bool LKE::read()
{
    if (passiveModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace passiveModels
} // End namespace Foam

// ************************************************************************* //
