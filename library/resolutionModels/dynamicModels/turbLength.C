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

#include "turbLength.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbLength, 0);
addToRunTimeSelectionTable(dynamicModel, turbLength, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbLength::turbLength
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    dynamicModel(typeName, meshL, meshR),

    RResolvedAvg_
    (
        meshL.lookupObject<const volSymmTensorField>("RResolvedAvg")
    ),
    
    kSgsAvg_
    (
        meshL.lookupObject<const volScalarField>("kSgsAvg")
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
            "volLengthWeight", 0.5
        )
    ),
    adjustCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "adjustCoeff", paramDict_, 5.0
        )
    ),
    lower_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lower", paramDict_, 0.9
        )
    ),
    upper_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "upper", paramDict_, 1.1
        )
    ),
    epsSmall_
    (
        "epsSmall", dimVelocity*dimVelocity/dimTime, SMALL*10.0
    )
{
    assertOrder(lower_.value(), upper_.value());    
    computeEquivDelta();

    printParams();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbLength::computeEquivDelta()
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
    deltaL_.internalField() /= 1.5;
    deltaL_.internalField() *= (1-beta);
    deltaL_.internalField() += beta * pow(meshL_.V(), 0.333333);
}

  //- Print bount for resolution quantities
  void turbLength::printResolutionBounds()
  {

      if(printBounds_)
      {
          Info << "Turbulence quantity range: " << endl;
          Info << "    SGS TKE             : " 
              << "min = " 
              << Foam::gMin(kSgsAvg_.internalField()) << ", " 
              << "max = " 
              << Foam::gMax(kSgsAvg_.internalField()) << endl;
          
          Info << "    resolved TKE        : "
              << "min = " 
              << Foam::gMin(scalar(0.5) * 
              (tr(RResolvedAvg_.internalField()) )) << ", "
              << "max = " 
              << Foam::gMax(scalar(0.5) * 
              (tr(RResolvedAvg_.internalField()) ) )
              << endl;
          
          Info << "    total TKE           : " 
              << "min = " 
              << Foam::gMin(kTotal_.internalField()) << ", "
              << "max = " 
              << Foam::gMax(kTotal_.internalField()) << endl;

          Info << "    turb. length scale  : " 
              << "min = " 
              << Foam::gMin(turbLengthScale_.internalField())  << ", "
              << "max = "
              << Foam::gMax(turbLengthScale_.internalField()) << ", "
              << "avg = " << Foam::gAverage(turbLengthScale_.internalField())
              << endl;
      }
  }



void turbLength::evaluateFlag()
{
    kTotal_ = scalar(0.5) * tr(RResolvedAvg_) + kSgsAvg_;
    
    turbLengthScale_ = 
        Foam::sqrt(pow3(kTotal_)) / Foam::max(epsilonAvg_, epsSmall_); 

    resIndicator_.internalField() = 
        turbLengthScale_.internalField() / deltaL_.internalField();
    
    // Update lesFlag field with latency strategy
    lesFlagL_.internalField() += 
        pos
        (
            (
                turbLengthScale_.internalField()
                * ( lower_.value() / adjustCoeff_.value() )
            )
            - deltaL_.internalField() 
        ) * (1.0 - lesFlagL_.internalField())
        -
        pos
        ( 
            deltaL_.internalField() 
            -
            (
                turbLengthScale_.internalField()
                * ( upper_.value() / adjustCoeff_.value() )
            )
        ) * lesFlagL_.internalField();
    
    interpLesFlagL2R();
    printResolutionBounds();
}


bool turbLength::read()
{
    if (dynamicModel::read())
    {
        lower_.readIfPresent(paramDict());
        upper_.readIfPresent(paramDict());

        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dynamicModels
} // End namespace Foam

// ************************************************************************* //
