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

#include "turbLengthAniso.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbLengthAniso, 0);
addToRunTimeSelectionTable(dynamicModel, turbLengthAniso, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbLengthAniso::turbLengthAniso
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
        dimensionedVector
        ("turbLengthScale", dimLength, vector::zero)
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
    
    kTotal_
    ( 
        IOobject
        (
            "kTotal",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedVector
        (
            "turbLengthScale", dimensionSet(0, 3, -3, 0, 0), vector::zero
        )
    ),

    resIndX_
    ( 
        IOobject
        (
            "resIndX",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar
        (
            "resx", dimless, 0.0
        )
    ),

    resIndY_
    ( 
        IOobject
        (
            "resIndY",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar
        (
            "resy", dimless, 0.0
        )
    ),

    resIndZ_
    ( 
        IOobject
        (
            "resIndZ",
            runTime_.timeName(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshL,
        dimensionedScalar
        (
            "resy", dimless, 0.0
        )
    ),
    
    adjustCoeff_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "adjustCoeff", paramDict_, 5.0
        )
    ),
    scaleTurbLengthBy_
    (
        dimensioned<vector>::lookupOrAddToDict
        (
            "scaleTurbLengthBy", paramDict_, vector(1.0, 1.0, 1.0)
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

	forAll(kTotal_, cellI)
	{
        kTotal_[cellI].component(vector::X) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::XX) + kSgsAvg_[cellI]/3.0));
        kTotal_[cellI].component(vector::Y) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::YY) + kSgsAvg_[cellI]/3.0));
        kTotal_[cellI].component(vector::Z) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::ZZ) + kSgsAvg_[cellI]/3.0));
	}

    assertOrder(lower_.value(), upper_.value());    
    computeEquivDelta();

    printParams();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbLengthAniso::computeEquivDelta()
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
        vector maxVec = vector::zero;
        forAll(vertexPos[cellI], pI) 
        {
            forAll(vertexPos[cellI], pJ)
            {
                if 
                    ( pI!=pJ && 
                    (tmpDist=mag(vertexPos[cellI][pI]-vertexPos[cellI][pJ])) > maxDist 
                    )
                {
                    maxDist = tmpDist;
                    maxVec = vertexPos[cellI][pI]-vertexPos[cellI][pJ];
                }
                   
            }
        }

        deltaL_.internalField()[cellI] = maxVec;
    }

}

  //- Print bount for resolution quantities
  void turbLengthAniso::printResolutionBounds()
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



void turbLengthAniso::evaluateFlag()
{

    forAll(kTotal_, cellI)
      {
        kTotal_[cellI].component(vector::X) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::XX) + kSgsAvg_[cellI]/3.0));
        kTotal_[cellI].component(vector::Y) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::YY) + kSgsAvg_[cellI]/3.0));
        kTotal_[cellI].component(vector::Z) =  Foam::sqrt(pow3( 0.5*RResolvedAvg_[cellI].component(symmTensor::ZZ) + kSgsAvg_[cellI]/3.0));
      }
  
    turbLengthScale_ = kTotal_ / Foam::max(epsilonAvg_/3.0, epsSmall_); 
  
    resIndX_.internalField() = scaleTurbLengthBy_.value().x() *
      turbLengthScale_.internalField().component(vector::X) / mag(deltaL_.internalField().component(vector::X));
  
    resIndY_.internalField() = scaleTurbLengthBy_.value().y() *
      turbLengthScale_.internalField().component(vector::Y) / mag(deltaL_.internalField().component(vector::Y));
  
    resIndZ_.internalField() = scaleTurbLengthBy_.value().z() *
      turbLengthScale_.internalField().component(vector::Z) / mag(deltaL_.internalField().component(vector::Z));

    forAll(resIndicator_, cellI)
      {
        resIndicator_[cellI] = min(min(resIndX_[cellI], resIndY_[cellI]), resIndZ_[cellI]);
      }

    // Update lesFlag field with latency strategy
    lesFlagL_.internalField() += 
        pos
        (
            (
                resIndicator_.internalField()
                * ( lower_.value() / adjustCoeff_.value() )
            )
            - 1.0
        ) * (1.0 - lesFlagL_.internalField())
        -
        pos
        ( 
            1.0
            -
            (
                resIndicator_.internalField()
                * ( upper_.value() / adjustCoeff_.value() )
            )
        ) * lesFlagL_.internalField();

    // Interpolate to lesFlagR
    interpLesFlagL2R();
    
    printResolutionBounds();
}


bool turbLengthAniso::read()
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
