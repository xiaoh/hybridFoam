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

#include "resolvedEPS.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(resolvedEPS, 0);
addToRunTimeSelectionTable(dynamicModel, resolvedEPS, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

resolvedEPS::resolvedEPS
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    dynamicModel(typeName, meshL,  meshR),

    UAvg_
    (
        meshL.lookupObject<const volVectorField>("UAvg")
    ),

    epsilonAvg_
    (
        meshL.lookupObject<const volScalarField>("epsilonAvg")
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

    runTime_(meshL.time()),

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
    
    lower_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "lower", paramDict_, 0.78
        )
    ),
    upper_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "upper", paramDict_, 0.82
        )
    ),
    epsSmall_
    (
        "epsSmall", dimVelocity*dimVelocity/dimTime, SMALL*10.0
    )
{
    assertOrder(lower_.value(), upper_.value());    

    printParams();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  //- Print bount for resolution quantities
  void resolvedEPS::printResolutionBounds()
  {

      if(printBounds_)
      {
          Info << "Turbulence quantity range: " << endl;
          Info << "Mean Motion Dissipation             : " 
              << "min = " 
              << Foam::gMin(epsilonMMAvg_.internalField()) << ", " 
              << "max = " 
              << Foam::gMax(epsilonMMAvg_.internalField()) << endl;
          
          Info << "Total dissipation             : "
               << "min = " 
               << Foam::gMin(epsilonAvg_.internalField()) << ", "
               << "max = " 
               << Foam::gMax(epsilonAvg_.internalField())
               << endl;
      }
  }



void resolvedEPS::evaluateFlag()
{
    dimensionedScalar existWt = 1.0/(1.0 + (runTime_.deltaT()/averagingTime_) );
    dimensionedScalar instantWt  = 1 - existWt;
    
    volSymmTensorField DAvg = dev(symm(fvc::grad(UAvg_)));
    
    epsilonMMAvg_ = existWt * epsilonMMAvg_
        + 2.0 * (nu_ + nuSgs_) * (magSqr(DAvg) * instantWt);

    resIndicator_.internalField() = 1.0 - 
        epsilonMMAvg_.internalField() / Foam::max(epsSmall_.value(), epsilonAvg_.internalField());
      
    // Update lesFlagEps_ using latency strategy.
    lesFlagL_.internalField() += 
        (
            pos(resIndicator_.internalField() - upper_.value()) * (1.0 - lesFlagL_.internalField()) 
            - pos(lower_.value() - resIndicator_.internalField()) * lesFlagL_.internalField() 
        ); 
    
    interpLesFlagL2R();
    
    printResolutionBounds();
}


bool resolvedEPS::read()
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
