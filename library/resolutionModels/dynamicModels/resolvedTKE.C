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

#include "resolvedTKE.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(resolvedTKE, 0);
addToRunTimeSelectionTable(dynamicModel, resolvedTKE, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

resolvedTKE::resolvedTKE
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    dynamicModel(typeName, meshL,  meshR),

    RResolvedAvg_
    (
        meshL.lookupObject<const volSymmTensorField>("RResolvedAvg")
    ),
    
    kSgsAvg_
    (
        meshL.lookupObject<const volScalarField>("kSgsAvg")
    ),
    
    kTotal_(scalar(0.5) * tr(RResolvedAvg_) + kSgsAvg_),
    
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
  void resolvedTKE::printResolutionBounds()
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
      }
  }



void resolvedTKE::evaluateFlag()
{
    kTotal_ = scalar(0.5) * tr(RResolvedAvg_) + kSgsAvg_;

    resIndicator_.internalField() = 
        1.0 - kSgsAvg_.internalField() / kTotal_.internalField();
    
    // Update lesFlag field with latency strategy
    lesFlagL_.internalField() +=
        (
            pos(resIndicator_.internalField() - upper_.value()) * (1.0 - lesFlagL_.internalField() ) 
            - pos(lower_.value() - resIndicator_.internalField()) * lesFlagL_.internalField()
        );  // Validated.
    
    interpLesFlagL2R();
    
    printResolutionBounds();
}


bool resolvedTKE::read()
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
